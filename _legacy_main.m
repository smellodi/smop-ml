% The ML process using DIFFERENTIAL EVOLUTION algorithm
% Described in Mariaana's manuscript.

%% NOTES
%
% IF WE WANT TO INCLUDE PID, THEN THE WHOLE SYSTEM NEEDS TO BE REDESIGNED
%
% Philipp Muller, January 2024  
%
% IMPORTANT: Recipe name is parsed in the SMOP app, so do not change the
% way it is constructed in this script.
%
% Oleg Spakov, March 2024

%% APP ENTRY
function main(varargin)
    args = argparser(varargin,nargin);
    if size(args) == 0
        return
    end

    % dbstop error
    % rng("shuffle")

    %% STEP 1: Establish a connection to SMOP and get config parameters.
    
    % create a TCP client and connect to SMOP software
    fprintf("Connecting to %s ...", args.ip);
    smopClient = SmopClient(args.ip);
    
    % continue only if the connection is OK
    if ~smopClient.isConnected
        fprintf("  Failed. Exiting...\n");
        pause(3);
        return;
    end
    
    fprintf("  Done, port = %d\n", smopClient.port);
    
    % continue only if config packet was received
    if ~smopClient.readConfig()
        fprintf("Exiting...\n");
        pause(2);
        return;
    end
    
    %% STEP 2: Define parameters
    
    gas = struct(); % stores gas info and constant parameters

    % Minimum and maximum flows of chemicals in the mixture (in sccm)
    gas.min = 0;
    gas.max = smopClient.getMaxFlow(50);
    
    % Maximum number of iterations
    args.mi = smopClient.getMaxIterationCount(args.mi);
    
    % Distance threshold for terminating iterative process
    args.th = smopClient.getThreshold(args.th);
    
    % Define which distance measure to use
    args.alg = smopClient.getAlgorithm(args.alg);
    
    % Size of a vector (number of chemicals)
    gas.n = smopClient.getChannelCount(1);
    
    fprintf("%-20s n = %d, min = %d ccm, max = %d ccm\n", ...
        "Gas parameters:", gas.n, gas.min, gas.max);
    fprintf("%-20s name = %s, cr = %.2f, f = %.2f\n", ...
        "Alg parameters:", args.alg, args.cr, args.f);
    fprintf("%-20s mi = %d, th = %.2f\n", ...
        "Search parameters:", args.mi, args.th);

    % Randomization offset. This is max offset applied to the originally 
    % generated flow value if this flow does not pass the validation process
    % (used in "validate" function)
    gas.df = 0.2 * (gas.max - gas.min);

    % Limits are used to detect if the generated flow values may cause
    % measurement sensor oversaturation. 
    % To be used in "limitVectors" function. Note that this function 
    % operates with 2-3 gases only, thus only 3 limits are used
    % (no limitation is applied if there are 4 or more gases)
    gas.limits = smopClient.getCriticalFlows();

    % Gas names often used to print out info
    gas.names = smopClient.gases;

    %% STEP 3: Read in message with DMS measurement from target scent
    
    initM = smopClient.getInitialMeasurement();

    if isstruct(initM.dms)
        initM.dms = initM.dms.data.positive;
        
        % Separation voltage used in DMS scope mode
        usv = smopClient.getUsv();

        fprintf("DMS scan parameters: ");
        if usv > 0
            fprintf("Us = %.1f V\n", usv);
        else
            fprintf("full scan\n");
        end
    else
        initM.dms = [];
    end

    if isstruct(initM.snt)
        initM.snt = initM.snt.data.resistances;
    else
        initM.snt = [];
    end

    %% STEP 4: Initialization
    
    % Start without any prior information.
    
    %% STEP 4a: Init vectors and consts
    % All measured vectors (flow rates)
    F = createInitialVectors(gas.n,gas.min,gas.max);

    % Both flows set to gas.max result in oversaturated PID, therefore 
    % we apply some limitations as described in limitVectors.
    F = limitVectors(F,gas.limits);
    F = roundTo(F,args.dec);

    % Population size: set of best vectors
    % Typically 10 times the dimension n of x, but this may be too long... 
    % ..keep it same as the number of initial vectors.
    X = F;
    
    XI = 1:size(X,2); % range of X and U
    
    %% STEP 4b: Collect measurements

    % All measurements
    M = arrayfun(@(x)struct("dms",0,"snt",0),XI);

    gm = 1e8;   % overall minimum distance (global minima)
    gmId = -1;  % M and F index of the global minima;
    cfm = 1e8;  % last search distance (measured trials only)
    
    fprintf("\nCollecting initial measurements:\n");
    for jj = XI
        pause(0.1); % simply, to make ML being not too fast in SMOP interface
    
        recipeName = sprintf("Reference #%d", jj);
        smopClient.sendRecipe(recipeName,F(:,jj),false,cfm);
        clear recipeName;

        fprintf("[%d] %s", jj, formatVector(gas.names,F,jj));
    
        M(jj) = getMeasurement(smopClient);

        cfm = getSimilarityMeasure(args.alg,initM,M(jj));
        if (cfm < gm)
            gm = cfm;
            gmId = jj;
        end
        fprintf(" DIST=%6.3f\n", cfm);
    end

    fprintf("GM: %.4f [%s]\n", gm, formatVector(gas.names,F,gmId));
    
    %% STEP 5: Iterative step

    iter = 1;   % iteration counter
    MI = XI;    % indexes of M corresponding to vectors consisting X

    isFinished = false;     % switch for terminating the iterative process

    while (~isFinished)
        fprintf("\nIteration #%d:\n", iter);

        im = 1e8;   % minimum distance for vectors tested in this iteration
        imId = -1;  % U index that corresponds to 'im' (iteration minima)
    
        %% STEP 5a: Differential evolution
        % [https://en.wikipedia.org/wiki/Differential_evolution]

        V = mutate(X,args.f);               % generate new vectors
        V = limitValues(V,gas.min,gas.max); % limit their values 
        U = crossover(X,V,args.cr);         % mix old and new vectors
        U = validate(U,gas.df);             % remove repetitions
        U = limitValues(U,gas.min,gas.max); % limit values after rndmzation
        U = limitVectors(U,gas.limits);     % avoid oversaturation
        U = roundTo(U,args.dec);            % round flow values

        fprintf("Flows to test:\n%s", formatVectorsAll(gas.names,U));
        clear V;
        
        %% STEP 5b: Selecting better vectors

        % Compute distances from the initial measurements and measureemnts 
        % of trial vectors and corresponding target vectors
    
        for jj = XI
            %% STEP 5b1: Compute distance of TARGET vector

            cfX = getSimilarityMeasure(args.alg,initM,M(MI(jj)));

            %% STEP 5b2: Compute distance of TRIAL vector
    
            % Find shuffled intertersection of indexes for which F and U
            % have same flow rates
            C = getCommonIndices(F,U,jj);
    
            if ~isempty(C)  % a measurement of this vector exists
                kk = C(1);  % already, lets take it from the database M
                fprintf("[%d] REPEAT  %s", jj, formatVector(gas.names,F,kk));
            else
                pause(0.1);     % emulate some heavy ML search :)
    
                % This vector was not measured yet, so lets do it
                recipeName = sprintf("Iteration #%d, Search #%d", iter, jj); 
                smopClient.sendRecipe(recipeName,U(:,jj),isFinished,cfm);
                clear recipeName;
    
                % Wait for new measurement and add it to the table of 
                % measured recipes
                fprintf("[%d] MEASURE %s", jj, formatVector(gas.names,U,jj));
                kk = length(M) + 1;
                M(kk) = getMeasurement(smopClient);
                F(:,kk) = U(:,jj);
            end

            cfU = getSimilarityMeasure(args.alg,initM,M(kk));
            fprintf(" DIST=%6.3f", cfU);

            if isempty(C)   % memorize cf of measured data for recipe
                cfm = cfU;
            end

            clear C;
            
            %% STEP 5b3: update distance minimas

            cf = min(cfU,cfX);
            info = ["  ", "  ", ""];

            % Global minima
            if cf < gm
                gm = cf;
                gmId = kk;
                info(1) = "GM";
            end

            % Minima of the tested vectors
            if cfU < im
                im = cfU;
                imId = jj;
                info(2) = "IM";
            end
    
            %% STEP 5b4: replace target vec with trial vec if it has lower cf

            if cfU < cfX
                X(:,jj) = U(:,jj);
                info(3) = sprintf("[%d >> %d]", MI(jj), kk);
                MI(jj) = kk;
            end
    
            fprintf(" %s\n", join(info));
            clear cfU cfX cf kk info;
        end
        
        fprintf("IM: %.4f [%s]\n", im, formatVector(gas.names,U,imId));
        fprintf("GM: %.4f [%s]\n", gm, formatVector(gas.names,F,gmId));
        
        %% STEP 5c: Make a decision about the proximity of the best guess
    
        if (gm < args.th)
            isFinished = true;
            recipeName = "Final recipe";
            fprintf("\n%s:\n", recipeName);
        elseif (iter >= args.mi)
            isFinished = true;
            recipeName = sprintf("Best after %d iterations", args.mi);
            fprintf("\n%s:\n", recipeName);
        end
    
        if isFinished
            % Send the final recipe
            flows = F(:,gmId);
            smopClient.sendRecipe(recipeName, flows, isFinished, cfm);
            fprintf("  %s, DIST = %.4f\n\nFinished\n\n", ...
                formatVector(gas.names,F,gmId), gm);
            clear recipeName flows;
        else
            fprintf("Continuing the search, best vectors are:\n%s", ...
                formatVectorsAll(gas.names,X));
        end
        
        iter = iter + 1;
    end
    
    %% STEP 6: Finalize
    
    clear smopClient;   % disconnect from the SMOP
    pause(3);           % short pause before closing the CMD window

end

%%% FUNCTIONS

% Waits for new measuremenst to arrive and return the extracted data
function m = getMeasurement(smopClient)
    m = struct("dms",[],"snt",[]);

    newM = smopClient.waitForMeasurement();
    if isstruct(newM.dms)
        m.dms = newM.dms.data.positive;
    end
    if isstruct(newM.snt)
        m.snt = newM.snt.data.resistances;
    end
end

% Creates initial vectors each of n size. The result includes all possible
% combinations of min_ and max_, plus the vector with central values.
function F = createInitialVectors(n, min_, max_)
    V = [];
    for jj = n:-1:0
        a = repmat(min_,1,jj);      % [0, 0, ...  (or empty list if jj == 0)
        b = repmat(max_,1,n-jj);    % ... 50, 50] (or empty list if jj == n)
        U = [a b];                  % unite both parts
        U = perms(U);               % save all permutations as rows
        V = [V;U];                  % add to the resulting list
    end
    V = unique(V,"rows")';          % remove duplications produced by
                                    % permutations, and transpose the matrix
    center = (max_ + min_) / 2;
    F = [V repmat(center,n,1)];     % add central point

    while size(F) < 4               % in case n = 1 that gives only 3 points
        F = [F randi(max_)];        % add a random point
    end
end

% Finds common indices and returns a permutated array of them.
function M = getCommonIndices(V, U, ucol)
    [n,c] = size(V);

    % The initial intersection array stores all indices
    M = 1:c;

    % Find indices of samples for which V and U have same flow rates
    for jj = 1:n
        ID = find(V(jj,:) == U(jj,ucol));
        M = intersect(ID,M);
        if isempty(M)
            return      % the result is empty, no need to continue
        end
    end
    
    % Shuffle the resulting array in random order
    M = M(randperm(numel(M)));
end

% A cost function that measures similarity (distance) between dispersion plot 
% of the mixture we want to recreate and the training dispesion plot. 
function distance = getSimilarityMeasure(alg, measrm1, measrm2)
    distance = 1e8;

    if (alg == "euclidean")
        % For Euclidean algorithm we use RMSE as the distance.
        if ~isempty(measrm2.dms)
            distance = sqrt(mean((measrm1.dms - measrm2.dms).^2));
        elseif ~isempty(measrm2.snt)
            distance = sqrt(mean((measrm1.snt - measrm2.snt).^2));
            if distance > 10000            % huge values observed with real SNT
                distance = distance / 100000;  % TODO: remove in future
            end
        end
    end
end

% Mutation
function V = mutate(X, f)
    [n,np] = size(X);
    V = nan(n,np);

    for jj = 1:np
        % Pick three distinct vectors from X that are different from jj
        ids = randperm(np);    % random permutation of integers 1 to np
        ids(ids == jj) = [];   % remove integer jj
        
        % Compute donor vector from first three vectors
        % pick three distinct dispersion plots (also PID, SNT) 
        % (need to be distinct from each other and from x)
        value = X(:,ids(1)) + f * (X(:,ids(2)) - X(:,ids(3)));

        V(:,jj) = value;
    end
end

% Uses binomial method for combining components from target and donor vectors
function U = crossover(X, V, cr)
    [n,np] = size(X);

    % Generate randomly chosen indices to ensure that at least one
    % component of the donor vector is included in the target vector
    rv = randi(n,1,np);

    % Random numbers in [0,1] for each component of each target vector
    rv2 = rand(n,np);

    % Combining target and donor vectors to get trial vectors
    U = nan(n,np);
    for jj = 1:np
        for kk = 1:n
            if (rv2(kk,jj) <= cr) || (rv(jj) == kk)
                U(kk,jj) = V(kk,jj);
            else                        % OLEG: extreamly rare if cr=0.8, but
                U(kk,jj) = X(kk,jj);    % happens in 10-30% cases when cr=0.5
            end
        end
    end
end

% Avoid the search algorithm to stuck with testing same or similar vectors
function V = validate(V, delta)
    [n,np] = size(V);
    interval = [-delta, delta];

    % Replace repeated pairs with random pairs
    for jj = 2:np
        for kk = 1:(jj-1)
            if V(:,jj) == V(:,kk)
                V(:,kk) = V(:,kk) + randi(interval,n,1);
            end
        end
    end

    % If all flow values of a certain gas are the same ..
    allSame = arrayfun(@(x) true, V(:,1));
    for jj = 2:np
        for kk = 1:n
            allSame(kk) = allSame(kk) && V(kk,jj-1) == V(kk,jj);
        end
    end

    % .. then randomize those values a bit
    for kk = 1:n
        if allSame(kk)
            V(kk,:) = V(kk,:) + randi(interval,1,np);
        end
    end
end

% Limits the values to stay within bounds
function V = limitValues(V, min_, max_)
    V(V < min_) = min_;
    V(V > max_) = max_;
end

% Avoids measurement sensor oversaturation by limiting a combined gas flow.
% We assume that limits are semi-axis (a,b) of an ellipse, and if the
% vector lies outside of this ellipse, then it is "pulled" toward the
% center to the ellipse edge.
% NOTE: will be applied only for 2 or 3 -dimensional vectors (2-3 gases)
function V = limitVectors(V, limits)
    [n,np] = size(V);
    if n == 2
        for jj = 1:np
            if (V(1,jj) > 0) && (V(2,jj) > 0)
                a = atan2(V(2,jj),V(1,jj));
                lf = [limits(1) * cos(a); limits(2) * sin(a)];
                V(:,jj) = min(lf,V(:,jj));
            end
        end
    elseif n == 3
        for jj = 1:np
            if (V(1,jj) > 0) && (V(2,jj) > 0) && (V(3,jj) > 0)
                a = atan2(V(2,jj),V(1,jj));
                b = atan2(V(3,jj),V(1,jj));
                lf = [limits(1) * cos(a) * cos(b); ...
                      limits(2) * sin(a); ...
                      limits(3) * sin(b)];
                V(:,jj) = min(lf,V(:,jj));
            end
        end
    end
end

% Rounds values in the vector.
% If dec < 0, then rounds the values to be divisible by |dec|.
function V = roundTo(V, dec)
    r = max(dec,0);
    V = round(V,r);

    if dec < 0
        [n,np] = size(V);
        dec = abs(dec) + 1;
        for kk = 1:n
            for jj = 1:np
                r = round(V(kk,jj) / dec);
                V(kk,jj) = dec * r;
            end
        end
    end
end

% Formats gas names N with their flows V as "GAS1=FLOW1 GAS2=FLOW2 ..."
function s = formatVector(N, V, vcol)
    s = "";
    for ii = 1:size(V,1)
        if ii ~= 1
            s = sprintf("%s ", s);
        end
        s = sprintf("%s%s=%4.1f", s, N(ii), V(ii,vcol));
    end
end

% Formats gas names N with a list of flows V as
%   GAS1 FLOW1 FLOW2 ..
%   GAS2 FLOW1 FLOW2 ..
%   ..
function s = formatVectorsAll(N, V)
    s = "";
    [n,np] = size(V);
    for jj = 1:n
        s = sprintf("%s  %-12s", s, N(jj));
        for kk = 1:np
            s = sprintf("%s %4.1f", s, V(jj,kk));
        end
        s = sprintf("%s\n", s);
    end
end
