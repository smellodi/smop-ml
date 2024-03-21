% The learning process uses DIFFERENTIAL EVOLUTION algorithm described in 
% Mariaana's manuscript.


%% NOTES

% SUGGESTION: USE FLOW RATES OF 20 OR HIGHER FOR TARGET SCENT TO SPEED
% UP RESPSONSE TIME/SHORTEN STABILIZATION TIME
% 
% SUGGESTION: MEASURE ONLY REGIONS WITH ALPHA CURVES FROM ANTON'S PAPER:
%   SV 450 to 700, CV -1 to 4
%
% IF WE WANT TO INCLUDE PID, THEN THE WHOLE SYSTEM NEEDS TO BE REDESIGNED
%
% Philipp Muller, January 2024  

% IMPORTANT: Recipe name is parsed in the SMOP app, so do not change the
% way it is constructed in this script.
%
% Oleg Spakov, March 2024

%% APP ENTRY
function main(varargin)
    args = argparser(varargin,nargin);

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
    
    % RMSE threshold for terminating iterative process
    args.th = smopClient.getThreshold(args.th);
    
    % Define which distance measure to use
    args.alg = smopClient.getAlgorithm(args.alg);
    
    % Size of a vector (number of chemicals)
    gas.n = smopClient.getChannelCount(2);
    
    fprintf("%-20s n = %d, min = %d ccm, max = %d ccm\n", ...
        "Gas params:", gas.n, gas.min, gas.max);
    fprintf("%-20s name = %s, cr = %.2f, f = %.2f\n", ...
        "Alg params:", args.alg, args.cr, args.f);
    fprintf("%-20s mi = %d, th = %.2f\n", ...
        "Search params:", args.mi, args.th);

    % Randomization offset. This is max offset applied to the originally 
    % generated flow value if this flow does not pass the validation process
    % (used in "validate" function)
    gas.df = 0.2 * (gas.max - gas.min);

    % Limits are used to detect if the generated flow values may cause
    % measurement sensor oversaturation.
    % To be used in "limitVectors" function.
    gas.limits = smopClient.getCriticalFlow([55; 60]);

    % Gas names often used to print out info
    gas.names = smopClient.gases;

    %% STEP 3: Read in message with DMS measurement from target scent
    
    % Separation voltage used in DMS scope mode (i.e. args.ssv == 1)
    usv = 0;
    
    initM = smopClient.getInitialMeasurement();

    if isstruct(initM.dms)
        if args.ssv
            [usv, initDMS] = getDmsSingleLine(initM.dms,0.5);
            fprintf("DMS scan parameters: Us = %.1f V\n", usv);
            initM.dms = initDMS;
            clear initDMS;
        else   % use full DMS scan
            usv = 0;
            initM.dms = initM.dms.data.positive;
    
            fprintf("DMS scan parameters: full scan\n");
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
    %   we apply some limitations as described in limitVectors.
    F = limitVectors(F,gas.limits);
    F = roundTo(F,args.dec);

    % Population size: set of best vectors
    % Typically 10 times the dimension n of x, but this may be too long... 
    %    keep it same as the number of initial vectors.
    X = F;
    
    XI = 1:size(X,2); % range of X and U
    
    %% STEP 4b: Collect measurements

    % All measurements
    M = arrayfun(@(x)struct("dms",0,"snt",0),XI);

    gm = 1e8;   % overall minimum RMSE (global minima)
    gmId = -1;  % M and F index of the global minima;
    cfm = 1e8;  % last search RMSE (measured trials only)
    
    fprintf("\nCollecting initial measurements:\n");
    for jj = XI
        pause(0.5); % simply, to make ML being not too fast in SMOP interface
    
        recipeName = sprintf("Reference #%d", jj);
        smopClient.sendRecipe(recipeName,F(:,jj),false,cfm,usv);
        clear recipeName;

        fprintf("[%d] %s", jj, formatVector(gas.names,F,jj));
    
        M(jj) = getMeasurement(smopClient);

        % THIS VERSION calculates RMSEs of differences between vector 
        % pairs as a measure of closeness to the initial dispersion plot. 
        cfm = getSimilarityMeasure(args.alg,initM,M(jj));
        if (cfm < gm)
            gm = cfm;
            gmId = jj;
        end
        fprintf(" RMSE=%6.3f\n", cfm);
    end

    fprintf("GM: %.4f [%s]\n", gm, formatVector(gas.names,F,gmId));
    
    %% STEP 5: Iterative step

    iter = 1;   % iteration counter
    MI = XI;    % indexes of M corresponding to vectors consisting X

    isFinished = false;     % switch for terminating iterative process

    while (~isFinished)
        fprintf("\nIteration #%d:\n", iter);

        im = 1e8;   % minimum RMSE for vectors tested in this iteration
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

        % Compute cost function values for trial vectors and corresponding
        %   target vectors
    
        for jj = XI
            %% STEP 5b1: Compute RMSE of TARGET vector

            cfX = getSimilarityMeasure(args.alg,initM,M(MI(jj)));

            %% STEP 5b2: Compute RMSE of TRIAL vector
    
            % Find shuffled intertersection of indexes for which F and U
            % have same flow rates
            C = getCommonIndices(F,U,jj);
    
            if ~isempty(C)  % measurement of this vector exists
                kk = C(1);  % already, lets take it from the database
                fprintf("[%d] REPEAT  %s", jj, formatVector(gas.names,F,kk));
            else
                pause(0.5);     % emulate some heavy ML search :)
    
                % This vector was not measured yet, so lets do it
                recipeName = sprintf("Iteration #%d, Search #%d", iter, jj); 
                smopClient.sendRecipe(recipeName,U(:,jj),isFinished,cfm,usv);
                clear recipeName;
    
                % Wait for new measurement and add it to the table of 
                % measured recipes
                fprintf("[%d] MEASURE %s", jj, formatVector(gas.names,U,jj));
                kk = length(M) + 1;
                M(kk) = getMeasurement(smopClient);
                F(:,kk) = U(:,jj);
            end

            cfU = getSimilarityMeasure(args.alg,initM,M(kk));
            fprintf(" RMSE=%6.3f", cfU);

            if isempty(C)   % memorize cf of measured data for recipe
                cfm = cfU;
            end

            clear C;
            
            %% STEP 5b3: update RMSE minimas

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
    
            %% STEP 5b4: replace target vec with trial vec if it has lower RMSE

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
            fprintf("  %s, RMSE = %.4f\n\nFinished\n\n", ...
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
end

% Finds common indices and returns a permutated array of them.
% Note: currently supports 2-gas mixtures only.
function idMix = getCommonIndices(V, U, ucol)
    % Find samples for which V and U have same flow rate for 
    %   the first gas (isopropanol) and ..
    idGas1 = find(V(1,:) == U(1,ucol));
    %   .. the second gas (ethanol or n-Butanol)
    idGas2 = find(V(2,:) == U(2,ucol));
    
    % Get intersection and shuffle it in random order
    % NOTE: some vectors might be missing from the data
    idMix = intersect(idGas1,idGas2);
    idMix = idMix(randperm(numel(idMix)));
end

% Computes a similarity (distance) measure.
% For now, we use RMSE to represent such a measure.
function rmse = getSimilarityMeasure(alg, measrm1, measrm2)
    rmse = 1e8;

    % COST FUNCTION WILL BE FUNCTION THAT MEASURES SIMILARITY OF DISPERSION
    % PLOT OF MIXTURE WE WANT TO RECREATE AND ONE TRAINING DISPERSION PLOT
    % NAA WILL TEST DIFFERENT MEASURES IN HER THESIS
    if (alg == "euclidean")
        if ~isempty(measrm2.dms)
            rmse = sqrt(mean((measrm1.dms - measrm2.dms).^2));
        elseif ~isempty(measrm2.snt)
            rmse = sqrt(mean((measrm1.snt - measrm2.snt).^2));
            if rmse > 1000              % huge values observed with real SNT
                rmse = rmse / 100000;   % something to remove in future
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
        % pick three distinct dispersion plots (also PID, SNT) (need to be
        % distinct from each other and from x)
        % QUESTION: Do they need to be from different concentration
        % combinations or could they be also from the same concentration
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
% NOTE: will be applied only for 2-dimensional vectors (2 gases)
function V = limitVectors(V, limits)
    [n,np] = size(V);
    if (n == 2)
        for jj = 1:np
            a = atan2(V(2,jj),V(1,jj));
            lf = [limits(1) * cos(a); limits(2) * sin(a)];
            V(:,jj) = min(lf,V(:,jj));
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

% Computes the separation voltage (usv) closest to the linePosition that
% should be from [0:1] interval (0 corresponds to dms.setup.usv.min 
% and 1 corresponds to dms.setup.usv.max)
function [usv, data] = getDmsSingleLine(dms, linePosition)
    % compute the Us line at a specific location in the full DMS scan:
    step = round(linePosition * dms.setup.usv.steps);
    interval = (dms.setup.usv.max - dms.setup.usv.min) / ...
        (dms.setup.usv.steps - 1);
    usv = dms.setup.usv.min + interval * (step - 1);

    start = dms.setup.ucv.steps * (step - 1);
    end_ = start + dms.setup.ucv.steps;
    data = dms.data.positive(start+1:end_);
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
