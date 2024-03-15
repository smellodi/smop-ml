%% NOTES

% TODO:
% - ensure each generated flow pair does not oversaturate DMS or PID
% - apply generated flow rounding also after extra randomization

% TIMO SUGGEST TO USE FLOW RATES OF 20 OR HIGHER FOR TARGET SCENT TO SPEED
% UP RESPSONSE TIME/SHORTEN STABILIZATION TIME
% USE LOWER RESOLUTION --> MEASURE ONLY REGIONS WITH ALPHA CURVES FROM
% ANTON'S PAPER: SV 450 to 700, CV -1 to 4
% for future implementation I would still suggest to use Atte's training
% base as starting point to speed up search

% The learning process uses DIFFERENTIAL EVOLUTION algorithm described in 
% Mariaana's manuscript.
%
% IF WE NOW WANT TO ALSO INCLUDE SMELL INSPECTOR AND PID SENSOR THE WHOLE
% SYSTEM WILL NEED TO BE REDESIGNED
%
% Philipp Muller, January 2024

% THIS VERSION uses full network message exchange protocol as explained in 
% /TG-SMELLODI/General/Software/algorithm_Smellodi.docx
%
% IMPORTANT: Recipe name is parsed in the SMOP app, so do not change the
% way it is constructed.
%
% Oleg Spakov, March 2024

%%% App entry

function smop_ml(varargin)
    args = argparser(varargin, nargin);

    % dbstop error
    % rng("shuffle")

    %% STEP 1: Establish a connection to SMOP and get config parameters.
    
    % create a TCP client and connect to SMOP software
    fprintf('Connecting to %s ...\n', args.ip);
    smopClient = SmopClient(args.ip);
    
    % continue only if the connection is OK
    if ~smopClient.isConnected
        pause(3);
        return;
    end
    
    fprintf('Connected on port %d\n', smopClient.port);
    
    % continue only if config packet was received
    if ~smopClient.readConfig()
        fprintf('Exiting...\n');
        pause(2);
        return;
    end
    
    gases = smopClient.gases;
    
    %% STEP 2: Define parameters etc.
    
    % - switch for terminating iterative process
    isFinished = false;
    
    % - number of chemicals in the mixture (concentration can be also 0 sccm)
    n = smopClient.getChannelCount(2);
    
    % - minimum and maximum flows of chemicals in the mixture (in sccm)
    minFlow = 0;
    maxFlow = smopClient.getMaxFlow(50);
    
    % - maximum number of iterations
    mi = smopClient.getMaxIterationCount(args.mi);
    
    % - threshold RMSE for terminating iterative process
    limRMSE = smopClient.getThreshold(args.th);
    
    % Define which distance measure to use. Options are:
    % - 'euclidean': uses Euclidean distance between two matrices
    % - 'alpha': extracts alpha curves or dispersion plots using Anton's
    %   algorithm (ISOEN paper) and computes Euclidean distance between the
    %   alpha curve locations
    % - 'frobenius': Frobenius norm
    % - 'cos': Cosine similarity
    % - 'stat': stat tests on diff between generated and original odors?
    % - what else could be tested by Naa???
    % THIS VERSION calculates RMSEs of differences between vector 
    % pairs as a measure of closeness to the initial dispersion plot. 
    % FOR NOW USE ONLY DMS DATA, no PID or SNT data used.

    %alg = smopClient.getAlgorithm(args.alg);
    
    % crossover probability CR: [0,1] (Wikipedia lists 0.9 as standard)
    % - cr maintains diversity of population
    % - for non-separable problems with correlated parameters cr in [0.5, 1]
    %   is recommended
    % - Mariaana used cr = 0.5 for n <= 3 and cr = 0.9 for n >= 4
    % - she pointed out that cr = 0.9 could have been used for any n (also
    cr = args.cr;
    
    % differential weight f: [0,2] (Wikipedia lists 0.8 as standard setting)
    % - f controls differential variation and has most severe impact on
    %   performance of differential evolution algorithm
    f = args.f;
    
    % - DMS scan parameter
    useSingleUsv = args.ssv;

    % Randomization offset. This is max offset applied to the originally 
    % generated flow value if this flow does not pass the validation process
    % (used in "validate" function)
    randDelta = 0.2*(maxFlow-minFlow);

    % Limits used to detect if the generated flow values may cause
    % measurement sensor oversaturation. To be used in "limitVectors" function.
    limits = smopClient.getCriticalFlow([55; 60]);

    fprintf('Gas params: n = %d, minFlow = %d, maxFlow = %d\n', ...
        n, minFlow, maxFlow);
    fprintf('Search params: mi = %d, limRMSE = %.2f cr = %.2f, f = %.2f\n', ...
        mi, limRMSE, cr, f);
    
    %% STEP 3: Read in message with DMS measurement from target scent
    
    dms = smopClient.getInitialDms();

    if useSingleUsv
        [usv, posDP] = getDmsSingleLine(dms, 0.5);
        fprintf('DMS scan parameters: Us = %.1f V\n', usv);
    else   % use full DMS scan
        usv = 0;
        posDP = dms.data.positive;

        fprintf('DMS scan parameters: full scan\n');
    end

    clear dms;

    %% STEP 4: Initialization
    
    % Start without any prior information and measure np mixtures.
    
    %% STEP 4a: Init vectors and consts
    % All measured vectors (flow rates)
    F = [minFlow maxFlow minFlow maxFlow (maxFlow + minFlow)/2;
         minFlow minFlow maxFlow maxFlow (maxFlow + minFlow)/2];

    % Both flows set to maxFlow result in oversaturated PID, therefore 
    %   we apply some limitations as described in limitVectors.
    F = limitVectors(F,limits);
    F = roundTo(F,args.dec);

    % Population size: set of np best vectors
    % Typically 10 times the dimension n of x, but 
    %   this may be too long... keep it 5 and not 20 for the demo.
    X = F;
    
    np = size(X,2);

    % Indexes of M corresponding to vectors consisting X
    p = nan(1,np);

    % following loop might be not necessary, but keep it here to allow 
    % various combinations of flows in X that may come in future.
    for jj = 1:np
        % find shuffled intersection of indexes for which F and X
        % have same flow rates
        idMix = getCommonIndices(F, X, jj);
        
        % use DMS measurement related to first index for initial population
        if ~isempty(idMix)
            p(jj) = idMix(1);
        end
        clear idMix
    end
    
    %% STEP 4b: Collect np measurements

    % All measurements
    M = arrayfun(@(x)struct('pos',1),1:np);

    % overall minimum RMSE
    GM = 1e8;
    cf = 1e5;   % minima
    
    fprintf('\nCollecting initial measurements:\n');
    for jj = 1:np
        pause(0.5); % simply, to make ML being not too fast in SMOP interface
    
        recipeName = sprintf('Reference #%d', jj);
        smopClient.sendRecipe(recipeName, F(:,jj), false, cf, usv);
        fprintf('[%d] %s', jj, formatVector(gases, F, jj));
    
        dms = smopClient.waitForDms();
        M(jj).pos = dms.data.positive;

        cf = sqrt(mean((posDP - dms.data.positive).^2));
        if (cf < GM)
            GM = cf;
            idGM = jj;
        end
        fprintf(' RMSE=%6.3f\n', cf);

        clear dms;
    end

    fprintf('GM: %.4f [%s]\n', GM, formatVector(gases, F, idGM));
    
    %% STEP 5: Iterative step

    iter = 1;       % iteration counter
    cfm = cf;       % cf of measured trials only
    
    while (~isFinished)
        fprintf('\nIteration #%d:\n', iter);

        % LM = 1e8;   % minimum RMSE for vectors in X
        UM = 1e8;   % minimum RMSE for vectors tested in this iteration
    
        %% STEP 5a: Differential evolution
        % [https://en.wikipedia.org/wiki/Differential_evolution]

        V = mutate(X, f);                       % generate new vectors
        V = limitValues(V, minFlow, maxFlow);   % limit values of new vectors
        U = crossover(X, V, cr);                % mix old and new vectors
        U = validate(U, randDelta);             % remove repetitions
        U = limitValues(U, minFlow, maxFlow);   % limit values after rndmzation
        U = limitVectors(U, limits);            % avoid oversaturation
        U = roundTo(U, args.dec);               % round flow values

        fprintf('Flows to test:\n%s', formatVectorsAll(gases, U));
        
        %% STEP 5b: Selectiing better vectors

        % - compute cost function values for trial vectors and corresponding
        %   target vectors
        % COST FUNCTION WILL BE FUNCTION THAT MEASURES SIMILARITY OF DISPERSION
        % PLOT OF MIXTURE WE WANT TO RECREATE AND ONE TRAINING DISPERSION PLOT
        % NAA WILL TEST DIFFERENT MEASURES IN HER THESIS
        % 
        % for NOW USE RMSE AS DISTANCE MEASURE
    
        for jj = 1:np
    
            if isnan(p(jj)) % OLEG: seems that this is never true
                continue;
            end
    
            %% STEP 5b1: RMSE of TARGET vector

            tv = M(p(jj)).pos;
           
            cf_X = sqrt(mean((posDP - tv).^2));
            clear tv

            %% STEP 5b2: RMSE of TRIAL vector
    
            % find shuffled intertersection of indexes for which F and U
            % have same flow rates
            idMix = getCommonIndices(F, U, jj);
    
            if ~isempty(idMix)              % measurement of this vector
                idPair = idMix(1);          % exists already, lets 
                tv = M(idPair).pos;      % take it from the database
                fprintf('[%d] REPEAT  %s', jj, formatVector(gases, F, idPair));
            else
                pause(0.5);     % emulate some heavy ML search :)
    
                % this vector was not measured yet, so lets do it
                recipeName = sprintf('Iteration #%d, Search #%d', iter, jj); 
                smopClient.sendRecipe(recipeName, U(:,jj), isFinished, cfm, usv);
                fprintf('[%d] MEASURE %s', jj, formatVector(gases, U, jj));
    
                % wait for new DMS measurement and add it to the table of 
                % measured recipes
                dms = smopClient.waitForDms();
                M(end + 1).pos = dms.data.positive;
    
                F = [F U(:,jj)];
                idPair = size(F,2);
    
                tv = dms.data.positive;
                clear dms;
            end
    
            cf_U = sqrt(mean((posDP - tv).^2));
            fprintf(' RMSE=%6.3f', cf_U);

            if isempty(idMix)   % memorize cf of measured data for recipe
                cfm = cf_U;
            end

            clear idMix tv;
            
            %% STEP 5b3: update RMSE minima

            cf = min(cf_U,cf_X);

            % storing global minima
            if cf < GM
                GM = cf;
                idGM = idPair;
                fprintf(' GM');
            end

            % storing local minima
            % if cf < LM
            %     LM = cf;
            %     idLM = jj;
            % end

            % minima of the tested vectors
            if cf_U < UM
                UM = cf_U;
                idUM = jj;
                fprintf(' UM');
            end
    
            %% STEP 5b4: replace target vec with trial vec if it has lower RMSE

            if cf_U < cf_X
                X(:,jj) = U(:,jj);
                fprintf(' p(%d)= %d >> %d', jj, p(jj), idPair);
                p(jj) = idPair;
            end
    
            fprintf('\n');
            clear cf_U cf_X;
        end
        
        fprintf('UM: %.4f [%s]\n', UM, formatVector(gases, U, idUM));
        fprintf('GM: %.4f [%s]\n', GM, formatVector(gases, F, idGM));
        
        %% STEP 5c: Make a decision about the proximity of the best guess
    
        if (GM < limRMSE)
            isFinished = true;
            recipeName = 'Final recipe';
            fprintf('\n%s:\n', recipeName);
        elseif (iter >= mi)
            isFinished = true;
            recipeName = sprintf('Best after %d iterations', mi);
            fprintf('\n%s:\n', recipeName);
        end
    
        if isFinished
            % send the final recipe
            flows = F(:,idGM);
            smopClient.sendRecipe(recipeName, flows, isFinished, cfm);
            fprintf('  %s, RMSE = %.4f\n\nFinished\n\n', ...
                formatVector(gases, F, idGM), GM);
        else
            fprintf('Continuing the search, best vectors are:\n%s', ...
                formatVectorsAll(gases, X));
        end
        
        iter = iter + 1;  % increase counter iterations by one
    end
    
    %% STEP 6: Finalize
    
    clear smopClient;   % disconnect from the SMOP
    pause(3);           % short pause before closing the CMD window

end

%%% FUNCTIONS

% Finds common indices and returns a permutated array of them.
% 'index' is the column to use from U matrix
function idMix = getCommonIndices(V, U, index)
    % find samples for which V and U have same flow rate for ..
    % .. isopropanol and ..
    idGas1 = find(V(1,:) == U(1,index));
    % .. ethanol or n-Butanol
    idGas2 = find(V(2,:) == U(2,index));
    
    % get intersection and shuffle it in random order
    % NOTE: some vectors might be missing from the data
    idMix = intersect(idGas1,idGas2);
    idMix = idMix(randperm(numel(idMix)));
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

    % generate randomly chosen indices to ensure that at least one
    % component of the donor vector is included in the target vector
    rv = randi(n,1,np);

    % random numbers in [0,1] for each component of each target vector
    rv2 = rand(n,np);

    % combining target and donor vectors to get trial vectors
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

    % If all flow values of a certain gas are the same...
    allSame = arrayfun(@(x) true, V(:,1));
    for jj = 2:np
        for kk = 1:n
            allSame(kk) = allSame(kk) && V(kk,jj-1) == V(kk,jj);
        end
    end

    % ...then randomize those values a bit
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
% We assume that limirs are semi-axis (a,b) of an ellipse, and if the
% vector lies outside of this ellipse, then it is "pulled" toward the
% center to the ellipse edge.
% NOTE: applied only for 2-dimensional vectors (2 gases)
function V = limitVectors(V, limits)
    [n,np] = size(V);
    if (n == 2)
        for jj = 1:np
            a = atan2(V(2,jj), V(1,jj));
            lf = [limits(1)*cos(a); limits(2)*sin(a)];
            V(:,jj) = min(lf,V(:,jj));
        end
    end
end

% Rounds values in the vector.
% If dec < 0, then round value to be divisible by |dec|.
function V = roundTo(V, dec)
    r = max(dec, 0);
    V = round(V, r);

    if dec < 0
        [n,np] = size(V);
        dec = abs(dec) + 1;
        for kk = 1:n
            for jj = 1:np
                r = mod(V(kk,jj),dec);
                V(kk,jj) = V(kk,jj) - r;
            end
        end
    end
end

% Computes the separation voltage (usv) closest to the linePosition that
% should be from [0:1] interval (0 corresponds to dms.setup.usv.min 
% and 1 corresponds to dms.setup.usv.max)
function [usv, posDP] = getDmsSingleLine(dms, linePosition)
    % compute the Us line at a specific location in the full DMS scan:
    step = round(linePosition * dms.setup.usv.steps);
    interval = (dms.setup.usv.max - dms.setup.usv.min) / ...
        (dms.setup.usv.steps - 1);
    usv = dms.setup.usv.min + interval * (step - 1);

    start = dms.setup.ucv.steps * (step - 1);
    end_ = start + dms.setup.ucv.steps;
    posDP = dms.data.positive(start+1:end_);
end

% Formats gas names with their flows as "GAS1=FLOW1 GAS2=FLOW2 ..."
function s = formatVector(names, V, index)
    s = '';
    c = size(V,1);
    for ii = 1:c
        if ii ~= 1
            s = sprintf('%s ', s);
        end
        s = sprintf('%s%s=%4.1f', s, names(ii), V(ii, index));
    end
end

% Formats gas names with a list of flows as "  GAS1 FLOW1 FLOW2 ..\n GAS2 .."
function s = formatVectorsAll(names, V)
    s = '';
    [n,np] = size(V);
    for jj = 1:n
        s = sprintf('%s  %-12s', s, names(jj));
        for kk = 1:np
            s = sprintf('%s %4.1f', s, V(jj,kk));
        end
        s = sprintf('%s\n', s);
    end
end