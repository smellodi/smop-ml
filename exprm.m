% The ML process using DIFFERENTIAL EVOLUTION algorithm
% Experimental script.

%% APP ENTRY
function exprm(varargin)
    args = argparser(varargin,nargin);

    %% STEP 1: Init consts

    flows = [22  3; ...
            3 26; ...
            26	47; ...
            47	22; ...
            11	12; ...
            37	38; ...
            10	42; ...
            39	8; ...
            18	29; ...
            30	20];

    flowIndex = 3;

    min_ = 6;
    max_ = 46;

    flow1 = flows(flowIndex,1);
    flow2 = flows(flowIndex,2);

    F = [flow1 min_ flow1 max_ flow1 min_ flow1 45 flow1 25 flow1 10 flow1 25 flow1 25 flow1 40 flow1 flow1; ...
         flow2 min_ flow2 min_ flow2 max_ flow2 45 flow2 25 flow2 25 flow2 10 flow2 40 flow2 25 flow2 flow2];
    F = [F F F];

    XI = 1:size(F,2); % range of X and U

    %% STEP 2: Establish a connection to SMOP and get config parameters.
    
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
    
    %% STEP 3: Define parameters
    
    gas = struct(); % stores gas info and constant parameters

    gas.min = 0;
    gas.max = smopClient.getMaxFlow(50);
    gas.n = smopClient.getChannelCount(1);
    gas.names = smopClient.gases;

    args.alg = smopClient.getAlgorithm(args.alg);

    %% STEP 4: Read in message with DMS measurement from target scent
    
    initM = smopClient.getInitialMeasurement();

    if isstruct(initM.dms)
        initM.dms = initM.dms.data.positive;
    else
        initM.dms = [];
    end

    if isstruct(initM.snt)
        initM.snt = initM.snt.data.resistances;
    else
        initM.snt = [];
    end

    %% STEP 5: Collect measurements

    cfm = 1e8;  % last search distance (measured trials only)
    
    fprintf("\nTest measurements\n");
    fprintf("0 %s\n", formatVector(gas.names,flows(flowIndex,:)',1));

    for jj = XI
        pause(0.1); % simply, to make ML being not too fast in SMOP interface
    
        recipeName = sprintf("Reference #%d", jj);
        smopClient.sendRecipe(recipeName,F(:,jj),false,cfm);
        clear recipeName;

        fprintf("%d %s", jj, formatVector(gas.names,F,jj));
    
        M = getMeasurement(smopClient);

        cfm = getSimilarityMeasure(args.alg,initM,M);
        fprintf(" DIST %.3f\n", cfm);
    end

    % Send the final recipe
    pause(0.1);
    smopClient.sendRecipe("Final recipe", F(:,1), true, cfm);
    fprintf("Finished\n\n");
    
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

% Formats gas names N with their flows V as "GAS1=FLOW1 GAS2=FLOW2 ..."
function s = formatVector(N, V, vcol)
    s = "";
    for ii = 1:size(V,1)
        if ii ~= 1
            s = sprintf("%s ", s);
        end
        s = sprintf("%s%s %4.1f", s, N(ii), V(ii,vcol));
    end
end
