% Suppemental script:
% 1 - to study repeatability of measurements
%       use arg param to set the target vector
% 2 - test PID level for a certain flow rate
% 3 - measure RMSE for various offsets (rough)
%       use arg param to set the test channel/odor
% 4 - measure RMSE for various offsets (fine)
% 5 - heating up the system by repeating weak flow multiple times
%       use arg param to set the target flow (same for both channels)

%% APP ENTRY
function exprm(varargin)
    args = argparser(varargin,nargin);
    if size(args) == 0
        return
    end

    %% STEP 1: Init consts

    if args.mode == 1
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
    
        flowIndex = args.arg;
    
        flow1 = flows(flowIndex,1);
        flow2 = flows(flowIndex,2);
    
        min_ = 6;
        max_ = 46;
    
        F = [flow1 min_ flow1 max_ flow1 min_ flow1 45 flow1 25 flow1 10 flow1 25 flow1 25 flow1 40 flow1 flow1; ...
             flow2 min_ flow2 min_ flow2 max_ flow2 45 flow2 25 flow2 25 flow2 10 flow2 40 flow2 25 flow2 flow2];
        F = [F F F];
    elseif args.mode == 2
        flow1 = 40;
        flow2 = 0;
        F = [flow1 flow1 flow1 flow2 flow2 flow2 flow1 flow2 flow1 flow2 flow1 flow2; ...
             flow2 flow2 flow2 flow1 flow1 flow1 flow2 flow1 flow2 flow1 flow2 flow1];
        F = [F F F];
    elseif args.mode == 3
        channel = args.arg;
        F(channel,:) = [5 10 15 20 25 30 35 40 45 50 50 45 40 35 30 25 20 15 10 5];
        F(3 - channel,:) = 0;
        flow1 = 0;
        flow2 = 0;
    elseif args.mode == 4
        flow1 = 5; % 5 15 25 35 45
        flow2 = 0; % 0 10 20 30 40 50
        channel = 2;
        F(1,:) = 0;
        F(2,:) = 0;
        id = 1;
        for jj = -18:3:18
            flow = flow1 + jj;
            if (flow > 0) && (flow <= 50)
                F(channel,id) = flow;
                F(3 - channel,id) = flow2;
                id = id + 1;
            end
        end
        if channel == 2
            t = flow2;
            flow2 = flow1;
            flow1 = t;
        end
    elseif args.mode == 5
        if (args.arg < 5)
            f = 5;
        else
            f = args.arg;
        end
        A = [f f f f f f f f f f f f f f f f f f f f];
        F = [A; A];
        F = [F F F F F];  % 100 times
        flow1 = 0;
        flow2 = 0;
    else
        fprintf("Unknown mode '%d'. Exiting...\n", args.mode);
        pause(3);
        return;
    end

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
    fprintf("0 %s\n", formatVector(gas.names,[flow1;flow2],1));

    for jj = 1:size(F,2)
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

function dict = argparser(args, n)
    algorithms = [
        % implemented
        "euclidean"     % uses Euclidean distance between two matrices

        % the rest are not implemented
        "alpha"         % extracts alpha curves or dispersion plots using 
                        %   Anton's algorithm (ISOEN paper) and computes 
                        %   Euclidean distance between thealpha curve locations
        "frobenius"     % Frobenius norm
        "cos"           % Cosine similarity
        "stat"          % stat tests on diff between generated and original odors?

        % - what else could be tested by Naa???
    ];

    dict = struct( ...
        "ip", "127.0.0.1", ...
        "mode", 0, ...
        "arg", 1, ...
        "alg", algorithms(1) ...
    );

    if n == 0
        fprintf("===============================\n");
        fprintf("SMOP Machine Learning back-end.\n");
        fprintf("Flow test routine\n");
        fprintf("===============================\n");
        fprintf("\n");
        fprintf("Usage:\n");
        fprintf("  exprm [<arg>=<value>]\n\n");
        fprintf("<arg> can be one of the following:\n");
        fprintf("    ip [str]      IP of SMOP UI app (default is %s)\n", dict.ip);
        fprintf("    mode [int]    (required) one of the following numbers:\n");
        fprintf("                  1: distractors + test flow\n");
        fprintf("                  2: PID level test\n");
        fprintf("                  3: RMSE dependency on a flow offset\n");
        fprintf("                  4: Weak flow multiple times (heating up)\n");
        fprintf("                  5: 100 pulses of a mixture of the same flow rates \n");
        fprintf("    arg [int]     an argument depending on the mode:\n");
        fprintf("                  mode=1: flow pair index, 1-10\n");
        fprintf("                  mode=2,4: ignored\n");
        fprintf("                  mode=3: test channel, 1/2\n");
        fprintf("                  mode=5: flow rate, >=5sccm (same for both channels)\n");
        fprintf("    alg [str]     algorithm to use (default is '%s')\n", dict.alg);
        fprintf("                      available values: [%s]\n", join(algorithms,", "));
        fprintf("\n");
        fprintf("Example:\n");
        fprintf("    exprm mode=1 arg=2\n");

        dict = {};
        return
    end

    for jj = 1:n
        p = split(string(args(jj)),"=");
        wasParsed = true;

        if (length(p) == 1)    % parse flags as boolean values set to "true"
            if strcmpi(p(1), "-")
                % no such parameters yet, the comparison above is a dummy
            else
                wasParsed = false;
            end
        else        % parse all arguments as a pair of "key"="value"
            try
                if strcmpi(p(1), "ip")
                    dict.ip = p(2);
                elseif strcmpi(p(1), "arg")
                    dict.arg = str2double(p(2));
                elseif strcmpi(p(1), "mode")
                    dict.mode = str2double(p(2));
                elseif strcmpi(p(1), "alg")
                    dict.alg = p(2);
                else
                    wasParsed = false;
                end
            catch ex
                fprintf("Failed to parse command-line argument '%s': %s\n", ...
                    args(jj), ex.message);
            end
        end

        if (~wasParsed)
            fprintf("Unknown parameter '%s'\n", p(1));
        end
    end
end