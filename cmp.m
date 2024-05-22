% Compares pairs of DMS measurements
% Should be used with data collected earlier and saved in JSON files, or
% for comparing any two DMS measurements save in files.

%% APP ENTRY
function cmp(varargin)

    args = argparser(varargin,nargin);
    if size(args) == 0
        return
    end
    
    if args.action == "2cmp"
        files = [args.folder + "/" + args.fn1 ...
                 args.folder + "/" + args.fn2]';
    elseif strncmp(args.action, "5", 1)
        files = dir(args.folder + "/*.json");
        files = [ files(2) ...
            files(3) files(5) files(7) files(9) files(11) ...
            files(13) files(15) files(17) files(19) files(21) files(22) ...
            files(23) files(25) files(27) files(29) files(31) ...
            files(33) files(35) files(37) files(39) files(41) files(42) ...
            files(43) files(45) files(47) files(49) files(51) ...
            files(53) files(55) files(57) files(59) files(61) files(62) ...
        ]';
        files = arrayfun(@(x) args.folder + "/" + x.name,files);
    elseif strncmp(args.action, "6", 1)
        files = dir(args.folder + "/*.json");
        if args.action == "6ipa"
            files = [
                files(2) files(3) files(4) files(8) files(10) files(12) ...
                files(14) files(15) files(16) files(20) files(22) files(24) ...
                files(26) files(27) files(28) files(32) files(34) files(36) ...
            ]';
        else
            files = [
                files(5) files(6) files(7) files(9) files(11) files(13) ...
                files(17) files(18) files(19) files(21) files(23) files(25) ...
                files(29) files(30) files(31) files(33) files(35) files(37) ...
            ]';
        end
        files = arrayfun(@(x) args.folder + "/" + x.name,files);
    end

    if args.action == "5fst"    % comparison against the first DMS
        dms1 = readDms(files(1));
        for jj = 2:size(files)
            dms2 = readDms(files(jj));
            if ~dms1 || ~dms2
                continue
            end
            cfm = getSimilarityMeasure(args.alg,dms1,dms2);
            fprintf("%.3f\n", cfm);
        end
    else                        % sequential comparison
        for jj = 2:size(files)
            dms1 = readDms(files(jj-1));
            dms2 = readDms(files(jj));
            if ~dms1 || ~dms2
                continue
            end
            cfm = getSimilarityMeasure(args.alg,dms1,dms2);
            fprintf("%.3f\n", cfm);
        end
    end
end

%%% FUNCTIONS

% Read IonVision data file, returns positive DMS data
function dms = readDms(filename)
    try
        text = fileread(filename);
        data = jsondecode(text);
        dms = data.MeasurementData.IntensityTop;
    catch ex
        dms = 0;
        fprintf("Failed to read '%s': %s\n", filename, ex.message);
    end
end

% A cost function that measures similarity (distance) between dispersion plot 
% of the mixture we want to recreate and the training dispesion plot. 
function distance = getSimilarityMeasure(alg, dms1, dms2)
    distance = 1e8;

    if (alg == "euclidean")
        % For Euclidean algorithm we use RMSE as the distance.
        distance = sqrt(mean((dms1 - dms2).^2));
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

    actions = [
        % folder contains any DMS data:
        "2cmp"   % compare two DMS

        % folder contains data from MOD_5_heatup test:
        "5seq"   % compute RMSE for pairs of target mixtures.
        "5fst"   % compute RMSE for target mixture DMS compared against the first mixture.

        % folder contains data from MOD_6 test:
        "6ipa"   % compute RMSE for pairs of IPA flows.
        "6eth"   % compute RMSE for pairs of ethanol flows.
    ];

    dict = struct( ...
        "fn1", '', ...
        "fn2", '', ...
        "folder", "data", ...
        "action", actions(1), ...
        "alg", algorithms(1), ...
        "isInfo", false ...
    );

    if n == 0
        fprintf("===============================\n");
        fprintf("SMOP DMS comparison tool.\n");
        fprintf("===============================\n");
        fprintf("\n");
        fprintf("Usage:\n");
        fprintf("  cmp [<arg>=<value>]\n\n");
        fprintf("<arg> can be one of the following:\n");
        fprintf("    folder [str]  DMS folder (default is '%s')\n", dict.folder);
        fprintf("    action [str]  action to perform (default is '%s')\n", dict.action);
        fprintf("                      available values: [%s]\n", join(actions,", "));
        fprintf("    fn1 [str]     DMS filename #1, 'fn1=' can be omitted\n");
        fprintf("                      used only in action '%s'\n", actions(1));
        fprintf("    fn2 [str]     DMS filename #2, 'fn2=' can be omitted\n");
        fprintf("                      used only in action '%s'\n", actions(1));
        fprintf("    alg [str]     algorithm to use (default is '%s')\n", dict.alg);
        fprintf("                      available values: [%s]\n", join(algorithms,", "));
        fprintf("\n");
        fprintf("Example:\n");
        fprintf("    cmp 1.json 2.json\n");
        fprintf("    cmp folder=6_heatup_PID action=6ipa\n");

        dict = {};
        return
    end

    for jj = 1:n
        p = split(string(args(jj)),"=");
        wasParsed = true;

        if (length(p) == 1)    % parse flags as boolean values set to "true"
            if isempty(dict.fn1)
                dict.fn1 = p(1);
            elseif isempty(dict.fn2)
                dict.fn2 = p(1);
            else
                wasParsed = false;
            end
        else        % parse all arguments as a pair of "key"="value"
            try
                if strcmpi(p(1), "fn1")
                    dict.fn1 = p(2);
                elseif strcmpi(p(1), "fn2")
                    dict.fn2 = p(2);
                elseif strcmpi(p(1), "folder")
                    dict.folder = p(2);
                elseif strcmpi(p(1), "action")
                    dict.action = p(2);
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