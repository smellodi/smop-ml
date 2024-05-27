% Compares pairs of DMS measurements
% Should be used with data collected earlier and saved in JSON files, 
% or for comparing any two DMS measurements save in files.

%% APP ENTRY
function cmp(varargin)

    args = argparser(varargin,nargin);
    if size(args) == 0
        return
    end
    args.folder = args.folder + "/";
    
    if args.action == "2cmp"
        files = [args.folder + args.fn1 ...
                 args.folder + args.fn2]';
    elseif strncmp(args.action, "5", 1)
        range = [2 3 5 7 9 11 13 15 17 19 21 22 23 25 27 29 31 33 35 ...
                 37 39 41 42 43 45 47 49 51 53 55 57 59 61 62]';
        files = dir(args.folder + "*.json");
        try
            files = arrayfun(@(x) args.folder + files(x).name,range);
        catch
            fprintf("Not enough files for this action");
            return
        end
    elseif strncmp(args.action, "6", 1)
        if args.action == "6ipa"
            range = [2 3 4 8 10 12 14 15 16 20 22 24 26 27 28 32 34 36]';
        else
            range = [5 6 7 9 11 13 17 18 19 21 23 25 29 30 31 33 35 37]';
        end
        files = dir(args.folder + "*.json");
        try
            files = arrayfun(@(x) args.folder + files(x).name,range);
        catch
            fprintf("Not enough files for this action");
            return
        end
    end

    if args.action == "5fst"    % comparison against the first DMS
        dms1 = readDms(files(1));
        for jj = 2:size(files)
            dms2 = readDms(files(jj));
            if (size(dms1,1) <= 1) || (size(dms2,1) <= 1)
                continue
            end
            cfm = getSimilarityMeasure(args.alg,dms1,dms2);
            fprintf("%.3f\n", cfm);
        end
    else                        % sequential comparison
        dms1 = readDms(files(1));
        for jj = 2:size(files)
            dms2 = readDms(files(jj));
            if (size(dms1,1) <= 1) || (size(dms2,1) <= 1)
                continue
            end
            cfm = getSimilarityMeasure(args.alg,dms1,dms2);
            fprintf("%.3f\n", cfm);

            dms1 = dms2;
        end
    end
end

%%% FUNCTIONS

% Reads IonVision data file, returns positive DMS data
% or 0, if cannot read file or the file is corrupted, etc.
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
    if (alg == "euclidean")
        % For Euclidean algorithm we use RMSE as the distance.
        distance = sqrt(mean((dms1 - dms2).^2));
    else
        distance = pdist2(dms1',dms2',"euclidean");
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
        "action", '', ...
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
        fprintf("    action [str]  (required) action to perform\n");
        fprintf("                      available values: [%s]\n", join(actions,", "));
        fprintf("    fn1 [str]     DMS filename #1, 'fn1=' can be omitted\n");
        fprintf("                      used only in action '%s'\n", actions(1));
        fprintf("    fn2 [str]     DMS filename #2, 'fn2=' can be omitted\n");
        fprintf("                      used only in action '%s'\n", actions(1));
        fprintf("    alg [str]     algorithm to use (default is '%s')\n", dict.alg);
        fprintf("                      available values: [%s]\n", join(algorithms,", "));
        fprintf("\n");
        fprintf("Example:\n");
        fprintf("    cmp action=2cmp 1.json 2.json\n");
        fprintf("    cmp action=6ipa folder=6_heatup_PID\n");

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