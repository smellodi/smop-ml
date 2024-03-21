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
        "cr", 0.7, ...
        "f", 0.8, ...
        "dec", 0, ...
        "mi", 8, ...
        "th", 0.5, ...
        "alg", algorithms(1), ...
        "ssv", false ...
    );

    fprintf("===============================\n");
    fprintf("SMOP Machine Learning back-end.\n");
    fprintf("===============================\n");
    fprintf("\n");
    fprintf("Usage:\n");
    fprintf("  smop_ml [<arg>=<value>]\n\n");
    fprintf("<arg> can be one of the following:\n");
    fprintf("    ip [str]      IP of SMOP UI app (default is %s)\n", dict.ip);
    fprintf("    cr [double]   crossover population (default is %.2f)\n", dict.cr);
    fprintf("    f  [double]   differential weight(default is %.2f)\n", dict.f);
    fprintf("    dec [int]     decimals left after rounding generated flows" + ...
        " (default is %d)\n", dict.dec);
    fprintf("                      note that negative values round flows " + ...
        "to be divisible by |dec|+1\n");
    fprintf("                      (thus, dec=-1 makes flow values to be even\n");
    fprintf("    mi [int]      max iteration count (default is %d)\n", dict.mi);
    fprintf("    th [double]   proximity threshold (default is %.2f)\n", dict.th);
    fprintf("    alg [str]     algorithm to use (default is '%s')\n", dict.alg);
    fprintf("                      available values: [%s]\n", join(algorithms,", "));
    fprintf("    ssv [bool]    flag that enabled single DMS separation voltage\n");
    fprintf("                      (default is %s, i.e. full scan)\n", ...
        string(dict.ssv));
    fprintf("\n");
    fprintf("Example:\n");
    fprintf("    smop_ml ip=192.168.1.4 f=0.6 ssv\n");   
    fprintf("    smop_ml cr=0.9 dec=-1 mi=12 ssv=true\n");   
    fprintf("\n");
    fprintf("Note how flags can be declared with or without the boolean value:" + ...
        "\n    use either 'ssv', or 'ssv=true' or 'ssv=1' to enable 'ssv'\n");
    fprintf("\n\n");

    for jj = 1:n
        p = split(string(args(jj)),"=");
        wasParsed = true;

        if (length(p) == 1)    % parse flags as boolean values set to "true"
            if strcmpi(p(1), "ssv")
                dict.ssv = true;
            else
                wasParsed = false;
            end
        else        % parse all arguments as a pair of "key"="value"
            try
                if strcmpi(p(1), "ip")
                    dict.ip = p(2);
    % Crossover probability CR: [0,1] (Wikipedia lists 0.9 as standard)
    % - cr maintains diversity of population
    % - for non-separable problems with correlated parameters cr in [0.5, 1]
    %   is recommended
    % - Mariaana used cr = 0.5 for n <= 3 and cr = 0.9 for n >= 4
    % - she pointed out that cr = 0.9 could have been used for any n (also
                elseif strcmpi(p(1), "cr")
                    dict.cr = str2double(p(2));
    % Differential weight f: [0,2] (Wikipedia lists 0.8 as standard setting)
    % - f controls differential variation and has most severe impact on
    %   performance of differential evolution algorithm
                elseif strcmpi(p(1), "f")
                    dict.f = str2double(p(2));
    % Decimals left after rounding generated flows. Note that negative 
    % values round flows to be divisible by |dec|+1
                elseif strcmpi(p(1), "dec")
                    dict.dec = str2double(p(2));
    % Maximum number of iterations
                elseif strcmpi(p(1), "mi")
                    dict.mi = str2double(p(2));
    % RMSE threshold for terminating iterative process
                elseif strcmpi(p(1), "th")
                    dict.th = str2double(p(2));
    % Define which distance measure to use.
                elseif strcmpi(p(1), "alg")
                    dict.alg = p(2);
    % Enabled single-separation-voltage mode for DMS measurements (scope mode).
                elseif strcmpi(p(1), "ssv")
                    dict.ssv = strcmpi(p(2), "true") || strcmp(p(2), "1");
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