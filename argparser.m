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

    fprintf('===============================\n');
    fprintf('SMOP Machine Learning back-end.\n');
    fprintf('===============================\n');
    fprintf('\n');
    fprintf('Usage:\n');
    fprintf('  smop_ml { <arg>=<value>}\n\n');
    fprintf('<arg> can be one of the following:\n');
    fprintf('    ip [str]      IP of SMOP UI app (default is %s)\n', dict.ip);
    fprintf('    cr [double]   crossover population (default is %.2f)\n', dict.cr);
    fprintf('    f  [double]   differential weight(default is %.2f)\n', dict.f);
    fprintf(['    dec [int]     decimals left after rounding generated flows' ...
        ' (default is %d)\n'], dict.dec);
    fprintf(['                      note that negative values round flows ' ...
        'to be divisible by |dec|+1\n']);
    fprintf('                      (thus, dec=-1 makes flow values to be even\n');
    fprintf('    mi [int]      max iteration count (default is %d)\n', dict.mi);
    fprintf('    th [double]   proximity threshold (default is %.2f)\n', dict.th);
    fprintf('    alg [str]     algorithm to use (default is "%s")\n', dict.alg);
    fprintf('                      available values: [%s]\n', join(algorithms,", "));
    fprintf('    ssv [bool]    flag that enabled single DMS separation voltage\n');
    fprintf('                      (default is %s, i.e. full scan)\n', ...
        string(dict.ssv));
    fprintf('\n');
    fprintf('Example:\n');
    fprintf('    smop_ml ip=192.168.1.4 f=0.6 ssv\n');   
    fprintf('    smop_ml cr=0.9 dec=-1 mi=12 ssv=true\n');   
    fprintf('\n');
    fprintf(['Note how flags can be declared with or without the boolean value:' ...
        '\n    use either "ssv", or "ssv=true" or "ssv=1" to enable "ssv"\n']);
    fprintf('\n\n');

    for jj = 1:n
        p = split(string(args(jj)),"=");

        if (length(p) == 1)    % parse flags as boolean values set to 'true'
            if strcmpi(p(1), "ssv")
                dict.ssv = true;
            end
        else        % parse all arguments as a pair of 'key'='value'
            try
                if strcmpi(p(1), "ip")
                    dict.ip = p(2);
                elseif strcmpi(p(1), "cr")
    % Crossover probability CR: [0,1] (Wikipedia lists 0.9 as standard)
    % - cr maintains diversity of population
    % - for non-separable problems with correlated parameters cr in [0.5, 1]
    %   is recommended
    % - Mariaana used cr = 0.5 for n <= 3 and cr = 0.9 for n >= 4
    % - she pointed out that cr = 0.9 could have been used for any n (also
                    dict.cr = str2double(p(2));
                elseif strcmpi(p(1), "f")
    % Differential weight f: [0,2] (Wikipedia lists 0.8 as standard setting)
    % - f controls differential variation and has most severe impact on
    %   performance of differential evolution algorithm
                    dict.f = str2double(p(2));
                elseif strcmpi(p(1), "dec")
    % Decimals left after rounding generated flows. Note that negative 
    % values round flows to be divisible by |dec|+1
                    dict.dec = str2double(p(2));
                elseif strcmpi(p(1), "mi")
    % Maximum number of iterations
                    dict.mi = str2double(p(2));
                elseif strcmpi(p(1), "th")
    % RMSE threshold for terminating iterative process
                    dict.th = str2double(p(2));
                elseif strcmpi(p(1), "alg")
    % Define which distance measure to use.
                    dict.alg = p(2);
                elseif strcmpi(p(1), "ssv")
                    dict.ssv = strcmpi(p(2), 'true') || strcmp(p(2), "1");
                end
            catch ex
                fprintf("Failed to parse command-line argument '%s': %s", ...
                    args(jj), ex.message);
            end
        end
    end
end