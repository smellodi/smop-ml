function dict = argparser(args, n)
    algorithms = ["euclidean"];

    dict = struct( ...
        "ip", "127.0.0.1", ...
        "cr", 0.7, ...
        "f", 0.8, ...
        "dec", 0, ...
        "mi", 8, ...
        "th", 0.5, ...
        "alg", algorithms(1), ...
        "usv", false ...
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
    fprintf('    dec [int]     decimals left after rounding generated flows (default is %d)\n', dict.dec);
    fprintf('                      note that -1 round to even value\n');
    fprintf('    mi [int]      max iteration count (default is %d)\n', dict.mi);
    fprintf('    th [double]   proximity threshold (default is %d)\n', dict.th);
    fprintf('    alg [str]     algorithm to use (default is %s)\n', dict.alg);
    fprintf('                      Available values: [%s]\n', join(algorithms,", "));
    fprintf('    usv [bool]    flag to use single DMS separation voltage\n');
    fprintf('                      (default is %s, i.e. full scan)\n', string(dict.usv));
    fprintf('\n');
    fprintf('  Example:\n');
    fprintf('    smop_ml ip=127.0.0.1 f=0.6 usv=true\n');   
    fprintf('\n\n');

    for jj = 1:n
        p = split(string(args(jj)),"=");
        if (length(p) ~= 2)
            continue;
        end

        try
            if strcmpi(p(1), "ip")
                dict.ip = p(2);
            elseif strcmpi(p(1), "cr")
                dict.cr = str2double(p(2));
            elseif strcmpi(p(1), "f")
                dict.f = str2double(p(2));
            elseif strcmpi(p(1), "dec")
                dict.dec = str2double(p(2));
            elseif strcmpi(p(1), "mi")
                dict.mi = str2double(p(2));
            elseif strcmpi(p(1), "th")
                dict.th = str2double(p(2));
            elseif strcmpi(p(1), "alg")
                dict.th = p(2);
            elseif strcmpi(p(1), "usv")
                dict.usv = strcmpi(p(2), 'true') || strcmp(p(2), "1");
            end
        catch ex
            fprintf("Failed to parse command-line argument '%s': %s", ...
                args(jj), ex.message);
        end
    end
end