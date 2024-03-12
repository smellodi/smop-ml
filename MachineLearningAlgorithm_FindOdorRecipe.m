
% TO DO: ADD ALSO GLOBAL MINIMUM TO MESSAGE SEND TO OLEG
% TIMO SUGGEST TO USE FLOW RATES OF 20 OR HIGHER FOR TARGET SCENT TO SPEED
% UP RESPSONSE TIME/SHORTEN STABILIZATION TIME
% USE LOWER RESOLUTION --> MEASURE ONLY REGIONS WITH ALPHA CURVES FROM
% ANTON'S PAPER: SV 450 to 700, CV -1 to 4
% for future implementation I would still suggest to use Atte's training
% base as starting point to speed up search


% Script for the differential evolution algorithm described in Mariaana's
% manuscript. It returns a two-dimensional vector x that contains the flow
% rates for isopropanol (1st component) and ethanol (2nd component).
% Accepted values for flow rates are 0, 10, 20, .., 50 sccm.
%
% IF WE NOW WANT TO ALSO INCLUDE SMELL INSPECTOR AND PID SENSOR THE WHOLE
% SYSTEM WILL NEED TO BE REDESIGNED
%
% See TG-Smellodi/Documents/General/Software/algorithm_Smellodi.docx for
% details on structure of input from and output to Oleg's software (will
% use 'tcpclient' and send JSON files).
%
% THIS VERSION uses simplified message from and to Oleg's software:
% - received messages: [SV_steps, CV_steps, vector_DMS_values]
% - sent messages: [final_switch, flowrate_channel1, flowrate_channel2]
%   with final_switch = 1 if odor recipe ready and 0 otherwise
%
% Philipp Muller, January 2024

% OLEG, March 2024:
%
% function [w, h, dms] = readDMS(socket)
%   dataBytes = [];
%   while (~length(dataBytes)) || (dataBytes(1) ~= 91)  % 91 is '[' symbol
%     pause(1)
%     dataBytes = read(socket);
%   end
%   json = native2unicode(dataBytes);
%   array = jsondecode(json);
%   w = array(1)
%   h = array(2)
%   dms = array(3:end);
% end
%
% Usage example:
% [~, ~, posDP] = readDMS(ComIn)
%
%
% function sendRecipe(socket, isFinished, flowIPA, flowNBut, LM)
%   recipe = [isFinished, flowIPA, flowNBut, LM]
%   json = jsonencode(recipe)
%   dataBytes = uint8(json);
%   dataBytes = [dataBytes 13 10];
%   write(socket, dataBytes) 
% end
%
% Usage example:
% sendRecipe(ComIn, isFinished, X(1,ii)*10, X(2,ii)*10, 1e5)

clear all
close all

dbstop error

%rng(0)

%% STEP 0: Define parameters etc.
% - IP and PORT for communication with Oleg's software
IP = '192.168.1.3';
PORT = 2339;
% - create a TCP connection object and connect to Oleg's software

% UNCOMMENT FOLLOWING LINE FOR DEMO
ComIn = tcpclient(IP, PORT);
% % send some test message to Oleg's software to check if he can receive
% a = [0, 1, 2];
% recipeTxt = jsonencode(a);
% recipeBytes = uint8(recipeTxt);
% recipeBytes = [recipeBytes 13 10];
% write(ComIn, recipeBytes) 

% - switch for terminating iterative process
isFinished = 0;

% - number of chemicals in the mixture (concentration can be also 0 sccm)
n = 2;

% - maximum flow of each chemical in the mixture (in sccm)
max_f = 50;         % default value for testing: 50

% - maximum number of iterations
mi = 5;

% - threshold RMSE for terminating iterative process
limRMSE = 5;%0.2;

% Define which distance measure to use. Options are:
% - 'euclidean': uses Euclidean distance between two matrices
% - 'alpha': extracts alpha curves or dispersion plots using Anton's
%   algorithm (ISOEN paper) and computes Euclidean distance between the
%   alpha curve locations
% - what else could be tested by Naa???

% population size: typically 10 times the dimension n of x but we know that
% there are 6 different flow rates for each chemical in the final mixtures,
% thus there are 6^2 options (choosing 1 out of 6 twice in a row, order
% matters, with repetition).
%np = 6^2;
% % matrix with concentrations of n-Butanol (1st row) and isopropanol
% C = [zeros(1,6) ones(1,6) 2*ones(1,6) 3*ones(1,6) 4*ones(1,6) ...
%     5*ones(1,6); repmat(0:1:5,1,6)];

% crossover population CR: [0,1]
% - Mariaana used cr = 0.5 for n <= 3 and cr = 0.9 for n >= 4
% - she pointed out that cr = 0.9 could have been used for any n (also
%   Wikipedia lists 0.9 as standard setting)
% - cr maintains diversity of population
% - for non-separable problems with correlated parameters cr in [0.5, 1] is
%   recommended
cr = 0.9;%0.5;

% differential weight f: [0,2] (Wikipedia lists 0.8 as standard setting)
% - f controls differential variation and has most severe impact on
%   performance of differential evolution algorithm
% 
f = 0.8;%0.5;

%% STEP 1: Load DMS database
% We need database with dispersion plots from n-Butanol and isopropanol as
% well as their mixtures measured at flow rates of 0, 10, 20, .., 50 sccm.

% FOR NOW USE THE MAT FILE FROM ATTE'S TESTS in 2023
%load('DMSmeasIPAandNBut')

%% STEP 2: Read in message with DMS measurement from target scent

hasData = 0;
dataBytes = [];
    
while ~hasData
	pause(1)
    dataBytes = read(ComIn);
    hasData = length(dataBytes);
end
TS = native2unicode(dataBytes);
TS = jsondecode(TS);
posDP = TS(3:end);

%% STEP 3: Initialization

% % OPTION A: Use previously measured dispersion plots that contain same
% % chemicals at same possible flow rates.
% % Initialize a population of np n-dimensional parameter vectors that
% % contain concentration for each of the n chemicals in the mixture.
% X = randi([0,5],n,np);
% 
% % Find indices of DMS measurements in 'allIMS' with the flow levels defined
% % in X.
% % - extract flow rates for IPA and n-Butanol from 'allIMS'
% cIPA = [allIMS.ConcIPA];
% cNBut = [allIMS.ConcNBut];
% p = nan(np,1);
% for jj = 1:np
%     % find samples for which X and 'allIMS' have same flow rate for ..
%     % .. isopropanol and ..
%     idIPA = find(cIPA==X(1,jj)*10);
%     % .. n-Butanol
%     idNBut = find(cNBut==X(2,jj)*10);
%     
%     % get intersection and shuffle it in random order
%     % NOTE: some pairs of flow rates might be missing from the data
%     idMix = intersect(idIPA,idNBut);
%     idMix = idMix(randperm(numel(idMix)));
%     
%     % use DMS measurement related to first index for initial population
%     if ~isempty(idMix)
%         p(jj) = idMix(1);
%     end
%     clear idIPA idNBut idMix
% end

% OPTION B: Implement this version and use it for the demo



%start without any prior information --> measure np mixtures
% % - where do we define which flow rates we measure for this case?

% is [0 0] a realistic option? or should we use only [max_f 0], [0 max_f],
% [max_f max_f], and [max_f/2 max_f/2]

% % - how many dispersion plots do we measure before we run iterative step?

%X = .1 * [0 max_f 0 max_f max_f/2; 0 0 max_f max_f max_f/2];
% use the line below for testing only --> measures 3 mixtures only
X = .1 * [0 max_f 0 .8*max_f max_f/2; 0 0 max_f .8*max_f max_f/2];


np = size(X,2);
cIPA = []; cEth = []; p = nan(np,1);
for ii = 1: np
    disp(['Measure reference dispersion plot ' num2str(ii) '.'])
    recipe = [isFinished, X(1,ii)*10, X(2,ii)*10, 1e5]
    recipeTxt = jsonencode(recipe)
    recipeBytes = uint8(recipeTxt);
    recipeBytes = [recipeBytes 13 10];
    write(ComIn, recipeBytes) 
    % wait for new DMS measurement and add it to the table of measured
    % recipes
    hasData = 0;
    dataBytes = 1; %[];
    while (~hasData) || (dataBytes(1)~=91)
    	pause(1)
        dataBytes = read(ComIn);
        hasData = length(dataBytes);
    end
    TS = native2unicode(dataBytes);
    TS = jsondecode(TS);
    allIMS(ii).MeasurementData.IntensityTop = TS(3:end);
    
    cIPA = [cIPA X(1,ii)*10];
    cEth = [cEth X(2,ii)*10];
    
    clear recipe recipeTxt recipeBytes
end    

% following loop might be not necessary
for jj = 1:np
    % find samples for which X and 'allIMS' have same flow rate for ..
    % .. isopropanol and ..
    idIPA = find(cIPA==X(1,jj)*10);
    % .. n-Butanol
    idEth = find(cEth==X(2,jj)*10);
    
    % get intersection and shuffle it in random order
    % NOTE: some pairs of flow rates might be missing from the data
    idMix = intersect(idIPA,idEth);
    idMix = idMix(randperm(numel(idMix)));
    
    % use DMS measurement related to first index for initial population
    if ~isempty(idMix)
        p(jj) = idMix(1);
    end
    clear idIPA idNBut idMix
end


%% STEP 4: Iterative step
GM = 1e8;       % overall minimum RMSE
cIt = 1;        % counter iterations

while(isFinished==0)
    % minimum RMSE of current iteration
    LM = 1e8;

    % mutation
    V = nan(n,np);
    for jj = 1:np
        % - pick three distinct vectors from X that are different from jj
        ids = randperm(np);     % random permutation of integers 1 to np
        ids(ids==jj)=[];        % remove integer jj
        
        % - compute donor vector from first three vectors
        V(:,jj) = X(:,ids(1)) + f * (X(:,ids(2)) - X(:,ids(3)));
    end
    % PROBLEM FOR NOW THAT WE HAVE ONLY FLOW RATES 0, 10, 20, 50
    % IF THE NEW CANDIDATE IS, E.G., 14 THEN ROUND IT TO NEAREST FLOW RATE,
    % I.E. 10
    % THE FOLLOWING LINE CAN BE REMOVED LATER ON, ONCE WE CAN FREELY
    % ADJUST THE CONCENTRATIONS
    
    
% Rounding of flow rates might be a good idea if we only use 0, 10, 20, .. 
% for the flow rates.
    %V = round(V);
    
    % CHECK IF SOME VALUE IS NEGATIVE (NOT POSSIBLE) or larger than 5 (too
    % high flow rate)
    V(V<0) = 0;
    V(V>.1*max_f) = .1*max_f;
        
    
    % crossover
    % use binomial method for combining components from target and donor
    % vectors
    %  - generate randomly chosen indices to ensure that at least one
    %    component of the donor vector is included in the target vector
    rv = randi(2,1,np);
    % - random numbers in [0,1] for each component of each target vector
    rv2 = rand(n,np);
    % - combining target and donor vectors to get trial vectors
    U = nan(n,np);
    for jj = 1:np
        for kk = 1:n
            if (rv2(kk,jj) <= cr) || (rv(jj) == kk)
%                disp(['Component ' num2str(kk) ' of vector ' num2str(jj) ' is changed.'])
                U(kk,jj) = V(kk,jj);
            else
                U(kk,jj) = X(kk,jj);
            end
        end
    end
    
    
    % selection
    % - compute cost function values for trial vectors and corresponding
    %   target vectors
    % COST FUNCTION WILL BE FUNCTION THAT MEASURES SIMILARITY OF DISPERSION
    % PLOT OF MIXTURE WE WANT TO RECREATE AND ONE TRAINING DISPERSION PLOT
    % NAA WILL TEST DIFFERENT MEASURES IN HER THESIS
    % 
    % for NOW USE RMSE AS DISTANCE MEASURE
       
    for jj = 1:np
        
        if ~isnan(p(jj))
            % calculate distances for pairs of dispersion plots, PID curves and
            % SmellInspector curves between candidate solution and original
            % measured scent
            % FOR NOW USE ONLY DMS DATA, i.e. calculate RMSEs of differences
            % between vector pairs

            % statistical tests on differences between generated and original
            % odors?

            % calculate distances for pairs of dispersion plots, PID curves and
            % SmellInspector curves between candidate solution and original
            % measured scent
            % FOR NOW USE ONLY DMS DATA, i.e. calculate RMSEs of differences
            % between vector pairs

            % statistical tests on differences between generated and original
            % odors?
        
        
            % calculate root mean square error of target scent and potential
            % artifical scent defined by target vector
            tv = allIMS(p(jj)).MeasurementData.IntensityTop;
%             tv = allIMS(jj).MeasurementData.IntensityTop;
           
            
%      % FOR TESTING ONLY
%      tv = tv(1:12);
            
            cf_X = sqrt(mean((posDP - tv).^2));
            clear tv
            
            % ALTERNATIVES to evaluate closeness of dispersion plots:
            % 1. Euclidean distance
            % 2. Frobenius norm
            % 3. Cosine similarity
            % 4. Extract alpha curves and calculate their similarities

            % calculate root mean square error of target scent and potential
            % artifical scent defined by trial vector
            % - find samples for which U and 'allIMS' have same flow rate for ..
            % .. isopropanol and ..
            idIPA = find(cIPA==U(1,jj)*10);
            % .. ethanol
            idEth = find(cEth==U(2,jj)*10);
            % - get intersection and shuffle it in random order
            idMix = intersect(idIPA,idEth);
            idMix = idMix(randperm(numel(idMix)));
            % - calculate RMSE
            if ~isempty(idMix)
                tv = allIMS(idMix(1)).MeasurementData.IntensityTop;

%      % FOR TESTING ONLY
%      tv = tv(1:12);
     
                cf_U = sqrt(mean((posDP - tv).^2));

                clear tv idIPA idNBut idMix
                
                % storing global and local minima
                if min(cf_U,cf_X) < GM
                    GM = min(cf_U,cf_X);
                    idGM = jj;
                end
                if min(cf_U,cf_X) < LM
                    LM = min(cf_U,cf_X);
                    idLM = jj;
                end
                % if trial vector yields lower RMSE than target vector then replace
                % target by trial vector
                if cf_U < cf_X
                    X(:,jj) = U(:,jj);
                    % no need for else clause, because it will not change X
                end
                clear cf_U
            else
                % if set of concentrations not measured yet then ask
                % IonVision to measure it
                recipe = [isFinished, U(1,jj)*10, U(2,jj)*10, LM]
                recipeTxt = jsonencode(recipe)
                recipeBytes = uint8(recipeTxt);
                recipeBytes = [recipeBytes 13 10];
                write(ComIn, recipeBytes) 
                % wait for new DMS measurement and add it to the table of measured
                % recipes
                hasData = 0;
                dataBytes = 1; %[];
                while (~hasData) || (dataBytes(1)~=91)
                    pause(1)
                    dataBytes = read(ComIn);
                    hasData = length(dataBytes);
                end
                TS = native2unicode(dataBytes);
                TS = jsondecode(TS);
                allIMS(length(allIMS)+1).MeasurementData.IntensityTop = TS(3:end);
                tv = TS(3:end);
                cf_U = sqrt(mean((posDP - tv).^2));

                clear tv idIPA idNBut idMix
                
                % storing global and local minima
                if min(cf_U,cf_X) < GM
                    GM = min(cf_U,cf_X);
                    idGM = jj;
                end
                if min(cf_U,cf_X) < LM
                    LM = min(cf_U,cf_X);
                    idLM = jj;
                end
                % if trial vector yields lower RMSE than target vector then replace
                % target by trial vector
                if cf_U < cf_X
                    X(:,jj) = U(:,jj);
                    % no need for else clause, because it will not change X
                end
                clear cf_U                
            end
            clear cf_X
        end
    end
    
    % pick three distinct dispersion plots from training set (need to be
    % distinct from each other and from x)
    % also pick corresponding PID and SmellInspector data
    % QUESTION: Do they need to be from different concentration
    % combinations or could they be also from the same concentration
    
    
    % How to decide which is the best recipe right now? store local and
    % global minima of RMSE between target scent and any potential
    % artificail scent as well as its column in X
    disp(['Minimum RMSE in iteration ' num2str(cIt) ': ' num2str(LM) ])
    disp(['- flow rate IPA: ' num2str(X(1,idLM)*10)]);
    disp(['- flow rate n-Butanol: ' num2str(X(2,idLM)*10)]);
    
    
    cIt = cIt + 1;  % increase counter iterations by one

    % terminate iterative process or not
    if (LM < limRMSE) || (cIt >= mi)
%    if (cIt >= mi)
        isFinished = 1;
        % send message
        disp('Final recipe found:')
        recipe = [isFinished, X(1,idLM)*10, X(2,idLM)*10, LM]
        % UNCOMMENT 3 ROWS BELOW FOR DEMO
        recipeTxt = jsonencode(recipe);
        recipeBytes = uint8(recipeTxt);
        recipeBytes = [recipeBytes 13 10];
        write(ComIn, recipeBytes)    
    else
        % send message
        disp('New recipe should be tested:')
        recipe = [isFinished, X(1,idLM)*10, X(2,idLM)*10, LM]
        recipeTxt = jsonencode(recipe)
        recipeBytes = uint8(recipeTxt);
        recipeBytes = [recipeBytes 13 10];
        write(ComIn, recipeBytes) 
        % wait for new DMS measurement and add it to the table of measured
        % recipes
        hasData = 0;
        dataBytes = 1; %[];
        while (~hasData) || (dataBytes(1)~=91)
            pause(1)
            dataBytes = read(ComIn);
            hasData = length(dataBytes);
        end
        TS = native2unicode(dataBytes);
        TS = jsondecode(TS);
        allIMS(length(allIMS)+1).MeasurementData.IntensityTop = TS(3:end);
        
        clear recipe recipeBytes recipeTxt
    end
    clear LM TS
end

%% Disconnect once recipe is ready
clear ComIn
