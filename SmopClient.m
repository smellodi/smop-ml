%% SMOP Client class

% Handles communication with the SMOP server
classdef SmopClient < handle
    
    properties (SetAccess = private)
        isConnected = false
        config = 0
        initialDMS = 0
        initialSNT = 0
        ip
        port
        gases
    end
    
    properties (Access = private)
        socket
    end
    
    methods
        % Constructor
        function obj = SmopClient(ip, port)
            arguments
                ip = "127.0.0.1"
                port = 2339
            end

            obj.ip = ip;
            obj.port = port;

            try
                obj.socket = tcpclient(ip, port, "Timeout", 3);
                obj.isConnected = true;
            catch ex
                fprintf("Cannot connect to %s:%d >> %s\n", ip, port, ex.message);
            end
        end
        
        % Reads config packet and, possibly, the initial DMS packet that may arrive
        % together with the config packet. Returns "true" if the config received,
        % "false" otherwise
        % Exceptions: network stream reading
        function isOk = readConfig(obj)
            isOk = false;

            str = obj.readSocket();
            lines = splitlines(str);

            line1 = trim(lines{1});
            if isempty(line1)
                return;
            end

            try
                obj.config = toConfig(line1);
                c = length(obj.config.printer.channels);
                obj.gases = strings(1,c);
                for ii = 1:c
                    obj.gases(ii) = obj.config.printer.channels(ii).odor;
                end
            catch ex
                fprintf("Failed to parse config: %s\n", ex.message);
                return;
            end
          
            % Lets check if the initial measurement was already received
            if (length(lines) > 1) && (~isempty(trim(lines{2})))
                % yes it was, lets memorize it
                mline = trim(lines{2});

                try
                    obj.initialDMS = toMeasurement(mline, "dms");
                catch ex
                    fprintf("Failed to parse the initial DMS: %s\n", ex.message);
                end

                try
                    obj.initialSNT = toMeasurement(mline, "snt");
                catch ex
                    fprintf("Failed to parse the initial SNT: %s\n", ex.message);
                end
            end

            isOk = true;
        end

        % Return the initial measurement received with the config packet, or 
        % waits until it is received if measurement packet wasn't received yet.
        function measrm = getInitialMeasurement(obj)
            measrm = struct( ...
                "dms", 0, ...
                "snt", 0 ...
            );

            if isstruct(obj.initialDMS)  % if we have the initial DMS received
                measrm.dms = obj.initialDMS; % with the config packet, then 
                return;                  % lets return it immediately
            end

            if isstruct(obj.initialSNT)  % if we have the initial SNT received
                measrm.snt = obj.initialSNT; % with the config packet, then 
                return;                  % lets return it immediately
            end

            measrm = obj.waitForMeasurement();  % otherwise lets read 
                                                % the stream again
            if isstruct(measrm.dms)
                obj.initialDMS = measrm.dms;
            end
            if isstruct(measrm.snt)
                obj.initialSNT = measrm.snt;
            end
        end
        
        % Reads the network stream, returns DMS, SNT and PID measurements if any.
        % Note that only one DMS measurement can be returned (or 0 if DMS was not 
        % received), while the SNT and PID measurements returned as arrays (can 
        % be empty).
        % Exceptions: network stream reading, JSON parsing, packet format
        function [dms, snts, pids] = readData(obj)
          str = obj.readSocket();
          lines = splitlines(str);
          dms = 0; snts = []; pids = [];
        
          for i = 1:length(lines)
            line = trim(lines{i});
            if ~isempty(line)
                dms = toMeasurement(line, "dms");
                if isstruct(dms)
                    continue;
                end
        
                snt = toMeasurement(line, "snt");
                if (isstruct(snt))
                    snts = [snts snt];
                    continue;
                end
        
                pid = toMeasurement(line, "pid");
                if (isstruct(pid))
                    pids = [pids pid];
                end
            end
          end
        end
        
        % Does not return until a measurement is received. 
        % Currently, only DMS are are awaited.
        function measrm = waitForMeasurement(obj)
          dms = 0;
          snt = 0;

          while ~isstruct(dms) && ~isstruct(snt)
            try
                [dms, snts, ~] = obj.readData();
                if ~isempty(snts)
                    snt = snts(1);
                end
            catch ex
                fprintf("Failed to parse the packet: %s\n", ex.message);
            end
          end

          measrm = struct( ...
              "dms", dms, ...
              "snt", snt ...
          );
        end
        
        % Sends the recipe to SMOP.
        % Exceptions: network stream reading
        function sendRecipe(obj, name, flows, isFinished, distance)
          arguments
            obj
            name
            flows
            isFinished = false
            distance = 1e5
          end
        
          [offsets, gains] = obj.getFlowTransformations();

          c = length(flows);
          channels(1:c) = struct;
          for i = 1:c
              channels(i).id = obj.config.printer.channels(i).id;
              channels(i).flow = offsets(i) + gains(i) * flows(i);
              channels(i).duration = -1;
          end
        
          recipe = struct( ...
              "type", "recipe", ...
              "content", struct( ...
                "name", name, ...
                "isFinal", isFinished, ...
                "distance", distance, ...
                "channels", channels ...
              ));

          if c == 1
              recipe.content.channels = {channels};
          end
        
          json = jsonencode(recipe);
          dataBytes = uint8(json);
          dataBytes = [dataBytes 13 10];
          write(obj.socket, dataBytes);
        end

        % Config parameters

        % As of Dec 2024, props of printer channels have the following fields:
        %   minFlow, maxFlow, criticalFlow, pidCheckLevel, 
        %   shortKnownName, fullKnownName, isKnownOdor

        function maxFlow = getMaxFlow(obj, default)
            maxFlow = default;
            
            % The previous way is not valid anymore: we want to keep 
            % individual ranges for each odor, and this info comes
            % as minFlow and maxFlow of each printer channel

            % if ~isstruct(obj.config)
            %     return
            % end
            % 
            % % We pretend that the max flow for all channels is same, so we
            % % take the first channel's max flow value and return it.
            % firstChannel = obj.config.printer.channels(1);
            % if isfield(firstChannel.props, "maxFlow")
            %     maxFlow = firstChannel.props.maxFlow;
            % end
        end

        function channelCount = getChannelCount(obj, minimum)
            channelCount = minimum;

            if isstruct(obj.config)
                channelCount = length(obj.config.printer.channels);
                if channelCount < minimum
                    msg = sprintf("There must be at least %d channels" + ...
                        " enabled in Odor Printer", minimum);
                    throw(MException("smop:client", msg));
                end
            end
        end

        function maxIterCount = getMaxIterationCount(obj, default)
            maxIterCount = default;
            if isstruct(obj.config) && (obj.config.maxIterationNumber > 0)
                maxIterCount = obj.config.maxIterationNumber;
            end
        end

        function threshold = getThreshold(obj, default)
            threshold = default;
            if isstruct(obj.config) && (obj.config.threshold > 0)
                threshold = obj.config.threshold;
            end
        end

        function algorithm = getAlgorithm(obj, default)
            algorithm = default;
            if isstruct(obj.config) && ~isempty(obj.config.algorithm)
                algorithm = obj.config.algorithm;
            end
        end

        function usv = getUsv(obj)
            usv = 0;
            if isstruct(obj.initialDMS) && ~isstruct(obj.initialDMS.setup.usv)
                usv = obj.initialDMS.setup.usv;
            end
        end

        function criticalFlows = getCriticalFlows(obj)
            channels = obj.config.printer.channels;
            c = length(channels);
            criticalFlows = ones(c,1) * 100;
            for jj = 1:c
                if (isfield(channels(jj).props, "criticalFlow"))
                    criticalFlows(jj) = channels(jj).props.criticalFlow;
                end
            end
        end
    end

    methods (Access = private)
        % Destructor
        function delete(obj)
            clear obj.socket;
        end

        % Reads network stream data, returns it as a string
        function data = readSocket(obj)
          dataBytes = [];
          while (isempty(dataBytes))
            pause(1);
            dataBytes = read(obj.socket);
            % remove leading zeros
            dataBytes = dataBytes(find(dataBytes ~= 0, 1, "first") : end);
          end
          data = native2unicode(dataBytes);
        end

        % Use it only for the study in Dec 2024 
        function [offsets, gains] = getFlowTransformations(obj)
            c = length(obj.config.printer.channels);
            offsets = zeros();
            gains = ones();

            for ii = 1:c
                props = obj.config.printer.channels(ii).props;
                
                offsets(ii) = props.minFlow;
                gains(ii) = (props.maxFlow - props.minFlow) / 50;

                % offset = 0;
                % gain = 1;
                
                % odorName = obj.config.printer.channels(ii).odor;
                % if odorName == "Cyclohexanone"  % 0 .. 2.8
                %     gain = 2.8 / 50;
                % elseif odorName == "Limonene"   % 5 .. 40
                %     offset = 5;
                %     gain = 40.0 / 50.0;
                % elseif odorName == "Citronellol"   % 40 .. 100
                %     offset = 50;
                %     gain = 60.0 / 50.0;
                % end

                % offsets(ii) = offset;
                % gains(ii) = gain;
            end
        end
    end
end

%% STATIC private function

% Removes leading blank characters from the string.
% The leading blank chars may appear in the packets, as the SMOP server sends
% \0 byte every second to check the connection status. These null-bytes are
% then converted into " " by native2unicode, and then jsondecode fails
% to parse the string if it has leading blanks.
function s = trim(text)
  s = text(find(text > 32, 1, "first") : end);
end

% Returns either a packet content if its type matches the desired type,
% or 0 is does not match, or throws an exception if the packet is not a
% valid packet.
function content = getPacketContent(packet, type)
    content = 0;
    isValidJson = all(isfield(packet, ["type", "content"]));
    if ~isValidJson
        throw(MException("smop:client", "invalid packet"))
    elseif strcmp(packet.type, type)
        content = packet.content;
    end
end

% Returns content of the config packet.
% Note how we define some propertes right HERE!
% (or should these properties be defined in SMOP?)
function config = toConfig(json)
  packet = jsondecode(json);
  config = getPacketContent(packet, "config");
end

% Returns content of a measurement packet, or 0 if this measurement source
% is not the one that is expected.
function measurement = toMeasurement(json, source)
  packet = jsondecode(json);
  measurement = getPacketContent(packet, "measurement");

  if ~isstruct(measurement) || ~strcmp(measurement.source, source)
      measurement = 0;
  end
end

%% NET Protocol

% Copied here for convenience.
% Please refer the online document to keep it up-to-date

% { "type": "config",
%   "content" : {
%     "sources": ["dms" | "snt" | "pid"],
%     "algorithm": string,       // algorithm name,  "" - default
%     "maxIterationNumber": int, // 0 - ignore
%     "threshold": float, // diff between iterations; 0 – ignore
%     "printer": {
%       "channels": [{   // an array of the device channels
%         "id": int,	 // channel ID
%         "odor": string, // name of the odor
%         "props": {      // list of properties
%            “maxflow”: float,		// max allowed flow in nccm
%            “criticalFlow”: float		// PID oversaturating flow
%         }
%       }]
%     }
% } }

% { "type": "measurement",
%   "content" : {
%     "source": "dms",
%     "setup": { // from IV.Param. MeasurementParameters.SteppingControl
%       "usv": { min, max, step }, OR <value> in single-Us mode
%       "ucv": { min, max, step }
%     },
%     "conditions": { // from IV.Scan.SystemData.Sample  - OPTIONAL
%       "flow": { avg, min max },
%       "temperature": { avg, min max },
%       "pressure": { avg, min max },
%       "humidity": { avg, min max },
%       "pumppwm": { avg, min max }
%     },
%     "data": {
%       "positive": [float], // IV.Scan.MeasurementData.IntensityTop
%       "negative": [float]  // IV.Scan.MeasurementData.IntensityBottom
%     }
% } }

% { "type": "measurement",
%   "content" : {
%     "source": "snt",
%     "data": {
%       "resistors": [float],	// 64 values
%       "temperature": float,
%       "humidity": float
%     }
% } }

% { "type": "measurement",  // CURRENTLY IGNORED
%   "content" : {
%     "source": "pid",
%     "data": float
% } }

% { "type": "recipe",
%   "content" : {
%     "name": string,     // whatever name to display in SMOP UI
%     "isFinal": bool,
%     "distance": float,
%     "humidity": float?,  // % - OPTIONAL
%     "channels": [{      // an array of each channel configuration
%       "id": int,        // 1..5
%       "flow": float,    // sccm
%       "duration": float,     // seconds, or -1, or 0	
%       "temperature": float?, // °C (will be available in future)
%     }]
% } }
