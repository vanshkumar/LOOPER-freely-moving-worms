function worms = load_atanas_data(split, varargin)
%LOAD_ATANAS_DATA Load per-worm JSON files from atanas-data/{baseline,heat}.
%
% Usage:
%   worms = load_atanas_data('baseline');
%   worms = load_atanas_data('heat', 'TraceKey', 'trace_original');
%
% Output:
%   worms is a struct array. Each entry includes:
%     - uid, path, dataset_type
%     - RawData (neurons x time) ready for LOOPER
%     - X (time x neurons) convenience transpose
%     - T, num_neurons, trace_key
%     - dt_sec, timestamp_confocal
%     - ranges (struct array with 1-based and 0-based indices)
%     - events (heat indices if present)
%     - behavior (selected fields), labeled, neuron_categorization

    p = inputParser;
    p.addRequired('split', @(s)ischar(s) || isstring(s));
    p.addParameter('TraceKey', 'trace_array', @(s)ischar(s) || isstring(s));
    p.parse(split, varargin{:});

    split = char(lower(p.Results.split));
    traceKey = char(p.Results.TraceKey);

    rootDir = fileparts(mfilename('fullpath'));
    switch split
        case 'baseline'
            dataDir = fullfile(rootDir, 'baseline');
        case 'heat'
            dataDir = fullfile(rootDir, 'heat');
        otherwise
            error('Unknown split: %s (expected baseline or heat)', split);
    end

    files = dir(fullfile(dataDir, '*.json'));
    files = sortfiles_(files);

    worms = repmat(struct, numel(files), 1);
    for i = 1:numel(files)
        path = fullfile(files(i).folder, files(i).name);
        raw = jsondecode(fileread(path));

        trace = raw.(traceKey);
        if isempty(trace)
            error('Missing trace key %s in %s', traceKey, path);
        end

        if iscell(trace)
            trace = cell2mat(trace);
        end

        numNeurons = size(trace, 1);
        T = size(trace, 2);

        worms(i).path = path;
        worms(i).uid = raw.uid;
        worms(i).dataset_type = raw.dataset_type;
        worms(i).trace_key = traceKey;
        worms(i).RawData = trace;          % neurons x time (LOOPER input)
        worms(i).X = trace.';              % time x neurons (convenience)
        worms(i).num_neurons = numNeurons;
        worms(i).T = T;
        worms(i).timestamp_confocal = raw.timestamp_confocal;
        worms(i).dt_sec = dt_seconds_(raw);
        worms(i).ranges_raw = raw.ranges;
        worms(i).ranges = ranges_to_segments_(raw.ranges);

        worms(i).events_raw = [];
        worms(i).events = struct;
        if isfield(raw, 'events')
            worms(i).events_raw = raw.events;
            if isfield(raw.events, 'heat') && ~isempty(raw.events.heat)
                heatIdx = raw.events.heat(1);
                worms(i).events.heat_1b = heatIdx;
                worms(i).events.heat_0b = heatIdx - 1;
            end
        end

        worms(i).behavior = struct;
        behaviorKeys = { ...
            'velocity', 'angular_velocity', 'head_curvature', 'body_curvature', ...
            'pumping', 'reversal_events', 'dorsalness', 'feedingness', ...
            'forwardness', 'rel_enc_str_v', 'rel_enc_str_θh', 'rel_enc_str_P', ...
            'tau_vals', 'encoding_changing_neurons' ...
        };
        for k = 1:numel(behaviorKeys)
            key = behaviorKeys{k};
            if isfield(raw, key)
                worms(i).behavior.(key) = raw.(key);
            end
        end

        if isfield(raw, 'labeled')
            worms(i).labeled = raw.labeled;
        else
            worms(i).labeled = [];
        end
        if isfield(raw, 'neuron_categorization')
            worms(i).neuron_categorization = raw.neuron_categorization;
        else
            worms(i).neuron_categorization = [];
        end
    end
end

function dt = dt_seconds_(raw)
    dt = [];
    if isfield(raw, 'timestamp_confocal') && numel(raw.timestamp_confocal) > 1
        ts = double(raw.timestamp_confocal(:));
        dt = mean(diff(ts));
    elseif isfield(raw, 'avg_timestep') && ~isempty(raw.avg_timestep)
        dt = double(raw.avg_timestep) * 60.0;
    end
end

function segments = ranges_to_segments_(ranges)
    segments = [];
    if isempty(ranges)
        return;
    end

    if isnumeric(ranges)
        % rows are segments
        for r = 1:size(ranges, 1)
            idx = ranges(r, :);
            idx = idx(:).';
            segments = [segments; segment_from_idx_(idx)]; %#ok<AGROW>
        end
    elseif iscell(ranges)
        for r = 1:numel(ranges)
            idx = ranges{r};
            if isempty(idx)
                continue;
            end
            segments = [segments; segment_from_idx_(idx(:).')]; %#ok<AGROW>
        end
    end
end

function seg = segment_from_idx_(idx)
    seg.start_1b = idx(1);
    seg.end_1b = idx(end);
    seg.start_0b = idx(1) - 1;
    seg.end_0b = idx(end) - 1;
    seg.length = numel(idx);
end

function files = sortfiles_(files)
    [~, idx] = sort({files.name});
    files = files(idx);
end
