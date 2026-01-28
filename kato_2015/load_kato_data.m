function worms = load_kato_data(matPath, maxWorms)
%LOAD_KATO_DATA Load Kato 2015 WT_NoStim dataset into per-worm structs.
%
% Usage:
%   worms = load_kato_data();
%   worms = load_kato_data('KATO_WT_NoStim.mat', 3);
%
% Output per worm:
%   - uid (worm index as string)
%   - RawData (neurons x time, LOOPER-ready)
%   - X (time x neurons)
%   - T, num_neurons
%   - dt_sec, time
%   - neuron_ids (if available)
%   - source_dataset

    if nargin < 1 || isempty(matPath)
        matPath = fullfile(fileparts(mfilename('fullpath')), 'KATO_WT_NoStim.mat');
    end
    if nargin < 2
        maxWorms = [];
    end

    S = load(matPath);
    S = S.WT_NoStim;

    tracesCell = normalize_traces_local_(field_collect_(S, 'deltaFOverF_bc'));
    nWorms = numel(tracesCell);

    timeCell = normalize_time_cell_(field_collect_(S, 'tv'), nWorms);
    idsCell = normalize_ids_cell_(field_collect_(S, 'NeuronNames'), nWorms);
    datasetCell = normalize_ids_cell_(field_collect_(S, 'dataset'), nWorms);
    fpsCell = normalize_ids_cell_(field_collect_(S, 'fps'), nWorms);

    wormList = 1:nWorms;
    if ~isempty(maxWorms)
        wormList = wormList(1:min(numel(wormList), maxWorms));
    end

    worms = repmat(struct, numel(wormList), 1);
    wi_out = 1;
    for wi = wormList
        Xraw = tracesCell{wi};
        tvec = timeCell{wi};
        fps_i = [];
        if ~isempty(fpsCell{wi})
            fps_i = double(fpsCell{wi});
        end

        [X_tn, tvec, dt] = coerce_trace_to_t_by_n_(Xraw, tvec, fps_i);

        worms(wi_out).uid = string(wi);
        worms(wi_out).X = X_tn;
        worms(wi_out).RawData = X_tn'; % neurons x time
        worms(wi_out).T = size(X_tn, 1);
        worms(wi_out).num_neurons = size(X_tn, 2);
        worms(wi_out).dt_sec = dt;
        worms(wi_out).time = tvec;
        worms(wi_out).neuron_ids = idsCell{wi};
        worms(wi_out).source_dataset = datasetCell{wi};
        worms(wi_out).mat_path = string(matPath);

        wi_out = wi_out + 1;
    end
end

% ---- local helpers ----
function tracesCell = normalize_traces_local_(traces)
    if iscell(traces)
        tracesCell = traces;
        return;
    end
    if isnumeric(traces)
        if ndims(traces) == 2
            tracesCell = {traces};
            return;
        end
        if ndims(traces) == 3
            nW = size(traces, 3);
            tracesCell = cell(1, nW);
            for i = 1:nW
                tracesCell{i} = traces(:, :, i);
            end
            return;
        end
    end
    error('Unsupported traces format. Expected cell or numeric array.');
end

function val = field_collect_(S, name)
    if ~isfield(S, name)
        val = [];
        return;
    end
    if numel(S) > 1
        val = {S.(name)};
    else
        val = S.(name);
    end
end

function timeCell = normalize_time_cell_(timeData, nWorms)
    timeCell = cell(1, nWorms);
    if isempty(timeData)
        return;
    end
    if iscell(timeData)
        timeCell(1:min(nWorms, numel(timeData))) = timeData(1:min(nWorms, numel(timeData)));
        return;
    end
    if isnumeric(timeData)
        if isvector(timeData)
            for i = 1:nWorms
                timeCell{i} = timeData(:);
            end
            return;
        end
        if size(timeData,2) == nWorms
            for i = 1:nWorms
                timeCell{i} = timeData(:, i);
            end
            return;
        end
        if size(timeData,1) == nWorms
            for i = 1:nWorms
                timeCell{i} = timeData(i, :)';
            end
            return;
        end
        for i = 1:nWorms
            timeCell{i} = timeData(:);
        end
    end
end

function idsCell = normalize_ids_cell_(idsData, nWorms)
    idsCell = cell(1, nWorms);
    if isempty(idsData)
        return;
    end
    if iscell(idsData)
        idsCell(1:min(nWorms, numel(idsData))) = idsData(1:min(nWorms, numel(idsData)));
        return;
    end
    for i = 1:nWorms
        idsCell{i} = idsData;
    end
end

function [X_tn, tvec, dt] = coerce_trace_to_t_by_n_(Xraw, tvecIn, fps)
    Xraw = double(Xraw);
    tvec = [];
    if ~isempty(tvecIn)
        tvec = double(tvecIn(:));
    end

    if ~isempty(tvec) && numel(tvec) == size(Xraw,1)
        X_tn = Xraw; % time x neuron
    elseif ~isempty(tvec) && numel(tvec) == size(Xraw,2)
        X_tn = Xraw'; % neuron x time -> time x neuron
    else
        if size(Xraw,1) >= size(Xraw,2)
            X_tn = Xraw;
        else
            X_tn = Xraw';
        end
    end

    if isempty(tvec)
        tlen = size(X_tn,1);
        if ~isempty(fps) && isfinite(fps) && fps > 0
            dt = 1 / fps;
        else
            dt = 1;
        end
        tvec = (0:tlen-1)' * dt;
    else
        tlen = min(numel(tvec), size(X_tn,1));
        tvec = tvec(1:tlen);
        X_tn = X_tn(1:tlen, :);
        dt = median(diff(tvec));
    end
end
