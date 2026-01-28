% Estimate a delay window from autocorrelation in Atanas baseline worms.
% eLife 2019 (Brennan & Proekt) choose delay t based on the *minimum across
% neurons* autocorrelation time constant (1/e crossing), for activity and
% its derivative.

addpath(fullfile(pwd, 'atanas-data'));
addpath(fullfile(pwd, 'LOOPER_github_2020', 'Functions', 'Library'));

MAX_LAG = 200; % frames (adjust if needed)

worms = load_atanas_data('baseline');

lag1e_act_all = nan(1, numel(worms));
lag1e_deriv_all = nan(1, numel(worms));

for wi = 1:numel(worms)
    worm = worms(wi);
    X = worm.RawData; % neurons x time
    if isempty(X)
        continue;
    end
    % Paper smoothing (sigma = 1); traces are already z-scored in trace_array.
    X = filterData(X, 1, 'gaussian', true, 1);
    X = zscore(X, 0, 2);

    maxLag = min(MAX_LAG, size(X,2) - 2);
    actLag1e = [];
    derivLag1e = [];

    for ni = 1:size(X,1)
        x = X(ni, :);
        x = x - mean(x);
        ac = xcorr(x, maxLag, 'coeff');
        acPos = ac(maxLag+1:end); % lags 0..maxLag
        idx1e = find(acPos <= exp(-1), 1, 'first');
        if ~isempty(idx1e)
            actLag1e(end+1) = idx1e - 1; %#ok<AGROW>
        end

        dx = [0 diff(x)];
        dx = dx - mean(dx);
        acd = xcorr(dx, maxLag, 'coeff');
        acdPos = acd(maxLag+1:end);
        idx1ed = find(acdPos <= exp(-1), 1, 'first');
        if ~isempty(idx1ed)
            derivLag1e(end+1) = idx1ed - 1; %#ok<AGROW>
        end
    end

    dt = worm.dt_sec;
    if isempty(dt) || ~isfinite(dt)
        dt = 1;
    end

    if ~isempty(actLag1e)
        lag1e_act_all(wi) = min(actLag1e);
    end
    if ~isempty(derivLag1e)
        lag1e_deriv_all(wi) = min(derivLag1e);
    end

    fprintf('%s: lag1e_act_min=%d (%.2fs), lag1e_deriv_min=%d (%.2fs)\n', ...
        worm.uid, lag1e_act_all(wi), lag1e_act_all(wi)*dt, ...
        lag1e_deriv_all(wi), lag1e_deriv_all(wi)*dt);
end

fprintf('\nSummary across worms:\n');
if any(~isnan(lag1e_act_all))
    fprintf('lag1e_act_min median=%d frames\n', round(median(lag1e_act_all, 'omitnan')));
end
if any(~isnan(lag1e_deriv_all))
    fprintf('lag1e_deriv_min median=%d frames\n', round(median(lag1e_deriv_all, 'omitnan')));
end
