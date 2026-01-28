% DEPRECATED: heat-pulse path not pursued (freely moving baseline did not recover loops).
% Kept for reference only; do not use for current experiments.
% Evaluate Experiment 1 single (heat) run: recovery metrics + LOOPER diagnostics.

addpath(fullfile(pwd, 'atanas-data'));
addpath(pwd); % experiment1_helpers, run_looper_diagnostics

% Ensure figures open as docked tabs (not separate windows).
set(0, 'DefaultFigureWindowStyle', 'docked');
set(0, 'DefaultFigureVisible', 'on');

% --- Settings (match run) ---
T_STAR_SEC = 10;

outDir = fullfile(pwd, 'results', 'experiment1_single');
matPath = fullfile(outDir, 'experiment1_single.mat');
if ~exist(matPath, 'file')
    error('No Experiment1 single results found. Run experiment1_single_run first.');
end

S = load(matPath);
if ~isfield(S, 'saveData') || ~isfield(S, 'worm')
    error('Expected worm and saveData in %s', matPath);
end

worm = S.worm;
saveData = S.saveData;

% Heat event index.
if ~isfield(worm, 'events') || ~isfield(worm.events, 'heat_1b')
    error('No heat event found in worm metadata. Re-run experiment1_single_run.');
end
t_on = worm.events.heat_1b;
dt = worm.dt_sec;

if isempty(worm.ranges_raw)
    error('Expected paper windows (ranges_raw) for Experiment 1.');
end
[preIdx, postIdx] = experiment1_helpers.normalize_ranges(worm.ranges_raw);
experiment1_helpers.warn_ranges_contiguous(preIdx, postIdx, worm.uid);

% Project full trace onto the learned scaffold (nearest emission).
rawFull = worm.RawData;
if isfield(saveData, 'Detrend') && isfield(saveData.Detrend, 'enabled') && saveData.Detrend.enabled
    rawFull = experiment1_helpers.detrend_apply(rawFull, saveData.Detrend);
end
[alpha, theta, ~, d] = experiment1_helpers.project_to_model(rawFull, saveData);

out = experiment1_helpers.compute_delta_d(d, dt, saveData, t_on, preIdx, postIdx);
tOnProc = out.tOnProc;
preIdxProc = out.preIdxProc;
postIdxProc = out.postIdxProc;
relTime = out.relTime;
baseline = out.baseline;
delta_d = out.delta_d;
idxEnd = postIdxProc(end);

outDir = fullfile(pwd, 'results', 'experiment1_single');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

figure(10); clf;
plot(relTime, delta_d, 'LineWidth', 1.5);
xlabel('Seconds from heat onset');
ylabel('\Delta d (baseline-normalized)');
title(sprintf('Recovery (worm %s)', worm.uid));
saveas(gcf, fullfile(outDir, 'recovery.png'));

% --- E1 metrics for this worm ---
t_off = t_on; % impulse (single index)
t_star = min(numel(d), tOnProc + max(1, round(T_STAR_SEC / dt)));

% Phase at stimulus (t_on^-), delay-embedded timebase.
t_on_prev = max(1, tOnProc - 1);
alpha0 = alpha(t_on_prev);
theta0 = theta(t_on_prev);
nBins = max(theta);
theta0_norm = theta0 / nBins;

% Pre-phase velocity (bins per frame, circular).
preWindow = preIdxProc;
v_pre = mean(experiment1_helpers.circ_diff_bins(theta(preWindow), nBins), 'omitnan');

% PRC-like phase shift.
theta_star = theta(t_star);
pred = theta0 + v_pre * (t_star - tOnProc);
dtheta_bins = experiment1_helpers.wrap_diff_bins(theta_star - pred, nBins);
dtheta_norm = dtheta_bins / nBins;

% Recovery metrics.
postWindow = tOnProc:min(numel(d), idxEnd);
delta_post = d(postWindow) - baseline;
d_peak = max(delta_post);
t_half = nan;
if ~isnan(d_peak) && d_peak > 0
    idx_half = find(delta_post <= 0.5 * d_peak, 1, 'first');
    if ~isempty(idx_half)
        t_half = idx_half * dt;
    end
end


% Save per-worm event row.
tabPath = fullfile(outDir, 'events.csv');
headers = {'uid','t_on','t_off','dt_sec','nBins','alpha0','theta0','theta0_norm', ...
           'd_peak','t_half_sec','dtheta_bins','dtheta_norm'};
row = {worm.uid, t_on, t_off, dt, nBins, alpha0, theta0, theta0_norm, ...
       d_peak, t_half, dtheta_bins, dtheta_norm};
T = cell2table(row, 'VariableNames', headers);
writetable(T, tabPath);

% Summary CSV (single-row, all-baseline-compatible fields).
summaryPath = fullfile(outDir, 'summary.csv');
summary = experiment1_helpers.summary_metrics( ...
    alpha, theta, d, relTime, delta_d, worm, t_on, dt, nBins);
Tsum = struct2table(summary);
writetable(Tsum, summaryPath);

% --- Diagnostics (LOOPER scripts) ---
diagDir = fullfile(outDir, 'diagnostics');
diagTag = sprintf('atanas_heat_%s', worm.uid);
run_looper_diagnostics(saveData, struct('tag', diagTag, 'outDir', diagDir, ...
    'dt_sec', dt, 'preIdxRaw', preIdx, 'postIdxRaw', postIdx, 'tOnRaw', t_on));

figPath = fullfile(outDir, 'final_stream_pca.fig');
if exist(figPath, 'file')
    openfig(figPath, 'new', 'visible');
end

fprintf('Diagnostics saved in: %s\n', diagDir);
