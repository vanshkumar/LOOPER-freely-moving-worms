% folder pointing to where the data are located
maindir = 'F:/MonkeyData/Working Memory/data/';

% output file name
outputFileName = 'data_romo_eLife.mat';

DISPLAY = 0;

% This script was written to preprocess several different datasets and is
% general enough to handle slightly varying data formats. Only PFC data 
% from monkey 14 and monkey 15 are part of the eLife publication. The
% following parameters set it up. Other sets of parameters are commented
% out and are not relevant.

dirs = {'rr014/prefront', 'rr015/prefront/left', 'rr015/prefront/right'};
monkeys = [1 2 2];     % monkey mask for each folder
areas   = [1 1 1];     % area mask for each folder (always PFC in this case)
switchAreaAt = NaN;    % only applies to monkeys 24-25 (where recordings from two areas are stored in the same data-file)
% f1stimulationset = [10 14 18 24 30 34];  % set of F1 stimulation frequencies
f1stimulationset = [10 24 34];  % set of F1 stimulation frequencies
f1onlyLower =  [];     % f1 which were always lower than f2 (none for monkeys 14-15)
f1onlyHigher = [];     % f1 which were always higher than f2 (none for monkeys 14-15)
electrodeNum = 7;      % number of electrodes
timeInt = [-2500 10000]; % interval around F1 onset to be cut as one trial (in ms)

% dirs = {'rr014/prefront', 'rr014/s2', 'rr015/prefront/left', 'rr015/prefront/right'};
% monkeys = [1 1 2 2];    % monkey mask for each folder
% areas   = [1 2 3 4];    % area mask for each folder
% switchAreaAt = NaN;     % only applies to monkeys 24-25 (where recordings from two areas are stored in the same data-file)
% f1stimulationset = [10 14 18 24 30 34]; % set of F1 stimulation frequencies
% f1onlyLower =  [];      % f1 which were always lower than f2 (none for monkeys 14-15)
% f1onlyHigher = [];      % f1 which were always higher than f2 (none for monkeys 14-15)
% electrodeNum = 7;       % number of electrodes
% timeInt = [-500 4500]; % interval around F1 onset to be cut as one trial (in ms)

% dirs = {'rr015/prefront/left', 'rr015/prefront/right'};
% monkeys = [1 1];    % monkey mask for each folder
% areas   = [1 2];    % area mask for each folder
% switchAreaAt = NaN;     % only applies to monkeys 24-25 (where recordings from two areas are stored in the same data-file)
% f1stimulationset = [10 14 18 24 30 34]; % set of F1 stimulation frequencies
% f1onlyLower =  [];      % f1 which were always lower than f2 (none for monkeys 14-15)
% f1onlyHigher = [];      % f1 which were always higher than f2 (none for monkeys 14-15)
% electrodeNum = 7;       % number of electrodes
% timeInt = [-500 7500]; % interval around F1 onset to be cut as one trial (in ms)

% dirs = {'rr013/m1', 'rr013/prefront/left', 'rr013/prefront/right', 'rr013/s2cont', 'rr013/sma'};
% monkeys = [1 1 1 1 1];
% areas   = [1 2 3 4 5];
% switchAreaAt = NaN;
% f1stimulationset = [10 14 18 22 26 30 34];
% f1onlyLower =  [1 2]; 
% f1onlyHigher = [6 7]; 
% electrodeNum = 7;     
% timeInt = [-500 4500]; % interval around F1 onset to be cut as one trial (in ms)

% dirs = {'R24/', 'R25/'};
% monkeys = [1 2];
% areas   = [1 2]; % NB: units 1-42 are from S1 and 43-84 from S2
% switchAreaAt = 42;
% f1stimulationset = [10 14 16 18 20 22 24 26 28 30 34];
% f1onlyLower =  [1 2 3 5];
% f1onlyHigher = [7 9 10 11];
% electrodeNum = 84;
% timeInt = [-500 4500]; % interval around F1 onset to be cut as one trial (in ms)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f1Both = setdiff(1:length(f1stimulationset), [f1onlyLower f1onlyHigher]);

unitCounter = 0;
rejectedUnits = 0;
usedUnits = 0;
rejectedSessions = 0;
okSessions = 0;
unstableNeurons = 0;

x = -1000:1000;
window = 200;
gaussKernel = 1/sqrt(2*pi)/window * exp(-x.^2/window^2/2);

time = timeInt(1)/1000 : 0.01 : timeInt(2)/1000;

% Initializing arrays with sizes larger than necessary. In the end the 
% script will reduce the sizes to the ones actually used. If not 
% initialized large enough, the script will work but will be slower.
firingRates = zeros(1500, 6, 2, length(time), 50);
firingRatesAverage = zeros(1500, 6, 2, length(time));
trialNum = [];

firingRatesMisses = zeros(1500, 6, 2, length(time), 50);
firingRatesAverageMisses = zeros(1500, 6, 2, length(time));
trialNumMisses = [];

spikeRasters = zeros(1500, 6, 2, length(time), 50);

monkeyMask = [];
areaMask = [];

totalMistakes = [];

minimumTrialCount = 10;
minimumMistakes = 0;
useMisses = 0;

% Loop over all folders
for d = 1:length(dirs)
    if DISPLAY
        display(['Processing folder ' dirs{d} ' (' num2str(d) ' out of ' num2str(length(dirs)) ')'])
    end
    
    files = what([maindir dirs{d}]);
    
    % Loop over all files in the folder
    for f = 1:length(files.mat)
        if DISPLAY
            fprintf(['Processing session ' num2str(f) ' out of ' num2str(length(files.mat))])
        end
        
        load ([maindir dirs{d} '/' files.mat{f}])
        
        freqs = cell2mat(result(2:end, 4));
        [freqset, ~, freqs] = unique(freqs);
        freqsToUse = find(ismember(freqset, f1stimulationset));
        
        % Checking that this session contains all frequencies
        if length(freqsToUse) ~= length(f1stimulationset)
            if DISPLAY
                fprintf([' - skipping, not all frequencies present: ' num2str(freqset') '\n'])
            end
            rejectedSessions = rejectedSessions + 1;
            continue
        end
        
        % sometimes some info is missing - omit these trials
        for tt = 2:size(result,1)
            if isempty(result{tt,9})
                omit(tt-1) = 1;
                if DISPLAY
                    fprintf(['omitting some trials in ' files.mat{f}])
                end
            else
                omit(tt-1) = 0;
            end
        end
        
        decisions = double(cell2mat(result(2:end, 4)) > cell2mat(result(2:end,5)));
        mistakes = ~cell2mat(result(2:end, 3));
        decisions(mistakes) = ~decisions(mistakes);
                
        decisions(omit==1) = nan;
        
        % Checking that this session contains all decisions for those
        % frequencies for which f2>f1 and f2<f1 were both present
        freqsToUseBothDec = find(ismember(freqset, f1stimulationset(f1Both)));
        subset = find(ismember(freqs,freqsToUseBothDec) & mistakes==0);
        if (size(unique([freqs(subset) decisions(subset)], 'rows'), 1) ~= length(f1Both)*2)
            if DISPLAY
                fprintf(' - skipping, not enough decisions for "main" frequencies\n')
            end
            rejectedSessions = rejectedSessions + 1;
            continue
        end
        
        % Checking that this session contains expected decisions for the
        % frequencies for which always f2<f1 or always f2<f1
        ifReject = false;
        for ii = f1onlyLower
            if(isempty(find(freqs==freqsToUse(ii) & decisions==0 & mistakes==0)))
                ifReject = true;
                break
            end
        end
        for ii = f1onlyHigher
            if(isempty(find(freqs==freqsToUse(ii) & decisions==1 & mistakes==0)))
                ifReject = true;
                break
            end
        end        
        if ifReject
            if DISPLAY
                fprintf([' - skipping: no right decision for frequency ' num2str(f1stimulationset(ii)) '\n'])
            end
            rejectedSessions = rejectedSessions + 1;
            continue
        end
        
        
        for ii = 1:length(freqsToUse)
            if(sum((freqs==freqsToUse(ii) & decisions==1 & mistakes==0)) < minimumTrialCount)
                ifReject = true;
                break
            end
        end
        for ii = 1:length(freqsToUse)
            if(sum((freqs==freqsToUse(ii) & decisions==0 & mistakes==0)) < minimumTrialCount)
                ifReject = true;
                break
            end
        end
        if ifReject
            if DISPLAY
                fprintf([' - skipping: not enough trials for frequency ' num2str(f1stimulationset(ii)) '\n'])
            end
            rejectedSessions = rejectedSessions + 1;
            continue
        end
        
        for ii = 1:length(freqsToUse)
            if(sum((freqs==freqsToUse(ii) & decisions==1 & mistakes==1)) < minimumMistakes)
                ifReject = true;
                break
            end
        end
        for ii = 1:length(freqsToUse)
            if(sum((freqs==freqsToUse(ii) & decisions==0 & mistakes==1)) < minimumMistakes)
                ifReject = true;
                break
            end
        end
        if ifReject
            if DISPLAY
                fprintf([' - skipping: not enough misses for frequency ' num2str(f1stimulationset(ii)) '\n'])
            end
            rejectedSessions = rejectedSessions + 1;
            continue
        end
        
        % If this point is reached, the session is fine
        okSessions = okSessions + 1;
        
        % Counting mistakes
        totalMistakes = [totalMistakes; mistakes];
        %continue
        
        % Selecting all units that have at least 1 spike in at least 1 trial
        % Nota bene: electrode #9 is recording stimulus! (for monkeys 13,14,15)
        unitsToUse = nan(1,electrodeNum);       %length(result{2,6})
        for tr=2:size(result,1)
            for u=1:electrodeNum 
                if ~isempty(result{tr,6}{u})
                    unitsToUse(u) = u;
                end
            end
        end
        
        rejectedUnits = rejectedUnits + length(isnan(unitsToUse));
        unitsToUse = find(~isnan(unitsToUse));
        
        % Loop over all responding units
        for unit = unitsToUse
            if DISPLAY
                fprintf('.')
            end
            
            unitCounter = unitCounter + 1;
                        
            % Filling the masks
            monkeyMask(unitCounter) = monkeys(d);
            if isnan(switchAreaAt)
                areaMask(unitCounter) = areas(d);
            else
                if unit > switchAreaAt
                    areaMask(unitCounter) = areas(2);
                else
                    areaMask(unitCounter) = areas(1);
                end
            end
            
            % Loop over all frequencies
            for fr = 1:length(freqsToUse)*2
                
                % Loop over all decisions
                for dec = 0:1
                    if fr > length(freqsToUse)
                        miss = 1;
                    else
                        miss = 0;
                    end
                    
                    realFR = mod(fr - 1, length(freqsToUse)) + 1;
                    
                    % For impossible/very rare frequency-decision
                    % combinations we fill the data with NaNs
                    if ~isempty(f1onlyLower) || ~isempty(f1onlyHigher)
                        if(ismember(realFR, f1onlyLower) && dec==1 || ismember(realFR, f1onlyHigher) && dec==0)
                            firingRatesAverage(unitCounter, realFR, dec+1, :) = nan(size(firingRatesAverageMisses,4),1);
                            trialNum(unitCounter, realFR, dec+1) = 0;
                            continue
                        end
                    end

                    % For possible combinations -- Loop over trials
                    trialsToUse = find(freqs == freqsToUse(realFR) & decisions == dec & mistakes == miss);
                    timepointCount = (timeInt(2)-timeInt(1))/10+1;

                    repetitions = zeros(length(trialsToUse), timepointCount);
                    rawSpikes = zeros(length(trialsToUse), timepointCount);
                    for tr = 1:length(trialsToUse)
                        sp = result{1+trialsToUse(tr),6};
                        so1 = round(result{1+trialsToUse(tr),9});

                        if so1 < -timeInt(1)
                            continue
                        end

                        spiketrain = round(sp{unit});
                        spiketrainFull = zeros(1,30000);
                        spiketrainFull(round(spiketrain(spiketrain>0))) = 1;
                        trialrate = conv(spiketrainFull, gaussKernel, 'same');
                        chunk = trialrate(so1 + timeInt(1) : so1 + timeInt(2));
                        repetitions(tr, :) = chunk(:, 1:10:end) * 1000;

                        trialrate = spiketrainFull;
                        chunk = trialrate(so1 + timeInt(1) : so1 + timeInt(2));
                        rawSpikes(tr, :) = chunk(:, 1:10:end) * 1000;
                    end

                    if isempty(trialsToUse)
                        trialNum(unitCounter, fr, dec+1) = 0;
                        continue
                    end         

                    % Averaging over trials
                    firingRatesAverage(unitCounter, fr, dec+1, 1:size(repetitions,2)) = mean(repetitions, 1);
                    trialNum(unitCounter, fr, dec+1) = size(repetitions, 1);

                    % Filling in single trial array
                    firingRates(unitCounter, fr, dec+1, 1:size(repetitions,2), 1:size(repetitions,1)) = repetitions';

                    spikeRasters(unitCounter, fr, dec+1, 1:size(rawSpikes,2), 1:size(rawSpikes,1)) = rawSpikes';
                        
                end
            end
        end
        if DISPLAY
            fprintf('\n')
        end
    end
end

useFrequencies = length(f1stimulationset);
if useMisses
    useFrequencies = useFrequencies * 2;
end

% drop nonused zeros
firingRates = firingRates(1:unitCounter, 1:useFrequencies, ...
                            1:2, 1:timepointCount, 1:max(trialNum(:)));
spikeRasters = spikeRasters(1:unitCounter, 1:useFrequencies, ...
                            1:2, 1:timepointCount, 1:max(trialNum(:)));
firingRatesAverage = firingRatesAverage(1:unitCounter, 1:useFrequencies, ...
                            1:2, 1:timepointCount);
trialNum = trialNum(1:unitCounter, 1:useFrequencies, :);
                        
firingRatesMisses = firingRatesMisses(1:unitCounter, 1:length(f1stimulationset), ...
                            1:2, 1:timepointCount, 1:max(trialNumMisses(:)));
firingRatesAverageMisses = firingRatesAverageMisses(1:unitCounter, 1:length(f1stimulationset), ...
                            1:2, 1:timepointCount);

if DISPLAY
    display(['Number of rejected sessions: ' num2str(rejectedSessions)])
    display(['Number of processed sessions: ' num2str(okSessions)])
    display(['Number of rejected (empty) units: ' num2str(rejectedUnits)])
    display(['Number of processed units: ' num2str(unitCounter)])
end
                        
% sparsifying 
firingRates_size = size(firingRates);
firingRates_sparse = sparse(firingRates(:));
firingRatesMisses_size = size(firingRates);
firingRatesMisses_sparse = sparse(firingRates(:));
timeEvents = [0 0.5 3.5 4.0 7.0];
timeEventsNames = {'F1start', 'F1stop', 'F2start', 'F2stop', 'Cue'};

% display('Saving...')
% save([maindir outputFileName], 'firingRates_sparse', 'firingRates_size', ...
%                      'firingRatesAverage', 'trialNum', ...
%                      'time', 'timeEvents', 'timeEventsNames', ...
%                      'monkeyMask', 'areaMask', 'firingRatesMisses_sparse', 'firingRatesAverageMisses')
%                  
display('Done')