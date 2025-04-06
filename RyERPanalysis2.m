%% ERP ANALYSIS – Modern Notch Filter, Bandpass, Polarity Fix, Trial Logging
clear; clc;

%% === 1. Load EEG File ===
filename = 'OpenBCI-RAW-2025-04-04_13-12-14.txt';

opts = detectImportOptions(filename, 'FileType', 'text');
opts.VariableNamesLine = 5;
opts.DataLine = 6;
EEGtable = readtable(filename, opts);

Fs = 250; % sampling rate
preStim = 0.2; postStim = 0.8;
samplesPre = round(preStim * Fs);
samplesPost = round(postStim * Fs);
timeVec = (-samplesPre:samplesPost-1)/Fs;

%% === 2. Extract EEG & Photodiode ===
analog0 = EEGtable.AnalogChannel0;
eegData = [EEGtable.EXGChannel0, EEGtable.EXGChannel1, EEGtable.EXGChannel2, ...
           EEGtable.EXGChannel3, EEGtable.EXGChannel4, EEGtable.EXGChannel5, ...
           EEGtable.EXGChannel6, EEGtable.EXGChannel7];

channelLabels = {
    'Anterior Frontal midline (AFz)', 
    'Frontal midline (Fz)', 
    'Central midline (Cz)', 
    'Centroparietal right (CP2)', 
    'Parietal midline (Pz)', 
    'Parietal left (P3)', 
    'Parietal right (P4)', 
    'Occipital left (O1)'};

%% === 3a. 60 Hz Notch Filter – Toolbox-Free ===
notchFreq = 60; 
bandwidth = 1;

lowStop = (notchFreq - bandwidth) / (Fs/2);
highStop = (notchFreq + bandwidth) / (Fs/2);

[bn, an] = butter(2, [lowStop highStop], 'stop');
eegNotched = filtfilt(bn, an, double(eegData));

%% === 3b. Bandpass Filter (0.5–30 Hz) ===
lowCut = 0.5; highCut = 30;
[b, a] = butter(4, [lowCut highCut]/(Fs/2), 'bandpass');
filteredEEG = filtfilt(b, a, eegNotched);

%% === 4. Polarity Correction for OpenBCI ===
correctedEEG = -filteredEEG;

%% === 5. Detect Stimulus Onsets from Photodiode Drops ===
threshold = 700;
belowThresh = analog0 < threshold;
stimOnsets = find(diff(belowThresh) == 1) + 1;

shortPulseSamples = 32;  % ~128ms
longPulseSamples  = 128; % ~512ms
shortTol = 5; longTol = 10;

targetEpochs = [];
nonTargetEpochs = [];

% Logging trial inclusion/exclusion
includedIndices = [];
excludedIndices = [];
excludedReasons = {};

for idx = 1:length(stimOnsets)
    onset = stimOnsets(idx);
    pulseSamples = 0;
    
    while (onset + pulseSamples <= length(analog0)) && (analog0(onset + pulseSamples) < threshold)
        pulseSamples = pulseSamples + 1;
    end

    epochStart = onset - samplesPre;
    epochEnd = onset + samplesPost - 1;

    if epochStart > 0 && epochEnd <= length(correctedEEG)
        epochEEG = correctedEEG(epochStart:epochEnd, :);

        if abs(pulseSamples - shortPulseSamples) <= shortTol
            nonTargetEpochs = cat(3, nonTargetEpochs, epochEEG);
            includedIndices(end+1) = idx;
        elseif abs(pulseSamples - longPulseSamples) <= longTol
            targetEpochs = cat(3, targetEpochs, epochEEG);
            includedIndices(end+1) = idx;
        else
            excludedIndices(end+1) = idx;
            excludedReasons{end+1} = sprintf('Unrecognized pulse duration: %d samples', pulseSamples);
        end
    else
        excludedIndices(end+1) = idx;
        excludedReasons{end+1} = 'Out-of-bounds epoch';
    end
end

%% === 6. Print Trial Summary and Exclusions ===
fprintf('\n--- Trial Summary ---\n');
fprintf('Total detected trials: %d\n', length(stimOnsets));
fprintf('Included trials: %d\n', length(includedIndices));
fprintf('Excluded trials: %d\n', length(excludedIndices));
fprintf('\nExcluded trial indices: ');
disp(excludedIndices);

for i = 1:length(excludedIndices)
    fprintf('Trial %d excluded – %s\n', excludedIndices(i), excludedReasons{i});
end

%% === 7. Compute ERP Averages ===
meanTargetERP = mean(targetEpochs, 3);
meanNonTargetERP = mean(nonTargetEpochs, 3);

%% === 8. Plot ERP per Channel ===
for ch = 1:size(correctedEEG,2)
    figure('Color','white');
    
    plot(timeVec, meanNonTargetERP(:,ch), 'Color', [0.929, 0.694, 0.125], 'LineWidth', 1.5); hold on;
    plot(timeVec, meanTargetERP(:,ch), 'Color', [0, 0.447, 0.741], 'LineWidth', 1.5);
    
    yLimits = ylim;
    xline(0, '--k', 'LineWidth', 1.2);
    xline(-0.2, ':k', 'LineWidth', 1);
    ylim(yLimits); xlim([-0.2, 0.8]);

    title(['ERP – ', channelLabels{ch}], 'FontSize', 14);
    xlabel('Time (s)'); ylabel('Amplitude (\muV)');
    legend('Non-Target','Target','Location','best');
    set(gca, 'FontSize', 12, 'LineWidth', 1.2);
    grid on;
end

%% === 9. Save Output ===
save('ERP_results_clean.mat', ...
     'meanTargetERP', 'meanNonTargetERP', 'timeVec', ...
     'channelLabels', 'targetEpochs', 'nonTargetEpochs', ...
     'includedIndices', 'excludedIndices', 'excludedReasons');
