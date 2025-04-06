%% ERP ANALYSIS WITH 60HZ NOTCH + BANDPASS FILTERING + POLARITY CORRECTION
clear; clc;

%% === 1. Load OpenBCI EEG File ===
filename = 'OpenBCI-RAW-2025-04-04_13-12-14.txt';

opts = detectImportOptions(filename, 'FileType', 'text');
opts.VariableNamesLine = 5;
opts.DataLine = 6;
EEGtable = readtable(filename, opts);

Fs = 250;                       % Sampling rate
preStim = 0.2; postStim = 0.8;
samplesPre = round(preStim * Fs);
samplesPost = round(postStim * Fs);
epochLength = samplesPre + samplesPost;
timeVec = (-samplesPre:samplesPost-1)/Fs;

%% === 2. Extract EEG & Photodiode Signal ===
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

%% === 3a. Apply 60 Hz Notch Filter ===
notchFreq = 60;    
bandwidth = 1;     % Narrow bandwidth

wo = notchFreq / (Fs/2);   % Normalize
bw = bandwidth / (Fs/2);   
[bn, an] = iirnotch(wo, bw);

eegNotched = filtfilt(bn, an, double(eegData));

%% === 3b. Apply Bandpass Filter (0.5–30 Hz) ===
lowCut = 0.5; highCut = 30;
[b, a] = butter(4, [lowCut highCut]/(Fs/2), 'bandpass');
filteredEEG = filtfilt(b, a, eegNotched);

%% === 4. Invert Polarity for OpenBCI EEG ===
correctedEEG = -filteredEEG;

%% === 5. Detect Stimulus Onsets with Photodiode Thresholding ===
threshold = 700;
belowThresh = analog0 < threshold;
stimOnsets = find(diff(belowThresh) == 1) + 1;

shortPulseSamples = 32;  % ~128 ms
longPulseSamples = 128;  % ~512 ms

targetEpochs = [];
nonTargetEpochs = [];

for idx = 1:length(stimOnsets)
    onset = stimOnsets(idx);

    % Pulse duration
    pulseSamples = 0;
    while (onset + pulseSamples <= length(analog0)) && analog0(onset + pulseSamples) < threshold
        pulseSamples = pulseSamples + 1;
    end

    % Extract epoch
    epochStart = onset - samplesPre;
    epochEnd = onset + samplesPost - 1;

    if epochStart > 0 && epochEnd <= length(correctedEEG)
        epochEEG = correctedEEG(epochStart:epochEnd, :);

        % Classify by duration
        if abs(pulseSamples - shortPulseSamples) <= 5
            nonTargetEpochs = cat(3, nonTargetEpochs, epochEEG);
        elseif abs(pulseSamples - longPulseSamples) <= 10
            targetEpochs = cat(3, targetEpochs, epochEEG);
        end
    end
end

fprintf('Extracted %d Target epochs, %d Non-Target epochs\n', ...
    size(targetEpochs,3), size(nonTargetEpochs,3));

%% === 6. Compute Average ERP ===
meanTargetERP = mean(targetEpochs, 3);
meanNonTargetERP = mean(nonTargetEpochs, 3);

%% === 7. Plot ERPs by Channel ===
for ch = 1:size(correctedEEG,2)
    figure('Color','white');
    
    plot(timeVec, meanNonTargetERP(:,ch), 'Color', [0.9290, 0.6940, 0.1250], 'LineWidth', 1.5); hold on;
    plot(timeVec, meanTargetERP(:,ch), 'Color', [0, 0.4470, 0.7410], 'LineWidth', 1.5);

    % Onset / baseline lines
    yLimits = ylim;
    line([0 0], yLimits, 'Color','k','LineStyle','--','LineWidth',1.2);
    line([-0.2 -0.2], yLimits, 'Color','k','LineStyle',':','LineWidth',1);
    ylim(yLimits); xlim([-0.2, 0.8]);

    title(['ERP – ', channelLabels{ch}], 'FontSize', 14);
    xlabel('Time (s)'); ylabel('Amplitude (\muV)');
    legend('Non-Target','Target','Location','best');
    set(gca, 'FontSize', 12, 'LineWidth', 1.2);
    grid off; box off;
end

%% === 8. Save Output Data ===
save('ERP_results_with_notch.mat', ...
     'meanTargetERP', 'meanNonTargetERP', 'timeVec', 'channelLabels', ...
     'targetEpochs', 'nonTargetEpochs');
