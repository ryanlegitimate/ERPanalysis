%% ======================================================================
%  OPENBCI ERP ANALYSIS
%  • Imports the raw *.txt* exactly the way MATLAB’s Import Tool does
%  • Extracts the D17 digital trigger (LOW = flash)
%  • Removes micro‑glitches (< 30 samples)
%  • Uses the two dominant pulse widths you observed:
%        – SHORT  ≈ 45 samples  (~176 ms)  →  non‑target   (80 %)
%        – LONG   ≈ 135 samples (~544 ms)  →  target       (20 %)
%  • Baseline‑corrects each epoch (‑200 ms..0 ms)
%  • Prints trial counts, plots ERPs, saves a .mat file
% ======================================================================

clear; clc

%% 1. Import  -----------------------------------------------------------
filename = 'C:\Users\ryleg\Desktop\OpenBCI-RAW-2025-06-11_10-37-35.txt';

opts = delimitedTextImportOptions("NumVariables",25);
opts.DataLines = [6 Inf];           % skip the 5 header lines
opts.Delimiter = ",";

opts.VariableNames = [ ...
    "SampleIndex","EXGChannel0","EXGChannel1","EXGChannel2","EXGChannel3", ...
    "EXGChannel4","EXGChannel5","EXGChannel6","EXGChannel7", ...
    "AccelChannel0","AccelChannel1","AccelChannel2","NotUsed", ...
    "DigitalChannel0D11","DigitalChannel1D12","DigitalChannel2D13", ...
    "DigitalChannel3D17","NotUsed1","DigitalChannel4D18", ...
    "AnalogChannel0","AnalogChannel1","AnalogChannel2", ...
    "Timestamp","MarkerChannel","TimestampFormatted" ];

opts.VariableTypes           = repmat({'double'},1,25);
opts.VariableTypes{25}       = 'datetime';
opts.ExtraColumnsRule        = "ignore";
opts.EmptyLineRule           = "read";
opts = setvaropts(opts,"TimestampFormatted", ...
                  "InputFormat","yyyy-MM-dd HH:mm:ss.SSS", ...
                  "DatetimeFormat","preserveinput");

EEGtable = readtable(filename, opts);
clear opts filename

%% 2. Constants  --------------------------------------------------------
Fs            = 250;        % Hz
preStim       = 0.2;        % s  (baseline window)
postStim      = 0.8;        % s
samplesPre    = round(preStim  * Fs);
samplesPost   = round(postStim * Fs);
timeVec       = (-samplesPre:samplesPost-1)/Fs;

% Pulse definitions (derived from pulse‑width overview)
shortPulseSamples = 45;     shortTol = 10;        
longPulseSamples  = 135;    longTol  = 10;        

channelLabels = {'AFz','Fz','Cz','CP2','Pz','P3','P4','O1'};

%% 3. EEG extraction & cleaning  ----------------------------------------
eeg = [EEGtable.EXGChannel0, EEGtable.EXGChannel1, EEGtable.EXGChannel2, ...
       EEGtable.EXGChannel3, EEGtable.EXGChannel4, EEGtable.EXGChannel5, ...
       EEGtable.EXGChannel6, EEGtable.EXGChannel7];

% 60 Hz notch
[bn,an] = butter(2,[(60-1) (60+1)]/(Fs/2),'stop');
eeg = filtfilt(bn,an,double(eeg));

% 0.5–30 Hz band‑pass
[bb,aa] = butter(4,[0.5 30]/(Fs/2),'bandpass');
eeg = filtfilt(bb,aa,eeg);

% Polarity fix
eeg = -eeg;

%% 4. Trigger from D17  --------------------------------------------------
byteD  = uint8(EEGtable.DigitalChannel3D17);   % 0‑255 byte value
pulse  = bitget(byteD,8)==0;                   % LOW (0) during flash

starts = find(diff([0; pulse])==1);
stops  = find(diff([pulse;0])==-1);
nPulses = min(numel(starts),numel(stops));
starts  = starts(1:nPulses);   stops = stops(1:nPulses);
pulseDur = stops - starts;                     % width in samples

% --- remove micro‑pulses (<30 samples) ------------------------------
keep = pulseDur >= 30;
starts   = starts(keep);
pulseDur = pulseDur(keep);
stimOnsets = starts;

%% 5. Epoch extraction, baseline correction, classification --------------
targetEpochs    = [];
nonTargetEpochs = [];
incIdx = [];  excIdx = [];  excWhy = {};

for k = 1:numel(stimOnsets)
    onset = stimOnsets(k);
    width = pulseDur(k);
    epochStart = onset - samplesPre;
    epochEnd   = onset + samplesPost - 1;

    if epochStart>0 && epochEnd<=size(eeg,1)
        ep = eeg(epochStart:epochEnd,:);
        ep = ep - mean(ep(1:samplesPre,:),1);     % baseline

        if abs(width-shortPulseSamples)<=shortTol
            nonTargetEpochs = cat(3,nonTargetEpochs,ep);  incIdx(end+1)=k;
        elseif abs(width-longPulseSamples)<=longTol
            targetEpochs    = cat(3,targetEpochs,ep);     incIdx(end+1)=k;
        else
            excIdx(end+1)=k;  excWhy{end+1}=sprintf('Width %d',width);
        end
    else
        excIdx(end+1)=k;  excWhy{end+1}='Out‑of‑bounds';
    end
end

%% 6. Print trial counts  ------------------------------------------------
fprintf('\n--- Trial Summary ---\n');
fprintf('Total pulses found : %d\n', numel(pulseDur)+numel(find(pulseDur<30)));
fprintf('Included trials    : %d\n', numel(incIdx));
fprintf('Target trials      : %d\n', size(targetEpochs,3));
fprintf('Non‑target trials  : %d\n\n', size(nonTargetEpochs,3));

%% 7. Grand averages & plots  -------------------------------------------
meanT  = mean(targetEpochs,3);
meanNT = mean(nonTargetEpochs,3);

if isempty(meanT)||isempty(meanNT)
    warning('No valid epochs — adjust pulse settings.');
else
    for ch = 1:8
        figure('Color','white');
        plot(timeVec,meanNT(:,ch),'LineWidth',1.4,'Color',[0.929 0.694 0.125]); hold on
        plot(timeVec,meanT(:,ch),'LineWidth',1.4,'Color',[0 0.447 0.741]);
        xline(0,'--k'); xlim([-0.2 0.8]);
        title(['ERP – ' channelLabels{ch}]);
        xlabel('Time (s)'); ylabel('Amplitude (\muV)');
        legend('Non‑Target','Target','Location','best'); grid on
    end
end

%% 8. Save ---------------------------------------------------------------
save('ERP_results_D17_baseline.mat', ...
     'meanT','meanNT','timeVec','channelLabels', ...
     'targetEpochs','nonTargetEpochs','incIdx','excIdx','excWhy');
