function DATA = PwaveDetection(Filebase,params)

%
% this function is to detect P-waves 
%
% created by Shuzo Sakata (shuzo.sakata@strath.ac.uk)
% for further information, please refer to https://www.biorxiv.org/content/10.1101/752683v1.abstract

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. defining threshold for pEEG (based on amp during the longest NREM epoch) 
% 2. detecting candidates for P-waves (based on 5-30 Hz bandpass-filtered pEEG) 
% 3. excluding artefacts (1. mean+3SD of EMG rms; 2. mean+3SD of T2P of P-waves^REM)
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% useage
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUTS
%   - Filebase
%   - params
%       .nChs -- the total number of channels
%       .pChs -- 2 elements vector (for pEEG)
%       .emgCh -- 1 or 2 elements vector (for EMG)
%       .ext -- extension ('.dat' or '.eeg')
%
% OUTPUT
%   - DATA
%       .PwavTime --- timing of P-waves (vector)
%       .PwavWav --- P-waves waveforms (n x time)
%


%% default parameters
Fs = 1000; % sampling rate
maxVolts = 5;
sampsPerVolt = double(intmax('int16')/maxVolts); % convert to get volts

Twin=10; % calculate RMS every 10 msec
[b,a]=butter(5,[5 30]/(Fs/2),'bandpass'); % filter param (5-30 Hz bandpass) %%IMPORTANT PARAMETERS!!

Epoch = 4; % analyze every 4 sec
Time = 50; % extracted time window (in msec)
skip = Time; % skip next 50 msec from Pwave peak 

% basic parameters
Ext = params.ext;
Pch1 = params.pChs(1);
Pch2 = params.pChs(2);
EMGch = params.emgCh;
nCh = params.nChs;

params.twin_rms = Twin;
params.epoch_win = Epoch;
params.pWavWin = Time;

% calculating of total time of recording
fileInfo = dir([Filebase, Ext]);
fileSize = fileInfo.bytes;
FinishSec=fileSize/2/nCh/Fs;
fprintf(['rec dur = ',num2str(FinishSec/60, '%.1f'),' min\n']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Set threshold %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1. choose the longest NREM sleep episode
load('.\sleepscoring_results.mat', 'Score');

n_NREM = find(Score==1); % find NREM wins
n_NREM_diff = diff(n_NREM); % 
n_trans = find(n_NREM_diff ~= 1); % find transitions
NREMoffidx = [n_NREM(n_trans); n_NREM(end)]; % NREM off
NREMonidx = [n_NREM(1); n_NREM(n_trans+1)]; % NREM on

NREMon = NREMonidx * Epoch - Epoch; % in sec from start
NREMoff = NREMoffidx * Epoch; % in sec
duration = NREMoff - NREMon; % NREM duration

[NREMmax,I] = max(duration); % choose the longest NREM
start = NREMon(I,1);
startsec = start + 2*Epoch; 
stopsec = startsec + 10*Epoch; %%%%% need to ensure that this period should be in NREM!!
fprintf(['max NREM duration: ', num2str(NREMmax), 'sec\n']);
if NREMmax < (stopsec - start)
    error('NREM duration is not enough to establish a threshold!! \n');
end

%% 2. Extract pEEG channels and calculate RMS of EEG every 10 msec
% subtract pEEG
MySignal = LoadDatSeg([Filebase,Ext],nCh,Fs*(startsec-1)+1, Fs*(stopsec-1))/sampsPerVolt;
pEEG = MySignal(Pch1,:) - MySignal(Pch2,:);
pEEG = filtfilt(b,a,pEEG); % bandpass filter

% Calculate RMS every Twin
nTwins = (stopsec-startsec)*Fs/Twin;
MyEEG = reshape(pEEG, Twin, nTwins);
RMS = rms(MyEEG);

%% 3. Set threshold
Tval = 5;
MeanValue = mean(RMS);
SD = std(RMS);
thresh = -MeanValue - Tval*SD; %%%%%%%%%% IMPORTANT!!!!!!

% display
figure(5); 
plot(pEEG,'k');hold on;
plot([1 length(pEEG)],thresh*ones(2,1),'r:');hold on;
axis([0 length(pEEG)+1 5*thresh -5*thresh]);
title('extracted pEEG during the longest NREM episode');

%% 4. Extract EMG channels and calculate RMS of EMG every 10 msec
% load EMG signals
MySignal = LoadDatSeg([Filebase,Ext],nCh,1, Fs*FinishSec)/sampsPerVolt;
if length(EMGch) == 1
    EMG = MySignal(EMGch,:);
else
    EMG = MySignal(EMGch(1),:) - MySignal(EMGch(2),:);
end

nTwins = floor(FinishSec*Fs/Twin);
EMG = EMG(1:Twin*nTwins);
MyEMG = reshape(EMG, Twin, nTwins);
RMS_EMG = rms(MyEMG);

%% 5. Set thereshold for EMG
Tval = 3;
MeanValueEMG = mean(RMS_EMG);
SDEMG = std(RMS_EMG);
threshEMG = MeanValueEMG + Tval*SDEMG;

% display
figure(6);
plot(RMS_EMG, 'k'); hold on;
plot([1 length(RMS_EMG)], threshEMG*ones(2,1),'r'); hold on;
xlim([0 length(RMS_EMG)]);
title('RMS of EMG');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% P-wave detection (initial screening)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. extract entire pEEG signals
% MySignal = LoadDatSeg([Filebase, Ext],nCh,1, Fs*FinishSec)/sampsPerVolt;
pEEG = MySignal(Pch1,:) - MySignal(Pch2, :);
pEEG = filtfilt(b,a,pEEG); % bandpass filter
clear MySignal

%% 2. detect candidates
PwavTime=[]; % timing vector for Pwave (make container)
PwavWav = []; % waveform matrix 
EMGpwr = []; % EMG power during P-wave generations
MyEndPoint = length(Score) * Epoch * Fs;

t = Time;
while t < MyEndPoint-Time*2
    Tmp = pEEG(t); % scan each point (a bit inefficient ...)
    if Tmp < thresh % below the threshold?
        % determine the peak points
        MyData = pEEG(1,t:t+Time);
        [~, I] = min(MyData); % negative peak
        MyP = t+I-1; % netative peak point (PwavTime)
        
        %% P-wave    
        MyWav = pEEG(MyP - Time:MyP + Time);
        
        %% EMG
        EMGpwr(end+1, 1) = rms(EMG(MyP - Time:MyP + Time));
        
        %% update
        PwavTime(end+1,1) = MyP;
        PwavWav(end+1, :) = MyWav;
        t = MyP + skip; 
    else % no P-wave
        t = t+1;
    end
end

%%%%%%%%%%%%
%% post-cleaning 1 -- based on EMG rms (< mean + 3SD)
Idx = find(EMGpwr >= threshEMG);
PwavTime(Idx) = [];
PwavWav(Idx,:) = [];
fprintf(['EMG thre...:', num2str(length(Idx)), ' events excluded\n']);

%% post-cleaning 2 -- based on P2T amp of P-waves during REM (<= mean + 3SD)
% labeling
PwavSec = ceil(PwavTime/Fs);
tmp = repmat(Score',Epoch,1);
StateVec = tmp(:);
PwavState = StateVec(PwavSec);

% P-waves during REM sleep
if ~isempty(find(PwavState == 2))
    REMpWav = PwavWav(PwavState == 2, :); % REM
    T2P = peak2peak(REMpWav,2);
    T2Pthresh = mean(T2P) + 3*std(T2P); %%% threshold
else % no REM episode
    NREMpWav = PwavWav(PwavState == 1, :); % NREM
    T2P = peak2peak(NREMpWav,2);
    T2Pthresh = mean(T2P) + 3*std(T2P); %%% threshold
end
AllT2P = peak2peak(PwavWav,2);
Idx = find(AllT2P >= T2Pthresh);
PwavTime(Idx) = [];
PwavWav(Idx, :) = [];
fprintf(['T2P thre...:', num2str(length(Idx)), ' events excluded\n']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% saving data
DATA.thresh = thresh;
DATA.threshEMG = threshEMG;
DATA.threshT2P = T2Pthresh;
DATA.PwavTime = PwavTime; % in msec
DATA.PwavWav = PwavWav;
DATA.EMGpwr = EMGpwr;

save(['PwaveInfo_',num2str(Pch1),'_',num2str(Pch2),'.mat'],'DATA', 'params');


