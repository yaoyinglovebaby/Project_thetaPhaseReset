%% Set path
addpath '\\analyse4.psy.gla.ac.uk\project0309\Ying_Phasereset\analyses\prediction\simulation'
addpath '\\analyse4.psy.gla.ac.uk\project0309\Ying_Phasereset\analyses\functions'

%% Simulate data -- Phase reset %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:100
    cfg_all(i).centerFreq = 5 + randi([-2, 30]);
    cfg_all(i).freqVariation = 0.5+ 0.1 * randn();
    cfg_all(i).freqDriftWindow = 200;
    cfg_all(i).spikeDensitySmoothing = 100;
    cfg_all(i).epochDuration = [-4 4];
    cfg_all(i).numTrials = 100;
    cfg_all(i).eventTime = 0.2 + 0.01 * randi([0, 10]);  % 加一点随机抖动
    cfg_all(i).eventJitter = 0.1+ 0.01 * randi([0, 10]);
    cfg_all(i).spikePhase = pi * (randi([0 1]));  % 0 或 pi    
    cfg_all(i).noise = 0.6 + 0.01 * randi([0, 10]);       % baseline noise 随机波动
    cfg_all(i).firingRate = 100+ 1 * randi([0, 10]);%  cfg_all(i).centerFreq;%
    cfg_all(i).responseDuration =round(1 + 0.1 * randi([0, 10]))*0.1; %  0;%
    cfg_all(i).plotFigures = 0;  
end
% 生成随机type标签
types = [repmat({'pr'}, 1, 50), repmat({'ERP'}, 1, 50)];
% 随机打乱顺序
types = types(randperm(100));

classify_by_prediction_with_roc(cfg_all,1,types)
%%
% Simulate phase reset data
[spikeTimes_phr, timeVec, ~, ~] = simulate_phase_reset(centerFreq, freqVariation, freqDriftWindow, ...
    spikeDensitySmoothing, epochDuration, numTrials, eventTime, eventJitter, spikePhase, plotFigures, noise);

% Simulate ERP data
[spikeTimes_erp, timeVec, ~, ~] = simulate_evoked_resp(centerFreq, freqVariation, freqDriftWindow, ...
    spikeDensitySmoothing, epochDuration, numTrials, eventTime, eventJitter, responseDuration, firingRate, spikePhase, plotFigures);

classify_and_evaluate_with_roc(spikeTimes_erp, spikeTimes_phr, timeVec)
%%
% 画ROC
% classify_by_prediction_with_roc(spikeTimes_erp, spikeTimes_phr, timeVec)

