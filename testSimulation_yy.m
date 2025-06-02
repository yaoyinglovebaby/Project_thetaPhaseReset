% Generate two spiking neurons where one demonstrates a phase-reset and the
% other demonstrates an evoked response after stimulus onset
%%
datapath='Z:\Luca\data\microLFP\mySpkInt';

sessionName={'Encoding','Retrieval'};
lengendName={'Hit','Miss'};
titleText={'encoding hit','retrieval hit','encoding miss','retrieval miss'};       
spkmat_hit={'spkd_e.trial(hitsIdx'')','spkd_r.trial(hitsIdx'')'};
spkmat_miss={'spkd_e.trial(~hitsIdx'')','spkd_r.trial(~hitsIdx'')'};
spkmat={'spkd_e.trial(:)','spkd_r.trial(:)'};
colorValue={[0 0 1],[1 0 0]};
tvec=-4:0.001:6; %time
% load data
cd(datapath)
%%
n=262; selectedBand = 'Theta';  % 'Theta', 'Alpha', 'Beta', 'lowGamma', 'highGamma'
tmpspks = round(spks2plot(n).spks);
dt = 0:1:max(tmpspks) + 1; % 时间间隔
[spkdens, ~] = hist(tmpspks, dt);

pi=1;
ntrls_e=size(spks2plot(n).encTrigger, 1);%trial number
trls_e = [spks2plot(n).encTrigger(:, 1) + time(1), spks2plot(n).encTrigger(:, 1) + time(end), ones(ntrls_e, 1) .* -4, spks2plot(n).hitsIdx];
[spkd_e, countE] = continuous2trls(trls_e, data);
cfg = [];
cfg.channel = spks2plot(n).wirename; % 指定通道名，支持单个或多个通道
selected_data = ft_selectdata(cfg, spkd_e);
osci=vertcat(selected_data.trial{spks2plot(n).hitsIdx});

spks.trial{1, 1} = spkdens;
spks.time{1, 1} = dt(2:end-1) ./ 1000;
spks.label = {[spks2plot(n).wirename '-' num2str(allSpks(n).su)]};
spks.fsample = 1000;
[spks1, count] = continuous2trls(trls_e, spks);

pi=2;
ntrls_r = size(spks2plot(n).retTrigger, 1);
trls_r = [spks2plot(n).retTrigger(:, 1) + time(1), spks2plot(n).retTrigger(:, 1) + time(end), ones(ntrls_r, 1) .* -4, spks2plot(n).hitsIdx];
[spkd_r, countR] = continuous2trls(trls_r, data);
cfg = [];
cfg.channel = spks2plot(n).wirename; % 指定通道名，支持单个或多个通道
selected_data = ft_selectdata(cfg, spkd_r);
osci=vertcat(selected_data.trial{spks2plot(n).hitsIdx});

spks.trial{1, 1} = spkdens;
spks.time{1, 1} = dt(2:end-1) ./ 1000;
spks.label = {[spks2plot(n).wirename '-' num2str(allSpks(n).su)]};
spks.fsample = 1000;
[spks1, count] = continuous2trls(trls_r, spks);

%%
 figure;
    subplot(10,2,1:6);imagesc(tvec,[1:size(osci,1)],osci);
    title('ERP image');
    subplot(10,2,7:10);plot(tvec,mean(osci,1));
    title('ERP');
    subplot(10,2,11:16);
    title('Rasterplot');
    hold on;
    cntrl=0;spk_phi=[];
    spk_dens = zeros(size(osci,1),numel(tvec));
    for trial = 1:ntrls_e
        spike_times = spks1.trial{trial}; % Get spike times for the current trial        
        if spks2plot(n).hitsIdx(trial)==1 % only hit trials
            cntrl=cntrl+1;
            spk_dens(cntrl,:)=smoothdata(spike_times,'Gaussian', 250); %smooth by trial
            anlyt_osci=hilbert(osci(cntrl,:));
            phs_osci=angle(anlyt_osci);            
            spike_times=find(spike_times == 1);
            spk_phi=[spk_phi phs_osci(spike_times)];
            
        for spike = 1:length(spike_times)
            % Plot each spike as a vertical line segment
            line([tvec(spike_times(spike)), tvec(spike_times(spike))], [cntrl-0.5, cntrl+0.5], 'Color', 'k', 'LineWidth', 1);
        end
        end

    end
    hold off;

    % Format the plot
    xlabel('Time (s)');
    ylabel('Trial');
    title('Raster Plot');
    ylim([0.5, size(osci,1)+0.5]); % Adjust y-axis limits to fit trials
    yticks(1:10:size(osci,1)); % Set y-ticks for each trial
    xlim([min(tvec) max(tvec)])
    subplot(10,2,17:20);
    plot(tvec,mean(spk_dens,1));
    title('Spike Density');
    set(gcf, 'Units', 'inches', 'Position', [1, 1, 12, 16]);
    set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0, 0, 12, 16]);
    savepath= ['Z:\Ying_Phasereset\analyses\simulation\testSimulation\',selectedBand];
    if ~exist (savepath)
        mkdir(savepath)
    end
    print(gcf, [savepath,slash,'No.',num2str(n),'-',sessionName{pi},'a.png'], '-dpng', '-r300'); % 600 DPI 
    close all
    figure;
    rose(spk_phi);title('Spike Phase');
    print(gcf, [savepath,slash,'No.',num2str(n),'-',sessionName{pi},'ph.png'], '-dpng', '-r300'); % 600 DPI 

    close all

%% 
% Let's now generate a rhythmic neurons at theta that doesn't reset it's
% phase but instead shows a constant increase in spike rate at 0.3 seconds

freq=5;% centre frequency of rhythm
freq_range=0.5; % no oscillation is perfectly stationary so we add a little variation within which the frequency will vary from moment to moment
sm_freqdrft=200; % this is the smoothing window (in ms) that is needed for the frequency drift; 
sm_spkdens=250; % this is the smoothing for calculating spike density (i.e. window length of the Gaussian, in ms);
epoch = [-4 4];% start and end time of epoch in seconds; stimulus onset is 0;
ntrials = 100;% how many trials do you want to simulate?
t_er = 0.3; % start time of evoked response
jitt_er = 0.05; % no physiological processs is perfect so we att a little jitter for evoked response in seconds (i.e. 0.025 = 25 ms)
dur_er = 0;%0.1;  % duration of evoked response in seconds
fr = 0;% 5;% firing rate of evoked response in Hz
fig = 1; % set this to 1 if you want a figure, or 0 if you don't want a figure; 0 speeds up things;
spk_phs = pi; % phase in radians [-pi to pi] to which spikes are locked; 0 = peak, pi = trough; cannot be negative numbers here 

[spk_ts,time,lfp, spk_phi] = simulate_evoked_resp(freq,freq_range,sm_freqdrft,sm_spkdens,epoch,ntrials,t_er,jitt_er,dur_er, fr, spk_phs,fig);

%%

