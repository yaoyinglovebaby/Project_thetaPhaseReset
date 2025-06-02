function [spk_ts,time,lfp,spk_phi] = simulate_phase_reset(freq,freq_range,sm_freqdrft,sm_spkdens,epoch,ntrials,t_pr,jitt_pr,spk_phs,fig,noise)

%% Simulation to illustrate phase-reset in rhythmically spiking single neurons
% Input:
% freq=4;% centre frequency of rhythm in Hz
% freq_range=0.5;% scale of frequency drift (i.e. 1 = frequency will drift between freq +/- 1 Hz)
% sm_freqdrft=200;% smoothing for frequency drift; the higher this value is the slower the drift
% sm_spkdens=1000/freq;% smoothing for spk density
% epoch=[-4 4];% start/end time of epoch in seconds
% ntrials=110;% number of trials
% t_pr=0.3; % time at which phase-reset occurs in seconds
% jitt_pr=0.025;% jitter of phase reset in seconds (i.e. 0.025 = 25 ms)
% fig = 1 if you like a figure, 0 if you don't want a figure
% Output:
% spk_ts = spike times in samples
% time = time vector
% lfp = simulated local field potential that drives the neuron (average to % get ERP)
% spk_phi = vector containing spike phase in radians, concentanated accross trials (n spikes * ntrials); 
% noise = value between 1 and 0; 0 = no noise is added; > 0 = white noise is
% added with an amplitude as per the entered value

sr=1000;% sampling rate
freqdrift=0.05; 
tvec=epoch(1):1/sr:epoch(2);
cycle_length=sr/freq;% get cycle length in samples (ms)
phs_shift=(cycle_length/(2*pi))*spk_phs;

% Generate FREQUENCY DRIFT for conductor and follower dipoles
% Employs a "random walk" type of approach
drift_sig = zeros(ntrials,numel(tvec));
tmp_drft=0;
spk_dens = zeros(ntrials,numel(tvec));

% Create drift
for trl = 1:ntrials % Loop through trials
    for i = 1:numel(tvec)-1 % And timepoints
       r = -freqdrift + (freqdrift+freqdrift)*rand(1,1); % Random drift
       tmp_drft(1,i+1) = tmp_drft(1,i)+r; % Add new drift to ongoing drift
    end
    tmp_sm = smoothdata(tmp_drft,'Gaussian',sm_freqdrft); % Smooth using 1-second window
    tmp_minmax=max(abs(tmp_sm));
    drift_sig(trl,:) = (smoothdata(tmp_sm,'Gaussian',sm_freqdrft)./tmp_minmax).*freq_range; % Smooth using 1-second window
end

%figure;plot(tvec,drift_sig);
tol=-0.1;
spk_phi=[];
for n=1:ntrials   
    reset_t=(rand(1,1))*jitt_pr; %reset时间加上随机大小jitter的时间，从下面的公式可以看出来这个jitter是加在每个trial上的
    t_pr2=round(t_pr*sr+reset_t*sr)/sr;% this is a bit more complicated than it has to be but it's all to do with getting discrete indices for phase reset that vary from trial to trial
    tvec_pre=tvec(1:(t_pr2*sr)+abs(epoch(1))*sr);
    tvec_post=tvec((t_pr2*sr)+abs(epoch(1))*sr+1:end);%将时间分成reset之前的和reset之后的
    drft_pre=drift_sig(n,1:numel(tvec_pre));
    drft_post=drift_sig(n,numel(tvec_pre)+1:end);% drift也按照时间分成前后两部分
    randphs=rand(1,1)*2*pi; %产生一个随机相位，是lfp开始的相位
    osci_pre=sin(2*pi*tvec_pre.*(drft_pre+freq)+randphs); %创造lfp，reset之前从随机相位开始
    osci_post=sin(2*pi*tvec_post.*(drft_post+freq)+pi); % reset后都从PI开始 
    osci(n,:)=[osci_pre osci_post];
    anlyt_osci=hilbert(osci(n,:));
    phs_osci=angle(anlyt_osci);
    % here we subtract the spike phase and flip the signal so we can easily
    % extract the time points where phs_osci visits spk_phs
    osci4spk = abs(phs_osci - spk_phs).*-1; %越接近设置的spike的相位值越接近于0 【-6 0】     
    [~,spks{1,n}]=findpeaks(osci4spk,'MinPeakHeight',tol); %设置阈值：只找大于等于 tol 的峰值，根据lfp产生spike
    %spkts=spkts+phs_shift;
    %spks{1,n}=spkts(find(spkts<numel(tvec)));
    % spk_amp{1,n}=osci(n,spks{1,n});
    spk_dens(n,spks{1,n})=1;
    spk_dens(n,:)=smoothdata(spk_dens(n,:),'Gaussian', sm_spkdens); %smooth by trial
  % %  plotting stuff; useful for debugging 
  % figure
    spk_phi=[spk_phi phs_osci(spks{1,n})];
  %   plot(tvec,osci(n,:));
  %   hold on
  %   scatter(tvec(spks{1,n}),spk_amp{1,n});
  %   hold off
end

if fig == 1
    figure;
    subplot(10,2,1:6);imagesc(tvec,[1:n],osci);
    title('ERP image');
    subplot(10,2,7:10);plot(tvec,mean(osci,1));
    title('ERP');
    subplot(10,2,11:16);
    title('Rasterplot');
    hold on;
    for trial = 1:n
        spike_times = spks{trial}; % Get spike times for the current trial
        for spike = 1:length(spike_times)
            % Plot each spike as a vertical line segment
            line([tvec(spike_times(spike)), tvec(spike_times(spike))], [trial-0.5, trial+0.5], 'Color', 'k', 'LineWidth', 1);
        end

    end
    hold off;

    % Format the plot
    xlabel('Time (s)');
    ylabel('Trial');
    title('Raster Plot');
    ylim([0.5, n+0.5]); % Adjust y-axis limits to fit trials
    yticks(1:10:n); % Set y-ticks for each trial
    xlim([min(tvec) max(tvec)])
    subplot(10,2,17:20);
    plot(tvec,mean(spk_dens,1));
    title('Spike Density');
    figure;
    rose(spk_phi);title('Spike Phase');
end

if noise > 0
    noise = rand(size(osci));
else
    noise = zeros(size(osci));
end

spk_ts=spks;
time=tvec;
lfp=osci+noise;
