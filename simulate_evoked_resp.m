function [spk_ts,time,lfp, spk_phi] = simulate_evoked_resp(freq,freq_range,sm_freqdrft,sm_spkdens,epoch,ntrials,t_er,jitt_er,dur_er, fr, spk_phs,fig)

%% Simulation to illustrate phase-reset in rhythmically spiking single neurons
% % Input:
% freq=4;% centre frequency of rhythm in Hz
% freq_range=0.5;% scale of frequency drift (i.e. 1 = frequency will drift between freq +/- 1 Hz)
% sm_freqdrft=200;% smoothing for frequency drift; the higher this value is the slower the drift
% sm_spkdens=1000/freq;% smoothing for spk density
% epoch=[-4 4];% start/end time of epoch in seconds
% ntrials=110;% number of trials
% t_er=0.3; % time at which evoked response occurs in seconds
% jitt_er=0.025;% jitter of evoked response in seconds (i.e. 0.025 = 25 ms)
% dur_er = 0.2; % duration of ERP in seconds
% fr = 2; % firing rate for evoked response in Hz
% fig = 1; % if you like a figure, 0 if you don't want a figure
% Output:
% spk_ts = spike times in samples
% time = time vector
% lfp = simulated local field potential that drives the neuron (average to % get ERP)
% spk_phi = vector containing spike phase in radians, concentanated accross trials (n spikes * ntrials); 


sr=1000;% sampling rate
freqdrift=0.05; 
tvec=epoch(1):1/sr:epoch(2);
%jitt_pr=(sr/freq)/jitt_pr/1000;

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
    randphs=rand(1,1)*2*pi;
    osci=sin(2*pi*tvec.*(drift_sig(n,:)+freq)+randphs);   % back ground oscillations goes on without phase reset 
    er_t=(rand(1,1))*jitt_er;
    t_pr2=round(t_er*sr+er_t*sr)/sr;% this is a bit more complicated than it has to be but it's all to do with getting discrete indices for phase reset that vary from trial to trial
    n_pre=(t_pr2*sr)+abs(epoch(1))*sr;
    n_post=numel(tvec)-n_pre;
    erp=[zeros(1,n_pre-dur_er*sr/2) hanning(dur_er*sr)' zeros(1,n_post-dur_er*sr/2)];
    lfp(n,:)=osci+erp;% build ERP by adding evoked response on top of background oscillation
    anlyt_lfp=hilbert(lfp(n,:));% extract phase for compound signal
    phs_lfp=angle(anlyt_lfp);
    anlyt_osci=hilbert(osci);% extract phase for oscillation only
    phs_osci=angle(anlyt_osci);
    % here we subtract the spike phase and flip the signal so we can easily
    % extract the time points where phs_osci visits spk_phs
    osci4spk = abs(phs_osci - spk_phs).*-1;
    [~,spks_osci]=findpeaks(osci4spk,'MinPeakHeight',tol);% spikes driven by background oscillation
    % add evoked response spikes
    er_idx=find(erp>0);
    n_spks=round((fr/sr)*numel(er_idx));
    spks_er=er_idx(randperm(numel(er_idx),n_spks));% put random spikes during evoked response window
    spks{1,n}=[spks_osci spks_er];% add evoked response to background oscillation
    %spks{1,n}=find(osci(n,:)>=threshld);
    spk_dens(n,spks{1,n})=1;
    spk_dens(n,:)=smoothdata(spk_dens(n,:),'Gaussian', sm_spkdens);
    spk_phi=[spk_phi phs_lfp(spks{1,n})];
end

if fig == 1
    figure;
    subplot(10,2,1:6);imagesc(tvec,[1:n],lfp);
    title('ERP image');
    subplot(10,2,7:10);plot(tvec,mean(lfp,1));
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

spk_ts=spks;
time=tvec;
% lfp=osci; %I think it shouldn't be? 2025.04.19yy