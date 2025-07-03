%% addpath & load data
clear

addpath('Z:\Ying_Phasereset\code\Scripts_BetaProject');
addpath 'Z:\Ying_Phasereset\code\TryOscore\func_internal';
addpath '\\analyse4.psy.gla.ac.uk\project0309\Ying_Phasereset\analyses';
addpath '\\analyse4.psy.gla.ac.uk\project0309\Ying_Phasereset\code';
datapath = 'Z:\Luca\ESN_code_data_19052023\data_ESN\experiment 2\';
addpath 'Z:\Ying_Phasereset\analyses\prediction\simulation'
addpath 'Z:\Ying_Phasereset\analyses\prediction'
load('Z:\Ying_Phasereset\analyses\oscore\allSpks_exp1_withESN_new.mat')

% allSpks=struct2table(allSpks);
% allSpks=allSpks(uniqVals,:);
% allSpks=table2struct(allSpks);
%% spike lfp time align

time=-4:0.001:6;
ifig=true;
ifprint=true;

savepath='Z:\Ying_Phasereset\analyses\prediction\realdata\step2_2.0\enc_1\';
if ~exist("savepath","dir")
   mkdir(savepath)
end
savepath_r='Z:\Ying_Phasereset\analyses\prediction\realdata\step2_2.0\ret_1\';
if ~exist("savepath_r","dir")
   mkdir(savepath_r)
end
%%
clearvars -regexp  ^spikeTimes
CNTE=0;cntE=0;CNTR=0;cntR=0;
for n=1:size(allSpks,1)
  %% prepare data   
  if allSpks(n).isOsci
     tmpspks = round(allSpks(n).spks);
     dt = 0:1:tmpspks(end)+1;%相当于tspan*fs
     [spkdens,~] = hist(tmpspks,dt);% transform spike timestamps into binary data (0=no spike,1=spike)
     spks.label={[allSpks(n).wirename '-' num2str(allSpks(n).su)]};
     spks.trial{1,1}=spkdens;
     spks.time{1,1}=dt./1000;
     spks.fsample = 1000;
     % Set up trial definition based on encoding and retrieval triggers
     ntrls_e=size(allSpks(n).encTrigger,1); 
     ntrls_r=size(allSpks(n).retTrigger,1); 
       
     trigger_e=allSpks(n).encTrigger(:,1);
     trigger_r=allSpks(n).retTrigger(:,1);

     trls_e=[trigger_e+time(1) trigger_e+time(end) ones(ntrls_e,1).*-4 allSpks(n).hitsIdx]; 
     trls_r=[trigger_r+time(1) trigger_r+time(end) ones(ntrls_r,1).*-4 allSpks(n).hitsIdx]; 
        
     [spkd_e] = continuous2trls(trls_e,spks);
     [spkd_r] = continuous2trls(trls_r,spks);

     for tri = 1:length(spkd_e.trial)
         spikeTimes_e{tri}=find(spkd_e.trial{tri} == 1);  
     end

     for tri = 1:length(spkd_r.trial)
         spikeTimes_r{tri}=find(spkd_r.trial{tri} == 1);  
     end
%% Classify using resetime-based prediction
    [smoothed_bin_avg_e, bin_starts, resetime_e,  ~, p_e,~, Tt_e, baseline_mean_e, baseline_std_e] = ...
        prediction_perm(spikeTimes_e, time, num2str(n), 10000, false);%uniqVals(n)
    [smoothed_bin_avg_r, ~, resetime_r, ~,p_r, ~, Tt_r, baseline_mean_r, baseline_std_r] = ...
        prediction_perm(spikeTimes_e, time,num2str(n), 10000, false);
    % close all;

    %% If classified as PR, run alignment and permutation test
     if p_e<=0.05 || p_r<=0.05
        nPerms = 10000;
        [alignmentOverTime_e, timeCenters] = itSpktimeAlignment(spikeTimes_e, time);
        [alignmentOverTime_r, ~] = itSpktimeAlignment(spikeTimes_r, time);

        %% Permutation-based significance test
        if resetime_e>0.5
            tmrang_e=[resetime_e-0.5 1];
        else
            tmrang_e=[resetime_e-0.5 resetime_e+0.5];
        end
        if resetime_r>0.5
            tmrang_r=[resetime_r-0.5 1];
        else
            tmrang_r=[resetime_r-0.5 resetime_r+0.5];
        end
        [permP_e, isSignificant_e] = perm_test_1d_diff(alignmentOverTime_e, ...
            timeCenters, resetime_e, tmrang_e, nPerms);
        [permP_r, isSignificant_r] = perm_test_1d_diff(alignmentOverTime_r, ...
            timeCenters, resetime_r, tmrang_r, nPerms);
        %% label
        if isSignificant_e && p_e<=0.05
            allSpks(n).alnEnc3=1;
            Tt_e=['PR Neuron - ',num2str(n)];
        else
            allSpks(n).alnEnc3=0;
            Tt_e=['Non-PR Neuron - ',num2str(n)];
        end
        if isSignificant_r&& p_r<=0.05
            allSpks(n).alnRet3=1;
            Tt_r=['PR Neuron - ',num2str(n)];
        else
            allSpks(n).alnRet3=0;
            Tt_r=['Non-PR Neuron - ',num2str(n)];
        end
        %% Plot results
        if ifig
            %------encoding ------------------------------------------------
            figure('Units','centimeters','Position',[5 5 14 18]);
            %% Subplot 1: Resetime error trace
            subplot(3,1,1)
            plot(bin_starts, smoothed_bin_avg_e, '-o', 'LineWidth', 1.5); hold on;

            % Add z-threshold line
            threshold = 1.64; % You can modify this threshold
            threshold_line = baseline_mean_e + threshold * baseline_std_e;
            yline(threshold_line, '--', ['z = ', num2str(threshold)], 'Color', [1, 0.1, 0]);

            % Annotate p-value
            x_range = xlim;
            y_range = ylim;
            % Replace "p_perm" with actual variable or permP
            text(x_range(1) + 0.05 * range(x_range), ...
                 y_range(2) - 0.05 * range(y_range), ...
                 sprintf('p (perm) = %.3f', p_e), ...
                 'FontSize', 12, 'Color', 'k', 'HorizontalAlignment', 'left');

            xlabel('Time (s)');
            ylabel('Mean Relative Error');
            title(Tt_e);

            %% Subplot 2: Spike alignment
            subplot(3,1,2)
            % plot(timeCenters, alignmentOverTime_e, 'k'); hold on;
            plot(timeCenters, smoothdata(alignmentOverTime_e, 'gaussian', 15), 'k'); hold on;
            xline(resetime_e, '--r'); hold on;
     
            if isSignificant_e
                
                title(sprintf('Significant alignment: p = %.4f', permP_e));
            else
                title(sprintf('Non-significant alignment: p = %.4f', permP_e));
            end
            xlabel('Time (s)'); ylabel('Spike alignment');
            
            %% Subplot 3: Rastor
            subplot(3,1,3)
            title('Rasterplot');
            hold on;
            for trial = 1:length(spikeTimes_e)
                spike_times = spikeTimes_e{trial}; % Get spike times for the current trial
                for spike = 1:length(spike_times)
                    % Plot each spike as a vertical line segment
                    line([ time(spike_times(spike)),  time(spike_times(spike))], [trial-0.5, trial+0.5], 'Color', 'k', 'LineWidth', 1);
                end

            end
            xlim([-1 1])
            xlabel('T0=enc1')
            
            % -----Save alignment figure ---------------------------
            if ifprint
                print(gcf, [savepath, Tt_e, '_cutime.png'], '-dpng', '-r600');
            end
            close all;
            %-------------------------------------------------------
            
            %------retrieval ---------------------------
           
            %% Subplot 1: Resetime error trace
            figure('Units','centimeters','Position',[5 5 14 18]);
            subplot(3,1,1)
            plot(bin_starts, smoothed_bin_avg_r, '-o', 'LineWidth', 1.5); hold on;

            % Add z-threshold line
            threshold = 1.64; % You can modify this threshold
            threshold_line = baseline_mean_r + threshold * baseline_std_r;
            yline(threshold_line, '--', ['z = ', num2str(threshold)], 'Color', [1, 0.1, 0]);

            % Annotate p-value
            x_range = xlim;
            y_range = ylim;
            % Replace "p_perm" with actual variable or permP
            text(x_range(1) + 0.05 * range(x_range), ...
                 y_range(2) - 0.05 * range(y_range), ...
                 sprintf('p (perm) = %.3f', p_r), ...
                 'FontSize', 12, 'Color', 'k', 'HorizontalAlignment', 'left');

            xlabel('Time (s)');
            ylabel('Mean Relative Error');
            title(Tt_r);

            %% Subplot 2: Spike alignment
            subplot(3,1,2)
            % plot(timeCenters, alignmentOverTime_r, 'k'); hold on;
            plot(timeCenters, smoothdata(alignmentOverTime_r, 'gaussian', 15), 'k'); hold on;
            
            xline(resetime_r, '--r'); hold on;
     
            if isSignificant_r
                
                title(sprintf('Significant alignment: p = %.4f', permP_r));
            else
                title(sprintf('Non-significant alignment: p = %.4f', permP_r));
            end
            xlabel('Time (s)'); ylabel('Spike alignment');
            
            %% Subplot 3: Rastor
            subplot(3,1,3)
            title('Rasterplot');
            hold on;
            for trial = 1:length(spikeTimes_r)
                spike_times = spikeTimes_r{trial}; % Get spike times for the current trial
                for spike = 1:length(spike_times)
                    % Plot each spike as a vertical line segment
                    line([ time(spike_times(spike)),  time(spike_times(spike))], [trial-0.5, trial+0.5], 'Color', 'k', 'LineWidth', 1);
                end

            end
            xlim([-1 1])
            xlabel('T0=ret1')

            % -----Save alignment figure ---------------------------
            if ifprint
               print(gcf, [savepath_r, Tt_r, '_cutime.png'], '-dpng', '-r600');
            end
            close all;
            %-------------------------------------------------------
        end
     end
  end
end



