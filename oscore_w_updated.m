%% addpath & load data
clear

addpath('Z:\Ying_Phasereset\code\Scripts_BetaProject');
addpath 'Z:\Ying_Phasereset\code\TryOscore\func_internal';
addpath '\\analyse4.psy.gla.ac.uk\project0309\Ying_Phasereset\analyses';
addpath '\\analyse4.psy.gla.ac.uk\project0309\Ying_Phasereset\code';
%datapath = 'Z:\Luca\ESN_code_data_19052023\data_ESN\experiment 1\';
datapath = 'Z:\Luca\ESN_code_data_19052023\data_ESN\experiment 2\';
% load([datapath,'allSpksFull.mat'],'allSpks');
load('\\analyse4.psy.gla.ac.uk\project0309\Luca\ESN_code_data_19052023\data_ESN\experiment 2\allSpks.mat');
% spks2plot=allSpks;
% cnt=0;
% for nii=1:585
%     cnt=cnt+1;
% spks2plot(nii).no=cnt;
% end  
allSpks=spks2plot;
%% set parameter
%
band_all={'Theta', 'Alpha', 'Beta', 'lowGamma' 'highGamma'};
for si=1:5
computecs=false;
fs = 1000;
plotsmooth=false;
plotfigure=false;
reduceNonspk=true;
selectedBand = band_all{si};  %  'Theta', 'Alpha', 'Beta', 'lowGamma' 'highGamma'
delta_f=0.25;
% 
time=-4:0.001:6;
switch selectedBand
    case 'Theta'
        Gwin = 125;  % Theta (4-8 Hz): highfrq border 60 ms
        FMin = 3;  FMax =8;  
        % w = 2^nextpow2(max(3 * fs / FMin, fs / 4)); % method on paper
        %1000 (3 CIRCLES OF LOWER FREQBAND)
    case 'Alpha'
        Gwin = 80;  % Alpha (8-12 Hz):  40 ms
        FMin = 8;  FMax =12;
       
        w =  2^nextpow2(max(3 * fs / FMin, fs / 4)); %375
    case 'Beta'
        Gwin = 35;  % Beta (13-30 Hz):  20 ms
        FMin = 13;  FMax =30;
        %230 CorrWind =w;
        w =  2^nextpow2(max(3 * fs / FMin, fs / 4));
    case 'lowGamma'%need rerun %OK
        Gwin = 20;  % Gamma (30-50 Hz):  20 ms
        FMin = 30;  FMax =50;
    case 'highGamma'
        Gwin = 10;  % Gamma (30-100 Hz):  10 ms
        FMin = 50;  FMax =100;
        % w =  2^nextpow2(max(3 * fs / FMin, fs / 4));
    otherwise
        error('The specified band is not found, select Theta, Alpha, Beta or Gamma.');
end

CorrWind= calculate_w(FMin, FMax, delta_f, fs);
   %%      
h = waitbar(0, '处理中...');  
for n=1:size(allSpks,2)
    waitbar(n/size(allSpks,2), h, sprintf('正在进行 %d/%d', n, size(allSpks,2)));   
  %% prepare data
        tmpspks = round(allSpks(n).spks);
        dt = 0:1:tmpspks(end)+1;%相当于tspan*fs
        [spkdens,~] = hist(tmpspks,dt);% transform spike timestamps into binary data (0=no spike,1=spike)
        spkdens_smoothed = conv(spkdens(2:end-1),gausswin(Gwin),'same');% compute smoothed firing rate using gaussian kernel

  %% oscore
        ach= xcorr(spkdens_smoothed,spkdens_smoothed,CorrWind);
        [OS, OFq, ach_nopeak, S,f,Fall,Sall] = OScore_ACH(ach, FMin, FMax, fs);
%       figure; plot(ach_nopeak);

  %% statistic
        num_permutations = 500;
        segment_length =3*length(-1:0.001:6);
        [oscore_permuted,OscFreq_permuted,z_perm,p_perm] = OscScoreStats(spkdens,OS,num_permutations,segment_length,FMin, FMax, fs,Gwin,CorrWind);

        if ~isempty(OFq) && p_perm<=0.05 && OS~=-1
            allSpks(n).isOsci=true;
            eval(['allSpks(n).',selectedBand,'= 1;']);
            eval(['allSpks(n).',selectedBand,'F= num2str(OFq);']);
            eval(['allSpks(n).',selectedBand,'_p= num2str(p_perm);']);
    
        else 
            allSpks(n).isOsci=false;
            eval(['allSpks(n).',selectedBand,'= 0;']);
            eval(['allSpks(n).',selectedBand,'F= 0;']);
            eval(['allSpks(n).',selectedBand,'_p= 0;']);
        end

    %% plot
     %    spks.label={[allSpks(n).wirename '-' num2str(allSpks(n).su)]};
     %    spks.trial{1,1}=spkdens;
     %    spks.time{1,1}=dt./1000;
     %    spks.fsample = 1000;
     %    % Set up trial definition based on encoding and retrieval triggers
     %    ntrls_e=size(allSpks(n).encTrigger,1);
     %    ntrls_r=size(allSpks(n).retTrigger,1);
     %    trls_e=[allSpks(n).encTrigger(:,1)+time(1) allSpks(n).encTrigger(:,1)+time(end) ones(ntrls_e,1).*-4 allSpks(n).hitsIdx]; 
     %    trls_r=[allSpks(n).retTrigger(:,1)+time(1) allSpks(n).retTrigger(:,1)+time(end) ones(ntrls_r,1).*-4 allSpks(n).hitsIdx]; 
     %    [spkd_e] = continuous2trls(trls_e,spks);
     %    [spkd_r] = continuous2trls(trls_r,spks);
     % 
     % 
     % if ~isempty(OFq)&&p_perm<=0.05
     %    screenSize = get(0, 'ScreenSize');  % [left, bottom, width, height]
     %    figureWidth = screenSize(3) * 0.4;   % 设定图形宽度为屏幕宽度的 60%
     %    figureHeight = figureWidth * 4/3.5;  % 高度根据 4:3 比例计算
     %    left = (screenSize(3) - figureWidth) / 2;
     %    bottom = (screenSize(4) - figureHeight) / 2;
     % 
     %    if allSpks(n).suaMua==1
     %        suaMuatmp='Sua';
     %    else
     %        suaMuatmp='Mua';
     %    end
     % 
     %    figure('Position', [left, bottom, figureWidth, figureHeight],'Visible', 'off'); %,'Visible', 'off'
     % 
     %    % 1 plot waveshape 
     %    wavf=allSpks(n).waveshapes;
     %    subplot('Position', [0.05, 0.7, 0.425, 0.25])
     %    density_plot(wavf,(1:size(wavf,1)), 1);
     %    title([suaMuatmp,' no.',num2str(allSpks(n).no),': ',[allSpks(n).wirename '-' num2str(allSpks(n).su)]],'FontSize',14);
     % 
     %    % 2 Plot spectrum of ach_nopeak
     % 
     %    % eval(['p_perm=spks2plot(n).',selectedBand,'_p;']);
     %    subplot('Position', [0.525, 0.7, 0.425, 0.25]);
     %    plot(Fall,Sall);
     %    axis ([min(f) max(f),0,([max(S)])*1.25]); 
     %    hold on;
     %    [maxVal, idx] = max(S);
     %    maxFreq = f(idx);
     %    plot(maxFreq, maxVal, 'ro');  % 红色圆圈标记最高点
     % 
     %    text(maxFreq, maxVal, sprintf('  %.2f Hz', maxFreq), 'VerticalAlignment', 'bottom', 'Color', 'red','FontSize',14, 'FontWeight', 'bold');
     %    % xlabel('Frequency [Hz]','FontSize',14);
     %    ylabel('Magnitude','FontSize',14);
     %    % title(['OFreq:',num2str(OFq),' OScore:',num2str(OS),' P(perm):',num2str(p_perm)],'FontSize',14)
     %    title('Spectrum of Autocorrelation Histogram')
     %    title(['OFreq:',num2str(OFq,'%.2f'),' OScore:',num2str(OS,'%.2f'),' P(perm):',num2str(p_perm,'%.2f')])
     % 
     %    % 3 raster plot
     %    spke_mat=vertcat(spkd_e.trial{:});
     %    if isfield(spkd_r,'trial')
     %    spkr_mat=vertcat(spkd_r.trial{:}); 
     %    spk_mat=[spke_mat;spkr_mat];
     %    else
     %    spk_mat=spke_mat;
     %    end
     %    spk_ts={find(spk_mat')}; % trials included in period from spike start to end
     %    xtix=[1 [1 2 3 4 5 6 7 8 9 10.001].*1000];
     %    xlabs=num2cell(round(time(xtix)));
     %    ntrls=ntrls_e+ntrls_r;
     % 
     %    subplot('Position', [0.05, 0.375, 0.9, 0.25]); 
     %    rasterplot(spk_ts,ntrls,numel(time),xtix,xlabs);   
     %    title('Rasterplot(Enco+R)','FontSize',14);
     % 
     % 
     %    % 4 plot spike density          
     %    subplot('Position', [0.05, 0.05, 0.9, 0.25]) 
     %    ISI = diff(tmpspks)';
     %    histogram(ISI, 'BinWidth', 0.001, 'Normalization', 'probability'); % 设置柱宽为 0.1 秒，归一化为概率密度
     %    xlabel('Inter-Spike Interval (ms)');
     %    ylabel('Proportion');
     %    title('Inter-Spike Interval (ISI) Histogram');
     %    grid on;
     % 
     %    %% save plot
     %    print(gcf, ['\\analyse4.psy.gla.ac.uk\project0309\Ying_Phasereset\analyses\oscore\rhythmic\',selectedBand,num2str(n,'%02d'),...
     %        '_w2048spksmth'], '-dpng', '-r600'); % 600 DPI 
     %    close all
     % end
end
end
% remember to change filename every time you store. donot use the same many
% times otherwise the file will be damaged and you will need to run again!
save('Z:\Ying_Phasereset\analyses\oscore\allSpks_exp1_withESN_new.mat', 'allSpks', '-v7.3');

