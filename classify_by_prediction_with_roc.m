function classify_by_prediction_with_roc(cfg, ifig, type)
addpath 'Z:\Ying_Phasereset\toolbox\Violinplot-Matlab-master'
    ori_labels = [];       % Ground truth: 1 = phase reset, 0 = ERP
    predictions = [];      % Classification results
    scorevalue = [];
    prop_sig = [];

    for i = 1:length(type)
        %% Extract parameters
        centerFreq = cfg(i).centerFreq;
        freqVariation = cfg(i).freqVariation;
        freqDriftWindow = cfg(i).freqDriftWindow;
        spikeDensitySmoothing = cfg(i).spikeDensitySmoothing;
        epochDuration = cfg(i).epochDuration;
        numTrials = cfg(i).numTrials;
        eventTime = cfg(i).eventTime;
        eventJitter = cfg(i).eventJitter;
        spikePhase = cfg(i).spikePhase;
        firingRate = cfg(i).firingRate;
        responseDuration = cfg(i).responseDuration;
        noise = cfg(i).noise;
        plotFigures = cfg(i).plotFigures;

        %% Simulate and classify
        switch type{i}
            case 'pr'
                [spikeTimes, timeVec, ~, ~] = simulate_phase_reset(centerFreq, freqVariation, freqDriftWindow, ...
                    spikeDensitySmoothing, epochDuration, numTrials, eventTime, eventJitter, spikePhase, plotFigures, noise);
                ori_labels(i) = 1;
                ori_label{1}='PhaseReset';

            case 'pRhy'
                [spikeTimes, timeVec, ~, ~] = simulate_evoked_resp(centerFreq, freqVariation, freqDriftWindow, ...
                    spikeDensitySmoothing, epochDuration, numTrials,  0, 0, 0, centerFreq, spikePhase, plotFigures);
                ori_labels(i) = 0;
                ori_label{2}='pure Rhythmic';
            case 'ERP'
                [spikeTimes, timeVec, ~, ~] = simulate_evoked_resp(centerFreq, freqVariation, freqDriftWindow, ...
                    spikeDensitySmoothing, epochDuration, numTrials, eventTime, eventJitter, responseDuration, firingRate, spikePhase, plotFigures);
                ori_labels(i) = 0;
                ori_label{2}='ERP';
        end
        % [fig, resetime, score_value, p_val, h_vec, p_vec] = prediction1_perm(spikeTimes, timeVec, 'auto');
        % [fig, ~, scorevalue(i), qmin, hvec, qvec] = prediction1(spikeTimes, timeVec, type{i});
        [fig, ~, scorevalue(i), ~,~,predictions(i), hvec,~] = prediction_perm(spikeTimes, timeVec, type{i});
        print(fig,['Z:\Ying_Phasereset\analyses\prediction\simulation\test_ERP_perm\',num2str(i),' shift_values for ',type{i},'.png'], '-dpng')
        % 
        % [~, ~, scorevalue(i), prop_sig(i)] = prediction1(spikeTimes, timeVec, 'auto');
        close all;

        % % Classification rule
        % if any(hvec) && scorevalue(i) > 0
        %     predictions(i) = 1;  % phase reset
        % else
        %     predictions(i) = 0;  % ERP or uncertain
        % end
    end

   

    %% Overall ROC
    [X, Y, ~, AUC] = perfcurve(ori_labels, scorevalue, 1);

    %% Confusion matrix (FDR-significant neurons only)
    % sig_idx = find(h == 1);
    % ori_sig = ori_labels(sig_idx);
    % pred_sig = predictions(sig_idx);

    TP = sum((predictions == 1) & (ori_labels == 1));
    TN = sum((predictions == 0) & (ori_labels == 0));
    FP = sum((predictions == 1) & (ori_labels == 0));
    FN = sum((predictions == 0) & (ori_labels == 1));

    fprintf('Confusion Matrix :\n');
    fprintf(['                Predicted ',ori_label{2},'   Predicted ',ori_label{1},'\n']);
    fprintf(['Actual ',ori_label{2},'       %10d        %10d\n'], TN, FP);
    fprintf(['Actual ',ori_label{1},'%10d        %10d\n'], FN, TP);
    fprintf('AUC (all neurons): %.4f\n', AUC);
    % fprintf('Significant neurons after FDR: %d / %d\n\n', length(sig_idx), length(h));

    %% Classification result bar plot
    % figure;
    % bar([TN, FP, FN, TP]);
    % set(gca, 'XTickLabel', {'True ERP', 'False PR', 'False ERP', 'True PR'});
    % ylabel('Neuron Count');
    % title('Classification Results ');

    %% ROC plot
    if ifig
        figure;
        plot(X, Y, 'LineWidth', 2);
        hold on;
        plot([0 1], [0 1], '--k');
        xlabel('False Positive Rate');
        ylabel('True Positive Rate');
        title('ROC Curve for Phase Reset Detection (All Neurons)');
        legend(sprintf('AUC = %.4f', AUC), 'Location', 'SouthEast');
        grid on;
    end


    %% Violin + Box Plot of shiftvalue (FDR-significant only)
% h 是 FDR 显著性逻辑向量

% %% 只保留 FDR 显著的单元
% sig_idx   = find(h);
% shift_sig = shiftvalue_sum(sig_idx);
% label_sig = ori_labels(sig_idx);

% 构造标签 cell 数组
labels = cell(size(ori_labels));
labels(ori_labels==0) = ori_label(2);
labels(ori_labels==1) = ori_label(1);

% 颜色定义为 cell，每个元素是 1x3 的 RGB 向量
colors = {
    [0.85 0.33 0.10];    % ERP: 橙
    [0.00 0.45 0.74];    % PR : 蓝
};


% 分组
group1 = scorevalue(strcmp(labels, ori_label{2}));
group2 = scorevalue(strcmp(labels, ori_label{1}));

% 统计检验（非参数 Mann-Whitney U test）
[p, h, stats] = ranksum(group1, group2);


figure;
% --- 子图 1：Violin Plot ---
subplot(1,2,1);
v = violinplot(scorevalue, labels);
for i = 1:length(v)
    % 将 colors{i} 再包成 cell：colors(i)
    v(i).ViolinColor = colors(i);
    v(i).EdgeColor   = colors{i};
    % 中位数线加粗
    v(i).MedianPlot.LineWidth = 2;
    % 去掉 boxplot 填充
    v(i).BoxPlot.FaceColor = 'w';
end
ylabel('Shift Value');
title('Shift Value Distribution (Violin Plot)');
set(gca, 'FontSize', 12);

% --- 子图 2：Box Plot with fill ---
subplot(1,2,2);
boxplot(scorevalue, labels, 'Colors', 'k', 'Symbol', '');
hold on;
hBoxes = findobj(gca, 'Tag','Box');
% findobj 返回对象按绘制顺序倒序排列，需要反向索引
for i = 1:length(hBoxes)
    idx = length(hBoxes) - i + 1;
    xData = get(hBoxes(idx), 'XData');
    yData = get(hBoxes(idx), 'YData');
    % 直接用 colors{i}，它是 RGB 向量
    patch(xData, yData, colors{i}, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
end
ylabel('Shift Value');
title('Shift Value Distribution (Box Plot)');
set(gca, 'FontSize', 12);

% 添加统计显著性标记（星号）
y_max = max([group1, group2]);
y_sig = y_max + 2;  % 星号位置
line([1 2], [y_sig y_sig], 'Color', 'k', 'LineWidth', 1.5);
text(1.5, y_sig + 1, getSigStar(p), 'FontSize', 16, 'HorizontalAlignment', 'center');

hold off;

% --- 辅助函数：根据 p 值返回星号表示 ---



end
