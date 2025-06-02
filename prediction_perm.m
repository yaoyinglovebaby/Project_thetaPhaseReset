function [fig, resetime, score_value, h_fdr, p_perm, prediction, h_perm_fdr, p_perm_all] = prediction_perm(spikeTimes, timeVec, type)
% PREDICTION_PERM   Detect phase-reset vs ERP in one neuron using binned error with optional permutation test
%
%   Inputs:
%     spikeTimes: cell array, each cell is a vector of spike times (s) for one trial
%     timeVec:    1xT vector of time points at which to compute prediction error
%     type:       string, used in the plot title
%     do_cluster_perm: boolean, whether to run cluster-based permutation
%
%   Outputs:
%     fig:         handle to diagnostic figure (empty if no significant effect)
%     resetime:    time (in seconds) of maximal significant deviation
%     score_value: signed magnitude of deviation (positive = phase reset, negative = ERP)
%     h_fdr:       binary vector indicating significant bins (after FDR on Wilcoxon p-values)
%     p_perm:      overall permutation-based significance (1 = significant, 0 = not)
%     prediction:  1 if significant, 0 if not
%     h_perm_fdr:  binary vector indicating bins significant via permutation + FDR
%     p_perm_all:  vector of permutation-based p-values per bin

    %% Parameters
    binSize = 0.1;             % 100 ms bin
    stepSize = 0.03;           % 20 ms step
    min_valid_trials = 2;
    alpha = 0.05;
    n_perm = 1000;
    threshold = 1.64;             % z-threshold for baseline comparison

    %% Initialization
    n_trials = numel(spikeTimes);
    n_times = numel(timeVec);
    all_err = nan(n_trials, n_times);

    %% Step 1: Construct error matrix
    for tr = 1:n_trials
        st = spikeTimes{tr};
        rel = nan(1, n_times);
        for k = 3:numel(st)
            t1 = st(k-2);
            t2 = st(k-1);
            t3 = st(k);
            T = t2 - t1;
            t_pred = t2 + T;
            if t3 > 0 && t3 <= n_times
                rel(t3) = (t3 - t_pred) / T;
            end
        end
        all_err(tr, :) = rel;
    end

    %% Step 2: Time binning
    mask = find((timeVec >= -1) & (timeVec <= 1));
    t_window = timeVec(mask);
    bin_starts = t_window(1):stepSize:(t_window(end) - binSize);
    n_bins = numel(bin_starts);

    %% Step 3: Mean error per bin per trial
    mean_err_bins = nan(n_trials, n_bins);
    for b = 1:n_bins
        bin_mask = (timeVec >= bin_starts(b)) & (timeVec < bin_starts(b) + binSize);
        trial_means = nanmean(all_err(:, bin_mask), 2);
        if sum(~isnan(trial_means)) >= min_valid_trials
            mean_err_bins(:, b) = trial_means;
        end
    end

    %% Step 4: Wilcoxon test
    p_vec = nan(1, n_bins);
    for b = 1:n_bins
        vals = mean_err_bins(:, b);
        if sum(~isnan(vals)) >= 5
            p_vec(b) = signrank(vals, 0);
        end
    end

    %% Step 4.1: FDR correction
    [h_fdr, ~, ~, ~] = fdr_bh(p_vec, alpha);

    %% Step 5: Baseline thresholding
    baseline_idx = bin_starts >= -1 & bin_starts < 0;
    task_idx = bin_starts >= 0 & bin_starts <= 1;

    baseline_vals = mean_err_bins(:, baseline_idx);
    baseline_mean = nanmean(baseline_vals(:));
    baseline_std  = nanstd(baseline_vals(:));
    bin_avg = nanmean(mean_err_bins, 1);
    z_value=(bin_avg - baseline_mean) / baseline_std;
    h_sig = z_value > threshold ;
    % h_sig = abs(bin_avg - baseline_mean) > threshold * baseline_std;
    sig_idx = find(h_sig & task_idx);

    if numel(sig_idx) < 1
        resetime = 0;
        score_value = 0;
        prediction = 0;
        h_perm_fdr = [];
        p_perm_all = [];
         p_perm=0;
         fig = figure;
    subplot(2,1,1);
    hold on;
    plot(bin_starts, bin_avg, '-o', 'LineWidth', 1.5);
    scatter(bin_starts(h_fdr == 1), bin_avg(h_fdr == 1), 80, 'r', 'filled');
    yline(0, '--k');
    xlabel('Time (s)');
    ylabel('Mean Relative Error');
    title(['Wilcoxon + FDR - ' type]);

    subplot(2,1,2);
    hold on;
    plot(bin_starts, bin_avg, '-o', 'LineWidth', 1.5);
    scatter(bin_starts(h_perm_fdr == 1), bin_avg(h_perm_fdr == 1), 80, 'r', 'filled');
    yline(0, '--k');
    xlabel('Time (s)');
    ylabel('Mean Relative Error');
    title(['Permutation + FDR - ' type]);
        return;
    end

    %% Step 6: Permutation test (shuffle all values)
    h_sig_perm = false(n_perm, n_bins);
    p_all = zeros(n_perm, 1);

    for p = 1:n_perm
        shuffled = mean_err_bins(:, randperm(n_bins));
        bin_avg_perm = nanmean(shuffled, 1);

        bl_perm_vals = shuffled(:, baseline_idx);
        bl_mean_perm = nanmean(bl_perm_vals(:));
        bl_std_perm = nanstd(bl_perm_vals(:));
        % if bl_std_perm < 1e-8
        %     continue;
        % end
        z_perm=(bin_avg_perm - bl_mean_perm) / bl_std_perm;
        h_sig_perm(p,:) = z_perm > threshold ;
        % h_sig_perm(p,:) = abs(bin_avg_perm - bl_mean_perm) > threshold * bl_std_perm;
        sig_idx_perm = find(h_sig_perm(p,:) & task_idx);
        if numel(sig_idx_perm) >= numel(sig_idx)
            p_all(p) = 1;
        end
    end

    p_perm = mean(p_all);
    p_perm_all = mean(h_sig_perm, 1);
    p_perm_all(p_perm_all == 0) = NaN;
    [h_perm_fdr, ~, ~, ~] = fdr_bh(p_perm_all, alpha);

    %% Step 7: Final scoring
    if p_perm < alpha
        [~, imax] = max(abs(bin_avg(sig_idx)));
        resetime = bin_starts(sig_idx(imax));
        score_value = bin_avg(sig_idx(imax));
        prediction = 1;
    else
        resetime = 0;
        score_value = 0;
        prediction = 0;
        
         fig = figure;
    subplot(2,1,1);
    hold on;
    plot(bin_starts, bin_avg, '-o', 'LineWidth', 1.5);
    scatter(bin_starts(h_fdr == 1), bin_avg(h_fdr == 1), 80, 'r', 'filled');
    yline(0, '--k');
    xlabel('Time (s)');
    ylabel('Mean Relative Error');
    title(['Wilcoxon + FDR - ' type]);

    subplot(2,1,2);
    hold on;
    plot(bin_starts, bin_avg, '-o', 'LineWidth', 1.5);
    scatter(bin_starts(h_perm_fdr == 1), bin_avg(h_perm_fdr == 1), 80, 'r', 'filled');
    yline(0, '--k');
    xlabel('Time (s)');
    ylabel('Mean Relative Error');
    title(['Permutation + FDR - ' type]);
        return;
    end

    %% Step 8: Plot
    fig = figure;
    subplot(2,1,1);
    hold on;
    plot(bin_starts, bin_avg, '-o', 'LineWidth', 1.5);
    scatter(bin_starts(h_fdr == 1), bin_avg(h_fdr == 1), 80, 'r', 'filled');
    yline(0, '--k');
    xlabel('Time (s)');
    ylabel('Mean Relative Error');
    title(['Wilcoxon + FDR - ' type]);

    subplot(2,1,2);
    hold on;
    plot(bin_starts, bin_avg, '-o', 'LineWidth', 1.5);
    scatter(bin_starts(h_perm_fdr == 1), bin_avg(h_perm_fdr == 1), 80, 'r', 'filled');
    yline(0, '--k');
    xlabel('Time (s)');
    ylabel('Mean Relative Error');
    title(['Permutation + FDR - ' type]);
end
