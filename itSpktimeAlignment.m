function [alignmentOverTime, timeCenters] = itSpktimeAlignment(spikeTimes, timeVec)

    %% Params
    dt = 0.001;
    timeAxis = -1:dt:1;
    winSize = 0.05;
    stepSize = 0.01;
    winSamples = round(winSize / dt);
    stepSamples = round(stepSize / dt);
    nTrials = length(spikeTimes);
    nSteps = floor((length(timeAxis) - winSamples) / stepSamples) + 1;
    alignmentOverTime = zeros(1, nSteps);
    timeCenters = zeros(1, nSteps);

    % gaussSigma = 0.02;
    % gaussWin = -0.1:dt:0.1;
    % gaussKernel = exp(-gaussWin.^2 / (2 * gaussSigma^2));
    % gaussKernel = gaussKernel / sum(gaussKernel);

    %% Compute alignment
    for s = 1:nSteps
        idxStart = (s - 1) * stepSamples + 1;
        idxEnd = idxStart + winSamples - 1;
        trialWindow = timeAxis(idxStart:idxEnd);
        timeCenters(s) = mean(trialWindow);

        smoothed = zeros(nTrials, winSamples);
        for ni = 1:nTrials
            aligned = (spikeTimes{ni} - find(timeVec == 0)) / 1000;
            spikeTrain = histcounts(aligned, [trialWindow, trialWindow(end) + dt]);
            smoothed(ni, :) = spikeTrain;%conv(spikeTrain, gaussKernel, 'same');
        end
        R = corr(smoothed');
        alignmentOverTime(s) = nanmean(R(triu(true(nTrials), 1)));
    end

      
end

