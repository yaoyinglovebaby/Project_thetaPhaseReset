function [permP, isSignificant] = perm_test_1d_diff(dataVec, timeVec, cutime, edgetime, nPerms)
    % dataVec: 1D timecourse, e.g. alignmentOverTime
    % timeVec: corresponding time axis
    % baselineWin: [start, end] of baseline
    % effectWin: [start, end] of test window
    % nPerms: number of permutations
    baselineWin = [edgetime(1) cutime];
    effectWin = [cutime edgetime(2)];  
    baseIdx = find(timeVec >= baselineWin(1) & timeVec <= baselineWin(2));
    effIdx  = find(timeVec >= effectWin(1) & timeVec <= effectWin(2));
    totalIdx = [baseIdx, effIdx];
    nBase = numel(baseIdx);

    baselineMean = mean(dataVec(baseIdx));
    effectMean   = mean(dataVec(effIdx));
    realDiff     = effectMean - baselineMean;

    nullDist = zeros(1, nPerms);
    for p = 1:nPerms
        shuffledIdx = totalIdx(randperm(length(totalIdx)));
        fakeBase = shuffledIdx(1:nBase);
        fakeEff  = shuffledIdx(nBase+1:end);
        nullDist(p) = mean(dataVec(fakeEff)) - mean(dataVec(fakeBase));
    end

    permP = (sum(nullDist >= realDiff) + 1) / (nPerms + 1);
    isSignificant = permP < 0.05;
end
