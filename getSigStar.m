function star = getSigStar(p)
    if p < 0.001
        star = '***';
    elseif p < 0.01
        star = '**';
    elseif p < 0.05
        star = '*';
    else
        star = 'n.s.';
    end
end