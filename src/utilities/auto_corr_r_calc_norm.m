% returns the average raw moment autocorrelation of the traces
function auto_corr_r = auto_corr_r_calc_norm(traces, max_delay)


% ----------calculates individual correaltions (with weights)--------------

    corrs = cell([1 length(traces)]);
    for i = 1:length(traces)
        limit = min(max_delay, length(traces{i}));
        len = length(traces{i});
        corr = zeros([1 limit]);
        for j = 1:limit
            corr(j) = traces{i}(1:len - j + 1) * traces{i}(j:len)';
        end
        corr = corr ./ (len:-1:(len - limit + 1));
        corrs{i} = corr / corr(1);
    end
    
% -----------------combines correaltions together------------------------    

    auto_corr_r = zeros([1 max_delay]);
    for i = 1:length(corrs)
        auto_corr_r(1:length(corrs{i})) = auto_corr_r(1:length(corrs{i})) ...
            + corrs{i};
    end
    auto_corr_r = auto_corr_r / auto_corr_r(1);

end


