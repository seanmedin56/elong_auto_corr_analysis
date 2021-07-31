% returns the average raw moment autocorrelation of the traces
function auto_corr_r = auto_corr_r_calc_norm(traces, max_delay)


% ----------calculates individual correaltions (with weights)--------------

    corrs = cell([1 length(traces)]);
    trace_lens = cell([1 length(traces)]);
    counts = zeros([1 max_delay]);
    for i = 1:length(traces)
        limit = min(max_delay, length(traces{i}));
        len = length(traces{i});
        corr = zeros([1 limit]);
        num_points = zeros([1 limit]);
        for j = 1:limit
            temp_mult = traces{i}(1:len - j + 1) .* traces{i}(j:len);
            corr(j) = nansum(temp_mult);
            num_points(j) = sum(~isnan(temp_mult));
        end
        counts = counts + num_points;
        corr = corr ./ num_points;
        corrs{i} = corr / corr(1);
        trace_lens{i} = len;
    end
    
% -----------------combines correaltions together------------------------    

    auto_corr_r = zeros([1 max_delay]);
    for i = 1:length(corrs)
        auto_corr_r(1:length(corrs{i})) = auto_corr_r(1:length(corrs{i})) ...
            + corrs{i} * trace_lens{i};
    end
    %auto_corr_r = auto_corr_r ./ counts;
    auto_corr_r = auto_corr_r / auto_corr_r(1);

end


