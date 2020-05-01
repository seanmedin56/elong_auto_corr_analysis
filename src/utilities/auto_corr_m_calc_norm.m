% returns the average autocovariance of the traces
function auto_corr_m = auto_corr_m_calc_norm(traces, max_delay)



% ----------------------calculates the global mean-------------------------

    global_mean = 0;
    count = 0;
    for i = 1:length(traces)
        global_mean = global_mean + nansum(traces{i});
        count = count + sum(~isnan(traces{i}));
    end
    global_mean = global_mean / count;
    
%--------------calculates autocovariance using auto_corr_r_calc------------
    new_traces = cell([1 length(traces)]);
    for i = 1:length(traces)
        new_traces{i} = traces{i} - global_mean;
    end
    auto_corr_m = auto_corr_r_calc_norm(new_traces, max_delay);

end


