function cross_corr_m = cross_corr_m_calc(trace1, trace2, max_delay)

% returns the average crosscovariance of the traces

% ----------------------calculates the global means-------------------------

    global_mean1 = 0;
    count = 0;
    for i = 1:length(trace1)
        global_mean1 = global_mean1 + sum(trace1{i});
        count = count + length(trace1{i});
    end
    global_mean1 = global_mean1 / count;
    
    global_mean2 = 0;
    count = 0;
    for i = 1:length(trace2)
        global_mean2 = global_mean2 + sum(trace2{i});
        count = count + length(trace2{i});
    end
    global_mean2 = global_mean2 / count;
    
%--------------calculates crosscovariance using cross_corr_r_calc----------

    new_traces1 = cell([1 length(trace1)]);
    for i = 1:length(trace1)
        new_traces1{i} = trace1{i} - global_mean1;
    end
    
    new_traces2 = cell([1 length(trace2)]);
    for i = 1:length(trace2)
        new_traces2{i} = trace2{i} - global_mean2;
    end
    cross_corr_m = cross_corr_r_calc(new_traces1, new_traces2, max_delay);

end