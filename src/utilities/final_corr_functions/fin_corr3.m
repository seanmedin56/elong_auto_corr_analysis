function auto = fin_corr3(trace1, trace2, max_delay)

% calculates central moment
% does normalize individual correlations
% does not weigh individual correlations by number of points per trace

% ----------------------calculates the global means-------------------------

    global_mean1 = 0;
    count = 0;
    for i = 1:length(trace1)
        global_mean1 = global_mean1 + nansum(trace1{i});
        count = count + sum(~isnan(trace1{i}));
    end
    global_mean1 = global_mean1 / count;
    
    global_mean2 = 0;
    count = 0;
    for i = 1:length(trace2)
        global_mean2 = global_mean2 + nansum(trace2{i});
        count = count + sum(~isnan(trace2{i}));
    end
    global_mean2 = global_mean2 / count;
    
%--------------subtracts means----------

    new_traces1 = cell([1 length(trace1)]);
    for i = 1:length(trace1)
        new_traces1{i} = trace1{i} - global_mean1;
    end
    
    new_traces2 = cell([1 length(trace2)]);
    for i = 1:length(trace2)
        new_traces2{i} = trace2{i} - global_mean2;
    end
    
    trace1 = new_traces1;
    trace2 = new_traces2;
   % ----------calculates individual correaltions (with weights)----------
   
    corrs = cell([1 length(trace1)]);
    counts = zeros([1 max_delay]);
    for i = 1:length(trace1)
        limit = min([max_delay, length(trace1{i}), length(trace2{i})]);
        len = min(length(trace1{i}), length(trace2{i}));
        corr = zeros([1 limit]);
        num_points = zeros([1 limit]);
        for j = 1:limit
            temp_mult = trace1{i}(1:len - j + 1) .* trace2{i}(j:len);
            corr(j) = nansum(temp_mult);
            num_points(j) = sum(~isnan(temp_mult));
        end
        counts = counts + num_points;
        corr = corr ./ num_points;
        corrs{i} = corr / corr(1);
    end
    
% -----------------combines correaltions together------------------------    

    cross_corr_r = zeros([1 max_delay]);
    for i = 1:length(corrs)
        cross_corr_r(1:length(corrs{i})) = cross_corr_r(1:length(corrs{i})) ...
            + corrs{i};
    end
    %cross_corr_r = cross_corr_r / counts;
    auto = cross_corr_r / cross_corr_r(1); 

end