function [stds,stds1,stds2] = corr_bootstraps(trace1, trace2, max_delay, num_times, type)

% takes random selections of traces and calculates the standard deviation
% of the autocorrelation and its derivatives
% can be done for cross correlation (trace1 != trace2) or for
% autocorrelation (trace1 == trace2)
%   max_delay: number of time delay points to take in the auto correlation
%   num_times: number of times to sample the traces
%   type: 'r' = raw moment, anything else = central moment

    vals = cell([1 max_delay]);
    first = cell([1 max_delay-1]);
    second = cell([1 max_delay-2]);
    
    %iterates through random subsamples of traces
    for i = 1:num_times
        sample_idx = randi([1 length(trace1)], 1, length(trace1));
        sample1 = cell([1 length(trace1)]);
        sample2 = cell([1 length(trace2)]);
        for j = 1:length(sample1)
            sample1{j} = trace1{sample_idx(j)};
            sample2{j} = trace2{sample_idx(j)};
        end

        if type == 'r'
            corr = cross_corr_r_calc(sample1, sample2, max_delay);
        else
            corr = cross_corr_m_calc(sample1, sample2, max_delay);
        end

        for j = 1:length(corr)
            vals{j}(i) = corr(j);
        end
        
        corr_1st_deriv = corr(2:max_delay) - corr(1:max_delay-1);
        for j = 1:length(corr_1st_deriv)
            first{j}(i) = corr_1st_deriv(j);
        end
        
        corr_2nd_deriv = corr_1st_deriv(2:max_delay-1) - corr_1st_deriv(1:max_delay-2);
        for j = 1:length(corr_2nd_deriv)
            second{j}(i) = corr_2nd_deriv(j);
        end
        
    end
    stds = zeros([1 length(vals)]);
    for i = 1:length(vals)
        stds(i) = std(vals{i});
    end
    
    stds1 = zeros([1 length(first)]);
    for i = 1:length(first)
        stds1(i) = std(first{i});
    end
    
    stds2 = zeros([1 length(second)]);
    for i = 1:length(second)
        stds2(i) = std(second{i});
    end
end