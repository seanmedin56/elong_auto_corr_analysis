function std_derivs = corr_bootstraps(trace1, trace2, max_delay, num_times, ...
    num_derivs, type)

% takes random selections of traces and calculates the standard deviation
% of the autocorrelation and its derivatives
% can be done for cross correlation (trace1 != trace2) or for
% autocorrelation (trace1 == trace2)
%   max_delay: number of time delay points to take in the auto correlation
%   num_times: number of times to sample the traces
%   type: 'r' = raw moment, anything else = central moment

    vals = zeros(num_times, max_delay);
    std_derivs = cell(1, num_derivs + 1);
    
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
        elseif type == 'm'
            corr = cross_corr_m_calc_norm(sample1, sample2, max_delay);
        elseif type == "r1"
            corr = cross_corr_r_calc2(sample1, sample2, max_delay);
        elseif type == "m1"
            corr = cross_corr_m_calc2(sample1, sample2, max_delay);
        elseif type == "r2"
            corr = cross_corr_r_calc3(sample1, sample2, max_delay);
        elseif type == "m2"
            corr = cross_corr_m_calc3(sample1, sample2, max_delay);
        elseif type == "r3"
            corr = cross_corr_r_calc4(sample1, sample2, max_delay);
        elseif type == "r4"
            corr = cross_corr_m_calc4(sample1, sample2, max_delay);
        elseif type == "c1"
            corr = fin_corr1(sample1, sample2, max_delay);
        elseif type == "c2"
            corr = fin_corr2(sample1, sample2, max_delay);
        elseif type == "c3"
            corr = fin_corr3(sample1, sample2, max_delay);
        else
            corr = fin_corr4(sample1, sample2, max_delay);
        end

        vals(i,:) = corr;
    end
    for deriv = 0:num_derivs
        new_vals = vals;
        for idx = 1:deriv
            new_vals = new_vals(:,2:end) - new_vals(:,1:end-1);
        end
        std_derivs{deriv + 1} = std(new_vals, 0, 1);
    end
end