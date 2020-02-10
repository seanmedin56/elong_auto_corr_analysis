function new_traces = curate_traces(traces, max_perc, num_per_group, max_delay)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    max_keep = length(traces) * max_perc;
    while length(traces) > max_keep
        traces = traces(randperm(numel(traces)));
        num_groups = floor(length(traces) / num_per_group);
        scores = zeros(1, num_groups);
        auto_cors2nd = cell(1, num_groups);
        for i = 1:num_groups - 1
            b = (i - 1)* num_per_group + 1;
            e = i * num_per_group;
            auto_cors2nd{i} = diff(diff(auto_corr_m_calc_norm(traces(b:e), max_delay)));
        end
        auto_cors2nd{num_groups} = diff(diff(auto_corr_m_calc_norm(traces(e+1:end),max_delay)));
        
        for i = 1:num_groups
            for j = i+1:num_groups
                dif = sum((auto_cors2nd{i} - auto_cors2nd{j}).^2);
                scores(i) = scores(i) + dif;
                scores(j) = scores(j) + dif;
            end
        end
        
        [~,worst] = max(scores);
        if worst == 1
            traces = traces(num_per_group + 1:end);
        elseif worst == num_groups
            traces = traces(1:num_per_group * (num_groups - 1));
        else
            b = (worst - 1) * num_per_group;
            e = worst * num_per_group + 1;
            traces = [traces(1:b) traces(e:end)];
        end
        
    end
    new_traces = traces;
end

