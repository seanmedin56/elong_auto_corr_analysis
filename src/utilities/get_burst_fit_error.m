function [tot_error,spread_error, tot_fit_error, spread_fit_error, ...
    on_error, off_error] = ...
    get_burst_fit_error(traces, arrival_traces,elong_time,rise_time, ...
    time_res, fluo_per_rna, on_thresh)
%GET_BURST_FIT_ERROR Summary of this function goes here
%   Detailed explanation goes here
    tot_error = 0;
    spread_error = 0;
    tot_fit_error = 0;
    spread_fit_error = 0;
    on_error = 0;
    off_error = 0;
    for j = 1:length(traces)
        trace = traces{j};
        arrival_trace = arrival_traces{j};
        [fitted_trace, x] = fit_bursts(trace,elong_time,rise_time);
        
        % calculates real bursts
        real_instants = zeros(1, length(trace));
        for i = 1:length(real_instants)
            t_start = (i - 1) * time_res;
            t_end = i * time_res;
            idxes = find(arrival_trace > t_start & arrival_trace < t_end);
            real_instants(i) = length(idxes) * fluo_per_rna;
        end
        
        % calculates differences and adds them to the errors
        tot_fit_error = tot_fit_error + mean((fitted_trace(1:end-elong_time)' ...
            - trace(1:end-elong_time)).^2);
        spread_fit_error = spread_fit_error + std((fitted_trace(1:end-elong_time)'...
            - trace(1:end-elong_time)).^2);
        tot_error = tot_error + mean((real_instants(1:end-elong_time) ...
            - x(1:end-elong_time)').^2);
        spread_error = spread_error + std((real_instants(1:end-elong_time) ...
            - x(1:end-elong_time)').^2);
        
        % calculates on/off errors
        real_on_idxes = find(real_instants(1:end-elong_time) > 0);
        real_off_idxes = find(real_instants(1:end-elong_time) == 0);
        fit_on_idxes = find(x(1:end-elong_time) > on_thresh)';
        fit_off_idxes = find(x(1:end-elong_time) <= on_thresh)';
        on_error = on_error + (1 - ...
            sum(ismember(fit_on_idxes, real_on_idxes)) / length(real_on_idxes));
        off_error = off_error + (1 - ...
            sum(ismember(fit_off_idxes,real_off_idxes)) / length(real_off_idxes));
        
    end
    tot_fit_error = sqrt(tot_fit_error / length(traces));
    spread_fit_error = sqrt(spread_fit_error / length(traces));
    tot_error = sqrt(tot_error / length(traces));
    spread_error = sqrt(spread_error / length(traces));
    on_error = on_error / length(traces);
    off_error = off_error / length(traces);
end
