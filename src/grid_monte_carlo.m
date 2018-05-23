function [best_elong,best_rise,all] = grid_monte_carlo(time_res, points_per_trace, ...
                            num_traces, num_states, trans_mat, ...
                            rna_per_sec, fluo_per_rna, ...
                            init_dist,auto_cor,elong_min,elong_max, ...
                            elong_step, rise_min, rise_max, rise_step)
%Calculates the most likely elongation time and rise time for the
%autocorrelation with all the other parameters given
%   time_res: how often the fluorescence level is measured
%   poins_per_trace: how many time points to take per synthetic trace
%   num_traces: number of synthetic traces to make
%   num_states: number of initiation rates
%   trans_mat: the matrix representing the transitions between states
%   rna_per_sec: how much rna is produced per second for each state
%   fluo_per_rna: fluorescence levels for each rna
%   init_dist: distribution of which state a trace starts at
%   auto_cor: auto correlation of the data we're trying to estimate
%   elong_min: estimated minimum possible elongation time for the real data
%   elong_max: estimated maximum possible elongation time for the real data
%   elong_step: distance between attempted elongation times
%   rise_min,rise_max,rise_step: same as three above except for the rise
%   time

    %calculates second derivative of given auto correlation
    auto_cor1st = auto_cor(2:end) - auto_cor(1:end-1);
    auto_cor2nd = auto_cor1st(2:end) - auto_cor1st(1:end-1);
    
    %things to iterate through
    elongs = elong_min:elong_step:elong_max;
    rises = rise_min:rise_step:rise_max;
    all = zeros(length(elongs),length(rises));
    
    best_err = 2^31;
    best_rise = rise_min;
    best_elong = elong_min;
    
    
    %iterates through all the options
    i = 1;
    for elong = elongs
        j = 1;
        for rise = rises
            traces = gen_data(elong,time_res,points_per_trace, ...
                            num_traces, num_states, trans_mat, ...
                            rna_per_sec, fluo_per_rna, rise, ...
                            init_dist, 0);
            corr = auto_corr_m_calc_norm(traces, length(auto_cor));
            corr1st = corr(2:end) - corr(1:end-1);
            corr2nd = corr1st(2:end) - corr1st(1:end-1);
            
            %this *might* make it more robust against different parameters
            %(will play around with it more later)
            adjust = (corr2nd * auto_cor2nd') / (corr2nd * corr2nd');
            corr2nd = corr2nd * adjust;
            
            %ignores first value
            err_arr = corr2nd(2:end) - auto_cor2nd(2:end);
            err = err_arr * err_arr';
            all(i,j) = err;
            if err < best_err
                best_err = err;
                best_rise = rise;
                best_elong = elong;
            end
            j = j + 1;
        end
        i = i + 1;
    end
end

