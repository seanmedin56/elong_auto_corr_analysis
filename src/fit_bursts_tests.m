addpath('utilities/');
% tests whether bursts in traces can be fit well for a variety of different
% conditions

results_struct = struct;
idx = 1;

% default variables
time_res = 20;
elong_time = 160;
time_per_trace = 2000; % seconds
num_traces = 20;
num_states = 3;
fluo_per_rna = 350;
on_thresh = fluo_per_rna;
noise = 0;
rna_per_sec = [.00001, .2, .4];
rise_time = 40;
init_dist = [1 0 0];
trans_mat = [-0.0115,  0.0095,  0.0000; ...
         0.0115, -0.0180,  0.0440; ...
         0.0000,  0.0085, -0.0440];

     
% analyzes how noise effects fit
noises = [0:0.5:3] * 350;
for n = noises
    points_per_trace = ceil(time_per_trace / time_res);

    [traces, arrival_traces] = gen_data(elong_time, time_res, points_per_trace, ...
    num_traces, num_states, trans_mat, rna_per_sec, ...
    fluo_per_rna, rise_time, init_dist, n);
    elong = floor(elong_time / time_res);
    rise = rise_time / time_res;
    [tot_error, spread_error, tot_fit_error, spread_fit_error, on_error, off_error] = ...
        get_burst_fit_error(traces, arrival_traces, elong, rise, ...
        time_res,fluo_per_rna,on_thresh);
    results_struct(idx).label = ['noise ' num2str(n)];
    results_struct(idx).tot_error = tot_error;
    results_struct(idx).spread_error = spread_error;
    results_struct(idx).tot_fit_error = tot_fit_error;
    results_struct(idx).spread_fit_error = spread_fit_error;
    results_struct(idx).on_error = on_error;
    results_struct(idx).off_error = off_error;
    mean_val = 0;
    for tr = traces
        mean_val = mean_val + mean(tr{1});
    end
    mean_val = mean_val / length(traces);
    results_struct(idx).mean_fluo_per = mean_val / elong;
    idx = idx + 1;
end

% analyzes how different bursting durations effects fit
mags = [.5, .1, .05 .01, .005, .001];
for mag = mags
    trans = [-mag, mag; mag, -mag];
    points_per_trace = ceil(time_per_trace / time_res);
    new_num_states = 2;
    new_dist = [0 1];
    new_rna_per_sec = [.000001 .3];
    [traces, arrival_traces] = gen_data(elong_time, time_res, points_per_trace, ...
        num_traces, new_num_states, trans, new_rna_per_sec, ...
        fluo_per_rna, rise_time, new_dist, noise);
    elong = floor(elong_time / time_res);
    rise = rise_time / time_res;
    [tot_error, spread_error, tot_fit_error, spread_fit_error,on_error,off_error] = ...
        get_burst_fit_error(traces, arrival_traces, elong, rise, ...
        time_res,fluo_per_rna,on_thresh);
    results_struct(idx).label = ['mag ' num2str(mag)];
    results_struct(idx).tot_error = tot_error;
    results_struct(idx).spread_error = spread_error;
    results_struct(idx).tot_fit_error = tot_fit_error;
    results_struct(idx).spread_fit_error = spread_fit_error;
    results_struct(idx).on_error = on_error;
    results_struct(idx).off_error = off_error;
    mean_val = 0;
    for tr = traces
        mean_val = mean_val + mean(tr{1});
    end
    mean_val = mean_val / length(traces);
    results_struct(idx).mean_fluo_per = mean_val / elong;
    idx = idx + 1;
end

% examines effect of variable falloff points on fit
fall_off_stds = [1:5];
fall_off_stds = [fall_off_stds [1:10] * 10];
for fall_off = fall_off_stds
    points_per_trace = ceil(time_per_trace / time_res);
    elong_vars = {'Gaussian', elong_time, fall_off};
    [traces, arrival_traces] = gen_data_var_len(elong_vars, time_res, points_per_trace, ...
    num_traces, num_states, trans_mat, rna_per_sec, ...
    fluo_per_rna, rise_time, init_dist, noise);
    elong = floor(elong_time / time_res);
    rise = rise_time / time_res;
    
    [tot_error, spread_error, tot_fit_error, spread_fit_error,on_error,off_error] = ...
        get_burst_fit_error(traces, arrival_traces, elong, rise, ...
        time_res,fluo_per_rna,on_thresh);
    results_struct(idx).label = ['fall off std ' num2str(fall_off)];
    results_struct(idx).tot_error = tot_error;
    results_struct(idx).spread_error = spread_error;
    results_struct(idx).tot_fit_error = tot_fit_error;
    results_struct(idx).spread_fit_error = spread_fit_error;
    results_struct(idx).on_error = on_error;
    results_struct(idx).off_error = off_error;
    mean_val = 0;
    for tr = traces
        mean_val = mean_val + mean(tr{1});
    end
    mean_val = mean_val / length(traces);
    results_struct(idx).mean_fluo_per = mean_val / elong;
    idx = idx + 1;
end
