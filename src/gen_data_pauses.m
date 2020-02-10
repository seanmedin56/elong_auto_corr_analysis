function traces = gen_data_pauses(elong_time,pause_locs,pause_means, ...
                            pause_sds,time_res,points_per_trace, ...
                            num_traces, num_states, trans_mat, ...
                            rna_per_sec, fluo_per_rna,MS2_rise_time, ...
                            init_dist, noise)
%This function returns traces with the features specified in the above
% elong_time: how long it would take to transcribe the gene without pauses
% pause_locs: array of where the pauses are
% pause_means: array of the mean times of how long each pause is
% pause_sds: array of the standard deviations of how long each pause is
% time_res: How often the fluorescence is measured
% points_per_trace: How many measurements are taken
% num_traces: How many traces are taken
% num_state: Number of distinct rates of transcription
% trans_mat: The transition flow matrix between distrinct rates of
%            transcription
% rna_per_sec: Vector of how many rnas are loaded per second for each state
% fluo_per_rna: How much fluorescence each rna produces when its fully
% loaded
% MS2_rise_time: What percentage along the gene must be traversed before
% all the MS2 has been transcribed
% init_dist: The initial probability distribution for being in a given
% state

addpath('utilities/');

traces = cell([1,num_traces]);
elongs = cell([1,num_traces]);
pause_funcs = cell(1,length(pause_locs));
[pause_locs,idxes] = sort(pause_locs);
pause_means = pause_means(idxes);
pause_sds = pause_sds(idxes);
for i = 1:length(pause_locs)
    pause_funcs{i} = @(naive_state,time) max(0,normrnd(pause_means(i),pause_sds(i)));
end
for i = 1:num_traces
    [trace,elong] = gillespie_pauses(elong_time, pause_funcs, pause_locs, ...
                              time_res, points_per_trace, ...
                              num_states, trans_mat, rna_per_sec, ...
                              fluo_per_rna, MS2_rise_time,init_dist,noise);
    traces{i} = trace;
    elongs{i} = elong;
end

% uncomment below to also plot spread
plot_elong_spread(elongs,20);
    
end


