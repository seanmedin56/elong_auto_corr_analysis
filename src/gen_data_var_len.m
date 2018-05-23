function traces = gen_data_var_len(elong_vars,time_res,points_per_trace, ...
                            num_traces, num_states, trans_mat, ...
                            rna_per_sec, fluo_per_rna,MS2_rise_time, ...
                            init_dist, noise)
%This function returns traces with the features specified in the above
% elong_vars: array with two possibilities 
% if the first element is "Gaussian" the next two are the mean and standard
% deviation, otherwise the next three are the start, end, and step for a
% uniform distribution
% time_res: How often the fluorescence is measured
% points_per_trace: How many measurements are taken
% num_traces: How many traces are taken
% num_state: Number of distinct rates of transcription
% trans_mat: The transition flow matrix between distrinct rates of
%            transcription
% rna_per_sec: Vector of how many rnas are loaded per second for each state
% fluo_per_rna: How much fluorescence each rna produces when its fully
% loaded
% MS2_rise_time: How long it takes for a polymerase to fully fluoresce
% init_dist: The initial probability distribution for being in a given
% state

addpath('utilities/');

traces = cell([1,num_traces]);
if strcmp(elong_vars{1},'Gaussian')
    elong_func = @(state,t) normrnd(elong_vars{2},elong_vars{3});
else
    elong_func = @(state,t) randsample(elong_vars{2}:elong_vars{4}:elong_vars{3},1);
end

for i = 1:num_traces
    traces{i} = gillespie_var_len(elong_func, time_res, points_per_trace, ...
                              num_states, trans_mat, rna_per_sec, ...
                              fluo_per_rna, MS2_rise_time,init_dist,noise);
end
    
end




