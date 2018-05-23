function traces = gen_data(elong_time,time_res,points_per_trace, ...
                            num_traces, num_states, trans_mat, ...
                            rna_per_sec, fluo_per_rna,MS2_rise_time, ...
                            init_dist, noise)
%This function returns traces with the features specified in the above
% elong_time: how long it takes an rna polymerase to traverse the gene
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
for i = 1:num_traces
    traces{i} = gillespie_gen(elong_time, time_res, points_per_trace, ...
                              num_states, trans_mat, rna_per_sec, ...
                              fluo_per_rna, MS2_rise_time,init_dist,noise);
end
    
end

