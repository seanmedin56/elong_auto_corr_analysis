function traces = get_correction(time_res,points_per_trace, ...
                            num_traces, num_states, trans_mat, ...
                            rna_per_sec, init_dist)
%This function returns traces with the features specified in the above
% time_res: How often the fluorescence is measured
% points_per_trace: How many measurements are taken
% num_traces: How many traces are taken
% num_state: Number of distinct rates of transcription
% trans_mat: The transition flow matrix between distrinct rates of
%            transcription
% rna_per_sec: Vector of how many rnas are loaded per second for each state
% init_dist: The initial probability distribution for being in a given
% state
addpath('utilities/');

traces = cell([1,num_traces]);
for i = 1:num_traces
    traces{i} = gillespie_get_init(time_res, points_per_trace, ...
                              num_states, trans_mat, rna_per_sec, ...
                              init_dist);
end
    
end

