function [trace,elongs] = gillespie_rate_spread(elong_func, time_res, points_per_trace, ...
                                 num_states, trans_mat, rna_per_sec, ...
                                 fluo_per_rna, MS2_rise_perc,init_dist,noise)
%Generates individual traces with the gillespie algorithm that have rnas
%with random elongation times picked with the elong-func function
% elong_func: function for generating how long it takes an rna polymerase 
% to traverse the gene (function of time and state)
% time_res: How often the fluorescence is measured
% points_per_trace: How many measurements are taken
% num_state: Number of distinct rates of transcription
% trans_mat: The transition flow matrix between distrinct rates of
%            transcription
% rna_per_sec: Vector of how many rnas are loaded per second for each state
% fluo_per_rna: How much fluorescence each rna produces when its fully
% loaded
% MS2_rise_perc: How far along the gene before the MS2 genes has been fully
% transcribed
% init_dist: The initial probability distribution for being in a given
% state

% the trace that will be returned
trace = zeros(1,points_per_trace);

%indexes polymerase number
idx_pol2 = 1;

%indexes naive state number
idx_naive = 1;

% preallocation size for arrays
pre_alloc = 100;

% uniformly distributed time points with resolution deltaT 
% and length seq_length
times_unif = (1:points_per_trace) * time_res;

%array for keeping track of polymerase arrivals
arrival_times = zeros(1,pre_alloc);

%keeps track of the intrinsic rates for each polymerase
int_rates = zeros(1,pre_alloc);

%keep track of the state when each polymerase came
naive_states = zeros(1,pre_alloc);

% duration of the simulated process
t_max = points_per_trace * time_res;

% first state obtained by sampling the initial state pmf
naive_states(1) = randsample(1:num_states, 1, true, init_dist);

% variable to keep track of the current reaction time
t = 0;

while (t < t_max)
    
    % determine when transition to another state occurs
    lambda = -trans_mat(naive_states(idx_naive),naive_states(idx_naive));
    dt = exprnd(1/lambda);
    t = t + dt;
    
    %determine probability of transitioning to each of the other states
    rates = trans_mat(:,naive_states(idx_naive));
    rates(naive_states(idx_naive)) = 0;
    probs = rates / lambda;
    
    %determine next state
    idx_naive = idx_naive + 1;
    naive_states(idx_naive) = randsample(1:num_states, 1, true, probs);
    
    %determine how many polymerases arrived before the transition occurred
    avg_time = 1 / rna_per_sec(naive_states(idx_naive-1));
    ddt = exprnd(avg_time);
    while ddt < dt
        arrival_times(idx_pol2) = t - dt + ddt;
        int_rates(idx_pol2) = elong_func(naive_states(end), t - dt + ddt);
        idx_pol2 = idx_pol2 + 1;
        next = exprnd(avg_time);
        ddt = ddt + next;
    end   
end

% distribution of elongation times
elongs = zeros(1,idx_pol2-1);

% uses the arrival times to generate the traces
bottleneck = 0;
crosses = [];
cross_locs = [];
for pol = 1:(idx_pol2 - 1)
    new_crosses = [];
    new_cross_locs = [];
    start = arrival_times(pol);
    start_idx_list = find(times_unif >= start);
    elongs(pol) = max(bottleneck - start, int_rates(pol));
    if isempty(start_idx_list)
        continue
    end
    first_idx = start_idx_list(1);
    end_idx_list = find(times_unif <= max(bottleneck,start + int_rates(pol)));
    
    if isempty(end_idx_list)
        continue
    end
    last_idx = end_idx_list(end);
    trace_idxes = first_idx:last_idx;
    for idx = trace_idxes
        end_time = times_unif(idx);
        if bottleneck > start + int_rates(pol) && ismember(end_time,crosses)
            rel_cross = find(crosses == end_time);
            relative_loc = min(cross_locs(rel_cross(1),(end_time - start) / int_rates(pol)));
        else
            relative_loc = (end_time - start) / int_rates(pol);
        end
        if relative_loc > MS2_rise_perc
            trace(idx) = trace(idx) + fluo_per_rna;
        else
            trace(idx) = trace(idx) + fluo_per_rna * (relative_loc / MS2_rise_perc);
        end
        new_crosses = [new_crosses end_time];
        new_cross_locs = [new_cross_locs relative_loc];
    end
    bottleneck = max(bottleneck, start + int_rates(pol));
end

% adds Gaussian noise
gauss_noise = normrnd(0,noise,1,points_per_trace);
trace = trace + gauss_noise;
trace(trace < 0) = 0;

end



