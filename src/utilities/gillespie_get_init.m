function [trace] = gillespie_get_init(time_res, points_per_trace, ...
                                 num_states, trans_mat, rna_per_sec, ...
                                 init_dist)
%Generates individual traces for possible initiation times with the 
% gillespie algorithm
% elong_time: how long it takes an rna polymerase to traverse the gene
% time_res: How often the fluorescence is measured
% points_per_trace: How many measurements are taken
% num_state: Number of distinct rates of transcription
% trans_mat: The transition flow matrix between distrinct rates of
%            transcription
% rna_per_sec: Vector of how many rnas are loaded per second for each state
% fluo_per_rna: How much fluorescence each rna produces when its fully
% loaded
% MS2_rise_time: How long it takes for a polymerase to fully fluoresce
% init_dist: The initial probability distribution for being in a given
% state

% the trace that will be returned
trace = zeros(1,points_per_trace);

% uniformly distributed time points with resolution deltaT 
% and length seq_length
times_unif = (1:points_per_trace) * time_res;

%keep track of the state when each polymerase came
naive_states = [];

%keep track of when each state was added
arrival_times = [0];

% duration of the simulated process
t_max = points_per_trace * time_res;

% first state obtained by sampling the initial state pmf
naive_states(1) = randsample(1:num_states, 1, true, init_dist);

% variable to keep track of the current reaction time
t = 0;

while (t < t_max)
    
    % determine when transition to another state occurs
    lambda = -trans_mat(naive_states(end),naive_states(end));
    dt = exprnd(1/lambda);
    t = t + dt;
    arrival_times = [arrival_times t + dt];
    
    %determine probability of transitioning to each of the other states
    rates = trans_mat(:,naive_states(end));
    rates(naive_states(end)) = 0;
    probs = rates / lambda;
    
    %determine next state
    naive_states = [naive_states randsample(1:num_states, 1, true, probs)];
    
    
end

% uses the arrival times to generate the traces
for k = 1:points_per_trace
    t_end = times_unif(k);

    ind_start_before = find(arrival_times > t_end - 20);
    
    if isempty(ind_start_before)
        i_start = length(arrival_times) - 1;
    else
        i_start = ind_start_before(1) - 1;
    end

    ind_end_before = find(arrival_times <= t_end);

    i_end = ind_end_before(end);

    times_window = arrival_times(i_start:i_end);
    states_window = naive_states(i_start:i_end);
    
    if length(times_window) > 1
        for i = 2:(length(times_window) - 1)

            t2 = times_window(i+1) - times_window(i);
            trace(k) = trace(k) + t2 * states_window(i);
        end
    end
    trace(k) = trace(k) + (t_end - times_window(end)) * states_window(end);
end


end

