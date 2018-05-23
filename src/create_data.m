% code for generating traces AND running inference on those traces
% liberally copy and pasted from mhmm repository

addpath('utilities/');

% create directories for outputing data and the results
dirs = {'../out/synth_dat/', '../dat/synth_dat/', ...
         '../fig/synth_dat/'};
     
for i = 1:numel(dirs)
    if (exist(dirs{i}, 'dir') ~= 7)
        mkdir(dirs{i});
    end
end

% number of promoter states
K = 3;

% number of groups
nGR = 1;

% time resolution
deltaT = {12.1, 20, 10};

% elongation time
t_elong = 12.1 * 4;

% system memory
w = 4;

% time it takes to transcribe the MS2 loops [sec]
t_MS2 = 32.5;

% alpha: length of the MS2 loop in time steps
alpha = cell(1,nGR);
for n = 1:nGR
    alpha{n} = t_MS2 / deltaT{n};
end

% transition rates [sec^-1]
R = cell([1, nGR]);
R{1} = [-0.0115,  0.0095,  0.0000; ...
         0.0115, -0.0180,  0.0440; ...
         0.0000,  0.0085, -0.0440];
     
R{2} = R{1};
R{3} = R{1};

% transition probabilities
A = cell([1, nGR]);
for i = 1:nGR
    A{i} = rate_to_prob(R{i}, deltaT{i});
end

% initial state pmf
pi0 = cell([1, nGR]);
pi0{1} = [0.05, 0.05, 0.9];
pi0{2} = pi0{1};
pi0{3} = pi0{1};

% emission rates [a.u. / sec]
r_emission = cell([1, nGR]);
r_emission{1} = [2, 65, 130];
r_emission{2} = r_emission{1};
r_emission{3} = r_emission{2};

% background noise [a.u.]
noise = cell([1, nGR]);
for i = 1:nGR
    noise{i} = 1000;
end

% state conversion (which states are active in the two loci)
conv = cell([1, nGR]);
conv{1} = [1, 0; 0, 2; 0, 3];
conv{2} = conv{1};
conv{3} = conv{1};

% fluorescence per rna [a.u. / rna]
fluo_per_rna = 350;

% time it takes for another polymerase to load [sec]
promoter_pause = 5;
loading_rate = 1;
wait_time = promoter_pause + loading_rate;

% mandatory wait time between arrivals
flor = 0;

% total number of time points in a pooled data set
n_points_total = [ 50000, 3000, 3000];

% average number of points per trace;
seq_length = [50, 50, 50];

% number of traces in each AP position
n_traces = round(n_points_total./seq_length);

% ---------------- Save parameters into a structure -------------------

% structure to store all synthetic parameter values
synthetic_parameters = struct;
synthetic_parameters.K = K;
synthetic_parameters.w = w;
synthetic_parameters.t_MS2 = t_MS2;
synthetic_parameters.deltaT = deltaT;
synthetic_parameters.alpha = alpha;

synthetic_parameters.R = R;
synthetic_parameters.A = A;
synthetic_parameters.pi0 = pi0;
synthetic_parameters.noise = noise;
synthetic_parameters.r_emission = r_emission;

synthetic_parameters.n_points_total = n_points_total;
synthetic_parameters.seq_length = seq_length;
synthetic_parameters.n_traces = n_traces;

% structure to store synthetic data sets
data = struct;

% index variable for the current generated trace
set_ind = 1;
n_sets = 1;

% maximum number of EM iterations
n_steps_max = 1000;

% tolerance parameter for inference convergence
eps = 10^(-4);

% matrix to store the time resolution and set number pairs
set_count_mat = zeros(n_sets, nGR);

for i = 1:nGR
    for set = 1:n_sets
        
        fluo_data = cell([1, n_traces(i)]);
        old_fluo_data = cell([1, n_traces(i)]);
        state_data = cell([1, n_traces(i)]);
        poll2_data = cell([1, n_traces(i)]);
        for tr = 1:n_traces(i)
            fluo_gill = synthetic_rate_gillespie_two(seq_length(i), alpha{i}, ...
                K, w, R{i}, deltaT{i}, r_emission{i}, noise{i}, pi0{i}, ...
                fluo_per_rna, flor, conv{i});
            fluo_data{tr} = fluo_gill.fluo_MS2;
        end

        data(set_ind).fluo_data = fluo_data;
        data(set_ind).deltaT = deltaT{i};
        data(set_ind).deltaT_ind = i;
        data(set_ind).n_traces = n_traces(i);
        data(set_ind).alpha = alpha{i};
        data(set_ind).set = set;
        
        set_count_mat(set, i) = set_ind;
        set_ind = set_ind + 1;
        
        
    end
end

%runs inference if run_inference is set to 1
run_inference = 0;

if run_inference == 1
    % structure array to store the analysis data
    outputs = struct;

    pool = parpool(str2num(getenv('SLURM_CPUS_ON_NODE')));
    for i = 1:nGR
        disp('very early');
        % set and time resolution indices

        % extract the synthetic data corresponding to the [set, i] pair
        fluo_data = data(i).fluo_data;
        disp('early');
        logL_max = -Inf;

        % random initialization of model parameters
        param_init2 = initialize_random (K, w, fluo_data);

        % approximate inference assuming iid data
        local_iid_out = local_em_iid (fluo_data, param_init2.v, ...
                param_init2.noise, K, w, alpha{i}, 1000, 1e-4);
        noise_iid = 1/sqrt(exp(local_iid_out.lambda_log));
        v_iid = exp(local_iid_out.v_logs);
        disp('kind of early');
        local_struct = struct;
        % localEM runs
        
        n_localEM = 50;
        parfor i_local = 1:n_localEM
            disp('not that early');
            % random initialization of model parameters
            param_init = initialize_random_with_priors(K, noise_iid, v_iid);

            pi0_log_init = log(param_init.pi0);
            A_log_init = log(param_init.A);
            v_init = param_init.v;
            noise_init = param_init.noise;

            % localEM call
            local_out = local_em_MS2 (fluo_data, ...
                v_init, noise_init, pi0_log_init', A_log_init, K, w, ...
                alpha{i}, n_steps_max, eps);

            
            local_struct(i_local).logL = local_out.logL;
            local_struct(i_local).A_log_inf = local_out.A_log;
            local_struct(i_local).A_inf = exp(local_out.A_log);
            local_struct(i_local).R_inf = prob_to_rate(local_out.A_log, deltaT{i});
            v_inf = exp(local_out.v_logs).*local_out.v_signs;
            local_struct(i_local).v_inf = v_inf;
            local_struct(i_local).r_inf = v_inf / deltaT{i};
            lambda_inf = exp(local_out.lambda_log);
            local_struct(i_local).noise_inf = 1/sqrt(lambda_inf);
            local_struct(i_local).pi0_log_inf = local_out.pi0_log;
                
        end
        disp('first marker');
        [logL, max_index] = max([local_struct.logL]);
        outputs(i).pi0 =exp(local_struct(max_index).pi0_log_inf);
        outputs(i).pi0_log = local_struct(max_index).pi0_log_inf;

        outputs(i).v = local_struct(max_index).v_inf(:);
        outputs(i).r = local_struct(max_index).r_inf(:);

        outputs(i).noise = local_struct(max_index).noise_inf;

        outputs(i).A = local_struct(max_index).A_inf(:);
        outputs(i).A_mat = local_struct(max_index).A_inf;
        outputs(i).A_log = local_struct(max_index).A_log_inf;

        outputs(i).R = local_struct(max_index).R_inf(:);
        outputs(i).R_mat = local_struct(max_index).R_inf;

        outputs(i).deltaT = deltaT{i};

        disp('second marker');
    end
    delete(pool);
end

% extract the current date in a string format
date_str = ['dT' num2str(deltaT{1}) 'w' num2str(w) 'n' num2str((noise{1} / 1000)) 'kal' ...
    num2str(alpha{1}) 'poisson'];

if exist(['../dat/synth_dat/' date_str '_data.mat'], 'file')
    warning('DO NOT SAVE OVER PAST DATA');
    return
end

% save the generated data into a '.mat' file
save(['../dat/synth_dat/' date_str '_data.mat'], 'data');
save(['../dat/synth_dat/' date_str '_params.mat'], 'synthetic_parameters');

% save the outputs into a '.mat' file (if inference is run)
if run_inference == 1
    save(['../out/synth_dat/' date_str '_outputs.mat'], 'outputs');
end