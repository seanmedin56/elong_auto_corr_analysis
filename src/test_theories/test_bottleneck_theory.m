addpath('../utilities/');
addpath('../');

num_states = 3;
time_res = 10; % seconds
points_per_trace = 240; % 40 minutes
max_delay = 100;
elongs = 40:25:200;
elong_prop_stds = zeros(1,7);
elong_prop_stds(1) = 0.1;
for i = 2:length(elong_prop_stds)
    elong_prop_stds(i) = elong_prop_stds(i-1) * 10^(1 / 6);
end
rise_times = 0:0.25:1;
trans_mat = [-0.0115,  0.0095,  0.0000; ...
             0.0115, -0.0180,  0.0440; ...
             0.0000,  0.0085, -0.0440];
rna_per_sec = [.0001, 0.2, 0.4];
init_dist = [0, 1, 0];
noise = 0;
num_traces = 300;

dists = struct;
idx = 1;
for T = elongs
    for a0 = rise_times
        for Tsig = elong_prop_stds
            elong_vars = {'Gaussian', T, T * Tsig};
            [traces,mean_elong,elong_std] = gen_data_rate_spread(elong_vars, ...
                time_res, points_per_trace, num_traces, num_states, trans_mat,...
                rna_per_sec, 1000, a0, init_dist, noise);
            auto_cor = auto_corr_m_calc_norm(traces, max_delay);
            deriv1 = diff(auto_cor);
            deriv2 = diff(deriv1);
            deriv3 = diff(deriv2);
            [~,peak2] = max(deriv2);
            idxes = find(deriv2 < 0);
            first_zero = idxes(find(idxes > peak2, 1, 'first'));
            [~,I] = min(deriv3(peak2:first_zero));
            dists(idx).T = T;
            dists(idx).Tsig = Tsig;
            dists(idx).mean_elong = mean_elong;
            dists(idx).a0 = a0;
            dists(idx).std_elong = elong_std;
            dists(idx).min = I + peak2 - 1;
            idx = idx + 1;
        end
    end
end

% plot to check results
Tsigs = [dists.Tsig];
errs = [dists.min] - [dists.mean_elong] / time_res;
figure();
scatter(Tsigs, errs);
