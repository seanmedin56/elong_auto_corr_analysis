% This script generates auto correlations for different elongation dynamics
% with different gene length and creates a graph comparing different select
% features (specified before any calculations are done)

addpath('utilities/');

%functions for features to be analyzed

features = {@extract_2nd_deriv_peak};
feature_names = {'Peak'};
num_cases = 4;
case_names = {'constant rate' , 'rate spread', 'var len', 'pausing'};
feature_plots = cell(num_cases, length(features));

%defines variables universal to all scripts

time_res = 10; 

points_per_trace = 200;

num_traces = 250;

num_states = 2;

trans_mat = [-.01,.01;.01,-.01];

rna_per_sec = [.00001,.2];

fluo_per_rna = 350;

init_dist = [.5,.5];

noise = 300;

construct_lengths = .95:.4:9.4; %in kb
base_length = 4.7;

base_time = 160; % 160 seconds to transribe 4.7kb

alpha_perc = .25; % percent of 4.7kb that is ms2 loops

max_delay = 100;

% generates data and extracts statistics for each case

case_num = 1;

%constant elongation time
elong_times = (base_time / base_length) * construct_lengths;
rise_time = alpha_perc * base_time;

for i = 1:length(elong_times)
    elong = elong_times(i);
    traces = gen_data(elong,time_res,points_per_trace,num_traces,num_states, ...
        trans_mat, rna_per_sec, fluo_per_rna, rise_time, init_dist,noise);
    
    corr = auto_corr_m_calc_norm(traces, max_delay);

    for f = 1:length(features)
        funct = features{f};
        feature = funct(corr);
        feature_plots{case_num, f}(i) = feature;
    end
     
end
case_num = case_num + 1;
% diff rates elongation time gaussian

base_width = 20;
elong_widths = (base_width / base_length) * construct_lengths;
alpha_perces = (alpha_perc / base_length) * construct_lengths;

for i = 1:length(elong_times)
    elong_var = {'Gaussian', elong_times(i), elong_widths(i)};
    traces = gen_data_rate_spread(elong_var, time_res, points_per_trace, ...
        num_traces, num_states, trans_mat, rna_per_sec, fluo_per_rna, ...
        alpha_perces(i), init_dist, noise);
    
    corr = auto_corr_m_calc_norm(traces, max_delay);
    
    for f = 1:length(features)
        funct = features{f};
        feature = funct(corr);
        feature_plots{case_num, f}(i) = feature;
    end
    
end
case_num = case_num + 1;

% different falloff points (no bottleneck)

for i = 1:length(elong_times)
    elong_var = {'Gaussian', elong_times(i), elong_widths(i)};
    traces = gen_data_var_len(elong_var, time_res, points_per_trace, ...
        num_traces, num_states, trans_mat, rna_per_sec, fluo_per_rna, ...
        rise_time, init_dist, noise);
    
    corr = auto_corr_m_calc_norm(traces, max_delay);
    
    for f = 1:length(features)
        funct = features{f};
        feature = funct(corr);
        feature_plots{case_num, f}(i) = feature;
    end
    
end
case_num = case_num + 1;

%sequence dependent pausing
pause_loc = .9;
mean_pause = 4;
sd_pause = 2;
for i = 1:length(elong_times)
    traces = gen_data_pauses(elong_times(i), pause_loc, mean_pause, ...
        sd_pause,time_res, points_per_trace, num_traces, num_states, ...
        trans_mat, rna_per_sec, fluo_per_rna, alpha_perc, init_dist, noise);
    
    corr = auto_corr_m_calc_norm(traces, max_delay);
    
    for f = 1:length(features)
        funct = features{f};
        feature = funct(corr);
        feature_plots{case_num, f}(i) = feature;
    end
    
end

case_num = case_num + 1;

for f = 1:length(features)
    fig = figure();
    
    for i = 1:num_cases
        plot(construct_lengths, feature_plots{i,f})
        hold on
    end
    xlabel('Size of gene (kb)');
    ylabel(strcat(feature_names{f}, ' (', num2str(time_res), ' second delays)'));
    title(strcat('Comparing ', feature_names{f}));
    legend(case_names);
end

