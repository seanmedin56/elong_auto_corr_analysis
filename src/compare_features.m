% This script generates auto correlations for different elongation dynamics
% with different gene length and creates a graph comparing different select
% features (specified before any calculations are done)

addpath('utilities/');

%functions for features to be analyzed

features = {@extract_width, @extract_center, @extract_2nd_deriv_peak};
feature_names = {'Widths', 'Centers', 'Peaks'};
num_cases = 1;
case_names = {'rate spread', 'rate spread fit'};

%case_names = {'constant rate', 'rate spread', ...
%     'var len', 'pausing'};

feature_plots = cell(num_cases, length(features));

% stuff for dealing with errorbars
use_errorbars = true;
feature_errors = cell(num_cases, length(features));
num_bootstraps = 10;

%defines variables universal to all scripts

time_res = 5; 

points_per_trace = 400;

num_traces = 250;

num_states = 2;

trans_mat = [-.01,.01;.01,-.01];

rna_per_sec = [.00001,.2];

fluo_per_rna = 350;

init_dist = [.5,.5];

noise = 300;

%construct_lengths = [1.9,2.65,3.4,4.7]; %in kb
construct_lengths = 1.7:.5:15.7;
base_length = 4.7;

base_time = 160; % 160 seconds to transribe 4.7kb

alpha_perc = .25; % percent of 4.7kb that is ms2 loops

max_delay = 200;

% generates data and extracts statistics for each case

case_num = 1;

% %constant elongation time
elong_times = (base_time / base_length) * construct_lengths;
rise_time = alpha_perc * base_time;
% 
% for i = 1:length(elong_times)
%     elong = elong_times(i);
%     traces = gen_data(elong,time_res,points_per_trace,num_traces,num_states, ...
%         trans_mat, rna_per_sec, fluo_per_rna, rise_time, init_dist,noise);
%     
%     [fps, fes] = extract_features(features,traces, ...
%         use_errorbars,max_delay,num_bootstraps);  
%     
%     for f = 1:length(features)
%         feature_plots{case_num, f}(i) = fps(f);
%         feature_errors{case_num, f}(i) = fes(f);
%     end
%      
% end
% case_num = case_num + 1;
% diff rates elongation time gaussian

base_width = 80;
elong_widths = (base_width / base_length) * construct_lengths;
alpha_perces = (alpha_perc / base_length) * construct_lengths;

for i = 1:length(elong_times)
    elong_var = {'Gaussian', elong_times(i), elong_widths(i)};
    traces = gen_data_rate_spread(elong_var, time_res, points_per_trace, ...
        num_traces, num_states, trans_mat, rna_per_sec, fluo_per_rna, ...
        alpha_perces(i), init_dist, noise);
    
    [fps, fes] = extract_features(features,traces, ...
        use_errorbars,max_delay,num_bootstraps);  
    
    for f = 1:length(features)
        feature_plots{case_num, f}(i) = fps(f);
        feature_errors{case_num, f}(i) = fes(f);
    end
    
end
case_num = case_num + 1;

% different falloff points (no bottleneck)

% for i = 1:length(elong_times)
%     elong_var = {'Gaussian', elong_times(i), elong_widths(i)};
%     traces = gen_data_var_len(elong_var, time_res, points_per_trace, ...
%         num_traces, num_states, trans_mat, rna_per_sec, fluo_per_rna, ...
%         rise_time, init_dist, noise);
%     
%     [fps, fes] = extract_features(features,traces, ...
%         use_errorbars,max_delay,num_bootstraps);  
%     
%     for f = 1:length(features)
%         feature_plots{case_num, f}(i) = fps(f);
%         feature_errors{case_num, f}(i) = fes(f);
%     end  
%     
% end
% case_num = case_num + 1;
% 
% %sequence dependent pausing
% pause_loc = .9;
% mean_pause = 4;
% sd_pause = 2;
% for i = 1:length(elong_times)
%     traces = gen_data_pauses(elong_times(i), pause_loc, mean_pause, ...
%         sd_pause,time_res, points_per_trace, num_traces, num_states, ...
%         trans_mat, rna_per_sec, fluo_per_rna, alpha_perc, init_dist, noise);
%     
%     [fps, fes] = extract_features(features,traces, ...
%         use_errorbars,max_delay,num_bootstraps);  
%     
%     for f = 1:length(features)
%         feature_plots{case_num, f}(i) = fps(f);
%         feature_errors{case_num, f}(i) = fes(f);
%     end
% end
% 
% case_num = case_num + 1;

for f = 1:length(features)
    fig = figure();
    
    for i = 1:num_cases
        if use_errorbars
            errorbar(construct_lengths, feature_plots{i,f}, feature_errors{i,f}, '-o');
        else
            plot(construct_lengths, feature_plots{i,f}, '-o')
        end
        hold on
        if num_cases <= 1
            p = polyfit(construct_lengths, feature_plots{i,f}, 1);
            plot(construct_lengths, polyval(p, construct_lengths));
            hold on
        end
    end
    xlabel('Size of gene (kb)');
    ylabel(strcat(feature_names{f}, ' (', num2str(time_res), ' second delays)'));
    title(['Comparing ' feature_names{f}]);
    legend(case_names);
end

