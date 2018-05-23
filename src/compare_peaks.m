% This script generates auto correlations for different elongation dynamics
% with different gene length and creates a graph comparing the widths and
% center of the peak of their second derivatives

addpath('utilities/');

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

construct_lengths = [4.7,3.4,2.65,1.9]; %in kb

base_time = 160; % 160 seconds to transribe 4.7kb

alpha_perc = .25; % percent of 4.7kb that is ms2 loops

max_delay = 45;

% generates data and extracts statistics for each case

center_fig = figure();
width_fig = figure();

%constant elongation time
elong_times = (base_time / construct_lengths(1)) * construct_lengths;
rise_time = alpha_perc * base_time;
centers = zeros(1,length(elong_times));
widths = zeros(1,length(elong_times));
for i = 1:length(elong_times)
    elong = elong_times(i);
    traces = gen_data(elong,time_res,points_per_trace,num_traces,num_states, ...
        trans_mat, rna_per_sec, fluo_per_rna, rise_time, init_dist,noise);
    
    corr = auto_corr_m_calc_norm(traces, max_delay);

    corr_1st = corr(2:max_delay) - corr(1:max_delay-1);

    corr_2nd = corr_1st(2:max_delay-1) - corr_1st(1:max_delay-2);
    
    [center,width] = get_center_and_width(corr_2nd);
    centers(i) = center;
    widths(i) = width;  
end
figure(center_fig);
plot(construct_lengths,centers, '-ob');
figure(width_fig);
plot(construct_lengths,widths, '-ob');

% diff rates elongation time gaussian

base_width = 20;
elong_widths = (base_width / construct_lengths(1)) * construct_lengths;
alpha_perces = (alpha_perc / construct_lengths(1)) * construct_lengths;

for i = 1:length(elong_times)
    elong_var = {'Gaussian', elong_times(i), elong_widths(i)};
    traces = gen_data_rate_spread(elong_var, time_res, points_per_trace, ...
        num_traces, num_states, trans_mat, rna_per_sec, fluo_per_rna, ...
        alpha_perces(i), init_dist, noise);
    
    corr = auto_corr_m_calc_norm(traces, max_delay);
    corr_1st = corr(2:max_delay) - corr(1:max_delay-1);
    corr_2nd = corr_1st(2:max_delay-1) - corr_1st(1:max_delay-2);
    
    [center,width] = get_center_and_width(corr_2nd);
    centers(i) = center;
    widths(i) = width; 
    
end

figure(center_fig);
hold on
plot(construct_lengths,centers, '-or');
figure(width_fig);
hold on
plot(construct_lengths,widths, '-or');

% different falloff points (no bottleneck)

for i = 1:length(elong_times)
    elong_var = {'Gaussian', elong_times(i), elong_widths(i)};
    traces = gen_data_var_len(elong_var, time_res, points_per_trace, ...
        num_traces, num_states, trans_mat, rna_per_sec, fluo_per_rna, ...
        rise_time, init_dist, noise);
    
    corr = auto_corr_m_calc_norm(traces, max_delay);
    corr_1st = corr(2:max_delay) - corr(1:max_delay-1);
    corr_2nd = corr_1st(2:max_delay-1) - corr_1st(1:max_delay-2);
    
    [center,width] = get_center_and_width(corr_2nd);
    centers(i) = center;
    widths(i) = width; 
    
end
figure(center_fig);
hold on
plot(construct_lengths,centers, '-og');
figure(width_fig);
hold on
plot(construct_lengths,widths, '-og');

%sequence dependent pausing
pause_loc = .9;
mean_pause = 4;
sd_pause = 2;
for i = 1:length(elong_times)
    traces = gen_data_pauses(elong_times(i), pause_loc, mean_pause, ...
        sd_pause,time_res, points_per_trace, num_traces, num_states, ...
        trans_mat, rna_per_sec, fluo_per_rna, alpha_perc, init_dist, noise);
    
    corr = auto_corr_m_calc_norm(traces, max_delay);
    corr_1st = corr(2:max_delay) - corr(1:max_delay-1);
    corr_2nd = corr_1st(2:max_delay-1) - corr_1st(1:max_delay-2);
    
    [center,width] = get_center_and_width(corr_2nd);
    centers(i) = center;
    widths(i) = width; 
    
end

figure(center_fig);
hold on
plot(construct_lengths,centers, '-om');
figure(width_fig);
hold on
plot(construct_lengths,widths, '-om');

legend('constant rate' , 'rate spread', 'var len', 'pausing');
title('Comparing Peak Widths');
xlabel('Size of gene (kb)');
ylabel('Width of peak (5 second delays)');

figure(center_fig);
legend('constant rate' , 'rate spread', 'var len', 'pausing');
title('Comparing Peak Centers');
xlabel('Size of gene (kb)');
ylabel('Center of peak (5 second delays)');


