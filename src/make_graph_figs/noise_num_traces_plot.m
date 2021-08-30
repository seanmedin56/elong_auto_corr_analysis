% makes SI plot showing how noise and interpolation effect simple
% autocorrelation model and derivatives

% key variables
elong_time = 80;
time_res = 10;
rise_time = 40;
points_per_trace = 240; % 240
rna_per_sec = 0.2; 
noises = [0,700,1400];
num_traces = 400;
fluo_per_rna = 350; % not that it matters
prob_start_on = 1; 
k_on = 0.01;
k_off = 0.05;
cut = 10;

max_delay = 18;

% generates simple poisson promoter plots
traces = cell(1,num_traces);
traces2 = cell(1,num_traces);
traces3 = cell(1,num_traces);
for i = 1:num_traces
    traces{i} = gillespie_gen(elong_time, time_res, points_per_trace, ...
                                 2, [-k_on, k_off; k_on, -k_off], [.001, rna_per_sec], ...
                                 fluo_per_rna, rise_time,[0,1],noises(1));
    traces{i} = traces{i}(1+cut:end);   
    
    traces2{i} = traces{i} + normrnd(0,noises(2),1,points_per_trace-cut);
    traces2{i}(traces2{i} < 0) = 0;
    traces3{i} = traces{i} + normrnd(0,noises(3),1,points_per_trace-cut);
    traces3{i}(traces3{i} < 0) = 0;
end


% simulated autocorrelation
auto = fin_corr4(traces, traces, max_delay);
auto2 = fin_corr4(traces2, traces2, max_delay);
auto3 = fin_corr4(traces3, traces3, max_delay);

% simulates autocorrelation bootstraps
auto_stds = corr_bootstraps(traces, traces, max_delay, 100, 0, "c4");
auto_stds2 = corr_bootstraps(traces2, traces2, max_delay, 100, 0, "c4");
auto_stds3 = corr_bootstraps(traces3, traces3, max_delay, 100, 0, "c4");

figure('DefaultAxesFontSize',10)
errorbar(0:length(auto)-1, auto, auto_stds{1}, 'o-');
hold on
errorbar(0:length(auto)-1, auto2, auto_stds2{1}, 'o-');
hold on
errorbar(0:length(auto)-1, auto3, auto_stds3{1}, 'o-');

xline(8, '-', 'Elongation Time');
xlabel('time delay');
ylabel(['correlation']);
grid on
legend({'No Noise', 'Medium Noise', 'Heavy Noise'});
%saveas(gcf, ['./../../fig/noise2_figs/auto_noise.svg']);
%close;

cur_corr = auto;
cur_corr2 = auto2;
cur_corr3 = auto3;
for i = 1:3
    cur_corr = diff(cur_corr); % calculates latest deriv
    cur_corr2 = diff(cur_corr2);
    cur_corr3 = diff(cur_corr3);

    corr_stds = corr_bootstraps(traces, traces, max_delay, 100, i, "c4");
    corr_stds2 = corr_bootstraps(traces2, traces2, max_delay, 100, i, "c4");
    corr_stds3 = corr_bootstraps(traces3, traces3, max_delay, 100, i, "c4");

    figure('DefaultAxesFontSize',10)
    errorbar(2:length(cur_corr), cur_corr(2:end), corr_stds{i+1}(2:end), 'o-');
    hold on
    errorbar(2:length(cur_corr), cur_corr2(2:end), corr_stds2{i+1}(2:end), 'o-');
    hold on
    errorbar(2:length(cur_corr), cur_corr3(2:end), corr_stds3{i+1}(2:end), 'o-');

    xline(8, '-', 'Elongation Time');
    xlabel('time delay');
    ylabel(['\Delta^' int2str(i) ' correlation']);
    grid on
    legend({'No Noise', 'Medium Noise', 'Heavy Noise'}, 'Location', 'southeast');
    %saveas(gcf, ['./../../fig/noise2_figs/deriv' int2str(i) '.svg']);
    %close;
end