% makes short/long rise time regimes for SI

% key variables
elong_time = 80;
time_res = 10;
rise_time = 40;
points_per_trace = 40 * 60 / time_res; % 240
rna_per_sec = 0.2; % ON state
noise = 0;
num_traces = 400;%300
fluo_per_rna = 350; % not that it matters
prob_start_on = 1; 
cut = 10;

k_on = 0.01; 
k_off = 0.05; 

time_std1 = 8;
time_std2 = 16;

max_delay = 20;

% generates simple poisson promoter plots
traces = cell(1,num_traces);
for i = 1:num_traces
    traces{i} = gillespie_gen(elong_time, time_res, points_per_trace, ...
                                 2, [-k_on, k_off; k_on, -k_off], [.001, rna_per_sec], ...
                                 fluo_per_rna, rise_time,[0,1],noise);
    traces{i} = traces{i}(1+cut:end);
end

% simulated autocorrelation
auto = fin_corr4(traces, traces, max_delay);
auto_stds = corr_bootstraps(traces, traces, max_delay, 100, 0, "c4");

% expected autocorrelation
auto_eq = zeros(1, max_delay);
[aes, bes] = decompose_matrix([-k_on, k_off; k_on, -k_off]*10, [.001, rna_per_sec]*10);
[M,I] = sort(abs(bes));
aes = aes(I(2:end));
bes = abs(bes(I(2:end)));
norm_auto = full_func_cor(elong_time / time_res, rise_time / elong_time, 0, aes, bes);
for i = 0:(max_delay-1)
    auto_eq(i + 1) = full_func_cor(elong_time / time_res, ...
        rise_time / elong_time,i,aes,bes) / norm_auto;
end

figure('DefaultAxesFontSize',10)
errorbar(0:length(auto)-1, auto, auto_stds{1}, 'o-');
hold on
plot(0:length(auto)-1, auto_eq);
hold on
%plot(8, auto(9), 'x', 'MarkerEdgeColor', 'Red', ...
%    'MarkerSize', 12, 'LineWidth', 3);
xline(4, '-', 'Rise Time');
xline(8, '-', 'Elongation Time');
xlabel('time delay (time steps)');
ylabel(['correlation']);
grid on
legend({'Simulated', 'Expected Value'});
saveas(gcf, ['./../../fig/gauss_plot_figs/gauss0_auto.svg']);
%close;

cur_corr = auto;
cur_corr_eq = auto_eq;
for i = 1:3
    cur_corr = diff(cur_corr); % calculates latest deriv
    corr_stds = corr_bootstraps(traces, traces, max_delay, 100, i, "c4");
    cur_corr_eq = diff(cur_corr_eq);   
    
    figure('DefaultAxesFontSize',10)
    errorbar(1:length(cur_corr), cur_corr, corr_stds{i+1}, 'o-');
    hold on
    plot(1:length(cur_corr), cur_corr_eq);
    hold on
    %plot(8, cur_corr(8), 'x', 'MarkerEdgeColor', 'Red', ...
    %'MarkerSize', 12, 'LineWidth', 3);
    xline(4, '-', 'Rise Time');
    xline(8, '-', 'Elongation Time');
    xlabel('time delay (time steps)');
    ylabel(['\Delta^' int2str(i) ' correlation']);
    grid on
    legend({'Simulated', 'Expected Value'});
    saveas(gcf, ['./../../fig/gauss_plot_figs/gauss0_deriv' int2str(i) '.svg']);
    %close;
end

% generates shorter std plots
traces = cell(1,num_traces);
for i = 1:num_traces
    traces{i} = gillespie_gen_gauss(elong_time, time_res, points_per_trace, ...
                                 2, [-k_on, k_off; k_on, -k_off], [.001, rna_per_sec], ...
                                 fluo_per_rna, rise_time,[0,1],noise,time_std1);
    traces{i} = traces{i}(1+cut:end);
end

% simulated autocorrelation
auto = fin_corr4(traces, traces, max_delay);
auto_stds = corr_bootstraps(traces, traces, max_delay, 100, 0, "c4");

% expected autocorrelation
auto_eq = zeros(1, max_delay);
[aes, bes] = decompose_matrix([-k_on, k_off; k_on, -k_off]*10, [.001, rna_per_sec]*10);
[M,I] = sort(abs(bes));
aes = aes(I(2:end));
bes = abs(bes(I(2:end)));
norm_auto = full_func_cor(elong_time / time_res, rise_time / elong_time, 0, aes, bes);
for i = 0:(max_delay-1)
    auto_eq(i + 1) = full_func_cor(elong_time / time_res, ...
        rise_time / elong_time,i,aes,bes) / norm_auto;
end

figure('DefaultAxesFontSize',10)
errorbar(0:length(auto)-1, auto, auto_stds{1}, 'o-');
hold on
plot(0:length(auto)-1, auto_eq);
hold on
%plot(8, auto(9), 'x', 'MarkerEdgeColor', 'Red', ...
%    'MarkerSize', 12, 'LineWidth', 3);
xline(4, '-', 'Rise Time');
xline(8, '-', 'Elongation Time');
xlabel('time delay (time steps)');
ylabel(['correlation']);
grid on
legend({'Simulated', 'Expected Value'});
saveas(gcf, ['./../../fig/gauss_plot_figs/gauss1_auto.svg']);
%close;

cur_corr = auto;
cur_corr_eq = auto_eq;
for i = 1:3
    cur_corr = diff(cur_corr); % calculates latest deriv
    corr_stds = corr_bootstraps(traces, traces, max_delay, 100, i, "c4");
    cur_corr_eq = diff(cur_corr_eq);   
    
    figure('DefaultAxesFontSize',10)
    errorbar(1:length(cur_corr), cur_corr, corr_stds{i+1}, 'o-');
    hold on
    plot(1:length(cur_corr), cur_corr_eq);
    hold on
    %plot(8, cur_corr(8), 'x', 'MarkerEdgeColor', 'Red', ...
    %'MarkerSize', 12, 'LineWidth', 3);
    xline(4, '-', 'Rise Time');
    xline(8, '-', 'Elongation Time');
    xlabel('time delay (time steps)');
    ylabel(['\Delta^' int2str(i) ' correlation']);
    grid on
    legend({'Simulated', 'Expected Value'});
    saveas(gcf, ['./../../fig/gauss_plot_figs/gauss1_deriv' int2str(i) '.svg']);
    %close;
end

% generates larger std plots
traces = cell(1,num_traces);
for i = 1:num_traces
    traces{i} = gillespie_gen_gauss(elong_time, time_res, points_per_trace, ...
                                 2, [-k_on, k_off; k_on, -k_off], [.001, rna_per_sec], ...
                                 fluo_per_rna, rise_time,[0,1],noise,time_std2);
    traces{i} = traces{i}(1+cut:end);
end

% simulated autocorrelation
auto = fin_corr4(traces, traces, max_delay);
auto_stds = corr_bootstraps(traces, traces, max_delay, 100, 0, "c4");

% expected autocorrelation
auto_eq = zeros(1, max_delay);
[aes, bes] = decompose_matrix([-k_on, k_off; k_on, -k_off]*10, [.001, rna_per_sec]*10);
[M,I] = sort(abs(bes));
aes = aes(I(2:end));
bes = abs(bes(I(2:end)));
norm_auto = full_func_cor(elong_time / time_res, rise_time / elong_time, 0, aes, bes);
for i = 0:(max_delay-1)
    auto_eq(i + 1) = full_func_cor(elong_time / time_res, ...
        rise_time / elong_time,i,aes,bes) / norm_auto;
end

figure('DefaultAxesFontSize',10)
errorbar(0:length(auto)-1, auto, auto_stds{1}, 'o-');
hold on
plot(0:length(auto)-1, auto_eq);
hold on
%plot([2,6,8], auto([3,7,9]), 'x', 'MarkerEdgeColor', 'Red', ...
%    'MarkerSize', 12, 'LineWidth', 3);
xline(4, '-', 'Rise Time');
xline(8, '-', 'Elongation Time');
xlabel('time delay (time steps)');
ylabel(['correlation']);
grid on
legend({'Simulated', 'Expected Value'});
saveas(gcf, './../../fig/gauss_plot_figs/gauss2_auto.svg');
%close;

cur_corr = auto;
cur_corr_eq = auto_eq;
for i = 1:3
    cur_corr = diff(cur_corr); % calculates latest deriv
    corr_stds = corr_bootstraps(traces, traces, max_delay, 100, i, "c4");
    cur_corr_eq = diff(cur_corr_eq);   
    
    figure('DefaultAxesFontSize',10)
    errorbar(1:length(cur_corr), cur_corr, corr_stds{i+1}, 'o-');
    hold on
    plot(1:length(cur_corr), cur_corr_eq);
    hold on
    %plot([2,6,8], cur_corr([2,6,8]), 'x', 'MarkerEdgeColor', 'Red', ...
    %    'MarkerSize', 12, 'LineWidth', 3);
    xline(4, '-', 'Rise Time');
    xline(8, '-', 'Elongation Time');
    xlabel('time delay (time steps)');
    ylabel(['\Delta^' int2str(i) ' correlation']);
    grid on
    legend({'Simulated', 'Expected Value'});
    saveas(gcf, ['./../../fig/gauss_plot_figs/gauss2_deriv' int2str(i) '.svg']);
    %close;
end

