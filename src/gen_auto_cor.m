function hs = gen_auto_cor(traces, auto, first, second, ...
                                    bootstraps, max_delay,cut)
% Generates an autocorrelation and/or the derivatives of an
% autocorrelation for the traces
% traces: The traces we are taking the autocorrelation of
% auto: Boolean saying whether or not to include the autocorrelation
% first: Boolean determining if to include 1st derivative
% second: Boolean determining if to include 2nd derivative
% bootstraps: Boolean determining if to include bootstraps
% max_delay: How many points to take of the autocorrelation
% returns handles for plots
addpath('utilities/');
hs = [];
%cuts the traces by amount cut
for i=1:length(traces)
    traces{i} = traces{i}(1 + cut:end);
end

corr = auto_corr_m_calc_norm(traces, max_delay);

corr_1st = corr(2:max_delay) - corr(1:max_delay-1);

corr_2nd = corr_1st(2:max_delay-1) - corr_1st(1:max_delay-2);

if bootstraps
    [stds, stds1, stds2] = corr_bootstraps(traces,traces, max_delay,100,'m');
end

if auto
    h = figure;
    if bootstraps
        errorbar(0:max_delay-1,corr,stds);
    else
        plot(0:max_delay-1, corr);
    end
    title('Central Moment', 'FontSize', 14);
    xlabel('time delay', 'FontSize', 14);
    grid on
    hs(1) = h;
end

if first
    h = figure;
    if bootstraps
        errorbar(corr_1st,stds1);
    else
        plot(corr_1st);
    end
    title('Central Moment 1st Derivative', 'FontSize', 14);
    xlabel('time delay', 'FontSize', 14);
    grid on
    hs(2) = h;
end

if second
    h = figure;
    if bootstraps
        errorbar(corr_2nd,stds2);
    else
        plot(corr_2nd);
    end
    title('Central Moment 2nd Derivative', 'FontSize', 14);
    xlabel('time delay', 'FontSize', 14);
    grid on
    hs(3) = h;
end

