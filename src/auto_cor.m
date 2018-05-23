%examines the autocorrelation statistics of a trace
%saves a plot of R and M for the autocorrelation
%also saves a plot for the first and second derivative

addpath('utilities/');

if exist('project') ~= 1
    warning('Define "project" var (str): project identifier');
    return
end

% cut determines how many of the first few points of each trace to discard
if exist('cut') ~= 1
    cut = 1;
end
project = [project 'cut' num2str(cut)];

% create directories for outputing plots
dirs = {['../out/auto_corr/' project '/']};

for i = 1:numel(dirs)
    if (exist(dirs{i}, 'dir') ~= 7)
        mkdir(dirs{i});
    end
end

%checks if traces exist in the workspace
if exist('traces') ~= 1
    warning('Define set of traces (cell array of arrays)');
    return
end

%checks if max_delay exist in the workspace
if exist('max_delay') ~= 1
    warning('Define how many time delay points to be analyzed (int)');
    return
end

for i = 1:length(traces)
    traces{i} = traces{i}(cut:end);
end

% -------------plots the autocovariance and the raw moment----------------

corr_m = auto_corr_m_calc_norm(traces, max_delay);

% sample code for generating bootstraps for the autocorrelation and the
% second derivative of the autocorrelation
%bootstrap_m = corr_bootstraps(traces, max_delay, 1000, 'm');
%boot_m_2nd = corr_2nd_deriv_bootstraps(traces, max_delay, 1000, 'm');

h = figure;
plot(0:max_delay-1, corr_m);
title('central moment');
xlabel('time delay');
grid on
savefig([dirs{1} 'central_moment.fig']);
close(h);

% -----------plots the first and second derivatives ----------------------
% -----------of the autocovariance and the raw moment --------------------

m_1st_deriv = corr_m(2:max_delay) - corr_m(1:max_delay - 1);

m_2nd_deriv = m_1st_deriv(2:max_delay - 1) - m_1st_deriv(1:max_delay - 2);

h = figure;
plot(m_1st_deriv);
title('first derivative of central moment');
xlabel('time delay');
savefig([dirs{1} 'central_moment_1st_deriv.fig']);
close(h);

h = figure;
plot(m_2nd_deriv);
title('second derivative of central moment');
xlabel('time delay');
savefig([dirs{1} 'central_moment_2nd_deriv.fig']);
close(h);
