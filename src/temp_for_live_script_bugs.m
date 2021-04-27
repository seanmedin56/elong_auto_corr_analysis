max_delay = 20;
set(gca,'FontSize',16)
% Finds longest cleaned up trace
longest = [];
for i = 1:length(tracesp)
    if length(tracesp{i}) > length(longest)
        longest = tracesp{i};
    end
end

figure();
plot(longest);
xlabel('time points');
ylabel('fluorescence (a.u)');
title('Sample Trace')

all = auto_corr_m_calc_norm(tracesp, max_delay);
std_derivs = corr_bootstraps(tracesp,tracesp, max_delay,100, 0,'m');
figure()
a = errorbar(0:(max_delay-1), all, std_derivs{1});


for i=1:length(tracesp)
    hold on
    if length(tracesp{i}) > 80
        single = auto_corr_m_calc_norm({tracesp{i}}, max_delay);
        a = plot(0:(max_delay-1), single);
        a.Color = [0,0,0,0.08];
    end
end



title('Autocorrelation of Traces')
xlabel('time delay');
ylabel('normalized autocorrelation')
legend({'All Traces', 'Individual Traces'})