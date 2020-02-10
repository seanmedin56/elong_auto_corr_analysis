function plot_expected_auto_cor(hs, elong, alph_perc, aes, bes, max_delay,...
    is_gauss, elong_spread, rise_spread)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

[M,I] = sort(abs(bes));
aes = aes(I(2:end));
bes = bes(I(2:end));

% calculates and normalizes expected autocorrelation
new_auto_cor = zeros(1, max_delay + 1);
if is_gauss
    for i = 0:max_delay
        new_auto_cor(i + 1) = full_func_cor_gauss(elong,elong_spread, ...
            alph_perc,rise_spread,i,aes,bes);
    end
else
    for i = 0:max_delay
        new_auto_cor(i + 1) = full_func_cor(elong,alph_perc,i,aes,bes);
    end
end
new_auto_cor = new_auto_cor / new_auto_cor(1);

% plots appropriate expressions
if ishandle(hs(1))
    set(0, 'CurrentFigure', hs(1));
    hold on
    plot(0:length(new_auto_cor)-1,new_auto_cor);
    legend('Simulated', 'Expected Value');
end
deriv1 = new_auto_cor(2:end) - new_auto_cor(1:end-1);
if length(hs) > 1 && ishandle(hs(2)) && strcmp(get(hs(2),'type'),'figure')
    set(0, 'CurrentFigure', hs(2));
    hold on
    plot(deriv1);
    legend('Simulated', 'Expected Value');
end
deriv2 = deriv1(2:end) - deriv1(1:end-1);
if length(hs) > 2 && ishandle(hs(3)) && strcmp(get(hs(3),'type'),'figure')   
    set(0, 'CurrentFigure', hs(3));
    hold on
    plot(deriv2);
    legend('Simulated', 'Expected Value');
end
deriv3 = diff(deriv2);
if length(hs) > 3 && ishandle(hs(4)) && strcmp(get(hs(4),'type'),'figure')
    set(0, 'CurrentFigure', hs(4));
    hold on
    plot(deriv3);
    legend('Simulated', 'Expected Value');
end

