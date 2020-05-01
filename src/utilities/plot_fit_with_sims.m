function plot_fit_with_sims(auto_cor, cor_fun, opt, num_eigs, delay_offset, hs)
%PLOT_FIT_WITH_SIMS Summary of this function goes here
%   Detailed explanation goes here

    approx = zeros(1,length(auto_cor));
    for i = 1:length(approx)
        approx(i) = cor_fun(opt(1),opt(2),i-1,opt(3:2+num_eigs), ...
            opt(3+num_eigs:end)) / cor_fun(opt(1),opt(2),0, ...
            opt(3:2+num_eigs),opt(3+num_eigs:end));
    end
    approx(1:delay_offset) = auto_cor(1:delay_offset);
    approx(delay_offset+1:end) = approx(delay_offset+1:end) * ...
        auto_cor(delay_offset+1) / approx(delay_offset+1);
    if ishandle(hs(1))
        set(0, 'CurrentFigure', hs(1));
        hold on
        plot(0:length(approx)-1,approx);
        legend('Original', 'Estimate');
    end
    deriv1 = approx(2:end) - approx(1:end-1);
    if length(hs) > 1 && ishandle(hs(2)) && strcmp(get(hs(2),'type'),'figure')
        set(0, 'CurrentFigure', hs(2));
        hold on
        plot(deriv1);
        legend('Original', 'Estimate');
    end
    deriv2 = deriv1(2:end) - deriv1(1:end-1);
    if length(hs) > 2 && ishandle(hs(3))
        set(0, 'CurrentFigure', hs(3));
        hold on
        plot(deriv2);
        legend('Original', 'Estimate');
    end
    deriv3 = deriv2(2:end) - deriv2(1:end-1);
    if length(hs) > 3 && ishandle(hs(4))
        set(0, 'CurrentFigure', hs(4));
        hold on
        plot(deriv3);
        legend('Original', 'Estimate');
    end
end

