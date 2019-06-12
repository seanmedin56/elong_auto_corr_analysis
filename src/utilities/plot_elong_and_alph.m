function plot_elong_and_alph(hs, elong, alph)
%PLOT_ELONG_AND_ALPH Summary of this function goes here
%   Detailed explanation goes here
    if ishandle(hs(1))
        set(0, 'CurrentFigure', hs(1));
        xline(elong, '-', 'Elongation Time');
        xline(alph, '-', 'Rise Time');
    end

    if length(hs) > 1 && ishandle(hs(2)) && strcmp(get(hs(2),'type'),'figure')
        set(0, 'CurrentFigure', hs(2));
        xline(elong, '-', 'Elongation Time');
        xline(alph, '-', 'Rise Time');
    end
    if length(hs) > 2 && ishandle(hs(3))
        set(0, 'CurrentFigure', hs(3));
        xline(elong, '-', 'Elongation Time');
        xline(alph, '-', 'Rise Time');
    end
end

