function plot_elong_and_alph(hs, elong, alph)
%PLOT_ELONG_AND_ALPH Summary of this function goes here
%   Detailed explanation goes here
    if ishandle(hs(1))
        set(0, 'CurrentFigure', hs(1));
        xline(elong, '-', 'Elongation Time', 'HandleVisibility', 'Off');
        xline(alph, '-', 'Rise Time', 'HandleVisibility', 'Off');
    end

    if length(hs) > 1 && ishandle(hs(2)) && strcmp(get(hs(2),'type'),'figure')
        set(0, 'CurrentFigure', hs(2));
        xline(elong, '-', 'Elongation Time', 'HandleVisibility', 'Off');
        xline(alph, '-', 'Rise Time', 'HandleVisibility', 'Off');
    end
    if length(hs) > 2 && ishandle(hs(3)) && strcmp(get(hs(3),'type'),'figure')
        set(0, 'CurrentFigure', hs(3));
        xline(elong, '-', 'Elongation Time', 'HandleVisibility', 'Off');
        xline(alph, '-', 'Rise Time', 'HandleVisibility', 'Off');
    end
    if length(hs) > 3 && ishandle(hs(4)) && strcmp(get(hs(4),'type'),'figure')
        set(0, 'CurrentFigure', hs(4));
        xline(elong, '-', 'Elongation Time', 'HandleVisibility', 'Off');
        xline(alph, '-', 'Rise Time', 'HandleVisibility', 'Off');
    end
end

