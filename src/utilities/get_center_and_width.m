function [center,width] = get_center_and_width(deriv2)
% Returns the center and the width of the 2nd derivative of the
% autocorrelation function (deriv2)

deriv2 = deriv2(2:end);

peak = max(deriv2);

idxes = find(deriv2 > peak / 2);

diff_begin_len_end = (idxes(end) - idxes(1) + 1) - length(idxes);

if diff_begin_len_end ~= 0
    disp(idxes);
    fixer = 0;
    if idxes(end) - idxes(end-1) > 1
        
        fixer = fixer + idxes(end) - idxes(end-1) - 1;
        idxes = idxes(1:end-1);
    end
    if idxes(2) - idxes(1) > 1
        
        fixer = fixer + idxes(2) - idxes(1) - 1;
        idxes = idxes(2:end);
    end
    if fixer == diff_begin_len_end
        width = length(idxes);
        center = (idxes(end) + idxes(1)) / 2;
    else
        center = -1;
        width = 0;
        return
    end
else
    width = length(idxes);
    center = (idxes(end) + idxes(1)) / 2;
end

end

