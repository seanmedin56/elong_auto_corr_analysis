% note: this is not used anywhere at the moment
function elong_cor = elong_term2(elong,alph,tau)
% Calculates the elongation correlation term exactly
%   Detailed explanation goes here
    elong_cor = 0;
    if elong - tau - alph > 0
        elong_cor = elong_cor + (elong - tau - alph);
    end
    if tau ~= 0 && tau < alph
        elong_cor = elong_cor + tau * (1 - (tau / 2 / alph));
        
        elong_cor = elong_cor + alph / 3 + tau^3 / (6 * alph^2) - tau / 2;
    end
    if tau > alph
        elong_cor = elong_cor + alph / 2;
    end
    
end

