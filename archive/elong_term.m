% note: this is an approximation of the analytical expression, for exact
% expression, see full_func_cor
function elong_cor = elong_term(elong,alph,tau)
%Calculates the elongation correlation term
%   Detailed explanation goes here
    elong_cor = 0;
    if tau < elong
        for i=1:(elong-tau)
            elong_cor = elong_cor + alpha(elong-i,alph)*alpha(elong-i-tau,alph);
        end
        if elong - tau < 1
            i = 0;
        end
        elong_cor = elong_cor + (elong - tau - i) * alpha(tau,alph)*alpha(0,alph);
    end
end

