function cor_deriv = full_func_cor_deriv_gauss(elong, elong_spread, ...
    alph,alph_spread, tau,a,b, deriv,cor_fun,normalizer)
% takes derivative of expectation value of autocorrelation at give point
%   Detailed explanation goes here
    
    if deriv == 0
        cor_deriv = cor_fun(elong, elong_spread, alph, alph_spread,tau,a,b);
    elseif deriv == 1
        cor_deriv = cor_fun(elong, elong_spread, alph, alph_spread,tau,a,b) -  ...
            cor_fun(elong, elong_spread, alph, alph_spread,tau-1,a,b);
    else
        first = cor_fun(elong, elong_spread, alph, alph_spread,tau - 1,a,b);
        second = cor_fun(elong, elong_spread, alph, alph_spread,tau,a,b);
        third = cor_fun(elong, elong_spread, alph, alph_spread,tau + 1,a, b);
        cor_deriv = third + first - 2 * second;
    end
    cor_deriv = cor_deriv / cor_fun(elong, elong_spread, alph, alph_spread,normalizer,a,b);
    
    
end

