function cor_deriv = full_func_cor_deriv(elong,alph_perc,tau,a,b,deriv,cor_fun,normalizer)
% takes derivative of expectation value of autocorrelation at give point
%   Detailed explanation goes here
    
    if deriv == 0
        cor_deriv = cor_fun(elong,alph_perc,tau,a,b);
    elseif deriv == 1
        cor_deriv = cor_fun(elong,alph_perc,tau,a,b) -  ...
            cor_fun(elong,alph_perc,tau - 1,a,b);
    else
        first = cor_fun(elong,alph_perc,tau - 1,a,b);
        second = cor_fun(elong,alph_perc,tau,a,b);
        third = cor_fun(elong,alph_perc,tau + 1,a, b);
        cor_deriv = third + first - 2 * second;
    end
    cor_deriv = cor_deriv / cor_fun(elong,alph_perc,normalizer,a,b);
    
    
end

