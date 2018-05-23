function cor_deriv = full_func_cor_deriv(elong,alph,tau,a,b,deriv,full_func_cor)
% takes derivative of expectation value of autocorrelation at give point
%   Detailed explanation goes here
    
    if deriv == 0
        cor_deriv = full_func_cor(elong,alph,tau,a,b);
    elseif deriv == 1
        cor_deriv = full_func_cor(elong,alph,tau+1,a,b) -  ...
            full_func_cor(elong,alph,tau,a,b);
    else
        first = full_func_cor(elong,alph,tau,a,b);
        if first == 0 && tau == 1
            debug = true;
        end
        second = full_func_cor(elong,alph,tau + 1,a,b);
        third = full_func_cor(elong,alph,tau + 2,a, b);
        cor_deriv = third + first - 2 * second;
    end
    cor_deriv = cor_deriv / full_func_cor(elong,alph,1,a,b);
    
    
end

