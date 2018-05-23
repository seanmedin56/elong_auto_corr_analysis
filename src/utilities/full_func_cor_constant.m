% EXPERIMENTAL WAY OF CALCULATING ANALYTICAL EXPRESSOIN, DOES NOT WORK,
% MIGHT TRY AGAIN LATER
function cor_tot = full_func_cor_constant(elong,alph,tau,aes,bes)
% calculates autocorrelation at time delay tau with a constant elongation
% time elong
%   Detailed explanation goes here
 
    cor_tot = 0;
    %calculates poisson term

    if elong > tau
        fun = @(t1) min(1,(t1 + tau) / alph) * min(1, tau / alph);
        cor_tot = integral(fun, 0, elong - tau);
    end
    
    %calculates dynamics terms
    for i = 1:length(aes)
        a = aes(i);
        b = bes(i);
        fun = @(t1,t2) min(1, t1 ./ alph) .* min(1,t2 / alph) .* ... 
            exp(-b .* abs(t2 + tau - t1));
        cor_tot = cor_tot + a * integral2(fun,0,elong,0,elong);
    end
    

end



