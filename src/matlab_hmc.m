function model = matlab_hmc(auto_cor, start_estimate)
%MATLAB_HMC Summary of this function goes here
%   Detailed explanation goes here
    
    num_vars = length(start_estimate);
    num_eigs = (num_vars - 2) / 2;
    num_deriv = 2;
    to_fit = auto_cor(2:end) / auto_cor(2);
    cor_fun = @(elong,alph,tau,aes,bes) full_func_cor(elong,alph,tau,aes,bes);

    logpdf = @(vars) [0]; %vars = [el,al,a1,a2,...b1,b2,...]
    for i = 1:length(to_fit)
        logpdf = @(vars) logpdf(vars) + log(abs(full_func_cor_deriv(vars(1),vars(2), ...
            i,vars(3:2+num_eigs),vars(3+num_eigs:end),num_deriv,cor_fun) - to_fit(i)));
    end  
    model = hmcSampler(logpdf, start_estimate);
    
end

