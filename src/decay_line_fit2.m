function opt = decay_line_fit2(auto_cor,initial_conditions, upper_limits,lower_limits,hs)
%Fits the elongation correlation term plus an exponential decay term for
%the dynamics correlation to the given auto correlation with a nonlinear
%least squares fitting
%  auto_cor: the auto correlation function that we are trying to fit
%  x0: intial values for [b,el,al]
%  el: estimated elongation time
%  al: estimated rise time
%  b: parameter representing the dynamics

    addpath('utilities/');
    num_vars = length(upper_limits);
    num_eigs = (num_vars - 2) / 2;
    num_derivs = [1,2];
    delay_offset = 2;
    max_delay = 100;
    cor_fun = @(elong,alph_perc,tau,aes,bes) full_func_cor(elong,alph_perc,tau,aes,bes);
    f = @(vars) [0]; %vars = [el,al,a1,a2,...b1,b2,...]

    for num_deriv = num_derivs
        to_fit = auto_cor / auto_cor(delay_offset + 1);

        %num_deriv determines what's being fit
        for i = 1:num_deriv
            to_fit = to_fit(2:end) - to_fit(1:end-1);
        end

        %generate function which is the expectation value of the
        %autocorrelation
        if num_deriv == 0
            del_offset = 1;
        else
            del_offset = 0;
        end
        for i = (1 + delay_offset):min(length(to_fit),max_delay)
            f = @(vars) [f(vars) (full_func_cor_deriv(vars(1),vars(2), ...
                i - del_offset,vars(3:2+num_eigs),vars(3+num_eigs:end),...
                num_deriv,cor_fun,delay_offset) - to_fit(i))];
        end  
    end
    % run non linear least squares on function with multiple random
    % starting points and choose the one with the lowest error
    low_err = 10000;
    opt = zeros(1,num_vars);
    options = optimoptions('lsqnonlin', 'FunctionTolerance', 1e-7, ...
        'StepTolerance', 1e-7, 'OptimalityTolerance', 1e-7, 'display', 'off');
    for j = 1:5

        if j == 1
            x0 = initial_conditions;
        else
            x0 = zeros(1,num_vars);
            for i =1:num_vars
                x0(i) = rand() * (upper_limits(i) - lower_limits(i)) + lower_limits(i);
            end
        end
        [x,err] = lsqnonlin(f,x0,lower_limits,upper_limits,options);
        if err < low_err
            low_err = err;
            opt = x;
            disp(low_err);
            disp(opt);
        end
    end
    display(low_err);
    % plots result
    if exist('hs', 'var')
        plot_fit_with_sims(auto_cor, cor_fun, opt, num_eigs, delay_offset,hs)
    end
end

