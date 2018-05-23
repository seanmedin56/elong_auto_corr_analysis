function opt = decay_line_fit2(auto_cor,upper_limits,lower_limits,hs)
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
    num_deriv = 2;
    to_fit = auto_cor(2:end) / auto_cor(2);
    cor_fun = @(elong,alph,tau,aes,bes) full_func_cor(elong,alph,tau,aes,bes);
    
    %num_deriv determines what's being fit
    for i = 1:num_deriv
        to_fit = to_fit(2:end) - to_fit(1:end-1);
    end
    
    %generate function which is the expectation value of the
    %autocorrelation
    f = @(vars) [0]; %vars = [el,al,a1,a2,...b1,b2,...]
    for i = 1:length(to_fit)
        f = @(vars) [f(vars) full_func_cor_deriv(vars(1),vars(2), ...
            i,vars(3:2+num_eigs),vars(3+num_eigs:end),num_deriv,cor_fun) - to_fit(i)];
    end  
    
    % run non linear least squares on function with multiple random
    % starting points and choose the one with the lowest error
    low_err = 10000;
    opt = zeros(1,num_vars);
    options = optimoptions('lsqnonlin', 'FunctionTolerance', 1e-7, ...
        'StepTolerance', 1e-7, 'OptimalityTolerance', 1e-7, 'display', 'off');
    for j = 1:1
        x0 = zeros(1,num_vars);
        for i =1:num_vars
            x0(i) = rand() * (upper_limits(i) - lower_limits(i)) + lower_limits(i);
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
    if exist('hs', 'var')
        %plot result
        approx = zeros(1,length(auto_cor));
        for i = 1:length(approx)
            approx(i) = cor_fun(opt(1),opt(2),i-1,opt(3:2+num_eigs), ...
                opt(3+num_eigs:end)) / cor_fun(opt(1),opt(2),0, ...
                opt(3:2+num_eigs),opt(3+num_eigs:end));
        end
        if ishandle(hs(1))
            set(0, 'CurrentFigure', hs(1));
            hold on
            plot(0:length(approx)-1,approx);
            legend('Original', 'Estimate');
        end
        deriv1 = approx(2:end) - approx(1:end-1);
        if length(hs) > 1 && ishandle(hs(2)) && strcmp(get(hs(2),'type'),'figure')
            set(0, 'CurrentFigure', hs(2));
            hold on
            plot(deriv1);
            legend('Original', 'Estimate');
        end
        if length(hs) > 2 && ishandle(hs(3))
            deriv2 = deriv1(2:end) - deriv1(1:end-1);
            set(0, 'CurrentFigure', hs(3));
            hold on
            plot(deriv2);
            legend('Original', 'Estimate');
        end
    end
end

