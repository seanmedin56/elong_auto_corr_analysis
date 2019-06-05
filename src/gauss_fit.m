function opt = gauss_fit(auto_cor,initial_guesses, upper_limits,lower_limits,hs)
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
    num_eigs = (num_vars - 4) / 2;
    num_derivs = [1,2];
    delay_offset = 2;
    max_delay = 100;
    cor_fun = @(elong,elong_spread, alph, alph_spread, tau,aes,bes) ...
        full_func_cor_gauss(elong,elong_spread, alph, alph_spread, tau,aes,bes);
    
    f = @(vars) [0]; %vars = [el,al,a1,a2,...b1,b2,..., el_sp, al_sp]
    for num_deriv = num_derivs
        to_fit = auto_cor / auto_cor(delay_offset + 1);
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
        for i = (1 + delay_offset):min(length(to_fit), max_delay)
            f = @(vars) [f(vars) full_func_cor_deriv_gauss(vars(1), vars(end - 1), ...
                vars(2), vars(end),i - del_offset,vars(3:2+num_eigs), ...
                vars(3+num_eigs:end - 2),num_deriv,cor_fun,delay_offset) - to_fit(i)];
        end  
    end
    % run non linear least squares on function with multiple random
    % starting points and choose the one with the lowest error
    low_err = 10000;
    opt = zeros(1,num_vars);
    options = optimoptions('lsqnonlin', 'FunctionTolerance', 1e-7, ...
        'StepTolerance', 1e-7, 'OptimalityTolerance', 1e-7, 'display', 'off');
    for j = 1:8
        if j == 1
            x0 = initial_guesses;
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
    if exist('hs', 'var')
        %plot result
        approx = zeros(1,length(auto_cor));
        for i = 1:length(approx)
            approx(i) = cor_fun(opt(1),opt(end-1),opt(2),opt(end),i-1,...
                opt(3:2+num_eigs), ...
                opt(3+num_eigs:end)) / cor_fun(opt(1),opt(end-1),opt(2), ...
                opt(end),0, opt(3:2+num_eigs),opt(3+num_eigs:end));
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

