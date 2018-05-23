% no longer in use (will delete if not used again next time I clean up
% code)
function opt = decay_line_fit(auto_cor,upper_limits,lower_limits)
%Fits the elongation correlation term plus an exponential decay term for
%the dynamics correlation to the given auto correlation with a nonlinear
%least squares fitting
%  auto_cor: the auto correlation function that we are trying to fit
%  x0: intial values for [a,b,c,el,al]
%  el: estimated elongation time
%  al: estimated rise time
%  a,b,c: miscellaneous other parameters

    addpath('utilities/');
    
    %generate function which is a combination of an expoential decay term
    %and an elongation term
    f = @(vars) [0]; %vars = [a,b,c,el,al]
    for i = 1:length(auto_cor)
        f = @(vars) [f(vars) (dynam_term(vars(4),vars(5),i-1,vars(2)) + ...
            vars(1)*elong_term2(vars(4),vars(5),i-1) - vars(3)) / ...
            (dynam_term(vars(4),vars(5),0,vars(2)) + ...
            vars(1)*elong_term2(vars(4),vars(5),0) - vars(3)) - auto_cor(i)];
    end
    
    % run non linear least squares on function with multiple random
    % starting points and choose the one with the lowest error
    low_err = 10000;
    opt = [0,0,0,0,0];
    options = optimoptions('lsqnonlin', 'FunctionTolerance', 1e-8, ...
        'StepTolerance', 1e-8, 'OptimalityTolerance', 1e-8, 'display', 'off');
    for j = 1:10
        x0 = zeros(1,5);
        for i =1:5
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
    %plot result
    approx = zeros(1,length(auto_cor));
    approx_elong = zeros(1, length(auto_cor));
    approx_dynam = zeros(1, length(auto_cor));
    for i = 1:length(approx)
        approx_elong(i) = opt(1)*elong_term(opt(4),opt(5),i-1);
        approx_dynam(i) = dynam_term(opt(4),opt(5),i-1,opt(2));
        approx(i) = (dynam_term(opt(4),opt(5),i-1,opt(2)) + ...
            opt(1)*elong_term2(opt(4),opt(5),i-1) - opt(3)) / ...
            (dynam_term(opt(4),opt(5),0,opt(2)) + ...
            opt(1)*elong_term2(opt(4),opt(5),0) - opt(3));
    end
    figure();
    grid on
    plot(0:length(approx)-1,approx);
    hold on
    plot(0:length(approx_elong)-1,approx_elong);
    hold on
    plot(0:length(approx_dynam)-1,approx_dynam);
    figure();
    grid on
    plot(0:length(approx)-1,approx);
end

