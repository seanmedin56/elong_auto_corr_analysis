function real_pars = matlab_hmc(auto_cor, start_estimate, priors, num_eigs, hs)
%MATLAB_HMC Summary of this function goes here
%   Detailed explanation goes here

    % defines infitissimal for variables
    elong_inf = .00001;
    rise_inf = .0001;
    a_inf = .001;
    b_inf = .000001;
    
    % defines priors
    elong_prior_mean = priors(1);
    elong_prior_std = priors(2);
    rise_prior_mean = priors(3);
    rise_prior_std = priors(4);
    a_prior_mean = priors(5);
    a_prior_std = priors(6);
    b_prior_mean = priors(7);
    b_prior_std = priors(8);
    log_noise_prior_mean = 0;
    log_noise_prior_std = 1;

    num_derivs = [1,2];
    cor_fun = @(elong,alph_perc,tau,aes,bes) full_func_cor(elong,alph_perc,...
        tau,aes,bes);
    delay_offset = 2;
    
    logpdf = @(Parameters) log_posterior(Parameters);
    
    % creates hmc sampler
    log_noise_start = 0;
    start_estimate = log(start_estimate);
    start_estimate(2) = sqrt(-start_estimate(2));
    smp = hmcSampler(logpdf,[start_estimate'; log_noise_start], 'NumSteps', 25);
    
    % estimates maximum a posteriori
    [MAPpars,fitInfo] = estimateMAP(smp,'VerbosityLevel',0);
    inter_MAPpars = MAPpars;
    inter_MAPpars(2) = -inter_MAPpars(2)^2;
    real_pars = exp(inter_MAPpars);
    real_pars = real_pars(1:end-1);
    
    % tunes sample parameters
    [smp,tuneinfo] = tuneSampler(smp,'Start',MAPpars);
     
    figure;
    plot(tuneinfo.StepSizeTuningInfo.StepSizeProfile);
    xlabel('Iteration');
    ylabel('Step size');

    accratio = tuneinfo.StepSizeTuningInfo.AcceptanceRatio
    
    % draw samples
    NumChains = 4;
    chains = cell(NumChains,1);
    Burnin = 500;
    NumSamples = 1000;
    for c = 1:NumChains
        if (c == 1)
            level = 1;
        else
            level = 0;
        end
        randomness = randn(size(MAPpars)) * 0.1;
        chains{c} = drawSamples(smp,'Start',MAPpars + randomness, ...
            'Burnin',Burnin,'NumSamples',NumSamples, ...
            'VerbosityLevel',level,'NumPrint',300);
    end 
    diags = diagnostics(smp,chains)
    
    % plots elongation and rise time statistics
    all_samples = vertcat(chains{:});
    all_samples(:,2) = - all_samples(:,2)  .^ 2;
    all_samples = exp(all_samples);
    figure();
    histogram(all_samples(:,1));
    figure();
    histogram(all_samples(:,2));
    figure();
    hist3(all_samples(:,1:2), 'Nbins', [20,20]);
    
    plot_fit_with_sims(auto_cor, cor_fun, real_pars, num_eigs, delay_offset, hs)
    
    function [log_pdf, grad_log_pdf] = log_posterior(Parameters)
        
        log_pdf = 0;
        grad_log_pdf = zeros(length(Parameters),1);
        for num_deriv = num_derivs
            
            to_fit = auto_cor / auto_cor(delay_offset + 1);
            for i = 1:num_deriv
                to_fit = to_fit(2:end) - to_fit(1:end-1);
            end
            
            elong = exp(Parameters(1));
            rise = exp(-Parameters(2)^2);
            aes = exp(Parameters(3:2 + num_eigs));
            bes = exp(Parameters(3 + num_eigs:end-1));
            log_noise = Parameters(end);
            sigma = sqrt(exp(log_noise));
            
            fit_auto_cor = zeros(1, length(to_fit));
            if num_deriv == 0
                del_offset = 0;
            else
                del_offset = 1;
            end
            for delay = 0:(length(to_fit) - 1)
                fit_auto_cor(delay + 1) = full_func_cor_deriv(elong,rise, ...
                delay + del_offset,aes,bes, num_deriv,cor_fun,delay_offset);
            end
            
            Z = (fit_auto_cor(1 + delay_offset:end) - ...
                to_fit(1 + delay_offset:end)) / sigma;
            
            % compute log liklihood and gradients
            loglik = sum(-log(sigma) - .5*log(2*pi) - .5*Z.^2);
            grad_elongs = zeros(1, length(Z));
            grad_rises = zeros(1, length(Z));
            grad_aes = zeros(num_eigs, length(Z));
            grad_bes = zeros(num_eigs, length(Z));
            for i = 1:length(Z)
                delay = i + delay_offset + del_offset - 1;
                grad_elongs(i) = (full_func_cor_deriv(elong+elong_inf,rise,...
                    delay, aes,  bes, num_deriv,cor_fun,delay_offset) - ...
                    fit_auto_cor(delay)) / elong_inf;
                grad_rises(i) = (full_func_cor_deriv(elong,rise+rise_inf,...
                    delay, aes,  bes, num_deriv,cor_fun,delay_offset) - ...
                    fit_auto_cor(delay)) / rise_inf;
                for j = 1:num_eigs
                    new_aes = aes;
                    new_aes(j) = aes(j) + a_inf;
                    new_bes = bes;
                    new_bes(j) = bes(j) + b_inf;
                    
                    grad_aes(j,i) = (full_func_cor_deriv(elong,rise,...
                    delay, new_aes,  bes, num_deriv,cor_fun,delay_offset) - ...
                    full_func_cor_deriv(elong,rise, delay,aes,bes, num_deriv, ...
                    cor_fun,delay_offset)) / a_inf;
                    
                    grad_bes(j,i) = (full_func_cor_deriv(elong,rise,...
                    delay, aes,  new_bes, num_deriv,cor_fun,delay_offset) - ...
                    full_func_cor_deriv(elong,rise, delay,aes,bes, num_deriv, ...
                    cor_fun,delay_offset)) / b_inf;
                end
            end
            grad_elong = -grad_elongs * Z' * elong / sigma;
            grad_rise = -grad_rises * Z' * (-2 * Parameters(2) * rise) / sigma;
            grad_a = -grad_aes .* aes * Z' / sigma;
            grad_b = -grad_bes .* bes * Z' / sigma;
            grad_noise = sum(-.5 + .5*(Z.^2));
            grad_log_pdf = grad_log_pdf + [grad_elong; grad_rise; grad_a; grad_b; grad_noise];
            log_pdf = loglik + log_pdf;        
        end
        
        % calculate priors
        [elong_prior, elong_grad_prior] = normalPrior(elong_prior_mean, ...
            elong, elong_prior_std);
        [rise_prior, rise_grad_prior] = normalPrior(rise_prior_mean, ...
            rise, rise_prior_std);
        [a_prior, a_grad_prior] = normalPrior(a_prior_mean, ...
            aes, a_prior_std);
        [b_prior, b_grad_prior] = normalPrior(b_prior_mean, ...
            bes, b_prior_std);
        [log_noise_prior, log_noise_grad_prior] = normalPrior(...
            log_noise_prior_mean, log_noise, log_noise_prior_std);

        log_prior = elong_prior + rise_prior + a_prior + b_prior + log_noise_prior;
        grad_prior = [elong_grad_prior; rise_grad_prior; ...
            a_grad_prior; b_grad_prior; log_noise_grad_prior];
        
        log_pdf = log_pdf + log_prior;
        grad_log_pdf = grad_log_pdf + grad_prior;
        
    end  

    % log prior probability function
     function [logpdf,gradlogpdf] = normalPrior(P,Mu,Sigma)
        Z          = (P - Mu)./Sigma;
        logpdf     = sum(-log(Sigma) - .5*log(2*pi) - .5*(Z.^2));
        gradlogpdf = Z./Sigma;
    end

end

