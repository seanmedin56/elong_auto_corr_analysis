function create_error_heatmap(auto_cor,upper_limits,lower_limits)
%CR Summary of this function goes here
%   Detailed explanation goes here
    addpath('utilities/');
    num_vars = length(upper_limits);
    num_eigs = (num_vars - 2) / 2;
    num_derivs = [1,2];
    delay_offset = 2;
    num_iters1 = 50;
    num_iters2 = 30;
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
    
    map = zeros(num_iters1, num_iters1);
    all_elongs = lower_limits(1):((upper_limits(1) - ...
        lower_limits(1)) / (num_iters1 - 1)):upper_limits(1);
    all_rises = lower_limits(2):((upper_limits(2) - ...
        lower_limits(2)) / (num_iters1 - 1)):upper_limits(2);
    eig_combos = {[]};
    for i = 1:num_eigs
        a_idx = 3 + (i - 1) * 2;
        b_idx = 4 + (i - 1) * 2;
        aes = lower_limits(a_idx):((upper_limits(a_idx) - ...
            lower_limits(a_idx)) / (num_iters2 - 1)):upper_limits(a_idx);
        bes = lower_limits(b_idx):((upper_limits(b_idx) - ...
            lower_limits(b_idx)) / (num_iters2 - 1)):upper_limits(b_idx);
        new_eig_combos = cell(1, num_iters2 ^ (2 * i));
        idx = 1;
        for a = aes
            for b = bes
                for j = 1:length(eig_combos)
                    new_eig_combos{idx} = [eig_combos{j} a b];
                    idx = idx + 1;
                end
            end
        end
        eig_combos = new_eig_combos;
    end
    
    for i = 1:length(all_elongs)
        for j = 1:length(all_rises)
            best = 10000;
            for k = 1:length(eig_combos)
                err = sum(f([all_elongs(i), all_rises(j), eig_combos{k}]).^2);
                best = min(err, best);
            end
            map(i,j) = best;
        end
    end
    
    figure();
    heatmap(all_elongs, all_rises, map);
    colormap jet
end

