function cor_tot = full_func_cor_gauss(elong, elong_spread, alph, ... 
    alph_spread, tau, aes, bes)
% calculates expected autocorrelation at time delay tau
%   elong: mean elongation time
%   elong_spread: standard deviation elongation time
%   alph: mean rise time
%   alph_spread: standard deviation rise time
%   tau: time delay
%   aes: coefficiencts for exponential terms
%   bes: decay values for exponential terms

    %how many points to take on either side of the mean for the gaussian
    num_points = 20; 
    left_dist_elong = max(2 * elong_spread, elong * 9 / 10);
    left_dist_alph = max(2 * alph_spread, alph * 9 / 10);
    elongs_to_take = elong - left_dist_elong: ... 
        left_dist_elong / num_points:elong + left_dist_elong;
    alphs_to_take = alph - left_dist_alph: ...
        left_dist_alph / num_points: alph + left_dist_alph;
    cor_tot = 0;
    prob_tot = 0;
    for i = 1:num_points * 2 + 1
        el = elongs_to_take(i);
        al = alphs_to_take(i);
        elong_prob = normpdf(el, elong, elong_spread);
        prob_tot = prob_tot + elong_prob;
        cor_tot = cor_tot + elong_prob * full_func_cor(el, al, tau, aes, bes);
    end
    cor_tot = cor_tot / prob_tot;     
end

