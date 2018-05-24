function peak_idx = extract_2nd_deriv_peak(corr)
% extracts peak of the 2nd derivative of the autocorrelaiton corr
%   Ignores 1st point in 2nd derivative when figuring out peak

    deriv1 = corr(2:end) - corr(1:end-1);
    deriv2 = deriv1(2:end) - deriv1(1:end-1);
    [~, peak_idx] = max(deriv2(2:end));
    peak_idx = peak_idx + 1;

end

