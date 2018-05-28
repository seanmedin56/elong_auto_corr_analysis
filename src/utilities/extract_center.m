function center = extract_center(corr)
% extracts the center of the peak of the 2nd derivative of the autocorrelaiton 
% corr
%   Ignores 1st point in 2nd derivative when figuring out peak

    deriv1 = corr(2:end) - corr(1:end-1);
    deriv2 = deriv1(2:end) - deriv1(1:end-1);
    thresh = max(deriv2(2:end)) / 2;
    idxes = find(deriv2 > thresh);
    if idxes(1) == 1
        center = (idxes(2) + idxes(end)) / 2;
    else
        center = (idxes(1) + idxes(end)) / 2;
    end
end



