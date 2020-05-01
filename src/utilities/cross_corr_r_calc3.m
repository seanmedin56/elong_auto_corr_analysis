function cross_corr_r = cross_corr_r_calc3(trace1, trace2, max_delay)

% returns the average raw moment crosscorrelation of the traces

% ----------calculates individual correaltions (with weights)--------------

    corr = zeros([1 max_delay]);
    counts = zeros([1 max_delay]);
    for i = 1:length(trace1)
        limit = min([max_delay, length(trace1{i}), length(trace2{i})]);
        len = min(length(trace1{i}), length(trace2{i}));
        weights1 = zeros(1, length(trace1{i}));
        weights2 = zeros(1, length(trace2{i}));
        for j = 1:limit
            temp_mult = trace1{i}(1:len - j + 1) .* trace2{i}(j:len);
            weights1(1:len-j+1) = weights1(1:len-j+1) + ~isnan(temp_mult);
            weights2(j:len) = weights2(j:len) + ~isnan(temp_mult);
        end
        for j = 1:limit
            temp_mult = trace1{i}(1:len - j + 1) .* trace2{i}(j:len);
            temp_mult = temp_mult .* weights1(1:len-j+1) .* weights2(j:len);
            corr(j) = corr(j) + nansum(temp_mult);
            counts(j) = counts(j) + sum(~isnan(temp_mult) .* ...
                weights1(1:len-j+1) .* weights2(j:len));
        end
    end
    
% -----------------combines correaltions together------------------------    
    cross_corr_r = corr ./ counts;
    cross_corr_r = cross_corr_r / cross_corr_r(1);

end

