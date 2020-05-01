data_set = load('../dat/3_4kb10_5/nucleus_struct.mat');
nucleus_struct = data_set.nucleus_struct;
subset = nucleus_struct([nucleus_struct.ncStart] == 14);

% clean up parameters that I'm not testing
max_delay = 20;
field = 'fluo3D2_interp';
skip = 1; % downsampling if greater than 1 (must be an integer > 1)

% clean up parameters
cut = 0; % how many beginning points should be cut
min_valid = 40; % minimum # of nonzero values each trace must have
min_val = 80000; % see above

% processes traces
traces = {};
idx = 1;
for s = subset
    s.(field)(isnan(s.(field))) = 0;
    if sum(isnan(s.(field))) == 0 && s.inference_flag
        len = length(s.(field));
        trace = s.(field);
        if len >= min_valid
            splits = [1 find(trace < min_val) len];
            for i = 1:length(splits) - 1
                if splits(i + 1) - splits(i) > min_valid
                    new_trace = trace(splits(i)+1:splits(i+1)-1);
                    new_trace = new_trace(1+cut:end-cut);
                    traces{idx} = new_trace;
                    idx = idx + 1;
                end
            end
        end
    end
end

% bootstrap autocorrelation
num_times = 100;
mins = zeros(1, num_times);
for time = 1:num_times
    sample_idx = randi([1 length(traces)], 1, length(traces));
    sample1 = cell([1 length(traces)]);
    sample2 = cell([1 length(traces)]);
    for j = 1:length(sample1)
        sample1{j} = traces{sample_idx(j)};
        sample2{j} = traces{sample_idx(j)};
    end
    corr = cross_corr_m_calc(sample1, sample2, max_delay);
    deriv1 = corr(2:end) - corr(1:end-1);
    deriv2 = deriv1(2:end) - deriv1(1:end-1);
    deriv3 = deriv2(2:end) - deriv2(1:end-1);
    
%     % determines min_elong
%     [~, I] = max(deriv2);
%     min_idx = I - 1;
%     
%     % determines max_elong
%     crossings = find(deriv3(min_idx:end) > 0) + min_idx - 1;
%     max_idx = length(deriv3);
%     for i = 1:length(crossings)
%         if deriv2(crossings(i) + 1) / deriv2(crossings(i)) > 1.5
%             max_idx = crossings(i);
%         end
%     end
%     [~, I] = min(deriv3(min_idx:max_idx));
%     
%     mins(time) = I + min_idx - 1;

    [~, I] = min(deriv3);
    mins(time) = I;
end

figure()
histogram(mins)