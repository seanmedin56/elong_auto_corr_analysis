% this script tests the consistency of the appearance of key features when
% I chance the clean up parameters

% imports data
data_set = load('../dat/1_9kb10_5/nucleus_struct.mat');
nucleus_struct = data_set.nucleus_struct;
subset = nucleus_struct([nucleus_struct.ncStart] == 14);

% clean up parameters that I'm not testing
max_delay = 20;
field = 'fluo3D2_interp';
skip = 1; % downsampling if greater than 1 (must be an integer > 1)

% clean up parameters to span
min_valids = 40; % minimum # of nonzero values each trace must have
min_vals = 0:100000:800000; % see above

num_times = 100;

dips = cell(length(min_valids), length(min_vals));
num_traces = cell(length(min_valids), length(min_vals));

for ii = 1:length(min_valids)
    min_valid = min_valids(ii);
    for ij = 1:length(min_vals)
        min_val = min_vals(ij);

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
            [~, I] = min(deriv3);
            mins(time) = I(1);
        end
        num_traces{ii,ij} = length(traces);
        dips{ii,ij} = mins;
    end
end
