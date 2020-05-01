addpath('../utilities/');

% this script tests the consistency of the appearance of key features when
% I chance the clean up parameters

% imports data
data_set = load('../dat/3_4kb_highthresh/nucleus_struct.mat');
nucleus_struct = data_set.nucleus_struct;
subset = nucleus_struct([nucleus_struct.ncStart] == 14);

% clean up parameters that I'm not testing
max_delay = 25;
field = 'fluo3D2_interp';
skip = 1; % downsampling if greater than 1 (must be an integer > 1)

% clean up parameters to span
cuts = 0; % how many beginning points should be cut
min_valids = 80; % minimum # of nonzero values each trace must have
min_kills = 20:5:40; % # points needed between values below min_val in order to keep segment
min_vals = 0:200000:800000; % see above

results = struct;
res_idx = 1;
for cut = cuts
    for min_valid = min_valids
        for min_kill = min_kills
            for min_val = min_vals
                results(res_idx).cut = cut;
                results(res_idx).min_valid = min_valid;
                results(res_idx).min_kill = min_kill;
                results(res_idx).min_val = min_val;
                
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
                                if splits(i + 1) - splits(i) <= min_kill
                                    trace(splits(i):splits(i+1)) = 0;
                                end
                            end
                            min_nonzero = find(trace > min_val, 1, 'first');
                            max_nonzero = find(trace > min_val, 1, 'last');
                            trace = trace(min_nonzero + cut:max_nonzero - cut);
                            %trace = trace - mean(trace);
                            %trace = trace(2:end) - trace(1:end- 1);
                            len = length(trace);
                            if length(find(trace > 0)) >= min_valid
                                for k = 1:skip
                                    traces{idx} = trace(k:skip:len);
                                    idx = idx + 1;
                                end
                            end
                        end
                    end
                end
                
                % finds autocorrelation, takes derivatives, extracts
                % features
                corr = auto_corr_m_calc_norm(traces, max_delay);
                deriv1 = corr(2:end) - corr(1:end-1);
                deriv2 = deriv1(2:end) - deriv1(1:end-1);
                deriv3 = deriv2(2:end) - deriv2(1:end-1);
                [~, I] = min(deriv1);
                results(res_idx).d1_min = I(1);
                [~, I] = max(deriv2);
                results(res_idx).d2_max  = I(1);
                [~, I] = min(deriv3);
                results(res_idx).d3_min = I(1);
                res_idx = res_idx + 1;
            end
        end
    end
end

% plots possibly useful information
figure()
histogram([results.d1_min]);
title('deriv1 min');

figure()
histogram([results.d2_max]);
title('deriv2 max');

figure()
histogram([results.d3_min]);
title('deriv3 min');


