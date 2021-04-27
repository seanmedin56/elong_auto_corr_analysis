% script for confirming the elongation time can still be successfully
% picked out with my data cleaning methods

% autocorrelaiton parameters
max_delay = 20;
cut = 0;

% range of cleaning/dirtying parameters
threshes = 0;
jump_percs = [99, 99.25, 99.5, 99.7, 99.9, 99.99];
min_valids = 25:5:80;
data_type = 'fluo3DN'; % options: 'fluo3D2', 'fluo3DN', 'fluo3DNRaw'

max_del_z = 3;
max_del_xy = 6;
f = 1.5;

% imports experimental data to analyze
%raw_data = load('../../dat/3_4kb_21_04_24/nucleus_struct.mat');
%raw_data = load('../../dat/2_65kb_21_04_24/nucleus_struct.mat');
raw_data = load('../../dat/1_9kb_21_04_24/nucleus_struct.mat');
%raw_data = load('../../dat/3_4kb_9_28_20/nucleus_struct.mat');
%raw_data = load('../../dat/2_65kb_9_28_20/nucleus_struct.mat');
%raw_data = load('../../dat/1_9kb_9_28_20/nucleus_struct.mat');

nc = 14;
nucleus_struct = raw_data.nucleus_struct;
nucleus_struct = nucleus_struct([nucleus_struct.setID] >= 0);
subset = nucleus_struct([nucleus_struct.ncStart] == nc);

% figures out time resolution
time_res = nanmedian(diff([subset.time]));

all_data = struct;

% loop through combinations
data_idx = 1;
% loops through combinations of cleaning/dirtying parameters
for thresh = threshes
    for jump_perc = jump_percs
        for min_valid = min_valids
            num_big = 0;
            % figures out threshhold for fluo change
            jump_thresh = prctile(diff([subset.(data_type)]), jump_perc);

            % interpolates and divides traces
            tracesp = {};
            tracesp2 = {};
            idx = 1;
            idx2 = 1;
            for i = 1:length(subset)
                fluo = subset(i).(data_type);
                time = subset(i).time;
                xPos = subset(i).xPosParticle;
                yPos = subset(i).yPosParticle;
                zPos = subset(i).zPosParticle;

                % figures out timing for interpolation
                t_start = time(1);
                t_end = round((time(end) - t_start) ...
                    / time_res) * time_res + t_start;
                time_interp = t_start:time_res:t_end;

                % nans too big movements
                del_xy = sqrt(diff(xPos).^2 + diff(yPos).^2);
                del_z = abs(diff(zPos));
                tr_del_xy = [0 max(del_xy(1:end-1), del_xy(2:end)) 0];
                tr_del_z1 = [0 del_z(1:end-1) 0];
                tr_del_z2 = [0 del_z(2:end) 0];
                fluo(tr_del_xy >= max_del_xy) = NaN;
                fluo(tr_del_z1 >= max_del_z & tr_del_z2 >= max_del_z) = NaN;   

                tr_dd1 = abs([0 diff(diff(fluo)) 0]);
                fluo(tr_dd1>jump_thresh) = NaN;
                num_big = num_big + length(find(tr_dd1 > jump_thresh));

                % interpolates
                non_nans = find(~isnan(fluo));
                if length(non_nans) < min_valid
                    continue
                end
                trace = interp1(time(non_nans), ...
                    fluo(non_nans), time_interp);

                % nans below thresh
                trace(trace < thresh) = NaN;

                % nans anything not close to known point
                new_times = time(non_nans);
                for j = 1:length(time_interp)
                    higher = find(new_times > time_interp(j));
                    lower = find(new_times < time_interp(j));
                    if isempty(higher) || isempty(lower)
                        continue
                    end
                    dist_high = new_times(higher(1)) - time_interp(j);
                    dist_low = time_interp(j) - new_times(lower(end));
                    if dist_low > f*time_res && dist_high > f*time_res
                        trace(j) = NaN;
                    end
                end

                % divide into segments
                bad_pts = [1 find(isnan(trace)) length(trace)];
                for j = 2:length(bad_pts)
                    if bad_pts(j) - bad_pts(j-1) > min_valid
                        new_tr = trace(bad_pts(j-1)+1:bad_pts(j)-1);
                        tracesp{idx} = new_tr;
                        idx = idx + 1;
                    end
                end
                
                % processes traces for the nan autocorrelation methods
                trace = interp1(time, fluo, time_interp);

                % nans below thresh
                trace(trace < thresh) = NaN;

                % nans anything not close to known point
                new_times = time(non_nans);
                for j = 1:length(time_interp)
                    higher = find(new_times > time_interp(j));
                    lower = find(new_times < time_interp(j));
                    if isempty(higher) || isempty(lower)
                        continue
                    end
                    dist_high = new_times(higher(1)) - time_interp(j);
                    dist_low = time_interp(j) - new_times(lower(end));
                    if dist_low > f*time_res && dist_high > f*time_res
                        trace(j) = NaN;
                    end
                end

                % divide into segments
                bad_pts = [1 find(isnan(trace)) length(trace)];
                for j = 2:length(bad_pts)
                    if bad_pts(j) - bad_pts(j-1) <= min_valid
                        trace(bad_pts(j-1)+1:bad_pts(j)-1) = nan;

                    end
                end
                if length(trace) >= min_valid
                    tracesp2{idx} = trace;
                    idx2 = idx2 + 1;
                end

            end

            %--------- runs first autocorrelation-------------------%
            auto = cross_corr_m_calc(tracesp, tracesp, max_delay);
            deriv3 = diff(diff(diff(auto)));

            % smoothing
            deriv3_smooth = deriv3;
            for i = 4:(length(deriv3) - 1)
                deriv3_smooth(i) = 0.25*deriv3(i-1) + ...
                    0.5*deriv3(i) + 0.25*deriv3(i+1);
            end
            deriv3_smooth(end) = (deriv3_smooth(end-1) + ...
                deriv3_smooth(end) * 2) / 3;

            % needs offset due to noise correlations
            % created by interpolation
            [~,r1_min_idx] = min(deriv3(3:end));
            r1_min_idx = r1_min_idx + 2;

            [~,r1_min_idx_s] = min(deriv3_smooth(3:end));
            r1_min_idx_s = r1_min_idx_s + 2;
            
            %--------- runs second autocorrelation-------------------%
            auto = cross_corr_m_calc2(tracesp2, tracesp2, max_delay);
            deriv3 = diff(diff(diff(auto)));

            % smoothing
            deriv3_smooth = deriv3;
            for i = 4:(length(deriv3) - 1)
                deriv3_smooth(i) = 0.25*deriv3(i-1) + ...
                    0.5*deriv3(i) + 0.25*deriv3(i+1);
            end
            deriv3_smooth(end) = (deriv3_smooth(end-1) + ...
                deriv3_smooth(end) * 2) / 3;

            % needs offset due to noise correlations
            % created by interpolation
            [~,r2_min_idx] = min(deriv3(3:end));
            r2_min_idx = r2_min_idx + 2;

            [~,r2_min_idx_s] = min(deriv3_smooth(3:end));
            r2_min_idx_s = r2_min_idx_s + 2;
            
            %--------- runs third autocorrelation-------------------%
            auto = cross_corr_m_calc3(tracesp2, tracesp2, max_delay);
            deriv3 = diff(diff(diff(auto)));

            % smoothing
            deriv3_smooth = deriv3;
            for i = 4:(length(deriv3) - 1)
                deriv3_smooth(i) = 0.25*deriv3(i-1) + ...
                    0.5*deriv3(i) + 0.25*deriv3(i+1);
            end
            deriv3_smooth(end) = (deriv3_smooth(end-1) + ...
                deriv3_smooth(end) * 2) / 3;

            % needs offset due to noise correlations
            % created by interpolation
            [~,r3_min_idx] = min(deriv3(3:end));
            r3_min_idx = r3_min_idx + 2;

            [~,r3_min_idx_s] = min(deriv3_smooth(3:end));
            r3_min_idx_s = r3_min_idx_s + 2;
            
            %--------- runs fourth autocorrelation-------------------%
            auto = cross_corr_m_calc4(tracesp, tracesp, max_delay);
            deriv3 = diff(diff(diff(auto)));

            % smoothing
            deriv3_smooth = deriv3;
            for i = 4:(length(deriv3) - 1)
                deriv3_smooth(i) = 0.25*deriv3(i-1) + ...
                    0.5*deriv3(i) + 0.25*deriv3(i+1);
            end
            deriv3_smooth(end) = (deriv3_smooth(end-1) + ...
                deriv3_smooth(end) * 2) / 3;

            % needs offset due to noise correlations
            % created by interpolation
            [~,r4_min_idx] = min(deriv3(3:end));
            r4_min_idx = r4_min_idx + 2;

            [~,r4_min_idx_s] = min(deriv3_smooth(3:end));
            r4_min_idx_s = r4_min_idx_s + 2;

            all_data(data_idx).time_res = time_res;
            all_data(data_idx).thresh = thresh;
            all_data(data_idx).jump_perc = jump_perc;
            all_data(data_idx).min_valid = min_valid;
            all_data(data_idx).r1_min3 = r1_min_idx;
            all_data(data_idx).r1_min3_s = r1_min_idx_s;
            all_data(data_idx).r2_min3 = r2_min_idx;
            all_data(data_idx).r2_min3_s = r2_min_idx_s;
            all_data(data_idx).r3_min3 = r3_min_idx;
            all_data(data_idx).r3_min3_s = r3_min_idx_s;
            all_data(data_idx).r4_min3 = r4_min_idx;
            all_data(data_idx).r4_min3_s = r4_min_idx_s;
            all_data(data_idx).num_traces = length(tracesp);
            all_data(data_idx).num_traces2 = length(tracesp2);
            data_idx = data_idx + 1;
        end
    end
end     

save('1_9kb_fluo3DN_04_24_21.mat', 'all_data');

