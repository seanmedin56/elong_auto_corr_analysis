% script for confirming the elongation time is consistent across
% bootstrapping

% autocorrelaiton parameters
max_delay = 20;
cut = 0;
num_boots = 100;

% range of cleaning/dirtying parameters
thresh = 0;
jump_perc = 99.5;
min_valid = 30;

max_del_z = 4;
max_del_xy = 30;
f = 1;

% imports experimental data to analyze
raw_data = load('../../dat/3_4kb10_5/nucleus_struct.mat');
%raw_data = load('../../dat/2_65kb10_5/nucleus_struct.mat');
%raw_data = load('../../dat/1_9kb10_5/nucleus_struct.mat');

nc = 14;
nucleus_struct = raw_data.nucleus_struct;
nucleus_struct = nucleus_struct([nucleus_struct.setID] >= 0);
subset = nucleus_struct([nucleus_struct.ncStart] == nc);

% figures out time resolution
time_res = nanmedian(diff([subset.time]));

%---------------- generates base set of traces----------------------------

% figures out threshhold for fluo change
jump_thresh = prctile(diff([subset.fluo3D2]), jump_perc);

% interpolates and divides traces
tracesp = {};
idx = 1;
num_big = 0;
for i = 1:length(subset)
    fluo = subset(i).fluo3D2;
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
    tr_del_xy = [0 del_xy(1:end-1) + del_xy(2:end) 0];
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

end

%----------bootstraps traces and records 3rd derivative minimums-----------

real_mins = [];
smooth_mins = [];
for boot = 1:num_boots
    
    sample_idx = randi([1 length(tracesp)], 1, length(tracesp));
    sample = cell([1 length(tracesp)]);
    for j = 1:length(sample)
        sample{j} = tracesp{sample_idx(j)};
    end

    % runs autocorrelation
    auto = cross_corr_m_calc3(sample, sample, max_delay);
    deriv1 = diff(auto);
    deriv2 = diff(deriv1);
    deriv3 = diff(deriv2);

    % smoothing
    deriv3_smooth = deriv3;
    deriv3_smooth(3) = (deriv3(3) * 2 + deriv3(3)) / 3;
    for i = 4:(length(deriv3) - 1)
        deriv3_smooth(i) = 0.25*deriv3(i-1) + ...
            0.5*deriv3(i) + 0.25*deriv3(i+1);
    end
    deriv3_smooth(end) = (deriv3_smooth(end-1) + ...
        deriv3_smooth(end) * 2) / 3;

    % needs offset due to noise correlations
    % created by interpolation
    [~,min_idx] = min(deriv3(3:end));
    min_idx = min_idx + 2;
    real_mins = [real_mins min_idx];

    [~,min_idx_s] = min(deriv3_smooth(3:end));
    min_idx_s = min_idx_s + 2;
    smooth_mins = [smooth_mins min_idx_s];
end

all_poss = sort(unique([real_mins smooth_mins]));
vals = zeros(1, length(all_poss));
vals_s = zeros(1, length(all_poss));
for i = 1:length(all_poss)
    vals(i) = length(find(real_mins == all_poss(i)));
    vals_s(i) = length(find(smooth_mins == all_poss(i)));
end

vals = vals / length(real_mins);
vals_s = vals_s / length(smooth_mins);

figure();
plot(all_poss, vals);
hold on
plot(all_poss, vals_s);
title('distribution of bootstrapping 3rd deriv mins');
xlabel('3rd deriv min');
ylabel('prob');
legend('normal', 'smoothed');
            