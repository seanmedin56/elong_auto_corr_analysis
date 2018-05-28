function [fps, fes] = extract_features(features, traces, use_errorbars, ...
    max_delay, num_bootstraps)
% Extracts features from traces and puts them into feature_plots (and puts
% errorbars into feature_errors if use_errorbars is true)
%   Detailed explanation goes here

corr = auto_corr_m_calc_norm(traces, max_delay);
    
    if use_errorbars
        alt_cors = cell(1, num_bootstraps);
        for j = 1:num_bootstraps
            subset_idxes = randi([1 length(traces)], 1, length(traces));
            subset = cell(1, length(traces));
            for k = 1:length(traces)
                subset{k} = traces{subset_idxes(k)};
            end
            alt_cors{j} = auto_corr_m_calc_norm(subset, max_delay);
        end
    end
    fps = zeros(1, length(features));
    fes = zeros(1, length(features));
    for f = 1:length(features)
        funct = features{f};
        feature = funct(corr);
        fps(f) = feature;
        
        if use_errorbars
            errors = zeros(1, num_bootstraps);
            for j = 1:num_bootstraps
                errors(j) = funct(alt_cors{j});
            end
            fes(f) = std(errors);
        end
    end
end

