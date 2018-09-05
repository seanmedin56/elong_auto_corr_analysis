function [errors] = time_res_function(trace_struct, dt)
%TIME_RES_FUNCTION Summary of this function goes here
%   Detailed explanation goes here
    
    errors = [];
    for i =  1:length(trace_struct)
        times = trace_struct(i).time;
        times_interp = ceil(times(1) / dt):dt:times(end);
        for time = times_interp
            errors = [errors min(abs(time - times))];
        end
    end
    errors = errors / sqrt(length(errors));
end

