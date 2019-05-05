function [fitted_trace, x] = fit_bursts(trace, elong_time, rise_time)
%FIT_BURSTS Finds best fit bursts using nonnegative least squares
%   trace: array of transcription spots at consecutive time steps
%   elong_time: how many time steps it takes for the polymerase to fall off
%   the gene (must be an integer)
%   rise_time: how many time steps it takes for the polymerase to fully
%   transcribe the MS2 (or PP7) loops (probably doesn't have to be an
%   integer)
    num_points = length(trace);
    A = zeros(num_points,num_points);
    for i = 1:(num_points)
        for j = 1:elong_time
            if i + j - 1 > num_points
                break
            end
            if j <= rise_time
                A(i + j - 1, i) = j /( 2 * rise_time);
            else
                A(i + j - 1, i) = 1;
            end
        end
    end
    x = lsqnonneg(A, trace');   
    fitted_trace  = A * x;
end

