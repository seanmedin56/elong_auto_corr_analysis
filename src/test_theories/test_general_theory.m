addpath('../utilities/');

% % establish cases to check
Ts = [4:0.2:60];
%Ts = 6:0.2:20;
%Ts = [4:100];
kons = [0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, ...
    10, 20, 50, 100];
koffs = [0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, ...
    10, 20, 50, 100];
koffs = [0.001];
aes = 0:0.02:1;
% ces = [0.1, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, 100];
% b_interv = 20;
% b_orders = [0.001, 0.01, 0.1, 1, 10, 100];
dists = struct;
idx = 1;

% check all cases for feature
for T = Ts
    for a0 = aes
        for kon = kons
            for koff = koffs
                auto = zeros(1, ceil(T) + 5);
                c = koff / sqrt(2 * (kon^2 + koff^2));
                b = kon + koff;
                for delay =  0:length(auto)-1
                    [cor_tot, p_term, d_term] = full_func_cor(T,a0,delay,1,b);
                    auto(delay + 1) = p_term + c * d_term;
                end

                deriv1 = diff(auto);
                deriv2 = diff(deriv1);
                deriv3 = diff(deriv2);

                dists(idx).T = T;
                dists(idx).a0 = a0;
                dists(idx).b = b;
                dists(idx).kon = kon;
                dists(idx).koff = koff;
                dists(idx).mins = find(islocalmin(deriv3));
                [~, real_min] = min(deriv3(4:end));
                dists(idx).abs_min = real_min + 3;
                if isempty(dists(idx).mins)
                    stop = true;
                end
                dists(idx).maxs = find(islocalmax(deriv3));
                idx = idx + 1;
            end
        end
    end
end

% plot results of check
xes = zeros(1, length(dists));
yes = zeros(1, length(dists));
yes2 = zeros(1, length(dists));
zes = zeros(1, length(dists));

for i = 1:length(dists)
    dist = dists(i);
    xes(i) = dist.T;
    if ~isempty(dist.mins)
        yes(i) = dist.mins(end) - dist.T;
        yes2(i) = dist.abs_min - dist.T;
    else
        yes(i) = 5;
    end
    if ~isempty(dist.maxs)
        zes(i) = dist.maxs(end) - dist.T;
    else
        zes(i) = nan;
    end
end

figure()
scatter(xes,yes);

figure()
scatter(xes,yes2);

figure()
scatter(xes, zes);
