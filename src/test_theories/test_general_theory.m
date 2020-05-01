addpath('../utilities/');

% establish cases to check
Ts = [4:0.2:20 21:100];
aes = 0:0.02:1;
ces = [0.1, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, 100];
ces = [1];
b_interv = 20;
b_orders = [0.001, 0.01, 0.1, 1, 10, 100];
dists = struct;
idx = 1;

% check all cases for feature
for T = Ts
    for a0 = aes
        for i = 2:length(b_orders)
            bes = linspace(b_orders(i-1),b_orders(i),b_interv);
            for b = bes
                for c = ces
                    auto = zeros(1, ceil(T) + 20);
                    for delay =  0:length(auto)-1
                        [cor_tot, p_term, d_term] = full_func_cor(T,a0,delay,1,b);
                        auto(delay + 1) = c * d_term;
                    end
                    deriv1 = diff(auto);
                    deriv2 = diff(deriv1);
                    deriv3 = diff(deriv2);
                    deriv4 = diff(deriv3);
                    to_break = false;
                    for j = length(deriv4) - 15:length(deriv4)
                        if deriv4(j) <= 0
                            to_break = true;
                            break;
                        end
                    end
                    if to_break
                        continue;
                    end
                    dists(idx).T = T;
                    dists(idx).a0 = a0;
                    dists(idx).b = b;
                    dists(idx).mins = find(islocalmin(deriv3));
                    dists(idx).maxs = find(islocalmax(deriv3));
                    idx = idx + 1;
                end
            end
        end
    end
end

% plot results of check
xes = zeros(1, length(dists));
yes = zeros(1, length(dists));
zes = zeros(1, length(dists));
for i = 1:length(dists)
    dist = dists(i);
    xes(i) = dist.T;
    yes(i) = dist.mins(end) - dist.T;
    if ~isempty(dist.maxs)
        zes(i) = dist.maxs(end) - dist.T;
    else
        zes(i) = nan;
    end
end
figure()
scatter(xes,yes);

figure()
scatter(xes, zes);
