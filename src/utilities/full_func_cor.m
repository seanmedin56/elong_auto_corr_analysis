function cor_tot = full_func_cor(elong,alph,tau,aes,bes)
% calculates autocorrelation at time delay tau
%   elong: elongatio time
%   alph: rise time
%   tau: time delay
%   aes: coefficiencts for exponential terms
%   bes: decay values for exponential terms

    alph = min(alph,elong);
    
    %calculates poisson term

    cor_tot = max(0,elong - tau - alph) + (min(alph,max(elong - tau, 0))^2 ...
        - max(alph - tau,0)^2) / 2 / alph + max(alph - tau,0)^3 ...
        / 3 / alph^2 + max(alph - tau,0)^2 * tau / 2 / alph^2;
    
    %calculates dynamics terms
    
    for i = 1:length(aes)
        a = aes(i);
        b = bes(i);
        cor = 0;
        if elong > alph

            %t > alph and t' > alph

            % t' + tau - t > 0
            % part 1
            if elong > alph + tau
                cor = cor + (elong - alph - tau) / b;
                cor = cor + 1 / b^2;
            else
                cor = cor + exp(-b * (alph + tau - elong)) / b^2;
            end
            if tau > alph
                if tau > elong
                    cor = cor + (exp(-b * tau) - exp(-b * (alph + tau - elong))) / b^2;
                else
                    cor = cor + (exp(-b * elong) - exp(-b * alph)) / b^2;
                end
            else
                cor = cor + (exp(-b * (elong + tau - alph)) - exp(-b * tau)) / b^2;
            end
            cor = cor - exp(-b * tau) / b^2;

            %part 2

            if tau > alph
                if elong > tau
                    cor = cor + (exp(-b* alph) - exp(-b * tau) + ...
                        exp(-b * (elong + tau - alph)) - exp(-b * elong)) / b^2;
                else
                    cor = cor + (exp(-b* (alph + tau - elong)) - exp(-b * tau) + ...
                        exp(-b * (elong + tau - alph)) - exp(-b * tau)) / b^2;
                end
            end

            %t' + tau - t < 0

            if alph + tau < elong
                cor = cor + (elong - alph - tau) / b;
                cor = cor + (exp(b * (alph + tau - elong)) - 1) / b^2;
            end

            %t > alph and t' < alph
            if alph > 0.1
                %t' + tau - t > 0


                piece1 = min(max(alph,tau),elong);
                piece2 = min(alph + tau,elong);

                if alph < tau
                    cor = cor + (exp(-b * (tau - piece1)) - exp(-b * (tau - alph))) ...
                    / (alph * b^3);
                end
                cor = cor + (exp(-b * tau) - exp(-b * (alph + tau - piece2))) ...
                    / b^2 + (exp(-b * tau) - ...
                    exp(-b * (alph + tau - piece2))) / b^3 / alph + (piece2 - piece1) / b^2 / alph + ...
                    (piece2^2 - piece1^2 + 2*tau * piece1 - 2*tau * piece2) / (2 * alph * b);

                % t' + tau - t < 0
                if elong > tau
                    cor = cor + (exp(b * (tau - piece1)) - exp(b * (tau - elong))) ...
                        / (alph * b^3) - (piece2 - piece1) / (alph * b^2) + ...
                        (piece2^2 - piece1^2) / (2 * alph * b) + tau * ...
                        (piece1 - piece2) / (alph * b);
                    if elong > tau + alph
                        cor = cor + (1/b^2 - 1 / alph / b^3) * ...
                            (exp(b * (alph + tau - piece2)) - exp(b * (alph + tau - elong)));
                    end
                end

                % t < alph and t' > alph

                % t' + tau - t > 0

                cor = cor + (exp(-b * tau) - exp(-b * (elong + tau - alph))) / b^2 ...
                    + (exp(-b * (elong + tau - alph)) + exp(-b * (alph + tau)) - ...
                    exp(-b * (elong + tau)) - exp(-b * tau)) / (alph * b^3);
            end
        end

        % t < alph and t' < alph

        if alph > 0.1 || alph >= elong
            alph = min(alph,elong);
            %t' + tau - t > 0

            piece3 = min(tau,alph);
            piece4 = max(tau - alph, 0);

            cor = cor + ((-1 / b^2) + (-1 / (b^3 * alph))) * exp(-b * tau) + ...
                ((1 / (b^3 * alph)) + (1 / (b^4 * alph^2))) * ...
                (exp(-b * tau) - exp(-b * (alph + tau))) + piece3 * ...
                exp(-b * piece4) / (b^3 * alph^2) + (exp(-b * tau) - ...
                exp(-b * piece4)) / (b^4 * alph^2) + 1 / (2 * b^2) - ...
                piece3^2 / (2 * b^2 * alph^2) + alph / (3 * b) - tau / (2 * b) ...
                + piece3^2 * tau / (2 * b * alph^2) - piece3^3 / (3 * b * alph^2);

            %t' + tau - t < 0
            if tau < alph
                cor = cor - exp(b * (tau - alph)) / b^3 / alph + piece3 * ...
                    exp(b * (tau - piece3)) / b^3 / alph^2 + (exp(b * (tau - piece3)) ...
                    - exp(b * (tau - alph))) / b^4 / alph^2 - 1 / 2 / b^2 + ...
                    piece3^2 / 2 / b^2 / alph^2 + alph / 3 / b - tau / 2 / b + ...
                    piece3^2 * tau / 2 / b / alph^2 - piece3^3 / 3 / b / alph^2;
            end

        end

        if (isnan(cor) || cor == inf || cor == -inf)
            debug = true;
        end
        cor_tot = cor_tot + a * cor;
    end

end

