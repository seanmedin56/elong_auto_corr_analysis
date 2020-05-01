function [a_array, b_array] = decompose_matrix(R, lambdas)
%DECOMPOSEMATRIX Summary of this function goes here
%   Detailed explanation goes here
    
    lambda_mat = diag(lambdas);
    [V,D] = eig(R);
    V_inv = inv(V);
    a_array = zeros(1, length(lambdas));
    b_array = zeros(1, length(lambdas));
    [~, idx0] = min(abs(diag(D)));
    %p = V(:,idx0(1)).^2; % old code that I'm pretty sure is wrong
    p = V(:,idx0(1)) / sum(V(:,idx0(1)));
    for i = 1:length(lambdas)
        a_array(i) = lambdas * V(:,i) * V_inv(i,:) * lambda_mat * p;
        a_array(i) = abs(a_array(i) / sum(p .* lambdas'));
        b_array(i) = D(i,i);
    end
end

