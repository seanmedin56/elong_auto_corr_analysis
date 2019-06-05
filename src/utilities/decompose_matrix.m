function [a_array, b_array] = decompose_matrix(R, lambdas)
%DECOMPOSEMATRIX Summary of this function goes here
%   Detailed explanation goes here
    
    lambda_mat = diag(lambdas);
    [V,D] = eig(R);
    V_inv = inv(V);
    a_array = zeros(1, length(lambdas));
    b_array = zeros(1, length(lambdas));
    [~, idx0] = min(abs(diag(D)));
<<<<<<< HEAD
    p = V(:,idx0(1)).^2;
    for i = 1:length(lambdas)
        a_array(i) = lambdas * V(:,i) * V_inv(i,:) * lambda_mat * p;
        a_array(i) = abs(a_array(i) / sum(p .* lambdas'));
=======
    p = V(:,idx0(1));
    for i = 1:length(lambdas)
        a_array(i) = lambdas * V(:,i) * V_inv(i,:) * lambda_mat * p;
>>>>>>> c1b127ab052d404dae55258065f42d213bab972e
        b_array(i) = D(i,i);
    end
end

