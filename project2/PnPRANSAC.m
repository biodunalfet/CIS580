function [C,R] = PnPRANSAC(X, x, K, threshold, max_iter)
N = size(X,1);
R = [];
C = [];
if N < 6
    return
end
X_aug = [X ones(N,1)];
n_best = 0;

median_err = zeros(max_iter,1);
for i = 1:max_iter
    r_idx = ceil(rand(6,1)*N);
    [C_temp, R_temp] = LinearPnP(X(r_idx,:),x(r_idx,:),K);

    P = K*R_temp*[eye(3) -C_temp];
    x_p = bsxfun(@rdivide,P(1:2,:)*X_aug',P(3,:)*X_aug')';
    
    mask = sqrt(sum((x-x_p).^2,2)) < threshold;
    median_err(i) = median(sqrt(sum((x-x_p).^2,2)));
    num_in = sum(mask);
    if n_best < num_in
        n_best = num_in;
        C = C_temp;
        R = R_temp;
        inliers = mask(:);
    end
end

[C, R] = LinearPnP(X(inliers,:),x(inliers,:),K);
%{
figure(200)
clf
hist(median_err,1000)
%}