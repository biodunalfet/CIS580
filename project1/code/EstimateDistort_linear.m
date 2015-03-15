function [ks] = EstimateDistort_linear(x, X, K, Rs, ts)
% Linear parameter estimation of k
%
%   x:  2D points. n x 2 x N matrix, where n is the number of corners in
%   a checkerboard and N is the number of calibration images
%       
%   X:  3D points. n x 2 matrix, where n is the number of corners in a
%   checkerboard and assumes the points are on the Z=0 plane
%
%   K: a camera calibration matrix. 3 x 3 matrix.
%
%   Rs: rotation matrices. 3 x 3 x N matrix, where N is the number of calibration images. 
%
%   ts: translation vectors. 3 x 1 x N matrix, where N is the number of calibration images. 
%
%   ks: radial distortion parameters. 2 x 1 matrix, where ks(1) = k_1 and
%   ks(2) = k_2.
%

%% Your code goes here
N = size(x,3);
n = size(x,1);
r_mat = zeros(2*N*n,2);
K_p = K;
K_p(1,2) = 0;
pix_ideal = zeros(2*N*n,1);
pix_img = pix_ideal;
for i = 1:N
    ab_h = [Rs(:,:,i) ts(:,:,i)]*[X';zeros(1,n);ones(1,n)];
    ab_ih = bsxfun(@rdivide,ab_h(1:2,:),ab_h(3,:));
    r = repmat(sqrt(sum(ab_ih.^2)),2,1);
    r_mat((i-1)*2*n+1:i*2*n,:) = [reshape(r.^2,2*n,1) reshape(r.^4,2*n,1)];
    
    pix_i = K_p*[ab_ih; ones(1,n)];
    pix_ideal((i-1)*2*n+1:i*2*n) = reshape(bsxfun(@rdivide,pix_i(1:2,:),pix_i(3,:)), 2*n,1);
    
    % store actual corner pixel coordinates
    pix_img((i-1)*2*n+1:i*2*n,:) = reshape(x(:,:,i)', 2*n, 1);
end

c_vec = repmat(K(1:2,3),N*n,1);
M_l = repmat(pix_ideal-c_vec,1,2).*r_mat;
M_r = pix_img-pix_ideal;
ks = M_l\M_r;
end