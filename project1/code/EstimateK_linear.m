function [K, Hs] = EstimateK_linear(x, X)
% Linear parameter estimation of K
%
%   x:  2D points. n x 2 x N matrix, where n is the number of corners in
%   a checkerboard and N is the number of calibration images
%       
%   X:  3D points. n x 2 matrix, where n is the number of corners in a
%   checkerboard and assumes the points are on the Z=0 plane
%
%   imgs: calibration images. N x 1 cell, where N is the number of calibration images
%
%   K: a camera calibration matrix. 3 x 3 matrix.
%
%   Hs: a homography from the world to images. 3 x 3 x N matrix, where N is 
%   the number of calibration images. You can use est_homography.m to
%   estimate homographies.
%

%% Your code goes here
N = size(x,3);
Hs = zeros(3,3,N);
v = zeros(2*N,6);
for i = 1:N
    Hs(:,:,i) = est_homography(x(:,1,i),x(:,2,i),X(:,1),X(:,2));
    v(2*i-1,:) = (v_ij(Hs(:,:,i),1,1)-v_ij(Hs(:,:,i),2,2))';
    v(2*i,:) = v_ij(Hs(:,:,i),1,2)';
end

[~,~,V] = svd(v);
B = V(:,end);
B_11 = B(1);
B_12 = B(2);
B_22 = B(3);
B_13 = B(4);
B_23 = B(5);
B_33 = B(6);

p_y = (B_12*B_13 - B_11*B_23)/(B_11*B_22-B_12^2);
c = B_33-(B_13^2+p_y*(B_12*B_13-B_11*B_23))/B_11;
f_y = sqrt(c*B_11/(B_11*B_22-B_12^2));
f_x = sqrt(c/B_11);
s = -B_12*f_x^2*f_y/c;
p_x = s*p_y/f_y-B_13*f_x^2/c;

K = [f_x s p_x; 0 f_y p_y; 0 0 1];

function vec = v_ij(H,i,j)
   vec = [H(1,i)*H(1,j);...
       H(1,i)*H(2,j)+H(2,i)*H(1,j);...
       H(2,i)*H(2,j);...
       H(3,i)*H(1,j)+H(1,i)*H(3,j);...
       H(3,i)*H(2,j)+H(2,i)*H(3,j);...
       H(3,i)*H(3,j)];
end

end