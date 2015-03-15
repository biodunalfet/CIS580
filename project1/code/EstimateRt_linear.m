function [ Rs, ts ] = EstimateRt_linear( Hs, K )
% Linear parameter estimation of R and t
%
%   K: a camera calibration matrix. 3 x 3 matrix.
%
%   Hs: a homography from the world to images. 3 x 3 x N matrix, where N is 
%   the number of calibration images. 
%
%   Rs: rotation matrices. 3 x 3 x N matrix, where N is the number of calibration images. 
%
%   ts: translation vectors. 3 x 1 x N matrix, where N is the number of calibration images. 
%

%% Your code goes here.
N = size(Hs,3);
Rs = zeros(size(Hs));
ts = zeros(3,1,N);
for i = 1:N
    z_p = norm(K\Hs(:,1,i));
    r_1 = 1/z_p*(K\Hs(:,1,i));
    r_2 = 1/z_p*(K\Hs(:,2,i));
    r_3 = cross(r_1,r_2);
    ts(:,:,i) = 1/z_p*(K\Hs(:,3,i));
    Rs(:,:,i) = [r_1 r_2 r_3];
end

