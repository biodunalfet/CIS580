clear all
close all

load('pcd.mat')

figure(1)
X = X3D(logical(ReconX),:);
C_data = Color(logical(ReconX),:);
mask = sqrt(sum(X.^2,2)) < 15 & X(:,3) > 0;
showPointCloud(X(mask,:)*[0 -1 0; 0 0 -1; 1 0 0], uint8(C_data(mask,:)))