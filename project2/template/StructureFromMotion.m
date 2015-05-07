%function StructureFromMotion
close all
addpath('..')
addpath('../SBA_example')
K = [568.996140852 0 643.21055941;
     0 568.988362396 477.982801038;
     0 0 1];
nImages = 6;

R_w=[0 -1 0; 0 0 -1; 1 0 0]';
R_c = [0 -1 0; 0 0 -1; 1 0 0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Load images
im = cell(nImages,1);
disp('Loading images...')
for iImage = 1 : nImages
    str = sprintf('../image%07d.bmp', iImage);
    im{iImage} = imread(str);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Load matching
Mx = []; My = []; M = []; Color = [];
disp('Loading matches...')

for iImage = 1 : nImages-1;
    str = sprintf('../matching%d.txt', iImage);
    [mx, my, m, color] = LoadMatching(str, iImage, nImages);
    Mx = [Mx;mx];
    My = [My;my];
    M = [M;m];
    Color = [Color; color];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Initialize 3D points, reconstruction index, and visibility matrix

%load('ransac_data.mat')
X3D = zeros(size(M,1), 3);
ReconX = zeros(size(M,1),1);
V = zeros(size(M,1), nImages);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Exclude outliers using F matrix
disp('Rejecting outliers...')
n_inliers = 0;
for iImage1 = 1 : nImages-1
    for iImage2 = iImage1+1 : nImages
        idx1 = find(M(:,iImage1)==1);
        idx2 = find(M(:,iImage2)==1);
        idx = intersect(idx1, idx2);
        
        x1 = [Mx(idx,iImage1) My(idx,iImage1)];
        x2 = [Mx(idx,iImage2) My(idx,iImage2)];
        if size(x1,1) < 8
            continue;
        end
        [~, ~, inlier] = GetInliersRANSAC(x1, x2,0.005,2000);
        n_inliers = n_inliers + sum(inlier);
        M(idx(~inlier),iImage1) = 0;
    end
end

figure(1)
clf
match_idx = all(logical(M(:,1:2)),2);
showMatchedFeatures(im{1},im{2}, x1,x2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Set initial two frames
initialframe1 = 1;
initialframe2 = 4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Get point index for two frames
idx1 = find(M(:,initialframe1)==1);
idx2 = find(M(:,initialframe2)==1);
idx = intersect(idx1, idx2);

x1 = [Mx(idx,initialframe1) My(idx,initialframe1)];
x2 = [Mx(idx,initialframe2) My(idx,initialframe2)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Get fundamental matrix and essential mtraix
disp('Estimating fundamental and essential matrices...')
F = EstimateFundamentalMatrix(x1, x2);
E = EssentialMatrixFromFundamentalMatrix(F,K);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Camera pose estimation
disp('Getting potential camera poses...')
[Cset, Rset] = ExtractCameraPose(E);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Triangulation and pose disambiguation
disp('Triangulating points for each potential camera pose...')
Xset = cell(4,1);
for i = 1 : 4
    Xset{i} = LinearTriangulation(K, zeros(3,1), eye(3), Cset{i}, Rset{i}, x1, x2);    
end
disp('Picking correct camera pose...')
[C, R, X] = DisambiguateCameraPose(Cset, Rset, Xset);

cam_verts = [0 0 0;1 0.75 0.5;1 0.75 -0.5;1 -0.75 -0.5;1 -0.75 0.5];

faces = [1 2 3 1;...
        2 3 4 5;...
        1 4 5 1;...
        1 4 3 1;...
        1 5 2 1];
    t_1 = -Rset{1}*Cset{1};
    t_2 = -Rset{2}*Cset{2};
    t_3 = -Rset{3}*Cset{3};
    t_4 = -Rset{4}*Cset{4};
   
%{d
figure(7)
clf
mask = sqrt(sum(X.^2,2)) < 30 & X(:,3) > 0;
c_mask = find(idx);
X_r = X(mask,:)*R_w';
plot3(X_r(:,1),X_r(:,2),X_r(:,3),'ro')
grid on
axis equal
hold on
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Nonlinear triangulation
X_lin = X;
disp('Starting nonlinear triangulation');
X = NonlinearTriangulation(K, zeros(3,1), eye(3), C, R, x1, x2, X);

%{d
mask = sqrt(sum(X.^2,2)) < 30 & X(:,3) > 0;
X_r = X(mask,:)*R_w';
plot3(X_r(:,1),X_r(:,2),X_r(:,3),'bo')
p_a = patch('Faces',faces,'Vertices',bsxfun(@plus,cam_verts*R_c',(-R*C)')*R*R_w','FaceColor',[0.5 0.5 0.5]);
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Set reconstructed frame
r_idx = [initialframe1, initialframe2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Set camera pose
Cr_set{1} = zeros(3,1);
Rr_set{1} = eye(3,3);
Cr_set{2} = C;
Rr_set{2} = R;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Set points and visibility matrix
X3D(idx,:) = X;
ReconX(idx) = 1;
V(idx, initialframe1) = 1;
V(idx, initialframe2) = 1;
%% here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Add images
color_scheme = [[0.5 0.5 0.5];
                [1 0 0];
                [0 1 0];
                [0 0 1];
                [1 0.5 0];
                [0 0 0]];
X_reconstructed = cell(5,1);
%{d
figure(10)
clf
patch('Faces',faces,'Vertices',bsxfun(@plus,cam_verts*R_c',(-Rr_set{1}*Cr_set{1})')*Rr_set{1}*R_w','FaceColor',color_scheme(1,:));
hold on
X_r = X*R_w';
mask = sqrt(sum(X_r.^2,2)) < 30 & X_r(:,3) > 0;
plot3(X_r(mask,1),X_r(mask,2),X_r(mask,3),'o','color',color_scheme(1+1,:))
patch('Faces',faces,'Vertices',bsxfun(@plus,cam_verts*R_c',(-Rr_set{2}*Cr_set{2})')*Rr_set{2}*R_w','FaceColor',color_scheme(1+1,:));

axis equal
grid on
drawnow
%}
for iImage = 1 : nImages
    if ~isempty(find(r_idx==iImage,1))
        continue;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
    % Get 2D-3D correspondences
    idx1 = find(ReconX==1);
    idx2 = find(M(:,iImage)==1);
    idx = intersect(idx1, idx2);
    if length(idx) < 6
        continue;
    end
    
    X = X3D(idx,:);
    x = [Mx(idx,iImage) My(idx,iImage)];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
    % Run PnP
    
    C_pnp=[];
    R_pnp=[];
    n_bytes = 0;
    while isempty(C_pnp) || isempty(R_pnp)
        fprintf(repmat('\b',1,n_bytes))
        n_bytes = fprintf('no pose estimate yet...');
        [C_pnp, R_pnp] = PnPRANSAC(X, x, K,10,2000);
    end
    fprintf('Got a pose estimate!\n')
    %{
    figure(11)
    clf
    imshow(im{iImage})
    P_i = K*R_pnp*[eye(3) -C_pnp];
    x_p = bsxfun(@rdivide,P_i(1:2,:)*[X ones(size(X,1),1)]',P_i(3,:)*[X ones(size(X,1),1)]')';
    hold on
    plot(x_p(:,1),x_p(:,2),'ro')
    %}
    disp('Nonlinear PnP');
    [C, R,e_total,residual] = NonlinearPnP(X, x, K, C_pnp, R_pnp);
    
    %{
    P_3_r = K*R*[eye(3) -C];
    x_p_r = bsxfun(@rdivide,P_3_r(1:2,:)*[X ones(size(X,1),1)]',P_3_r(3,:)*[X ones(size(X,1),1)]')';
    plot(x_p_r(:,1),x_p(:,2),'bo')
    plot(x(:,1),x(:,2),'*','color',[1 0.5 0])
    axis equal
    grid on
    xlabel('x')
    ylabel('y')
    zlabel('z')
    legend('Linear PnP Estimate','Non-linear PnP Estimate','Original Features')
    return
    %}
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
    % Set camera poses and reconstructed frame index
    Cr_set{end+1} = C;
    Rr_set{end+1} = R;
    r_idx(end+1) = iImage;
    V(idx, iImage) = 1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
    % Triangulation for additional points
    disp('Adding more points');
    for iImage1 = 1 : length(r_idx)-1
        idx1 = find(ReconX~=1);
        idx2 = find(M(:,r_idx(iImage1))==1);
        idx3 = find(M(:,iImage)==1);
        idx = intersect(idx1, idx2);
        idx = intersect(idx, idx3);
        x1 = [Mx(idx,r_idx(iImage1)) My(idx,r_idx(iImage1))];
        x2 = [Mx(idx,iImage) My(idx,iImage)];
        X = LinearTriangulation(K, Cr_set{iImage1}, Rr_set{iImage1}, C, R, x1, x2);
        X = NonlinearTriangulation(K, Cr_set{iImage1}, Rr_set{iImage1}, C, R, x1, x2, X);
        
        X3D(idx,:) = X;
        ReconX(idx) = 1;
    end
    %{d
    X_r = X3D(logical(V(:,iImage)),:)*R_w';
    mask = sqrt(sum(X_r.^2,2)) < 30 & X_r(:,1) > 0;
    plot3(X_r(mask,1),X_r(mask,2),X_r(mask,3),'o','color',color_scheme(r_idx==iImage,:))
    patch('Faces',faces,'Vertices',bsxfun(@plus,cam_verts*R_c',(-R*C)')*R*R_w','FaceColor',color_scheme(r_idx==iImage,:));
    drawnow
    
    P_temp = K*R*[eye(3) -C];
    X_aug_temp = [X ones(size(X,1),1)];
    x_p = bsxfun(@rdivide,P_temp(1:2,:)*X_aug_temp',P_temp(3,:)*X_aug_temp')';
    %}

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
    % Set visibiltiy and measurements for bundle adjustment
    V_bundle = V(:,r_idx);
    Mx_bundle = Mx(:,r_idx);
    My_bundle = My(:,r_idx);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
    % Run bundle adjustment
    disp('Bundle adjustment');
    [Cr_set, Rr_set, X3D] = BundleAdjustment(K, Cr_set, Rr_set, X3D, ReconX, V_bundle, Mx_bundle, My_bundle);
end

for i = 1:6
    figure(i+18)
    clf
    imshow(im{i})
    P = K*Rr_set{i}*[eye(3) -Cr_set{i}];
    mask = logical(ReconX(:));
    x_p = bsxfun(@rdivide,P(1:2,:)*[X3D(mask,:) ones(size(X3D(mask,:),1),1)]',P(3,:)*[X3D(mask,:) ones(size(X3D(mask,:),1),1)]')';
    mask2 = mask;
    x = [Mx(mask2,i) My(mask2,i)];
    hold on
    plot(x_p(:,1),x_p(:,2),'ro')
    plot(x(:,1),x(:,2),'go')
    legend('Reprojected points','Original features')
end


figure(16)
clf
R_w=[0 -1 0; 0 0 -1; 1 0 0]';
R_c = [0 -1 0; 0 0 -1; 1 0 0];
X3D_sub = X3D(logical(ReconX),:);
mask = sqrt(sum(X3D_sub.^2,2)) < 30 & X3D_sub(:,3) > 0;
Color_recon = Color(logical(ReconX),:);
showPointCloud(X3D_sub(mask,:)*R_w',uint8(Color_recon(mask,:)))
xlabel('x')
ylabel('y')
zlabel('z')
hold on
p_a = patch('Faces',faces,'Vertices',bsxfun(@plus,cam_verts*R_c',(-Rr_set{1}*Cr_set{1})')*Rr_set{1}*R_w','FaceColor',color_scheme(1,:));
p_b = patch('Faces',faces,'Vertices',bsxfun(@plus,cam_verts*R_c',(-Rr_set{2}*Cr_set{2})')*Rr_set{2}*R_w','FaceColor',color_scheme(2,:));
p_c = patch('Faces',faces,'Vertices',bsxfun(@plus,cam_verts*R_c',(-Rr_set{3}*Cr_set{3})')*Rr_set{3}*R_w','FaceColor',color_scheme(3,:));
p_d = patch('Faces',faces,'Vertices',bsxfun(@plus,cam_verts*R_c',(-Rr_set{4}*Cr_set{4})')*Rr_set{4}*R_w','FaceColor',color_scheme(4,:));
p_e = patch('Faces',faces,'Vertices',bsxfun(@plus,cam_verts*R_c',(-Rr_set{5}*Cr_set{5})')*Rr_set{5}*R_w','FaceColor',color_scheme(5,:));
P_f = patch('Faces',faces,'Vertices',bsxfun(@plus,cam_verts*R_c',(-Rr_set{6}*Cr_set{6})')*Rr_set{6}*R_w','FaceColor',color_scheme(6,:));

%save('pcd.mat','X3D','Rr_set','Cr_set','ReconX','Color')