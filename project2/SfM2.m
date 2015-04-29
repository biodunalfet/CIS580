%clear all
close all

%% load data
imgPath = './';
imgType = '*.bmp';
images = dir([imgPath imgType]);
for idx = 1:length(images)
    im{idx} = imread([imgPath images(idx).name]);
end

K = [568.996140852 0 643.21055941;
     0 568.988362396 477.982801038;
     0 0 1];

matches = cell(5,6);
match_files = dir('matching*');
for i = 1:length(match_files)
    fid = fopen(match_files(i).name);
    first_line = textscan(fgetl(fid),'%s %f');
    num_features = first_line{2};
    data = textscan(fid,repmat('%f ',1,21));
    % loop through each feature in source image
    for j = 1:length(data{1})
        % loop through each match
        for k = 1:data{1}(j) - 1
            % append match to cell array
            target_idx = data{3*(k-1)+7}(j);
            matches{i,target_idx} = [matches{i,target_idx};...
                                    data{5}(j) data{6}(j) data{3*(k-1)+8}(j) data{3*(k-1)+9}(j) data{2}(j) data{3}(j) data{4}(j) j];
        end
    end
    fclose(fid);
end

%% run SfM
match_inliers = cell(size(matches));
ransac_indices = cell(size(matches));
tic
n_bytes = 0;

%% run RANSAC
for i = 1:5
    for j = 1:6
        % skip the dataset if it's empty
        if isempty(matches{i,j})
            continue
        end
        
        % pull features and color data
        feature_1 = matches{i,j}(:,1:2);
        feature_2 = matches{i,j}(:,3:4);
        color_data = matches{i,j}(:,5:7);
        source_index = matches{i,j}(:,8);
        
        % RANSAC outlier rejection
        t_s = toc;
        [x1, x2, idx] = GetInliersRANSAC(feature_1,feature_2,0.05,10000);
        t_e = toc;
        n_bytes = fprintf('time to RANSAC: %6.6f seconds\n',t_e-t_s);
        color_inliers = color_data(idx,:);
        match_inliers{i,j} = [x1 x2 color_inliers];
        ransac_indices{i,j} = source_index(idx);
    end
end
fprintf('\n')

%% get first set of 3d points
x1 = match_inliers{1,2}(:,1:2);
x2 = match_inliers{1,2}(:,3:4);
colors = match_inliers{1,2}(:,5:end);

% calculate fundamental matrix
F = EstimateFundamentalMatrix(x1,x2);

% calculate essential matrix
E = EssentialMatrixFromFundamentalMatrix(F,K);

% get all four possible camera poses
[Cset, Rset] = ExtractCameraPose(E);

% triangulate features
Xset = cell(4,1);
for i = 1:4
    Xset{i} = LinearTriangulation(K, zeros(3,1), eye(3), Cset{i}, Rset{i}, x1, x2);
end

% get correct camera pose
[C, R, X0] = DisambiguateCameraPose(Cset,Rset,Xset);

% run non-linear optimization
X = NonlinearTriangulation(K, zeros(3,1), eye(3), C, R, x1, x2, X0);

%{
% plot triangulated features and their 3d points
figure
showMatchedFeatures(im{1},im{2},x1,x2)
drawnow

figure
mask = sqrt(sum(X.^2,2))<5000 & X(:,3) > 0;
showPointCloud(X(mask,:)*[0 -1 0; 0 0 -1; 1 0 0], uint8(colors(mask,:)))
axis equal
xlabel('x')
ylabel('y')
zlabel('z')
%}

C_set = {C};
R_set = {R};
X_set = X;
for i = 1:1
    for j = 1:6
        if ~isempty(matches{i,j})
            n_empty_idx = j;
            break
        end
    end
    for j = i:3
        if isempty(matches{i,j}) || (i == 1 && j == 2)
            continue
        end
        x2 = match_inliers{i,j}(:,3:4);
        [~,idx1,idx2] = intersect(ransac_indices{1,2},ransac_indices{1,3});

        %{
        figure
        clf
        imshow(im{1})
        hold on
        plot(match_inliers{i,n_empty_idx}(idx1,1),match_inliers{i,n_empty_idx}(idx1,2),'r*')

        figure
        clf
        imshow(im{3})
        hold on
        plot(match_inliers{i,j}(idx2,3),match_inliers{i,j}(idx2,4),'b*')
        %}
        
        [Cnew, Rnew] = PnPRANSAC(X(idx1,:),x2(idx2,:),K,0.1,10000);
        [Cnew, Rnew] = NonlinearPnP(X(idx1,:),x2(idx2,:),K,Cnew,Rnew);
        
        %{d
        figure
        clf
        imshow(im{j})
        hold on
        plot(x2(idx2,1),x2(idx2,2),'b*')
        P = K*Rnew*[eye(3) -Cnew];
        X_aug = [X(idx1,:) ones(length(idx1),1)]';
        x_p = bsxfun(@rdivide,P(1:2,:)*X_aug,P(3,:)*X_aug)';
        plot(x_p(:,1),x_p(:,2),'r*')
        %}
        %{
        x1 = match_inliers{i,j}(:,1:2);
        Xnew = LinearTriangulation(K,C,R,Cnew,Rnew,x1,x2);
        Xnew = NonlinearTriangulation(K,C,R,Cnew,Rnew,x1,x2,Xnew);
        %}
        %{
        % plot triangulated features and their 3d points
        figure
        showMatchedFeatures(im{i},im{j},x1,x2)
        drawnow

        figure
        mask = sqrt(sum(Xnew.^2,2))<500 & Xnew(:,3) > 0;
        showPointCloud(Xnew(mask,:)*[0 -1 0; 0 0 -1; 1 0 0], uint8(colors(mask,:)))
        axis equal
        xlabel('x')
        ylabel('y')
        zlabel('z')
        %}
    end
end