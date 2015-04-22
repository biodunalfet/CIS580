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
                                    data{5}(j) data{6}(j) data{3*(k-1)+8}(j) data{3*(k-1)+9}(j) data{2}(j) data{3}(j) data{4}(j)];
        end
    end
    fclose(fid);
end

%% run SfM
X_triangulated = cell(size(matches));
R_measured = cell(size(matches));
C_measured = cell(size(matches));

for i = 1:5
    for j = 1:6
        if isempty(matches{i,j})
            continue
        end
        
        %% pull features and color data
        feature_1 = matches{i,j}(:,1:2);
        feature_2 = matches{i,j}(:,3:4);
        color_data = matches{i,j}(:,5:end);
        
        %% adaptive non-maximal suppression
        %{
        N = size(matches{i,j},1);
        [feature_1,feature_2,anms_idx] = anms(feature_1,feature_2,floor(N*0.9));
        color_data = color_data(anms_idx,:);
        %}
        
        %% RANSAC outlier rejection
        [x1, x2, idx] = GetInliersRANSAC(feature_1,feature_2,0.4,100000);
        color_inliers = color_data(idx,:);

        %% estimate fundamental and essential matrices
        F = EstimateFundamentalMatrix(x1,x2);
        E = EssentialMatrixFromFundamentalMatrix(F,K);
        
        %% estimate possible camera poses and triangulate points
        [Cset, Rset] = ExtractCameraPose(E);
        Xset = cell(4,1);
        for k = 1:4
            Xset{k} = LinearTriangulation(K, zeros(3,1), eye(3), Cset{k},Rset{k},x1,x2);
        end
        
        %% extract proper camera pose
        [C, R, X0] = DisambiguateCameraPose(Cset, Rset, Xset);
        
        %% filter out points behind the camera
        x_norm_mask = sum(X0.^2,2)<1000 & X0(:,3) > 0;
        
        X0 = X0(x_norm_mask,:);
        x1 = x1(x_norm_mask,:);
        x2 = x2(x_norm_mask,:);
        color_inliers = color_inliers(x_norm_mask,:);
        
        %% non-linear triangulation refinement
        X = X0;%NonlinearTriangulation(K, zeros(3,1), eye(3), C, R, x1, x2, X0);
        
        X_triangulated{i,j} = X;
        R_measured{i,j} = R;
        C_measured{i,j} = C;
        
        %% some plotting        
        %{
        % concatenated images
        imcat = [im{1} im{2}];
        figure
        imshow(imcat);
        hold on
        for k= 1:size(x1,1)
            plot([x1(k,1) x2(k,1)+1280],[x1(k,2) x2(k,2)],'r*--')
        end
        %}
        %{
        % images with overlayed features
        figure
        imshow(im{1})
        hold on
        plot(x1(:,1),x1(:,2),'r*')

        figure
        imshow(im{2})
        hold on
        plot(x2(:,1),x2(:,2),'r*')
        %}
        
        %{
        % 3d point cloud
        figure
        showPointCloud(X0*[0 -1 0; 0 0 -1; 1 0 0],uint8(color_inliers))
        axis equal
        xlabel('x')
        ylabel('y')
        zlabel('z')
        %}
        
        %{
        % projected 3d points onto camera view
        projection_1 = K*[eye(3) zeros(3,1)]*[X0 ones(length(X0),1)]';
        scaled_projection_1 = bsxfun(@rdivide,projection_1(1:2,:), projection_1(3,:));
        figure
        plot(scaled_projection_1(1,:),-scaled_projection_1(2,:),'r*')
        axis equal
        xlim([0 1280])
        ylim([-960 0])
        
        projection_2 = K*R*[eye(3) -C]*[X0 ones(length(X0),1)]';
        scaled_projection_2 = bsxfun(@rdivide,projection_2(1:2,:), projection_2(3,:));
        figure
        plot(scaled_projection_2(1,:),-scaled_projection_2(2,:),'r*')
        axis equal
        xlim([0 1280])
        ylim([-960 0])
        %}
        
        %{
        % points tuned with nonlinear optimization
        figure
        showPointCloud(X*[0 -1 0; 0 0 -1; 1 0 0],uint8(color_inliers))
        axis equal
        xlabel('x')
        ylabel('y')
        zlabel('z')
        drawnow
        %}
    end
end

%% move points wrt origin and create point cloud
X_final = [];
for i = 1:5
    for j = 1:6
        if isempty(X_triangulated{i,j})
            continue
        end
        
        if i ~= 1
            X_temp = [X_triangulated{i,j} ones(length(X_triangulated{i,j}),1)]';
            H = eye(4);
            for k = 1:j
                if isempty(R_measured{i,k})
                    continue
                end
                H = H/[R_measured{i,k} -R_measured{i,k}*C_measured{i,k};...
                        0 0 0 1];
            end
            X_moved = H*X_temp;
            X_final = [X_final; X_moved(1:3,:)'];
        else
            X_final = [X_final; X_triangulated{i,j}];
        end
    end
end
figure
showPointCloud(X_final*[0 -1 0; 0 0 -1; 1 0 0])
axis equal
xlabel('x')
ylabel('y')
zlabel('z')