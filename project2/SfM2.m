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
match_inliers = cell(size(matches));
tic
n_bytes = 0;
for i = 1:5
    for j = 1:6
        if isempty(matches{i,j})
            continue
        end
        
        %% pull features and color data
        feature_1 = matches{i,j}(:,1:2);
        feature_2 = matches{i,j}(:,3:4);
        color_data = matches{i,j}(:,5:end);
        
        %% RANSAC outlier rejection
        t_s = toc;
        [x1, x2, idx] = GetInliersRANSAC(feature_1,feature_2,0.5,10000);
        t_e = toc;
        %fprintf(repmat('\b',1,n_bytes));
        n_bytes = fprintf('time to RANSAC: %6.6f seconds\n',t_e-t_s);
        color_inliers = color_data(idx,:);
        match_inliers{i,j} = [x1 x2 color_inliers];
    end
end
fprintf('\n')
x1 = match_inliers{1,2}(:,1:2);
x2 = match_inliers{1,2}(:,3:4);
colors = match_inliers{1,2}(:,5:end);
F = EstimateFundamentalMatrix(x1,x2);
E = EssentialMatrixFromFundamentalMatrix(F,K);
[Cset, Rset] = ExtractCameraPose(E);
Xset = cell(4,1);
for i = 1:4
    Xset{i} = LinearTriangulation(K, zeros(3,1), eye(3), Cset{i}, Rset{i}, x1, x2);
end
[C, R, X0] = DisambiguateCameraPose(Cset,Rset,Xset);

%{
mask = X0(:,3) > 0 & sqrt(sum(X0.^2,2)) < 30;
X0 = X0(mask,:);
colors = colors(mask,:);
x1 = x1(mask,:);
x2 = x2(mask,:);
%}

t_s = toc;
X = NonlinearTriangulation(K, zeros(3,1), eye(3), C, R, x1, x2, X0);
t_e = toc;
fprintf('Time to triangulate: %6.6f seconds\n', t_e-t_s);
%{
C_f_set{1} = C;
R_f_set{1} = R;
for i = 3:6
    
    [Cnew Rnew] = PnPRANSAC(X, 
end
%}
% images with overlayed features
figure
showMatchedFeatures(im{1},im{2},x1,x2,'montage')
%{d
figure
imshow(im{1})
hold on
plot(x1(:,1),x1(:,2),'r*')
%}
%{
figure
imshow(im{2})
hold on
plot(x2(:,1),x2(:,2),'r*')
%}
%{d
mask = X(:,3) > 0 & sqrt(sum(X.^2,2)) < 30;
X = X(mask,:);
colors = colors(mask,:);
%}
figure
showPointCloud(X*[0 -1 0; 0 0 -1; 1 0 0], uint8(colors))
axis equal
xlabel('x')
ylabel('y')
zlabel('z')