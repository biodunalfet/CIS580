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
for i = 1:5
    for j = 1:6
        if isempty(matches{i,j})
            continue
        end
        feature_1 = matches{i,j}(:,1:2);
        feature_2 = matches{i,j}(:,3:4);
        color_data = matches{i,j}(:,5:end);
        
        [x1, x2, idx] = GetInliersRANSAC(feature_1,feature_2);
        color_inliers = color_data(idx,:);
        
        %{
        imcat = [im{1} im{2}];
        figure
        imshow(imcat);
        hold on
        for k= 1:size(x1,1)
            plot([x1(k,1) x2(k,1)+1280],[x1(k,2) x2(k,2)],'r*--')
        end
        %}
        %{d
        figure
        imshow(im{1})
        hold on
        plot(x1(:,1),x1(:,2),'r*')

        figure
        imshow(im{2})
        hold on
        plot(x2(:,1),x2(:,2),'r*')
        %}
        F = EstimateFundamentalMatrix(x1,x2);
        E = EssentialMatrixFromFundamentalMatrix(F,K);
        [Cset, Rset] = ExtractCameraPose(E);
        Xset = cell(4,1);
        for k = 1:4
            Xset{k} = LinearTriangulation(K, zeros(3,1), eye(3), Cset{k},Rset{k},x1,x2);
        end
        [C, R, X0] = DisambiguateCameraPose(Cset, Rset, Xset);
        x_norm_mask = sum(X0.^2,2)<100;
        figure
        showPointCloud(([1 0 0; 0 0 -1; 0 1 0]*X0(x_norm_mask,:)')',uint8(color_inliers(x_norm_mask,:)))
        axis equal
        xlabel('x')
        ylabel('y')
        zlabel('z')
        %X = NonlinearTriangulation(K, zeros(3,1), eye(3), C, R, x1, x2, X0);
        
        
    end
end
%{

%}