clear all

matches = cell(5,6);

match_files = dir('matching*');
for i = 1:length(match_files)
    fid = fopen(match_files(i).name);
    first_line = textscan(fgetl(fid),'%s %f');
    num_features = first_line{2};
    %data = textscan(fid,'%d %d %d %d %f %f %d %f %f %d %f %f %d %f %f %d %f %f %d %f %f');
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
