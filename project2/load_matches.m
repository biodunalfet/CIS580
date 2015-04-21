matchPath = './';
matchType = 'matching*';
match_files = dir([imgPath imgType]);
for idx = 1:length(match_files)
    match_data{idx} = fileread([matchPath match_files(idx).name]);
    split_texts = strsplit(text,'\n');
    
    num_features
end

text = fileread(filename);

% split text by line
split_texts = strsplit(text,'\n');

% find uncommented and otherwise valid phrases
valid_boundary = regexp(split_texts,'^\s*boundary(\s*[-]?\d+\.?\d*){6}\s*$');
valid_blocks = regexp(split_texts,'^\s*block(\s*[-]?\d+\.?\d*){9}\s*$');

% get indices of valid phrases
blocks_idx = find(~cellfun(@isempty, valid_blocks));

% extract boundary and obstacle values
boundary_text = split_texts{~cellfun(@isempty,valid_boundary)};
boundary = cell2mat(textscan(boundary_text(length('boundary')+1:end), '%f'))';
blocks = zeros(9,length(blocks_idx));
for i = 1:length(blocks_idx)
    block_text = split_texts{blocks_idx(i)};
    blocks(:,i) = cell2mat(textscan(block_text(length('block')+1:end), '%f'));
end
blocks = blocks';
