function [x1_nms, x2_nms, idx] = anms(x1,x2,num_feat)
    per_list = floor(num_feat/2);
    sup_list_1 = zeros(length(x1),1);
    sup_list_2 = zeros(length(x2),1);
    for i = 1:length(x1)
        [sup_list_1(i,1), sup_list_1(i,2)] = min(sum(bsxfun(@minus,x1(i,:),x1).^2,2));
        [sup_list_2(i,1), sup_list_2(i,2)] = min(sum(bsxfun(@minus,x2(i,:),x2).^2,2));
    end
    [~,sort_idx_1] = sort(sup_list_1(:,1));
    [~,sort_idx_2] = sort(sup_list_2(:,1));
    
    picked_indices_1 = sup_list_1(sort_idx_1,2);
    picked_indices_2 = sup_list_2(sort_idx_2,2);
    idx = unique([picked_indices_1(1:per_list);picked_indices_2(1:per_list)]);
    x1_nms = x1(idx,:);
    x2_nms = x2(idx,:);
end