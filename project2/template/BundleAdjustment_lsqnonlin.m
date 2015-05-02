function [Cset, Rset, Xset] = BundleAdjustment_lsqnonlin(K, Cr_set, Rr_set, X3D, ReconX, V_bundle, Mx_bundle, My_bundle)
opts = optimoptions(@lsqnonlin, 'Algorithm', 'levenberg-marquardt',...
                                'TolX', 1e-64,...
                                'TolFun', 1e-64,...
                                'MaxFunEvals', 1e64,...
                                'MaxIter', 1e64,...
                                'Display','iter');
n_frames = length(Rr_set);
mask = logical(ReconX);
X = X3D(mask,:);
n_features = size(X,1);
Mx = Mx_bundle(mask,:);
My = My_bundle(mask,:);
V = V_bundle(mask,:);

params0 = zeros(7*n_frames+numel(X),1);
for i = 1:n_frames
    q = R2q(Rr_set{i});
    params0((i-1)*7+1:i*7) = [q;Cr_set{i}];
end
params0(7*n_frames+1:end) = reshape(X,numel(X),1);

[params,~] = lsqnonlin(@error,params0,[],[],opts,Mx,My,V,K,n_frames, n_features);

Cset = cell(size(Cr_set));
Rset = cell(size(Rr_set));
for i = 1:n_frames
    Rset{i}=q2R(params(i-1)*7+1:i*7-3);
    Cset{i}=params((i-1)*7+5:i*7);
end
Xset = reshape(params(7*n_frames+1:end),n_features,3);

    function err = error(param, Mx, My, V, K, n_frames, n_features)
        e_p = cell(1,n_frames);
        X_aug = [reshape(param(7*n_frames+1:end),n_features,3) ones(n_features,1)];
        for iter = 1:n_frames
            P = K*q2R(param((iter-1)*7+1:iter*7-3))*[eye(3) -param((iter-1)*7+5:iter*7)];
            e_p{iter} = bsxfun(@times,bsxfun(@rdivide, P(1:2,:)*X_aug',P(3,:)*X_aug')' - [Mx(:,iter) My(:,iter)],V(:,iter));
        end
        err = cell2mat(e_p);
        err = err(:);
    end
end