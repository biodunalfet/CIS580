function [C, R] = NonlinearPnP(X,x,K,C0,R0)

opts = optimoptions(@lsqnonlin, 'Algorithm', 'levenberg-marquardt',...
                                'TolX', 1e-200,...
                                'TolFun', 1e-200,...
                                'MaxFunEvals', 1e64,...
                                'MaxIter', 1e64,...
                                'Display','none');

[params,~] = lsqnonlin(@error,[R2q(R0);C0],[],[],opts,X,x,K);
%[params,~] = lsqnonlin(@error,[R0(:);C0],[],[],opts,X,x,K);


R = q2R(params(1:4));
%R = reshape(params(1:end-3),3,3);
C = params(end-2:end);

    function err = error(param0, X, x, K)
        X_aug = [X ones(size(X,1),1)];
        P_n = K*q2R(param0(1:4))*[eye(3) -param0(5:end)];
        %P_n = K*reshape(param0(1:end-3),3,3)*[eye(3) -param0(end-2:end)];
        x_p = bsxfun(@rdivide,P_n(1:2,:)*X_aug',P_n(3,:)*X_aug')';
        err = x-x_p(:,1:2);
    end
end