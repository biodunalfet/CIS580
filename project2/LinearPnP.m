function [C, R] = LinearPnP(X, x, K)
N = length(X);

M = zeros(3*N,12);
for i=1:N
    M(3*(i-1)+1:3*i,:) = cross_mat([X(i,:) 1],x(i,:));
    %{
    M(3*(i-1)+1:3*i,:) = [zeros(1,4) -X(i,:) x(i,2)*X(i,:);...
                            X(i,:) zeros(1,4) -x(i,1)*X;...
                            -x(i,2)*X(i,:) x(i,1)*X(i,:) zeros(1,4)];
    %}
end

[~,~,V] = svd(M);
P = reshape(V(:,end),3,4);
H = K\P;
R = H(:,1:3);
C = -R'*H(:,4);
%{d
    function A = cross_mat(X,x)
        A = [zeros(1,4) -X x(2)*X;...
            X zeros(1,4) -x(1)*X;...
            -x(2)*X x(1)*X zeros(1,4)];
    end
%}
end