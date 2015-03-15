clear all
close all
clc

%% 2.1
H = [0.0002 0 -0.1482;...
    0 0 0.2983;...
    0 -0.0002 0.9429];
r1_unscaled = H(:,1);
lambda2 = 1/(sqrt(r1_unscaled'*r1_unscaled))

p_1 = [1201 2127];
p_2 = [1539 2127];
p_3 = [939.5 1799];
p_4 = [2484 1805];

X_1 = [0 0 1]';
X_2 = [510 0 1]';
X_p_1 = lambda2*H*X_1
X_p_2 = lambda2*H*X_2
lambda1 = (X_p_1(3)+X_p_2(3))/2
[X_p_1(1)  X_p_1(3) 0;...
          X_p_1(2)  0  X_p_1(3);...
          X_p_2(1)  X_p_2(3) 0]
params = [X_p_1(1)  X_p_1(3) 0;...
          X_p_1(2)  0  X_p_1(3);...
          X_p_2(1)  X_p_2(3) 0]\(lambda1*[p_1(1) p_1(2) p_2(1)]')

K = [params(1) 0 params(2);...
    0 params(1) params(3);...
    0 0 1]

%% 2.2
% a
p_c = [3264 2448]/2;
fx = (lambda1*p_1(1)-p_c(1)*X_p_1(3))/X_p_1(1)
fy = (lambda1*p_1(2)-p_c(2)*X_p_1(3))/X_p_1(2)

% b
ux = (lambda1*p_2(1)-X_p_2(1)*fx)/X_p_2(3)
uy = (lambda1*p_2(2)-X_p_2(2)*fy)/X_p_2(3)

% c
fx_c = (lambda1*p_1(1)-ux*X_p_1(3))/X_p_1(1)
fy_c = (lambda1*p_1(2)-uy*X_p_1(3))/X_p_1(2)

% d
ux_d = (lambda1*p_2(1)-X_p_2(1)*fx_c)/X_p_2(3)
uy_d = (lambda1*p_2(2)-X_p_2(2)*fy_c)/X_p_2(3)

%% 2.3


%% 2.4
K_new = [fx_c 0 ux_d; 0 fy_c uy_d; 0 0 1]
P_3 = (K_new*lambda2*H)\(lambda1*[p_3 1]');
P_3 = P_3/P_3(3)
P_4 = (K_new*lambda2*H)\(lambda1*[p_4 1]');
P_4 = P_4/P_4(3)
norm(P_3/P_3(3)-P_4/P_4(3))