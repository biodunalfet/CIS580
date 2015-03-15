function H = dh_to_homogeneous(alpha,a,d,theta)

T_d = [1 0 0 0; 0 1 0 0; 0 0 1 d; 0 0 0 1];
R_theta = [cos(theta) -sin(theta) 0 0; sin(theta) cos(theta) 0 0; 0 0 1 0; 0 0 0 1];
T_a = [1 0 0 a; 0 1 0 0; 0 0 1 0; 0 0 0 1];
R_alpha = [1 0 0 0; 0 cos(alpha) -sin(alpha) 0; 0 sin(alpha) cos(alpha) 0; 0 0 0 1];
H = T_d*R_theta*T_a*R_alpha;