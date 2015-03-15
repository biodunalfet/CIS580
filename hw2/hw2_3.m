clear all
close all
clc
addpath('fp')
im1 = imread('image_a_front.jpg');
im2 = imread('image_b_pan.jpg');
im3 = imread('image_c_tilt.jpg');

K = [2810.8515 0 1643.16492;...
    0 2803.50997 1240.43791;...
    0 0 1];

M_a = [0.9514 0.3148 -102.4771;...
    0.0159 0.0098 910.3956;...
    0.3076 -1.0897 2677.0886];

M_b = [0.9834 0.1730 277.4586;...
    0.0147 0.0102 920.2367;...
    0.1806 -1.0812 2655.7909];

M_c = [0.9276 0.3454 -607.849;...
    -0.1191 0.3374 260.2264;...
    0.354 -0.8441 1788.7189];

figure(1)
clf
%subplot(3,1,1);
im_plot_1 = imshow(im1);
axis on
hold on
%points = ginput
%{d
vx1 = load('vx1.mat');
vy1 = load('vy1.mat');
plot(vx1.points(1:2,1),vx1.points(1:2,2),'y-*')
plot(vx1.points(3:4,1),vx1.points(3:4,2),'y-*')
plot(vy1.points(1:2,1),vy1.points(1:2,2),'g-*')
plot(vy1.points(3:4,1),vy1.points(3:4,2),'g-*')

% x vanishing point
l12 = cross([vx1.points(1,:) 1],[vx1.points(2,:) 1]);
l34 = cross([vx1.points(3,:) 1],[vx1.points(4,:) 1]);
vx_1 = cross(l12,l34);
vx_1 = vx_1/vx_1(3);

% y vanishing point
l12 = cross([vy1.points(1,:) 1],[vy1.points(2,:) 1]);
l34 = cross([vy1.points(3,:) 1],[vy1.points(4,:) 1]);
vy_1 = cross(l12,l34);
vy_1 = vy_1/vy_1(3);
plot(vy_1(1),vy_1(2),'r*')

% plot manually calculated horizon line
plot([vx_1(1) vy_1(1)],[vx_1(2) vy_1(2)],'r-*')

% plot actual horizon
lp_1 = (K*M_a)'\[0 0 1]';
p1_1 = [0 -lp_1(3)/lp_1(2)];
p1_2 = [size(im1,2) (-lp_1(3)-lp_1(1)*size(im1,2))/lp_1(2)];
plot([p1_1(1) p1_2(1)],[p1_1(2) p1_2(2)],'b*-')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
clf
%subplot(3,1,2);
im_plot_2 = imshow(im2);
axis on
hold on
%points = ginput
%{d
vx2 = load('vx2.mat');
vy2 = load('vy2.mat');
plot(vx2.points(1:2,1),vx2.points(1:2,2),'y-*')
plot(vx2.points(3:4,1),vx2.points(3:4,2),'y-*')
plot(vy2.points(1:2,1),vy2.points(1:2,2),'g-*')
plot(vy2.points(3:4,1),vy2.points(3:4,2),'g-*')

% x vanishing point
l12 = cross([vx2.points(1,:) 1],[vx2.points(2,:) 1]);
l34 = cross([vx2.points(3,:) 1],[vx2.points(4,:) 1]);
vx_2 = cross(l12,l34);
vx_2 = vx_2/vx_2(3);

% y vanishing point
l12 = cross([vy2.points(1,:) 1],[vy2.points(2,:) 1]);
l34 = cross([vy2.points(3,:) 1],[vy2.points(4,:) 1]);
vy_2 = cross(l12,l34);
vy_2 = vy_2/vy_2(3);
plot(vy_2(1),vy_2(2),'r*')

% plot manually calculated horizon line
plot([vx_2(1) vy_2(1)],[vx_2(2) vy_2(2)],'r-*')

% plot actual horizon
lp_2 = (K*M_b)'\[0 0 1]';
p2_1 = [0 -lp_2(3)/lp_2(2)];
p2_2 = [size(im2,2) (-lp_2(3)-lp_2(1)*size(im2,2))/lp_2(2)];
plot([p2_1(1) p2_2(1)],[p2_1(2) p2_2(2)],'b*-')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3)
clf
%subplot(3,1,3);
im_plot_3 = imshow(im3);
axis on
hold on
%points = ginput;
vx3 = load('vx3.mat');
vy3 = load('vy3.mat');
plot(vx3.points(1:2,1),vx3.points(1:2,2),'y-*')
plot(vx3.points(3:4,1),vx3.points(3:4,2),'y-*')
plot(vy3.points(1:2,1),vy3.points(1:2,2),'g-*')
plot(vy3.points(3:4,1),vy3.points(3:4,2),'g-*')

% x vanishing point
l12 = cross([vx3.points(1,:) 1],[vx3.points(2,:) 1]);
l34 = cross([vx3.points(3,:) 1],[vx3.points(4,:) 1]);
vx_3 = cross(l12,l34);
vx_3 = vx_3/vx_3(3);

% y vanishing point
l12 = cross([vy3.points(1,:) 1],[vy3.points(2,:) 1]);
l34 = cross([vy3.points(3,:) 1],[vy3.points(4,:) 1]);
vy_3 = cross(l12,l34);
vy_3 = vy_3/vy_3(3);
plot(vy_3(1),vy_3(2),'r*')

% plot manually calculated horizon line
plot([vx_3(1) vy_3(1)],[vx_3(2) vy_3(2)],'r-*')

% plot actual horizon
lp_3 = (K*M_c)'\[0 0 1]';
p3_1 = [0 -lp_3(3)/lp_3(2)];
p3_2 = [size(im3,2) (-lp_3(3)-lp_3(1)*size(im3,2))/lp_3(2)];
plot([p3_1(1) p3_2(1)],[p3_1(2) p3_2(2)],'b*-')
%}