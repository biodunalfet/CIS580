%% initialize
clear all
close all
clc

%% display image
im = imread('IMG_0220.JPG');
figure(1)
clf
picture = imshow(im);
hold on
axis on
axis equal
%set(gca,'YDir','normal')

%% get points
%points = ginput;

%% get joe
%joe = ginput;

%% get door
%door = ginput;

%% get points for z vanishing point lines
%points = ginput;
%{d
%% calculate and draw lines
points1 = load('points.mat');
points2 = load('points2.mat');
points3 = load('points3.mat');

plot([points1.p1(1) points1.p2(1)],[points1.p1(2) points1.p2(2)],'r*-')
plot([points1.p3(1) points1.p4(1)],[points1.p3(2) points1.p4(2)],'r*-')
plot([points2.p1(1) points2.p2(1)],[points2.p1(2) points2.p2(2)],'r*-')
plot([points2.p3(1) points2.p4(1)],[points2.p3(2) points2.p4(2)],'r*-')
plot(points3.points(1:2,1), points3.points(1:2,2),'r*-')
plot(points3.points(3:4,1), points3.points(3:4,2),'r*-')

% get y vanishing point
l12 = cross([points1.p1 1],[points1.p2 1]);
l12 = l12/norm(l12(1:2));

l34 = cross([points1.p3 1],[points1.p4 1]);
l34 = l34/norm(l34(1:2));

vy = cross(l12,l34);
vy = vy/vy(3);
plot(vy(1),vy(2),'r*')

% get x vanishing point
l12 = cross([points2.p1 1],[points2.p2 1]);
l12 = l12/norm(l12(1:2));

l34 = cross([points2.p3 1],[points2.p4 1]);
l34 = l34/norm(l34(1:2));

vx = cross(l12,l34);
vx = vx/vx(3);
plot(vx(1),vx(2),'r*')

% get z vanishing opint
l12 = cross([points3.points(1,:) 1],[points3.points(2,:) 1]);
l12 = l12/norm(l12(1:2));

l34 = cross([points3.points(3,:) 1],[points3.points(4,:) 1]);
l34 = l34/norm(l34(1:2));

vz = cross(l12,l34);
vz = vz/vz(3);

plot(vz(1),vz(2),'r*')

plot([vx(1) vy(1)],[vx(2) vy(2)],'r-')
xlim([vy(1) vz(1)])
ylim([0 size(im,1)])
%}
%{d
%% calculate stuff
% line from joe's feet to bottom of door
door = load('door.mat');
joe = load('joe.mat');
joe_top = joe.joe(2,:);
door_top = door.door(2,:);
joe_bot = joe.joe(1,:);
door_bot = door.door(1,:);

% line between bottom of joe and bottom of door
door_joe_bot_line = cross([joe_bot 1],[door_bot 1]);
%plot([joe_bot(1) door_bot(1)],[joe_bot(2) door_bot(2)],'r*-')

vanishing = cross(vx,vy);

% intersection of bottom of joe and door to horizon
bot_vanish_int = cross(vanishing,door_joe_bot_line);
bot_vanish_int = bot_vanish_int/bot_vanish_int(3);

% line connecting top of joe to vanishing point intersection with bottom
% line
top_vanish_line = cross(bot_vanish_int, [joe_top 1]);
plot([joe_bot(1) bot_vanish_int(1)],[joe_bot(2) bot_vanish_int(2)],'r*-')
plot([bot_vanish_int(1) joe_top(1)],[bot_vanish_int(2) joe_top(2)],'r*-')
plot([door_bot(1) door_top(1)],[door_bot(2) door_top(2)],'r*-')

% line connecting door bottom to door top
door_bot_top_line = cross([door_bot 1],[door_top 1]);

% intersection of door vertical and joe top to horizon
top_vanish_int = cross(top_vanish_line, door_bot_top_line);
top_vanish_int = top_vanish_int/top_vanish_int(3);
plot(top_vanish_int(1),top_vanish_int(2),'r*')
plot(bot_vanish_int(1),bot_vanish_int(2),'r*')

% line connecting joe bottom to joe top
joe_bot_top_line = cross([joe_bot 1],[joe_top 1]);

% line connecting top of door to vanishing point intersection with bottom
% line
top_door_vanish_line = cross(bot_vanish_int, [door_top 1]);

% intersection of joe vertical and door top to horizon
joe_vert_door_top_int = cross(top_door_vanish_line,joe_bot_top_line);
joe_vert_door_top_int = joe_vert_door_top_int/joe_vert_door_top_int(3);

plot(joe_vert_door_top_int(1),joe_vert_door_top_int(2),'r*')

joe_height = 1917.7; % joe's height in meters

A = joe_bot;
B = joe_top;
C = joe_vert_door_top_int(1:2);
D = vz(1:2);

d = norm(C-A)*norm(D-B)/(norm(B-A)*norm(D-C))*joe_height

% intersection of horizon with joe vertical
joe_horiz_int = cross(vanishing,joe_bot_top_line);
joe_horiz_int = joe_horiz_int/joe_horiz_int(3);
plot(joe_horiz_int(1),joe_horiz_int(2),'r*')

A_p = joe_bot;
B_p = joe_horiz_int(1:2);
C_p = joe_top;
D_p = vz(1:2);

h = (norm(B_p-A_p)*norm(D_p-C_p))/(norm(C_p-A_p)*norm(D_p-B_p))*joe_height
%}