clear all
H = [cos(pi/4) -sin(pi/4) 0 0;...
    sin(pi/4) cos(pi/4) 0 0;...
    0 0 1 1;...
    0 0 0 1];

V = eye(4);

V_p = H*V;

q = 1;
elbowOffsetX = 5;
lowerArmLength = 5;
DH = [-pi/2, 0, 0, q;...
    pi/2, 0, 0, q+pi/2;...
    pi/2, 0, upperArmLength, q+pi/2;...
    pi/2, elbowOffsetX, 0, q;...
    -pi/2, -elbowOffsetX, lowerArmLength, q-pi/2;...
    -pi/2, 0, 0, q;...
    pi/2, 0, 0, q;...
    -pi/2, 0, 0, -pi/2];


figure(1)
clf
plot3([V(1,1) 0]+V(1,4),[V(2,1) 0]+V(2,4),[V(3,1) 0]+V(3,4),'r-*')
hold on
grid on
axis equal
plot3([V(1,2) 0]+V(1,4),[V(2,2) 0]+V(2,4),[V(3,2) 0]+V(3,4),'g-*')
plot3([V(1,3) 0]+V(1,4),[V(2,3) 0]+V(2,4),[V(3,3) 0]+V(3,4),'b-*')


plot3([V_p(1,1) 0]+V_p(1,4),[V_p(2,1) 0]+V_p(2,4),[V_p(3,1) 0]+V_p(3,4),'r-*')
plot3([V_p(1,2) 0]+V_p(1,4),[V_p(2,2) 0]+V_p(2,4),[V_p(3,2) 0]+V_p(3,4),'g-*')
plot3([V_p(1,3) 0]+V_p(1,4),[V_p(2,3) 0]+V_p(2,4),[V_p(3,3) 0]+V_p(3,4),'b-*')