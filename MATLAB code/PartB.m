clear
clc
% Darian Irani
% 30 Nov 2022
% Part B

%assigning variables
x1=3.72;
x2=2.55;
x3=2.54;
y1=-0.82;
y2=-0.27;
y3=-1.70;
E=215*10^9;
fx=115.4;
fy=158.5;
v=0.33;

%% Question 2

%setting matrices for K
delta=[(x2*y3-x3*y2)+(x3*y1-x1*y3)+(x1*y2-x2*y1)]*0.5;
B=(1/(2*delta))*[y2-y3 0 y3-y1 0 y1-y2 0;0 x3-x2 0 x1-x3 0 x2-x1;x3-x2 y2-y3 x1-x3 y3-y1 x2-x1 y1-y2];
D=E/(1-v-v^2)*[(1-0.33) 0.33 0;0.33 (1-0.33) 0;0 0 (0.5-0.33)];

%obtain stiffness matrix
K=delta*transpose(B)*D*B;

u = [0.765;3.98;0;0;0;0]*(10^-3);
f = K*u

%% Question 3

%find u1x and u1y
k=[1.56 -0.01; -0.01 0.40]*10^11;
f=[115.4;158.5]*10^6;
u=(inv(k))*f;

%% Rough work

u_2=[0.765;3.98]*(10^-3);
u_2

