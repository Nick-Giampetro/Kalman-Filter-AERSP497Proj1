
clc;
clear;
%----------------------------problem 2-------------------------------------
A_c = [0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1; -2 1 0 -1.6 0.8 0; 1 -2 1 0.8 -1.6 0.8; 0 1 -1 0 0.8 -0.8];
B_c = [0 0 0; 0 0 0; 0 0 0; 1 0 0; 0 1 0; 0 0 1];
C_c = [0 1 0 0 0 0 ;0 0 0 1 0 0 ; 0 0 0 0 0 1];

dT = 0.1;

A_d = expm(A_c * dT);
B_d = (expm(A_c * dT) - eye(6))*inv(A_c)*B_c;

sigma_u = [0.2 0;0 0.2];
sigma_v = [0.1 0;0 0.1];
%-------------------discrete time system simulation------------------------
x0 =[1 2 3 0 0 0]';
x_k = x0;
x_k1 = x_k;
u=[0 0 0]';
n=0;
xp = zeros(6,100);
xp(:,1) = x_k;
while(n<100)
x_k = x_k1;
x_k1 = A_d * x_k + B_d*u;
xp(:,n+1) = x_k1;
n=n+1;
end
x2_k = x0;
x2_k1 = x2_k;
n2=0;
xp2 = zeros(6,100);
xp2(:,1) = x2_k;
while(n2<100)
x2_k = x2_k1;
urand = normrnd(0,0.2*eye(3));
uk=[urand(1,1);urand(2,2);urand(3,3)];
x2_k1 = A_d * x2_k + B_d*uk;
xp2(:,n2+1) = x2_k1;
n2=n2+1;
end
t=0:0.1:10-0.1;
h1 = plot(t, xp(1,:), 'r', t, xp(2,:),'g',t, xp(3,:), 'b', t, xp2(1,:), 'c', t, xp2(2,:),'m', t, xp2(3,:),'y');
%h1 = plot(t, xp(1,:), '-', t, xp(2,:),'-.');
xlabel('Time(sec)');
ylabel('masses position');
legend('mass 1','mass 2','mass 3','mass 1 u~(0,0.2I)','mass 2 u~(0,0.2I)','mass 3 u~(0,0.2I)');
%-----------------------continous system simulation------------------------
x0 = [1 2 3 0 0 0]' ;
t0 = 0: 0.1: 10;
zz = my_rk4(@spring_sys,t0,x0);
zz2 = my_rk4(@spring_sys2,t0,x0);
figure
h2 = plot(t0,zz(1,:),'r',t0,zz(2,:),'g',t0,zz(3,:),'b',t0,zz2(1,:),'c',t0,zz2(2,:),'m',t0,zz2(3,:),'y' );

xlabel('Time(sec)');
ylabel('masses position');
legend('mass 1','mass 2','mass 3','mass 1 u~(0,0.2I)','mass 2 u~(0,0.2I)','mass 3 u~(0,0.2I)');
% system model with u=[0 0]'
function dx= spring_sys(t,x)
u=[0 0 0]' ;
dx=[1 2 3 0 0 0]';
A_c = [0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1; -2 1 0 -1.6 0.8 0; 1 -2 1 0.8 -1.6 0.8; 0 1 -1 0 0.8 -0.8];
B_c = [0 0 0; 0 0 0; 0 0 0; 1 0 0; 0 1 0; 0 0 1];
dx = A_c * x +B_c *u;
end
% system model with u ~ N(0, 0.2*eye(2))
function dx= spring_sys2(t,x)
urand = normrnd(0,0.2*eye(3));
u=[urand(1,1);urand(2,2);urand(3,3)];
dx=[1 2 3 0 0 0]';
A_c = [0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1; -2 1 0 -1.6 0.8 0; 1 -2 1 0.8 -1.6 0.8; 0 1 -1 0 0.8 -0.8];
B_c = [0 0 0; 0 0 0; 0 0 0; 1 0 0; 0 1 0; 0 0 1];
dx = A_c * x +B_c *u;
end
function YY = my_rk4(f,t0,x0)
dt = t0(2)-t0(1);
YY = zeros(length(x0),length(t0));
xn = x0;
YY(:,1)= xn;
dx1 = f(dt, xn)*dt;
xh1 = xn + dx1/2;
dx2 = f(dt, xh1)*dt;
xh2 = xn + dx2/2;
dx3 = f(dt, xh2)*dt;
xh3 = xn + dx3;
dx4 = f(dt, xh3)*dt;
xn1 = xn+(dx1+2*dx2+2*dx3+dx4)/6;
flag = 1;
while(flag<length(t0))
xn = xn1;
YY(:,flag+1) = xn;
dx1 = f(dt, xn)*dt;
xh1 = xn + dx1/2;
dx2 = f(dt,xh1)*dt;
xh2 = xn + dx2/2;
dx3 = f(dt,xh2)*dt;
xh3 = xn + dx3;
dx4 = f(dt,xh3)*dt;
xn1 = xn +(dx1+2*dx2+2*dx3+dx4)/6;
flag = flag +1;
end
end
