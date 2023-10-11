
clc;
clear;
%----------------------------problem 2-------------------------------------
A_c = [0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1; -2 1 0 -1.6 0.8 0; 1 -2 1 0.8 -1.6 0.8; 0 1 -1 0 0.8 -0.8];
B_c = [0 0 0; 0 0 0; 0 0 0; 1 0 0; 0 1 0; 0 0 1];
C_c = [1 0 0 0 0 0 ;0 1 0 0 0 0 ; 0 0 1 0 0 0];

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

for n = 0:99
x_k = x_k1;
x_k1 = A_d * x_k + B_d*u;
xp(:,n+1) = x_k1;
end

x2_k = x0;
x2_k1 = x2_k;
n2=0;
xp2 = zeros(6,100);
xp2(:,1) = x2_k;

for n2 = 0:99
urandv = normrnd(0,0.1*eye(3));
ukv=[urandv(1,1);urandv(2,2);urandv(3,3)];
y(:,n2+1) = C_c * x2_k + ukv ;

x2_k = x2_k1;
urandu = normrnd(0,0.2*eye(3));
uku=[urandu(1,1);urandu(2,2);urandu(3,3)];
x2_k1 = A_d*x2_k + B_d*uku;
xp2(:,n2+1) = x2_k1;
end

t=0:dT:10-0.1;


figure
h0 = plot(t, y(1,:),'r',t, y(2,:),'g',t, y(3,:),'b');
xlabel('Time(sec)');
ylabel('masses position according to sensor');
legend('mass 1 u~(0,0.1I)','mass 2 u~(0,0.1I)','mass 3 u~(0,0.1I)');

figure
h1 = plot(t, xp(1,:), '--', t, xp(2,:),'--',t, xp(3,:), '--', t, xp2(1,:), 'r', t, xp2(2,:),'g', t, xp2(3,:),'b');
xlabel('Time(sec)');
ylabel('masses position aacording to dynamics');
legend('mass 1 (ideal)','mass 2 (ideal)','mass 3 (ideal)','mass 1 u~(0,0.2I)','mass 2 u~(0,0.2I)','mass 3 u~(0,0.2I)');

