
clc;
clear;
close all;


A_c = [0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1; -2 1 0 -1.6 0.8 0; 1 -2 1 0.8 -1.6 0.8; 0 1 -1 0 0.8 -0.8];
B_c = [0 0 0; 0 0 0; 0 0 0; 1 0 0; 0 1 0; 0 0 1];
C = [1 0 0 0 0 0 ;0 1 0 0 0 0 ; 0 0 1 0 0 0];
%C = [0 1 0 0 0 0 ;0 0 0 1 0 0 ; 0 0 0 0 0 1];

dT = 0.1;
A_d=expm(A_c*dT);
B_d=(A_d-eye(6))*inv(A_c)*B_c ;

x0 =[1 2 3 0 0 0]';
x_k = x0;
x_k1 = x_k;

x_hat = zeros(6,100);
u_bar = 0 ;

S_u = 0.2 * eye(3) ;
S_v = 0.1 * eye(3) ;
Q = B_d*S_u*B_d' ;
R = eye(3)*S_v*eye(3)' ;

P0 = blkdiag(eye(3),0.1*eye(3)) ;
% P0 = zeros(6,6)' ;
P_hat(:,:,1) = P0 ;
xhat0 = x0 ;
x_hat(:,1)= xhat0 ;
x_tru(:,1)= x0 ;


for k = 2:100
    urandu = normrnd(u_bar,S_u*eye(3)) ;
    u_true(:,k-1)=[urandu(1,1);urandu(2,2);urandu(3,3)]; % gaussian random noise with mean S_u
    x_tru(:,k)=A_d*x_tru(:,k-1) + B_d*u_true(:,k-1); % true state change
    
    %C = [1 0 0 0 0 0 ;0 1 0 0 0 0 ; 0 0 1 0 0 0];
  
    urandv = normrnd(u_bar,S_v*eye(3)) ;
    u_sensor(:,k-1)= [urandv(1,1);urandv(2,2);urandv(3,3)] ;
    y(:,k-1) = C*x_tru(:,k) + sqrtm(S_v)*u_sensor(:,k-1);

    x_hat(:,k) = A_d*x_hat(:,k-1); % the model assumes ubar=0
    P_hat(:,:,k) = A_d*P_hat(:,:,k-1)*A_d' + Q;
    
    %C = [0 1 1 1 1 1 ; 1 0 1 1 1 1 ; 1 1 0 1 1 1];

    % measurement update
    P_hat(:,:,k) = inv(inv(P_hat(:,:,k)) + C'*inv(R)*C) ;
    K(:,:,k) = P_hat(:,:,k)*C'*inv(R) ;
    x_hat(:,k)=x_hat(:,k) + K(:,:,k)*(y(k) - C*x_hat(:,k)) ;
    
    % now compute estimate error
    e(:,k)=x_tru(:,k)-x_hat(:,k);
    trP(k)=trace(P_hat(:,:,k));
    Sx(:,k)=sqrt(diag(P_hat(:,:,k)));
end

urandv = normrnd(u_bar,S_v*eye(3)) ;
u_sensor(:,100)= [urandv(1,1);urandv(2,2);urandv(3,3)] ;
y(:,100) = C*x_tru(:,k) + sqrtm(S_v)*u_sensor(:,k-1);

sig=sqrt([squeeze(P_hat(1,1,:))';squeeze(P_hat(2,2,:))']);

t=0:dT:10-dT;

figure
h0 = plot(t, y(1,:),'r',t, y(2,:),'g',t, y(3,:),'b');
xlabel('Time(sec)');
ylabel('Mass positions according to sensor');
legend('mass 1 u~(0,0.1I)','mass 2 u~(0,0.1I)','mass 3 u~(0,0.1I)');

figure
h1 = plot(t, x_tru(1,:), '--', t, x_tru(2,:),'--',t, x_tru(3,:), '--', t, x_hat(1,:),'r',t, x_hat(2,:),'g',t, x_hat(3,:),'b');
xlabel('Time(sec)');
ylabel('Mass positions according to dynamics');
legend('mass 1 u~(0,0.2I)','mass 2 u~(0,0.2I)','mass 3 u~(0,0.2I)','mass 1 (KF)','mass 2 (KF)','mass 3 (KF)');

figure
h2 = plot(t,e(1,:),'r',t,e(2,:),'g',t,e(3,:),'b') ;
xlabel('Time(sec)');
ylabel('Error');
legend('mass 1','mass 2','mass 3');

