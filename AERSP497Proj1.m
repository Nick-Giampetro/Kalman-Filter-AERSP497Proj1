
clc;
clear;
close all;


A_c = [0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1; -2 1 0 -1.6 0.8 0; 1 -2 1 0.8 -1.6 0.8; 0 1 -1 0 0.8 -0.8];
B_c = [0 0 0; 0 0 0; 0 0 0; 1 0 0; 0 1 0; 0 0 1];
C_c = [1 0 0 0 0 0 ;0 1 0 0 0 0 ; 0 0 1 0 0 0];

dT = 0.1;
A_d=expm(A_c*dT);
B_d=(A_d-eye(6))*inv(A_c)*B_c ;


x0 =[1 2 3 0 0 0]';
x_k = x0;
x_k1 = x_k;
u=[0 0 0]';
xp = zeros(6,100);
xp(:,1) = x_k;

for k = 0:99
    x_k = x_k1;
    x_k1 = A_d * x_k + B_d*u;
    xp(:,k+1) = x_k1;
end

x2_k = x0;
x2_k1 = x2_k;
xp2 = zeros(6,100);
xp2(:,1) = x2_k;

y = zeros(3,100);

for k = 0:99
    urandv = normrnd(0,0.1*eye(3));
    ukv=[urandv(1,1);urandv(2,2);urandv(3,3)];
    y(:,k+1) = C_c * x2_k + ukv ;
    
    x2_k = x2_k1;
    urandu = normrnd(0,0.2*eye(3));
    uku=[urandu(1,1);urandu(2,2);urandu(3,3)];
    x2_k1 = A_d*x2_k + B_d*uku;
    xp2(:,k+1) = x2_k1;
end



x_hat = zeros(6,100);
u_bar=0 ;

S_u = 0.2 * eye(3) ;
S_v = 0.1 * eye(3) ;
Q_d=B_d*S_u*B_d' ;

P0 = blkdiag(eye(3), 0.1*eye(3)) ;
P_hat(:,:,1) = P0 ;
xhat0 = x0 ;
x_hat(:,1)= xhat0 ;
x_tru(:,1)= x0 ;


for k = 2:100
    urandu = normrnd(u_bar,S_u*eye(3)) ;
    u_true(:,k-1)=[urandu(1,1);urandu(2,2);urandu(3,3)]; % gaussian random noise with mean S_u
    x_tru(:,k)=A_d*x_tru(:,k-1) + B_d*u_true(:,k-1); % true state change
    
    urandv = normrnd(u_bar,S_v*eye(3)) ;
    u_sensor(:,k-1)= [urandv(1,1);urandv(2,2);urandv(3,3)] ;
    y1(:,k-1) = C_c*x_tru(:,k) + sqrtm(S_v)*u_sensor(:,k-1);
    
    R = eye(3)*S_v*eye(3)' ;

    x_hat(:,k) = A_d*x_hat(:,k-1); % the model assumes ubar=0
    P_hat(:,:,k) = A_d*P_hat(:,:,k-1)*A_d' + Q_d;
    
    % measurement update
    K(:,:,k)=P_hat(:,:,k)*C_c'*inv(C_c*P_hat(:,:,k)*C_c' + R);
    x_hat(:,k)=x_hat(:,k) + K(:,:,k)*(y1(k) - C_c*x_hat(:,k));
    P_hat(:,:,k)=(eye(6) - K(:,:,k)*C_c)*P_hat(:,:,k);
    
    % now compute estimate error
    e(:,k)=x_tru(:,k)-x_hat(:,k);
    trP(k)=trace(P_hat(:,:,k));
    Sx(:,k)=sqrt(diag(P_hat(:,:,k)));
end

urandv = normrnd(u_bar,S_v*eye(3)) ;
u_sensor(:,100)= [urandv(1,1);urandv(2,2);urandv(3,3)] ;
y1(:,100) = C_c*x_tru(:,k) + sqrtm(S_v)*u_sensor(:,k-1);

sig=sqrt([squeeze(P_hat(1,1,:))';squeeze(P_hat(2,2,:))']);

t=0:dT:10-dT;

figure
h1 = plot(t, y1(1,:),'r',t, y1(2,:),'g',t, y1(3,:),'b');
xlabel('Time(sec)');
ylabel('masses position according to sensor');
legend('mass 1 u~(0,0.1I)','mass 2 u~(0,0.1I)','mass 3 u~(0,0.1I)');

figure
h0 = plot(t, x_tru(1,:), '--', t, x_tru(2,:),'--',t, x_tru(3,:), '--', t, x_hat(1,:),'r',t, x_hat(2,:),'g',t, x_hat(3,:),'b');
xlabel('Time(sec)');
ylabel('masses position according to dynamics');
legend('mass 1 u~(0,0.2I)','mass 2 u~(0,0.2I)','mass 3 u~(0,0.2I)','mass 1 (KF)','mass 2 (KF)','mass 3 (KF)');


