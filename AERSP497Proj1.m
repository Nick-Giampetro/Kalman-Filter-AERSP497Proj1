
clc;
clear;


A_c = [0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1; -2 1 0 -1.6 0.8 0; 1 -2 1 0.8 -1.6 0.8; 0 1 -1 0 0.8 -0.8];
B_c = [0 0 0; 0 0 0; 0 0 0; 1 0 0; 0 1 0; 0 0 1];
C_c = [1 0 0 0 0 0 ;0 1 0 0 0 0 ; 0 0 1 0 0 0];

dT = 0.1;
A_d = expm(A_c * dT);
B_d = (expm(A_c * dT) - eye(6))*inv(A_c)*B_c;

sigma_u = [0.2 0;0 0.2];
sigma_v = [0.1 0;0 0.1];

<<<<<<< HEAD

=======
P0 = blkdiag(eye(3), 0.1*eye(3));
%-------------------discrete time system simulation------------------------
>>>>>>> 4527af9d3ddb1189d43791d63c7e29481c0de2c4
x0 =[1 2 3 0 0 0]';
x_k = x0;
x_k1 = x_k;
u=[0 0 0]';
xp = zeros(6,100);
xp(:,1) = x_k;

for n = 0:99
    x_k = x_k1;
    x_k1 = A_d * x_k + B_d*u;
    xp(:,n+1) = x_k1;
end

x2_k = x0;
x2_k1 = x2_k;
xp2 = zeros(6,100);
xp2(:,1) = x2_k;

y = zeros(3,100);

for n = 0:99
    urandv = normrnd(0,0.1*eye(3));
    ukv=[urandv(1,1);urandv(2,2);urandv(3,3)];
    y(:,n+1) = C_c * x2_k + ukv ;
    
    x2_k = x2_k1;
    urandu = normrnd(0,0.2*eye(3));
    uku=[urandu(1,1);urandu(2,2);urandu(3,3)];
    x2_k1 = A_d*x2_k + B_d*uku;
    xp2(:,n+1) = x2_k1;
end

P0 = blkdiag(eye(3), 0.1*eye(3)) ;

for n = 1:100
    x_hat(:,n)=A_model*x_hat(:,n-1); % the model assumes ubar=0
    P_hat(:,:,n)=A_model*P_hat(:,:,n-1)*A_model' + Q_model + Q_extra;
    
    % measurement update
    Sr(:,:,n)=C*P_hat(:,:,n)*C' + R; % innovations covariance
    K(:,n)=P_hat(:,:,n)*C'*inv(Sr(:,:,n));
    ry(:,n)=y(n)-C*x_hat(:,n); % innovations
    x_hat(:,n)=x_hat(:,n) + K(:,n)*ry(:,n);
    P_hat(:,:,n)=(eye(2) - K(:,n)*C)*P_hat(:,:,n);
    eta(n)=ry(:,n)'*inv(Sr(:,:,n))*ry(:,n);
    
    Sx(:,n)=sqrt(diag(P_hat(:,:,n)));
end
e=x_tru-x_hat; % compute estimate error
sig=sqrt([squeeze(P_hat(1,1,:))';squeeze(P_hat(2,2,:))']);

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

