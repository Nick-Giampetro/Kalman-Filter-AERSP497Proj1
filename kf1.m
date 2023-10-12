% KF 1
close all
clear all

% time parameters
dT=0.1;
t0=0;
tf=20;

%define filter parameters
xhat0=[0;0;0];
P0=eye(3);
Q_extra=0*eye(3);   % fictitious process noise
R=0.1; % assumed measurement noise

% define true initial conditions
x0=[2;0;1];

% define model parameters
m=1;
b=0.4;
k=1;

% define noise parameters
S_u=0.2;
S_v=0.1;

% measurement of position only
C=[1 0 0];

% continuous time dynamics
A_c=[0 1 0;-(k)/m -b/m 1/m;0 0 0];     % assumed dynamics

B_c=[0;1/m;0];

sys_c=ss(A_c,B_c,C,0);
[A,B,C,D]=ssdata(c2d(sys_c,dT,'zoh'));

% Q=B'Su*B
Q=B*S_u*B';

time=t0:dT:tf;
MAXK=length(time);

x_tru(:,1)=x0;      % initialize truth
x_hat(:,1)=xhat0;   % initialize estimate
P_hat(:,:,1)=P0;
e(:,1)=x_tru(:,1)-x_hat(:,1);
trP(1)=trace(P_hat(:,:,1));
Sx(:,1)=sqrt(diag(P_hat(:,:,1)));
for k=2:MAXK
    % simulate model dynamics by one time step
    u_true(k-1)=sqrtm(S_u)*randn(1,1); % gaussian random noise with mean u_bar
    x_tru(:,k)=A*x_tru(:,k-1) + B*u_true(k-1); % true state change
    y(k)=C*x_tru(:,k) + sqrtm(S_v)*randn(1,1); % get a measurement (noisy)
    
    % time update the estimate
    x_hat(:,k)=A*x_hat(:,k-1); 
    P_hat(:,:,k)=A*P_hat(:,:,k-1)*A' + Q;
    
    % measurement update
    K(:,k)=P_hat(:,:,k)*C'*inv(C*P_hat(:,:,k)*C' + R);
    x_hat(:,k)=x_hat(:,k) + K(:,k)*(y(k) - C*x_hat(:,k));
    P_hat(:,:,k)=(eye(3) - K(:,k)*C)*P_hat(:,:,k);
    
    % now compute estimate error
    e(:,k)=x_tru(:,k)-x_hat(:,k);
    trP(k)=trace(P_hat(:,:,k));
    Sx(:,k)=sqrt(diag(P_hat(:,:,k)));
end

h0 = plot(time, x_tru(1,:), '--', time, x_hat(1,:));
xlabel('Time(sec)');
ylabel('masses position according to dynamics');
legend('mass 1 u~(0,0.2I)','mass 1 (KF)');

figure
subplot(3,1,1)
    plot(time,e(1,:),'b',time,Sx(1,:),'b:',time,-Sx(1,:),'b:')
    xlim([t0 tf])
    ylabel('x_1')
subplot(3,1,2)
    plot(time,e(2,:),'b',time,Sx(2,:),'b:',time,-Sx(2,:),'b:')
    xlim([t0 tf])
    ylabel('x_2')
subplot(3,1,3)
    plot(time,e(3,:),'b',time,Sx(3,:),'b:',time,-Sx(3,:),'b:')
    xlim([t0 tf])
    ylabel('x_3')
    