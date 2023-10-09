% KF 2
close all
clear all

% time parameters
dT=0.1;
t0=0;
tf=20;

%define filter parameters
xhat0=[0;0];
P0=eye(2);
Q_extra=0*[1 0;0 1];   % fictitious process noise
R=0.1; % assumed measurement noise

% define true initial conditions
x0=[2;0];

% define model parameters
m=1;
b=0.4;
k=1;

% define noise parameters
u_bar=0;
S_u=0.2;
S_v=0.1;

% measurement of position only
C=[1 0];

dk=0; % potential for unmodelled plant error...
db=0; % -0.3;

% continuous time dynamics
A_ct=[0 1;-(k+dk)/m -(b+db)/m];	% true dynamics
A_cm=[0 1;-(k)/m -b/m];     % assumed dynamics

B_c=[0;1/m];

% convert to discrete time using ZOH
A_true=expm(A_ct*dT);
A_model=expm(A_cm*dT);
B_true=(A_true-eye(2))*inv(A_ct)*B_c;
B_model=(A_model-eye(2))*inv(A_cm)*B_c;

% Q=B'Su*B
Q_model=B_model*S_u*B_model';

time=t0:dT:tf;
MAXK=length(time);

x_tru(:,1)=x0;      % initialize truth
x_hat(:,1)=xhat0;   % initialize estimate
P_hat(:,:,1)=P0;
Sx(:,1)=sqrt(diag(P_hat(:,:,1)));

for k=2:MAXK
    % simulate model dynamics by one time step
    u_true(k-1)=u_bar + sqrtm(S_u)*randn(1,1); % gaussian random noise with mean u_bar
    x_tru(:,k)=A_true*x_tru(:,k-1) + B_true*u_true(k-1); % true state change
    y(k)=C*x_tru(:,k) + sqrtm(S_v)*randn(1,1); % get a measurement (noisy)
    
    % time update the estimate
    x_hat(:,k)=A_model*x_hat(:,k-1); % the model assumes ubar=0
    P_hat(:,:,k)=A_model*P_hat(:,:,k-1)*A_model' + Q_model + Q_extra;
    
    % measurement update
    Sr(:,:,k)=C*P_hat(:,:,k)*C' + R; % innovations covariance
    K(:,k)=P_hat(:,:,k)*C'*inv(Sr(:,:,k));
    ry(:,k)=y(k)-C*x_hat(:,k); % innovations
    x_hat(:,k)=x_hat(:,k) + K(:,k)*ry(:,k);
    P_hat(:,:,k)=(eye(2) - K(:,k)*C)*P_hat(:,:,k);
    eta(k)=ry(:,k)'*inv(Sr(:,:,k))*ry(:,k);
    
    Sx(:,k)=sqrt(diag(P_hat(:,:,k)));
end
e=x_tru-x_hat; % compute estimate error
sig=sqrt([squeeze(P_hat(1,1,:))';squeeze(P_hat(2,2,:))']); % std dev of estimate error

Sr=squeeze(Sr); % remove singleton dimensions
% plot results

figure
subplot(2,2,1)
    plot(time,e(1,:),'b',time,Sx(1,:),'b:',time,-Sx(1,:),'b:')
    xlim([t0 tf])
    ylabel('x_1')
    
subplot(2,2,2)
    plot(time,e(2,:),'b',time,Sx(2,:),'b:',time,-Sx(2,:),'b:')
    xlim([t0 tf])
    ylabel('x_2')

subplot(2,2,3)
    plot(1:MAXK,ry,'b',1:MAXK,sqrt(Sr),'b:',1:MAXK,-sqrt(Sr),'b:');
    xlim([0 200])
    ylim([-2 2])
    ylabel('innovations')
    text(10,1,'k=knom, b=bnom')
subplot(2,2,4)
    plot(1:MAXK,eta)
    xlim([0 200])
    ylim([0 50])
    ylabel('eta')
    text(10,40,sprintf('mean(eta): %g',mean(eta)))
    