
clc;
clear;
close all;


A_c = [0 0 0 1 0 0 ; 
       0 0 0 0 1 0 ; 
       0 0 0 0 0 1 ; 
       -2 1 0 -1.6 0.8 0 ; 
       1 -2 1 0.8 -1.6 0.8 ; 
       0 1 -1 0 0.8 -0.8] ;
B_c = [0 0 0 ; 
       0 0 0 ; 
       0 0 0 ; 
       1 0 0 ;  
       0 1 0 ; 
       0 0 1] ;
C = [1 0 0 0 0 0 ;0 1 0 0 0 0 ; 0 0 1 0 0 0];

dT = 0.1;
A_d=expm(A_c*dT);
B_d=(A_d-eye(6))*inv(A_c)*B_c ;

x0 =[1 2 3 0 0 0]';
x_k = x0;
x_k1 = x_k;

x_hat = zeros(6,100);
u_bar = 0 ;

S_u = 0.2 ;
S_v = 0.1 ;

Q = B_d*S_u*B_d' ;
R = 0.1*eye(3) ;

P0 = blkdiag(eye(3),0.1*eye(3)) ;
% P0 = zeros(6,6)' ;
P(:,:,1) = P0 ;
xhat0 = x0 ;
x_hat(:,1)= xhat0 ;
x_tru(:,1)= x0 ;

for rn = 1:50
    for k = 2:100
        urandu = normrnd(u_bar,S_u*eye(3)) ;
        u_true(:,k-1)=[urandu(1,1);urandu(2,2);urandu(3,3)]; % gaussian random noise with mean S_u
        x_tru(:,k)=A_d*x_tru(:,k-1) + B_d*u_true(:,k-1); % true state change
      
        urandv = normrnd(u_bar,S_v*eye(3)) ;
        u_sensor(:,k-1)= [urandv(1,1);urandv(2,2);urandv(3,3)] ;
        y(:,k-1) = C*x_tru(:,k) + sqrtm(S_v)*u_sensor(:,k-1);
    
        x_hat(:,k) = A_d*x_hat(:,k-1) ;         % the model assumes ubar=0
        P(:,:,k) = A_d*P(:,:,k-1)*A_d' + Q;
    
        % measurement update
        K(:,:,k)=P(:,:,k)*C'*inv(C*P(:,:,k)*C' + R);
        x_hat(:,k)=x_hat(:,k) + K(:,:,k)*(y(:,k-1) - C*x_hat(:,k));
        P(:,:,k)=(eye(6) - K(:,:,k)*C)*P(:,:,k);
        
        % now compute estimate error
        e(:,k)=x_tru(:,k)-x_hat(:,k);
        e_norm(rn,k) = norm(e(:,k));
        P_trace(rn,k)=sqrt(trace(P(:,:,k)));
    end
    
    urandv = normrnd(u_bar,S_v*eye(3)) ;
    u_sensor(:,100)= [urandv(1,1);urandv(2,2);urandv(3,3)] ;
    y(:,100) = C*x_tru(:,k) + sqrtm(S_v)*u_sensor(:,k-1) ;
end

t=0:dT:10-dT;

figure
h0 = plot(t,mean(e_norm(:,:)),t,mean(P_trace(:,:)));
xlabel('Time(sec)');
ylabel('Average of errors at each step in time');
legend('Avg error','P Trace');