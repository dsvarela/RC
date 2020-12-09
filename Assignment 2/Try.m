s=tf('s');
G=1/(0.2*s+1)/(s+1)*[ 1 1; 1+2*s 2];
wB1=0.25; % desired closed-loop bandwidth
wB2=0.25; % desired closed-loop bandwidth
A=1/1000; % desired disturbance attenuation inside bandwidth
M=1.5 ; % desired bound on hinfnorm(S)
Wp=[(s/M+wB1)/(s+wB1*A) 0; 0 (s/M+wB2)/(s+wB2*A)]; % Sensitivity weight
Wu=eye(2); % Control weight
Wt=[]; % Empty weight


[K,N,GAM,INFO]=mixsyn(G,Wp,Wu,Wt);
K = minreal(K);
K = tf(K);
T = G*K/(eye(2) + G*K);
S = inv(eye(2) + G*K);

%Uses the robust control toolbox
systemnames ='G Wp Wu'; % Define systems
inputvar ='[r(2); u(2)]'; % Input generalized plant
input_to_G= '[u]';
input_to_Wu= '[u]';
input_to_Wp= '[r-G]';
outputvar= '[Wp; Wu; r-G]'; % Output generalized plant
cleanupsysic = 'yes';
sysoutname='P';
sysic;


[K2,CL2,GAM2,INFO2] = hinfsyn(P,2,2); % Hinf design
K2 = minreal(K2);
K2 = tf(K2);
T2 = G*K2/(eye(2) + G*K2);
S2 = inv(eye(2) + G*K2);