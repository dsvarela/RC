% Uses the Mu-toolbox
G=nd2sys(1,conv( [10 1],conv([0.05 1],[0.05 1])),200); %Plant is G.
M=1.5; wb=10; A=1.e-4; Wp = nd2sys([1/M wb], [1 wb*A]); Wu = 1; % Weights.
%
% Generalized plant P is found with function sysic:
%
s = tf('s');
wp2 = (s/M^(1/2)+ wb)^2/(s + A^(1/2)*wb)^2;
wp2 = ss(wp2);
wp2 = minreal(wp2);
Wp = pck(wp2.A,wp2.B,wp2.C,wp2.D);
systemnames = 'G Wp Wu';
inputvar = '[r(1); u(1)]';
outputvar = '[Wp; Wu; r-G]';
input_to_G = '[u]';
input_to_Wp = '[r-G]';
input_to_Wu = '[u]';
sysoutname = 'P';
cleanupsysic = 'yes';
sysic;
%
% Find H-infinity optimal controller:
%
nmeas=1; nu=1; gmn=0.5; gmx=20; tol=0.001;

% [K,CL,GAM] = hinfsyn(P,NMEAS,NCON,GAMTRY);
[K,CL,GAM] = hinfsyn(P,nmeas,nu,gmn,gmx,tol);
[a,b,c,d] = unpck(K);
[num,den]=ss2tf(a,b,c,d);
K = tf(num,den);

G = 200/(10*s+1)*1/(0.05*s+1)^2;
Wp = (s/M^(1/2)+ wb)^2/(s + A^(1/2)*wb)^2;
Wu = 1;

p11 = [-Wp ; 0];
p12 =  [-Wp*G ; -Wu];
p21 = -1;
p22 = -G;
P = [p11 , p12; p21, p22];
[K2,CL2,GAM2] = hinfsyn(P,nmeas,nu);

K2 = minreal(K2)
[num,den]=ss2tf(K2.A,K2.B,K2.C,K2.D);
K2 = tf(num,den);