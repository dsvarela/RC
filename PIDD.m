s = tf('s');
Gp = (-0.07988 *s^4 - 0.003315 *s^3 - 0.8677* s^2 + 0.006493* s - 0.03458)/...
    (s^5 + 0.5979* s^4 + 10.98 *s^3 + 4.709 *s^2 + 0.5421* s + 0.1827);

Gv = (0.2204 *s^4 + 0.02348 *s^3 + 2.394* s^2 + 0.003981 *s + 0.09541)/...
    (s^5 + 0.5979 *s^4 + 10.98 *s^3 + 4.709 *s^2 + 0.5421 *s + 0.1827);

%%
% C= -2.759 s^2 - 0.003532 s - 0.1104 / s^3 + 2 s^2 + s
% beta = 10;
% T = 1/1.97;
% LagC = beta*(T*s+1)/(beta*T*s+1);
% Gc = -.023/s*LagC;
% PIDap = K/Ti * (1+s*Ti)*(1+s*Td)/s/(1+s*Tf);
% Gc = K/Ti * (1+ s*(Ti+Tf) + Ti*(Td+Tf)*s^2)/s/(1+s*Tf);
% LP2 = w^2/(s^2 + 2*w*b*s + w^2)
% N = (s^2 + 2*w*b1*s + w^2)(s^2 + 2*w*b2*s + w^2)
%I;
% 
% w = .21;
% b = .25; % 0.01
% K = -.33;
% Gc = K*(s^2 + 2*w*b*s + w^2)/s/(s^2 + 2*w*.5*s + w^2);

Gc0 = -.18/s;

beta = 10;
T = 1/1.97;
LagC = beta*(T*s+1)/(beta*T*s+1);
Gc1 = -.0443/s/(s+.197);

w = .21;
b = .25; % 0.01
K = -.33;
Gc = K*(s^2 + 2*w*b*s + w^2)/s/(s^2 + 2*w*.5*s + w^2);

L = Gc * Gp;
close all
figure('Position', [0 40 960 960]);
subplot(2,1,1);
bode(Gc, Gc1);
subplot(2,1,2)
bode(Gc*Gp, Gc1*Gp);

figure('Position', [960 40 960 960]);
CL = feedback(L,1);
DL = (Gv/(1+L));

subplot(2,1,1);
step(CL);
hold on
step(feedback(Gc1*Gp,1));
% step(feedback(Gc1*Gp,1));
hold off
subplot(2,1,2)
step(DL);
hold on
%step((Gv/(1+Gc1*Gp)));
% step((Gv/(1+Gc2*Gp)));
hold off
stepinfo(CL).Overshoot,stepinfo(CL).SettlingTime

%% DR
close all
%P + Lag Filter;
%Best for Reference Tracking
w = 1.97;
b1 = .1;
b2 = 10;
N = (s^2 + 2*w*b1*s+w^2)/(s^2 + 2*w*b2*s+w^2);
Gcd = -.2  *(s^2 + 2*w*b1*s+w^2)/(s^2 + 2*w*b2*s+w^2)/s;

L = Gcd*Gp;
CL = Gcd*Gp/(1+Gcd*Gp);
S = 1/(1+Gcd*Gp);   
figure;
bode(Gp, Gcd, L,S)
figure;
step(Gv*S)
[stepinfo(CL).Overshoot, stepinfo(CL).SettlingTime];
