%% Load Data
close all
clc
s = tf('s'); % For transfer functions

load('Assignment_Data_SC42145-mat')

[num,den] = ss2tf(A,B,C,D,1);
Gp = tf(num(1,:),den); % TF - Blade Picth to Rotational Velocity

[num,den] = ss2tf(A,B,C,D,3);
Gv = tf(num(1,:),den); % TF - Wind Disturbance to Rotational Velocity

clear num den
%% Controllers
Gc1 = -.26/s; % Integrator
L1 = Gp*Gc1; % Open Loop
CLr1 = feedback(L1,1); % Closed Loop for Reference Tracking
CLd1 = Gv/(1+L1);  % Closed Loop for Disturbance Rejection
 
Gc2 = -.225/s/(5*s+1); % Integrator + Filter
L2 = Gp*Gc2; % Open Loop
CLr2 = feedback(L2,1); % Closed Loop for Reference Tracking
CLd2 = Gv/(1+L2);  % Closed Loop for Disturbance Rejection

% Integrator + Notch
w = sqrt(.04101);
b2 = .295;
K = -.3128;
Gc3 = K*(s^2 + 0.02109*s + .04101)/s/(s^2 + 2*w*b2*s + .04101);  
L3 = Gp*Gc3; % Open Loop
CLr3 = feedback(L3,1); % Closed Loop for Reference Tracking
CLd3 = Gv/(1+L3);  % Closed Loop for Disturbance Rejection

%% Plot Stuff
figure('Position', [0 40 960*2 960]);
subplot(2,2,1);
bode(L1, L2, L3, Gp);
grid on;
legend('$L1$','$L2$','$L3$','$Gp$','Interpreter','Latex',...
    'Fontsize', 16, 'Location', 'southeast');
title('Bode Diagram of Open Loop TFs and Uncontrolled Plant', ...
'Fontsize', 18,'FontWeight','bold');
subplot(2,2,3);
grid on;
bode(1/(1+L1), 1/(1+L2), 1/(1+L3));
legend('$S1$','$S2$','$S3$','Interpreter','Latex',...
    'Fontsize', 16, 'Location', 'southeast');
title('Bode Diagram of Sensitivity Functions', ...
'Fontsize', 18,'FontWeight','bold');
subplot(2,2,2);
grid on;
step(CLr1,CLr2,CLr3);
legend('$T1$','$T2$','$T3$','Interpreter','Latex',...
    'Fontsize', 16, 'Location', 'southeast');
title('Step Response - Reference Tracking', ...
'Fontsize', 18,'FontWeight','bold');

subplot(2,2,4);
grid on;
step(CLd1,CLd2,CLd3);
legend('$GvS1$','$GvS2$','$GvS3$','Interpreter','Latex',...
    'Fontsize', 16, 'Location', 'northeast');
title('Step Response - Disturbance Rejection', ...
'Fontsize', 18,'FontWeight','bold');
hold off