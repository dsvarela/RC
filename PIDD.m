s = tf('s');
Gp = (-0.07988 *s^4 - 0.003315 *s^3 - 0.8677* s^2 + 0.006493* s - 0.03458)/...
    (s^5 + 0.5979* s^4 + 10.98 *s^3 + 4.709 *s^2 + 0.5421* s + 0.1827);

Gv = (0.2204 *s^4 + 0.02348 *s^3 + 2.394* s^2 + 0.003981 *s + 0.09541)/...
    (s^5 + 0.5979 *s^4 + 10.98 *s^3 + 4.709 *s^2 + 0.5421 *s + 0.1827);

Gc1 = -.26/s;

Gc2 = -.225/s/(5*s+1);

w = sqrt(.04101);
b2 = .295;
K = -.3128;
Gc3 = K*(s^2 + 0.02109*s + .04101)/s/(s^2 + 2*w*b2*s + .04101);


close all
figure('Position', [0 40 960 960]);
subplot(2,1,1);
bode(Gc1, Gc2, Gc3);
subplot(2,1,2)
bode(Gc1*Gp, Gc2*Gp, Gc3*Gp);

figure('Position', [960 40 960 960]);

subplot(2,1,1);
step(feedback(Gc1*Gp,1));
hold on
step(feedback(Gc2*Gp,1));
step(feedback(Gc3*Gp,1));
subplot(2,1,2);
step(Gv / (1+ Gc1*Gp));
hold on
step(Gv / (1+ Gc2*Gp));
step(Gv / (1+ Gc3*Gp));
hold off

%%
Gc  = Gc;
close all

figure('Position', [0 40 960*2 960]);
step(feedback(-16.5*Gp,1));
set(gca, 'Fontsize', 16);
grid on;
xlabel('Time','Fontsize', 20,'FontWeight','bold' );
ylabel('Amplitude','Fontsize', 20,'FontWeight','bold')
legend('$Gp$','Interpreter','Latex', 'Fontsize', 24, 'Location', 'southeast');
title('Step Response for the Closed Loop Plant with K=-16.5', ...
'Fontsize', 24,'FontWeight','bold','Interpreter','Latex')

%%
figure('Position', [0 40 960*2 960]);
bode(-16.5*Gp);
grid on;
legend('$Gc1$','$Gc2$','$Gc3$','Interpreter','Latex', 'Fontsize', 24, 'Location', 'southeast')
title('Bode Diagram of $\mathbf{L(j\omega)}$ (Comparison)', ...
'Fontsize', 24,'FontWeight','bold','Interpreter','Latex')