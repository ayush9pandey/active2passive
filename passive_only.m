s = tf('s');
tau0 = 1;
C0 = 1e-5;
R = tau0/C0;
R1a = 1e9;
R1b = 1e4;
a = (1/tau0) + (1/(R1a*C0));
tfVoutVref = (1/tau0) / (s + a);
tfVoutId = (R/tau0) / (s + a);
[Y1,T] = step(tfVoutVref);
Y2 = impulse(tfVoutId,T);

a = (1/tau0) + (1/(R1b*C0));
tfVoutVrefb = (1/tau0) / (s + a);
tfVoutIdb = (R/tau0) / (s + a);
Y1b = step(tfVoutVrefb,T);
Y2b = impulse(tfVoutIdb,T);
%% Plots 
figure
subplot(2,2,1);
hold on
plot(T,Y1,'LineWidth',3)
plot(T,zeros(1,size(Y1b,1)) + 1,'LineWidth',3)
hold off
title('Step Reference Tracking (High R)')
xlabel('Time','FontSize',48)
ylabel('Capacitor voltage','FontSize',48)
t1 = annotation('textbox',[.2 .8 0.1 0.01],'String','Good step tracking','FitBoxToText','on');
t1.FontSize = 48;
xt = get(gca, 'XTick');
set(gca, 'FontSize', 48)
grid on
subplot(2,2,3); 
plot(T,Y2,'g','LineWidth',3)
title('Impulse Disturbance (High R)')
xlabel('Time','FontSize',48)
ylabel('Capacitor voltage','FontSize',48)
t2 = annotation('textbox',[.2 .2 0.1 0.01],'String','Poor disturbance rejection','FitBoxToText','on');
t2.FontSize = 48;
xt = get(gca, 'XTick');
set(gca, 'FontSize', 44)
grid on
subplot(2,2,2); 
hold on
plot(T,Y1b,'LineWidth',3)
plot(T,zeros(1,size(Y1b,1)) + 1,'LineWidth',3)
hold off
title('Step Reference Tracking (Low R)')
xlabel('Time','FontSize',48)
ylabel('Capacitor voltage','FontSize',48)
t3 = annotation('textbox',[.6 .8 0.1 0.01],'String','Poor step tracking','FitBoxToText','on');
t3.FontSize = 48;
xt = get(gca, 'XTick');
set(gca, 'FontSize', 48)
grid on
subplot(2,2,4); 
plot(T,Y2b,'g','LineWidth',3)
title('Impulse Disturbance (Low R)')
xlabel('Time','FontSize',48)
ylabel('Capacitor voltage','FontSize',48)
t4 = annotation('textbox',[.6 .2 0.1 0.01],'String','Good disturbance rejection','FitBoxToText','on');
t4.FontSize = 48;
xt = get(gca, 'XTick');
set(gca, 'FontSize', 44)
grid on