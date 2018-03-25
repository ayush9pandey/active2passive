s = tf('s');
C0 = 1e-6;
L = 2e-1;
R1a = 1e7;
R1b = 1e-1;
A = [-1/(R1b*C0) 1/C0;-1/L 0];
B = [0 1/C0;1/L 0];
C = eye(2);
D = 0;
G = ss(A,B,C,D,'StateName',{'Vc' 'Il'},'InputName',{'Vref','id'},'OutputName',{'Vc','Il'});

%% Plots

subplot(2,2,1);
step(G(1,1))
subplot(2,2,2);
impulse(G(1,2))
% subplot(2,2,3);
% step(G(2,2))
% subplot(2,2,4);
% impulse(G(2,1))