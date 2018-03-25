clear all
Tspan = 0:0.01:2.5;
tauv = 0.4:0.6:1;
I = zeros(length(Tspan),length(tauv));
EI = zeros(length(Tspan),length(tauv));
for j = 1:length(tauv)
    tau = tauv(j);
    C = 1e-6;
    R = tau / C;
    beta = 1/C;
    alpha = 1/R;
    v0 = 5;
    k = 1.38064852e-23;
    T = 300;
    for i = 1:length(Tspan)
        t = Tspan(i);
        I(i,j) = (alpha^2 * k * T * beta) + (alpha^2 * k * T * beta * exp(-2 * alpha * beta * t) ) + (alpha^2 * v0^2 * exp(-alpha * beta * t));
        EI(i,j) = alpha * exp(-alpha * beta * t) * v0;
    end
    I(1,j) = I(1) + 2 * k * T * alpha;
end


%% plots
hold on
% plot(I*1e5)
plot(Tspan,I(:,1)*1e5,'r','LineWidth',3)
plot(Tspan,EI(:,1),'b','LineWidth',3)
plot(Tspan,I(:,2)*1e5,'r:','LineWidth',3)
plot(Tspan,EI(:,2),'b:','LineWidth',3)
xlabel('Time \rightarrow','FontSize',48);
ylabel('Current variance and transient speed ','FontSize',48);
title('C = 1e-6, v0 = 5V, T = 300K','FontSize',48)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 48)
legend({'\tau = 0.4, Variance (Accuracy)','\tau = 0.4, Speed','\tau = 1, Variance (Accuracy)','\tau = 1, Speed'},'FontSize',48);
