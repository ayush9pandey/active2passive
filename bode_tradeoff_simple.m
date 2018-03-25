s = tf('s');
t = 1; %plant parameter t = RC
T1 = 1:1:5; %controller parameter t1 = R1C1
P = 1 / (t*s + 1);
for i = 1:length(T1)
    t1 = T1(i);
    C = 1 / (t1*s);
    Ls = loopsens(P,C);
    S = Ls.Si;
    [H,W] = freqresp(S);
    hold on
    plot(W,log(abs(reshape(H,size(H,3),size(H,2)))),'LineWidth',2);
end
plot(W,zeros(size(W)),'r-','LineWidth',2)
xlabel('Frequency','FontSize',40)
ylabel('log(|S(jw)|','FontSize',40)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 28)

