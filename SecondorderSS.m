t = 1; t1 = 1; E0v = 1:1:100;
A = [0 1;-1/(t1*t) -1/t];
B = [0; 1/(t*t1)];
C = [1 0];
D = 1;

vd = 1;

tspan = 0:0.1:10;
iniCon = [0;0];
for k = 1:length(E0v)
    E0 = E0v(k);
    [i, x] = ode45(@(i,x)sys(i,x,A,B,vd), tspan, iniCon);
    y = C*x' + D*vd;

    iniConh = [0;0;sqrt(2*E0)];
    [i, xh] = ode45(@(i,xh)sysh(i,xh,A,B,C,D,vd,E0), tspan, iniConh);
    yh = (xh(3)/sqrt(2*E0)) * (C * [xh(1);xh(2)] + vd);

    xhat = xh(:,1:2);
    error(k) = norm(x-xhat);
end

plot(E0v,error,'r','LineWidth',3)
xlabel('Energy stored initially','FontSize',40)
ylabel('2-Norm of Error in Approximation','FontSize',40)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 28)
hold on

function dx = sys(i,x,A,B,vd)
    dx = A*x + B*vd;
end

function dxh = sysh(i,xh,A,B,C,D,vd,E0)
    dxv = (xh(3)*A*[xh(1);xh(2)]) / sqrt(2*E0) + xh(3)*B*vd/sqrt(2*E0);
    dxh(1) = dxv(1); dxh(2) = dxv(2);
    dxh(3) = (1/sqrt(2*E0)) * ([xh(1)' xh(2)']*C' + D'*vd) * vd - (1/sqrt(2*E0)) * [xh(1)' xh(2)']* (A*[xh(1);xh(2)] + B*vd);
    dxh = dxh';
end
    