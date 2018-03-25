%%% With inifinite Av
s = tf('s');
tau0 =0.5;
C0 = 1e-5;

K1 = 10;
tau1 = 1/K1;
C1 = 1e-6;

R1 = tau1 / C1;
R = tau0/C0;

tau = C1 * (R + R1);

a = 1/tau0 + 1/(R1*C0);
% 
A = [-a -1/tau0; 1/tau1 0];
Bref = [a ; -1/tau1];
% Bid = [R/tau0;0];
C = [1 0];
D = 0;

sys = ss(A,Bref,C,D);
% [zv,pv,kv] = ss2zp(A,Bref,C,D);
% [zi,pi,ki] = ss2zp(A,Bid,C,D);
% tfVoutVref = (s*tau + 1)/(tau1*tau0*s^2 + s *tau + 1);
% tfVoutid = (s*R*tau1)/(tau1*tau0*s^2 + s*tau + 1);
% zv1 = zero(tfVoutVref);
% pv1 = pole(tfVoutVref);
% 
% zi1 = zero(tfVoutid);
% pi1 = pole(tfVoutid);



%%% With finite Av

Av = 40000;
b = (1/(1-Av)) * ( (Av/tau0) + (1/(R1*C0)) );
Aa = [-a b;1/tau1 1/((1-Av)*tau1)];
% Brefa = [-( (Av/tau0) + ( Av^2 / ( (1-Av)*tau0 ) ) + ( Av/(R1 *(1-Av)) ) ); Av/ ((1-Av)*tau1) ];
Brefa = [(-Av*a)/(1-Av);Av/( (1-Av)*tau1 )];
% Bida = [R/tau0;0];
Ca = [1 0];
Da = 0;

sys1 = ss(Aa,Brefa,Ca,Da);
[zva,pva,kva] = ss2zp(Aa,Brefa,Ca,Da);
% den = (s^2 * tau0 * tau1) + s*tau + 1 - (s*tau0)/(Av - 1);
% den = (s^2 * tau0 * tau1) + s*tau - (s*tau0)/(Av - 1) - (s * ( R*tau0 + R1*tau0*Av ) / (1-Av) ) + ((R1 + R) / (R1*(1-Av)) ) ;
% tfVoutVrefa = (Av * (s*tau + 1)) / ( (Av - 1) * den);
% sys2 = ss(Aa,Bida,Ca,Da);
% [zia,pia,kia] = ss2zp(Aa,Bida,Ca,Da);
% tfVoutida = ((s*tau1*R) + ( R / (1-Av) )  ) / den;

% zv1a = zero(tfVoutVrefa);
% pv1a = pole(tfVoutVrefa);

% zi1a = zero(tfVoutida);
% pi1a = pole(tfVoutida);

%%% plots 
% hold on
% grid
% % step(sys)
% % step(sys1)
% plot(step(tfVoutVref),'r')
% plot(step(tfVoutVrefa),'b')


%%%  solve state space equations
tspan = 0:0.01:4;
iniCon = [0;0];
Vref = 1;
Id = 0;
[i, x] = ode45(@(i,x)sys_solve(i,x,Aa,Brefa,Vref), tspan, iniCon);
% [i, x] = ode45(@(i,x)sys_solve(i,x,Aa,[Brefa Bida],[Vref;Id]), tspan, iniCon);
% y = C*x' + D*Vref;
Vc = x(:,1);
Vc1 = x(:,2);
Vc = Vc';
Vc1 = Vc1';
Vrefn = zeros(size(Vc)) + Vref;
Idn = zeros(size(Vc)) + Id;
% Generate step function
% Vrefn = [zeros(1, floor(size(Vc,2)/2)) Vref + zeros(1, floor(size(Vc,2)/2)) 1];
V0 = (Av/(1-Av)) * (Vc1 - Vrefn);
% V0 = Vrefn - Vc1;
V1 = (1/(1-Av)) * (Vc1 - Av * Vrefn);
% V1 = Vrefn;
Ir = ( (V0 - Vc) / R );
I = ( (Vc - V1) / R1 );
I0 = Ir  - I;
PowerT = V0 .* I0;
Powersat = PowerT/1;
V0sat = Powersat ./ I0;
Vc1sat = ( (1-Av)/(Av) ) .* V0sat + Vrefn;
V1sat = (V0sat ./ Av) + Vrefn;
Vcsat = ( ( (R1/R) .* V0sat ) - (R1 .* I0) + V1sat) ./ ( (R1/R) + 1);
% supplyV = 1:0.01:10;
supplyV = 1.01;
for k = 1:length(supplyV)
    supply = supplyV(k);
    MAX = max(PowerT)/supply;
    count = 0; 
    for i = 1:length(PowerT)
        if PowerT(i) > MAX
            Powersat(i) = MAX;
            count = count + 1;
        else
            Powersat(i) = PowerT(i);
        end
        V0sat(i) = Powersat(i) / I0(i);
        Vc1sat(i) = ( (1-Av)/(Av) ) * V0sat(i) + Vrefn(i);
        V1sat(i) = (V0sat(i) / Av) + Vrefn(i);
    %     Vcsat(i) = ( ( (R1/R) * V0sat(i) ) - (R1 * I0(i)) + V1sat(i) ) / ( (R1/R) + 1);
        k1 = R1 / (R1 + R);
        k2 = (R1 * R) / (R1 + R);
        Vcsat(i) = V1sat(i) - k1 * Vc1sat(i) - k2 * I0(i);
    %     Vc(i)
    %     Vcsat(i)
    end
    error(k) = norm(Vc - Vcsat);
end

%%%% V0 saturate
% PowerT = V0 .* I0;
% Powersat = PowerT/1;
% V0sat = Powersat ./ I0;
% Vc1sat = ( (1-Av)/(Av) ) .* V0sat + Vrefn;
% V1sat = (V0sat ./ Av) + Vrefn;
% Vcsat = ( ( (R1/R) .* V0sat ) - (R1 .* I0) + V1sat) ./ ( (R1/R) + 1);
% supplyV = 1:0.01:10;
% supplyV = 2;
% for k = 1:length(supplyV)
%     supply = supplyV(k);
%     MAX = max(V0)/supply;
%     count = 0; 
%     for i = 1:length(V0)
%         if V0(i) > MAX
%             V0sat(i) = MAX;
%             count = count + 1;
%         else
%             V0sat(i) = V0(i);
%         end
% %         V0sat(i) = Powersat(i) / I0(i);
%         Vc1sat(i) = ( (1-Av)/(Av) ) * V0sat(i) + Vrefn(i);
%         V1sat(i) = (V0sat(i) / Av) + Vrefn(i);
%     %     Vcsat(i) = ( ( (R1/R) * V0sat(i) ) - (R1 * I0(i)) + V1sat(i) ) / ( (R1/R) + 1);
%         k1 = R1 / (R1 + R);
%         k2 = (R1 * R) / (R1 + R);
%         Vcsat(i) = V1sat(i) - k1 * Vc1sat(i) - k2 * I0(i);
%     %     Vc(i)
%     %     Vcsat(i)
%     end
%     error(k) = norm(Vc - Vcsat);
% end
% 
kk = 0.2;
sysidp = idtf([kva (kva*zva - Av*kva*kk)/supplyV Av*kva*kk*zva/supplyV],[1 -(pva(1) + pva(2)) pva(1)*pva(2)]);
dat = iddata(Vcsat',Vrefn',0.01);
syss = tfest(dat,sysidp);
zero_est = zero(syss)
zva
pole_est = pole(syss)
pva
Yy1 = step(sys1,tspan);
Yyest = step(syss,tspan);
errr = norm(Yy1 - Yyest)


% % 
%% 
hold on 
% plot(Yy1)
% plot(Yyest,'r')
% plot(Vcsat,'g--')
grid on
plot(tspan,Vrefn,'LineWidth',4)
plot(tspan,Vc,'k','LineWidth',4)
% plot(tspan,Vcsat,'r--','LineWidth',4)
% title('Step voltage reference tracking','FontSize',36)
xlabel('Time','FontSize',48)
ylabel('Capacitor voltage','FontSize',48)
% legend({'reference','active','passive'},'FontSize',42);
legend({'reference','Vc'},'FontSize',48);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 48)

t1 = annotation('textbox',[.15 .48 0.1 0.01],'String','-','FitBoxToText','on','EdgeColor','none');
t1.FontSize = 48;
t2 = annotation('textbox',[.25 .8 0.1 0.01],'String','+','FitBoxToText','on','EdgeColor','none');
t2.FontSize = 48;
t3 = annotation('textbox',[.39 .62 0.1 0.01],'String','-','FitBoxToText','on','EdgeColor','none');
t3.FontSize = 48;
% t3 = annotation('textbox',[.22 .68 0.1 0.01],'String','+','FitBoxToText','on','EdgeColor','none');
% t3.FontSize = 28;

% plot(supplyV,error,'LineWidth',3);
% xlabel('Power supply','FontSize',36);
% ylabel('2-norm error in step tracking','FontSize',36);
%%

% % [j,Vcsat] = ode45(@(j,Vcsat)newsys(j,Vcsat,tau1,Av,Vref,Vc1sat),tspan,iniconVc);
% B = (1-Av)*s*tau1;
% Vrefn = zeros(size(Vc)) + Vref;
% % Vrefn = Vrefn';
% V0 = (Vc/(B+1)) + ( (-Av*s*tau1 * Vrefn)/(B+1) ) + (Vc1);
% I0 = (V0/R) - (Vc/R) - (s*C1*Vc1);
% % 
% iniConh = [0;0;sqrt(2*E0)];
% [i, xh] = ode45(@(i,xh)sysh_solve(i,xh,A,B,C,D,vd,E0), tspan, iniConh);
% yh = (xh(3)/sqrt(2*E0)) * (C * [xh(1);xh(2)] + vd);
% 
% xhat = xh(:,1:2);
% error(k) = norm(x-xhat);


function dx = sys_solve(i,x,Aa,Brefa,Vref)
    dx = Aa*x + Brefa*Vref;
end

% % function dx = newsys(j,Vcsat,tau1,Av,Vref)
% %     dx = V
% % end

% function dxh = sysh_solve(i,xh,A,B,C,D,vd,E0)
%     dxv = (xh(3)*A*[xh(1);xh(2)]) / sqrt(2*E0) + xh(3)*B*vd/sqrt(2*E0);
%     dxh(1) = dxv(1); dxh(2) = dxv(2);
%     dxh(3) = (1/sqrt(2*E0)) * ([xh(1)' xh(2)']*C' + D'*vd) * vd - (1/sqrt(2*E0)) * [xh(1)' xh(2)']* (A*[xh(1);xh(2)] + B*vd);
%     dxh = dxh';
% end
    
