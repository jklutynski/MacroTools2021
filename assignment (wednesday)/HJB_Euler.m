clear;

% Declare parameters

delta = 0.05;
rho = 0.03;
gamma = 2;
alpha = 0.3;

% Grid

N = 101;
kss = (alpha/(rho+delta))^(1/(1-alpha));
kmin = 0.5*kss;
kmax = 1.5*kss;
k = linspace(kmin,kmax,N)';
dk = k(2)-k(1);
y = k.^(alpha)-delta*k;

% First solve the problem using iterations on the HJB equation.

% Initial guess

v0 =  (((y).^(1-gamma))./(1-gamma))./rho;

metric = 1;
Delta = 0.2*dk;
tic
while metric>1e-9
    
    dv0 = gradient(v0)./dk;
    c0 = dv0.^(-1/gamma);
    u0 = (c0.^(1-gamma))/(1-gamma);
    v1 = Delta*(u0+dv0.*(y-c0)-rho*v0)+v0;
    metric = max(abs(v0-v1));
    v0 = v1;
    
end
toc

c0_VFI = c0;

% Now solve the problem using the Euler equation.

% Initial guess

R = alpha*k.^(alpha-1)-delta-rho;
c0 = 5*y;
dc0 = gradient(c0)./dk;
metric = 1;
Delta = 0.001;

tic
while metric>1e-6
    
    c1 = (dc0.*y)./(1/gamma*R+dc0);   
    metric = max(abs(c1-c0));
    
    dc0 = Delta*gradient(c1)./dk+(1-Delta)*dc0;
    c0 = Delta*c1+(1-Delta)*c0;
    
end
toc

% Plot the results

figure;
y1 = plot(k,(y-c0_VFI),'-','LineWidth',2);
hold on
plot(k,zeros(1,N),'--','LineWidth',2,'color',[0 0 0]);
y2 = plot(k,(y-c0),'ro','markersize',8);
grid
xlabel('$k_t$','interpreter','latex')
ylabel('$\dot{k_t}$','interpreter','latex')
xlim([kmin kmax])
legend1 = legend([y1 y2],'Value function iteration','Euler equation iteration');
set(legend1,'Location','southwest','interpreter','latex');
set(gca,'FontSize',8,'TickLabelInterpreter','latex')

xSize = 11; 
ySize = 7.5;

set(gcf,'Units','centimeters','Position',[0 0 xSize ySize],'PaperUnits','centimeters' ...
     ,'PaperPosition',[0 0 xSize ySize],'PaperSize',[xSize ySize],'PaperPositionMode','auto')
 
print -dpdf policy_f.pdf

% Plot an impulse response

% The function will be used for the ODE solver.
fun1 = @(t,x) pchip(k,y-c0,x);
fun2 = @(t,x) pchip(k,y-c0_VFI,x);

% Compute impulse responses for 12 quareters.
[T1,Y1] = ode45(fun1,[0 100],kss*0.95);
[T2,Y2] = ode45(fun2,[0 100],kss*0.95);

% Plot the results.
figure;
y1 = plot(T1,Y1,'r','linewidth',2);
hold on
y2 = plot(T2,Y2,'b','linewidth',2);
grid;
title('Comparison','interpreter','latex')
xlabel('Time (quarters)','interpreter','latex')
ylabel('Capital','interpreter','latex')
legend1 = legend([y1 y2],'Euler equation iteration','Value function iteration');
set(legend1,'Location','southeast','interpreter','latex');
set(gca,'FontSize',8,'TickLabelInterpreter','latex')

set(gcf,'Units','centimeters','Position',[0 0 xSize ySize],'PaperUnits','centimeters' ...
     ,'PaperPosition',[0 0 xSize ySize],'PaperSize',[xSize ySize],'PaperPositionMode','auto')
 
print -dpdf impulse.pdf


