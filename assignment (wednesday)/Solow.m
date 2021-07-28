clear;

% Declare parameters

s = 0.1;
alpha = 1/3;
g = 0.02;
eta = 0.01;
delta = 0.3;

kss = ((g+eta+delta)/s)^(1/(alpha-1));

% Shock to parameters

sp = 0.4;
alphap = alpha;
gp = g;
etap = eta;
deltap = delta;

kssp = ((gp+etap+deltap)/sp)^(1/(alphap-1));

% Solve first using ODE

fun = @(t,k) sp*k^alphap-(gp+etap+deltap)*k;
options = odeset('RelTol',1e-6,'AbsTol',[1e-6])
[time,capital] = ode45(fun,[0 100],kss,options);

% Done. Now assume that the economy was in the steady state initially.

T = 40;

timem = linspace(-T,-1,T-1)';
capitalm = kss*ones(size(timem));

% Producitivy before and after change

Am = exp(g*timem);
Ap = exp(gp*time);
A = [Am;Ap];

% Full time path

time = [timem;time];

% Full capital path.

capital = [capitalm;capital];

% Log capital and output per worker

log_Capitalm = log(capital(1:T-1))+log(A(1:T-1));
log_Capitalp = log(capital(T:end))+log(A(T:end));
log_Capital = [log_Capitalm;log_Capitalp];

log_Outputm = log(A(1:T-1))+alpha*log(capital(1:T-1));
log_Outputp = log(A(T:end))+alphap*log(capital(T:end));
log_Output = [log_Outputm;log_Outputp];

% Log consumption and investment per worker

Cm = log_Output(1:39)+log(1-s);
Cp = log_Output(40:end)+log(1-sp);

log_Consumption = [Cm;Cp];

Im = log_Output(1:39)+log(s);
Ip = log_Output(40:end)+log(sp);

log_Investment = [Im;Ip];

% Plot the results
figure;

subplot(2,2,1)
plot(time,log_Output)
ylabel('Output','interpreter','latex')
set(gca,'FontSize',10,'TickLabelInterpreter','latex')
axis tight

subplot(2,2,2)
plot(time,log_Capital)
ylabel('Capital','interpreter','latex')
set(gca,'FontSize',10,'TickLabelInterpreter','latex')
axis tight

subplot(2,2,3)
plot(time,log_Investment)
ylabel('Investment','interpreter','latex')
xlabel('Time (years)','interpreter','latex')
set(gca,'FontSize',10,'TickLabelInterpreter','latex')
axis tight

subplot(2,2,4)
plot(time,log_Consumption)
ylabel('Consumption','interpreter','latex')
xlabel('Time (years)','interpreter','latex')
set(gca,'FontSize',10,'TickLabelInterpreter','latex')
axis tight

set(gcf,'NextPlot','add');
axes;
h = title('\textbf{Solow growth model: Saving like China}','interpreter','latex','FontSize',10);
set(gca,'Visible','off');
set(h,'Visible','on');

xSize = 13; 
ySize = 8.5;

set(gcf,'Units','centimeters','Position',[0 0 xSize ySize],'PaperUnits','centimeters' ...
     ,'PaperPosition',[0 0 xSize ySize],'PaperSize',[xSize-1 ySize],'PaperPositionMode','auto')
 
print -dpdf solow.pdf