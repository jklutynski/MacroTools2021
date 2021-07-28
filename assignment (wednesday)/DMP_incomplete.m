clear;

% Declare parameters

delta = 0.05;               % Separation rate
rho = 0.03;                 % Discount rate
gamma = 0 %5;                % Inverse of EIS
eta = 0.5;                  % Elasticity of matching wrt vacancies
pi = (rho+delta);                  % Profits
n_ss = 0.94;                % Steady state employment
theta_ss = 1;               % Steady state labor market tightness
f = n_ss*delta/(1-n_ss);    % Steady state job finding rate
psi = f/(theta_ss.^(eta));  % Matching efficiency
Jss = pi/(rho + delta);               % Steady state firm value
kappa = psi*Jss;            % Vacancy posting cost
up = 0.9;

% Create the grid

N = 21;                     % This is pretty low, but my own (limited) experience is that these methods work better for coarse grids.
nmin = 0.85;
nmax = 0.97;
n = linspace(nmin,nmax,N)';
dn = n(2)-n(1);

% Initial guesses

ndot = (1-n)*f-delta*n;
c = n*(1+kappa)-kappa;
dc = 1+kappa;
dJ = 0;
J = Jss;

% Set convergence numbers

Delta = 0.01;
metric = 1;

% A counter to keep track of the number of iterations

i = 0;

% And let's do some time keeping.

tic

% Let's rock'n'roll

while metric>1e-6
    
%     i = i+1
%     %dc = gradient(c).*ndot/dn;
%     J1 = (pi - gamma.*dc./c.*J + dJ)/(delta + rho);  %1/(delta + rho + dc/c*gamma)*(pi + dJ);
%     metric = max(abs(J-J1));
%     
%     J = up* J + (1 - up)*J1;
%     theta = (kappa./(psi*J)).^(1/(eta-1));
%     ndot = (1 - n).*psi*theta.^eta - n*delta;
%     c = n - kappa*theta.*(1 - n);
%     dc = gradient(c).*ndot/dn;
%     dJ = gradient(J).*ndot/dn;

    i = i+1;
    
    J1 = 1/(rho+delta)*(pi+dJ.*ndot-gamma*dc.*ndot./c.*J);
    metric = max(abs(J-J1));
    
    J = Delta*(pi+dJ.*ndot-gamma*dc.*ndot./c.*J-(rho+delta)*J)+J;
    theta = (kappa./(J*psi)).^(1/(eta-1));
    ndot = (1-n)*psi.*theta.^(eta)-delta*n;
    c = n.*(1+kappa*theta)-kappa*theta;
    dc = gradient(c)/dn;
    dJ = gradient(J)/dn;

    
end

toc

% Ok done. That was fast. You can compare the solution to one with 251
% gridpoints if you set Delta to 0.0001 or lower. That will be very slow
% and the difference appears almost indiscernible.

% Now let's calculate an impulse response function from an unexpected job
% destruction shock that bring employment to 0.85.

% The function will be used for the ODE solver.
fun2 = @(t,x) (1 - x)*f - x*delta;
% Let's compare with the infinite EIS solution:
fun1 = @(t,x) pchip(n,ndot,x);

% Compute impulse responses for 12 quarters.
[T1,Y1] = ode45(fun1, [0 12], 0.85);
[T2,Y2] = ode45(fun2, [0 12], 0.85);

% Plot the results.

% y1 = plot(T1,Y1,'linewidth',2);
% hold on
% y2 = plot(T2,Y2,'linewidth',2);
% grid;
% title('Comparison')
% xlabel('Time (quarters)')
% ylabel('Employment')
% legend1 = legend([y1 y2],'Finite EIS','Infinite EIS');
% set(legend1,'Location','southeast');
% set(gca,'FontSize',14)
xSize = 19; 
ySize = 6.5;

% Plot the results.
figure;
subplot(1,2,1);
plot(n,J,'linewidth',2);
grid;
hold on
plot(n,Jss*ones(size(n)),'linewidth',2);
xlabel('Employment')
ylabel('Firm value')
axis tight
set(gca,'FontSize',10)

subplot(1,2,2);
y1 = plot(T1,Y1,'linewidth',2);
hold on
y2 = plot(T2,Y2,'linewidth',2);
grid;
xlabel('Time (quarters)')
ylabel('Employment')
legend1 = legend([y1 y2],'Some EIS','Another EIS');
set(legend1,'Location','southeast');
axis tight
set(gca,'FontSize',10)

set(gcf,'Units','centimeters','Position',[0 0 xSize ySize],'PaperUnits','centimeters' ...
     ,'PaperPosition',[0 0 xSize ySize],'PaperSize',[xSize-2.7 ySize+0.5],'PaperPositionMode','auto')

