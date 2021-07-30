clear;

% Parameters

sigma = 2;              % Coefficient of relative risk-aversion.
beta = 1.05^(-1/4);     % Discount factor.
alpha = 0.3;            % Capital share of income.
delta = 0.025;          % Depreciation rate

% Transition matrix: Symmetric just for simplicity.

rho = 0.95;

T = [rho, 1-rho;1-rho,rho];

[a,b]=eigs(T',1);

TT = a./sum(a); 

% Productivity

Z = [1.02;0.98];

% Steady state capital.

kss = ((1/beta+delta-1)/alpha)^(1/(alpha-1));

% Set up the two Euler equations symbolically. Looks messy, but it's
% actually quite repetitive.

syms km k kpg kpb

Euler_g = [-(Z(1)*km^(alpha)+(1-delta)*km-k)^(-sigma)+beta*( ...
    T(1,1)*(1+alpha*Z(1)*k^(alpha-1)-delta)*(Z(1)*k^(alpha)+(1-delta)*k-kpg)^(-sigma)+ ...
    T(1,2)*(1+alpha*Z(2)*k^(alpha-1)-delta)*(Z(2)*k^(alpha)+(1-delta)*k-kpb)^(-sigma))];

Euler_b = [-(Z(2)*km^(alpha)+(1-delta)*km-k)^(-sigma)+beta*( ...
    T(2,1)*(1+alpha*Z(1)*k^(alpha-1)-delta)*(Z(1)*k^(alpha)+(1-delta)*k-kpg)^(-sigma)+ ...
    T(2,2)*(1+alpha*Z(2)*k^(alpha-1)-delta)*(Z(2)*k^(alpha)+(1-delta)*k-kpb)^(-sigma))];

% Linearize them. We will end up with a system like A*(k-kss)+B*(kp-kss)+C_g*(kpp_g-kss)+C_b*(kpp_b-kss)+D = 0.
% Let's start with Euler_g

Xss = [kss kss kss kss];
Xvar = [km k kpg kpb];

A_g = jacobian(Euler_g,km); A_g = double(subs(A_g,Xvar,Xss));
B_g = jacobian(Euler_g,k); B_g = double(subs(B_g,Xvar,Xss));
C_gg = jacobian(Euler_g,kpg); C_gg = double(subs(C_gg,Xvar,Xss));
C_gb = jacobian(Euler_g,kpb); C_gb = double(subs(C_gb,Xvar,Xss));

D_g = double(subs(Euler_g,Xvar,Xss));

% Done. Let's do the same with Euler_b.

A_b = jacobian(Euler_b,km); A_b = double(subs(A_b,Xvar,Xss));
B_b = jacobian(Euler_b,k); B_b = double(subs(B_b,Xvar,Xss));
C_bg = jacobian(Euler_b,kpg); C_bg = double(subs(C_bg,Xvar,Xss));
C_bb = jacobian(Euler_b,kpb); C_bb = double(subs(C_bb,Xvar,Xss));

D_b = double(subs(Euler_b,Xvar,Xss));

% Ok, very good. Now the policy functions will look like (kp-kss) =
% E_{i}+F_{i}*(k-kss) for i=g,b. Let's for simplicity call (k_t-kss) for
% u_{t-1}. Then our system is a little more compactly written as
% A*u_{t-1}+B*u_{t}+C_g*u_{t+1,g}+C_b*u_{t+1,b}+D = 0. And our policy
% function is u_t=E+F*u_{t-1}. If we insert the policy function in the system
% and solve, we find that:

% E_i = (B_i+C_{i,g}F_g+C_{i,b}F_b)^(-1)*(-(C_{i,g}E_g+C_{i,b}E_b+D))
% F_i = (B_i+C_{i,g}F_g+C_{i,b}F_b)^(-1)*(-A_i)

% This is a system we can iterate on. Very similar to time iteration. Very
% similar to solving a Riccati equation.

% Initial guesses

E_g = 0;
E_b = 0;

F_g = 0;
F_b = 0;

metric = 1;

% Let's rock'n'roll

tic

while metric>1e-13
    
    % First calculate new guesses given the previous:

    E_g1 = inv(B_g+C_gg*F_g+C_gb*F_b)*(-(C_gg*E_g+C_gb*E_b+D_g));
    E_b1 = inv(B_b+C_bg*F_g+C_bb*F_b)*(-(C_bg*E_g+C_bb*E_b+D_b));

    F_g1 = inv(B_g+C_gg*F_g+C_gb*F_b)*(-A_g);
    F_b1 = inv(B_b+C_bg*F_g+C_bb*F_b)*(-A_b);
    
    % Then update the guesses for the next iteration:

    E_g = E_g1;
    E_b = E_b1;
    F_g = F_g1;
    F_b = F_b1;
    
    % Check if the Euler equation is close to zero:

    metric = max(abs([A_g+(B_g+C_gg*F_g+C_gb*F_b)*(E_g+F_g)+(C_gg*E_g+C_gb*E_b+D_g) A_b+(B_b+C_bg*F_g+C_bb*F_b)*(E_b+F_b)+(C_bg*E_g+C_bb*E_b+D_b)]));

end

timer = toc;

% Done! You probably have noticed how fast this is. If not, just type
% 'timer' into the command window and you'll see how long time the iterations
% took (measured in seconds). This is, apart from some stability properties,
% the true beauty of linearized models (at least according to me).

% Now, using these policy functions we can do a lot of stuff, such as stochastic
% simulations and impulse responses. I will focus on impulse responses. As
% discussed in class, I will calculate both the impulse responses from a
% range of predetermined durations of a shock, and for the expected value.
% Surprisingly, it is often more tedious to work out the impulse responses
% than it is to actually solve the model. So from hereon the code is all
% about computing and plotting results.

% Set the length of the impulse response.

T_imp = 50;

% Find initial state. It's not obvious that this is actually kss.

Ex_g = 0;

Ex_b = 0;

v = [0,1];

Ex = v(1)*Ex_g+v(2)*Ex_b;

metric = 1;

while metric>1e-13
    
    vn = v*T;
    
    Ex_gn = T(1,1)*v(1)/vn(1)*(E_g+F_g*Ex_g)+T(2,1)*v(2)/vn(1)*(E_g+F_g*Ex_b);
    
    Ex_bn = T(1,2)*v(1)/vn(2)*(E_b+F_b*Ex_g)+T(2,2)*v(2)/vn(2)*(E_b+F_b*Ex_b);
    
    metric = abs(vn(1)*Ex_gn+vn(2)*Ex_bn-(v(1)*Ex_g+v(2)*Ex_b));
    
    Ex_g = Ex_gn;
    
    Ex_b = Ex_bn;
    
    v = vn;
    
end

% Initial steady state is uss

uss = vn(1)*Ex_gn+vn(2)*Ex_bn;

% Now calculate the impulse response

Ex_g = uss;

Ex_b = uss;

Ex_z = 1;

v = [0,1];

Ex = v(1)*Ex_g+v(2)*Ex_b;

for i = 1:T_imp
    
    v(i+1,:) = v(i,:)*T;
    
    Ex_g(:,i+1) = T(1,1)*v(i,1)/v(i+1,1)*(E_g+F_g*Ex_g(:,i))+T(2,1)*v(i,2)/v(i+1,1)*(E_g+F_g*Ex_b(:,i));
    
    Ex_b(:,i+1) = T(1,2)*v(i,1)/v(i+1,2)*(E_b+F_b*Ex_g(:,i))+T(2,2)*v(i,2)/v(i+1,2)*(E_b+F_b*Ex_b(:,i));
    
    Ex_z(:,i+1) = v(i,1)*Z(1)+v(i,2)*Z(2);
    
    Ex(:,i+1) = v(i,1)*Ex_g(:,i+1)+v(i,2)*Ex_b(:,i+1);
    
end

Ex = [uss,Ex];

Ex_capital = Ex+kss;

y_ss = (kss+uss)^(alpha);

Ex_output = y_ss*Ex_z+alpha*(kss+uss)^(alpha-1)*(Ex_capital(1:end-1)-kss-uss);

Ex_investment = kss+Ex(2:end)-(1-delta)*(kss+Ex(1:end-1));

Ex_consumption = Ex_output-Ex_investment;

% Let's simulate a bunch of paths

No = 40;

rng(1979)
e = rand(T_imp,No);
s = 2*ones(T_imp,No);

E = [E_g,E_b];
F = [F_g,F_b];

x = ones(T_imp,No)*uss;

S = ones(1,No);

for j = 1:No

    for i = 1:T_imp-1

        x(i+1,j) = E(s(i,j))+F(s(i,j))*x(i,j);

        if e(i,j)<T(s(i,j),1)
            s(i+1,j) = 1;
        else
            s(i+1,j) = 2;
        end
        
        S(j) = S(j)*T(s(i,j),s(i+1,j));

    end

end

Sp = S;
Sp = Sp/sum(Sp);
S = S/(T(2,2)^(T_imp-1));
S = -log(S);
maxS = max(S);
minS = min(S);
dif = maxS-minS;
S = (S-minS)./dif;


capital = kss+x;
z = Z(s);
output = y_ss*z+alpha*(kss+uss)^(alpha-1)*(capital-kss-uss);
investment = capital(2:end,:)-(1-delta)*capital(1:end-1,:);
consumption = output(1:end-1,:)-investment;

y_ss = (kss+uss)^(alpha);
i_ss = kss+uss-(1-delta)*(kss+uss);
c_ss = y_ss - i_ss;


figure;
subplot(2,2,1);
for j = 1:No
    plot(100*(log(output(:,j))-log(y_ss)),'color',0.1+min(S(j)*2,1-0.1)*[1 1 1]);
    hold on
end
plot(100*(log(Ex_output(1:end-1))-log(y_ss)),'LineWidth',1.8,'color',[     0    0.4470    0.7410]);
hold off
title('Output','FontSize',12,'fontname','times')
set(gca,'FontSize',12,'fontname','times')
set(gca,'xlim',[0 50])

subplot(2,2,2);
for j = 1:No
    plot(100*(log(consumption(:,j))-log(c_ss)),'color',0.05+min(S(j)*1.5,1-0.05)*[1 1 1]);
    hold on
end
plot(100*(log(Ex_consumption(1:end-1))-log(c_ss)),'LineWidth',1.8,'color',[     0    0.4470    0.7410]);
hold off
title('Consumption','FontSize',12,'fontname','times')
set(gca,'FontSize',12,'fontname','times')
set(gca,'xlim',[0 50])

subplot(2,2,3);
for j = 1:No
    plot(100*(log(investment(:,j))-log(i_ss)),'color',0.05+min(S(j)*1.5,1-0.05)*[1 1 1]);
    hold on
end
plot(100*(log(Ex_investment(1:end-1))-log(i_ss)),'LineWidth',1.8,'color',[     0    0.4470    0.7410]);
hold off
title('Investment','FontSize',12,'fontname','times')
xlabel('Time (quarters)','FontSize',12,'fontname','times')
set(gca,'FontSize',12,'fontname','times')
set(gca,'xlim',[0 50])

subplot(2,2,4);
for j = 1:No
    plot(100*(log(z(:,j))),'color',0.05+min(S(j)*1.5,1-0.05)*[1 1 1]);
    hold on
end
plot(100*(log(Ex_z(1:end-1))),'LineWidth',1.8,'color',[     0    0.4470    0.7410]);
hold off
title('Productivity','FontSize',12,'fontname','times')
xlabel('Time (quarters)','FontSize',12,'fontname','times')
set(gca,'FontSize',12,'fontname','times')
set(gca,'xlim',[0 50])

xSize = 18; 
ySize = 13;

set(gcf,'Units','centimeters','Position',[0 0 xSize ySize],'PaperUnits','centimeters' ...
     ,'PaperPosition',[0 0 xSize ySize],'PaperSize',[xSize-2 ySize],'PaperPositionMode','auto')

 
print -dpdf Ramsey.pdf

