%**************************************************************************
%                   LSE Macroeconomics Summer Program
%                   Part II: Heterogeneous Agents
%                   Instructor: Wouter J. Den Haan
%
%                   use of this program in any fee-based program requires
%                   explicit permission (wjdenhaan@gmail.com)
%**************************************************************************
%
% Solving a model with both aggregate and invividual uncertainty using the
% algorithm of Krusell and Smith (1998).
%
%**************************************************************************
%
% notes: 
% 1) The problem of the individuals is solved using the dynare program "model_agg_uncertainty.mod" (second-order).
% 2) Simulated variables are all in absolute deviations from their own steady state levels.
% 3) To run this program, "disp_dr.m" is needed.
%
%**************************************************************************

warning off all
close all

% settings
T=10000;            % number of time periods in simulation step
N=1000;             % number of agents in simulation step
crit=0.000001;      % convergence criterion    

% parameter values
alpha = 0.25;
nu    = 1;
delta = 0.025;
beta  = 0.99;
rho   = 0.95;
k_ss  = (beta*alpha/(1-beta*(1-delta)))^(1/(1-alpha));
zeta0 = 0.01;
zeta1 = 0.3;
zeta2 = zeta1*exp(-zeta0*k_ss);
sig_e1   = 0.1;             
sig_e2   = 0.007;
up = 0.7;

s = RandStream('mt19937ar','seed',20140806);
RandStream.setGlobalStream(s);


% initial values for coefficients law of motion for aggregate capital 
b_z=0.95;                            
b_0=1.7;
b_K=0.9;
save parametervalues k_ss beta nu delta zeta0 zeta1 zeta2 alpha rho b_0 b_K b_z sig_e1 sig_e2

coef_old = [b_0 b_K b_z]';

% reserve memory
k_sim=zeros(N,T);                   % reserve memory for simulated individual capital stocks
k_sim(:,1)=0;                       % all agents start the simulation with the steady state level of capital
Ka_sim=zeros(T,1);                  % reserve memory for simulated aggregate capital stock
Ka_sim(1)=mean(k_sim(:,1));
Ka_sim_temp = zeros(T,1); 
z_sim=zeros(T,1);                   % reserve memory for aggregate productivity shock
z_sim(1)=0;
error=100;                          % initial value for error value in regression step

% draw shocks and simulate productivity
e1_sim=sig_e1*randn(N,T);               % draw innovations idiosyncratic productivity shock
e2_sim=sig_e2*randn(T,1);               % draw innovations aggregate productivity shock
for t=2:T
   z_sim(t) = 1 - rho + rho*z_sim(t-1) + e2_sim(t);
end

while error > crit
    % step 1: solve problem individual given a perceived law of motion for aggregate capital stock Ka
    save aggregatelaw b_0 b_K b_z sig_e1
    dynare model_agg_uncertainty.mod noclearall    % this generates the policy rule
    load dynarerocks                               % this uploads the decision rule into the matrix called decision

    % step 2: simulate economy, consisting of N agents, over T periods
    for t=2:T
        % compute individual capital using decision rules from dynare

        k_sim(:,t)=[ones(N,1), ...
            k_sim(:,t-1), ...
            ones(N,1)*(z_sim(t-1) - 1), ...
            ones(N, 1)*Ka_sim(t-1), ...
            e1_sim(:, t), ...
            ones(N, 1)*e2_sim(t), ...
            k_sim(:,t-1).^2, ...
            (z_sim(t-1) - 1)*k_sim(:, t-1), ...
            ones(N,1)*(z_sim(t-1) - 1).^2, ...
            Ka_sim(t-1)*k_sim(:, t-1), ...
            ones(N, 1)*Ka_sim(t-1)*(z_sim(t-1) - 1), ...
            ones(N, 1)*Ka_sim(t-1).^2, ...
            e1_sim(:, t).^2, ...
            e1_sim(:, t)*e2_sim(t), ...
            ones(N, 1)*e2_sim(t).^2, ...
            k_sim(:,t-1).*e1_sim(:, t), ...
            k_sim(:,t-1).*e2_sim(t), ...
            (z_sim(t-1) - 1).*e1_sim(:, t), ...
            ones(N, 1)*(z_sim(t-1) - 1).*e2_sim(t), ...
            Ka_sim(t-1).*e1_sim(:, t), ...
            ones(N, 1)*Ka_sim(t-1)*e2_sim(t)]* ...
            %decision(2:end, 1);
            [decision(1, 1) ; decision(3:end,1)];    

        % compute aggregate capital, taking account of the way dynare writes
        % decision rules.
        %Ka_sim(t) = mean(k_sim(:, t), 'all');         %b_0 + b_K*Ka_sim(t-1) + b_z*(z_sim(t))
        Ka_sim(t)=mean(k_sim(:,t));
        k_sim(:, t) = k_sim(:, t) - (decision(1, 1)-decision(2,1));
    end    
    
    % step 3: run a regression and update law of motion for Ka
    Ka_sim_temp(2:end) = Ka_sim(2:end) - b_0/(1 - b_K);
    y = Ka_sim(201:T);
    x0 = ones(T-200, 1);
    x1 = Ka_sim(200:T-1);
    x2 = z_sim(201:T) - 1;
    X = [x0 x1 x2];
    [B,~] = regress(y, X);
    error = norm(B - [b_0 b_K b_z])
    b_0 = b_0*up + (1-up)*B(1);
    b_K = b_K*up + (1-up)*B(2);
    b_z = b_z*up + (1-up)*B(3);
    Ka_sim(2:end) = Ka_sim_temp(2:end);
    %k_sim(2, 2)
end

             