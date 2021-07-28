function g2 = dynamic_g2(T, y, x, params, steady_state, it_, T_flag)
% function g2 = dynamic_g2(T, y, x, params, steady_state, it_, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T             [#temp variables by 1]     double   vector of temporary terms to be filled by function
%   y             [#dynamic variables by 1]  double   vector of endogenous variables in the order stored
%                                                     in M_.lead_lag_incidence; see the Manual
%   x             [nperiods by M_.exo_nbr]   double   matrix of exogenous variables (in declaration order)
%                                                     for all simulation periods
%   steady_state  [M_.endo_nbr by 1]         double   vector of steady state values
%   params        [M_.param_nbr by 1]        double   vector of parameter values in declaration order
%   it_           scalar                     double   time period for exogenous variables for which
%                                                     to evaluate the model
%   T_flag        boolean                    boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   g2
%

if T_flag
    T = model_agg_uncertainty.dynamic_g2_tt(T, y, x, params, steady_state, it_);
end
v2 = zeros(17,3);
v2(1,1)=1;
v2(2,1)=1;
v2(3,1)=2;
v2(4,1)=2;
v2(5,1)=2;
v2(6,1)=3;
v2(7,1)=3;
v2(8,1)=3;
v2(9,1)=4;
v2(10,1)=4;
v2(11,1)=4;
v2(12,1)=4;
v2(13,1)=4;
v2(14,1)=5;
v2(15,1)=5;
v2(16,1)=5;
v2(17,1)=5;
v2(1,2)=29;
v2(2,2)=85;
v2(3,2)=68;
v2(4,2)=32;
v2(5,2)=29;
v2(6,2)=68;
v2(7,2)=32;
v2(8,2)=29;
v2(9,2)=43;
v2(10,2)=127;
v2(11,2)=128;
v2(12,2)=140;
v2(13,2)=57;
v2(14,2)=8;
v2(15,2)=92;
v2(16,2)=116;
v2(17,2)=152;
v2(1,3)=(-(params(11)*T(5)));
v2(2,3)=getPowerDeriv(y(7),params(8),2);
v2(3,3)=(-(params(8)*T(4)));
v2(4,3)=v2(3,3);
v2(5,3)=(-(params(8)*y(6)*getPowerDeriv(y(3),params(8)-1,2)));
v2(6,3)=(-((1-params(8))*T(3)));
v2(7,3)=v2(6,3);
v2(8,3)=(-(y(6)*(1-params(8))*T(5)));
v2(9,3)=(y(4)+y(4))/(y(4)*y(4)*y(4)*y(4));
v2(10,3)=(-((1+y(11)-params(4))*(-((-params(2))*(y(10)+y(10))))/(y(10)*y(10)*y(10)*y(10))));
v2(11,3)=(-((-params(2))/(y(10)*y(10))));
v2(12,3)=v2(11,3);
v2(13,3)=(-(params(6)*(-params(5))*(-params(5))*exp((-params(5))*y(5))));
v2(14,3)=(-1);
v2(15,3)=v2(14,3);
v2(16,3)=(-1);
v2(17,3)=v2(16,3);
g2 = sparse(v2(:,1),v2(:,2),v2(:,3),6,169);
end
