function g1 = dynamic_g1(T, y, x, params, steady_state, it_, T_flag)
% function g1 = dynamic_g1(T, y, x, params, steady_state, it_, T_flag)
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
%   g1
%

if T_flag
    T = model_agg_uncertainty.dynamic_g1_tt(T, y, x, params, steady_state, it_);
end
g1 = zeros(6, 13);
g1(1,6)=(-params(12));
g1(1,3)=(-(params(11)*T(3)));
g1(1,7)=getPowerDeriv(y(7),params(8),1);
g1(2,6)=(-(params(8)*T(2)));
g1(2,3)=(-(params(8)*y(6)*T(4)));
g1(2,8)=1;
g1(3,6)=(-(T(1)*(1-params(8))));
g1(3,3)=(-(y(6)*(1-params(8))*T(3)));
g1(3,9)=1;
g1(4,4)=(-1)/(y(4)*y(4));
g1(4,10)=(-((1+y(11)-params(4))*(-params(2))/(y(10)*y(10))));
g1(4,5)=(-(params(6)*(-params(5))*exp((-params(5))*y(5))));
g1(4,11)=(-(params(2)/y(10)));
g1(5,4)=1;
g1(5,1)=(-(y(8)+1-params(4)));
g1(5,5)=1;
g1(5,8)=(-y(1));
g1(5,9)=(-(1+x(it_, 1)));
g1(5,12)=(-y(9));
g1(6,2)=(-params(9));
g1(6,6)=1;
g1(6,13)=(-1);

end
