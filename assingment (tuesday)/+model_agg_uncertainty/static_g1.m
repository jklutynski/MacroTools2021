function g1 = static_g1(T, y, x, params, T_flag)
% function g1 = static_g1(T, y, x, params, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T         [#temp variables by 1]  double   vector of temporary terms to be filled by function
%   y         [M_.endo_nbr by 1]      double   vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1]       double   vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1]     double   vector of parameter values in declaration order
%                                              to evaluate the model
%   T_flag    boolean                 boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   g1
%

if T_flag
    T = model_agg_uncertainty.static_g1_tt(T, y, x, params);
end
g1 = zeros(6, 6);
g1(1,3)=(-params(12));
g1(1,4)=T(3)-params(11)*T(3);
g1(2,3)=(-(params(8)*T(2)));
g1(2,4)=(-(params(8)*y(3)*getPowerDeriv(y(4),params(8)-1,1)));
g1(2,5)=1;
g1(3,3)=(-(T(1)*(1-params(8))));
g1(3,4)=(-(y(3)*(1-params(8))*T(3)));
g1(3,6)=1;
g1(4,1)=(-1)/(y(1)*y(1))-(1+y(5)-params(4))*(-params(2))/(y(1)*y(1));
g1(4,2)=(-(params(6)*(-params(5))*exp((-params(5))*y(2))));
g1(4,5)=(-(params(2)/y(1)));
g1(5,1)=1;
g1(5,2)=1-(y(5)+1-params(4));
g1(5,5)=(-y(2));
g1(5,6)=(-(1+x(1)));
g1(6,3)=1-params(9);
if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
end
end
