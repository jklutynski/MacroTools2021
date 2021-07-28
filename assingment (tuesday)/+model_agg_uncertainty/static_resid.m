function residual = static_resid(T, y, x, params, T_flag)
% function residual = static_resid(T, y, x, params, T_flag)
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
%   residual
%

if T_flag
    T = model_agg_uncertainty.static_resid_tt(T, y, x, params);
end
residual = zeros(6, 1);
lhs = T(1);
rhs = params(10)+T(1)*params(11)+params(12)*(y(3)-1);
residual(1) = lhs - rhs;
lhs = y(5);
rhs = params(8)*y(3)*T(2);
residual(2) = lhs - rhs;
lhs = y(6);
rhs = T(1)*y(3)*(1-params(8));
residual(3) = lhs - rhs;
lhs = 1/y(1)+params(7)-params(6)*exp((-params(5))*y(2));
rhs = params(2)/y(1)*(1+y(5)-params(4));
residual(4) = lhs - rhs;
lhs = y(1)+y(2);
rhs = y(5)*y(2)+y(6)*(1+x(1))+y(2)*(1-params(4));
residual(5) = lhs - rhs;
lhs = y(3);
rhs = 1-params(9)+y(3)*params(9)+x(2);
residual(6) = lhs - rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
end
