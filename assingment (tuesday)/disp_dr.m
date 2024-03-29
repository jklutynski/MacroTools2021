function disp_dr(dr,order,var_list)

% MODIFIED FROM THE ORIGINAL DYNARE VERSION BY 
% 1. COEFFICIENTS ARE ALWAYS WRITTEN TO SCREEN (no matter tolerance)
%    (this is done at the bottom where FLAG is defined
% 2. PUTTING WHAT IS WRITTEN TO
%    THE SCREEN IN MATRIX CALLED DECISION AND WRITING THIS MATRIX TO A FILE
%    CALLED dynarerocks.mat

% Display the decision rules
%
% INPUTS
% - dr         [struct]      decision rules.
% - order      [integer]     order of approximation.
% - var_list   [cell]        list of endogenous variables for which the decision rules should be printed.
%
% OUTPUTS1
% none

% Copyright (C) 2001-2019 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

global M_ options_

if order~=1 && M_.hessian_eq_zero
    order = 1;
    warning('disp_dr: using order = 1 because Hessian is equal to zero');
end

nx =size(dr.ghx,2);
nu =size(dr.ghu,2);
if options_.block
    k = M_.nspred;
    k1 = 1:M_.endo_nbr;
else
    k = find(dr.kstate(:,2) <= M_.maximum_lag+1);
    klag = dr.kstate(k,[1 2]);
    k1 = dr.order_var;
end

if isempty(var_list)
    var_list = M_.endo_names(1:M_.orig_endo_nbr);
end

nvar = length(var_list);

ivar=zeros(nvar,1);
for i=1:nvar
    i_tmp = strmatch(var_list{i}, M_.endo_names(k1), 'exact');
    if isempty(i_tmp)
        disp(var_list{i});
        error ('One of the variable specified does not exist') ;
    else
        ivar(i) = i_tmp;
    end
end

% get length of display strings
header_label_length=16; %default
for ii=1:length(ivar)
    header_label_length = max(header_label_length,length(M_.endo_names{k1(ivar(ii))})+2);
end
header_label_format  = sprintf('%%%ds',header_label_length);
value_format_float  = sprintf('%%%d.6f',header_label_length);
value_format_zero  = sprintf('%%%dd',header_label_length);

% account for additional characters introduced by auxiliary variables
if ~isempty(M_.aux_vars)
    aux_vars_type = [M_.aux_vars.type];
    if any(aux_vars_type==4)
        aux_var_additional_characters=14;
    else
        aux_var_additional_characters=3;
    end
else
    aux_var_additional_characters=0;
end

var_name_width = max([cellofchararraymaxlength(M_.endo_names(k1(ivar))), cellofchararraymaxlength(M_.exo_names)]);

%deal with covariances
if order > 1
    var_name_width=max(2*(var_name_width+aux_var_additional_characters)+2,20); %account for covariances, separated by comma
else
    var_name_width=max(var_name_width+aux_var_additional_characters,20);
end
label_format  = sprintf('%%-%ds',var_name_width);


%% start displayimg
disp('POLICY AND TRANSITION FUNCTIONS')
% variable names
str = char(32*ones(1,var_name_width));
for i=1:nvar
    str = [str sprintf(header_label_format, M_.endo_names{k1(ivar(i))})];
end

decision = [];

disp(str);
%
% constant
%
str=sprintf(label_format,'Constant');
flag = 0;
for i=1:nvar
    x = dr.ys(k1(ivar(i)));
    if order > 1
        x = x + dr.ghs2(ivar(i))/2;
    end
    decision = [decision x];
    [str,flag]=get_print_string(str,x,value_format_zero,value_format_float,flag,options_);
end
if flag
    disp(str)
end
if order > 1
    decisiontemp = [];
    str = sprintf(label_format,'(correction)');
    flag = 0;
    for i=1:nvar
        x = dr.ghs2(ivar(i))/2;
        [str,flag]=get_print_string(str,x,value_format_zero,value_format_float,flag,options_);
        decisiontemp = [decisiontemp x];
    end
    if flag
        disp(str)
    end
    decision = [decision;decisiontemp];
end
%
% ghx
%
for k=1:nx
    decisiontemp = [];
    flag = 0;
    if isfield(dr,'state_var')
        str1 = subst_auxvar(dr.state_var(k),-1);
    else
        str1 = subst_auxvar(k1(klag(k,1)),klag(k,2)-M_.maximum_lag-2);
    end
    str = sprintf(label_format,str1);
    for i=1:nvar
        x = dr.ghx(ivar(i),k);
        [str,flag]=get_print_string(str,x,value_format_zero,value_format_float,flag,options_);
        decisiontemp = [decisiontemp x];
    end
    if flag
        disp(str)
    end
    decision = [decision;decisiontemp];
end
%
% ghu
%
for k=1:nu
    decisiontemp = [];
    flag = 0;
    str = sprintf(label_format, M_.exo_names{k});
    for i=1:nvar
        x = dr.ghu(ivar(i),k);
        [str,flag] = get_print_string(str, x, value_format_zero, value_format_float, flag, options_);
        decisiontemp = [decisiontemp x];
    end
    if flag
        disp(str)
    end
    decision = [decision;decisiontemp];
end

if order > 1
    % ghxx
    for k = 1:nx
        for j = 1:k
            decisiontemp = [];
            flag = 0;
            str1 = sprintf('%s,%s',subst_auxvar(k1(klag(k,1)),klag(k,2)-M_.maximum_lag-2), ...
                           subst_auxvar(k1(klag(j,1)),klag(j,2)-M_.maximum_lag-2));
            str = sprintf(label_format, str1);
            for i=1:nvar
                if k == j
                    x = dr.ghxx(ivar(i),(k-1)*nx+j)/2;
                else
                    x = dr.ghxx(ivar(i),(k-1)*nx+j);
                end
                [str,flag]=get_print_string(str,x,value_format_zero,value_format_float,flag,options_);
                decisiontemp = [decisiontemp x];
            end
            if flag
                disp(str)
            end
        decision = [decision;decisiontemp];
        end
    end
    %
    % ghuu
    %
    for k = 1:nu
        for j = 1:k
            decisiontemp = [];
            flag = 0;
            str = sprintf(label_format, [M_.exo_names{k} ',' M_.exo_names{j}]);
            for i=1:nvar
                if k == j
                    x = dr.ghuu(ivar(i),(k-1)*nu+j)/2;
                else
                    x = dr.ghuu(ivar(i),(k-1)*nu+j);
                end
                [str,flag]=get_print_string(str, x, value_format_zero, value_format_float, flag, options_);
                decisiontemp = [decisiontemp x];
            end
            if flag
                disp(str)
            end
        decision = [decision;decisiontemp];
        end
    end
    %
    % ghxu
    %
    for k = 1:nx
        for j = 1:nu
            decisiontemp = [];
            flag = 0;
            str1 = sprintf('%s,%s',subst_auxvar(k1(klag(k,1)),klag(k,2)-M_.maximum_lag-2), M_.exo_names{j});
            str = sprintf(label_format,str1);
            for i=1:nvar
                x = dr.ghxu(ivar(i),(k-1)*nu+j);
                [str,flag]=get_print_string(str,x,value_format_zero,value_format_float,flag,options_);
                decisiontemp = [decisiontemp x];
            end
            if flag
                disp(str)
            end
        decision = [decision;decisiontemp];
        end
    end
end

save dynarerocks decision
end

% Given the index of an endogenous (possibly an auxiliary var), and a
% lead/lag, creates a string of the form "x(lag)".
% In the case of auxiliary vars for lags, replace by the original variable
% name, and compute the lead/lag accordingly.
function str = subst_auxvar(aux_index, aux_lead_lag)
global M_

if aux_index <= M_.orig_endo_nbr
    str = sprintf('%s(%d)', M_.endo_names{aux_index}, aux_lead_lag);
    return
end
for i = 1:length(M_.aux_vars)
    if M_.aux_vars(i).endo_index == aux_index
        switch M_.aux_vars(i).type
          case 0
            % endo leads >= 2
            str = sprintf('%s(%d)', M_.endo_names{aux_index}, aux_lead_lag);
            return
          case 1
            % endo lags >= 2
            orig_name = M_.endo_names{M_.aux_vars(i).orig_index};
          case 2
            % exo leads >= 1
            orig_name = M_.exo_names{M_.aux_vars(i).orig_index};
          case 3
            % exo lags >= 1
            orig_name = M_.exo_names{M_.aux_vars(i).orig_index};
          case 4
            % Expectation operator
            str = sprintf('EXPECTATION(%d)(...)', aux_lead_lag);
            return
          case 6
            % Ramsey's multipliers
            str = sprintf('%s(%d)', M_.endo_names{M_.aux_vars(i).endo_index}, aux_lead_lag);
            return
          case 8
            % Diff operator
            str = sprintf('diff(%s)', M_.endo_names{M_.aux_vars(i).orig_index});
            return
          case 9
            % Lagged diff
            lags = 0;
            j = i;
            while M_.aux_vars(j).type==9
                j = M_.aux_vars(j).orig_index;
                lags = lags+1;
            end
            str = sprintf('diff(%s(-%u))', M_.endo_names{M_.aux_vars(j).orig_index}, lags);
            return
          case 11
            % Leaded diff
            leads = 0;
            j = i;
            while M_.aux_vars(j).type==11
                j = M_.aux_vars(j).orig_index;
                lags = lags+1;
            end
            str = sprintf('diff(%s(%u))', M_.endo_names{M_.aux_vars(j).orig_index}, leads);
            return
          otherwise
            error('Invalid auxiliary type: %s', M_.endo_names{aux_index})
        end
        str = sprintf('%s(%d)', orig_name, M_.aux_vars(i).orig_lead_lag+aux_lead_lag);
        return
    end
end
error('Could not find aux var: %s', M_.endo_names{aux_index})
end

function [str,flag]=get_print_string(str, x, value_format_zero, value_format_float, flag, options_)
if abs(x) >= options_.dr_display_tol
    flag = 1;
    str = [str sprintf(value_format_float, x)];
else
    flag = 1;
    str = [str sprintf(value_format_float, x)];
end
end
