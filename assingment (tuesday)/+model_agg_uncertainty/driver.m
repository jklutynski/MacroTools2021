%
% Status : main Dynare file
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

tic0 = tic;
% Define global variables.
global M_ options_ oo_ estim_params_ bayestopt_ dataset_ dataset_info estimation_info ys0_ ex0_
options_ = [];
M_.fname = 'model_agg_uncertainty';
M_.dynare_version = '4.6.4';
oo_.dynare_version = '4.6.4';
options_.dynare_version = '4.6.4';
%
% Some global variables initialization
%
global_initialization;
diary off;
diary('model_agg_uncertainty.log');
M_.exo_names = cell(2,1);
M_.exo_names_tex = cell(2,1);
M_.exo_names_long = cell(2,1);
M_.exo_names(1) = {'e1'};
M_.exo_names_tex(1) = {'e1'};
M_.exo_names_long(1) = {'e1'};
M_.exo_names(2) = {'e2'};
M_.exo_names_tex(2) = {'e2'};
M_.exo_names_long(2) = {'e2'};
M_.endo_names = cell(6,1);
M_.endo_names_tex = cell(6,1);
M_.endo_names_long = cell(6,1);
M_.endo_names(1) = {'c'};
M_.endo_names_tex(1) = {'c'};
M_.endo_names_long(1) = {'c'};
M_.endo_names(2) = {'k'};
M_.endo_names_tex(2) = {'k'};
M_.endo_names_long(2) = {'k'};
M_.endo_names(3) = {'z'};
M_.endo_names_tex(3) = {'z'};
M_.endo_names_long(3) = {'z'};
M_.endo_names(4) = {'Ka'};
M_.endo_names_tex(4) = {'Ka'};
M_.endo_names_long(4) = {'Ka'};
M_.endo_names(5) = {'r'};
M_.endo_names_tex(5) = {'r'};
M_.endo_names_long(5) = {'r'};
M_.endo_names(6) = {'w'};
M_.endo_names_tex(6) = {'w'};
M_.endo_names_long(6) = {'w'};
M_.endo_partitions = struct();
M_.param_names = cell(14,1);
M_.param_names_tex = cell(14,1);
M_.param_names_long = cell(14,1);
M_.param_names(1) = {'k_ss'};
M_.param_names_tex(1) = {'k\_ss'};
M_.param_names_long(1) = {'k_ss'};
M_.param_names(2) = {'beta'};
M_.param_names_tex(2) = {'beta'};
M_.param_names_long(2) = {'beta'};
M_.param_names(3) = {'nu'};
M_.param_names_tex(3) = {'nu'};
M_.param_names_long(3) = {'nu'};
M_.param_names(4) = {'delta'};
M_.param_names_tex(4) = {'delta'};
M_.param_names_long(4) = {'delta'};
M_.param_names(5) = {'zeta0'};
M_.param_names_tex(5) = {'zeta0'};
M_.param_names_long(5) = {'zeta0'};
M_.param_names(6) = {'zeta1'};
M_.param_names_tex(6) = {'zeta1'};
M_.param_names_long(6) = {'zeta1'};
M_.param_names(7) = {'zeta2'};
M_.param_names_tex(7) = {'zeta2'};
M_.param_names_long(7) = {'zeta2'};
M_.param_names(8) = {'alpha'};
M_.param_names_tex(8) = {'alpha'};
M_.param_names_long(8) = {'alpha'};
M_.param_names(9) = {'rho'};
M_.param_names_tex(9) = {'rho'};
M_.param_names_long(9) = {'rho'};
M_.param_names(10) = {'b_0'};
M_.param_names_tex(10) = {'b\_0'};
M_.param_names_long(10) = {'b_0'};
M_.param_names(11) = {'b_K'};
M_.param_names_tex(11) = {'b\_K'};
M_.param_names_long(11) = {'b_K'};
M_.param_names(12) = {'b_z'};
M_.param_names_tex(12) = {'b\_z'};
M_.param_names_long(12) = {'b_z'};
M_.param_names(13) = {'sig_e1'};
M_.param_names_tex(13) = {'sig\_e1'};
M_.param_names_long(13) = {'sig_e1'};
M_.param_names(14) = {'sig_e2'};
M_.param_names_tex(14) = {'sig\_e2'};
M_.param_names_long(14) = {'sig_e2'};
M_.param_partitions = struct();
M_.exo_det_nbr = 0;
M_.exo_nbr = 2;
M_.endo_nbr = 6;
M_.param_nbr = 14;
M_.orig_endo_nbr = 6;
M_.aux_vars = [];
M_.Sigma_e = zeros(2, 2);
M_.Correlation_matrix = eye(2, 2);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
M_.sigma_e_is_diagonal = true;
M_.det_shocks = [];
options_.linear = false;
options_.block = false;
options_.bytecode = false;
options_.use_dll = false;
options_.linear_decomposition = false;
M_.nonzero_hessian_eqs = [1 2 3 4 5];
M_.hessian_eq_zero = isempty(M_.nonzero_hessian_eqs);
M_.orig_eq_nbr = 6;
M_.eq_nbr = 6;
M_.ramsey_eq_nbr = 0;
M_.set_auxiliary_variables = exist(['./+' M_.fname '/set_auxiliary_variables.m'], 'file') == 2;
M_.epilogue_names = {};
M_.epilogue_var_list_ = {};
M_.orig_maximum_endo_lag = 1;
M_.orig_maximum_endo_lead = 1;
M_.orig_maximum_exo_lag = 0;
M_.orig_maximum_exo_lead = 0;
M_.orig_maximum_exo_det_lag = 0;
M_.orig_maximum_exo_det_lead = 0;
M_.orig_maximum_lag = 1;
M_.orig_maximum_lead = 1;
M_.orig_maximum_lag_with_diffs_expanded = 1;
M_.lead_lag_incidence = [
 0 4 10;
 1 5 0;
 2 6 0;
 3 7 0;
 0 8 11;
 0 9 0;]';
M_.nstatic = 1;
M_.nfwrd   = 2;
M_.npred   = 3;
M_.nboth   = 0;
M_.nsfwrd   = 2;
M_.nspred   = 3;
M_.ndynamic   = 5;
M_.dynamic_tmp_nbr = [2; 2; 1; 0; ];
M_.model_local_variables_dynamic_tt_idxs = {
};
M_.equations_tags = {
  1 , 'name' , '1' ;
  2 , 'name' , 'r' ;
  3 , 'name' , 'w' ;
  4 , 'name' , '4' ;
  5 , 'name' , '5' ;
  6 , 'name' , 'z' ;
};
M_.mapping.c.eqidx = [4 5 ];
M_.mapping.k.eqidx = [4 5 ];
M_.mapping.z.eqidx = [1 2 3 6 ];
M_.mapping.Ka.eqidx = [1 2 3 ];
M_.mapping.r.eqidx = [2 4 5 ];
M_.mapping.w.eqidx = [3 5 ];
M_.mapping.e1.eqidx = [5 ];
M_.mapping.e2.eqidx = [6 ];
M_.static_and_dynamic_models_differ = false;
M_.has_external_function = false;
M_.state_var = [2 3 4 ];
M_.exo_names_orig_ord = [1:2];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(6, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(2, 1);
M_.params = NaN(14, 1);
M_.endo_trends = struct('deflator', cell(6, 1), 'log_deflator', cell(6, 1), 'growth_factor', cell(6, 1), 'log_growth_factor', cell(6, 1));
M_.NNZDerivatives = [22; 17; -1; ];
M_.static_tmp_nbr = [2; 1; 0; 0; ];
M_.model_local_variables_static_tt_idxs = {
};
load parametervalues;
set_param_value('alpha',alpha)
set_param_value('nu',nu)
set_param_value('delta',delta)
set_param_value('zeta0',zeta0)
set_param_value('beta',beta)
set_param_value('rho',rho)
set_param_value('zeta1',zeta1)
set_param_value('zeta2',zeta2)
set_param_value('k_ss',k_ss)
set_param_value('b_0',b_0)
set_param_value('b_K',b_K)
set_param_value('b_z',b_z)
set_param_value('sig_e1',sig_e1)
set_param_value('sig_e2',sig_e2)
%
% INITVAL instructions
%
options_.initval_file = false;
oo_.steady_state(2) = M_.params(1);
oo_.steady_state(1) = oo_.steady_state(6)-M_.params(4)*oo_.steady_state(2);
oo_.steady_state(4) = M_.params(1);
oo_.steady_state(5) = M_.params(8)*oo_.steady_state(4)^(M_.params(8)-1);
oo_.steady_state(6) = (1-M_.params(8))*oo_.steady_state(4)^M_.params(8);
oo_.steady_state(3) = 1;
if M_.exo_nbr > 0
	oo_.exo_simul = ones(M_.maximum_lag,1)*oo_.exo_steady_state';
end
if M_.exo_det_nbr > 0
	oo_.exo_det_simul = ones(M_.maximum_lag,1)*oo_.exo_det_steady_state';
end
%
% SHOCKS instructions
%
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = (M_.params(13))^2;
M_.Sigma_e(2, 2) = (M_.params(14))^2;
options_.irf = 0;
options_.nocorr = true;
options_.nomoments = true;
options_.order = 2;
var_list_ = {'k';'Ka'};
[info, oo_, options_, M_] = stoch_simul(M_, options_, oo_, var_list_);
save('model_agg_uncertainty_results.mat', 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save('model_agg_uncertainty_results.mat', 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save('model_agg_uncertainty_results.mat', 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save('model_agg_uncertainty_results.mat', 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save('model_agg_uncertainty_results.mat', 'estimation_info', '-append');
end
if exist('dataset_info', 'var') == 1
  save('model_agg_uncertainty_results.mat', 'dataset_info', '-append');
end
if exist('oo_recursive_', 'var') == 1
  save('model_agg_uncertainty_results.mat', 'oo_recursive_', '-append');
end


disp(['Total computing time : ' dynsec2hms(toc(tic0)) ]);
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
diary off
