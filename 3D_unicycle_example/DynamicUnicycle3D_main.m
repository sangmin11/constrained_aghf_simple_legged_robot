%% 3D dynamic unicycle main (quaternion)
clear;
clc;
restoredefaultpath;
addpath('func');
addpath('HF_models');
addpath('generated_files');
addpath('plot');


param_model = [];


extra_params = containers.Map; % a map for extra parameters

%-------------------heat flow parameters--------------------------------%

smax = 0.002; % max value of s
tgrids = 100; % number of grids for t
sgrids = 100; % number of grids for s
intgrids = 5000; % number of grids for integration
k_sym_flag = 1; % flag for making k a variable
alpha_sym_flag = 1; % flag for making alpha a variable
% if the flags are set to 0, when changing k or alpha, one has to "clear
% all" and rerun the script.
k =  50000;
alpha = 1;
T = 1; % total time

params = [k alpha];
flags = [k_sym_flag, alpha_sym_flag];

ax = [0 0 1];
ang = pi/6;
%ang = 0;
q0 = axang2quat([ax ang])';
X0 = [0; 0; 0;   q0;   0; 0; 0;  0; 0; 0];
ax = [0 0 1];
ang = pi/6+2*pi;
%ang = 0;
qf = axang2quat([ax ang])';
Xf = [0; 1; 1;   q0;   0; 0; 0;  0; 0; 0];

Xinit = [X0 Xf];

free_ind = [0 0 0   0 0 0 0   0 0 0   0 0 0;  %first column: initial value
            0 0 0   0 0 0 0   0 0 0   0 0 0]'; %second column: final value

admissible_control = [8 12 13]; %[vxd wyd wzd]
extra_params('admissible_control') = admissible_control;
extra_params('Quat_ind0') = 4;
extra_params('w_ind0') = 11;

%-------------------solve heat flow pde without drift------------------------------% 
extra_params('AbsTol') = 1;
extra_params('RelTol') = 1;


[p, dp, xfull, xdrift, bufull, budrift, sol, Xs, Xf, Xint, U, s, t, tint, cost, cost_all, cost_u] = solve_HF(@DynamicUnicycle3D_quat_model,param_model,smax,tgrids,intgrids,sgrids,Xinit,T,params,flags,free_ind,extra_params);
%-------------------------------------------------------------------------%

%% plot
run('DynamicUnicycle3D_quat_plot.m');
