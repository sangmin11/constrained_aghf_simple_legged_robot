%% diver midair
clear;
clc;
restoredefaultpath;
addpath('func');
addpath('HF_models');
mkdir('generated_files');
addpath('generated_files');
addpath('plot');



m = 1; m1 = 0.2; m2 = 1;
l = 0.8; l1 = 1; l2 = 1.5;
I = 1/12 * m*l^2; I1 = 1/12 * m1*l1^2; I2 = 1/12 * m2*l2^2;
Imax = I+I1+I2+m1*(l/2+l1/2)^2 + m2*(l/2+l2/2)^2;
Imin = I1+I2+I;
thd_avg = 2*pi;
th2d_max = Imax*thd_avg/(I2+m2*(l/2+l2/2)^2);
th2d_min = Imin*thd_avg/(I2+m2*(l/2+l2/2)^2);
beta = 0.8;
th2d0 = th2d_min*(1-beta)+th2d_max*beta;
%th2d0 = 1.5*th2d_max;

g = 9.8;
param_model = [m m1 m2 I I1 I2 l l1 l2 g];
extra_params = containers.Map; % a map for extra parameters

%-------------------heat flow parameters--------------------------------%

smax = 0.05; % without G_inverse
tgrids = 100; % number of grids for t
sgrids = 50; % number of grids for s
intgrids = 10000; % number of grids for integration using extracted control
k_sym_flag = 1; % flag for making k a variable
alpha_sym_flag = 1; % flag for making alpha a variable
% if the flags are set to 0, when changing k or alpha, one has to "clear
% all" and rerun the script.
flags = [k_sym_flag, alpha_sym_flag];
k = 1000;
alpha = 1;
params = [k alpha];
T = 1; % total time

X0 = [0;0;0;   0;0;th2d0]; % initial value
Xf = [2*pi;0;0;  0;0;0]; % final value
Xinit = [X0 Xf];
% set free/fixed bc, 1 for free, 0 for fixed
%case 1
%free_ind = [0 0 0 0 0 0;  %first column: initial value
%            0 0 0 1 1 1]'; %second column: final value
%case 2
free_ind = [0 0 0 0 0 1;  %first column: initial value
            0 0 0 0 1 0]'; %second column: final value
%free_ind = [0 0 0 0 0 0;  %first column: initial value
%            0 1 0 1 1 1]'; %second column: final value

%-------------------solve heat flow pde without drift------------------------------%
extra_params('AbsTol') = 1;
extra_params('RelTol') = 0.1;
[p, dp, xfull, xdrift, bufull, budrift, sol, Xs, Xf, Xint, U, s, t, tint, cost, cost_all, cost_u, T_solve] = solve_HF(@diver_midair_model,param_model,smax,tgrids,intgrids,sgrids,Xinit,T,params,flags,free_ind,extra_params);
%-------------------------------------------------------------------------%

%%
run('diver_midair_plot.m');