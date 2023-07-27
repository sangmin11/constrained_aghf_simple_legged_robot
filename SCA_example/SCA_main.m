%% soft continuum arm main
%clear;
%clc;
restoredefaultpath;
addpath('func');
addpath('HF_models');
addpath('generated_files');
addpath('plot');


param_model = [];


extra_params = containers.Map; % a map for extra parameters

%-------------------heat flow parameters--------------------------------%

smax = 1e10;
tgrids = 200; % number of grids for t
sgrids = 100; % number of grids for s
intgrids = 2000; % number of grids for integration
k_sym_flag = 1; % flag for making k a variable
alpha_sym_flag = 1; % flag for making alpha a variable
% if the flags are set to 0, when changing k or alpha, one has to "clear
% all" and rerun the script.
k =  100000;
alpha = 1;
T = 1; % total time

% initial guess is the shape defined by omega*
w_star = [0 0 pi/2]';
t_guess = linspace(0,T,50);
Xinit = zeros(7,length(t_guess));
speed = pi/2;
w = [0 0 speed];
r = (2*pi)/(speed)/(2*pi);
for i = 1:length(t_guess)
    th =  speed * t_guess(i);
    Xinit(4:7,i) = axang2quat([w th])';
    Xinit(1:3,i) = [r*sin(th) r*(1-cos(th)) 0]';
end



% for generating tip position
free_ind = [0 0 0 0 0 0 0;  %first column: initial value
            1 1 1 1 1 1 1]'; %second column: final value

        
% change Xinit(1:3,end) to desired tip position and fixe it:
%Xinit(1:3,end) = Xsol(1:3,end);
% for fixing tip position
%free_ind = [0 0 0 0 0 0 0;  %first column: initial value
%            0 0 0 1 1 1 1]'; %second column: final value

% specify force (generating tip position)
force = [-200 400 600];
% zero force (for reconstructing pose with tip position)
%force = [0 0 0];

params = [k alpha force];
flags = [k_sym_flag, alpha_sym_flag, [1 1 1]];

fixed_points = [];
fixed_length = [];
extra_params('fixed_points') = fixed_points;
extra_params('fixed_length') = fixed_length;

admissible_control = [6 7]; %[wx wy wz]
extra_params('admissible_control') = admissible_control;
extra_params('Quat_ind0') = 4;

%-------------------solve heat flow pde without drift------------------------------%
extra_params('AbsTol') = 1;

extra_params('RelTol') = 1;

[Xsol, dXsol, xfull, xdrift, bufull, budrift, sol, Xs, Xf, Xint, U, s, t, tint, cost, cost_all, cost_u] = solve_HF(@SCA_model,param_model,smax,tgrids,intgrids,sgrids,Xinit,T,params,flags,free_ind,extra_params);
%-------------------------------------------------------------------------%
%% plot
run('SCA_plot.m');