%% multi legs with fixed switch time, mutiple steps, main
clear;
clc;
addpath('func');
addpath('HF_models');
addpath('generated_files');
addpath('plot');


m = 2; %mass
J = 1; %inertia
g = 9.8; % gravity
k_leg = 2; % # of legs
param_model = [m J g k_leg];
Rmax_sq = 1; % max radius the feet can reach
mu = 1; % friction coefficient

extra_params = containers.Map; % a map for extra parameters
extra_params('Rmax_sq') = Rmax_sq;
extra_params('mu') = mu;


%-------------------heat flow parameters--------------------------------%
smax = 0.0005; % max value of s
tgrids = 100; % number of grids for t
sgrids = 100; % number of grids for s
intgrids = 100000; % number of grids for integration
k_sym_flag = 1; % flag for making k a variable
alpha_sym_flag = 1; % flag for making alpha a variable
% if the flags are set to 0, when changing k or alpha, one has to "clear
% all" and rerun the script.
k =  500000;
alpha = 1; % scaler of the drift term
T = 2; % total time

%---------------------predefined gait---------------------------
steps = 3;
dT = T/(2*steps+1); %for same stance and flight phase duration
r_stance = 0.5; %stance phase ratio, for different stance and flight phase duration
dT_1step = T/(steps+r_stance);
dT_stance = dT_1step * r_stance;
dT_flight = dT_1step * (1-r_stance);

t_sym_flag = ones(1,2*steps);
t_switch1 = [];
for i = 1:2*steps
    if(mod(i,2)==1)
        t_switch1 = [t_switch1 (i-1)/2*dT_1step+dT_stance];
    else
        t_switch1 = [t_switch1 i/2*dT_1step];
    end
end
% leg 2 has a timing delay 0.05s than leg 1
t_switch2 = t_switch1-0.05;  %need to rebuild model if timing is changed


T_switch = {t_switch1 t_switch2};
extra_params('T_switch') = T_switch;
params = [k alpha];
flags = [k_sym_flag, alpha_sym_flag];

%boudanry values
% state X = [px py theta \dot(px) \dot(py) \dot(theta) F1x F1y p1x p1y F2x F2y p2x p2y]
X0 = [0; 0.75; 0;   0.5; 0; 0;    1;1;0;0;   1;1;0;0]; % py = 0, flat terrian, 2 legs
Xf = [1.5; 0.75; 0;   0.5; 0; 0;    1;1;1.5;0;   1;1;1.5;0];


Xinit = [X0 Xf];
% flags for free boudanry values
free_ind = [0 0 0   0 0 0   1 1 0 0  1 1 0 0;  %first column: initial value
            0 0 0   0 0 0   1 1 0 0  1 1 0 0]'; %second column: final value



%-------------------solve heat flow pde without drift------------------------------% 
extra_params('AbsTol') = 1;
extra_params('RelTol') = 1;
[p, dp, xfull, xdrift, bufull, budrift, sol, Xs, Xf, Xint, U, s, t, tint, cost, cost_all, cost_u] = solve_HF(@multileg_multiple_fixed_switch_model,param_model,smax,tgrids,intgrids,sgrids,Xinit,T,params,flags,free_ind,extra_params);
%-------------------------------------------------------------------------%
%% plot and animation
w_robot = 0.4;
h_robot = 0.3;
l1=0.6; l2=0.6;
force_scale = 0.008;
%force_scale = 0.002;
run('multileg_multisteps_plot.m');