% one massless leg model
function [get_f, get_F, get_H, get_L, N, M] = SCA_model(flags, params, param_model, extra_params)
%-------------------this portion should be edited------------------

% symbolic states
N = 7; % # of states
M = length(extra_params('admissible_control')); % # of inputs
k = params(1);
alpha = params(2);

force = sym('force',[3 1],'real');

in = [];
if (flags(1))
    syms k real;
    in = [in; k];
end

if (flags(2))
    syms alpha real;
    in = [in; alpha];
end

in = [in; force];

syms x y z xd yd zd q0 q1 q2 q3 q0d q1d q2d q3d real
X = [x y z q0 q1 q2 q3]';
Xd = [xd yd zd q0d q1d q2d q3d]';


% HF formulation
f = sym(zeros(N,1));

Fp_xyz = [q0^2 + q1^2 - q2^2 - q3^2,          2*q1*q2 - 2*q0*q3,          2*q0*q2 + 2*q1*q3;
    2*q0*q3 + 2*q1*q2,  q0^2 - q1^2 + q2^2 - q3^2,          2*q2*q3 - 2*q0*q1;
    2*q1*q3 - 2*q0*q2,          2*q0*q1 + 2*q2*q3,  q0^2 - q1^2 - q2^2 + q3^2];

Fq_123 = [-q1/2, -q2/2, -q3/2;
    q0/2, -q3/2,  q2/2;
    q3/2,  q0/2, -q1/2;
    -q2/2,  q1/2,  q0/2];

f(1:3) = Fp_xyz * [1 0 0]';

%w_star = [0 2*pi 0]';
w_star = [0 0 pi/2]';
%w_star = [0 0 0]';
f(4:7) =  Fq_123 * w_star;

fq0 = [q0/2 q1/2 q2/2 q3/2]';

F_original = blkdiag(eye(3), [fq0 Fq_123]);

admissible_control = extra_params('admissible_control');

Fc = F_original;
Fc(:,admissible_control) = [];
F_free = F_original(:,admissible_control);

F = sym( [Fc F_free] ); %[vx vy vz w0 wx wy wz]

K = diag( [k*ones(1,N-M) 100*ones(1,M)] );

syms tt real

step = @(tt,tt1) heaviside(tt-tt1);

k_step = 100;
heaviside_appr = @(tt) 1/(1+exp(-2*k_step*tt));
a_step = 0.01;
dirac_appr = @(tt) exp(-tt^2/(2*a_step^2));

fixed_points = extra_params('fixed_points');
fixed_length = extra_params('fixed_length');

[~,mp] = size(fixed_points);

b = 0;
k_marker = 200*k;
for i = 1:mp
    b = b + k_marker*(x-fixed_points(1,i))^2*dirac_appr(tt-fixed_length(i));
    b = b + k_marker*(y-fixed_points(2,i))^2*dirac_appr(tt-fixed_length(i));
    b = b + k_marker*(z-fixed_points(3,i))^2*dirac_appr(tt-fixed_length(i));
end

density = 100;
%g = 9.8;
g = 0; % no gravity
b = b + 1*density*g*z; % potential of gravity

pd = [xd yd zd]';
b = b - pd'*force;

H = simplify( (F.')^(-1)*K*F^(-1) );
G = H;
L = simplify( (Xd - f*alpha).' * G * (Xd - f*alpha) + b );

%---------------------------------------------------------------------%

%---------the following part are generic, don't need to edit----------%

pLx = sym('pLx', [N 1],'real');
pLxd = sym('pLxd', [N 1],'real');
for i=1:N
    pLx(i) =  diff(L,X(i));
    pLxd(i) = diff(L,Xd(i));
end

in2 = [X;in];
%get_EL = matlabFunction([pLx, pLxd],'vars',[tt,X',Xd',in'],'File','generated_files\get_EL');
%get_G = matlabFunction(G,'vars', {tt,X,in},'File','generated_files\get_G');
get_f = matlabFunction(f,'vars',{tt,X},'File','generated_files\get_f');
get_F = matlabFunction(F,'vars',{tt,X});
get_H = matlabFunction(H,'vars',{tt,in2});
get_L = matlabFunction(L,'vars',{tt,X,Xd,in});


end