% one massless leg model
function [get_f, get_F, get_H, get_L, N, M] = DynamicUnicycle3D_quat_model(flags, params, param_model, extra_params)
global get_EL get_G
persistent get_EL1 get_G1 get_f1 get_F1 get_H1 get_L1 N1 M1
if( isempty(get_G1) )
    
    %-------------------this portion should be edited------------------
    
    % symbolic states
    N = 7+6; % # of states
    M = 3; % # of inputs
    k = params(1);
    alpha = params(2);
    
    in = [];
    if (flags(1))
        syms k real;
        in = [in; k];
    end
    
    if (flags(2))
        syms alpha real;
        in = [in; alpha];
    end
    
    syms x y z vx vy vz q0 q1 q2 q3 wx wy wz real
    syms xd yd zd vxd vyd vzd q0d q1d q2d q3d wxd wyd wzd real
    X = [x y z q0 q1 q2 q3 vx vy vz wx wy wz]';
    Xd = [xd yd zd q0d q1d q2d q3d vxd vyd vzd wxd wyd wzd]';
    
    v = [vx vy vz]';
    w = [wx wy wz]';
    
    
    % HF formulation
    
    Fp_xyz = [q0^2 + q1^2 - q2^2 - q3^2,          2*q1*q2 - 2*q0*q3,          2*q0*q2 + 2*q1*q3;
                  2*q0*q3 + 2*q1*q2,  q0^2 - q1^2 + q2^2 - q3^2,          2*q2*q3 - 2*q0*q1;
                  2*q1*q3 - 2*q0*q2,          2*q0*q1 + 2*q2*q3,  q0^2 - q1^2 - q2^2 + q3^2];
              
    
    Fq = [-q1/2, -q2/2, -q3/2;
           q0/2, -q3/2,  q2/2;
           q3/2,  q0/2, -q1/2;
          -q2/2,  q1/2,  q0/2];
    
    f = [ Fp_xyz*v; Fq*w; zeros(6,1)];
    
    F_original = eye(N); % [vx vy vz q0d q1d q2d q3d vxd vyd vzd wxd wyd wzd]
    
    admissible_control = extra_params('admissible_control');
    
    Fc = F_original;
    Fc(:,admissible_control) = [];
    F_free = F_original(:,admissible_control);
    
    F = sym( [Fc F_free] );
    
    K = diag( [k*ones(1,N-M) ones(1,M)] );
    
    b = 0;
    H = simplify( (F.')^(-1)*K*F^(-1) );
    G = H;
    %L =  (Xd - f).' * G * (Xd - f) ;
    L = simplify( (Xd - f*alpha).' * G * (Xd - f*alpha) + b );
    
    syms tt real
    
    %---------------------------------------------------------------------%
    
    %---------the following part are generic, don't need to edit----------%
    
    pLx = sym('pLx', [N 1],'real');
    pLxd = sym('pLxd', [N 1],'real');
    for i=1:N
        pLx(i) =  diff(L,X(i));
        pLxd(i) = diff(L,Xd(i));
    end
    
    in2 = [X;in];
    get_EL = matlabFunction([pLx, pLxd],'vars',[tt,X',Xd',in']);
    get_G = matlabFunction(G,'vars', {tt,X,in});
    get_f = matlabFunction(f,'vars',{tt,X});
    get_F = matlabFunction(F,'vars',{tt,X});
    get_H = matlabFunction(H,'vars',{tt,in2});
    get_L = matlabFunction(L,'vars',{tt,X,Xd,in});
    
    
    get_EL1 = get_EL;
    get_G1 = get_G;
    get_f1 = get_f;
    get_F1 = get_F;
    get_H1 = get_H;
    get_L1 = get_L;
    N1 = N;
    M1 = M;
    
else
    
    get_EL = get_EL1;
    get_G = get_G1;
    get_f = get_f1;
    get_F = get_F1;
    get_H = get_H1;
    get_L = get_L1;
    N = N1;
    M = M1;
    
end




end