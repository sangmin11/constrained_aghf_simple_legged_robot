% 2-link diver model
function [get_f, get_F, get_H, get_L, N, M] = diver_midair_model(flags, params,  param_model, extra_params)
global get_EL get_G
persistent get_EL1 get_G1 get_f1 get_F1 get_H1 get_L1 N1 M1
if( isempty(get_G1) )
    
    %-------------------this portion should be edited------------------%
    
    m = param_model(1); m1 = param_model(2); m2 = param_model(3);
    I = param_model(4); I1 = param_model(5); I2 = param_model(6);
    l = param_model(7); l1 = param_model(8); l2 = param_model(9);
    g = param_model(10);
    
    % symbolic states
    N = 6; % # of states
    M = 2; % # of inputs
    X = sym('X', [N 1],'real');
    Xd = sym('Xd', [N 1],'real');
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
    
    % HF formulation
    th=X(1); th1=X(2); th2=X(3);
    thd=X(4); th1d=X(5); th2d=X(6);
    
    %q = [x y th th1 th2].';
    qd = [thd th1d th2d].';
    
    D_dynamics = [ (4*I1*m + 4*I1*m1 + 4*I2*m + 4*I1*m2 + 4*I2*m1 + 4*I2*m2 + 4*I*m + 4*I*m1 + 4*I*m2 + l^2*m*m1 + l^2*m*m2 + l1^2*m*m1 + 4*l^2*m1*m2 + l1^2*m1*m2 + l2^2*m*m2 + l2^2*m1*m2 + 2*l1*l2*m1*m2*cos(th1 - th2) + 2*l*l1*m*m1*cos(th1) + 4*l*l1*m1*m2*cos(th1) + 2*l*l2*m*m2*cos(th2) + 4*l*l2*m1*m2*cos(th2))/(4*(m + m1 + m2)), (4*I1*m + 4*I1*m1 + 4*I1*m2 + l1^2*m*m1 + l1^2*m1*m2 + l1*l2*m1*m2*cos(th1 - th2) + l*l1*m*m1*cos(th1) + 2*l*l1*m1*m2*cos(th1))/(4*(m + m1 + m2)), (4*I2*m + 4*I2*m1 + 4*I2*m2 + l2^2*m*m2 + l2^2*m1*m2 + l1*l2*m1*m2*cos(th1 - th2) + l*l2*m*m2*cos(th2) + 2*l*l2*m1*m2*cos(th2))/(4*(m + m1 + m2));
        (4*I1*m + 4*I1*m1 + 4*I1*m2 + l1^2*m*m1 + l1^2*m1*m2 + l1*l2*m1*m2*cos(th1 - th2) + l*l1*m*m1*cos(th1) + 2*l*l1*m1*m2*cos(th1))/(4*(m + m1 + m2)),                                                                           (4*I1*m + 4*I1*m1 + 4*I1*m2 + l1^2*m*m1 + l1^2*m1*m2)/(4*(m + m1 + m2)),                                                                                                    (l1*l2*m1*m2*cos(th1 - th2))/(4*(m + m1 + m2));
        (4*I2*m + 4*I2*m1 + 4*I2*m2 + l2^2*m*m2 + l2^2*m1*m2 + l1*l2*m1*m2*cos(th1 - th2) + l*l2*m*m2*cos(th2) + 2*l*l2*m1*m2*cos(th2))/(4*(m + m1 + m2)),                                                                                                    (l1*l2*m1*m2*cos(th1 - th2))/(4*(m + m1 + m2)),                                                                           (4*I2*m + 4*I2*m1 + 4*I2*m2 + l2^2*m*m2 + l2^2*m1*m2)/(4*(m + m1 + m2))];
    
    
    C_dynamics = [ -(l*l1*m*m1*th1d*sin(th1) + 2*l*l1*m1*m2*th1d*sin(th1) + l*l2*m*m2*th2d*sin(th2) + 2*l*l2*m1*m2*th2d*sin(th2) + l1*l2*m1*m2*th1d*sin(th1 - th2) - l1*l2*m1*m2*th2d*sin(th1 - th2))/(4*(m + m1 + m2)), -(l1*m1*(th1d + thd)*(l*m*sin(th1) + 2*l*m2*sin(th1) + l2*m2*sin(th1 - th2)))/(4*(m + m1 + m2)), -(l2*m2*(th2d + thd)*(l*m*sin(th2) + 2*l*m1*sin(th2) - l1*m1*sin(th1 - th2)))/(4*(m + m1 + m2));
        (l1*m1*(l*m*thd*sin(th1) + 2*l*m2*thd*sin(th1) + l2*m2*th2d*sin(th1 - th2) + l2*m2*thd*sin(th1 - th2)))/(4*(m + m1 + m2)),                                                                                               0,                                     (l1*l2*m1*m2*sin(th1 - th2)*(th2d + thd))/(4*(m + m1 + m2));
        (l2*m2*(l*m*thd*sin(th2) + 2*l*m1*thd*sin(th2) - l1*m1*th1d*sin(th1 - th2) - l1*m1*thd*sin(th1 - th2)))/(4*(m + m1 + m2)),                                    -(l1*l2*m1*m2*sin(th1 - th2)*(th1d + thd))/(4*(m + m1 + m2)),                                                                                               0];
    
    
    G_dynamics = [0;0;0];
    
    Dinv = simplify(D_dynamics^(-1));
    
    f = [thd; th1d; th2d; -Dinv*(C_dynamics*qd+G_dynamics)];
    F = [eye(3) zeros(3,3); zeros(3,3) Dinv];
    F_inv = [eye(3) zeros(3,3); zeros(3,3) D_dynamics];
    

    
    %kb = 0.0001; % barrier function gain
    %pb = 1; % order of barrier function
    %b = kb/(th2-pi/2)^pb + kb/(th2+pi/2)^pb;
    
    K = diag([k k k k 1 1]);
    %H = simplify( (F.')^(-1)*K*F^(-1) );
    H = F_inv*K*F_inv; % since D is symmetric
    %G = (b+1)*H;
    G = H;
    %L =  (Xd - f).' * G * (Xd - f) ;
    k_step = 50;
    step = @(tt,tt1) 1/(1+exp(-2*k_step*(tt-tt1)));
    
    th2_max = 1.9;
    b = 0;
    %b = 10*( k*(th2-th2_max)^2*step(th2,th2_max) + k*(th2+th2_max)^2*(1-step(th2,-th2_max)) );
    
    %b = 10*k*(th+th1)^2;
    %b = 0;
    %L =  (Xd - f*alpha).' * G * (Xd - f*alpha);
    L =  (Xd - f*alpha).' * G * (Xd - f*alpha) + b;
    
    
    
    %---------------------------------------------------------------------%
    
    %---------the following part are generic, don't need to edit----------%
    
    syms tt real
    
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

rehash;


end