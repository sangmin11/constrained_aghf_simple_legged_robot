% one massless leg model
function [get_f, get_F, get_H, get_L, N, M] = multileg_multiple_fixed_switch_model(flags, params, param_model, extra_params)
global get_EL get_G get_G_using
persistent get_EL1 get_G1 get_f1 get_F1 get_H1 get_L1 N1 M1 get_G_using1
if( isempty(get_G1) )
    
    %-------------------this portion should be edited------------------
    
    m = param_model(1); % mass
    J = param_model(2); % inertia
    g = param_model(3); % gravity acc
    k_leg = param_model(4); % number of legs
    
    % symbolic states
    N = 6+4*k_leg; % # of states
    M = 4*k_leg; % # of inputs
    D = 6; % # of duals
    k = params(1); % k is lambda
    alpha = params(2);
    k_dynamics = params(3);
    T_switch = extra_params('T_switch');
    Rmax_sq = extra_params('Rmax_sq');
    mu = extra_params('mu');
    
    X = [];
    Xd = [];
    for i=1:N
        eval(['syms x' num2str(i) ' real']);
        eval(['syms x' num2str(i) 'd' ' real']);
        X = [X; eval(['x' num2str(i)])];
        Xd = [Xd; eval(['x' num2str(i) 'd'])];
    end

    W = [];
    for i = 1:D
        eval(['syms w' num2str(i) ' real']);
        W = [W; eval(['w' num2str(i)])];
    end
    
    Force = [];
    Pos = [];
    
    for i = 1:k_leg
        Force = [Force [X(7+4*(i-1));X(7+4*(i-1)+1)] ];
        Pos = [Pos [X(7+4*(i-1)+2);X(7+4*(i-1)+3)] ];
    end
    
    syms tt real
    
    in = [];
    if (flags(1))
        syms k real;
        in = [in; k];
    end
    
    if (flags(2))
        syms alpha real;
        in = [in; alpha];
    end

    if (flags(3))
        syms k_dynamics real;
        in = [in; k_dynamics];
    end
    
    % HF formulation
    
    % all control directions
    F = eye(N);
    F = sym(F);
    
    % totol torque on CoM
    Tor_sum = 0;
    for i=1:k_leg
        Tor_sum = Tor_sum + det([Pos(:,i)-[x1;x2] Force(:,i)]);
    end
    % drift term
    f = [ [x4 x5 x6]' ;  sum(Force,2)/m - [0;g] ; Tor_sum/J; zeros(M,1)];
    
    % cutomized step function
    step = @(tt,tt1) heaviside(tt-tt1);
    
    % approximation of heaviside function
    k_step = 100;
    %k_step = 50;
    heaviside_appr = @(tt) 1/(1+exp(-2*k_step*tt));

    % approximation of heaviside function derivative
    a_step = 0.04;
    %a_step = 0.08;
    dirac_appr = @(tt) 1/((2*pi)^0.5*a_step)*exp(-tt^2/(2*a_step^2));
    
    % activation functions for each leg
    Activation = [];
    for j = 1:k_leg
        activation = 1-step(tt, T_switch{j}(1) );
        for i = 2:2:length(T_switch{j})-1
            temp = (step(tt,T_switch{j}(i))-step(tt, T_switch{j}(i+1) ));
            activation = activation + temp;
        end
        activation = activation + step(tt, T_switch{j}(end) );
        Activation = [Activation; activation];
    end
    
    % lambda for each state
    K = diag( [k_dynamics k_dynamics k_dynamics k_dynamics k_dynamics k_dynamics repmat([1 1 1 1],1,k_leg)] );
    for i=1:k_leg
        K(9+4*(i-1),9+4*(i-1)) = 1 + k*Activation(i);
        K(9+4*(i-1)+1,9+4*(i-1)+1) = 1 + k*Activation(i);
    end
    
    % state constraints
    DpyDpx = [];
    Tan = []; % normal direction at contact point
    Nor = []; % tangent direction at contact point
    FNor = []; % contact force normal
    FTan = []; % contact force tangent
    R_sq = []; % contact point to CoM distance square
    dterrain_func = @(x) diff(terrain_func(x),x);
    for i = 1:k_leg
        dpydpx = dterrain_func(Pos(1,i));
        DpyDpx = [DpyDpx dpydpx];
        tan = [1 dpydpx]';
        Tan = [Tan tan];
        nor = [-dpydpx 1]';
        Nor = [Nor nor];
        FNor = [FNor Force(:,i)'*nor];
        FTan = [FTan Force(:,i)'*tan];
        R_sq = [R_sq (Pos(1,i)-x1)^2+(Pos(2,i)-x2)^2];
    end
    
    % min CoM height
    ymin = 0.3;
    
    %initialize the state constraints cost, b
    b = k*(x2-ymin)^2*(1-step(x2,ymin));

    % predefine max contanct force
    ForceMax = 4*m*g/k_leg;
    
    % sum up all state constraint for all legs
    for i = 1:k_leg
        
        pknee = find_knee(X(1:2),Pos(:,i),0.65,0.65,1);
        FnormSqr = Force(1,i)^2+Force(2,i)^2;
        
        % the constraints include:
        % zero force in flight phase, contact point to CoM distance,
        % non-negative contact normal force when in stance phase
        % friction cone
        % foot on terrain in stance phase
        % foot above terrain in flght phase
        % max contact force
        % knee height constraint
        b = b + ( k*Force(1,i)^2+k*Force(2,i)^2 )*(1-Activation(i)) + 10*k*(R_sq(i)-Rmax_sq)^2*step(R_sq(i),Rmax_sq) + ...
         k*FNor(i)^2*(1-step(FNor(i),0))*Activation(i) + ...
         k*(FTan(i)+mu*FNor(i))^2*(1-step(FTan(i)+mu*FNor(i),0))*Activation(i) + ...
         k*(FTan(i)-mu*FNor(i))^2*step(FTan(i)-mu*FNor(i),0)*Activation(i) + ...
         5*k*(terrain_func(Pos(1,i))-Pos(2,i))^2*Activation(i) + ...
         k*(terrain_func(Pos(1,i))-Pos(2,i))^2*step(terrain_func(Pos(1,i))-Pos(2,i),0)*(1-Activation(i)) +...
         k*FnormSqr*step(FnormSqr,ForceMax^2) +...
         5*k*(pknee(2)-0.00)^2*(1-step(pknee(2)-0.00,0));
    end
    
    % calculate the overall Lagrangian L
    H = (F.')^(-1)*K*F^(-1);
    G = H;
    L = (Xd - f*alpha).' * G * (Xd - f*alpha) + b + 2*(Xd(1:D) - f(1:D)*alpha).' * G(1:D,1:D) * W;
    
    %---------------------------------------------------------------------%
    
    %---------the following part are generic, don't need to edit----------%
    
    % get Euler Lagrange euqation symbolically
    % get dL/dx and dL/ddx
    pLx = sym('pLx', [N+D 1],'real');
    pLxd = sym('pLxd', [N+D 1],'real');
    
    for i=1:N
        pLx(i) =  diff(L,X(i));
        pLxd(i) = diff(L,Xd(i));
    end

    for i=1:D
        pLx(N+i) = -diff(L,W(i));
        pLxd(N+i) = 0;
    end
    
    % replace the heaviside funtion with customized heaviside function
    for i=1:N+D    
        temp_cell = arrayfun(@char, pLx(i), 'uniform', 0);
        str_temp = strrep(temp_cell{1},'heaviside','heaviside_appr');
        str_temp = strrep(str_temp,'dirac','dirac_appr');
        pLx(i) =  eval( evalin(symengine,str_temp) );
        
        temp_cell = arrayfun(@char, pLxd(i), 'uniform', 0);
        str_temp = strrep(temp_cell{1},'heaviside','heaviside_appr');
        str_temp = strrep(str_temp,'dirac','dirac_appr');
        pLxd(i) =  eval( evalin(symengine,str_temp) );
    end
    
    % generate functions
    in2 = [X;in];
    get_EL = matlabFunction([pLx, pLxd],'vars',[tt,X',W',Xd',in']);
    get_G = matlabFunction(G,'vars', {tt,X,in});
    get_G_using = matlabFunction(G,'vars', [tt,X',in']);
    get_f = matlabFunction(f,'vars',{tt,X});
    get_F = matlabFunction(F,'vars',{tt,X});
    get_H = matlabFunction(H,'vars',{tt,in2});
    get_L = matlabFunction(L,'vars',{tt,X,W,Xd,in});
    
    
    get_EL1 = get_EL;
    get_G1 = get_G;
    get_G_using1 = get_G_using;    
    get_f1 = get_f;
    get_F1 = get_F;
    get_H1 = get_H;
    get_L1 = get_L;
    N1 = N;
    M1 = M;
    
else
    
    get_EL = get_EL1;
    get_G = get_G1;
    get_G_using = get_G_using1;    
    get_f = get_f1;
    get_F = get_F1;
    get_H = get_H1;
    get_L = get_L1;
    N = N1;
    M = M1;
    
end




end