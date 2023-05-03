% one massless leg model
function [get_f, get_F, get_H, get_L, N, M] = multileg_multiple_fixed_switch_model(flags, params, param_model, extra_params)
global get_EL get_G
persistent get_EL1 get_G1 get_f1 get_F1 get_H1 get_L1 N1 M1
if( isempty(get_G1) )
    
    %-------------------this portion should be edited------------------
    
    m = param_model(1);
    J = param_model(2);
    g = param_model(3);
    k_leg = param_model(4);
    
    % symbolic states
    N = 6+4*k_leg; % # of states
    M = 4*k_leg; % # of inputs
    k = params(1);
    alpha = params(2);
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
    
    % HF formulation
    
    F = eye(N);
    F = sym(F);
    
    
    Tor_sum = 0;
    for i=1:k_leg
        Tor_sum = Tor_sum + det([Pos(:,i)-[x1;x2] Force(:,i)]);
    end
    f = [ [x4 x5 x6]' ;  sum(Force,2)/m - [0;g] ; Tor_sum/J; zeros(M,1)];
    
    step = @(tt,tt1) heaviside(tt-tt1);
 
    k_step = 100;
    %k_step = 50;
    heaviside_appr = @(tt) 1/(1+exp(-2*k_step*tt));
    a_step = 0.04;
    %a_step = 0.08;
    dirac_appr = @(tt) 1/((2*pi)^0.5*a_step)*exp(-tt^2/(2*a_step^2));
    
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
    
    K = diag( [k k k k k k repmat([1 1 1 1],1,k_leg)] );
    for i=1:k_leg
        K(9+4*(i-1),9+4*(i-1)) = 1 + k*Activation(i);
        K(9+4*(i-1)+1,9+4*(i-1)+1) = 1 + k*Activation(i);
    end
    
    
    DpyDpx = [];
    Tan = [];
    Nor = [];
    FNor = [];
    FTan = [];
    R_sq = [];
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
    
    ymin = 0.3;
    
    
    b = k*(x2-ymin)^2*(1-step(x2,ymin));
    ForceMax = 4*m*g/k_leg;
    
    for i = 1:k_leg
        
        pknee = find_knee(X(1:2),Pos(:,i),0.65,0.65,1);
        FnormSqr = Force(1,i)^2+Force(2,i)^2;
        
        b = b + ( k*Force(1,i)^2+k*Force(2,i)^2 )*(1-Activation(i)) + 10*k*(R_sq(i)-Rmax_sq)^2*step(R_sq(i),Rmax_sq) + ...
         k*FNor(i)^2*(1-step(FNor(i),0))*Activation(i) + ...
         k*(FTan(i)+mu*FNor(i))^2*(1-step(FTan(i)+mu*FNor(i),0))*Activation(i) + ...
         k*(FTan(i)-mu*FNor(i))^2*step(FTan(i)-mu*FNor(i),0)*Activation(i) + ...
         5*k*(terrain_func(Pos(1,i))-Pos(2,i))^2*Activation(i) + ...
         k*(terrain_func(Pos(1,i))-Pos(2,i))^2*step(terrain_func(Pos(1,i))-Pos(2,i),0)*(1-Activation(i)) +...
         k*FnormSqr*step(FnormSqr,ForceMax^2) +...
         5*k*(pknee(2)-0.00)^2*(1-step(pknee(2)-0.00,0));
    end
    
    H = (F.')^(-1)*K*F^(-1);
    G = H;
    L = (Xd - f*alpha).' * G * (Xd - f*alpha) + b;
    
    %---------------------------------------------------------------------%
    
    %---------the following part are generic, don't need to edit----------%
    
    pLx = sym('pLx', [N 1],'real');
    pLxd = sym('pLxd', [N 1],'real');
    for i=1:N
        pLx(i) =  diff(L,X(i));
        pLxd(i) = diff(L,Xd(i));
    end
    
    for i=1:N    
        temp_cell = arrayfun(@char, pLx(i), 'uniform', 0);
        str_temp = strrep(temp_cell{1},'heaviside','heaviside_appr');
        str_temp = strrep(str_temp,'dirac','dirac_appr');
        pLx(i) =  eval( evalin(symengine,str_temp) );
        
        temp_cell = arrayfun(@char, pLxd(i), 'uniform', 0);
        str_temp = strrep(temp_cell{1},'heaviside','heaviside_appr');
        str_temp = strrep(str_temp,'dirac','dirac_appr');
        pLxd(i) =  eval( evalin(symengine,str_temp) );
    end
    
    
    in2 = [X;in];
    %get_EL = matlabFunction([pLx, pLxd],'vars',[tt,X',Xd',in'],'File','generated_files\get_EL');
    %get_G = matlabFunction(G,'vars', {tt,X,in},'File','generated_files\get_G');
    %get_f = matlabFunction(f,'vars',{tt,X},'File','generated_files\get_f');
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