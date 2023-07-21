

function [p, dp, xfull, xdrift, bufull, budrift, sol, Xs, Xf, Xint, U, t, x, xnew, cost, cost_all, cost_u, T_solve] = solve_HF(model,param_model,tmax,xpoints,intgrids,tpoints,X,T,params,flags,free_ind,extra_params)

m = 0;
x = linspace(0,T,xpoints);                  % discretization of the curve
%t = [0 logspace(-4,log10(tmax),tpoints-1)]; % discretization of time interval
t = [0 logspace(-6,log10(tmax),tpoints-1)];

%---------------generate functions for HF (runtime is short if the functions already existed and not modified)----------------%
tic;
display('generating heat flow model...');
%[get_EL, get_G, get_f, get_F, get_H, get_L, N, M] = model(flags, params, param_model);
[get_f, get_F, get_H, get_L, N, M] = model(flags, params, param_model, extra_params);
toc;

tic;
display('writing PDE functions...');
make_pde_functions(N,flags,T,free_ind);
toc;


display('solving...');
tic;
profile on
%--------------solve HF pde--------------------------%
%opts = odeset('AbsTol',0.002);
opts = odeset('AbsTol',extra_params('AbsTol'),'RelTol',extra_params('RelTol'));
%sol = pdepe(m,@(x,t,u,DuDx) mypdexpde_generated(x,t,u,DuDx,params,get_EL),@(x) mypdexic_generated(x, X),@(xl,ul,xr,ur,t) mypdexbc_generated(xl,ul,xr,ur,t,X),x,t,opts);
sol = pdepe_AGHF(x,t,X,params,opts);
% The solution is of form sol(t,x,i)
profile off
T_solve = toc;
str = sprintf('solving time: %f', T_solve);
disp(str);

%-----------------------------Extract Controls-------------------------%
tic;
display('solution found, generatating path and other outputs...');
xpoints_new = intgrids;
budrift = zeros(M,xpoints_new);
bufull = zeros(N,xpoints_new);
xnew = linspace(0,T,xpoints_new);
p = zeros(N,xpoints_new);
for i = 1:N
    p(i,:) = spline(x,real(sol(end,:,i)),xnew);
end

t0 = 0;
drift0 = get_f(t0,p(:,1));
dp = [drift0 diff(p.').'/(T/(xpoints_new-1))];

for i = 1 : xpoints_new
    tt = t0 + (i-1)*(T/(xpoints_new-1));
    drift = get_f(tt,p(:,i));
    F = get_F(tt,p(:,i));
    bufull(:,i) = F^(-1) * ( dp(:,i) - drift );
    budrift(:,i)= bufull((end-M+1):end,i);
end


%Find x^* by integrating ODE
xfull   =   zeros(N,xpoints_new);
xfull(:,1)  =   p(:,1); %Initial state
xdrift   =   zeros(N,xpoints_new);
xdrift(:,1)  =   p(:,1); %Initial state
ind_nonSO3 = 1:N;
ind_SO3 = [];
ind_w = [];
if (isKey(extra_params,'Rotm_ind0'))
    ind_nonSO3(extra_params('Rotm_ind0'):extra_params('Rotm_ind0')+8) = [];
    ind_SO3 = extra_params('Rotm_ind0'):extra_params('Rotm_ind0')+8;
end
if (isKey(extra_params,'Quat_ind0'))
    ind_nonSO3(extra_params('Quat_ind0'):extra_params('Quat_ind0')+3) = [];
    ind_SO3 = extra_params('Quat_ind0'):extra_params('Quat_ind0')+3;
end
if (isKey(extra_params,'w_ind0'))
    ind_w = extra_params('w_ind0'):extra_params('w_ind0')+2;
end
admissible_control = [];
if (isKey(extra_params,'admissible_control'))
    admissible_control = extra_params('admissible_control');
end

Rtemp = eye(3);
for i =1: xpoints_new-1
    tt = t0 + i*(T/(xpoints_new-1));
    cBfull = get_F(tt,xfull(:,i));
    cBdrift_full = get_F(tt,xdrift(:,i));
    cBdrift = cBdrift_full(:,(end-M+1):end);
    drift = get_f(tt,xdrift(:,i));
    drift_full = get_f(tt,xfull(:,i));
    %xdrift(:,i+1) = xdrift(:,i) + ( drift + cBdrift*budrift(:,i) )*(xnew(i+1)-xnew(i));
    %xfull(:,i+1) = xfull(:,i) + ( drift_full + cBfull*bufull(:,i) )*(xnew(i+1)-xnew(i));
    xdrift(ind_nonSO3,i+1) = xdrift(ind_nonSO3,i) + ( drift(ind_nonSO3) + cBdrift(ind_nonSO3,:)*budrift(:,i) )*(xnew(i+1)-xnew(i));
    xfull(:,i+1) = xfull(:,i) + ( drift_full + cBfull*bufull(:,i) )*(xnew(i+1)-xnew(i));
    if (isKey(extra_params,'Rotm_ind0'))
        Rtemp = reshape(xdrift(ind_SO3,i),3,3);
        xdrift(ind_SO3,i+1) = ToColVec( Rtemp*expm((xnew(i+1)-xnew(i))*skew(xdrift(ind_w,i))) );
        xfull(:,i+1) = xfull(:,i) + ( drift_full + cBfull*bufull(:,i) )*(xnew(i+1)-xnew(i));
    end
    omega = [0 0 0]';
    if (isKey(extra_params,'Quat_ind0'))
        qtemp = xdrift(ind_SO3,i);
        if (isKey(extra_params,'w_ind0'))
            omega = xdrift(ind_w,i);
        else
            %w_star = [0 2*pi 0]';
            %w_star = [0 0 2*pi]';
            w_star = [0 0 0]';
            omega(admissible_control-4) = budrift(:,i) + w_star(admissible_control-4);
        end
        xdrift(ind_SO3,i+1) =  quatprod( qtemp, quatexp( (xnew(i+1)-xnew(i))*[0;omega]/2, 10) ) ;
        xfull(:,i+1) = xfull(:,i) + ( drift_full + cBfull*bufull(:,i) )*(xnew(i+1)-xnew(i));
    end
    
end

%-----------------calculate cost------------------------------%
X_temp = zeros(N,xpoints);
dX_temp = zeros(N,xpoints);
cost = zeros(1,tpoints); %total cost(includes barrier and full control)
cost_all = zeros(1,tpoints);% cost using full control only
cost_u = zeros(1,tpoints); % cost using only admissible control
for kk = 1:tpoints
    for j = 1:N
    [X_temp(j,:), dX_temp(j,:)]=pdeval(m,x,sol(kk,:,j),x);
    end
    for i = 1:xpoints
        tt = t0 + (i-1)*(T/(xpoints-1));
        X = X_temp(:,i);
        dX = dX_temp(:,i);
        f = get_f(tt,X);
        F = get_F(tt,X);
        G = get_G(tt,X,params');
        H = get_H(tt,[X;params']);
        L = get_L(tt,X,dX,params');
        ufull = F^(-1) * ( dX - f );
        udrift= ufull((end-M+1));
        cost_all(kk) = cost_all(kk) + (dX-f)'*H*(dX-f)*(T/(xpoints-1));
        cost_u(kk) = cost_u(kk) + udrift'*udrift*(T/(xpoints-1));
        %cost(kk) = cost(kk) + (dX-f)'*G*(dX-f)*(T/(xpoints-1));
        cost(kk) = cost(kk) + L*(T/(xpoints-1));
        
    end
end

% interpolate to original grids

Xs = zeros(N,xpoints); % HF solution
Xf = zeros(N,xpoints); % integreated path using full control
Xint = zeros(N,xpoints); % integreated path using admissible control
U = zeros(N,xpoints); % full control
for i = 1:N
    Xs(i,:) = interp1(xnew,p(i,:),x);
    Xf(i,:) = interp1(xnew,xfull(i,:),x);
    Xint(i,:) = interp1(xnew,xdrift(i,:),x);
    U(i,:) = interp1(xnew,bufull(i,:),x);
end

toc;
end


function make_pde_functions(N,flags,T,free_ind)    % Define PDE; Evaluate right-hand-side of GHF

%----------------generate pde function-------------------------------%

%id_mypdexpde = fopen('mypdexpde_generated.m','w');
%st = sprintf('%s\n','function [c,f,s] = mypdexpde_generated( x,t,u,DuDx,params,get_EL )');
id_mypdexpde = fopen('generated_files\mypdexpde_generated.m','w');
st = sprintf('%s\n','function [c,f,s] = mypdexpde_generated( x,t,u,DuDx,params)');
st = [st sprintf('%s\n','%this function is automatically generated...')];
st = [st sprintf('%s\n','N = length(u);')];

st1 = 'EL = get_EL(x,u(1)';
for i = 2:N
    st1 = [st1 ',u(' num2str(i) ')'];
end
for i = 1:N
    st1 = [st1 ',DuDx(' num2str(i) ')'];
end
for i = 1:length(flags)
   if(flags(i)==1)
       st1 = [st1 ',params(' num2str(i) ')'];
   end
end
st1 = [st1 ');'];
st = [st sprintf('%s\n',st1)];

st = [st sprintf('%s\n','pLx = EL(:,1);')];
st = [st sprintf('%s\n','pLxd = EL(:,2);')];
%st = [st sprintf('%s\n','G = get_G(x,u,transpose(params));')];
st = [st sprintf('%s\n','c = ones(N,1);')];
st = [st sprintf('%s\n','f = pLxd;')];
st = [st sprintf('%s\n','s = -pLx;')];
%st = [st sprintf('%s\n','f = G^(-1)*pLxd;')];
%st = [st sprintf('%s\n','s = -G^(-1)*pLx;')];
st = [st sprintf('%s\n','end')];

fprintf(id_mypdexpde,'%s', st);
fclose(id_mypdexpde);

%------------------generate ic function------------------------------%

id_mypdexic = fopen('generated_files\mypdexic_generated.m','w');
st = sprintf('%s\n','function u0 = mypdexic_generated(x, X)');
st = [st sprintf('%s\n','%this function is automatically generated...')];
st = [st sprintf('%s\n','X0 = X(:,1);')];
st = [st sprintf('%s\n','Xf = X(:,end);')];
st = [st sprintf('%s\n','N = length(X(1,:));')];
st = [st sprintf('%s\n','n = length(X(:,1));')];
st = [st sprintf('%s\n','u0 = zeros(n,1);')];

st1 = ['T=' num2str(T) ';'];
st = [st sprintf('%s\n',st1)];

st = [st sprintf('%s\n','if ( N <= 2) % if only boundary condition are specified, use a stright line')];
st = [st sprintf('%s\n','    u0=X0+(Xf-X0)*(x/T) + 0.00*sin(x/T*2*pi)*ones(n,1);')];
st = [st sprintf('%s\n','else % else if a initial guess curve is specified')];
st = [st sprintf('%s\n','    for i=1:n')];
st = [st sprintf('%s\n','        u0(i)=interp1( linspace(0,1,N), X(i,:), x/T);')];
st = [st sprintf('%s\n','    end')];
st = [st sprintf('%s\n','end')];
st = [st sprintf('%s\n','end')];

fprintf(id_mypdexic,'%s', st);
fclose(id_mypdexic);


%------------------generate bc function------------------------------%

id_mypdexbc = fopen('generated_files\mypdexbc_generated.m','w');
st = sprintf('%s\n','function [pl,ql,pr,qr] = mypdexbc_generated(xl,ul,xr,ur,t,X) % Boundary condition');
st = [st sprintf('%s\n','%this function is automatically generated...')];
st = [st sprintf('%s\n','X0 = X(:,1);')];
st = [st sprintf('%s\n','Xf = X(:,end);')];

[m,n] = size(free_ind);
st1 = 'free_ind=[';
for i = 1:m
    st1 = [st1 num2str(free_ind(i,1))];
    for j = 2:n
        st1 = [st1 ',' num2str(free_ind(i,j))];
    end
    if(i~=m) st1 = [st1 ';']; end
end
st1 = [st1 '];'];

st = [st sprintf('%s\n',st1)];
st = [st sprintf('%s\n','%intial states')];
st = [st sprintf('%s\n','pl = ul-X0;')];
st = [st sprintf('%s\n','pl = pl .* (~free_ind(:,1));')];
st = [st sprintf('%s\n','ql = free_ind(:,1);')];
st = [st sprintf('%s\n','% final states')];
st = [st sprintf('%s\n','pr = ur-Xf;')];
st = [st sprintf('%s\n','pr = pr .* (~free_ind(:,2));')];
st = [st sprintf('%s\n','qr = free_ind(:,2);')];
st = [st sprintf('%s\n','end')];

fprintf(id_mypdexbc,'%s', st);
fclose(id_mypdexbc);



end
