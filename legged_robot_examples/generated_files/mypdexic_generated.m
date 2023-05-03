function u0 = mypdexic_generated(x, X)
%this function is automatically generated...
X0 = X(:,1);
Xf = X(:,end);
N = length(X(1,:));
n = length(X(:,1));
u0 = zeros(n,1);
T=2;
if ( N <= 2) % if only boundary condition are specified, use a stright line
    u0=X0+(Xf-X0)*(x/T) + 0.00*sin(x/T*2*pi)*ones(n,1);
else % else if a initial guess curve is specified
    for i=1:n
        u0(i)=interp1( linspace(0,1,N), X(i,:), x/T);
    end
end
end
