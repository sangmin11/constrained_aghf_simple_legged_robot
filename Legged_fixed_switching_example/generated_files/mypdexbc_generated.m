function [pl,ql,pr,qr] = mypdexbc_generated(xl,ul,xr,ur,t,X) % Boundary condition
%this function is automatically generated...
X0 = X(:,1);
Xf = X(:,end);
free_ind=[0,0;0,0;0,0;0,0;0,0;0,0;1,1;1,1;0,0;0,0;1,1;1,1;0,0;0,0];
%intial states
pl = ul-X0;
pl = pl .* (~free_ind(:,1));
ql = free_ind(:,1);
% final states
pr = ur-Xf;
pr = pr .* (~free_ind(:,2));
qr = free_ind(:,2);
end
