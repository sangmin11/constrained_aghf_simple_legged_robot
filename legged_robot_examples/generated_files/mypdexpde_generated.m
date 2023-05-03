function [c,f,s] = mypdexpde_generated( x,t,u,DuDx,params)
%this function is automatically generated...
global get_EL get_G
N = length(u);
EL = get_EL(x,u(1),u(2),u(3),u(4),u(5),u(6),u(7),u(8),u(9),u(10),u(11),u(12),u(13),u(14),DuDx(1),DuDx(2),DuDx(3),DuDx(4),DuDx(5),DuDx(6),DuDx(7),DuDx(8),DuDx(9),DuDx(10),DuDx(11),DuDx(12),DuDx(13),DuDx(14),params(1),params(2));
pLx = EL(:,1);
pLxd = EL(:,2);
c = ones(N,1);
f = pLxd;
s = -pLx;
end
