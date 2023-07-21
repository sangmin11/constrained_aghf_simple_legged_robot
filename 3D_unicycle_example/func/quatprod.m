function out = quatprod( q1,q2 )
%QUATPROD implement quaternion product

v1 = q1(2:end);
v2 = q2(2:end);
out = [q1(1)*q2(1)-v1'*v2; q1(1)*v2+q2(1)*v1+cross(v1,v2)];

end

