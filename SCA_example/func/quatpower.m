function out = quatpower( q, n)
% self implemented quaternion power, input quaternion must be col vector
[l1,l2] = size(q);
out = zeros(l1,l2);
out(1) = 1;

for i = 1:n
   out = quatprod(out,q);
end

end

