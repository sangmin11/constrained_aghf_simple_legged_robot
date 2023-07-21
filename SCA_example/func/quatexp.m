function out = quatexp( q , order)
%self implemented version of quaternion exponential, input quaternion must
%be col vector

out = [1 0 0 0]';

for i = 1:order
   out = out +  quatpower(q,i)/(factorial(i));
end


end

