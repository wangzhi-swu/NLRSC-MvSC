function [E] = solve_l12l2(W,lambda)
n = size(W,2);
E = W;
for i=1:n
    E(:,i) = solve_l2(W(:,i),lambda);
end
end

function [x] = solve_l2(w,lambda)
% min lambda |x|_2 + |x-w|_2^2
nw = norm(w);
if nw > 54.0^(1/3)/4*lambda^(2/3)
    x = (func(nw,lambda))*w/nw;
else
    x = zeros(length(w),1);
end
end

function [R] = func(sigma, lambda)
temp = lambda/8*(abs(sigma)/3)^(-2/3);
if lambda/8*(abs(sigma)/3)^(-2/3) > 1
    temp = 1;
elseif lambda/8*(abs(sigma)/3)^(-2/3) < -1
    temp = -1;
else
    temp = temp;
end

R_temp1 = acos(temp);  %fai
R_temp2 = 1+cos(2*pi/3-2/3*R_temp1);
R = 2/3*sigma*R_temp2;
end