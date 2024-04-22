function [E] = solve_l23l2(W,lambda)
n = size(W,2);
E = W;
for i=1:n
    E(:,i) = solve_23(W(:,i),lambda);
end

end

function [x] = solve_23(w,lambda)
% min lambda |x|_2 + |x-w|_2^2
nw = norm(w);
if nw > 48^(1/4)/3*lambda^(3/4)
    x = (func(nw,lambda))*w/nw;
else
    x = zeros(length(w),1);
end
end

function [R] = func(sigma, lambda) 
if lambda == 0
    lambda = lambda + 1e-16;
end
x = 27*sigma*sigma/(16*lambda^(3/2));
fai = acosh(x);
Fai=(2/3^(1/2))*lambda^(1/4)*(cosh(fai/3))^(1/2);
h = abs(Fai)+(2*abs(sigma)/abs(Fai)-abs(Fai)*abs(Fai))^(1/2);
ht = h^3/8;
R=sign(sigma)*ht;
end