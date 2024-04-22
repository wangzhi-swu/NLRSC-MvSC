function [Q_new] = L23_norm_matrix(Q,lambdaW)


[m,n] = size(Q);
for i = 1:m
    for j = 1:n
        sigma = Q(i,j);
        lambda = lambdaW(i,j);
        if abs(sigma) > 48^(1/4)/3*lambda^(3/4)
            new_sigma = func(sigma, lambda);
        else 
            new_sigma = 0;
        end
        Q_new(i,j) = new_sigma;
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


end