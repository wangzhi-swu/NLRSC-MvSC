function [Q_new] = Half_norm_matrix(Q,lambdaW)

% paper: 1/2 regularized low-rank representation for hyperspectral imagery classification
[m,n] = size(Q);
for i = 1:m
    for j = 1:n
        sigma = Q(i,j);
        lambda = lambdaW(i,j);
        if abs(sigma) > 54.0^(1/3)/4*lambda^(2/3)
            new_sigma = func(sigma, lambda);
        else 
            new_sigma = 0;
        end
        Q_new(i,j) = new_sigma;
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

end