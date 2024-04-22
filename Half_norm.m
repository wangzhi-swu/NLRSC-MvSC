function [U,new_sigma,V] = Half_norm(W,lambda)

% paper: 1/2 regularized low-rank representation for hyperspectral imagery classification
[U,S,V] = svd(W,'econ');

sigma = diag(S);

len = length(sigma);

new_sigma = [];
for i = 1:len
    if abs(sigma(i)) > 54.0^(1/3)/4*lambda^(2/3)
        new_sigma(i) = func(sigma(i), lambda);
    else 
        new_sigma(i) = 0;
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