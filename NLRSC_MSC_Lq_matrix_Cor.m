function [Coe,Conver_iter,iter] = NLRSC_MSC_Lq_matrix_Cor(X,W,mu,rho,epsilon,lambda,beta,q)
nV = length(X);
n = size(X{1},2);  
for v = 1:nV
    C{v} = zeros(n,n); 
    P{v} = zeros(n,n); 
    S{v} = zeros(n,n); 
    M1{v} = zeros(n,n); 
    M2{v} = zeros(n,n); 
    dV = size(X{v},1);
    M3{v} = zeros(dV,n); 
    E{v} = zeros(dV,n);
end
Coe = zeros(n,n); 
mu_max = 1e10;
maxIter = 500;
iter = 0;

while iter<maxIter
%     finish = 0;
    iter = iter + 1;
    for v = 1:nV
        % update P
        A = C{v} + M1{v}./mu;
        if q == 2/3
            [U,sigma,V] = Two_thirds_norm(A,2/mu);
            P{v} = U*diag(sigma)*V';
        elseif q == 1/2
            [U,sigma,V] = Half_norm(A,2/mu);
            P{v} = U*diag(sigma)*V';
        end
        
        % update S
        B = C{v}+M2{v}./mu;
        tao = ((2*lambda)/mu).*W{v};
        if q == 2/3
            S{v} = L23_norm_matrix(B,tao);
        elseif q == 1/2
            S{v} = Half_norm_matrix(B,tao);
        end
        
        % update C
        Chat = X{v}'*(mu.*X{v}-mu.*E{v}+M3{v})+mu.*P{v}-M1{v}+mu.*S{v}-M2{v};
        C{v} = inv((2*mu).*eye(n)+mu.*(X{v}'*X{v}))*Chat;
        
        % update E
        tempE = X{v} - X{v}*C{v} + M3{v}./mu;
        if q == 2/3
            new_E = solve_l1l2(tempE,beta/mu);
            beta_W = ones(size(E{v},1),n).*((2*beta)/mu);
            E{v} = L23_norm_matrix(new_E,beta_W);
        elseif q == 1/2
            new_E = solve_l1l2(tempE,beta/mu);
            beta_W = ones(size(E{v},1),n).*((2*beta)/mu);
            E{v} = Half_norm_matrix(new_E,beta_W);
        end
   
        % Check
        stopC = max([max(max(abs(C{v} - P{v}))),max(max(abs(C{v} - S{v}))),max(max(abs(X{v} - X{v}*C{v}-E{v})))]); %两个约束项
        stopC_Views(1,v) = stopC;
   
        % update M1 M2
        M1{v} = M1{v} + mu.*(C{v} - P{v});
        M2{v} = M2{v} + mu.*(C{v} - S{v});
        M3{v} = M3{v} + mu.*(X{v} - X{v}*C{v}-E{v});
        mu = min(rho*mu, mu_max);
         
    end
    
    Conver_max = max(stopC_Views);  
    Conver_iter(iter,1) = Conver_max;
    
    finish = 1;
    for v = 1:nV
        if stopC_Views(1,v) >= epsilon 
            finish = 0;
            break;
        end
    end
    if finish == 1
        for v = 1:nV
            Coe = Coe + sqrt(C{v}.^2); 
        end
        break; 
    end
end
end