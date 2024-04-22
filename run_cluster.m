clc;
close all;
clear all;

addpath('tool');
load('100leaves.mat');


if exist('Y','var')
    gt = Y;  
end
if exist('y','var')
    gt = y;  
end
if exist('truelabel','var')
    gt = truelabel{1};
end
if exist('data','var')
    X = data;  
end

nV = length(X); 
class_num = length(unique(gt)); 

W = cell(1,nV); 
[data,vecX] = Con_Diagdata(X,nV);
vecX = NormalizeData(vecX); 

L = vecX'*vecX;
sigma = mean(mean(1 - L));
Wgt = ones(size(L)) - exp(-(ones(size(L))-L)/sigma);
for v=1:nV  
    X{v} = NormalizeData(X{v});
    W{v} = Wgt;
end


% 100 10 1 1e-1 1e-2 1e-3 1e-4 1e-5 1e-6 1e-7 1e-8
for beta = [1e-2] 
for lambda = [1e-1] 
mu = 1e-2; 
epsilon = 1e-6;
rho = 2;
q = 1/2;
% q = 2/3;
rng(5,'twister');
tic;
[Coe,Conver_iter,iter] = NLRSC_MSC_Lq_matrix_Cor(X,W,mu,rho,epsilon,lambda,beta,q);
affi = (Coe + Coe')./nV;
time = toc;

Predicted = SpectralClustering(affi, class_num);
result =  ClusteringMeasure(gt, Predicted);


fid=fopen('log_leaves-12-tt.txt','a');


fprintf(fid,'beta: %.1e\n', beta);
fprintf(fid,'lambda: %.1e\n', lambda);
fprintf(fid,'mu: %.1e\n', mu);
fprintf(fid,'eps: %.1e\n', epsilon);
fprintf(fid,'rho: %.1f\n', rho);
fprintf(fid,'q: %.1f\n', q);
fprintf(fid,'result: %.4g %.4g %.4g %.4g %.4g %.4g %.4g \n', result);
fprintf(fid,'time: %.2f\n', time);
fprintf(fid,'\n');

end
end
fprintf(fid,'--------------------------------------------------\n');
fid=fclose('all');

