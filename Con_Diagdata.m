function [ X,vecX ] = Con_Diagdata( data,nV )
    D = 0;
    for i = 1:nV
        D_set(i) = size(data{i},1); 
        D = D + D_set(i);
    end
    N = size(data{i},2); 
    dim_X = [D,N,nV]; 
    X = zeros(D,N,nV);
    vecX = zeros(D,N);
    for i = 1:N  
        B = zeros(D,nV);
        star_row = 1;
        for j = 1:nV  
            end_row = star_row + D_set(j)-1;
            B(star_row:end_row,j) = data{j}(:,i);
            vecX(star_row:end_row,i) = data{j}(:,i);
            star_row = end_row + 1;
        end
        X(:,i,:) = B; 
    end
end

