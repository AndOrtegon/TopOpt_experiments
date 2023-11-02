function xG = xGraph(x,n)
    M = zeros(n,n,4) ;
    
    M(:,:,1) = (1:2:2*n)'.*(2*n-1:-2:1)/((2*n-1)^2+3) ;
            
    M(:,:,2) = flip(M(:,:,1),2) ;
    M(:,:,3) = flip(M(:,:,2),1) ;
    M(:,:,4) = flip(M(:,:,3),2) ;
    
    xG = zeros(n*size(x,2),n*size(x,3)) ;
    
    for i = 1:size(x,2) 
        for j = 1:size(x,3)
            for k = 1:4
                xG(n*(i-1)+1:n*i,n*(j-1)+1:n*j) = ...
                    xG(n*(i-1)+1:n*i,n*(j-1)+1:n*j) + x(k,i,j)*M(:,:,k) ;
            end
        end
    end    
end