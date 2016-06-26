function [X_matrix, autocorrX, mu] = formX( filterSize, index, data, N )
% This function construct shifted versions of the data 


X_matrix = data(index:-1:index- filterSize+1)' ; 
 mu = 1/(X_matrix*X_matrix') ; 
for j=1:N-1 
dataToAdd = data(index+j:-1:index -filterSize+1+j) ; 
 mu = [mu, 1/(dataToAdd'*dataToAdd)] ; 
X_matrix = [X_matrix ;dataToAdd']   ; 
end

autocorrX = (1/N).* X_matrix'*X_matrix ;   
end

