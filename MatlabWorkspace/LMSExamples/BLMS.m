function [ error , MSerror, timeOfConvergence , w] = BLMS( referenceSignal, filtredSignal, filterSize,blockSize,  numberOfIterations )
% This function implements The Block LMS

x = referenceSignal ; 
d = filtredSignal ; 

 i = filterSize  ; 
 alpha = 0.9 ; 
 w = zeros(filterSize, 1) ;
 error = 0 ; 
 MSerror = 0; 
 epsilon = 10^-7 ; 
 timeOfConvergence =0 ; 
 tic
 while (i<numberOfIterations)
     [x_matrix, autocorrX]=formX(filterSize, i, x, blockSize) ;
     
             E=eig(autocorrX)   ; 
             mu =alpha* 1/max(E) ;
             e = d(i:i + blockSize -1)-x_matrix*w   ;   
             error = [error ; e] ; 
%              w = w + 2*mu*(1/blockSize)*x_matrix'*e  ;  
              w = w + 2*mu*(1/blockSize)*(x_matrix'*d(i:i + blockSize -1)-blockSize*autocorrX*w)  ; 
             MSE = (error'*error)/length(error) ; 
             MSerror =  [MSerror ;repmat(MSE,blockSize,1) ] ; 
               if(MSE < epsilon)
                  timeOfConvergence = toc ;  
               end
     
     
     i = i+ blockSize ; 
 end

if(timeOfConvergence==0)
    timeOfConvergence = toc ; 
end
           

end


function [X_matrix, autocorrX] = formX( filterSize, index, data, N )
% This function construct shifted versions of the data 


X_matrix = data(index:-1:index- filterSize+1)' ; 
for j=1:N-1 
dataToAdd = data(index+j:-1:index -filterSize+1+j) ; 
X_matrix = [X_matrix ;dataToAdd']   ; 
end

autocorrX = (1/N).* X_matrix'*X_matrix ;   
end

