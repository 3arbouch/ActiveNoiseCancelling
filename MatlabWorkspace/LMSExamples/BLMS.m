function [ error , MSerror, timeOfConvergence , timeOfComputation, w] = BLMS( referenceSignal, filtredSignal, filterSize,blockSize,  numberOfIterations, convergenceThreshold )
% This function implements The Block LMS
tic
x = referenceSignal ; 
d = filtredSignal ; 

 i = filterSize  ; 
 alpha = 0.05 ; 
 w = zeros(filterSize, 1) ;
 error = 0 ; 
 MSerror = 0; 
 
tstart = tic;
converged = 0 ;
 while (i<numberOfIterations)
            [x_matrix, autocorrX]=formX(filterSize, i, x, blockSize) ;
            % E=eig(autocorrX)  
           %  F = fft(autocorrX);
            mu =alpha* (blockSize/(trace(autocorrX))) ;
             %mu = 1/max(max((abs(F))));
             e = d(i:i + blockSize -1)-x_matrix*w   ;   
             error = [error ; e] ; 
%             w = w + 2*mu*(1/blockSize)*x_matrix'*e  ;  
             w = w + 2*mu*((1/blockSize)*x_matrix'*d(i:i + blockSize -1)-autocorrX*w)  ; 
             MSE = (error'*error)/length(error) ; 
             MSerror =  [MSerror ;repmat(MSE,blockSize,1) ] ; 
               if(~converged &&  MSE < convergenceThreshold && i> filterSize +1000)
                 timeOfConvergence = toc(tstart) ;  
                  converged = 1;
                  w_convergence = w ; 
               end
     
     
     i = i+ blockSize ; 
 end
 
timeOfComputation = toc ; 
 if(~converged)
      timeOfConvergence = toc(tstart) ;  
      
 else 
     w = w_convergence ;
 
 end
           

end


function [X_matrix, autocorrX] = formX( filterSize, index, data, N )
% This function construct shifted versions of the data 

x_matrix = zeros(N, filterSize);

for j=0:N-1 
dataToAdd = data(index+j:-1:index -filterSize+1+j) ; 
x_matrix(j+1,1:filterSize) = dataToAdd ;
end
X_matrix = x_matrix ;
autocorrX = (1/N).* X_matrix'*X_matrix ;   
end


