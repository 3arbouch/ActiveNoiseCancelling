function [ error , MSerror, timeOfConvergence , timeOfComputation, w] = BLMS( referenceSignal, filtredSignal, filterSize,blockSize,  numberOfIterations, convergenceThreshold )
% This function implements The Block LMS algorithm
% Estimate the filter of size filterSize between the reference signal and the filtred signal.
% Block size represents how much data a block contains 
% Returns:
% 1) The error computed at each sample
% 2) The mean sqaure error
% 3) Time of convergence: the time it takes to the MSE to be below the
%     convergence threshold
% 4) Time of computation: the time the algorithm takes for the whole
%   computaion
% 5) The estimated filter: the one that reaches the convergence threshold.
%    Otherwise, the filter estimated at the end of the computation

x = referenceSignal ; 
d = filtredSignal ; 

 i = filterSize  ; 
 alpha = 0.05 ; 
 w = zeros(filterSize, 1) ;
 error = 0 ; 
 MSerror = 0; 
 
tstart = tic; % Used to compute the convergence time 
converged = 0 ; % Used to indicate if the algorithm converged or not 
tic % used to compute the computation time
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
               if(~converged &&  MSE < convergenceThreshold && i> 2*filterSize)
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


