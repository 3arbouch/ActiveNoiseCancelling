function [ error , MSerror, timeOfConvergence, timeOfComputation , w ] = NLMS( referenceSignal, filtredSignal, filterSize, numberofIterations, convergenceThreshold )
% This function implement the NLMS algorithm.
% Estimate the filter of size filterSize between the reference signal and the filtred signal.
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


i = 1 ;

% initial guess of the filter 
w = zeros(filterSize, 1) ; 
% Vector that contain save the Mean square error
MSerror = zeros(numberofIterations, 1) ;

error = 0 ;
tstart = tic; % Used to compute the convergence time
converged = 0 ;
tic % used to compute the computation time 
while(i<numberofIterations)
               x_n = x(i+filterSize-1:-1:i) ; 
               e = d(i+filterSize-1) - x_n'*w  ;  
               mu = 1/(x_n'*x_n)   ; 
               w = w + mu.* x_n*e ; 
              error = [error ; e] ; 
              MSerror(i) = (error'*error)/length(error) ;
              % Test if the mean square error is below the convergence
              % threshold 
               if(~converged &&  MSerror(i) < convergenceThreshold && i > filterSize)
                  timeOfConvergence = toc(tstart) ; 
                  converged = 1;
                  w_convergence = w ;
               end
              
              
               
  i = i+1  ;             
end
timeOfComputation = toc ; 
 if(~converged)
      timeOfConvergence = toc(tstart) ;  
 else 
     w = w_convergence ;
 end

end

