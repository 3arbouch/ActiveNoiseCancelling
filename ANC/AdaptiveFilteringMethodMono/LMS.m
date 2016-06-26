function [ error , MSerror, timeOfConvergence, timeOfComputation , w, filters ] = LMS( referenceSignal, filtredSignal, filterSize, convergenceThreshold )
% This function implement the basic version of the LMS Algorithm.
% Estimate the filter of size filterSize between the reference signal and the filtred signal.
% Convergence Threshold indicates when the adaotation stops after Mean
% Square error < convergence Threshold
% Returns:
% 1) The error computed at each sample
% 2) The mean sqaure error
% 3) Time of convergence: the time it takes to the MSE to be below the
%     convergence threshold
% 4) Time of computation: the time the algorithm takes for the whole
%   computaion
% 5) The estimated filter: the one that reaches the convergence threshold.
%    Otherwise, the filter estimated at the end of the computation
% 6) Matrix containing all the estimated filters at each iteration

x = referenceSignal ; 
d = filtredSignal ; 

% Append reference Signal with zeros of length frame Size 
x = [x; zeros(filterSize,1)] ; 
d = [d; zeros(filterSize,1)] ; 

%Define the number of iterations
numberofIterations = min(length(x), length(d))-filterSize +1 ; 


% Define the frame Size.
frameSize = 1024 ; 

% This represents the step size of the LMS Algorithm
mu = 0.001;


% initial guess of the filter 
w = zeros(filterSize, 1) ; 

% Vector that contain save the Mean square error. Line Vector
MSerror = zeros(1,floor(numberofIterations/frameSize)) ;


%Matrix that contains the estimated filters each frameSize
filters = zeros(floor(numberofIterations/frameSize),filterSize) ; 

    
% initialize the error to zero. Line Vector.
error = zeros(1,numberofIterations) ;
    
    
    
% Used to compute the convergence time
 tstart = tic; 
    
% Indicates if the algorithm has converged or not 
converged = 0 ;

    % used to compute the computation time 
    tic 
    
    
    i = 1 ;   
while(i<numberofIterations)
    
                % get the reference and the error signal
               x_n = x(i+filterSize-1:-1:i) ; 
               e = d(i+filterSize-1) - x_n'*w  ;  
               
               error(i)= e ;  
               
               % Filter update
               w = w + mu.* x_n*e ;
               
               % Save the estimated filter if recahed frameSize iteration
               if(mod(i,frameSize)==0)    
                  m = i/frameSize ;
                  filters(m,:)=w;
                  MSerror(m) = (error*error')/i ;
                  
                    % Test if the mean square error is below the convergence
                    % threshold
                    if(~converged &&  MSerror(m) < convergenceThreshold && m > filterSize)
                        timeOfConvergence = toc(tstart) ; 
                        converged = 1;
                        w_convergence = w ;
                    end
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

