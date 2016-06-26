function [ error, MSerror, timeOfConvergence,timeOfComputation,  w1, w2, filters1, filters2] = NLMSSterio( referenceSignal1,referenceSignal2, filtredSignal, fs, filterSize, convergenceThreshold)
% This function implements the Sterio version of NLMS algorithm.
% Estimate the filter of size filterSize between the reference signal and the filtred signal.
% Returns:
% 1) The error computed at each sample
% 2) The mean sqaure error
% 3) Time of convergence: the time it takes to the MSE to be below the
%     convergence threshold
% 4) Time of computation: the time the algorithm takes for the whole
%   computaion
% 5) The estimated filters: the one that reaches the convergence threshold.
%    Otherwise, the filter estimated at the end of the computation
% 6) Matrix containing all the estimated filters at each iteration



x1 = referenceSignal1 ;
x2 = referenceSignal2 ; 
d = filtredSignal ; 

% Append reference Signal with zeros of length frame Size 
x1 = [x1; zeros(filterSize,1)] ; 
x2 = [x2; zeros(filterSize,1)] ; 
d = [d; zeros(filterSize,1)] ; 
%Define the number of iterations
numberofIterations = min(length(x1), length(d))-filterSize +1 ; 

i = 1 ;

% initial guess of the filter 
w1 = zeros(filterSize, 1) ; 
w2 = zeros(filterSize, 1) ; 


frameSize=1024 ; 

% Vector that contain save the Mean square error. Line Vector
MSerror = zeros(1,floor(numberofIterations/frameSize)) ;

%Matrix that contains the estimated filters each frameSize
filters1 = zeros(floor(numberofIterations/frameSize),filterSize) ; 
filters2 = zeros(floor(numberofIterations/frameSize),filterSize) ; 



% initialize the error to zero. Line Vector.
error = zeros(1,numberofIterations) ;

%Initialize the gradient norm. Line Vector.


% Used to compute the convergence time
tstart = tic; 

% Indicates if the algorithm has converged or not 
converged = 0 ;

% used to compute the computation time 
tic 


while(i<numberofIterations)
    
               % Get two blocks of signals
               x_n1 = x1(i+filterSize-1:-1:i) ; 
               x_n2 = x2(i+filterSize-1:-1:i) ; 
                 
               x_n = [x_n1' x_n2'] ; 
               mu= 1/(x_n*x_n')  ;
               
               
               
         
               % Calculate the error
               e = d(i+filterSize-1) - x_n1'*w1 - x_n2'*w2 ;  
               error(i)= e ;  
                
               % calculate the adaptation
               w1 = w1 + mu.* x_n1*e ;
               w2 = w2 + mu.*x_n2*e ; 


               % Save the estimated filter if recahed frameSize iteration
               if(mod(i,frameSize)==0)
                  m = i/frameSize ;
                  
                  filters1(m,:)=w1;
                  filters2(m,:)=w2;
                  MSerror(m) = (error*error')/i ;
                  
                           
%               Test if the mean square error is below the convergence
%               threshold
               if(~converged &&  MSerror(m) < convergenceThreshold && m > filterSize)
                  timeOfConvergence = toc(tstart) ; 
                  converged = 1;
                  w_convergence1 = w1 ;
                  w_convergence2 = w2 ;
               end
        
               end

  i = i+1  ;     

end
timeOfComputation = toc ; 
 if(~converged)
      timeOfConvergence = toc(tstart) ;  
 else 
     w1 = w_convergence1 ;
     w2 = w_convergence2 ; 
 end

end

