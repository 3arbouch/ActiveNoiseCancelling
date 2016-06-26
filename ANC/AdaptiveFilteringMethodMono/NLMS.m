function [ error , MSerror, timeOfConvergence, timeOfComputation , w, filters ] = NLMS( referenceSignal, filtredSignal, filterSize, convergenceThreshold, sparseCondition )
% This function implement two versions of the NLMS algorithm. the classic
% non constraint one and the L0 sparse conditionning one ( For sparse
% filtering estimation).
% Sparse condition =0,1 (0 for non sparse, 1 for L0 norm sparse solution)
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
% 6) Matrix containing all the estimated filters at each iteration

x = referenceSignal ; 
d = filtredSignal ; 

% Append reference Signal with zeros of length frame Size 
x = [x; zeros(filterSize,1)] ; 
d = [d; zeros(filterSize,1)] ; 

%Define the number of iterations
numberofIterations = min(length(x), length(d))-filterSize +1 ; 

i = 1 ;

% initial guess of the filter 
w = zeros(filterSize, 1) ; 

% Vector that contain save the Mean square error. Line Vector
MSerror = zeros(1,numberofIterations) ;

% Define the frame Size.
frameSize = 1024 ; 

%Matrix that contains the estimated filters each frameSize
filters = zeros(floor(numberofIterations/frameSize),filterSize) ; 

%initialize parameters for the NLMS Algorithm

    %Determine the shape of the exponantial function. The greater, the
    %higher is the peak. Please refer to documentation for more details 
    %Value fixed to 15: deduced from expeiments to
    %obtain the most stable solution.
    beta = 15 ;
    
    %Determine how much the non sparse solution is penalized: Penalysing
    %factor. PLease refer to documentation for more details.
    %Value fixed to 10^-6: deduced from experiemtns to obtain the
    %most stable solution.
    k = 10^-6; 
    
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
                
               % Get the reference microphone signal and the error signal
               x_n = x(i+filterSize-1:-1:i) ; 
               e = d(i+filterSize-1) - x_n'*w  ;  
               
               % calculate the step Size
               mu = 1/(x_n'*x_n)   ;
               
               error(i)= e ;  
             
               if(sparseCondition == 0)
                   % Non Sparse adaptation
                   w = w + mu.* x_n*e ;
                   MSerror(i) = (error*error')/i ;
               else
                   % L0 Norm adaptation
                  w = w + mu.* x_n*e +k*ExponentialAprrox(beta,w) ;
                  MSerror(i) =(error*error')/i ;
               end
               
               % Save the estimated filter if recahed frameSize iteration
               if(mod(i,frameSize)==0)
                  filters(i/frameSize,:)=w; 
               end
               
         
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

