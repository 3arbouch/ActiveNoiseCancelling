function [ error, MSerror, timeOfConvergence,timeOfComputation,  w1, w2, filters1, filters2, weights1, weights2, variances1, variances2] = NLMSSterioL1Weighted( referenceSignal1,referenceSignal2, filtredSignal, fs, filterSize, convergenceThreshold, mu_k1, mu_k2, interval1, interval2, upperLimit, memoryInSeconds )
% This function implements the L1 weighted Sterio version of NLMS algorithm.
% mu_k1, mu_k2 are the mean of the respective weighting function 
% interval1, interval2 represent the prior allwed error interval in
% samples: Usually set to filterSize
% Upper limit represents the upper limit of the weighting function 
% MemoryInSeconds represents how many previous weighting functions are
% taken into account fot the mean and variance calculation
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

weights1 = zeros(floor(numberofIterations/frameSize),filterSize) ; 
weights2 = zeros(floor(numberofIterations/frameSize),filterSize) ; 

variances1 = zeros(floor(numberofIterations/frameSize),1) ; 
variances2 = zeros(floor(numberofIterations/frameSize),1) ;

% initialize the error to zero. Line Vector.
error = zeros(1,numberofIterations) ;

%Initialize the gradient norm. Line Vector.


% Used to compute the convergence time
tstart = tic; 

% Indicates if the algorithm has converged or not 
converged = 0 ;

% used to compute the computation time 
tic 


%% Form k as a gaussian distribution

sigma_k1 = interval1/2 ; 
sigma_k2 = interval2/2 ; 

x = [1:filterSize]' ;
norm1 = normpdf(x,mu_k1,sigma_k1);
norm2 = normpdf(x,mu_k2,sigma_k2);


memory = memoryInSeconds;
memorySamples = floor((memory*fs)/frameSize) ; 
mu_buffer1 = [((interval1/2).*randn(memorySamples,1) + mu_k1); zeros(floor(numberofIterations/frameSize),1)];
mu_buffer2 = [((interval2/2).*randn(memorySamples,1) + mu_k2); zeros(floor(numberofIterations/frameSize),1)];



k1 = (max(norm1) - norm1)*(upperLimit/max(norm1));
k2 = (max(norm2) - norm2)*(upperLimit/max(norm2));


samplesNoise = 100 ;
sigma = samplesNoise/2 ; 
noise = sigma*randn(memorySamples + 1,1);





while(i<numberofIterations)
%      i
               % Get two blocks of signals
               x_n1 = x1(i+filterSize-1:-1:i) ; 
               x_n2 = x2(i+filterSize-1:-1:i) ; 
                
               x_n = [x_n1' x_n2'] ; 
               mu= 1/(x_n*x_n')  ;
               
               
                
               
               % Calculate the error
               e = d(i+filterSize-1) - x_n1'*w1 - x_n2'*w2 ;  
               error(i)= e ;  
                
               % calculate the adaptation
               w1 = w1 + mu*x_n1*e - k1.*sign(w1) ;
               w2 = w2 + mu*x_n2*e - k2.*sign(w2) ;

               
               
               
        
               
               % Save the estimated filter if recahed frameSize iteration
               if(mod(i,frameSize)==0)
                  m = i/frameSize ;
                  
                  filters1(m,:)=w1;
                  filters2(m,:)=w2;
                  MSerror(m) = (error*error')/i ;
                  
                  
                  % Update k1 and k2
               ignoreSamples = 70 ; 
               [~,mu_k1 ] = max(abs(w1(ignoreSamples:end-ignoreSamples)));
               [~,mu_k2 ] = max(abs(w2(ignoreSamples:end-ignoreSamples))) ;
               

               mu_k1 = mu_k1 + ignoreSamples ; 
               mu_k2 = mu_k2 +ignoreSamples ; 
                
                
               mu_buffer1(m+memorySamples) = mu_k1 ; 
               mu_buffer2(m+memorySamples) = mu_k2 ; 
       
               muVector1 = mu_buffer1(m:m+memorySamples) + noise  ; 
               muVector2 = mu_buffer2(m:m+memorySamples) + noise ; 
               
               mu_k1 = mean(muVector1)  ;
               sigma_k1 = std(muVector1) ; 
               
               
               mu_k2 = mean(muVector2)  ;
               sigma_k2 = std(muVector2) ;
        
               norm1 = normpdf(x,mu_k1,sigma_k1);
               norm2 = normpdf(x,mu_k2,sigma_k2);
               
               k1 = (max(norm1) - norm1)*(upperLimit/max(norm1));
               k2 = (max(norm2) - norm2)*(upperLimit/max(norm2));
               
               weights1(m,:) = k1 ; 
               weights2(m,:) = k2 ; 
               
               variances1(m) = sigma_k1 ; 
               variances2(m) = sigma_k2 ; 
        
               end
               
     
         
%               Test if the mean square error is below the convergence
%               threshold
%                if(~converged &&  MSerror(i) < convergenceThreshold && i > filterSize)
%                   timeOfConvergence = toc(tstart) ; 
%                   converged = 1;
%                   w_convergence1 = w1 ;
%                   w_convergence2 = w2 ;
%                end
              
              
               
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

