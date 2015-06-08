% This script is used to simulate the real time setting:
% The adaptation and the filtering opeerations made in real time 


close all 
clear all 
% define the reference Signal 
fs = 8000;
endTime = 10 ;
x=0.6*randn(1,endTime*fs) ;


% Define the filters: S2(z) and P(z)


filterSize = 256 ;
delayPerSample = 2*filterSize ;
[S2,~,~]=fircband(filterSize,[0 0.4 0.5 1], [1 1 0 0], [1 0.2],{'w' 'c'});
S2=1;
P = (1/filterSize)* ones(1,filterSize) ; 
y = filter (S2,1,x) ;
y = filter(P,1,y) ; 

% Divide by the lenght of P 
P1 = [zeros(1,delayPerSample),P] ;

d = filter(S2,1,x);
d = filter(P1,1,d);

filterToEstimate = filter(S2,1,P); 
figure ; 
plot(y) ; 
hold on 
plot(d(delayPerSample+1:end))

% Define the vector to store past signals 
past = zeros(1,4*filterSize);

numberOfIterations = (length(x)/(2*filterSize)) -1 ;

 % Define parameters of FDAF algorithm
 
  W0 = zeros(2*filterSize, 1) ;
  w = zeros(filterSize,1) ;
  P_k= ones(2*filterSize,1)  ; 
  g=[eye(filterSize),zeros(filterSize,filterSize); zeros(filterSize,filterSize), zeros(filterSize,filterSize)] ; 
  k = [zeros(filterSize, filterSize), eye(filterSize)] ; 
   past1 = zeros (1,filterSize);
   past2 = zeros (1,filterSize) ;
   past3 = zeros (1,filterSize) ; 
   past4 = zeros (1,filterSize) ;
   
   
   error = 0 ;
  
 i= 0 ;
%  filterSize = 3*filterSize ;
while( i < numberOfIterations -1)
    
    if(i==0)
    x_n1 = x(i*filterSize+1:(i+1)*filterSize);
    T = filter (S2,1,w) ;
    
%     temporary = filter(S2,1,[past3, past4,x_n1]);
%     temporary = temporary(filterSize+1:end);
%     cancellerPrevious1 = filter(w,1,temporary);
    cancellerPrevious1 = filter (w,1,[past4,x_n1]);
    cancellerPrevious1 = cancellerPrevious1 (filterSize+1:end);
    
    x_n2 = x((i+1)*filterSize+1:(i+2)*filterSize);
%     temporary = filter(S2,1,[past4,x_n1,x_n2]);
%     temporary = temporary(filterSize+1:end);
%     cancellerPrevious2 = filter(w,1,temporary);
    cancellerPrevious2 = filter(w,1,[x_n1,x_n2]) ; 
    cancellerPrevious2 = cancellerPrevious2 (filterSize+1:end);
    
    past3 =x_n1 ;
    past4 =x_n2 ;

    
    
    
    else 
      %% Acquire first chunk   
      x_n1 = x(i*filterSize+1:(i+1)*filterSize);
      
      %Adaptation
      d_n_1 = d(i*filterSize+1:(i+1)*filterSize);
      T = filter(S2,1,w) ;
      c= filter(T,1,[past2, past3]) ;
      e_n_1 = d_n_1 - c(filterSize+1:end);
%       e_n_1 = d_n_1 - cancellerPrevious1 ; 

      [ W, w ] = FDAFOSMOneIteration( [past2, past3]', e_n_1', filterSize, P_k, g,k ,W0);
      error = [error, e_n_1* e_n_1'] ; 
      
      

      % Construct the cancellerSignal
%      temporary = filter(S2,1,[past3,past4,x_n1]);
%      temporary = temporary(filterSize+1:end);
%      cancellerPrevious1 = filter(w,1,temporary);
%        T = filter (S2,1,w) ;
%        cancellerPrevious1 = filter(T,1,[past4,x_n1]) ;
% 
%      cancellerPrevious1 = cancellerPrevious1 (filterSize+1:end);

      %% Acquire second chunk 
      
      x_n2 = x((i+1)*filterSize+1:(i+2)*filterSize);
      
      % Adaptation
      d_n_2 = d((i+1)*filterSize+1:(i+2)*filterSize);
%       e_n_2 = d_n_2 - cancellerPrevious2 ;
        T = filter(S2,1,w) ;

        c = filter(T,1,[past3, past4]);
        e_n_2 = d_n_2 - c(filterSize+1:end);
      
      error = [error, e_n_2* e_n_2'] ; 

      W0 = W ;
      [ W, w ] = FDAFOSMOneIteration( [past3, past4]', e_n_2', filterSize, P_k, g,k ,W0);
     

      
       % Construct the cancellerSignal
%      temporary = filter(S2,1,[past4,x_n1,x_n2]);
%      temporary = temporary(filterSize+1:end);
%      cancellerPrevious2 = filter(w,1,temporary);
%      T = filter(S2,1,w) ; 
%      cancellerPrevious2 = filter(T,1,[x_n1,x_n2]);
%      cancellerPrevious2 = cancellerPrevious2 (filterSize+1:end);

     W0= W ;
        
     past2 = past4 ;
     past3 = x_n1 ;
     past4 = x_n2 ;
        
    end
    
   i   = i+2  ;
    
    
end




blockSize = filterSize ; 
filterSize = 3*filterSize ;
convergenceThreshold = 10^-10 ;
%  d = d(delayPerSample+1:end);
[errorFDAF, MSerrorFDAF, timeOfConvergenceNLMSFDAF,timeOfComputationFDAF, w_FDAF]=FDAFOSM(x',d', filterSize,blockSize, numberOfIterations/2, convergenceThreshold)  ;

figure ;
semilogy(MSerrorFDAF) ; 

figure ;
subplot(3,1,1)
stem(filterToEstimate);
subplot(3,1,2)
stem(w_FDAF) ; 
subplot(3,1,3)
stem(w);


y = filter(w_FDAF,1,x);


figure ;
plot(y(1:length(d))) ; 
hold on 
plot(d) ;


corrData = xcorr(y,d) ; 

[~, index]= max(corrData) ;
index-length(d)