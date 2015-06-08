close all





%%  Generate signals of interest
data = load('experiment50.mat');
filters = load ('experiment50Filters(1024)');

f0= 10000 ;
fs = 16000;
sineTime = 2 ;
filterSize  = 2048 ;


time = 0:(1/fs):sineTime - (1/fs) ;
sine= sin(2*pi*f0.*time);
noise = 0.5*randn(1,sineTime*fs) ; 




x = data.referenceMicPrime ;
%  x = x-mean(x) ;


% h = (1/filterSize)* ones(1,filterSize);
% Filter we want to estimate 
h= filters.P;

%Add additional delay to h and to the filtredSignal 
%   delay = zeros(2048,1) ; 
%   h = [delay;h];
 [~,index]= max(h) ; 
%  
 filterSize = floor(index/2) 
% h = h(index:end) ; 
% 
% % Add additional delay to the filter P(z): h
% delay = filterSize*2 ; 
% 
% h = [zeros(delay+10,1);h] ; 


% filterSize = filterSize/2 ; 
 %h = (1/filterSize)* ones(1,filterSize);
% This filter is estimated offline
h1 = filters.S2 ; 


 d = filter(h1,1,x);
% filtredSignal = filter(h,1,d);

%  filtredSignal = [delay',data.errorMicPrime ]; 

filtredSignal = data.errorMicPrime ;

numberOfIterations = floor(length(time)/filterSize);

filtredSignalBlocks = zeros(1,numberOfIterations*filterSize);
%% Parameters of the FDAFOSM algorithm
   W0 = zeros(2*filterSize, 1) ;
   w = zeros(filterSize,1) ;
  P_k= ones(2*filterSize,1)  ; 
  g=[eye(filterSize),zeros(filterSize,filterSize); zeros(filterSize,filterSize), zeros(filterSize,filterSize)] ; 
  k = [zeros(filterSize, filterSize), eye(filterSize)] ; 

  numberOfIterations = floor(length(x)/filterSize) - 1;
  completeMSerror = 0 ;
  
  
  for i =0:numberOfIterations-3
      
      if(i==0)
          % First take two past blocks and the current block 
         portion = [zeros(1,filterSize),zeros(1,filterSize),x(1:filterSize)];
         
         % Filter with the first known filter  
         temp= filter(h1,1,portion) ; 
         
         x_n = temp(filterSize+1:end) ; 
         
      elseif(i ==1)     
          % First take two past blocks and the current block 
         portion = [zeros(1,filterSize),x(1:filterSize),x(filterSize+1:2*filterSize)];
         
         % Filter with the first known filter  
         temp= filter(h1,1,portion) ; 
         
         x_n = temp(filterSize+1:end) ; 
         
      else 
         portion = x((i-2)*filterSize+1:1:(i+1)*filterSize) ;
          temp= filter(h1,1,portion) ; 
         
         x_n = temp(filterSize+1:end) ; 
      end
      
      
     
      d_k = filtredSignal((i+2)*filterSize+1:(i+3)*filterSize) ;
    
      y_k = filter(w,1,x_n) ; 
     
      e_k = d_k - y_k(filterSize+1:end) ;
      
     
    
      [ W, w ] = FDAFOSMOneIteration( x_n', e_k', filterSize, P_k, g,k ,W0);
     
   
        completeMSerror = [completeMSerror,e_k*e_k'] ;
      W0 = W ;
  end


%   Parform the adapttaion with the hole signal
filterSize =1024;
blockSize = filterSize ; 
convergenceThreshold = 10^-20 ;
filtredSignal = filtredSignal(index:end);
[errorFDAF, MSerrorFDAF, timeOfConvergenceNLMSFDAF,timeOfComputationFDAF, w_FDAF]=FDAFOSM(d',filtredSignal', filterSize,blockSize, numberOfIterations/8, convergenceThreshold)  ;




  figure ;
  
  semilogy(completeMSerror);
  title ('Error with adaptation block by block')
  
  
  figure ;
  semilogy(MSerrorFDAF); 
  title ('Error with adaptation with hole signal ')

 

tfplot(h,fs,'h',' original filter') ; 

tfplot(w,fs,'w','block adaptation');

tfplot(w_FDAF,fs,'w_FDAF','whole signal adaptation');



figure ; 

plot(filtredSignal); 
hold on 
plot(filter(h,1,d)) ;


figure ;

plot(filtredSignal)
hold on 
temp =filter(w_FDAF,1,d) ; 
plot(temp(index:end));


figure ; 
plot(filtredSignal)
hold on 
temp = filter(w,1,d) ; 
plot(temp(index:end));








