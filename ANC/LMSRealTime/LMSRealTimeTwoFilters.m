close all





%%  Generate signals of interest 
data = load ('experiment31Filters(2048)');

f0= 15000 ;
fs = 44100;
sineTime = 2 ;
filterSize  = 2048 ;


time = 0:(1/fs):sineTime - (1/fs) ;
sine= sin(2*pi*f0.*time);
noise = 0.5*randn(1,sineTime*fs) ; 




x = noise;
%  x = x-mean(x) ;


% h = (1/filterSize)* ones(1,filterSize);
% Filter we want to estimate 
h= data.P;
 %h = (1/filterSize)* ones(1,filterSize);
% This filter is estimated offline
h1 = data.S2 ; 
d = filter(h1,1,x);
filtredSignal = filter(h,1,d);


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
  for i =0:numberOfIterations-1
      
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
      
      
     
      d_k = filtredSignal(i*filterSize+1:(i+1)*filterSize) ;
    
      y_k = filter(w,1,x_n) ; 
     
      e_k = d_k - y_k(filterSize+1:end) ;
      
     
    
      [ W, w ] = FDAFOSMOneIteration( x_n', e_k', filterSize, P_k, g,k ,W0);
     
   
        completeMSerror = [completeMSerror,e_k*e_k'] ;
      W0 = W ;
  end


%   Parform the adapttaion with the hole signal
blockSize = filterSize ; 
convergenceThreshold = 10^-20 ;
[errorFDAF, MSerrorFDAF, timeOfConvergenceNLMSFDAF,timeOfComputationFDAF, w_FDAF]=FDAFOSM(d',filtredSignal', filterSize,blockSize, numberOfIterations, convergenceThreshold)  ;




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
plot(filter(w_FDAF,1,d));


figure ; 
plot(filtredSignal)
hold on 
plot(filter(w,1,d));








