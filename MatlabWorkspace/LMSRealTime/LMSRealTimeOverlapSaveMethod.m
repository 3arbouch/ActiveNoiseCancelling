function [] = LMSRealTimeOverlapSaveMethod
%% This script shows how filtering by blocks operates using the overlap save method  
close all 

data = load ('experiment90Filters(2048, 32khz)');


f0= 200 ;
fs = 32000;
sineTime = 1 ;
filterSize  = 2048 ;
time = 0:(1/fs):sineTime - (1/fs) ;

% Define the signal 
%x= sin(2*pi*f0.*time);
 x = randn(1,fs*sineTime) ; 
% mean(x(1:filterSize));
 
%  x= 0.5*ones(1,fs*sineTime) ;
% Define a moving average filter 
h = (1/filterSize)*ones(1,filterSize) ; 
h = data.P ; 
% [h,~,~]=fircband(filterSize,[0 0.4 0.5 1], [1 1 0 0], [1 0.2],{'w' 'c'});
% h = ones(1,filterSize) ; 



y = conv(h,x) ; 


tfplot(x,fs,'x','Original Signal') ; 


tfplot(y,fs,'y','Filtred Signal') ; 

 tfplot(h,fs,'h','Original Filter') ; 

% Simulate time of flight per samples
samplesTimeOfFlight = 250 ;
sampleDelay = samplesTimeOfFlight/2 ; 

% Take the frame Size as the nearest power of two to sample Delay in the
% lower band

power = (nextpow2(sampleDelay) -1) ; 
frameSize = 2^power    

numberOfFrameBeforeAdaptation  = filterSize/frameSize  ; 

lengthOutput = max([length(x)+length(h)-1,length(x),length(h)]) ; 

% Vector that will contain the output of the filtred blocks
output = zeros(1, lengthOutput) ; 

% Define the past and actual vectors that will be used for adaptation 

pastRefMic = zeros(1, filterSize) ; 
actualRefMic = zeros(1, filterSize) ; 

%Define the parameters of FDAF 
  W0 = zeros(2*filterSize, 1) ;
  w = zeros(filterSize,1) ;
  P_k= ones(2*filterSize,1)  ; 
  g=[eye(filterSize),zeros(filterSize,filterSize); zeros(filterSize,filterSize), zeros(filterSize,filterSize)] ; 
  k = [zeros(filterSize, filterSize), eye(filterSize)] ; 


numberOfOfteration = floor(length(x)/frameSize) ; 

% Used to determine when we reach the number of blocks needed to perform adaptation 
j =0 ; 

MSerror = 0 ;
%% Loop for adaptation 
for i =0: numberOfOfteration -1
    %Extract the block of the signal 
    x_n = x(i*frameSize+1:(i+1)*frameSize) ;
    actualRefMic(j*frameSize+1:(j+1)*frameSize) = x_n ; 
    
    % Perform the convolution
    y_n = conv(w,x_n) ; 
    %update the output vector 
    output =  updateOutput (output,y_n',i*frameSize+1) ; 
    
    j = j+1 ; 
    
    % Test if we reach the number of blocks nedded to perform the adaptation
    if(j == numberOfFrameBeforeAdaptation)
        % construct the errorSignal 
        j 
        i+1
      y_adpt = output((i+1-j)*frameSize+1:(i+1)*frameSize)  ; 
      e_adpt = y((i-j+1)*frameSize+1:(i+1)*frameSize) + y_adpt ;
       
     
      
     MSerror = [MSerror,e_adpt] ; 
      
      reference = [pastRefMic , actualRefMic] ; 
      
      [ W, w ] = FDAFOSMOneIteration( reference', e_adpt', filterSize, P_k, g,k ,W0);
      
      pastRefMic = actualRefMic ; 

      W0= W ;  
      j = 0 ; 
    end
    
    
end
    

figure ; 
semilogy(abs(MSerror)) ; 
title('MSE error') ; 


figure ; 
plot(output) ; 
hold on ; 
plot(y)
title ('Output of the overlap save method VS filtred signal ') ; 

figure ; 
subplot(1,2,1);
stem(w) ; 
title('Estimated Filter') ; 
subplot(1,2,2) ;
stem(h) ; 
title ('original filter') ; 


end


%% update the output vector 
function [outputNew] = updateOutput( output , y, index )
output(index: index+ length(y) -1 ) = output(index:index+ length(y) -1 ) - y ;  
 outputNew = output ; 
 
end




