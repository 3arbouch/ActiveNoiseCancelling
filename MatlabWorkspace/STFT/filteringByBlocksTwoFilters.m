close all 
clear all 


data = load ('experiment31Filters(2048)');

f0= 800 ;
fs = 44100;
sineTime = 20 ;
filterSize  = 2048 ;


time = 0:(1/fs):sineTime - (1/fs) ;
sine= sin(2*pi*f0.*time);

% h = (1/filterSize)* ones(1,filterSize);
h= data.P;
h1 = data.S2 ; 

filtredSignal = filter(h1,1,sine);
filtredSignal = filter(h,1,filtredSignal);


numberOfIterations = floor(length(time)/filterSize);

filtredSignalBlocks = zeros(1,numberOfIterations*filterSize);

 k = [zeros(filterSize, filterSize), eye(filterSize)] ; 
% Filtering in Time Domain
tic
for i=0:numberOfIterations-1
    if (i==0)
        % First take two past blocks and the current block 
         portion = [zeros(1,filterSize),zeros(1,filterSize),sine(1:filterSize)];
         
         % Filter with the first filter 
         temp= filter(h1,1,portion) ; 
         
         portion = temp(filterSize+1:end) ; 
         
    elseif (i==1)
        portion = [zeros(1,filterSize),sine(1:filterSize),sine(filterSize+1:2*filterSize)];
         
         temp= filter(h1,1,portion) ; 
         
         portion = temp(filterSize+1:end) ; 

         
    else  
         portion = sine((i-2)*filterSize + 1 : (i+1)*filterSize);
         
         temp = filter(h1, 1,portion) ;
         portion = temp(filterSize+1:end) ; 


    end
        output = filter (h,1,portion);
         output= output(filterSize+1:end);
        filtredSignalBlocks((i)*filterSize + 1:(i+1)*filterSize) = output ;
  
    
end

computationTiemConvolutaion = toc 
figure; 
plot(filtredSignal);
hold on 
plot(filtredSignalBlocks);


