close all 
clear all 


data = load ('experiment31Filters(2048)');

f0= 800 ;
fs = 44100;
sineTime = 20 ;
filterSize  = 1024 ;


time = 0:(1/fs):sineTime - (1/fs) ;
sine= sin(2*pi*f0.*time);

% h = (1/filterSize)* ones(1,filterSize);
h= data.P;

filtredSignal = filter(h,1,sine);


numberOfIterations = floor(length(time)/filterSize);

filtredSignalBlocks = zeros(1,numberOfIterations*filterSize);

 k = [zeros(filterSize, filterSize), eye(filterSize)] ; 
% Filtering in Time Domain
tic
for i=0:numberOfIterations-1
    if (i==0)
         portion = [zeros(1,filterSize),sine(1:filterSize)];
         
    else
         portion = sine((i-1)*filterSize + 1 : (i+1)*filterSize);
           
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


% Filtering in frequency domain 
 k = [zeros(filterSize, filterSize), eye(filterSize)] ; 
 H = fft([h',zeros(1,filterSize)]);

tic
for i=0:numberOfIterations-1
    if (i==0)
         portion = [zeros(1,filterSize),sine(1:filterSize)];
         
    else
         portion = sine((i-1)*filterSize + 1 : (i+1)*filterSize);
           
    end
       
        OUTPUT = fft(portion);
        OUTPUT = OUTPUT.*H;
        output = k*ifft(OUTPUT') ; 
        filtredSignalBlocks((i)*filterSize + 1:(i+1)*filterSize) = output ;
  
    
end

computationTimeFourrier = toc 

figure; 
plot(filtredSignal);
hold on 
plot(filtredSignalBlocks);