% This example illustrates two ways of seeing filtering: This is useful for
% the real time filtering when updating the buffer of the loudspeacker 




f0= 800 ;
fs = 44100;
sineTime = 20 ;
filterSize  = 1024 ;

% design a moving average filter 
h = (1/filterSize)* ones(1,filterSize) ;

x = ones(1,100) ;

y = filter(x,1,h) ;
plot(y);

% time = 0:(1/fs):sineTime - (1/fs) ;
% sine= sin(2*pi*f0.*time);
% 
% 
% filtredSignal = filter(h,1,sine);
% 
% 
% numberOfIterations = floor(length(time)/filterSize);
% 
% filtredSignalBlocks = zeros(1,numberOfIterations*filterSize);
% 
%  k = [zeros(filterSize, filterSize), eye(filterSize)] ; 
% % Filtering in Time Domain
% tic
% for i=0:numberOfIterations-1
%     if (i==0)
%          portion = [zeros(1,filterSize),sine(1:filterSize)];
%          
%     else
%          portion = sine((i-1)*filterSize + 1 : (i+1)*filterSize);
%            
%     end
%         output = filter (h,1,portion);
%          output= output(filterSize+1:end);
%         filtredSignalBlocks((i)*filterSize + 1:(i+1)*filterSize) = output ;
%   
%     
% end
% 
% 
% plot(filtredSignalBlocks)
% hold on 
% plot(filtredSignal)
