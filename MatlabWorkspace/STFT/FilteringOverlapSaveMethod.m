function [out] = FilteringOverlapSaveMethod
%% This script shows how filtering by blocks operates using the overlap save method  
close all 

f0= 800 ;
fs = 44100;
sineTime = 2 ;
filterSize  = 512 ;
time = 0:(1/fs):sineTime - (1/fs) ;

% Define the signal 
x= sin(2*pi*f0.*time);

%  x = 0.25*ones(1,400) ; 
  x = randn(1,10000) ; 
% Define a moving average filter 
h = (1/filterSize)*ones(1,filterSize) ; 

% h = 0.5*ones(1,filterSize) ; 
y = conv(h,x) ; 

% Vector that will contain the output of the filtred blocks 
frameSize = 128 ; 

lengthOutput = max([length(x)+length(h)-1,length(x),length(h)]) ; 
output = zeros(1, lengthOutput) ; 
size(output)
numberOfOfteration = floor(length(x)/frameSize) ; 


for i =0: numberOfOfteration -1
    %Extract the block of the signal 
    x_n = x(i*frameSize+1:(i+1)*frameSize) ; 
    y_n = conv(h,x_n) ; 
    
    size(y_n)
    %update the output vector 
    output =  updateOutput (output,y_n,i*frameSize + 1) ; 
    
    
end
    

figure ; 
plot(output) ; 
hold on ; 
plot(y)

figure ;
plot(output-y) ; 

datacorr=xcorr(output,y) ; 
[value, index]= max(datacorr) ; 

shift = index -length(y)   

end


%% update the output vector 
function [outputNew] = updateOutput( output , y, index )

output(index:index+ length(y) -1 ) = output(index:index+ length(y) -1 ) + y ;  
 outputNew = output ; 
 
end




