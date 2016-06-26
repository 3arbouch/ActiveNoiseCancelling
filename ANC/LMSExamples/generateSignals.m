function [ x,d,h ] = generateSignals( signalLength, originalFilterSize, optionReferenceSignal, optionFilter, colored )
% This function is used to generate signals for testing
% The reference signal is either a white noise or a sine tone accordig the
% optionalReferenceSignal:
%   1) 'wn' for white noise
%   2)  'st' for sine tone
% If the option colored is active , the generated reference signal will be a
% colored white gaussian Noise.
% colored equal:
%   1) 1 generate a colored noise
%   2) 0 generate a white noise
% The filter is either a low pass filter or a moving average with number of taps originalFilterSize
% The optionFilter is equal:
%   1) 'lp' to generate a low pass filter
%   2) 'ma' to generate a moving average filter
% return the reference Signal x, the filtred version of x d and the filter

if(strcmp(optionFilter,'lp'))
    [h,~,~]=fircband(originalFilterSize,[0 0.4 0.5 1], [1 1 0 0], [1 0.2],{'w' 'c'});
elseif(strcmp(optionFilter,'ma'))
     h = ones(originalFilterSize,1)/originalFilterSize ;
end

if(strcmp(optionReferenceSignal,'wn'))
     x = 0.25*randn(signalLength,1);
     
elseif (strcmp(optionReferenceSignal,'st'))
    f0 = 800; % frequencz of the sine 
    fs = 8000; % sampling freqeuncy 
    sineTime = signalLength/fs ; 
    x= 0.5*sin(2*pi*f0/fs.*(1:sineTime*fs)');

end

if(colored && strcmp(optionReferenceSignal,'wn'))
    A = [1 -0.9];
    x = filter(1,A,x);
end



x = x-mean(x) ;
d = filter(h,1,x);



end

