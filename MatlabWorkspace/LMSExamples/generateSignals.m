function [ x,d,h ] = generateSignals( signalLength, originalFilterSize )
% This functio is used to generate signals for testing
% The reference signal is a white noise 
% The filter is a low pass filter with number of taps originalFilterSize
% return the reference Signal x, the filtred version of x d and the filter
% h
[h,~,~]=fircband(originalFilterSize,[0 0.4 0.5 1], [1 1 0 0], [1 0.2],{'w' 'c'});
 x = 0.25*randn(signalLength,1);
 x = x-mean(x) ; 
 d = filter(h,1,x);



end

