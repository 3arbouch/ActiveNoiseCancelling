% TFPLOT Time and frequency plot
%    TFPLOT(S, FS, NAME, TITLE) displays a figure window with two subplots.
%    Above, the signal S is plotted in time domain; below, the signal is plotted
%    in frequency domain. NAME is the "name" of the signal, e.g., if NAME is
%    's', then the labels on the y-axes will be 's(t)' and '|s_F(f)|',
%    respectively.  TITLE is the title that will appear above the two plots. 

function [] = tfplot(s, fs, name, plottitle)

% Note: Since TITLE is the name of a built-in Matlab function, 
% you cannot use it as the name of a function argument. 
% In the following code, we called the argument PLOTTITLE instead 

% Compute the time and frequency scales
t = linspace(0, (length(s)-1) / fs, length(s));

NFFT = 2^nextpow2(length(s)); 
f = linspace(-fs/2, fs/2-fs/NFFT, NFFT);

% compute the FFT
s_f = fft(s,NFFT);

figure;
% First plot: time
subplot(3,1,1); plot(t, s);
xlabel('t [s]'); ylabel(sprintf('%s(t)', name));
title(plottitle);

%Second plot: Sample plot
subplot(3,1,2);
stem(s); 
xlabel('samples'); ylabel(sprintf('%s[n])', name));

% Third plot: frequency
% We use fftshift to move the coefficients for negative frequencies to the left
subplot(3,1,3); plot(f, fftshift(abs(s_f)));
xlabel('f [Hz]'); ylabel(sprintf('|%s_F(f)|', name));
