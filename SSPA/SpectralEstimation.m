%First, estimate the spectrum with periodogram
% Very high variation

rrSignal = rrSignal - mean(rrSignal);
pr = periodogram(rrSignal);
figure;
plot(pr);

% pwelch: let's reduce the variation of the periodogram 
% take blocks of N segmets of the segment: for each segement, compute the
% periodogra -Z Take the average of the periodograms 
% NP: Instead of using a sqaure window: use a cosine or a guassian window :
% reduce the power f certain samples: end ones 
% Use overlapping sections: 50 % overlapping 

prwelch= pwelch(rrSignal);
figure ;
plot (prwelch) ;

% Model the signal as autoregressive 

Pxx = pyulear(rrSignal,8);
figure ;
plot(Pxx);