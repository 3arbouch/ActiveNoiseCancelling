close all 
clear all 
dataExperiment = load('ExperimentWithMusicData30Sine.mat');
dataFilters = load('experiment30Filters(4096).mat');
fs = 44100 ;
% Cancel the echo from the headphones 
portion = dataExperiment.portion ;
cancellerSignal = filter(dataFilters.S1,1,portion);
cancelledSignal = dataExperiment.completeSignalRefMic - cancellerSignal ;

% Cancel the sine inside the headphones

temporary = filter (dataFilters.S2,1,cancelledSignal);
sineCancellerSignal = filter(dataFilters.P,1,temporary);
denoisedSignal = dataExperiment.completeSignalErrorMic-sineCancellerSignal ;

figure; 
stem(dataFilters.P);
title('Headphones Filter: P(z)')

figure; 
stem(dataFilters.S1);
title('Echo canceller Filter: S1(z)')

figure; 
stem(dataFilters.S2);
title('Intern headphone filter S2(z)')



figure ;
plot(dataExperiment.completeSignalErrorMic);
hold on 
plot(sineCancellerSignal)

       
    figure;
    subplot(2, 1, 1);
    [STFT1,F1,T1] = spectrogram(dataExperiment.completeSignalErrorMic,1024,512,1024,fs);
    imagesc(T1,F1,log10(abs(STFT1)))
    title('Listened Music with noise')
    subplot(2,1,2)
    [STFT3,F3,T3] = spectrogram(denoisedSignal,1024,512,1024,fs);
    imagesc(T3,F3,log10(abs(STFT3)))
    title('Cancelled Signal: Music witout noise')
