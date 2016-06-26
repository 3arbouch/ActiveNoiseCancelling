close all 
clear all 
dataExperiment = load('ExperimentWithMusicData23.mat');
dataFilters = load('experiment23Filters.mat');
fs = 44100 ;
% Cancel the echo from the headphones 
portion = dataExperiment.portion ;
cancellerSignal = filter(dataFilters.S1,1,portion);
cancelledSignal = dataExperiment.completeSignalRefMic - cancellerSignal ;

% Cancel the sine inside the headphones

temporary = filter (dataFilters.S2,1,cancelledSignal);
sineCancellerSignal = filter(dataFilters.P,1,temporary);
denoisedSignal = dataExperiment.completeSignalErrorMic-sineCancellerSignal ;

% Sound with music offline

%sound(dataExperiment.completeSignalErrorMic, fs) ; 
% sound(denoisedSignal, fs) ; 


%%  Sine Sound online 
data = load('RealTimeCancellingExperimentReportV3') ;

sound(data.completeErrorMic(1:15*data.fs),data.fs);
