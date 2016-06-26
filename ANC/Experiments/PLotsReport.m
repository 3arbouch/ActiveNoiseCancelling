close all 
clear all 
load('RealTimeCancellingExperimentReportV3') ; 


fs = 8000;

figure;
   
    subplot(2,1,1)
    
    [STFT1,F1,T1] = spectrogram(completeReferenceMic(1:20*fs),1024,512,1024,fs);
    imagesc(T1,F1,log10(abs(STFT1)))
    title('Signal captured at reference microphone')
    set(gca, 'FontSize', 20) ; 
    
        subplot(2,1,2)

    [STFT1,F1,T1] = spectrogram(completeErrorMic(1:20*fs),1024,512,1024,fs);
    imagesc(T1,F1,log10(abs(STFT1)))
    title('Signal captured at error microphone')
    set(gca, 'FontSize', 20) ; 


% figure; 
% stem(P); 
% xlabel('samples');
% ylabel('P[n]');
% set(gca, 'FontSize', 20) ; 
% title('Estimated headphone filter')
% 
% 
% figure; 
% stem(S1); 
% xlabel('samples');
% ylabel('S1[n]');
% set(gca, 'FontSize', 20) ; 
% title('Estimated headphone reference microphone filter')
% 
% 
% 
% figure; 
% stem(S2); 
% xlabel('samples');
% ylabel('S2[n]');
% set(gca, 'FontSize', 20) ; 
% title('Estimated intern headphone filter')
