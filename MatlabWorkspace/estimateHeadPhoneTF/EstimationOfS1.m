clear all
close all 
endTime=50;
f0= 6000 ; 
frameSize=1024;

%% Generate signals
%noise= dsp.SignalSource(0.25*randn(endTime*samplingRate, 1), frameSize);

song=dsp.AudioFileReader('song.mp3');
sineSweep = dsp.AudioFileReader('Sweep_10Hz_15000Hz_-3dBFS_10s.wav')   ; 
samplingRate= 8000 ; 
sine= dsp.SignalSource(0.5*sin(2*pi*f0/samplingRate.*(1:endTime*samplingRate)'), frameSize);

noise=dsp.SignalSource(0.1*randn(endTime*samplingRate, 1), frameSize);
%noise = sine ; 
%noise = song  ; 
%% Setup MOTU
H=dsp.AudioPlayer('DeviceName','MOTU 896','SampleRate', samplingRate, 'ChannelMapping', 1:3);

AudioInput = dsp.AudioRecorder(...
            'DeviceName', 'MOTU 896', ...
            'SampleRate', samplingRate, ...
            'NumChannels', 3,...
            'OutputDataType','double',...
            'QueueDuration', 2,...
            'SamplesPerFrame', frameSize);


null_signal=zeros(frameSize,1);
recordings=zeros(1,3);
noSound = zeros(frameSize,1) ; 
noiseData = 0 ; 

 

 
load('cropped_S1_hat') ; 
tic
while (~isDone(noise))
    audio1 = step(noise);
    
    %Speacker, headphoneLeft, headPhoneRight 
    step(H, [audio1(:,1),noSound, noSound]);
    
    
    noiseData = [noiseData;  audio1(:,1)] ; 
    multichannelAudioFrame = step(AudioInput);
    recordings=[recordings;multichannelAudioFrame];
    
  

    
    
end
toc

%% Release all the objects
release(noise); % release the input file
release(AudioInput);
release(H);  % release the audio output device
D = recordings(:,1) ;


%% Estimate the filter : LMS Algo


startPoint = 0 ; 
maxIter = samplingRate*18 ;
%maxIter = length(noiseData)-2 ; 
 
filterSize = 2000 ;   
i = samplingRate+5000; 

D= D(2+startPoint:end) ; 


noiseData_cropped= noiseData(2+startPoint:end) ; 
x = zeros(filterSize,1) ;
x = [x;noiseData_cropped] ; 


S1_hat = zeros(filterSize,1) ; 
 

epsilon = 10^-8 ; 
error =0 ;  
i=3 ; 
while(i< 15)
             
%       
%              x_n = x(i+filterSize:-1:i+1) ; 
%              e = D(i)- x_n'*S1_hat ; 
%              mu = 1/(x_n'*x_n) ; 
%              S1_hat = S1_hat + mu.* x_n*e ; 
%              error = [error; e] ;
%             
%            
%              
%              
%              if(abs(e)< epsilon)
%                  break ; 
%              end
             
             [x_matrix, mu ]=formX(filterSize, i, x) ;
           
             e = D(i*filterSize:(i+1)*filterSize-1)-x_matrix*S1_hat  ; 
             error = [error ; e] ; 
             S1_hat = S1_hat + 2*mu* x_matrix'*e ;
             
i  = i + 1
end

figure ; 
plot(error) ; 
title('error')  ; 

figure ; 
stem(S1_hat, 'b')
title ('Estimated Impulse Response')


figure ; 
plot(abs(fft(S1_hat))) ; 
title('Magnitude of the Fourier the estimated impulse response ')

%% Plot an STFT of all mics
% Reference mic
[STFT1,F1,T1] = spectrogram(recordings(:,1),1024,512,1024,samplingRate);

figure ; 
imagesc(T1,F1,log10(abs(STFT1)))

figure  ; 
plot(recordings(:,1)) ; 

figure ; 
plot(error) ; 
title ('Error')
%% Estimation of the Right transfer function of S1(z)
 
% 
% XFFT = fft(noiseData_cropped) ; 
%  
% 
% 
% figure ; 
% plot(D) ; 
% title ('Recorded data cropped')
% DFFT = fft(D)  ; 
% 
% S1FFT = DFFT./(XFFT) ; 
% 
% 
% figure ; 
% 
% plot(abs(S1FFT)) ; 
% title ('Fourrier Transform of the estimated Filter ')
% 
% 
% 
% figure  ; 
% S1_hat  = ifft(S1FFT) ; 
% stem(S1_hat) ; 
% title ('Estimated impulse response')
% 
% 
% 
% %% Plot the error
% f = filter(S1_hat, 1, noiseData) ; 
% figure ; 
% plot(f-recordings(:,1)) ; 
% title('error') ; 