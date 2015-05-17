f0= 800 ;
fs = 44100;
sineTime = 15 ;
sine= 0.5*sin(2*pi*f0/fs.*(1:sineTime*fs)');


sound(sine,fs);