    close all ; 
    clear all ; 
    
    % Load the data 
    load('ExperimentMovingZ1');
    if withNoise
       fprintf('SNR level is %s \n', SNR)   
    end
    
    if moving
      fprintf('The microphone is moving \n')   

    else 
       fprintf('The microphone is not moving \n')   

    end

    
 % Get the signals 
 referenceSignal1 = recBuffer(:,1) ; 
 referenceSignal2 =  recBuffer(:,2) ;
 movingMicSignal = recBuffer (:,3) ; 
        

 
 % Get the distances
 distanceRefSpeck = [distances.distanceRef1Speack1, distances.distanceRef2Speack2] ; 
 distancesMovSpeck = [distances.distanceMovSpeck1; distances.distanceMovSpeck2 ] ; 
    

 % Set up parameters for all Algos
    
    fs = 48000 ; 
    speedOfSound = 340.29 ; 
    filterSize = 1200 ; 
    convergenceThreshold = 10^-12 ; 

    sound(movingMicSignal(1:10*fs),fs);
% Set up parameters for some algos 
  
  %LMS
  mu1 = 10^-3 ; 
  mu2 = 10^-3 ; 
    
  lambda = 10^-3 ;
  
  
  
 mu_k1 = filterSize ; 
 mu_k2 = filterSize ; 
 interval1 = filterSize*2 ; 
 interval2 = filterSize*2 ; 
 upperLimit = 10^-3 ; 
 memoryInSeconds = 2 ; 

 noiseSamples = 100 ; 
%  mu_k1 = ((distancesMovSpeck(1) - distanceRefSpeck(1))/speedOfSound)*fs + ceil(abs(noiseSamples*randn(1,1))) ;
%  mu_k2 = ((distancesMovSpeck(2) - distanceRefSpeck(2))/speedOfSound)*fs + ceil(abs(noiseSamples*randn(1,1))) ;
%   
%  interval1 = 100 ; 
%  interval2 = 100 ; 
%  
 
weighted = 0 ; 

% Index of the algo to run: 
% 1 --> LMSSterio
    AlgosToRun = [13] ; 
    
    
  
for algo = AlgosToRun
    
   switch (algo)
        
       case 1
           [ error, MSerror, ~,~,  w1, w2, filters1, filters2] = LMSSterio( referenceSignal1,referenceSignal2, movingMicSignal, fs, filterSize, convergenceThreshold, mu1, mu2 );
         
       case 2 
           [ error, MSerror, timeOfConvergence,timeOfComputation,  w1, w2, filters1, filters2] = NLMSSterio( referenceSignal1,referenceSignal2, movingMicSignal, fs, filterSize, convergenceThreshold) ; 
    
           
       case 3 
           [ error, MSerror, timeOfConvergence,timeOfComputation,  w1, w2, filters1, filters2] = LMSSterioL0( referenceSignal1,referenceSignal2, movingMicSignal, fs, filterSize, convergenceThreshold, mu1, mu2, lambda );
           
           
       case 4
           [ error, MSerror, timeOfConvergence,timeOfComputation,  w1, w2, filters1, filters2] = NLMSSterioL0( referenceSignal1,referenceSignal2, movingMicSignal, fs, filterSize, convergenceThreshold, lambda ) ; 
       
       case 5
           
           [ error, MSerror, timeOfConvergence,timeOfComputation,  w1, w2, filters1, filters2] = LMSSterioL1( referenceSignal1,referenceSignal2, movingMicSignal, fs, filterSize, convergenceThreshold, mu1, mu2, lambda ) ; 
       case 6
           
           [ error, MSerror, timeOfConvergence,timeOfComputation,  w1, w2, filters1, filters2] = NLMSSterioL1( referenceSignal1,referenceSignal2, movingMicSignal, fs, filterSize, convergenceThreshold, lambda ) ;     
           
       case 7 
           
           [ error, MSerror, timeOfConvergence,timeOfComputation,  w1, w2, filters1, filters2] = FDAFSterio( referenceSignal1, referenceSignal2, movingMicSignal, filterSize,  convergenceThreshold ) ; 
           
           
       case 8 

           [ error, MSerror, timeOfConvergence,timeOfComputation,  w1, w2, filters1, filters2] = FDAFSterioL0( referenceSignal1, referenceSignal2, movingMicSignal, filterSize,  convergenceThreshold, lambda ) ; 
           
       case 9
           
           [ error, MSerror, timeOfConvergence,timeOfComputation,  w1, w2, filters1, filters2] = FDAFSterioL1( referenceSignal1, referenceSignal2, movingMicSignal, filterSize,  convergenceThreshold, lambda ) ; 

       case 10 
        
           [ error, MSerror, timeOfConvergence,timeOfComputation,  w1, w2, filters1, filters2, weights1, weights2, variances1, variances2] = NLMSSterioL0Weighted( referenceSignal1,referenceSignal2, movingMicSignal, fs, filterSize, convergenceThreshold, mu_k1, mu_k2, interval1, interval2, upperLimit, memoryInSeconds ) ; 
            weighted = 1 ; 
           
       case 11 
        
           [ error, MSerror, timeOfConvergence,timeOfComputation,  w1, w2, filters1, filters2, weights1, weights2, variances1, variances2] = NLMSSterioL1Weighted( referenceSignal1,referenceSignal2, movingMicSignal, fs, filterSize, convergenceThreshold, mu_k1, mu_k2, interval1, interval2, upperLimit, memoryInSeconds ) ;
            weighted = 1 ; 
       case 12 
           
           [ error, MSerror, timeOfConvergence,timeOfComputation,  w1, w2, filters1, filters2,weights1, weights2, variances1, variances2] = FDAFSterioL0Weighted( referenceSignal1, referenceSignal2, movingMicSignal, filterSize,  convergenceThreshold, mu_k1, mu_k2, interval1, interval2, upperLimit, memoryInSeconds );
            weighted = 1 ; 
           
        case 13 
           
           [ error, MSerror, timeOfConvergence,timeOfComputation,  w1, w2, filters1, filters2,weights1, weights2, variances1, variances2] = FDAFSterioL1Weighted( referenceSignal1, referenceSignal2, movingMicSignal, filterSize,  convergenceThreshold, mu_k1, mu_k2, interval1, interval2, upperLimit, memoryInSeconds );
           weighted = 1 ; 
   end
    
end

         
%% Plots 
close all
    distanceMovSpeck1 = [3.68 3.19 2.82 2.68 3.38] ; 
    distanceMovSpeck2 = [2.31 1.92 1.92 2.37 3.17] ; 

    numberRefPoints = length(distanceMovSpeck1) ; 

    a=5*ones(1,1000);
    sRef = 50*ones(1,numberRefPoints) ; 
    siz =[200,200,100,sRef,a];
    colorCirle = zeros(1000,3);
    colorCirle(:,2)=1;
    colorRefSpeck = [1 0.6471 0] ; 
    colorMovingMic = [1 0 0] ; 
    colorRefPoints = zeros(numberRefPoints,3) ; 
    colorRefPoints(:,3) = 1 ; 
    c= [colorRefSpeck;colorRefSpeck; colorMovingMic ;colorRefPoints;colorCirle];



% calculate the time axis    
  t = (filterSize*(0:size(filters1,1)-1))/fs ;

% calculate the distances
  distance1 = ones(1,length(t)) ; 
  distance2 = ones(1,length(t)) ;
    if(moving==1)

        for i=1:length(timesIntervals)-1
            lowerIndex = max(find(t<=timesIntervals(i)));
            maxIndex = max(find(t<=timesIntervals(i+1)));
            distance1(lowerIndex:maxIndex) = distances.distanceMovSpeck1(i);
            distance2(lowerIndex:maxIndex) = distances.distanceMovSpeck2(i);

        end


    else 
        distance1 = distances.distanceMovSpeck1(1)*distance1 ; 
        distance2 = distances.distanceMovSpeck2(1)*distance2 ; 
    end       

    
    % Mean sqaure error
    figure ;
    semilogy(t,MSerror);

    % Estimated filters
    figure ;
    plot(w1) ;
    hold on 
    plot(w2);



     
    % Calculate the distances from the delays
     numberOfFrames = size(filters1,1) ;
     numberOfSamplesToIgnore = 50 ;
     [~, delay1] = max(abs(filters1(:,numberOfSamplesToIgnore:filterSize-numberOfSamplesToIgnore)),[],2) ; 
     [~, delay2] = max(abs(filters2(:,numberOfSamplesToIgnore:filterSize-numberOfSamplesToIgnore)),[],2) ; 

     delay1 = delay1 +numberOfSamplesToIgnore ;
     delay2 = delay2 +numberOfSamplesToIgnore ;

     interval =10; 
     convergedDelay1 = CalculateConvergedDelays( delay1, interval ) ; 
     convergedDelay2 = CalculateConvergedDelays( delay2, interval ) ; 

     convergedDistance1 = (convergedDelay1/fs)*speedOfSound + distanceRefSpeck(1) ;
     convergedDistance2 = (convergedDelay2/fs)*speedOfSound + distanceRefSpeck(2) ;

    % Plot the calculated and the converged distances
     figure ;
     plot(t, convergedDistance1') ;
     hold on
     plot(t, distance1') ;

     rms1=sqrt(sum((convergedDistance1(:)-distance1(:)).^2)/numel(convergedDistance1))

     figure;
     plot(t,convergedDistance2');
     hold on
     plot(t, distance2') ;

     rms2=sqrt(sum((convergedDistance2(:)-distance2(:)).^2)/numel(convergedDistance2))


     
      if(~weighted)
      f = figure(5) ; 
       s = uicontrol(f,'Style','slider',...
            'Min',1,'Max',numberOfFrames,'Value',1,...
            'SliderStep',[1/(numberOfFrames-1) 10*(1/(numberOfFrames-1))],...
            'Position',[2 2 150 30], 'callback',{@slidercb,filters1, 5});


      f = figure(6) ; 
       s = uicontrol(f,'Style','slider',...
            'Min',1,'Max',numberOfFrames,'Value',1,...
            'SliderStep',[1/(numberOfFrames-1) 10*(1/(numberOfFrames-1))],...
            'Position',[2 2 150 30], 'callback',{@slidercb,filters2, 6});



            
      
           
           
       else
        f = figure(5) ; 
       s = uicontrol(f,'Style','slider',...
            'Min',1,'Max',numberOfFrames,'Value',1,...
            'SliderStep',[1/(numberOfFrames-1) 10*(1/(numberOfFrames-1))],...
            'Position',[2 2 150 30], 'callback',{@slidercb2,filters1,weights1, 5});


      f = figure(6) ; 
       s = uicontrol(f,'Style','slider',...
            'Min',1,'Max',numberOfFrames,'Value',1,...
            'SliderStep',[1/(numberOfFrames-1) 10*(1/(numberOfFrames-1))],...
            'Position',[2 2 150 30], 'callback',{@slidercb2,filters2,weights2, 6});
%         
        
      end
       
      % Concert varainces to centimeter uncertainty
      variances1 = 100*((variances1./fs).*speedOfSound) ;
      variances2 = 100*((variances2./fs).*speedOfSound) ; 
      
      %%
      
        positionSpeck1 = [120,-0] ; 
        positionSpeck2 = [-100,-120] ; 
        [positionRefPoints] = getPositionRefPoints (positionSpeck1,positionSpeck2,[distanceMovSpeck1; distanceMovSpeck2]);
        positions = [positionSpeck1;positionSpeck2;positionRefPoints] ; 
        
          f = figure(7) ; 
       s = uicontrol(f,'Style','slider',...
            'Min',1,'Max',numberOfFrames,'Value',1,...
            'SliderStep',[1/(numberOfFrames-1) 10*(1/(numberOfFrames-1))],...
            'Position',[2 2 150 30], 'callback',{@slidercb3,[convergedDistance1,convergedDistance2],[variances1, variances2],positions,c,siz, 7});


        
        
        %% Save Images
        
%              CaptureImagesFromPlot( [convergedDistance1,convergedDistance2],[variances1, variances2],positions,c,siz);
   
            
