function rc = PsychCamSettings(cmd, grabber, varargin)
% rc = PsychCamSettings(cmd, grabber [, arg0, arg1, ...])
%
% Setup tool for video sources for use with Psychtoolbox
% video capture functions. This function can mostly only operate
% on IIDC/DCAM machine vision standard compliant camera connected
% via IEEE1394-Firewire bus or USB bus + IIDC protocol. In other
% words, it only operates on cameras controlled via the libdc1394
% firewire video capture engine on Linux and OSX, not on standard
% consumer/prosumer class cameras controlled by the GStreamer engine,
% e.g., webcams, DV cameras etc.
%
% PsychCamSettings is used to query or change settings of such
% a video source that is supported by the Screen() subfunctions
% for video capture. The first parameter, a 'cmd' command string
% specifies the subfunction to execute. The second parameter
% 'grabber' is the device handle returned by Screen('OpenVideoCapture').
% All following parameters are dependent on the selected subfunction.
%
% Subfunctions of form 'XXX' change the setting of a parameter XXX.
% Subfunctions of form 'XXX' query the current setting of a parameter XXX.
% Subfunctions of form 'AutoXXX' try to switch parameter XXX into automatic
% control, if supported by the camera or video source.
%
% The type of parameters supported is highly dependent on the specific video
% source. Unsupported parameters are usually "no operations" - silently ignored.
%
% For known camera models, the tool will try to map physically meaningful
% values into camera settings and return camera settings as meaningful properties.
% E.g., exposure time is expected and returned in milliseconds, if the mapping
% for the specific camera is known. Otherwise it is set and returned in arbitrary
% units.
%
% Currently known cameras: Basler A602f grayscale firewire camera. Unibrain Fire-i
% firewire camera. Basler A312fc, AVT Marlin F033.
%
% Parameters and their meaning:
%
% curval = PsychCamSettings('ExposureTime', grabber, val)
% -- Set and/or return current exposure time in milliseconds for supported cams, and
% in raw camera specific system units for unknown cameras.
%
% curval = PsychCamSettings('Brightness', grabber, val)
% -- Set/Return brightness setting in arbitrary units. Brightness is the DC offset
% added to the CCD sensor signal before amplification and A/D conversion.
%
% curval = PsychCamSettings('Gain', grabber, val)
% -- Set/Return gain setting in arbitrary units. Gain is the multiplier
% applied to the CCD sensor signal during amplification and before A/D conversion.
%
% curval = PsychCamSettings('Gamma', grabber, val)
% -- Set/Return gamma setting. Gamma is used to influence or set gamma correction.
%
% curval = PsychCamSettings('Sharpness', grabber, val)
% -- Set/Return sharpness setting in arbitrary units. Manipulates digital post-
% processing of images in the camera.
%
% curval = PsychCamSettings('WhiteBalance', grabber, val)
% -- Set/Return white-balance setting in arbitrary units. Only meaningful on color cams.
%
% curval = PsychCamSettings('Saturation', grabber, val)
% -- Set/Return color saturation setting in arbitrary units. Only meaningful on color cams.
%
% Similar to above are queries and (auto-)settings for: Hue, WhiteShading, Iris, Focus, Pan, Tilt,
% Zoom, CaptureQuality, CaptureSize, Temperature, FrameRate, OpticalFilter and TriggerDelay.
%
% curval = PsychCamSettings('BacklightCompensation', grabber, val)
% -- Set/Return setting for backlight compensation mode. Backlight compensation is active
% if control of exposure time, gain and brightness are switched to automatic. It defines
% the algorithm to use for computing the overall image brightness. This is currently
% only supported on the Unibrain Fire-i camera and has the following meaning:
%
% 0 = Off. Just average across image.
% 1 = Use a disc in the center of the image.
% 2 = Use some weighted mix of a disc in the image center and the area outside the disc.
% 3 = Use some portrait mode for optimal exposure of a person sitting in front of the cam.
% 4 = Use upper third of image.
% 5 = Use middle third of image.
% 6 = Use lower third of image.
%
% Special commands and their meaning:
%
% vendor = PsychCamSettings('GetVendor', grabber) 
% -- Return camera vendor name string.
%
% model = PsychCamSettings('GetModel', grabber) 
% -- Return camera model name string.
%
% known = PsychCamSettings('IsKnownCamera', grabber)
% -- Return 1, if this camera is known to PsychCamSettings,
% so it can accept and return meaningful physical properties,
% instead of unknown device units.
%
% settings = PsychCamSettings('AutomateAllSettings', grabber)
% -- Switch all settings into automatic control mode, where possible,
% and return the current settings in a struct.
%
% settings = PsychCamSettings('GetAllSettings', grabber)
% -- Return all known settings in a struct.
%
% oldsettings = PsychCamSettings('SetAllSettings', grabber, settings)
% -- Set all settings from a struct 'settings'.
%
% latency = PsychCamSettings('EstimateLatency', grabber [, fps])
% -- This command analyses the current camera settings and, based on
% the specification of the camera, tries to estimate the latency (in
% seconds) between start of exposure of a video frame and arrival of
% the video frame in the computers video buffer.
%
% If you subtract this 'latency' value from the 'capturetimestamp'
% returned by the Screen('GetCapturedImage') function, you should
% get an estimate of when (in system time) the image was
% actually acquired that corresponds to the captured image.
%
% Please note: In theory Basler firewire cameras with "SFF cycle timer support"
% should be able to report the start time of image exposure directly via
% the returned capture timestamps if the command ...
% Screen('SetVideoCaptureParameter', grabber, 'BaslerFrameTimestampEnable');
% ... was executed before start of capture, therefore making this 'EstimateLatency'
% dance redundant. In practice though, the two tested Basler cameras which should
% support this feature didn't work properly, ie., they returned completely bogus
% capture timestamps when 'BaslerFrameTimestampEnable' was called. Your mileage
% may therefore vary, but it is at least worth a try if you possess a Basler camera.
% 
% Latency is defined as:
%
% Duration of sensor-exposure + sensor readout delay + transmission
% onset delay + time needed for data transfer over Firewire bus or
% other device.
%
% The latency values computed here are based on the official camera
% specs. If the spec is wrong or inaccurate, then this value will
% be wrong or inaccurate as well, so use with caution!
%

% History:
% 20.06.2006 Written (MK).
% 25.12.2013 Updated for support of Basler A312fc and more parameters.
% 07.02.2014 Add support for AVT Marlin F033.

if nargin < 1 || isempty(cmd)
   error('No subcommand specified!');
end

if nargin < 2 || isempty(grabber)
   error('No grabber handle provided!');
end

% Query model and vendor, match it against our camlist:
vendor = Screen('SetVideoCaptureParameter', grabber, 'GetVendorname');
model  = Screen('SetVideoCaptureParameter', grabber, 'GetModelname');
rc = model;

% Check for Unibrain Fire-i...
if ~isempty(strfind(vendor, 'Unibrain')) && ~isempty(strfind(model, 'Fire-i'))
   IsFirei = 1;
else
   IsFirei = 0;
end

% Check for Basler A602f...
if ~isempty(strfind(vendor, 'Basler')) && ~isempty(strfind(model, 'A602f'))
   IsA602f = 1;
else
   IsA602f = 0;
end

% Check for Basler A312fc...
if ~isempty(strfind(vendor, 'Basler')) && ~isempty(strfind(model, 'A312f'))
   IsA312fc = 1;
else
   IsA312fc = 0;
end

% Check for AVT Marlin F033...
if ~isempty(strfind(vendor, 'AVT')) && ~isempty(strfind(model, 'Marlin F033'))
   IsMarlinF033 = 1;
else
   IsMarlinF033 = 0;
end

if strcmp(cmd, 'GetVendor')
     rc = vendor;
     return;
end

if strcmp(cmd, 'GetModel')
     rc = model;
     return;
end

if strcmp(cmd, 'IsKnownCamera')
     % Return 1, if this camera is known to the script.
     rc = IsFirei || IsA602f || IsA312fc || IsMarlinF033;
     return;
end

if strcmp(cmd, 'EstimateLatency')
     % Estimate the latency from start of exposure to the point in time,
     % when the captured frame has finally fully arrived in system memory.
     % One can compute real start of frame acquisition by subtracting this
     % latency value from the capturetimestamp that is returned by
     % Screen('GetCapturedImage').

     % Latency is the sum of exposure time + transmission start delay +
     % transmission time.

     % Transmission time is defined by the IEEE-1394 firewire standard as
     % 1.0 secs / captureframerate:
     if nargin < 3
        % No hypothetical framerate specified. Query it from system.
        t_transmission = 1.0 / Screen('SetVideoCaptureParameter', grabber, 'GetFramerate');
     else 
        % Use provided value instead of querying real setting.
        t_transmission = 1.0 / varargin{1};
     end

     % Transmission onset delay and exposuretime are camera specific and
     % need to handled therefore by camera.
     
     % If we do not know the camera, then we do not know the real exposure time.
     % We just assume zero exposure time and a generic transmission delay of
     % 1 firewire bus cycle at 400 Mbit, which is a good guess:
     t_transdelay = 0.000125;
     t_exposure = 0;

     if IsFirei
        % We do not know the real transmission delay of Fire-i for sure, but
        % as the transmission time of Fire-i is at least 33.33 ms at maximum framerate,
        % and sensor readout time is probably much less, even at maximum ROI, it is
        % reasonable to assume that we have the standard delay of 125 microseconds.
        % Therefore, we leave transdelay at its default and just query exposure time.
        t_exposure = DoFireICamExposureTime(grabber, []) / 1000.0;
     end

     if IsA602f
        % Exposure can be queried.
        t_exposure = DoA602fCamExposureTime(grabber, []) / 1000.0;

        % Transmission delay depends on relative duration of transmission time and
        % sensor readout time. It is calculated as (ROIheight+3)*15.28 microseconds.
        readout_time = (RectHeight(Screen('SetVideoCaptureParameter', grabber, 'GetROI')) + 3);
        readout_time = readout_time * 0.00001528;

        if t_transmission > readout_time
            % Transmission delay is 125 microseconds.
            t_transdelay = 0.000125;
        else
            % Transmission delay is this.
            t_transdelay = (readout_time - t_transmission) + 0.000125;
        end
     end

     if IsA312fc
        % Exposure can be queried.
        t_exposure = DoA312fcCamExposureTime(grabber, []) / 1000.0;

        % Transmission delay depends on relative duration of transmission time and
        % sensor readout time. It is calculated as ROIheight * 28.27 microseconds
        % + 2178 usecs.
        readout_time = RectHeight(Screen('SetVideoCaptureParameter', grabber, 'GetROI'));
        readout_time = readout_time * 0.00002827 + 0.002178;

        if t_transmission > readout_time
            % Transmission start delay is 125 microseconds.
            t_transdelay = 0.000125;
        else
            % Transmission delay is this.
            t_transdelay = (readout_time - t_transmission) + 0.000125;
        end
     end

     if IsMarlinF033
        % Exposure can be queried.
        t_exposure = DoMarlinF033CamExposureTime(grabber, []) / 1000.0;

        % Transmission delay depends on relative duration of transmission time and
        % sensor readout time. It is calculated as ROIheight * 27.1 microseconds
        % + 98.5 usecs + 3.45 usecs * nonROIheight:
        aoiheight = RectHeight(Screen('SetVideoCaptureParameter', grabber, 'GetROI'));
        readout_time = aoiheight * 0.0000271 + 0.0000985;
        readout_time = readout_time + 0.00000345 * (494 - aoiheight);

        % Transmission start delay is 500 microseconds +/- 62.5 usecs.
        if t_transmission > readout_time
            % Transmission start delay is 500 microseconds.
            t_transdelay = 0.000500;
        else
            % Transmission delay is this.
            t_transdelay = (readout_time - t_transmission) + 0.000500;
        end
     end

     % Ok, final result:
     tdelay = t_exposure + t_transdelay + t_transmission;

     % Return value.
     rc = tdelay;
     return;
end

if strcmp(cmd, 'AutomateAllSettings')
     % Switch all settings to Auto control mode if possible, return current settings.
     rc.brightness = Screen('SetVideoCaptureParameter', grabber, 'AutoBrightness');
     rc.gain = Screen('SetVideoCaptureParameter', grabber, 'AutoGain');
     rc.shutter = Screen('SetVideoCaptureParameter', grabber, 'AutoShutter');
     rc.exposure = Screen('SetVideoCaptureParameter', grabber, 'AutoExposure');
     rc.sharpness = Screen('SetVideoCaptureParameter', grabber, 'AutoSharpness');
     [rc.whitebalance1, rc.whitebalance2] = Screen('SetVideoCaptureParameter', grabber, 'AutoWhiteBalance');
     [rc.whiteshadingR, rc.whiteshadingG, rc.whiteshadingB] = Screen('SetVideoCaptureParameter', grabber, 'AutoWhiteShading');
     rc.saturation = Screen('SetVideoCaptureParameter', grabber, 'AutoSaturation');
     rc.gamma = Screen('SetVideoCaptureParameter', grabber, 'AutoGamma');
     rc.vendor = Screen('SetVideoCaptureParameter', grabber, 'GetVendorname');
     rc.model = Screen('SetVideoCaptureParameter', grabber, 'GetModelname');
     rc.hue = Screen('SetVideoCaptureParameter', grabber, 'AutoHue');
     rc.iris = Screen('SetVideoCaptureParameter', grabber, 'AutoIris');
     rc.Focus = Screen('SetVideoCaptureParameter', grabber, 'AutoFocus');
     rc.Zoom = Screen('SetVideoCaptureParameter', grabber, 'AutoZoom');
     rc.Pan = Screen('SetVideoCaptureParameter', grabber, 'AutoPan');
     rc.Tilt = Screen('SetVideoCaptureParameter', grabber, 'AutoTilt');
     rc.OpticalFilter = Screen('SetVideoCaptureParameter', grabber, 'AutoOpticalFilter');
     rc.CaptureSize = Screen('SetVideoCaptureParameter', grabber, 'AutoCaptureSize');
     rc.CaptureQuality = Screen('SetVideoCaptureParameter', grabber, 'AutoCaptureQuality');
     rc.FrameRate = Screen('SetVideoCaptureParameter', grabber, 'FrameRate'); % Auto makes no sense.
     rc.TriggerDelay = Screen('SetVideoCaptureParameter', grabber, 'TriggerDelay'); % Auto makes no sense.
     rc.Temperature = Screen('SetVideoCaptureParameter', grabber, 'AutoTemperature');
     return;
end

if strcmp(cmd, 'GetAllSettings')
     % Read out all known settings and store them in returned struct:
     rc.brightness = Screen('SetVideoCaptureParameter', grabber, 'Brightness');
     rc.gain = Screen('SetVideoCaptureParameter', grabber, 'Gain');
     rc.shutter = Screen('SetVideoCaptureParameter', grabber, 'Shutter');
     rc.exposure = Screen('SetVideoCaptureParameter', grabber, 'Exposure');
     rc.sharpness = Screen('SetVideoCaptureParameter', grabber, 'Sharpness');
     [rc.whitebalance1, rc.whitebalance2] = Screen('SetVideoCaptureParameter', grabber, 'WhiteBalance');
     [rc.whiteshadingR, rc.whiteshadingG, rc.whiteshadingB] = Screen('SetVideoCaptureParameter', grabber, 'WhiteShading');
     rc.saturation = Screen('SetVideoCaptureParameter', grabber, 'Saturation');
     rc.gamma = Screen('SetVideoCaptureParameter', grabber, 'Gamma');
     rc.vendor = Screen('SetVideoCaptureParameter', grabber, 'GetVendorname');
     rc.model = Screen('SetVideoCaptureParameter', grabber, 'GetModelname');
     rc.hue = Screen('SetVideoCaptureParameter', grabber, 'Hue');
     rc.iris = Screen('SetVideoCaptureParameter', grabber, 'Iris');
     rc.Focus = Screen('SetVideoCaptureParameter', grabber, 'Focus');
     rc.Zoom = Screen('SetVideoCaptureParameter', grabber, 'Zoom');
     rc.Pan = Screen('SetVideoCaptureParameter', grabber, 'Pan');
     rc.Tilt = Screen('SetVideoCaptureParameter', grabber, 'Tilt');
     rc.OpticalFilter = Screen('SetVideoCaptureParameter', grabber, 'OpticalFilter');
     rc.CaptureSize = Screen('SetVideoCaptureParameter', grabber, 'CaptureSize');
     rc.CaptureQuality = Screen('SetVideoCaptureParameter', grabber, 'CaptureQuality');
     rc.FrameRate = Screen('SetVideoCaptureParameter', grabber, 'FrameRate');
     rc.TriggerDelay = Screen('SetVideoCaptureParameter', grabber, 'TriggerDelay');
     rc.Temperature = Screen('SetVideoCaptureParameter', grabber, 'Temperature');
     return;
end

if strcmp(cmd, 'SetAllSettings')
     % Set all known settings from passed struct:
     if nargin < 3
        error('You need to supply a conformant camera settings struct to SetAllSettings!');
     end

     % Retrieve new settings from rc.
     rc = varargin{1};

     % Set new settings from rc, assign old settings to returned struct rc:
     rc.brightness = Screen('SetVideoCaptureParameter', grabber, 'Brightness', rc.brightness);
     rc.gain = Screen('SetVideoCaptureParameter', grabber, 'Gain', rc.gain);
     rc.shutter = Screen('SetVideoCaptureParameter', grabber, 'Shutter', rc.shutter);
     rc.exposure = Screen('SetVideoCaptureParameter', grabber, 'Exposure', rc.exposure);
     rc.sharpness = Screen('SetVideoCaptureParameter', grabber, 'Sharpness', rc.sharpness);
     [rc.whitebalance1, rc.whitebalance2] = Screen('SetVideoCaptureParameter', grabber, 'WhiteBalance', rc.whitebalance1, rc.whitebalance2);
     [rc.whiteshadingR, rc.whiteshadingG, rc.whiteshadingB] = Screen('SetVideoCaptureParameter', grabber, 'WhiteShading', rc.whiteshadingR, rc.whiteshadingG, rc.whiteshadingB);
     rc.saturation = Screen('SetVideoCaptureParameter', grabber, 'Saturation', rc.saturation);
     rc.gamma = Screen('SetVideoCaptureParameter', grabber, 'Gamma', rc.gamma);
     rc.vendor = Screen('SetVideoCaptureParameter', grabber, 'GetVendorname');
     rc.model = Screen('SetVideoCaptureParameter', grabber, 'GetModelname');
     rc.hue = Screen('SetVideoCaptureParameter', grabber, 'Hue', rc.hue);
     rc.iris = Screen('SetVideoCaptureParameter', grabber, 'Iris', rc.iris);
     rc.Focus = Screen('SetVideoCaptureParameter', grabber, 'Focus', rc.Focus);
     rc.Zoom = Screen('SetVideoCaptureParameter', grabber, 'Zoom', rc.Zoom);
     rc.Pan = Screen('SetVideoCaptureParameter', grabber, 'Pan', rc.Pan);
     rc.Tilt = Screen('SetVideoCaptureParameter', grabber, 'Tilt', rc.Tilt);
     rc.OpticalFilter = Screen('SetVideoCaptureParameter', grabber, 'OpticalFilter', rc.OpticalFilter);
     rc.CaptureSize = Screen('SetVideoCaptureParameter', grabber, 'CaptureSize', rc.CaptureSize);
     rc.CaptureQuality = Screen('SetVideoCaptureParameter', grabber, 'CaptureQuality', rc.CaptureQuality);
     rc.FrameRate = Screen('SetVideoCaptureParameter', grabber, 'FrameRate', rc.FrameRate);
     rc.TriggerDelay = Screen('SetVideoCaptureParameter', grabber, 'TriggerDelay', rc.TriggerDelay);
     rc.Temperature = Screen('SetVideoCaptureParameter', grabber, 'Temperature', rc.Temperature);
     return;
end

if strcmp(cmd, 'BacklightCompensation')
     % This is only meaningful on the Unibrain Fire-i:
     % It maps to the 'Shutter' parameter...
     if IsFirei
        if nargin > 2
           rc = Screen('SetVideoCaptureParameter', grabber, 'Shutter', varargin{1});
        else
           rc = Screen('SetVideoCaptureParameter', grabber, 'Shutter');
        end
     else
        % Ignore it.
        rc = Inf;
     end

     return;
end

if ~isempty(strfind(cmd, 'ExposureTime'))
     if nargin > 2
        inval = varargin{1};
     else
        inval = [];
     end

     if IsA602f || IsFirei || IsA312fc || IsMarlinF033
        % Call cam specific routine.
        if IsA602f
           rc = DoA602fCamExposureTime(grabber, inval);
        end
        
        if IsA312fc
           rc = DoA312fcCamExposureTime(grabber, inval);
        end

        if IsFirei
           rc = DoFireICamExposureTime(grabber, inval);
        end

	if IsMarlinF033
	   rc = DoMarlinF033CamExposureTime(grabber, inval);
	end

        return;
     else
        % Unknown camera. We just fall-through, mapping exposure time to
        % Exposure setting.
        cmd = 'Exposure';
     end
end

if (IsA602f || IsA312fc || IsFirei || IsMarlinF033) && (~isempty(strfind(cmd, 'Exposure')) || ~isempty(strfind(cmd, 'Shutter')))
     % These settings are ignored on the Unibrain Fire-i, Marlin F033 and Basler A602f and A312fc.
     % They either control ExposureTime - for which we have a nice function -
     % or they control special settings...     
     rc = Inf;
     return;
end

% None above matched. We just assume we've got one of the
% generic parameter names and pass it to the Screen - function.
% This will handle the Set- , Get-, and Auto- case.
if nargin > 2
   rc = Screen('SetVideoCaptureParameter', grabber, cmd, varargin{1});
else
   rc = Screen('SetVideoCaptureParameter', grabber, cmd);
end

return;

% This routine maps exposure times in secs to register values and vice versa.
% It works for cameras from the Basler A60xf series. Formulas based on their
% online specs.
function rc = DoA602fCamExposureTime(grabber, inval)
   if isempty(inval)
      % Retrieve and map current setting. The A602f stores the value in Shutter.
      rawval = Screen('SetVideoCaptureParameter', grabber, 'Shutter');
      % Exposure time in seconds is just shuttervalue * timebase (fixed at 20 microseconds.)
      rc = rawval * 0.000020;
   else
      % Map and set new setting. Maximum exposuretime is also limited by 1/fps with fps = capture
      % framerate.
      rawval = min(inval/1000.0, (1.0/Screen('SetVideoCaptureParameter', grabber, 'GetFramerate'))-0.0005);
      rawval = round(rawval / 0.000020);
      % Clamp to valid range of cam: 20 microseconds to 33.3 milliseconds.
      rawval = min(rawval, 1665);
      rawval = max(rawval, 1);
      % Write new settings, return old settings:
      rawval = Screen('SetVideoCaptureParameter', grabber, 'Shutter', rawval);
      % Exposure time in seconds is just shuttervalue * timebase (fixed at 20 microseconds.)
      rc = rawval * 0.000020;
   end
 
   % Convert return value from seconds to milliseconds.
   rc = rc * 1000;
return;

% This routine maps exposure times in secs to register values and vice versa.
% It works for cameras from the Basler A312fc series. Formulas based on their
% online specs.
function rc = DoA312fcCamExposureTime(grabber, inval)
   if isempty(inval)
      % Retrieve and map current setting. The A312fc stores the value in Shutter.
      rawval = Screen('SetVideoCaptureParameter', grabber, 'Shutter');
      % Exposure time in seconds is just shuttervalue * timebase (fixed at 20 microseconds.)
      rc = rawval * 0.000020;
   else
      % Map and set new setting. Maximum exposuretime is also limited by 1/fps with fps = capture
      % framerate - 94.9 usecs, for a given target framerate:
      rawval = min(inval / 1000.0, (1.0 / Screen('SetVideoCaptureParameter', grabber, 'GetFramerate')) - 0.0000949);
      rawval = round(rawval / 0.000020);
      
      % Clamp to valid range of cam: 20 microseconds to 81.9 milliseconds.
      rawval = min(rawval, 4095);
      rawval = max(rawval, 1);
      
      % Write new settings, return old settings:
      rawval = Screen('SetVideoCaptureParameter', grabber, 'Shutter', rawval);
      
      % Exposure time in seconds is just shuttervalue * timebase (fixed at 20 microseconds.)
      rc = rawval * 0.000020;
   end

   % Convert return value from seconds to milliseconds.
   rc = rc * 1000;
return;

% This routine maps exposure times in secs to register values and vice versa.
% It works for cameras from the Unibrain Fire-i and Fire-400 series. Formulas based on their
% online lookup sheet:
function rc = DoFireICamExposureTime(grabber, inval)
   if ~isempty(inval)
      % Map and set new setting. Maximum exposuretime is also limited by 1/fps with fps = capture
      % framerate.
      expval = min(inval/1000.0, (1.0/Screen('SetVideoCaptureParameter', grabber, 'GetFramerate'))-0.0005);
      rawval = 0;

      % Mapping of exposure times to rawvalues is done via this hard-coded lookup table:
      if expval >= 0
         rawval = round(expval / 0.000016);
      end

      if expval > (128 * 0.000016)
         rawval = expval - (128 * 0.000016);
         rawval = 128 + round(rawval / 0.000032);
      end

      if expval > (128 * 0.000016 + (192-128) * 0.000032)
         rawval = expval - (128 * 0.000016 + (192-128) * 0.000032);
         rawval = 192 + round(rawval / 0.000064);
      end

      if expval > (128 * 0.000016 + (192-128) * 0.000032 + (384-192) * 0.000064)
         rawval = expval - (128 * 0.000016 + (192-128) * 0.000032 + (384-192) * 0.000064);
         rawval = 384 + round(rawval / 0.000128);
      end

      % Clamp to valid range of cam: 20 microseconds to 81.9 milliseconds.
      rawval = min(rawval, 511);
      rawval = max(rawval, 0);
      % Write new settings, return old settings:
      rawval = Screen('SetVideoCaptureParameter', grabber, 'Exposure', rawval);
   else
      % Retrieve current setting:
      rawval = Screen('SetVideoCaptureParameter', grabber, 'Exposure');
   end

   % Need to map returned rawval into meaningful value:
   if rawval >= 0 && rawval <= 128
      rc = rawval * 0.000016; % Stepping of 16 microsecs.
   end
   
   if rawval > 128 && rawval <= 192
      rc = 128 * 0.000016 + (rawval-128) * 0.000032; % Stepping of 32 microsecs.
   end

   if rawval > 192 && rawval <= 384
      rc = 128 * 0.000016 + (192-128) * 0.000032 + (rawval-192) * 0.000064; % Stepping of 64 microsecs.
   end

   if rawval > 384
      rc = 128 * 0.000016 + (192-128) * 0.000032 + (384-192) * 0.000064 + (rawval-384) * 0.000128; % Stepping of 128 microsecs.
   end

   % Convert rc from seconds to milliseconds.
   rc = rc * 1000.0;
return;

function rc = DoMarlinF033CamExposureTime(grabber, inval)
   if isempty(inval)
      % Retrieve and map current setting. The Marlin stores the value in Shutter.
      rawval = Screen('SetVideoCaptureParameter', grabber, 'Shutter');
      % Exposure time in seconds is just shuttervalue * timebase (fixed at 20 microseconds.) + 12 usecs offset:
      rc = rawval * 0.000020 + 0.000012;
   else
      % Map and set new setting. Maximum exposuretime is also limited by 1/fps with fps = capture
      % framerate, for a given target framerate:
      rawval = min(inval / 1000.0, (1.0 / Screen('SetVideoCaptureParameter', grabber, 'GetFramerate')));
      rawval = round((rawval - 0.000012) / 0.000020);
      
      % Clamp to valid range of cam: 32 microseconds to 81.912 milliseconds.
      rawval = min(rawval, 4095);
      rawval = max(rawval, 1);
      
      % Write new settings, return old settings:
      rawval = Screen('SetVideoCaptureParameter', grabber, 'Shutter', rawval);
      
      % Exposure time in seconds is just shuttervalue * timebase (fixed at 20 microseconds.) + 12 usecs offset:
      rc = rawval * 0.000020 + 0.000012;
   end

   % Convert return value from seconds to milliseconds.
   rc = rc * 1000;
return;
