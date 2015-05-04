function oldRes=SetResolution(screenNumber,width,height,hz,pixelSize)
% oldRes=SetResolution(screenNumber,width,height,[hz],[pixelSize])
% oldRes=SetResolution(screenNumber,res)
% 
% Set the resolution of the screen.  This is intended to be used in
% programs that run psychophysical experiments, so SetResolution is
% strict. It issues a fatal error unless it can provide exactly the
% resolution you specified. (For lenience, look at NearestResolution.)
% "screenNumber" is the screen number.
% "width" and "height" are the desired dimensions in pixels.  
% "hz" is the desired refresh rate (default is current frame rate).  
% "pixelSize" 8, 16, 24, or 32 bits and defaults to current pixelSize.
% Returns the current resolution as "oldRes".  
% 
%   oldRes=SetResolution(0,1024,768,75);
%   w=Screen(0,'OpenWindow');
%   Screen(w,'PutImage',image);
%   Screen(w,'Close');
%   SetResolution(0, oldRes);
% 
% To display a list of available resolutions, try ResolutionTest. Also see
% NearestResolution, ResolutionTest, and Screen('Resolution')
% and Screen('Resolutions').
% 
% NOTE: Apple has all the new LCD screens return a frame rate of 0, so
% we treat that value as a special case. A request for "hz" of NaN will
% match only with a frame rate of 0.
%
% NOTE: On Linux this only works as you'd expect if there is only one
% video output connected to a Psychtoolbox screen / X-Screen. On a
% setup with multiple outputs per screen, this will only change the
% size of the framebuffer, but not the resolution of the actual displays!
% Use Screen('ConfigureDisplay', 'Scanout', ...); to change settings on
% such a multi-display per screen setup on Linux.
% 
% Originally written by Sabina Wolfson.

% HISTORY:
% 1/27/00 SSW Wrote it.
% 1/28/00 dgp Cosmetic editing. Made screenNumber first argument. Renamed from set_resolution to SetResolution.
% 4/9/02  dgp Check the width and height arguments.
% 4/13/02 dgp Cosmetic.
% 4/29/02 dgp Screen Resolutions is only available on Mac.
% 6/6/02  dgp Accept res instead of parameter list.
% 6/7/02  dgp Hz value of NaN matches NaN.
% 9/23/07 mk  Adapted for Psychtoolbox-3.
% 4/15/14 mk  Adapted for multi-output per x-screen setups on Linux.

if nargin<2 || nargin>5
    error(sprintf('%s\n%s','USAGE: oldRes=SetResolution(screenNumber,width,height,[hz],[pixelSize])',...
                           '       oldRes=SetResolution(screenNumber,res)'));
end

if ~ismember(screenNumber, Screen('Screens'))
    error(sprintf('Invalid screenNumber %d.',screenNumber));
end

oldRes = Screen('Resolution', screenNumber);

if nargin<5 || isempty(pixelSize)
    pixelSize=oldRes.pixelSize;
end

if nargin<4 || isempty(hz)
    hz=oldRes.hz;
end

if nargin==2 && isa(width,'struct')
    res=width;
    if isfield(res,'width')
        width=res.width;
    end
    
    if isfield(res,'height')
        height=res.height;
    end
    
    if isfield(res,'hz')
        hz=res.hz;
    end
    
    if isfield(res,'pixelSize')
        pixelSize=res.pixelSize;
    end
end

if ~exist('height','var')  || isempty(height) || ~isfinite(height)
    error('height (in pixels) must be specified.');
end

if ~exist('width','var') || isempty(width) || ~isfinite(width)
    error('width (in pixels) must be specified.');
end

if isnan(hz)
    hz = 0;
end

% One display output per screen? Then we can validate the settings.
if Screen('ConfigureDisplay', 'NumberOutputs', screenNumber) < 2
    res = Screen('Resolutions', screenNumber, 0); % get a list of available resolutions

    resIndex=find([res.width]==width & [res.height]==height & [res.hz]==hz);
    if isempty(resIndex)
        error(sprintf('Can''t find a resolution of %g x %g x %g Hz. Resolution remains unchanged.',width,height,hz));
    end
    
    nres = res(resIndex(1));
else
    % Multiple outputs per screen. Can't validate, but take settings at face value and hope for the best:
    nres.width = width;
    nres.height = height;
    nres.hz = hz;
    nres.pixelSize = pixelSize;
    warning('SetResolution: Only changing virtual resolution of multi-display screen %i!\n');
end

if isempty(pixelSize) || ~isfinite(pixelSize)
    nres.pixelSize = [];
else
    nres.pixelSize = pixelSize;
end

oldRes = Screen('Resolution', screenNumber, nres.width, nres.height, nres.hz, nres.pixelSize);

return;
