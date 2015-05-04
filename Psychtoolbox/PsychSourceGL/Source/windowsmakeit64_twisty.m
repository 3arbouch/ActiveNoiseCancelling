function windowsmakeit64_twisty(what, onoctave)
% Builds the 64-Bit Psychtoolbox on MS-Windows for Octave and Matlab.
% This script is customized for MK's build machine "twisty", building
% against a VirtualBox VM running MS-Windows 7 64-Bit.
%

if ~(IsWin(1) || (IsWin && IsOctave))
    error('%s must be run on MS-Windows within 32-Bit Octave or 64-Bit Matlab!', mfilename);
end

if nargin < 1
    what = 0;
end

if nargin < 2
    onoctave = IsOctave;
end

% Rebuild all request?
if what == -1
    % Yes: Call ourselves recursively on all plugins/modes to rebuild
    % everything:
    tic;
    for what = 0:13
        windowsmakeit64_twisty(what);
    end
    elapsedsecs = toc;
    fprintf('Total rebuild time for all mex files was %f seconds. Bye.\n\n', elapsedsecs);
    return;
end

try

% Matlab or Octave build?
if onoctave == 0
    % Matlab build:
    if what == 0
        % Default: Build Screen with GStreamer support: Needs the
        % www.gstreamer.com GStreamer SDK for 64-Bit Windows.
        % Use this for verbose linker output: /VERBOSE:LIB 
        % MSVCRT 2010:
        if 1
            % GStreamer 1.0:
            mex -v -outdir ..\Projects\Windows\build -output Screen -DPTBMODULE_Screen -largeArrayDims -DPTB_USE_GSTREAMER -DGLEW_STATIC -LC:\gstreamer\1.0\x86_64\lib -IC:\gstreamer\1.0\x86_64\include -IC:\gstreamer\1.0\x86_64\include\gstreamer-1.0 -IC:\gstreamer\1.0\x86_64\include\glib-2.0 -IC:\gstreamer\1.0\x86_64\include\glib-2.0\include -IC:\gstreamer\1.0\x86_64\lib\glib-2.0\include -IC:\gstreamer\1.0\x86_64\include\libxml2 -I"C:\Program Files\Microsoft SDKs\Windows\v7.1\Include" -ICommon\Base -ICommon\Screen -IWindows\Base -IWindows\Screen Windows\Screen\*.c Windows\Base\*.c Common\Base\*.c Common\Screen\*.c kernel32.lib user32.lib gdi32.lib advapi32.lib glu32.lib opengl32.lib winmm.lib delayimp.lib -lgobject-2.0 -lglib-2.0 -lgstreamer-1.0 -lgstapp-1.0 -lgstvideo-1.0 -lgstpbutils-1.0 LINKFLAGS="$LINKFLAGS /DELAYLOAD:libgobject-2.0-0.dll /DELAYLOAD:libglib-2.0-0.dll /DELAYLOAD:libgstreamer-1.0-0.dll /DELAYLOAD:libgstapp-1.0-0.dll /DELAYLOAD:libgstvideo-1.0-0.dll /DELAYLOAD:libgstpbutils-1.0-0.dll"
        else
            % Legacy 0.10 GStreamer:
            mex -v -outdir ..\Projects\Windows\build -output Screen -DPTBMODULE_Screen -largeArrayDims -DPTB_USE_GSTREAMER -DGLEW_STATIC -LC:\gstreamer-sdk\0.10\x86_64\lib -IC:\gstreamer-sdk\0.10\x86_64\include -IC:\gstreamer-sdk\0.10\x86_64\include\gstreamer-0.10 -IC:\gstreamer-sdk\0.10\x86_64\include\glib-2.0 -IC:\gstreamer-sdk\0.10\x86_64\include\glib-2.0\include -IC:\gstreamer-sdk\0.10\x86_64\lib\glib-2.0\include -IC:\gstreamer-sdk\0.10\x86_64\include\libxml2 -I"C:\Program Files\Microsoft SDKs\Windows\v7.1\Include" -ICommon\Base -ICommon\Screen -IWindows\Base -IWindows\Screen Windows\Screen\*.c Windows\Base\*.c Common\Base\*.c Common\Screen\*.c kernel32.lib user32.lib gdi32.lib advapi32.lib glu32.lib opengl32.lib winmm.lib delayimp.lib -lgobject-2.0 -lgthread-2.0 -lglib-2.0 -lgstreamer-0.10 -lgstapp-0.10 -lgstinterfaces-0.10 -lgstpbutils-0.10 LINKFLAGS="$LINKFLAGS /DELAYLOAD:libgobject-2.0-0.dll /DELAYLOAD:libgthread-2.0-0.dll /DELAYLOAD:libglib-2.0-0.dll /DELAYLOAD:libgstreamer-0.10-0.dll /DELAYLOAD:libgstapp-0.10-0.dll /DELAYLOAD:libgstinterfaces-0.10-0.dll /DELAYLOAD:libgstpbutils-0.10-0.dll"
        end

        % Windows CRT build against MSCRT.dll in the Windows DDK: Didn't
        % help in any way fixing Matlab-64 Bit + Java JVM + GStreamer
        % crashes on some few machines:
        %         mex -v -outdir ..\Projects\Windows\build -output Screen -DPTBMODULE_Screen -largeArrayDims -DPTB_USE_GSTREAMER -DGLEW_STATIC ...
        %             COMPFLAGS="$COMPFLAGS " ...
        %             -DMBCS=1 -DWHISTLER_DDK=1 -D_DLL=1 -D_MT=1 -D__NO_MINGW_LFS -D__MSVCRT_VERSION__=0x0601 -D__MSVCRT__ -D_CRT_NONSTDC_FORCE_DEPRECATE ...
        %             -DWIN32_LEAN_AND_MEAN=1 -D_CRT_SECURE_FORCE_DEPRECATE -D_CRT_SECURE_NO_WARNINGS -D_CRT_SECURE_NO_DEPRECATE -D_CRT_NONSTDC_NO_DEPRECATE ...
        %             -IC:\WinDDK\7600.16385.1\inc\crt -IC:\WinDDK\7600.16385.1\inc\atl71 ...
        %             -LC:\gstreamer-sdk\0.10\x86_64\lib -IC:\gstreamer-sdk\0.10\x86_64\include -IC:\gstreamer-sdk\0.10\x86_64\include\gstreamer-0.10 ...
        %             -IC:\gstreamer-sdk\0.10\x86_64\include\glib-2.0 -IC:\gstreamer-sdk\0.10\x86_64\include\glib-2.0\include -IC:\gstreamer-sdk\0.10\x86_64\lib\glib-2.0\include ...
        %             -IC:\gstreamer-sdk\0.10\x86_64\include\libxml2 -I"C:\Program Files\Microsoft SDKs\Windows\v7.1\Include" ...
        %             -ICommon\Base -ICommon\Screen -IWindows\Base -IWindows\Screen ...
        %             Windows\Screen\*.c Windows\Base\*.c Common\Base\*.c Common\Screen\*.c ...
        %             oldnames.lib msvcrt.lib ...
        %             kernel32.lib user32.lib gdi32.lib advapi32.lib glu32.lib opengl32.lib winmm.lib delayimp.lib ...
        %             -lgobject-2.0 -lgthread-2.0 -lglib-2.0 -lgstreamer-0.10 -lgstapp-0.10 -lgstinterfaces-0.10 -lgstpbutils-0.10 ...
        %             LINKFLAGS="/LIBPATH:C:\WinDDK\7600.16385.1\lib\crt\amd64 /LIBPATH:C:\WinDDK\7600.16385.1\lib\win7\amd64 $LINKFLAGS ...
        %             /NODEFAULTLIB /OPT:NOREF ...
        %             /DELAYLOAD:libgobject-2.0-0.dll /DELAYLOAD:libgthread-2.0-0.dll /DELAYLOAD:libglib-2.0-0.dll /DELAYLOAD:libgstreamer-0.10-0.dll /DELAYLOAD:libgstapp-0.10-0.dll /DELAYLOAD:libgstinterfaces-0.10-0.dll /DELAYLOAD:libgstpbutils-0.10-0.dll"

        movefile(['..\Projects\Windows\build\Screen.' mexext], [PsychtoolboxRoot 'PsychBasic\MatlabWindowsFilesR2007a\']);
    end

    if what == 1
        % Build WaitSecs
        mex -v -outdir ..\Projects\Windows\build -output WaitSecs -DPTBMODULE_WaitSecs -largeArrayDims  -I"C:\Program Files\Microsoft SDKs\Windows\v7.1\Include" -ICommon\Base -ICommon\WaitSecs -IWindows\Base Windows\Base\*.c Common\Base\*.c Common\WaitSecs\*.c kernel32.lib user32.lib winmm.lib
        movefile(['..\Projects\Windows\build\WaitSecs.' mexext], [PsychtoolboxRoot 'PsychBasic\MatlabWindowsFilesR2007a\']);
    end

    if what == 2
        % Build PsychPortAudio
        mex -v -outdir ..\Projects\Windows\build -output PsychPortAudio -DPTBMODULE_PsychPortAudio -largeArrayDims -L..\Cohorts\PortAudio -I"C:\Program Files\Microsoft SDKs\Windows\v7.1\Include" -ICommon\Base -ICommon\PsychPortAudio -IWindows\Base Windows\Base\*.c Common\Base\*.c Common\PsychPortAudio\*.c kernel32.lib user32.lib winmm.lib delayimp.lib -lportaudio_x64 LINKFLAGS="$LINKFLAGS /DELAYLOAD:portaudio_x64.dll"
        movefile(['..\Projects\Windows\build\PsychPortAudio.' mexext], [PsychtoolboxRoot 'PsychBasic\MatlabWindowsFilesR2007a\']);
    end

    if what == 3
        % Build GetSecs
        mex -v -outdir ..\Projects\Windows\build -output GetSecs -DPTBMODULE_GetSecs -largeArrayDims  -I"C:\Program Files\Microsoft SDKs\Windows\v7.1\Include" -ICommon\Base -ICommon\GetSecs -IWindows\Base Windows\Base\*.c Common\Base\*.c Common\GetSecs\*.c kernel32.lib user32.lib winmm.lib
        movefile(['..\Projects\Windows\build\GetSecs.' mexext], [PsychtoolboxRoot 'PsychBasic\MatlabWindowsFilesR2007a\']);
    end

    if what == 4
        % Build IOPort
        mex -v -outdir ..\Projects\Windows\build -output IOPort -DPTBMODULE_IOPort -largeArrayDims  -I"C:\Program Files\Microsoft SDKs\Windows\v7.1\Include" -ICommon\Base -ICommon\IOPort -IWindows\Base -IWindows\IOPort Windows\Base\*.c Common\Base\*.c Common\IOPort\*.c Windows\IOPort\*.c kernel32.lib user32.lib winmm.lib
        movefile(['..\Projects\Windows\build\IOPort.' mexext], [PsychtoolboxRoot 'PsychBasic\MatlabWindowsFilesR2007a\']);
    end

    if what == 5
        % Build PsychCV
        % Disabled for now: As long as it contains 3rd party code, we can't
        % really distribute it in a precompiled version for Matlab under our
        % new license. Distribution of compiled mex files for octave would
        % be possible, but see below...
        % mex -v -outdir ..\Projects\Windows\build -output PsychCV -DPTBMODULE_PsychCV -largeArrayDims  -ID:\install\QuickTimeSDK\CIncludes -I"C:\Program Files\Microsoft SDKs\Windows\v7.1\Include" -ICommon\Base -ICommon\PsychCV -IWindows\Base -I..\Cohorts\ARToolkit\include Windows\Base\*.c Common\Base\*.c Common\PsychCV\*.c kernel32.lib user32.lib gdi32.lib advapi32.lib glu32.lib opengl32.lib winmm.lib delayimp.lib libARvideo.lib libARgsub.lib libARgsub_lite.lib libARgsubUtil.lib libARMulti.lib libAR.lib
        % movefile(['..\Projects\Windows\build\PsychCV.' mexext], [PsychtoolboxRoot 'PsychBasic\MatlabWindowsFilesR2007a\']);
    end

    if what == 7
        % Build Eyelink:
        mex -v -outdir ..\Projects\Windows\build -output Eyelink -DPTBMODULE_Eyelink -largeArrayDims -I"C:\Program Files\Microsoft SDKs\Windows\v7.1\Include" -I"C:\Program Files (x86)\SR Research\EyeLink\Includes\eyelink" -ICommon\Base -ICommon\Eyelink -IWindows\Base Windows\Base\*.c Common\Base\*.c Common\Eyelink\*.c user32.lib gdi32.lib advapi32.lib winmm.lib "C:\Program Files (x86)\SR Research\EyeLink\libs\x64\eyelink_core64.lib" "C:\Program Files (x86)\SR Research\EyeLink\libs\x64\eyelink_w32_comp64.lib"
        movefile(['..\Projects\Windows\build\Eyelink.' mexext], [PsychtoolboxRoot 'PsychBasic\MatlabWindowsFilesR2007a\']);
    end

    if what == 8
        % Build PsychKinectCore:
        mex -v -outdir ..\Projects\Windows\build -output PsychKinectCore -DPTBMODULE_PsychKinectCore -largeArrayDims -DWIN32 -L..\Cohorts\Kinect-v16-withsource\libusb\lib\msvc_x64 -I..\Cohorts\Kinect-v16-withsource -I"C:\Program Files\Microsoft SDKs\Windows\v7.1\Include" -ICommon\Base -IWindows\Base -ICommon\PsychKinect Windows\Base\*.c Common\Base\*.c Common\PsychKinect\*.c ..\Cohorts\Kinect-v16-withsource\*.cpp kernel32.lib user32.lib winmm.lib -lusb
        movefile(['..\Projects\Windows\build\PsychKinectCore.' mexext], [PsychtoolboxRoot 'PsychBasic\MatlabWindowsFilesR2007a\']);
    end

    if what == 9
        % Build PsychHID:
        mex -v -outdir ..\Projects\Windows\build -output PsychHID -DPTBMODULE_PsychHID -largeArrayDims -DWIN32 -I"C:\Program Files\Microsoft SDKs\Windows\v7.1\Include" -L..\Cohorts\libusb1-win32\MS64\dll -I..\Cohorts\libusb1-win32\include\libusb-1.0 -ICommon\Base -IWindows\Base -ICommon\PsychHID Windows\PsychHID\*.cpp Windows\PsychHID\*.c Windows\Base\*.c Common\Base\*.c Common\PsychHID\*.c -ldinput8 kernel32.lib user32.lib winmm.lib -lusb-1.0 setupapi.lib
        movefile(['..\Projects\Windows\build\PsychHID.' mexext], [PsychtoolboxRoot 'PsychBasic\MatlabWindowsFilesR2007a\']);
    end

else
    % Octave-4 build:
    if what == 0
        % Default: Build Screen.mex
        % TODO FIXME - x86 GStreamer SDK needed
        mexoctave -g -v --output T:\projects\OpenGLPsychtoolbox\Psychtoolbox-3\PsychSourceGL\Projects\Windows\build\Screen.mex -DPTBMODULE_Screen -DGLEW_STATIC -I"C:\Program Files (x86)\Microsoft SDKs\Windows\v7.0A\Include" -ICommon\Base -ICommon\Screen -IWindows\Base -IWindows\Screen Windows\Screen\*.c Windows\Base\*.c Common\Base\*.c Common\Screen\*.c kernel32.lib user32.lib gdi32.lib advapi32.lib glu32.lib opengl32.lib winmm.lib delayimp.lib
        system('copy T:\projects\OpenGLPsychtoolbox\Psychtoolbox-3\PsychSourceGL\Projects\Windows\build\Screen.mex T:\projects\OpenGLPsychtoolbox\Psychtoolbox-3\Psychtoolbox\PsychBasic\Octave3WindowsFiles\');
    end

    if what == 1
        % Build WaitSecs.mex
        mexoctave -g -v --output T:\projects\OpenGLPsychtoolbox\Psychtoolbox-3\PsychSourceGL\Projects\Windows\build\WaitSecs.mex -DPTBMODULE_WaitSecs  -I"C:\Program Files (x86)\Microsoft SDKs\Windows\v7.0A\Include" -ICommon\Base -ICommon\WaitSecs -IWindows\Base Windows\Base\*.c Common\Base\*.c Common\WaitSecs\*.c kernel32.lib user32.lib winmm.lib
        system('copy T:\projects\OpenGLPsychtoolbox\Psychtoolbox-3\PsychSourceGL\Projects\Windows\build\WaitSecs.mex T:\projects\OpenGLPsychtoolbox\Psychtoolbox-3\Psychtoolbox\PsychBasic\Octave3WindowsFiles\');
    end

    if what == 2
        % Build PsychPortAudio.mex
        % TODO FIXME - libportaudio_x86 compatible with mingw32 / gcc needed
        mexoctave -g -v --output T:\projects\OpenGLPsychtoolbox\Psychtoolbox-3\PsychSourceGL\Projects\Windows\build\PsychPortAudio.mex -DPTBMODULE_PsychPortAudio -L..\Cohorts\PortAudio -I"C:\Program Files (x86)\Microsoft SDKs\Windows\v7.0A\Include" -ICommon\Base -ICommon\PsychPortAudio -IWindows\Base Windows\Base\*.c Common\Base\*.c Common\PsychPortAudio\*.c kernel32.lib user32.lib winmm.lib delayimp.lib -lportaudio_x86
        system('copy T:\projects\OpenGLPsychtoolbox\Psychtoolbox-3\PsychSourceGL\Projects\Windows\build\PsychPortAudio.mex T:\projects\OpenGLPsychtoolbox\Psychtoolbox-3\Psychtoolbox\PsychBasic\Octave3WindowsFiles\');
    end

    if what == 3
        % Build GetSecs.mex
        mexoctave -g -v --output T:\projects\OpenGLPsychtoolbox\Psychtoolbox-3\PsychSourceGL\Projects\Windows\build\GetSecs.mex -DPTBMODULE_GetSecs  -I"C:\Program Files (x86)\Microsoft SDKs\Windows\v7.0A\Include" -ICommon\Base -ICommon\GetSecs -IWindows\Base Windows\Base\*.c Common\Base\*.c Common\GetSecs\*.c kernel32.lib user32.lib winmm.lib
        system('copy T:\projects\OpenGLPsychtoolbox\Psychtoolbox-3\PsychSourceGL\Projects\Windows\build\GetSecs.mex T:\projects\OpenGLPsychtoolbox\Psychtoolbox-3\Psychtoolbox\PsychBasic\Octave3WindowsFiles\');
    end

    if what == 4
        % Build IOPort.mex
        mexoctave -g -v --output T:\projects\OpenGLPsychtoolbox\Psychtoolbox-3\PsychSourceGL\Projects\Windows\build\IOPort.mex -DPTBMODULE_IOPort  -I"C:\Program Files (x86)\Microsoft SDKs\Windows\v7.0A\Include" -ICommon\Base -ICommon\IOPort -IWindows\Base -IWindows\IOPort Windows\Base\*.c Common\Base\*.c Common\IOPort\*.c Windows\IOPort\*.c kernel32.lib user32.lib winmm.lib
        system('copy T:\projects\OpenGLPsychtoolbox\Psychtoolbox-3\PsychSourceGL\Projects\Windows\build\IOPort.mex T:\projects\OpenGLPsychtoolbox\Psychtoolbox-3\Psychtoolbox\PsychBasic\Octave3WindowsFiles\');
    end

    if what == 5
        % Build PsychCV.mex
        % NOTE: Link is currently broken. Also we don't build and
        % distribute PsychCV.mex at the moment. Let's see if anybody
        % actually misses this mex file...
        % mexoctave -g -v --output T:\projects\OpenGLPsychtoolbox\Psychtoolbox-3\PsychSourceGL\Projects\Windows\build\PsychCV.mex -DPTBMODULE_PsychCV  -ID:\install\QuickTimeSDK\CIncludes -I"C:\Program Files (x86)\Microsoft SDKs\Windows\v7.0A\Include" -ID:\MicrosoftDirectXSDK\Include -ICommon\Base -ICommon\PsychCV -IWindows\Base -I..\Cohorts\ARToolkit\include Windows\Base\*.c Common\Base\*.c Common\PsychCV\*.c kernel32.lib user32.lib gdi32.lib advapi32.lib glu32.lib opengl32.lib winmm.lib delayimp.lib libARvideo.lib libARgsub.lib libARgsub_lite.lib libARgsubUtil.lib libARMulti.lib libAR.lib
        % system('copy T:\projects\OpenGLPsychtoolbox\Psychtoolbox-3\PsychSourceGL\Projects\Windows\build\PsychCV.mex T:\projects\OpenGLPsychtoolbox\Psychtoolbox-3\Psychtoolbox\PsychBasic\Octave3WindowsFiles\');
    end

    if what == 6
        % Build moglcore.mex
        % TODO FIXME
        mexoctave -g -v --output T:\projects\OpenGLPsychtoolbox\Psychtoolbox-3\PsychSourceGL\Projects\Windows\build\moglcore.mex  -I"C:\Program Files (x86)\Microsoft SDKs\Windows\v7.0A\Include" -IU:\projects\OpenGLPsychtoolbox\Psychtoolbox-3\Psychtoolbox\PsychOpenGL\MOGL\source -DWINDOWS -DGLEW_STATIC windowhacks.c gl_auto.c gl_manual.c mogl_rebinder.c moglcore.c glew.c user32.lib gdi32.lib advapi32.lib glu32.lib opengl32.lib glut32.lib
        system('copy T:\projects\OpenGLPsychtoolbox\Psychtoolbox-3\PsychSourceGL\Projects\Windows\build\moglcore.mex T:\projects\OpenGLPsychtoolbox\Psychtoolbox-3\Psychtoolbox\PsychBasic\Octave3WindowsFiles\');
    end

    if what == 7
        % Build Eyelink.mex
        % TODO FIXME 32-Bit EyeLink SDK needed.
        mexoctave -g -v --output T:\projects\OpenGLPsychtoolbox\Psychtoolbox-3\PsychSourceGL\Projects\Windows\build\Eyelink.mex  -LD:\SRResearch\EyeLink\libs -ID:\SRResearch\EyeLink\Includes\eyelink -I"C:\Program Files (x86)\Microsoft SDKs\Windows\v7.0A\Include" -ICommon\Base -ICommon\Eyelink -IWindows\Base Windows\Base\*.c Common\Base\*.c Common\Eyelink\*.c user32.lib gdi32.lib advapi32.lib winmm.lib eyelink_core.lib eyelink_w32_comp.lib eyelink_exptkit20.lib
        system('copy T:\projects\OpenGLPsychtoolbox\Psychtoolbox-3\PsychSourceGL\Projects\Windows\build\Eyelink.mex T:\projects\OpenGLPsychtoolbox\Psychtoolbox-3\Psychtoolbox\PsychBasic\Octave3WindowsFiles\');
    end

    if what == 8
        % Build PsychKinectCore.mex:
        % TODO FIXME 32-Bit libUSB-1.0 needed.
        mexoctave -g -v --output T:\projects\OpenGLPsychtoolbox\Psychtoolbox-3\PsychSourceGL\Projects\Windows\build\PsychKinectCore.mex -DPTBMODULE_PsychKinectCore -L..\Cohorts\Kinect-v16-withsource\libusb\lib\msvc -I..\Cohorts\Kinect-v16-withsource -ICommon\Base -IWindows\Base -ICommon\PsychKinect Windows\Base\*.c Common\Base\*.c Common\PsychKinect\*.c ..\Cohorts\Kinect-v16-withsource\*.cpp kernel32.lib user32.lib winmm.lib -lusb
        system('copy T:\projects\OpenGLPsychtoolbox\Psychtoolbox-3\PsychSourceGL\Projects\Windows\build\PsychKinectCore.mex T:\projects\OpenGLPsychtoolbox\Psychtoolbox-3\Psychtoolbox\PsychBasic\Octave3WindowsFiles\');
    end;

    if what == 9
        % Build PsychHID.mex:
        % TODO FIXME 32-Bit libUSB-1.0 needed.
        mexoctave -g -v --output T:\projects\OpenGLPsychtoolbox\Psychtoolbox-3\PsychSourceGL\Projects\Windows\build\PsychHID.mex -DPTBMODULE_PsychHID -L..\Cohorts\libusb1-win32\MS32\dll -I..\Cohorts\libusb1-win32\include\libusb-1.0 -ICommon\Base -IWindows\Base -ICommon\PsychHID Windows\PsychHID\*.cpp Windows\PsychHID\*.c Windows\Base\*.c Common\Base\*.c Common\PsychHID\*.c dinput8.lib kernel32.lib user32.lib winmm.lib libusb-1.0.lib setupapi.lib
        system('copy T:\projects\OpenGLPsychtoolbox\Psychtoolbox-3\PsychSourceGL\Projects\Windows\build\PsychHID.mex T:\projects\OpenGLPsychtoolbox\Psychtoolbox-3\Psychtoolbox\PsychBasic\Octave3WindowsFiles\');
    end
end

catch %#ok<*CTCH>
    % Empty. We just want to make sure the delete() call below is executed
    % in both success and failure case.
end

return;

% Special mex wrapper for Octave compile on Windows:
function mexoctave(varargin)
debugme = 0;
callmex = 1;

if (debugme), fprintf('nargin = %i\n', nargin); end;
myvararg = cell();
myvararg(end+1) = '-DPTBOCTAVE3MEX';
myvararg(end+1) = '-Wno-multichar';
myvararg(end+1) = '-Wno-unknown-pragmas';

outarg = '';
quoted = 0;
emitarg = 0;
for i=1:nargin
    curarg = char(varargin(i));

    if debugme, fprintf('Preparse Arg %i: %s\n', i, curarg); end
    if (length(strfind(curarg, '"')) > 0)
        if ~quoted
            % Start of quoted string:
            quoted = 1;
            outarg = [curarg];
            emitarg = 0;
        else
            % End of quoted string: Emit!
            quoted = 0;
            outarg = [ outarg curarg ];
            emitarg = 1;
        end
    else
        % Not start or end string of a quoted piece:
        if quoted
            % Within a quoted segment!
            outarg = [ outarg curarg ];
            emitarg = 0;
        else
            % Outside a quoted segment:
            % Expansion needed?
            ppos = strfind(curarg, '*.c');
            if length(ppos) > 0
                prefix  = curarg(1:ppos(1)-1);
                allfiles = dir(curarg);
                for j=1:length(allfiles)
                    expandedfiles = [ prefix allfiles(j).name ];
                    emitarg = 0;
                    myvararg(end+1) = expandedfiles;
                    if (debugme), fprintf('Emitted Arg %i : %s\n', length(myvararg), char(myvararg(end))); end;
                end
            else
                % Regular chunk: Emit it
                ppos = strfind(curarg, '.lib');
                if ~isempty(ppos)
                    curarg = [ '-l' curarg(1:ppos(1)-1) ];
                end

                outarg = curarg;
                emitarg = 1;
            end
        end
    end

    if emitarg
        emitarg = 0;
        myvararg(end+1) = outarg;
        if (debugme), fprintf('Emitted Arg %i : %s\n', length(myvararg), char(myvararg(end))); end;
    end
end

final = myvararg

if (debugme), outargtype = class(myvararg), end;

if (callmex)
    mex(myvararg{:});
end

return;
