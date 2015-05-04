/*
 * PortAudio Portable Real-Time Audio Library
 * Macintosh Core Audio specific extensions
 * portaudio.h should be included before this file.
 *
 * Copyright (c) 2005-2006 Bjorn Roche
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files
 * (the "Software"), to deal in the Software without restriction,
 * including without limitation the rights to use, copy, modify, merge,
 * publish, distribute, sublicense, and/or sell copies of the Software,
 * and to permit persons to whom the Software is furnished to do so,
 * subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR
 * ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
 * CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

/*
 * The text above constitutes the entire PortAudio license; however, 
 * the PortAudio community also makes the following non-binding requests:
 *
 * Any person wishing to distribute modifications to the Software is
 * requested to send the modifications to the original developer so that
 * they can be incorporated into the canonical version. It is also 
 * requested that these non-binding requests be included along with the 
 * license above.
 */

#ifdef __cplusplus
extern "C" {
#endif


/*
 * A pointer to a paMacCoreStreamInfo may be passed as
 * the hostApiSpecificStreamInfo in the PaStreamParameters struct
 * when opening a stream or querying the format. Use NULL, for the
 * defaults. Note that for duplex streams, both infos should be the
 * same or behaviour is undefined.
 */
typedef struct paMacCoreStreamInfo
{
    unsigned long size;         /**size of whole structure including this header */
    PaHostApiTypeId hostApiType;/**host API for which this data is intended */
    unsigned long version;      /**structure version */
    unsigned long flags;        /* flags to modify behaviour */
} paMacCoreStreamInfo;

/* MK: Disabled function body here in PsychPortAudio's private version of pa_mac_core.h
 * to avoid conflicting duplicate declaration of this function and symbols in
 * PsychPortAudio.o and the Portaudio library. However, we need to keep/define the
 * function prototype to avoid compiler warnings.
 */
/* Use this function to initialize a paMacCoreStreamInfo struct using the requested flags. */
void paSetupMacCoreStreamInfo( paMacCoreStreamInfo *data, unsigned long flags );
/*
{
   bzero( data, sizeof( paMacCoreStreamInfo ) );
   data->size = sizeof( paMacCoreStreamInfo );
   data->hostApiType = paCoreAudio;
   data->version = 0x01;
   data->flags = flags;
}
*/

/*
 * The following flags alter the behaviour of PA on the mac platform.
 * they can be ORed together. These should work both for opening and
 * checking a device.
 */
/* Allows PortAudio to change things like the device's frame size,
 * which allows for much lower latency, but might disrupt the device
 * if other programs are using it, even when you are just Querying
 * the device. */
#define paMacCore_ChangeDeviceParameters       0x01

/* In combination with the above flag,
 * causes the stream opening to fail, unless the exact sample rates
 * are supported by the device. */
#define paMacCore_FailIfConversionRequired     0x02

/* These flags set the SR conversion quality, if required. The wierd ordering
 * allows Maximum Quality to be the default.*/
#define paMacCore_ConversionQualityMin     0x0100
#define paMacCore_ConversionQualityMedium  0x0200
#define paMacCore_ConversionQualityLow     0x0300
#define paMacCore_ConversionQualityHigh    0x0400
#define paMacCore_ConversionQualityMax     0x0000

/*
 * Here are some "preset" combinations of flags (above) to get to some
 * common configurations. THIS IS OVERKILL, but if more flags are added
 * it won't be.
 */
/*This is the default setting: do as much sample rate conversion as possible
 * and as little mucking with the device as possible. */
#define paMacCorePlayNice  0x00
/*This setting is tuned for pro audio apps. It allows SR conversion on input
  and output, but it tries to set the appropriate SR on the device.*/
#define paMacCorePro       0x01
/*This is a setting to minimize CPU usage and still play nice.*/
#define paMacCoreMinimizeCPUButPlayNice  0x0100
/*This is a setting to minimize CPU usage, even if that means interrupting the device. */
#define paMacCoreMinimizeCPU  0x0101


#ifdef __cplusplus
}
#endif /* __cplusplus */