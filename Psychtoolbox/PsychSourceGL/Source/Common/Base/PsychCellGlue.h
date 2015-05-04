/*
  PsychToolbox2/Source/Common/PsychCellGlue.h		
  
  AUTHORS:
  Allen.Ingling@nyu.edu		awi 
  
  PLATFORMS: All
  
  PROJECTS:
  12/31/02	awi		Screen on OS X
   

  HISTORY:
  12/17/03  awi		wrote it.  
  
  DESCRIPTION:
  
	PsychStructGlue defines abstracted functions to create cell arrays passed 
	between the calling environment and the PsychToolbox. 
  
  NOTES:
  
  
  
*/



//begin include once 
#ifndef PSYCH_IS_INCLUDED_CellGlue
#define PSYCH_IS_INCLUDED_CellGlue

#include "Psych.h"

psych_bool PsychAllocOutCellVector(int position, PsychArgRequirementType isRequired, int numElements,  PsychGenericScriptType **pCell);
void PsychSetCellVectorStringElement(int index, const char *text, PsychGenericScriptType *cellVector);
void PsychSetCellVectorDoubleElement(int index, double value, PsychGenericScriptType *cellVector);
void PsychSetCellVectorNativeElement(int index, PsychGenericScriptType *pNativeElement,  PsychGenericScriptType *cellVector);
psych_bool PsychAllocInNativeCellVector(int position, PsychArgRequirementType isRequired, const PsychGenericScriptType **cellVector);
psych_bool PsychAllocInNativeString(int position, PsychArgRequirementType isRequired, const PsychGenericScriptType **nativeString);
void PsychConvertNativeCellArrayToNativeString(const PsychGenericScriptType **nativeCellArray, PsychGenericScriptType **nativeString);
         
                                                                             

//end include once
#endif


