# HLTMacros
Generic Analysis Macros for CMS Triggers
Works in CMSSW_7_4_0_pre7

Contains:
  compileTrigMacro.sh: Compiles both C++ macros
  matchTrigTree_HI.C: Computes trigger rates, creates turnon curves
  plotTrigTurnOn_HI.C: plots turn on curves (up to five per trigger "type") 
  inputTrigFile.txt: example of input trigger list 

Use:
  Need a file with TTree of HLT triggers to be tested and HiForest of matching offline objects

To create exectuable of histogram maker:
  sh compileTrigMacro.sh matchTrigTree_HI.C

To execute histogram maker:
  ./matchTrigTree_HI.exe inputHLTFile inputMatchingHiForestFile inputTriggerTextFile

To create executable of plotting macro:
  sh compileTrigMacro.sh plotTrigTurnOn_HI.C

To execute plotter:
  ./plotTrigTurnOn.exe inputHistFile inputTriggerTextFile

Notes: 
  * Compiler automatically makes directory "pdfDir" since plotting macro writes pdf to said dir
  * plotter is limited to 5 turn on curves automatically since more than 5 stuff gets crowded. You can change if you like, but I advise against. The first 5 listed are chosen
  * More complicated "types" i.e. offline objects are possible. Should be straightforward to add with a couple of tweaks to code