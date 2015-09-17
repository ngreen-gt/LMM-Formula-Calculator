Overview/Updates/Quick Start
============================================================================================
LMM formula calculatorsCHOFIT3 contains a fast algorithm for calculating the formula from exact mass for mass spectrometry (< 1 ppm error).  It is based of low-mass moiety (LMM).  Further background below.  Future project with he LMM algorithm involve incorporating it into mass spectra analysis software.     this code involves incorporation into other Example programs      Projects releated to the low-mass based formula assignment algorithm.  Projects in analysis of mass spectrometry for high/ultrahigh mass resolution.   formula assignment.    formulae assignment.  More background below.  Formula assignment Manage the various formula assingment tasks here.  Background about it is below.  For work related to formula assignment.  Include formula calculators, heuristics, analysis, This is a repository for all things low-mass moiety and formulae assingnment.  Right now it is home to simple formula calculator.  

2014-11-26 - User-NG: Created repository
2015-08-15 - ngreen-gt: Added CHOFIT3_min.pas, a stripped down version of CHOFIT3.  It is a supplement to the Analytical Chemistry paper.

CHOFIT3_min.pas (Version 1)
***************
Description: Formula calculator. *This version has limited user settings. 
ex. Mass List (ex. 590.232 Da, 592.343 Da)  >>> Formulae (ex. C20 H28 O16, C24 H28)
1) Compiled in FreePascal (newest compiler at FreePascal.org)
2) Non-GUI version.  Syntax for command prompt (Win) "CHOFIT3_min <input name> <output name> <low mass> <high mass> <max nitrogen> <max sulfur> <max phosphorus> <max carbon-13>"  *The input file must be .dat and output files are automatically .out
3) Input text file contain one mass per line.  Mass must be ion mass (charge = -1).
3) Limited user settings (ex. error = 0.4 ppm, charge = -1). Further details are in .pas file.


Background
======================================================================

