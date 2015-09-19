>hr2 folder has another readme for hr2[n-3].  
>For now, it is the best place to start.
>
>CHOFIT3_min is a supplement to my [paper](https://www.researchgate.net/publication/274720348_Fast_Graphically_Inspired_Algorithm_for_Assignment_of_Molecular_Formulae_in_Ultrahigh_Resolution_Mass_Spectrometry).
>
>Anal chem is romeo white (with 1 year embargo), so I haven't uploaded to ResearchGate yet.


**Details**

CHOFIT3 contains a faster algorithm for calculating the formula from exact mass for 
ultrahigh resolution mass spectrometry (< 1 ppm error).  

hr2 is similiar program, but written in c++ instead of freepascal.

2014-11-26 - User-NG: Created repository
2015-08-15 - ngreen-gt: Added CHOFIT3_min.pas, a stripped down version of CHOFIT3.  It is a supplement to the Analytical Chemistry paper.
2015-09-19 - ngreen-gt: Added hr2 folder.


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

