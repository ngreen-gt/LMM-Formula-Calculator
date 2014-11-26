LMM Formula Calculator
======================

Formula Calculator for accurate masses using a low-mass moiety algorithm.  Developed in C/C++, the source code was originally the HR2 program used in the Seven Golden Rules paper.  

2014-11-26  -  NG: Created repository.  I'm new on github, and learning.  I'll update here with news for now.  


Known Issues
===============
The charge correction to mass is turned off, and the electron won't be added/subtracted.  The current working version only uses exact masses, not ion/radical masses.  This is an integration issue with mass list, that can be fixed/customized, but is a low priority.

In the core calculation, the program stops searching for CHO formulae once the first one is found.  This may cause problems when large tolerances are used, and not find all formulae.  The minimum difference between two CHO formulae is 1.1 mDa (or ppm equivalent), which is when problems would begin.  The best current idea to fix this would using a less optimized program at higher tolerance, which searches through all potential CHO formula, and can return two CHO formulae.
