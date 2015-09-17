@echo off
rem Change directory and drive to batch file location
cd /d %~dp0

rem Run HR2 programs sequentially and output elapsed time

rem HR2 [n]
	rem Get start time:
	for /F "tokens=1-4 delims=:.," %%a in ("%time%") do (
	   set /A "start=(((%%a*60)+1%%b %% 100)*60+1%%c %% 100)*100+1%%d %% 100"
	)

	rem Any process here...
	programs\Test_HR2_n.exe -t 0.2 -C 1-84 -H 0-144 -O 0-37 -N 0-10 -S 0-6 -P 0-4 -1 0-1 mass_list\allCHO_HR2_short.dat > output\allCHO_HR2_short_n.out
	
	rem Get end time:
	for /F "tokens=1-4 delims=:.," %%a in ("%time%") do (
	   set /A "end=(((%%a*60)+1%%b %% 100)*60+1%%c %% 100)*100+1%%d %% 100"
	)

	rem Get elapsed time:
	set /A elapsed=end-start

	rem Show elapsed time:
	set /A hh=elapsed/(60*60*100), rest=elapsed%%(60*60*100), mm=rest/(60*100), rest%%=60*100, ss=rest/100, cc=rest%%100
	if %mm% lss 10 set mm=0%mm%
	if %ss% lss 10 set ss=0%ss%
	if %cc% lss 10 set cc=0%cc%
	echo HR2 [n] ELAPSED TIME: %hh%:%mm%:%ss%,%cc% > output\HR2_results.txt

rem HR2 [n-3]
	rem Get start time:
	for /F "tokens=1-4 delims=:.," %%a in ("%time%") do (
	   set /A "start=(((%%a*60)+1%%b %% 100)*60+1%%c %% 100)*100+1%%d %% 100"
	)

	rem Any process here...
	programs\Test_HR2_n3.exe -t 0.2 -C 1-84 -H 0-144 -O 0-37 -N 0-10 -S 0-6 -P 0-4 -1 0-1 mass_list\allCHO_HR2_short.dat > output\allCHO_HR2_short_n3.out
	
	rem Get end time:
	for /F "tokens=1-4 delims=:.," %%a in ("%time%") do (
	   set /A "end=(((%%a*60)+1%%b %% 100)*60+1%%c %% 100)*100+1%%d %% 100"
	)

	rem Get elapsed time:
	set /A elapsed=end-start

	rem Show elapsed time:
	set /A hh=elapsed/(60*60*100), rest=elapsed%%(60*60*100), mm=rest/(60*100), rest%%=60*100, ss=rest/100, cc=rest%%100
	if %mm% lss 10 set mm=0%mm%
	if %ss% lss 10 set ss=0%ss%
	if %cc% lss 10 set cc=0%cc%
	echo HR2 [n-3] ELAPSED TIME: %hh%:%mm%:%ss%,%cc% >> output\HR2_results.txt


rem View elapsed times in CHOFIT_results.txt 	
rem End of batch file

	