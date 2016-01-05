@echo off
setlocal

REM run:
REM scripts/build-test-win.bat [noecho] [cleanall] [notest]

REM commands
set CAT="c:\gnuwin32\bin\cat"
set DATE="c:\gnuwin32\bin\date"
set ECHO="c:\gnuwin32\bin\echo"
set TEE="c:\gnuwin32\bin\tee"
set LS="c:\gnuwin32\bin\ls"
set MAKE="c:\gnuwin32\bin\make"
set RM="c:\gnuwin32\bin\rm"
set UNAME="c:\gnuwin32\bin\uname"
set SCP="c:\msys\bin\scp"
set GCC="c:\mingw64\bin\gcc"
set GCXX="c:\mingw64\bin\g++"
set GFC="c:\mingw64\bin\gfortran"
set SVN="c:\Program Files\TortoiseSVN\bin\svn"

REM settings
set SSHRSA="c:/users/mad/.ssh/id_rsa"
set MACDIR="mad@macserv15865.cern.ch:Projects/madX"

if "%1"=="dont-redirect" shift & goto next
%rm% -f build-test-win.out
if "%1"=="noecho" goto noecho
call scripts\build-test-win.bat dont-redirect %* 2>&1 | %tee% build-test-win.out
goto doscp
:noecho
call scripts\build-test-win.bat dont-redirect %* > build-test-win.out 2>&1
:doscp
%scp% -q -i %sshrsa% build-test-win.out *-win32.exe *-win64.exe %macdir%
exit /B

:next
if "%1"=="noecho" shift

%echo% -e "\n===== Start of build and tests ====="
%date%
%uname% -m -n -r -s

%echo% -e "\n===== SVN update ====="
%svn% update
if ERRORLEVEL 1 (
	%echo% -e "\n===== SVN cleanup & update ====="
	%svn% cleanup
	%svn% update
	%echo% -e "\nERROR: svn update failed" && exit /B 1
)

%echo% -e "\n===== Release number ====="
%cat% VERSION

%echo% -e "\n===== Clean build ====="
if "%1"=="cleanall" (
REM cleanall not supported on windows (relies on find)
  shift
  %make% cleanbuild
) else (
  %echo% "Skipped (no explicit request)."
)

%echo% -e "\n===== Gnu build ====="
%gcc%  --version
%gcxx% --version 
%gfc%  --version
%make% all-win64-gnu
if ERRORLEVEL 1 %echo% -e "\nERROR: make all-win64-gnu failed" && exit /B 1

%echo% -e "\n===== Intel build ====="
call "C:\Program Files (x86)\IntelSWTools\compilers_and_libraries\windows\bin\compilervars.bat" ia32 vs2013
%make% all-win32-intel all-win32
if ERRORLEVEL 1 %echo% -e "\nERROR: make all-win32-intel failed" && exit /B 1

call "C:\Program Files (x86)\IntelSWTools\compilers_and_libraries\windows\bin\compilervars.bat" intel64 vs2013
%make% all-win64-intel all-win64
if ERRORLEVEL 1 %echo% -e "\nERROR: make all-win64-intel failed" && exit /B 1

%echo% -e "\n===== Binaries dependencies ====="
%make% infobindep
if ERRORLEVEL 1 %echo% -e "\nERROR: make infobindep failed" && exit /B 1

%echo% -e "\n===== Tests pointless files ====="
%make% cleantest && %make% infotestdep
if ERRORLEVEL 1 %echo% -e "\nERROR: make infotestdep failed" && exit /B 1

%echo% -e "\n===== Running tests (long) ====="
if "%1"=="notest" (
	shift
    %echo% "Skipped (explicit request)."
) else (
	%echo% ""

	%echo% -e "\n===== Testing madx-win64-intel ====="
	%make% madx-win64-intel && %ls% -l madx-win64-intel.exe madx64.exe && %make% cleantest && %make% tests-all COMP=intel ARCH=64 NOCOLOR=yes
	if ERRORLEVEL 1 %echo% -e "\nERROR: make tests-all for madx-win64-intel failed" && exit /B 1

	%echo% -e "\n===== Testing madx-win32-intel ====="
	%make% madx-win32-intel && %ls% -l madx-win32-intel.exe madx32.exe && %make% cleantest && %make% tests-all COMP=intel ARCH=32 NOCOLOR=yes
	if ERRORLEVEL 1 %echo% -e "\nERROR: make tests-all for madx-win32-intel failed" && exit /B 1

	%echo% -e "\n===== Testing madx-win64-gnu ====="
	set GFORTRAN_UNBUFFERED_PRECONNECTED=y
	%make% madx-win64-gnu && %ls% -l madx-win64-gnu.exe madx64.exe && %make% cleantest && %make% tests-all COMP=gnu ARCH=64 NOCOLOR=yes
	if ERRORLEVEL 1 %echo% -e "\nERROR: make tests-all for madx-win64-intel failed" && exit /B 1

	REM restore the default version
	%make% madx-win32 > tmp.out && %make% madx-win64 > tmp.out && %rm% -f tmp.out
	if ERRORLEVEL 1 %echo% -e "\nERROR: unable to restore the default version" && exit /B 1
)

REM date & end marker
%echo% -e "\nFinish:"
%date%
%echo% -e "\n===== End of build and tests ====="
