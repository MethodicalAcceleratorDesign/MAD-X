@echo off
setlocal

REM run:
REM scripts/build-test-win.bat [cleanall]

set CAT=c:\gnuwin32\bin\cat
set DATE=c:\gnuwin32\bin\date
set ECHO=c:\gnuwin32\bin\echo
set LS=c:\gnuwin32\bin\ls
set MAKE=c:\gnuwin32\bin\make
set RM=c:\gnuwin32\bin\rm
set UNAME=c:\gnuwin32\bin\uname

%rm% -f build-test-win.out
%uname% -n > build-test-win.run

%echo% -e "\n===== Start of build and tests ====="
%date%
%uname% -m -n -r -s

%echo% -e "\n===== SVN update ====="
"c:\Program Files\TortoiseSVN\bin\svn" update
if ERRORLEVEL 1 %echo% "ERROR: svn update failed"

%echo% -e "\n===== Release number ====="
%cat% VERSION

%echo% -e "\n===== Clean build ====="
if "%1"=="cleanall" (
   %make% cleanall && %make% cleanall ARCH=32
   if ERRORLEVEL 1 %echo% "ERROR: make cleanall failed"
) else (
   %echo% "Skipped (no explicit request)."
)

%echo% -e "\n===== Intel build ====="
call "C:\Program Files (x86)\Intel\Composer XE 2013 SP1\bin\ipsxe-comp-vars.bat" ia32 vs2012
%make% all-win32-intel all-win32
if ERRORLEVEL 1 %echo% "ERROR: make all-win32-intel failed"

call "C:\Program Files (x86)\Intel\Composer XE 2013 SP1\bin\ipsxe-comp-vars.bat" intel64 vs2012
%make% all-win64-intel all-win64
if ERRORLEVEL 1 %echo% "ERROR: make all-win64-intel failed"

%echo% -e "\n===== Binaries dependencies ====="
%make% infobindep
if ERRORLEVEL 1 %echo% "ERROR: make infobindep failed"

%echo% -e "\n===== Tests pointless files ====="
%make% cleantest && %make% infotestdep
if ERRORLEVEL 1 %echo% "ERROR: make infotestdep failed"

%echo% -e "\n===== Intel tests (64 bit) ====="
%make% madx-win64-intel && %ls% -l madx64.exe && %make% cleantest && %make% tests-all ARCH=64 NOCOLOR=yes
if ERRORLEVEL 1 %echo% "ERROR: make tests-all for madx-win64-intel failed"

%echo% -e "\n===== Intel tests (32 bit) ====="
%make% madx-win32-intel && %ls% -l madx32.exe && %make% cleantest && %make% tests-all ARCH=32 NOCOLOR=yes
if ERRORLEVEL 1 %echo% "ERROR: make tests-all for madx-win32-intel failed"

%echo% -e "\n===== End of build and tests ====="
%date%

%echo% "finished" > build-test-win.run
