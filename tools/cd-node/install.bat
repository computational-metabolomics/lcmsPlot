@echo off
rem Installs the Compounds Chromatograms node into Compound Discoverer.
rem Double-click this file; it asks for administrator rights, which are needed
rem to write into the Compound Discoverer folder under Program Files.
setlocal

net session >nul 2>&1
if %errorlevel% neq 0 (
    if "%~1"=="" (
        powershell -NoProfile -Command "Start-Process -FilePath '%~f0' -Verb RunAs"
    ) else (
        powershell -NoProfile -Command "Start-Process -FilePath '%~f0' -ArgumentList '%*' -Verb RunAs"
    )
    exit /b
)

powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0install.ps1" %*
set result=%errorlevel%
echo.
if %result% neq 0 (
    echo The installation did not complete. See install.log in this folder.
) else (
    echo The installation completed.
)
pause
exit /b %result%
