@echo off
cls

REM .paket\paket.bootstrapper.exe prerelease
REM if errorlevel 1 (
REM   exit /b %errorlevel%
REM )

.paket\paket.exe restore
if errorlevel 1 (
  exit /b %errorlevel%
)

dotnet fake run build.fsx -t %*