@echo off
call setEnv.bat
"C:\Program Files\MATLAB\R2020a\toolbox\shared\coder\ninja\win64\ninja.exe" -v %*
