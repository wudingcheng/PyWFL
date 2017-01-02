set if_inc=C:\WFL\intel\include\intel64
set if_a=C:\WFL\intel\lib\intel64\whigg.a
set if_lib=C:\WFL\intel\lib\intel64
rem python f2py.py -c --fcompiler=intelvem --compiler=msvc -m _pywfl pywfl.f90 -I%if_inc% %if_a% -L%if_lib% --build-dir build
python f2py.py -c --fcompiler=intelvem --compiler=msvc -m _statistics statistics.f90 -I%if_inc% %if_a% -L%if_lib% --build-dir build

echo Delete build directory
rd /s build