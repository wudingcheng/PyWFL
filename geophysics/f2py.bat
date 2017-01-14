rem set if_inc=C:\WFL\intel\include\intel64
rem set if_a=C:\WFL\intel\lib\intel64\whigg.a
python f2py.py -c --fcompiler=intelvem --compiler=msvc -m _geophysics geophysics.f90 -I%if_inc% %if_a% --build-dir build --quiet

echo Delete build directory
rd /s build