from numpy.distutils.core import Extension, setup
import glob
# ex = Extension(name='test',
#                sources=['filters.f90'] + glob.glob(r'C:\WFL\intel\include\intel64\*.f90'),
#                include_dirs=[r'C:\WFL\intel\include\intel64'],
#                library_dirs=[r'C:\WFL\intel\lib\intel64'],
#                libraries=['WHIGG']
#                )
ex = Extension(name='test',
               sources=['filters.f90'],
               include_dirs=[r'C:\WFL\intel\include\intel64'],
               library_dirs=[r'C:\WFL\intel\lib\intel64'],
               extra_f90_compile_args=[r'C:\WFL\intel\lib\intel64\whigg.a']
               )
if __name__ == '__main__':
    setup(name='example',
          ext_modules=[ex],)
