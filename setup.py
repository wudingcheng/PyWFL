from numpy.distutils.core import Extension
import glob
import os

wfl_include = 'C:\WFL\intel\include\intel64'
wfl_libarary = 'C:\WFL\intel\lib\intel64'


def get_submodules():
    f90 = glob.glob("*/*.f90")
    exts = []
    for f in f90:
        if "test" in f:
            continue
        # print f
        # exts.append(Extension('_%s'.format(os.path.split(f)[:-3])))
        name = "_{}".format(os.path.split(f)[1][:-4])
        exts.append(Extension(name, sources=[f],
                              library_dirs=[wfl_libarary],
                              include_dirs=[wfl_include],
                              f2py_options=['--quiet']))
    return exts


pywfl = Extension('_pywfl', sources=['pywfl.f90'],
                  library_dirs=[wfl_libarary],
                  include_dirs=[wfl_include],
                  f2py_options=['--quiet'])

# exts = get_submodules()

exts = [pywfl]

if __name__ == '__main__':
    from numpy.distutils.core import setup

    setup(name='PyWFL',
          ext_modules=exts)
