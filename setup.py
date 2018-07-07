import sys
from setuptools.command.bdist_egg import bdist_egg
from setuptools.extension import Extension
from skbuild import setup
import re


WHEELHOUSE_UPLOADER_COMMANDS = set(['fetch_artifacts', 'upload_all'])

cmdclass = {}

if WHEELHOUSE_UPLOADER_COMMANDS.intersection(sys.argv):

    import wheelhouse_uploader.cmd
    cmdclass.update(vars(wheelhouse_uploader.cmd))


setup(name='ad3',
      version="2.3.dev0",
      author="Andre Martins",
      description='Alternating Directions Dual Decomposition',
      url="http://www.ark.cs.cmu.edu/AD3",
      author_email="afm@cs.cmu.edu",
      package_dir={
          'ad3': 'python/ad3',
          'ad3.tests': 'python/ad3/tests'
      },
      package_data={"ad3": ["*.lib", "*.a", "*.pxd"]},
      packages=['ad3', 'ad3.tests'],
      include_package_data=True,
      zip_safe=False
)
