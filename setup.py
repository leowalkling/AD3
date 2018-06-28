import sys
from setuptools.command.bdist_egg import bdist_egg
from setuptools.extension import Extension
from skbuild import setup
import re

# Optionally retrieve the PyTorch installation info
try:
    import torch
    from torch.utils import cpp_extension

    def _s(s):
        match = re.fullmatch(r"([^:]*\:)(.+)", s)
        if " " in s and match is not None:
            return match[1] + '"' + match[2] + '"'
        else:
            return s

    cmake_args = [
        "-DPYTORCH_INCLUDES:PATH=\"" + ";".join(cpp_extension.include_paths(torch.cuda.is_available())) + '"',
        "-DPYTORCH_LDFLAGS:STRING=" + " ".join(map(_s, cpp_extension._prepare_ldflags([],
                                                                                      torch.cuda.is_available(),
                                                                                      False))),
    ]
except ImportError:
    cmake_args = []


AD3_COMPILE_ARGS = [
    '-fPIC',
    '-O3',
    '-c',
    '-fmessage-length=0'
]


libad3 = ('ad3', {
    'sources': ['ad3/FactorGraph.cpp',
                'ad3/GenericFactor.cpp',
                'ad3/Factor.cpp',
                'ad3/Utils.cpp',
                'examples/cpp/parsing/FactorTree.cpp'
                ],
    'include_dirs': ['.',
                     './ad3',
                     './Eigen',
                     './examples/cpp/parsing'
                     ],
    'extra_compile_args': AD3_COMPILE_ARGS
})

# this is a backport of a workaround for a problem in distutils.
# install_lib doesn't call build_clib
#
#
# class bdist_egg_fix(bdist_egg):
#     def run(self):
#         self.call_command('build_clib')
#         bdist_egg.run(self)


WHEELHOUSE_UPLOADER_COMMANDS = set(['fetch_artifacts', 'upload_all'])

cmdclass = {}
# cmdclass = {'bdist_egg': bdist_egg_fix}

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
      packages=['ad3', 'ad3.tests'],
      # libraries=[libad3],
      # cmdclass=cmdclass,
      include_package_data=True,
      cmake_args=cmake_args
)
