#! /usr/bin/env python

"""The **Wedring** pipeline setup script"""

from distutils.core import setup

from wedring import __version__

setup(name="Wedring",
      version=__version__,
      description="The Wedring pipeline for gene differential expression.",
      author="Rafael Covre",
      author_email="covrebio06@gmail.com",
      url="https://github.com/rcovre/Wedring",
      py_modules=["wedring.gffutils",
                  "wedring.indexstage",
                  "wedring.manager",
                  "wedring.mapstage",
                  "wedring.sysutils",
                  "wedring.tableutils",
                  "wedring.wedrerror"],
      data_files=[("bin", ["wedr", "diffExprStage.R"]),
                  ("share/doc/wedring", ["wedring.cfg"])])
