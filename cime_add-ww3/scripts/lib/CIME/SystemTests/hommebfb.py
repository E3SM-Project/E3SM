from CIME.SystemTests.hommebaseclass import HommeBase

import shutil
from distutils import dir_util

class HOMMEBFB(HommeBase):

    def __init__(self,case):
       HommeBase.__init__(self,case)
       self.cmakesuffix='-bfb'



