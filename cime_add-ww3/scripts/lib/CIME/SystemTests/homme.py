from CIME.SystemTests.hommebaseclass import HommeBase

import shutil
from distutils import dir_util

class HOMME(HommeBase):

    def __init__(self,case):
       HommeBase.__init__(self,case)
       self.cmakesuffix=''



