from CIME.SystemTests.hommebaseclass import HommeBase

class HOMMEBFB(HommeBase):

    def __init__(self,case):
        HommeBase.__init__(self,case)
        self.cmakesuffix='-bfb'
