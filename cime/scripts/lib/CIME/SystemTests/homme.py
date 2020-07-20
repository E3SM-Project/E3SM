from CIME.SystemTests.hommebaseclass import HommeBase

class HOMME(HommeBase):

    def __init__(self,case):
        HommeBase.__init__(self,case)
        self.cmakesuffix=''
