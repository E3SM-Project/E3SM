# Copyright (C) 2004, 2005, 2009, 2010
# Glimmer-CISM contributors - see AUTHORS file for list of contributors
#
# This file is part of Glimmer-CISM.
#
# Glimmer-CISM is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or (at
# your option) any later version.
#
# Glimmer-CISM is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Glimmer-CISM.  If not, see <http://www.gnu.org/licenses/>.
#
# Glimmer-CISM is hosted on BerliOS.de:
# https://developer.berlios.de/projects/glimmer-cism/

__all__ = ['GCEpoch']

import numpy, math, os.path
from GC_readlines import GCreadlines

class GCEpoch(object):
    """Handle epochs."""

    def __init__(self,fname=None,sharedir='data'):
        """Initialise epoch data from file.

        fname: name of file containing epoch data.

        file format: comments start with #
        each row contains 4 columns separate by commas: name, start time, end time and GMT RGB string."""

        if fname=='None':
            fn = os.path.join(sharedir,'stages')
        else:
            fn = fname

        self.data = []
        self.timescale = 0.001
        f = open(fn,'r')
        for l in GCreadlines(f):
            l = l.split(',')
            colour = l[3].strip().split('/')
            assert len(colour)==3
            for i in range(0,3):
                colour[i] = float(colour[i])/255.
            self.data.append({'name' : l[0], 'start':float(l[1]), 'end':float(l[2]),'colour':colour})
        f.close()

        self.__current = 0

    def get_epoch(self,t):
        """Return the name of the epoch given a time.

        t: time in ka"""

        time = t/self.timescale

        # first try if current points to the right epoch
        if self.data[self.__current]['start']<=time and self.data[self.__current]['end']>=time:
            return self.__current

        for c in range(0,len(self.data)):
            if self.data[c]['start']<=time and self.data[c]['end']>=time:
                self.__current = c
                return self.__current

        return None

    def get_colour(self,t):
        """Return the RGB string associated with time.

        t: time in ka"""

        c = self.get_epoch(t)
        if c != None:
            return self.data[c]['colour']
        else:
            return [1.,1.,1.]
