# Copyright (C) 2010
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

def GCreadlines(fobject, comment='#'):
    """Strip files from comments.

    fobject: file object to be read
    comment: string indicating comment.
    """

    lines = []
    for l in fobject.readlines():
        # ignore comments and empty lines
        l = l.strip()
        pos = l.find(comment)
        if pos>-1:
            l = l[:pos]
        if len(l)==0:
            continue
        lines.append(l)
    return lines
