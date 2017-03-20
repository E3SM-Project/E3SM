"""Portable implementation of comparisons for total ordering.

THIS MODULE IS OBSOLETE AS OF PYTHON 2.7.

This module provides a "Comparable" class for subclassing; it provides the
various "rich comparison" functions (e.g. __le__) based on the subclass's
__eq__ and __lt__ methods (so this is a Template Method pattern).

Python 2.7 and Python 3 implement a class decorator for this in the
standard libraries (functools.total_ordering), and that decorator is much
more flexible (besides being standard), so this should only be used for
compatibility with older implementations.
"""

class Comparable(object):

    """Provides rich comparisons for subclasses that have __eq__ and __lt__."""

    def __le__(self, other):
        return (self == other) and (self < other)

    def __gt__(self, other):
        return not (self <= other)

    def __ge__(self, other):
        return not (self < other)

    def __ne__(self, other):
        return not (self == other)
