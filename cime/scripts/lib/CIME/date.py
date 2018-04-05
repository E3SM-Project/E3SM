import re
from CIME.XML.standard_module_setup import *
logger = logging.getLogger(__name__)
###############################################################################
def get_file_date(filename):
###############################################################################
    """
    Returns the date associated with the filename as a date object representing the correct date
    Formats supported:
    "%Y-%m-%d_%h.%M.%s
    "%Y-%m-%d_%05s"
    "%Y-%m-%d-%05s"
    "%Y-%m-%d"
    "%Y-%m"
    "%Y.%m"

    >>> get_file_date("./ne4np4_oQU240.cam.r.0001-01-06-00435.nc")
    date(1, 1, 6, 0, 7, 15)
    >>> get_file_date("./ne4np4_oQU240.cam.r.0010-1-06_00435.nc")
    date(10, 1, 6, 0, 7, 15)
    >>> get_file_date("./ne4np4_oQU240.cam.r.0010-10.nc")
    date(10, 10, 1, 0, 0, 0)
    >>> get_file_date("0064-3-8_10.20.30.nc")
    date(64, 3, 8, 10, 20, 30)
    >>> get_file_date("0140-3-5")
    date(140, 3, 5, 0, 0, 0)
    >>> get_file_date("0140-3")
    date(140, 3, 1, 0, 0, 0)
    >>> get_file_date("0140.3")
    date(140, 3, 1, 0, 0, 0)
    """

    #
    # TODO: Add these to config_archive.xml, instead of here
    # Note these must be in order of most specific to least
    # so that lesser specificities aren't used to parse greater ones
    re_formats = [r"[0-9]*[0-9]{4}-[0-9]{1,2}-[0-9]{1,2}_[0-9]{1,2}\.[0-9]{1,2}\.[0-9]{1,2}", # [yy...]yyyy-mm-dd_hh.MM.ss
                  r"[0-9]*[0-9]{4}-[0-9]{1,2}-[0-9]{1,2}[\-_][0-9]{1,5}",                     # [yy...]yyyy-mm-dd_sssss
                  r"[0-9]*[0-9]{4}-[0-9]{1,2}-[0-9]{1,2}",                                    # [yy...]yyyy-mm-dd
                  r"[0-9]*[0-9]{4}[\-\.][0-9]{1,2}",                                          # [yy...]yyyy-mm
    ]

    for re_str in re_formats:
        match = re.search(re_str, filename)
        if match is None:
            continue
        date_str = match.group()
        date_tuple = [int(unit) for unit in re.split(r"-|_|\.", date_str)]
        year = date_tuple[0]
        month = date_tuple[1]
        day = 1
        second = 0
        if len(date_tuple) > 2:
            day = date_tuple[2]
            if len(date_tuple) == 4:
                second = date_tuple[3]
            elif len(date_tuple) == 6:
                # Create a date object with arbitrary year, month, day, but the correct time of day
                # Then use _get_day_second to get the time of day in seconds
                second = date.hms_to_second(hour = date_tuple[3],
                                            minute = date_tuple[4],
                                            second = date_tuple[5])
        return date(year, month, day, 0, 0, second)

    # Not a valid filename date format
    logger.debug("{} is a filename without a supported date!".format(filename))
    return None

class date:
    """
    Simple struct for holding dates and the time of day and performing comparisons

    Difference in Hour, Minute, or Second
    >>> date(4, 5, 6, 9) == date(4, 5, 6, 8)
    False
    >>> date(4, 5, 6, 9) != date(4, 5, 6, 8)
    True
    >>> date(4, 5, 6, 9) < date(4, 5, 6, 8)
    False
    >>> date(4, 5, 6, 9) <= date(4, 5, 6, 8)
    False
    >>> date(4, 5, 6, 9) >= date(4, 5, 6, 8)
    True
    >>> date(4, 5, 6, 9) > date(4, 5, 6, 8)
    True

    >>> date(4, 5, 6, 4) == date(4, 5, 6, 8)
    False
    >>> date(4, 5, 6, 4) != date(4, 5, 6, 8)
    True
    >>> date(4, 5, 6, 4) < date(4, 5, 6, 8)
    True
    >>> date(4, 5, 6, 4) <= date(4, 5, 6, 8)
    True
    >>> date(4, 5, 6, 4) >= date(4, 5, 6, 8)
    False
    >>> date(4, 5, 6, 4) > date(4, 5, 6, 8)
    False

    Difference in Day
    >>> date(4, 5, 8, 8) == date(4, 5, 6, 8)
    False
    >>> date(4, 5, 8, 8) != date(4, 5, 6, 8)
    True
    >>> date(4, 5, 8, 8) < date(4, 5, 6, 8)
    False
    >>> date(4, 5, 8, 8) <= date(4, 5, 6, 8)
    False
    >>> date(4, 5, 8, 8) >= date(4, 5, 6, 8)
    True
    >>> date(4, 5, 8, 8) > date(4, 5, 6, 8)
    True

    >>> date(4, 5, 5, 8) == date(4, 5, 6, 8)
    False
    >>> date(4, 5, 5, 8) != date(4, 5, 6, 8)
    True
    >>> date(4, 5, 5, 8) < date(4, 5, 6, 8)
    True
    >>> date(4, 5, 5, 8) <= date(4, 5, 6, 8)
    True
    >>> date(4, 5, 5, 8) >= date(4, 5, 6, 8)
    False
    >>> date(4, 5, 5, 8) > date(4, 5, 6, 8)
    False

    Difference in Month
    >>> date(4, 6, 6, 8) == date(4, 5, 6, 8)
    False
    >>> date(4, 6, 6, 8) != date(4, 5, 6, 8)
    True
    >>> date(4, 6, 6, 8) < date(4, 5, 6, 8)
    False
    >>> date(4, 6, 6, 8) <= date(4, 5, 6, 8)
    False
    >>> date(4, 6, 6, 8) >= date(4, 5, 6, 8)
    True
    >>> date(4, 6, 6, 8) > date(4, 5, 6, 8)
    True

    >>> date(4, 4, 6, 8) == date(4, 5, 6, 8)
    False
    >>> date(4, 4, 6, 8) != date(4, 5, 6, 8)
    True
    >>> date(4, 4, 6, 8) < date(4, 5, 6, 8)
    True
    >>> date(4, 4, 6, 8) <= date(4, 5, 6, 8)
    True
    >>> date(4, 4, 6, 8) >= date(4, 5, 6, 8)
    False
    >>> date(4, 4, 6, 8) > date(4, 5, 6, 8)
    False

    Difference in Year
    >>> date(5, 5, 6, 8) == date(4, 5, 6, 8)
    False
    >>> date(5, 5, 6, 8) != date(4, 5, 6, 8)
    True
    >>> date(5, 5, 6, 8) < date(4, 5, 6, 8)
    False
    >>> date(5, 5, 6, 8) <= date(4, 5, 6, 8)
    False
    >>> date(5, 5, 6, 8) >= date(4, 5, 6, 8)
    True
    >>> date(5, 5, 6, 8) > date(4, 5, 6, 8)
    True

    >>> date(3, 5, 6, 8) == date(4, 5, 6, 8)
    False
    >>> date(3, 5, 6, 8) != date(4, 5, 6, 8)
    True
    >>> date(3, 5, 6, 8) < date(4, 5, 6, 8)
    True
    >>> date(3, 5, 6, 8) <= date(4, 5, 6, 8)
    True
    >>> date(3, 5, 6, 8) >= date(4, 5, 6, 8)
    False
    >>> date(3, 5, 6, 8) > date(4, 5, 6, 8)
    False
    """
    @staticmethod
    def hms_to_second(hour, minute, second):
        _SECONDS_PER_HOUR = 3600
        _SECONDS_PER_MINUTE = 60
        return (hour * _SECONDS_PER_HOUR + minute * _SECONDS_PER_MINUTE
                + second)

    @staticmethod
    def second_to_hms(second):
        _SECONDS_PER_HOUR = 3600
        _SECONDS_PER_MINUTE = 60
        return { 'hour': second // _SECONDS_PER_HOUR,
                 'minute': (second % _SECONDS_PER_HOUR) // _SECONDS_PER_MINUTE,
                 'second': second % _SECONDS_PER_MINUTE
        }

    def __init__(self, year=1, month=1, day=1, hour=0, minute=0, second=0):
        self._year = year
        self._month = month
        self._day = day
        self._second = self.hms_to_second(hour, minute, second)

    def __str__(self):
        """
        >>> str(date(4, 5, 7, second=64))
        'date(4, 5, 7, 0, 1, 4)'
        """
        fmt_str = "date({year:d}, {month:d}, {day:d}, {hour:d}, {minute:d}, {second:d})"
        return fmt_str.format(year = self.year(),
                              month = self.month(),
                              day = self.day(),
                              hour = self.hour(),
                              minute = self.minute(),
                              second = self.second())

    def year(self):
        return self._year

    def month(self):
        return self._month

    def day(self):
        return self._day

    def hour(self):
        return self.second_to_hms(self._second)['hour']

    def minute(self):
        return self.second_to_hms(self._second)['minute']

    def second(self):
        return self.second_to_hms(self._second)['second']

    def second_of_day(self):
        return self._second

    def __repr__(self):
        return str(self)

    def __eq__(self, other):
        return ((self.year() == other.year()) and (self.month() == other.month())
                and (self.day() == other.day())
                and (self.second_of_day() == other.second_of_day()))

    def __ne__(self, other):
        return not (self == other)

    def __lt__(self, other):
        if self.year() < other.year():
            return True
        elif self.year() > other.year():
            return False
        # self.year == other.year
        if self.month() < other.month():
            return True
        elif self.month() > other.month():
            return False
        # self.month = other.month
        if self.day() < other.day():
            return True
        elif self.day() > other.day():
            return False
        # self.day = other.day
        if self.second_of_day() < other.second_of_day():
            return True
        else:
            # the dates are equal
            return False

    def __le__(self, other):
        return ((self < other) or (self == other))

    def __ge__(self, other):
        return not (self < other)

    def __gt__(self, other):
        return not (self <= other)
