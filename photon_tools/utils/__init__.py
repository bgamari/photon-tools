import numpy as np

def parse_int_list(s):
    """
    Parse a simple list of integers and ranges.

    :type s: :class:`str`
    :param s: A textual list of integers, e.g. ``1,3,4,5-9,17``
    """

    def parse_range(r):
        i = r.find('-')
        if i != -1:
            a = int(r[:i])
            b = int(r[i+1:])
            if a > b:
                raise RuntimeError('lower range bound must be less than or equal to upper bound')
            return range(a,b+1)
        else:
            return [int(r)]

    return [i for r in s.split(',') for i in parse_range(r)]

def parse_intervals(s):
    """
    Parse a list of intervals.

    :type s: :class:`str`
    :param s: A textual list of intervals, e.g. ``-3,5-9,15,20-``
    :returns: An array of (start,end) bounds (both inclusive)
    """
    def parse_interval(s):
        parts = s.split('-')
        if len(parts) == 1:
            start = end = float(parts[0])
        elif len(parts) == 2:
            start, end = parts
            start = None if start == '' else float(start)
            end = None if end == '' else float(end)
        else:
            raise RuntimeError('Expected interval (e.g. "5-10" or "80-")')
        return (start, end)

    return [parse_interval(a) for a in s.split(',')]

def in_interval(bounds, arr):
    """
    Return an array identifying those array elements that fall within the given
    interval.

    :param bounds: a tuple of lower and upper bounds; ``None`` indicates unbounded.
    :param arr: an array of times
    :returns: a boolean array indicating which elements fall in the given time range.
    """
    lower, upper = bounds
    mask = True
    if lower is not None:
        mask = np.logical_and(mask, arr >= lower)
    if upper is not None:
        mask = np.logical_and(mask, arr <= upper)
    return mask

def in_intervals(bounds, arr):
    """
    Return an array identifying those array elements that fall within one or more
    of the given intervals.

    :type bounds: array of tuples
    :param bounds: See :func:`in_interval`
    :param arr: an array of times
    :returns: a boolean array indicating which elements fall in one or more of the given intervals
    """
    a = np.zeros_like(arr, dtype=bool)
    for b in bounds:
        np.logical_or(a, in_interval(b, arr), a)
    return a
