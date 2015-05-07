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
