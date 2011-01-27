suffixes = {
        'T': 1e+12,
        'G': 1e+9,
        'M': 1e+6,
        'k': 1e+3,
        'c': 1e-2,
        'm': 1e-3,
        'u': 1e-6,
        'n': 1e-9,
        'p': 1e-12
}

def parse(v):
        """ Return a float from a suffixed string """
        if isinstance(v, str) and v[-1] in suffixes:
                return float(v[:-1]) * suffixes[v[-1]]
        return float(v)

def format(v, prec=8):
        """ Format a floating point number formatted with a suffix """
        # Dumb brute force
        for suf,val in suffixes.items():
                if v/val > 1 and v/val < 1000:
                        fmt = '%%1.%df%s' % (prec, suf)
                        return fmt % (v/val)

        # Out of range of our known suffixes, just use scientific notation
        fmt = '%%1.%de' % prec
        return fmt % v

