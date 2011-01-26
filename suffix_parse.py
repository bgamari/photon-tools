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
        if v.__class__ is str and v[-1] in suffixes:
                return float(v[:-1]) * suffixes[v[-1]]
        return float(v)

