#!/usr/bin/python

import json
import os

patterns = [
        '.%s.metadata',
        '%s.metadata',
        '.%s.meta',
        '%s.meta',
]

def get_metadata(f):
        """ Returns the parameters associated with the given file """
        fname = f if isinstance(f, str) else f.name
        for p in patterns:
                name = p % fname
                if os.path.isfile(name):
                        metad = json.load(open(name))
                        return metad

        return None

