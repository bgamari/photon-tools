import numpy as np
from photon_tools import timetag_parse, pt2_parse, metadata

def determine_filetype(fname):
    if fname.endswith('pt2'):       return 'pt2'
    elif fname.endswith('pt3'):     return 'pt2'
    elif fname.endswith('timetag'): return 'timetag'
    elif fname.endswith('times'):   return 'raw'
    elif fname.endswith('timech'):     return 'timech'
    raise RuntimeError("Unrecognized file type")
    
def verify_monotonic(times):
    """ Verify that timestamps are monotonically increasing """
    negatives = times[1:] <= times[:-1]
    if np.count_nonzero(negatives) > 0:
        indices = np.nonzero(negatives)
        raise RuntimeError('Found %d non-monotonic timestamps: photon indices %s' %
                           (np.count_nonzero(negatives), indices))

def verify_continuity(times, gap_factor=1000):
    """ Search for improbably long gaps in the photon stream """
    tau = (times[-1] - times[0]) / len(times)
    gaps = (times[1:] - times[:-1]) > gap_factor*tau
    if np.count_nonzero(gaps) > 0:
        print('Found %d large gaps:' % np.count_nonzero(gaps))
        gap_starts = np.nonzero(gaps)
        for s in gap_starts:
            print('    starting at %10d, ending at %10d, lasting %10d' %
                  (times[s], times[s+1], times[s+1] - times[s]))

class TimestampFile(object):
    """ A portable interface for reading photon timestamp files """
    def __init__(self, fname, channel, ftype=None):
        """ TimestampFile(filename, channel)

            Channel number is zero-based """
        self.jiffy = None
        self.metadata = None

        if ftype is None:
            ftype = determine_filetype(fname)

        if ftype == 'pt2':
            self.jiffy = 4e-12
            self.data = pt2_parse.read_pt2(fname, channel)

        elif ftype == 'timetag':
            if channel not in range(4):
                raise RuntimeError(".timetag files only have channels 0..3")
            self.metadata = metadata.get_metadata(fname)
            if self.metadata is not None:
                self.jiffy = 1. / self.metadata['clockrate']
            open(fname) # Ensure a reasonable error is generated
            self.data = timetag_parse.get_strobe_events(fname, 1<<channel)['t']

        elif ftype == 'raw':
            if channel != 0:
                raise RuntimeError("Raw timetag files have only channel 0")
            self.data = np.fromfile(fname, dtype='u8')

        elif ftype == 'timech':
            if channel > 255:
                raise RuntimeError("Raw timetag files have only 256 channels")
            d = np.fromfile(fname, dtype='u8,u1', names='time,chan')
	    self.data = d[d['chan'] == channel]['time']

        else:
            raise RuntimeError("Unknown file type")

        verify_monotonic(self.data)
        verify_continuity(self.data)
            
