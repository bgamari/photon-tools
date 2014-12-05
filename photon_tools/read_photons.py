import numpy as np
from photon_tools import timetag_parse, pt2_parse, metadata

def verify_monotonic(times):
    """ Verify that timestamps are monotonically increasing """
    if len(times) == 0: return
    negatives = times[1:] <= times[:-1]
    if np.count_nonzero(negatives) > 0:
        indices = np.nonzero(negatives)
        raise RuntimeError('Found %d non-monotonic timestamps: photon indices %s' %
                           (np.count_nonzero(negatives), indices))

def verify_continuity(times, gap_factor=1000):
    """ Search for improbably long gaps in the photon stream """
    if len(times) == 0: return
    tau = (times[-1] - times[0]) / len(times)
    gaps = (times[1:] - times[:-1]) > gap_factor*tau
    if np.count_nonzero(gaps) > 0:
        print('Found %d large gaps:' % np.count_nonzero(gaps))
        gap_starts, = np.nonzero(gaps)
        for s in gap_starts:
            print('    starting at %10d, ending at %10d, lasting %10d' %
                  (times[s], times[s+1], times[s+1] - times[s]))

class TimestampReader(object):
    extensions = []
    def __init__(self, fname, channel):
        raise Unimplemented

class Pt2File(TimestampReader):
    """ Read Picoquant PT2 timestamp files """
    extensions = ['pt2']
    def __init__(self, fname, channel):
        self.jiffy = 4e-12 # FIXME
        self.data = pt2_parse.read_pt2(fname, channel)

class TimetagFile(TimestampReader):
    """ Read Goldner FPGA timetagger files """
    extensions = ['timetag']
    def __init__(self, fname, channel):
        if channel not in range(4):
            raise RuntimeError(".timetag files only have channels 0..3")
        self.metadata = metadata.get_metadata(fname)
        if self.metadata is not None:
            self.jiffy = 1. / self.metadata['clockrate']
        open(fname) # Ensure a reasonable error is generated
        self.data = timetag_parse.get_strobe_events(fname, 1<<channel)['t']

class RawFile(TimestampReader):
    """ Read raw unsigned 64-bit timestamps """
    extensions = ['raw']
    def __init__(self, fname, channel):
        if channel != 0:
            raise RuntimeError("Raw timetag files have only channel 0")
        self.data = np.fromfile(fname, dtype='u8')

class RawChFile(TimestampReader):
    """ Read raw unsigned 64-bit timestamps, followed by 8-bit channel number """
    extensions = 'timech'
    def __init__(self, fname, channel):
        if channel > 255:
            raise RuntimeError("Raw timetag files have only 256 channels")
        d = np.fromfile(fname, dtype='u8,u1', names='time,chan')
        self.data = d[d['chan'] == channel]['time']

readers = [
    Pt2Reader,
    TimetagReader,
    RawReader,
    RawTimeChReader,
]

def reader_extensions():
    """
    Construct a map from supported file extensions to their
    associated reader.
    """
    extensions = {}
    for reader in readers:
        for ext in reader.extensions:
            extensions[ext] = reader
    return extensions

def find_reader(fname):
    exts = reader_extensions()
    root,ext = os.path.splitext(fname)
    return exts.get(ext)

class TimestampFile(object):
    """
    A format-agnostic interface for reading photon timestamp files.
    """

    @staticmethod
    def open(fname, channel, reader=None):
        """
        TimestampFile(filename, channel)

        Channel number is zero-based.
        """
        self.jiffy = None
        self.metadata = None

        if reader is None:
            reader = determine_filetype(fname)
        if reader is None:
            raise RuntimeError("Unknown file type")

        reader.read(fname, channel)
        verify_monotonic(self.data)
        verify_continuity(self.data)
