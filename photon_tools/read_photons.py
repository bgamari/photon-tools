import os
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

class InvalidChannel(RuntimeError):
    def __init__(self, requested_channel, valid_channels=[]):
        self.requested_channel = requested_channel
        self.valid_channels = valid_channels

    def __str__(self):
        return "Channel %s was requested but this file type only supports channels %s." \
            % (self.requested_channel, self.valid_channels)

class TimestampReader(object):
    """ An abstract reader of timestamp data """
    extensions = []
    def __init__(self):
        self.jiffy = None

class PicoquantFile(TimestampReader):
    """ Read Picoquant PT2 and PT3 timestamp files """
    extensions = ['pt2', 'pt3']
    def __init__(self, fname, channel):
        TimestampReader.__init__(self)
        self.jiffy = 4e-12 # FIXME
        self.data = pt2_parse.read_pt2(fname, channel)

class TimetagFile(TimestampReader):
    """ Read Goldner FPGA timetagger files """
    extensions = ['timetag']
    def __init__(self, fname, channel):
        TimestampReader.__init__(self)
        channels = range(4)
        if channel not in channels:
            raise InvalidChannel(channel, channels)
        self.metadata = metadata.get_metadata(fname)
        if self.metadata is not None:
            self.jiffy = 1. / self.metadata['clockrate']
        if not os.path.isfile(fname):
            raise IOError("File %s does not exist" % fname)
        self.data = timetag_parse.get_strobe_events(fname, 1<<channel)['t']

class RawFile(TimestampReader):
    """ Read raw unsigned 64-bit timestamps """
    extensions = ['times']
    def __init__(self, fname, channel):
        TimestampReader.__init__(self)
        if channel != 0:
            raise InvalidChannel(channel, [0])
        self.data = np.fromfile(fname, dtype='u8')

class RawChFile(TimestampReader):
    """ Read raw unsigned 64-bit timestamps, followed by 8-bit channel number """
    extensions = 'timech'
    def __init__(self, fname, channel):
        TimestampReader.__init__(self)
        if channel > 255:
            raise InvalidChannel(channel, range(256))
        d = np.fromfile(fname, dtype='u8,u1', names='time,chan')
        self.data = d[d['chan'] == channel]['time']

readers = [
    PicoquantFile,
    TimetagFile,
    RawFile,
    RawChFile,
]

def supported_extensions():
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
    exts = supported_extensions()
    root,ext = os.path.splitext(fname)
    return exts.get(ext[1:])

def open(fname, channel, reader=None):
    """
    open(filename, channel)

    Read a timestamp file. Channel number is zero-based.
    """
    if reader is None:
        reader = find_reader(fname)
    if reader is None:
        raise RuntimeError("Unknown file type")

    f = reader(fname, channel)
    verify_monotonic(f.data)
    verify_continuity(f.data)
    return f
