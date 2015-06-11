import os
import warnings
import numpy as np
from photon_tools.io import timetag_parse, pt2_parse, metadata

time_ch_dtype = np.dtype([('time', 'u8'), ('chan', 'u1')])

def verify_monotonic(times, filename):
    """ Verify that timestamps are monotonically increasing """
    if len(times) == 0: return
    negatives = times[1:] <= times[:-1]
    if np.count_nonzero(negatives) > 0:
        indices = np.nonzero(negatives)
        warnings.warn('%s: Found %d non-monotonic timestamps: photon indices %s' %
                      (filename, np.count_nonzero(negatives), indices))

def verify_continuity(times, filename, gap_factor=1000):
    """ Search for improbably long gaps in the photon stream """
    if len(times) == 0: return
    tau = (times[-1] - times[0]) / len(times)
    gaps = (times[1:] - times[:-1]) > gap_factor*tau
    if np.count_nonzero(gaps) > 0:
        msg = '%s: Found %d large gaps:\n' % (filename, np.count_nonzero(gaps))
        gap_starts, = np.nonzero(gaps)
        for s in gap_starts:
            msg += '    starting at %10d, ending at %10d, lasting %10d' % \
                   (times[s], times[s+1], times[s+1] - times[s])
        warnings.warn(msg)

class InvalidChannel(ValueError):
    def __init__(self, requested_channel, valid_channels=[]):
        self.requested_channel = requested_channel
        self.valid_channels = valid_channels

    def __str__(self):
        return "Channel %s was requested but this file type only supports channels %s." \
            % (self.requested_channel, self.valid_channels)

class TimestampFile(object):
    """
    Represents a timestamp file.

    A timestamp file is a file containing a sequential set of integer
    photon arrival timestamps taken in one or more channels.
    """
    def __init__(self, filename, jiffy, valid_channels=None):
        if valid_channels is None:
            valid_channels = self.__class__.valid_channels
        self._valid_channels = valid_channels
        self._fname = filename
        self._jiffy = jiffy

    @classmethod
    def extensions(self):
        """
        A list of supported file extensions
        """
        return []

    @property
    def jiffy(self):
        """
        The timestamp resolution in seconds or ``None`` is unknown.

        :returns: float
        """
        return self._jiffy

    @property
    def valid_channels(self):
        """
        The names of the channels of the file.
        Note that not all of these will have timestamps.

        :returns: list
        """
        return self._valid_channels

    @property
    def metadata(self):
        """
        Metadata describing the data set.

        :returns: dictionary mapping string metadata names to values
        """
        return self._metadata

    @property
    def name(self):
        """ File name of timestamp file """
        return self._fname

    def timestamps(self):
        """
        Read the timestamp data for all channels of the file

        :returns: An array of dtype :var:`time_ch_dtype` containing monotonically
        increasing timestamps annotated with channel numbers.
        """
        data = self._read_all()
        verify_monotonic(data['time'], self._fname)
        verify_continuity(data['time'], self._fname)
        return data

    def channel(self, channel):
        """
        Read photon data for a channel of the file

        :type channel: A valid channel name from :func:`valid_channels`.
        :returns: An array of ``u8`` timestamps.
        """
        self._validate_channel(channel)
        data = self._read_channel(channel)
        verify_monotonic(data, self._fname)
        verify_continuity(data, self._fname)
        return data

    def _read_all(self):
        """ Read the timestamps for all channels """
        raise NotImplementedError()

    def _read_channel(self, channel):
        """ Actually read the data of a channel """
        raise NotImplementedError()

    def _validate_channel(self, channel):
        """ A utility for implementations """
        if channel not in self.valid_channels:
            raise InvalidChannel(channel, self.valid_channels)

class PicoquantFile(TimestampFile):
    """ A Picoquant PT2 and PT3 timestamp file """
    valid_channels = [0,1,2,3]
    extensions = ['pt2', 'pt3']
    def __init__(self, fname):
        # TODO: Read metadata
        TimestampFile.__init__(self, fname, jiffy=4e-12)

    def _read_all(self):
        raise NotImplementedError()

    def _read_channel(self, channel):
        return pt2_parse.read_pt2(self._fname, channel)

class TimetagFile(TimestampFile):
    """ A timestamp file from the Goldner lab FPGA timetagger """
    extensions = ['timetag']
    valid_channels = [0,1,2,3]
    def __init__(self, fname):
        TimestampFile.__init__(self, fname, jiffy = None)
        self._metadata = metadata.get_metadata(fname)
        if self.metadata is not None:
            self._jiffy = 1. / self.metadata['clockrate']
        if not os.path.isfile(fname):
            raise IOError("File %s does not exist" % fname)

    def _read_all(self):
        res = timetag_parse.get_strobe_events(self._fname, 0xf)
        res.dtype.names = time_ch_dtype.names
        return res

    def _read_channel(self, channel):
        return timetag_parse.get_strobe_events(self._fname, 1<<channel)['t']

class RawFile(TimestampFile):
    """ Raw unsigned 64-bit timestamps """
    extensions = ['times']
    valid_channels = [0]
    def __init__(self, fname):
        TimestampFile.__init__(self, fname, jiffy = None)

    def _read_all(self):
        timestamps = np.fromfile(fname, dtype='u8')
        return np.from_records([timestamps, np.zeros_like(timestamps, dtype='u1')],
                               dtype=time_ch_dtype)

    def _read_channel(self, channel):
        return np.fromfile(fname, dtype='u8')

class RawChFile(TimestampFile):
    """ Raw unsigned 64-bit timestamps, followed by 8-bit channel number """
    extensions = ['timech']
    valid_channels = range(256)
    def __init__(self, fname, channel):
        TimestampFile.__init__(self, fname, jiffy = None)

    def _read_all(self):
        return np.fromfile(fname, dtype=time_ch_dtype)

    def _read_channel(self, channel):
        d = self._read_all()
        return d[d['chan'] == channel]['time']

readers = [
    PicoquantFile,
    TimetagFile,
    RawFile,
    RawChFile,
]

def supported_extensions():
    """
    A dictionary mapping supported file extensions to their associated
    :class:`TimestampFile` class.
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

def open(fname, reader=None):
    """
    Read a timestamp file.

    :type fname: ``str``
    :param fname: Path of timestamp file
    :type reader: class inheriting ``TimestampFile``, optional
    :param reader: Specify the file format explicitly
    """
    if reader is None:
        reader = find_reader(fname)
    if reader is None:
        raise RuntimeError("Unknown file type")

    return reader(fname)
