from codecs import (BOM_UTF8,
                    BOM_UTF16_BE, BOM_UTF16_LE,
                    BOM_UTF32_BE, BOM_UTF32_LE)
from datetime import datetime as dt
from glob import glob
import itertools
import logging
import matplotlib as mpl
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
import numpy as np
import os
from obspy import read, Stream
from obspy.signal.filter import envelope
import pandas as pd
import pygimli as pg
from pygimli.physics import Refraction, TravelTimeManager
from scipy import fftpack
from scipy.signal import boxcar, find_peaks
import sys
from typing import NamedTuple
import warnings

class nt_procmodes(NamedTuple):
    inactive: int = 0
    pick: int = 1
    trace_mute: int = 2
    trace_reverse: int = 3
    vel_est: int = 4
    
class nt_plotmodes(NamedTuple):
    seismograms: int = 0
    var_dens: int = 1
    
    def min(self):
        return self.seismograms
        
    def max(self):
        return self.var_dens

class nt_scrollmodes(NamedTuple):
    zoom: int = 0
    amplitude: int = 1
    
COMPUTE_KEYWORDS = ["cos", "autopick", "geometry"]
FILTER_KEYWORDS = ["lp", "hp", "bp", "bs", "3pf",
                   "hold", "remove"]
PLOT_KEYWORDS = ["traveltimes", "pseudosection", "pickperc", "spectrum"]
PICKSET_KEYWORDS = ["create", "load", "unload", "rename", "use", "delete",
                    "copy", "rename", 
                    "export", "import"]
SELECT_KEYWORDS = ["sin", "rin", "aoffset", "midpoint"]
CREATE_KEYWORDS = ["waveforms", "traveltimes", "mesh", "wavelet", "velmod"]
LOAD_KEYWORDS = ["config", "mesh", "scheme"]

PROC_MODES = nt_procmodes()
PLOT_MODES = nt_plotmodes()
SCROLL_MODES = nt_scrollmodes()

CONSOLECOLORS = {
    'DEBUG': 34,
    'INFO': 32,
    'WARNING': 93,
    'ERROR': 31,
    'CRITICAL': 31
}

BOMS = [BOM_UTF8,
        BOM_UTF16_BE, BOM_UTF16_LE,
        BOM_UTF32_BE, BOM_UTF32_LE]

def check_encoding(path):
    """Check the encoding of the file."""
    
    try:
        with open(path) as f:
            f.readline()
            return True
    except:
        return False

def check_bom(path):
    """Check for byte order mark in file."""
    
    hasbom = False
    
    with open(path) as f:
        # read first line of file and convert to bytes
        l = bytes(f.readline(), f.encoding)
        
        for bom in BOMS:
            if l.startswith(bom):
                hasbom = True
    
    return hasbom

def print_progressbar(current, total,
                     prefix='Progress', suffix='completed', length=30):
    """Print a progress bar in the terminal.
    
    Parameters
    ----------
    current: int
        Current state, e.g., number of currently processed items.
        
    total: int
        Total number, e.g., total number of items to process.
        
    prefix: str
        Leading text, i.e., text to be put in front of the progress bar.
        
    suffix: str
        Trailing text, i.e., text to be put behind the progress bar.
        
    length: int
        Total length of the progress bar in number of characters.
    """
    filled = int(length * current // total)
    bar = '\033[32m' + '=' * filled + \
          '\033[31m' + '-' * (length - filled) + '\033[0;0m'
    print('\r%s <%s> %4.1f%% %s' % 
        (prefix, bar, (current / total) * 10**2, suffix), end='\r')

def slidingwindow(tr, gate_length):
    """Move a window along the trace data.
    
    Parameters
    ----------
    tr: obspy Trace
        The Trace object containing the data.
        
    gate_length: float
        Gate length of the sliding window (%/100).
        
    Returns
    -------
        iters: list of iterators
            List of sliding windows.
    """
    
    win_len = int((tr.stats.npts - 1) * gate_length)

    iters = itertools.tee(np.arange(tr.stats.npts), win_len)
    
    for it, skipped in zip(iters, itertools.count()):
        for _ in np.arange(skipped):
            next(it, None)
    
    return zip(*iters)

def threepointfilter(data, n=1):
    """Recursively apply a three-point filter to the data.
    
    Parameters
    ----------
    data: numpy array
        Unfiltered data.
        
    n: int, default 1
        How many times to apply the filter
        
    Returns
    -------
    datafilt: numpy array
        Filtered data.
    """
    
    datafilt = data.copy()
    for i in np.arange(data.shape[0] - 2):
        datafilt[i] = .5 * data[i] + .25 * (data[i-1] + data[i+1])
    
    if n == 1:
        return datafilt
    else:
        return threepointfilter(datafilt, n-1)

def compute_perpendicular_vector(v):
    """Compute the perpendicular vector.
    
    Parameters
    ----------
    v: numpy array
        Given vector.
        
    Returns
    -------
    p: numpy array
        The perpendicular vector.
    """
    
    p = np.array([-v[1], v[0]])
    return p

def compute_line_intersection(sa, ea, sb, eb):
    """Compute the intersection of two lines by using vectors.
    
    Parameters
    ----------
    sa: numpy array
        Start point of line a.
    
    ea: numpy array
        End point of line a.
        
    sb: numpy array
        Start point of line b.
        
    eb: numpy array
        End point of line b.
        
    Returns
    -------
    isec: numpy array
        Intersection point or None if there is no intersection between
        the given lines.
    """
    
    # compute line vectors
    da = ea - sa
    db = eb - sb
    dp = sa - sb
    dap = compute_perpendicular_vector(da)
    denom = np.dot(dap, db)
    num = np.dot(dap, dp)
    
    # check if denominator is non-zero
    if denom != 0:
        isec = (num / denom) * db + sb
    else:
        isec = None
    
    return isec

def _loginput(self, message, *args, **kwargs):
    """Custom logging of console input.
    
    Parameters
    ----------
    message: string
        The logging message
    """
    
    if self.isEnabledFor(15):
        self._log(15, message, args, **kwargs)

def _logauto(self, message, *args, **kwargs):
    """Logging of automatically executed commands.
    
    Parameters
    ----------
    message: string
        The logging message
    """
    
    if self.isEnabledFor(16):
        self._log(16, message, args, **kwargs)
        
def _logautoinfo(self, message, *args, **kwargs):
    """Logging of output related to automatically executed commands.
    
    Parameters
    ----------
    message: string
        The logging message
    """
    
    if self.isEnabledFor(17):
        self._log(17, message, args, **kwargs)
    
def create_logger(name, logfile):
    """Create a logger.
    
    Parameters
    ----------
    name: string
        Name of the logger.
        
    logfile: string
        Path to the logfile inlcuding the file name.
        
    Returns
    -------
    logger: logger object
        The configured logger.
    
    """
    
    if name in logging.root.manager.loggerDict.keys():
        return logging.getLogger(name)
    else:
        logging.addLevelName(15, 'INPUT')
        logging.addLevelName(16, 'AUTO')
        logging.addLevelName(17, 'INFO')
        logging.Logger.input = _loginput
        logging.Logger.auto = _logauto
        logging.Logger.autoinfo = _logautoinfo
        
        # create logger
        logger = logging.getLogger(name)
        logger.setLevel(logging.DEBUG)
        logger.propagate = False
        
        # create file handler 
        fh = logging.FileHandler(logfile)
        fh.setLevel(logging.DEBUG)
        formatter = logging.Formatter('%(asctime)s %(levelname)-8s %(message)s',
                                      datefmt='%Y/%m/%d %H:%M:%S')
        fh.setFormatter(formatter)
        logger.addHandler(fh)
        
        # create console handler
        ch = logging.StreamHandler(sys.stdout)
        ch.setLevel(logging.INFO)
        formatter = logging.Formatter(
            '\033[1m%(levelname)-8s\033[0m: %(message)s')
        ch.setFormatter(formatter)
        logger.addHandler(ch)
        
        return logger

def get_filelist(directory, method='seismic'):
    """Check for valid files in the directory and return the 
    corresponding search string.
    
    Parameters
    ----------
    directory: str
        Path to the directory to check.
    
    method: str, default 'seismic'
        Define for which (geophysical) data to look for. 
    
    Returns
    -------
    search: str
        Regular expression to list all valid files found in the data
        directory.
    """
    
    if method == 'seismic':
        search = os.path.join(directory, '*' + '.sg2')
        if len(glob(search)) > 0:
            # ~ self._logger.info('Raw data with supported file ' \
                # ~ 'extension found (dmt summit)')
            return search, 'SEG2', False
        
        search = os.path.join(directory, '*' + '.dat')
        if len(glob(search)) > 0:
            # ~ self._logger.info('Raw data with supported file ' \
                # ~ 'extension found (geometric geode)')
            return search, 'SEG2', False
                
        search = os.path.join(directory, '*' + '.syn')
        if len(glob(search)) > 0:
            # ~ self._ftype = 'MSEED'
            # ~ self._syndata = True
            # ~ self._logger.info('Raw data with supported file ' \
                # ~ 'extension found (synthetic)')
            return search, 'MSEED', True
    
    return '', '', False
    # ~ self._logger.warning('No supported data found')
