from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from . import MethodManager
from formikoj.utils import *

class SeismicRefractionManager(MethodManager):
    
    'Class for managing the processing of seismic refraction data'
    
    def __init__(self, workdir):
        """Create an instance of the SeismicRefraction class.
        
        Parameters
        ----------
        workdir: str
            Path of the project directory.
        """
        
        if not MethodManager.__init__(self, workdir):
            return
        self._init_refraction_attributes()
        if 'left' in mpl.rcParams['keymap.back']:
            mpl.rcParams['keymap.back'].remove('left')
        if 'right' in mpl.rcParams['keymap.forward']:
            mpl.rcParams['keymap.forward'].remove('right')
        # ~ self._logger.info('Created instance of SeismicRefractionManager')
        
        # check whether the working directory contains shot files and start in
        # data preview mode
        filelist, self._ftype, self._syndata = get_filelist(self._workdir)
        if filelist != '':
            # ~ self._logger = create_logger('useis.SeismicRefractionManager',
            self._logger = create_logger(self._workdir,
                             os.path.join(self._workdir,
                                          self._logfile))
            self._pvd = 1
            self._preview_data()
        else:
            if not self._check_directory_structure():
                print('\033[31mDIRECTORY STRUCTURE INVALID\033[0;0m')
                return
            
            # ~ self._logger = create_logger('useis.SeismicRefractionManager',
            self._logger = create_logger(self._workdir,
                                         os.path.join(self._procdir,
                                                      self._logfile))
            
            if self._check_database():
                self._load_project()
            else:
                cgf = self._check_geometryfile(os.path.join(self._geomdir,
                                                            self._geomfile))
                if cgf == 1:
                    self.slh = SQLiteHandler(self._dbfile)
                    self._create_database()
                    self._logger.info('Project database created')
                    if self._read_geometry():
                        self._apply_geometry()
                        self._read_data()
                    else:
                        os.remove(self._dbfile)
                        self._logger.info('Project database deleted')
                elif cgf == 0:
                    self._pvd = 2
                    self._preview_data()
    
    def _check_database(self):
        """Check if the database is ready for processing."""
        
        if os.path.exists(self._dbfile):
            self.slh = SQLiteHandler(self._dbfile)
            
            cmd = """SELECT name
                     FROM sqlite_master
                     WHERE name == \'applied_geometry\'"""
            if self.slh.read_data(cmd).empty:
                ret = 0
            else:
                ret = 1
            self.slh.close_connection()
        else:
            ret = 0
        
        return ret

    def _create_database(self):
        """Create the project database."""
        
        self.slh.create_table(
            'geometry',
            ['station_id',
             'x', 'y', 'z',
             'geophone', 'shot',
             'first_geophone', 'num_geophones'],
            ['INTEGER',
             'FLOAT', 'FLOAT', 'FLOAT',
             'BOOL', 'INTEGER',
             'INTEGER', 'INTEGER'],
            ['station_id']
        )
        
        self.slh.create_table(
            'shots',
            ['shot_index_number', 'station_id'],
            ['INTEGER', 'INTEGER'],
            ['shot_index_number'],
            [('station_id', 'geometry', 'station_id')]
        )
        
        self.slh.create_table(
            'receivers',
            ['receiver_index_number', 'station_id', 'polarity'],
            ['INTEGER', 'INTEGER', 'INTEGER'],
            ['receiver_index_number'],
            [('station_id', 'geometry', 'station_id')]
        )
    
        self.slh.create_table(
            'applied_geometry',
            ['shot_index_number', 'receiver_index_number',
             'absolute_offset', 'midpoint_x', 'midpoint_y'],
            ['INTEGER', 'INTEGER',
             'FLOAT', 'FLOAT', 'FLOAT'],
            ['shot_index_number', 'receiver_index_number'],
            [('shot_index_number', 'shots', 'shot_index_number'), 
             ('receiver_index_number', 'receivers', 'receiver_index_number')]
        )
            
        self.slh.create_table(
            'fbpicks',
            ['pickset', 'shot_index_number', 'receiver_index_number', 
             'traveltime'],
            ['TEXT', 'INTEGER', 'INTEGER', 'FLOAT'],
            ['pickset', 'shot_index_number', 'receiver_index_number'],
            [('shot_index_number', 'shots', 'shot_index_number'),
             ('receiver_index_number', 'receivers', 'receiver_index_number')]
        )
    
    def _init_refraction_attributes(self):
        """Create and initialize attributes relevant for processing 
        seismic refraction data.
        """
        
        # data handling
        self._st = None
        self._pst = None
        self._cos = False
        self._cosst = None
        self._selected = ''
        self._actssns = np.array([])
        self._actrsns = np.array([])
        
        # data processing
        self._procmode = PROC_MODES.inactive
        self._scrollmode = SCROLL_MODES.zoom
        self._plotmode = PLOT_MODES.seismograms
        self._filtered = ''
        self._filterhold = False
        self._scaling = 1
        self._ylim = (-9999,-9999)
        self._autoplot = False
        self._npts = 0
        
        # plotting related
        self._title = ''
        self._mode = ''
        self._xlabel = 'SIN'
        self._fontsize = 10
        self._marker = 'x'
        
        self._fig = None
        self._ax = None
        self._ttfig = None
        self._ttax = None
        self._ppfig = None
        self._ppax = None
        self._vafig = None
        self._vaax = None
        self._plotset = []
        
        self._psg_sg = None
        self._psg_fb = []
        # ~ self._psg_fb = None
        self._psg_isvd = None

        self._pressed = False
        self._moved = False

        self._picklineplot = None
        self._vellineplot = None
        self._vellines = None
        self._velline = None
        self._vellinetext = []
        
        # set names auxiliary files
        self._geomfile = 'geometry.csv'
        self._logfile = 'history.log'

        # data handling
        self._syndata = False
        self._ftype = 'SEG2'
        self._fext = '.sg2'
        self._data = None
        self._pvd = 0
        self._st = None
        
        # data processing
        self._geomload = False
        self._geomapp = False
        self._3d = False
        
        # picking related
        self._picksets = []
        self._activeps = None
        self._picks = {}        
        self._numpicks = -9999.
        self._pickscatter = None

    def _check_plotting_ready(self):
        """Check if the project is ready for data plotting."""
        
        if self._data == None:
            self._logger.error('Data not loaded')
            return 0
        
        if self._st == None or len(self._st) == 0:
            self._logger.error('No data selected')
            return 0
            
        return 1

    def _check_processing_ready(self):
        """Check if the project is ready for data processing."""
        
        if self._data == None:
            self._logger.error('Data not loaded')
            return 0
        
        if not self._geomapp and not self._pvd:
            self._logger.error('Geometry not applied')
            return 0
        
        return 1
    
    def _load_project(self):
        """Load project information from database."""
        
        self.slh = SQLiteHandler(self._dbfile)
        
        geom = self.slh.read_data("""SELECT *
                                     FROM geometry""")
        if len(geom) > 0: 
            self._geomload = True
            self._3d = not np.allclose(geom.y, 0.)
        
        appgeom = self.slh.read_data("""SELECT COUNT(*)
                                        FROM applied_geometry""")
        self._numpicks = int(appgeom.iloc[0])                                
        if self._numpicks > 0: self._geomapp = True
        
        self._logger.info('Project information loaded')
        
        if self._geomapp:
            self._read_data()
            self._activate_pickset(['picks'])
            # ~ self._load_picks()

    def _copy_header_info(self, f, idx):
        """Copy relevant geometry information from the project file to the
        header of the file objects.
        
        Parameters
        ----------
        filename: str
            String containing the path to the file.
        
        idx: int
            Index of the file in the data dictionary of the MethodManager
            object.
        """
        
        shts = self.slh.read_data("""SELECT s.shot_index_number sin,
                                            g.first_geophone fg,
                                            g.shot ssn
                                     FROM shots s, geometry g
                                     ON s.station_id == g.station_id""")
                                                        
        for j, tr in enumerate(self._data[idx]):
            if self._syndata:
                ssn = int(tr.stats.station)
                rsn = int(tr.stats.channel)
            else:
                ssn = int(tr.stats.seg2['SHOT_SEQUENCE_NUMBER'])
                rsn = int(tr.stats.seg2['CHANNEL_NUMBER'])
            
            # set channel number
            tr.stats.channel = str(rsn)
            
            # set shot index number
            tr.stats.SIN = int(shts[shts.ssn==ssn].sin)
            
            # determine channel offset
            fg = int(shts[shts.ssn==ssn].fg)
            
            # set receiver index number
            tr.stats.RIN = j + fg
            
            # set file name
            tr.stats.file = f

    def _read_data(self):
        """Read data from disk."""
        
        self._check_directory_structure()
        
        self._data = {}
        self._fnames = {}

        if self._geomload:
            fnrs = self.slh.read_data("""SELECT shot
                                         FROM geometry
                                         WHERE shot > -1""")

        # ~ search = self._get_filelist()
        search, self._ftype, self._syndata = get_filelist(
            self._workdir if self._pvd == 1 else self._datadir)
        numfiles = float(len(glob(search)))
        for i, f in enumerate(sorted(glob(search))):
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                st = read(f, format=self._ftype)
            
            if not self._pvd and self._geomload:
                if self._syndata:
                    ssn = float(st[0].stats.station)
                else:
                    ssn = float(st[0].stats.seg2['SHOT_SEQUENCE_NUMBER'])
                
                # discard files without geometry information
                if ssn in fnrs.values:
                    self._data[i] = st.copy()
                    if self._npts == 0:
                        self._npts = st[0].stats.npts
                    
                    # copy information from seg2 header to stream header                        
                    self._copy_header_info(f, i)
            else:
                self._fnames[i] = f.split('/')[-1]
                self._data[i] = st.copy()
                
            # progress bar
            print_progressbar(i + 1, numfiles)
        print('')
        
        mk = min(list(self._data.keys()))
        self._ylim = (self._data[mk][0].stats.npts / \
            self._data[mk][0].stats.sampling_rate, 0)
        
        self._logger.info('Read %d files' % (len(self._data)))

    def _apply_geometry(self):
        """Apply geometry information and save it in project file."""
        
        if not self._geomload:
            self._logger.critical('Geometry information incomplete')
        
        # get shot geometry information from database
        cmd = """SELECT s.shot_index_number, gs.x sx, gs.y sy, gs.z sz, 
                        gs.first_geophone fg, gs.num_geophones ng
                 FROM geometry gs
                 JOIN shots s ON gs.station_id==s.station_id"""
        sgeom = self.slh.read_data(cmd)
        
        tmp = np.empty((0, 8))
        for i, sht in sgeom.iterrows():
            
            # get receiver geometry from database
            cond = ', '.join([str(e) for e in np.arange(sht.fg, 
                                                        sht.fg+sht.ng)])
            cmd = """SELECT r.receiver_index_number, gr.x rx, gr.y ry, gr.z rz
                     FROM geometry gr
                     JOIN receivers r ON gr.station_id==r.station_id
                     WHERE r.receiver_index_number in (%s)""" % (cond)
            rgeom = self.slh.read_data(cmd).to_numpy()
            
            tmp = np.vstack((tmp,
                             np.hstack((np.tile(sht[:-2], (rgeom.shape[0], 1)), 
                                        rgeom))))
        
        geom = pd.DataFrame(tmp, columns=['shot_index_number',
                                          'sx', 'sy', 'sz',
                                          'receiver_index_number',
                                          'rx', 'ry', 'rz'])
                                          
        # compute absolute offsets
        geom['absolute_offset'] = (((geom.sx-geom.rx)**2 +
                                    (geom.sy-geom.ry)**2 +
                                    (geom.sz-geom.rz)**2)**(1/2)).round(0)
        
        # compute midpoints
        geom['midpoint_x'] = (geom.sx + geom.rx)/2
        geom['midpoint_y'] = (geom.sy + geom.ry)/2
        
        # write to database
        self.slh.write_data('applied_geometry', 
                            geom[['shot_index_number', 'receiver_index_number', 
                                  'absolute_offset', 
                                  'midpoint_x', 'midpoint_y']])
        self._geomapp = True
        self._numpicks = geom.shape[0]
        self._logger.info('Applied geometry')
        
        # create default pickset
        fbpicks = geom.loc[:, ['shot_index_number', 'receiver_index_number']]
        fbpicks['traveltime'] = -1
        fbpicks['pickset'] = 'picks'
        self.slh.write_data('fbpicks', fbpicks)
        self._logger.info('Standard pickset \'picks\' created')
        self._activate_pickset(['picks'])

    def _check_geometryfile(self, path):
        """Check if the geometry exists and if the number of columns is
        sufficient.
        """
        
        if not os.path.exists(path):
            self._logger.warning('Geometry file not found')
            return 0
            
        if not check_encoding(path):
            self._logger.critical('Geometry file is not utf-8 encoded')
            return -1
            
        if check_bom(path):
            self._logger.critical('Geometry file contains byte order mark')
            return -1
            
        raw = pd.read_csv(path, sep=None, engine='python')
            
        if raw.shape[1] != 7:
            self._logger.critical(
                'Insufficient geometry information (number of columns)')
            return -1
            
        return 1

    def _read_geometry(self):
        """Read the geometry file."""
        
        # read geometry information from file
        path = os.path.join(self._geomdir, self._geomfile)
        geom = pd.read_csv(path,
                           sep=None, engine='python',
                           names=['x', 'y', 'z', 
                                  'geophone', 'shot', 
                                  'first_geophone', 'num_geophones'])
        geom.insert(0, 'station_id', np.arange(len(geom)) + 1)
        
        # write geometry information to database
        self.slh.write_data('geometry', geom)
        self._geomload = True
        self._3d = not np.allclose(geom.y, 0.)
        self._logger.info('Read geometry information from file')
        
        # extract receivers
        # ~ recs = pd.DataFrame(
            # ~ geom[geom.geophone > 0].sort_values(by=['x'])['station_id'])
        recs = pd.DataFrame(geom[geom.geophone > 0]['station_id'])
        recs.insert(0, 'receiver_index_number', np.arange(len(recs)) + 1)
        recs['polarity'] = 1
        self.slh.write_data('receivers', recs)
        self._logger.info('Extracted receiver geometry')

        # extract shots
        shots = pd.DataFrame(
            geom[geom.shot > 0].sort_values(
                # ~ by=['first_geophone', 'x'])['station_id'])
                by=['first_geophone', 'station_id'])['station_id'])
        
        search, _, _ = get_filelist(
            self._workdir if self._pvd else self._datadir)
        
        if shots.shape[0] > float(len(glob(search))):
            self._logger.error('Number of shots in geometry file exceeds ' + \
                               'the number of shot files')
            return 0
        else:
            shots.insert(0, 'shot_index_number', np.arange(len(shots)) + 1)
            self.slh.write_data('shots', shots)
            self._logger.info('Extracted shot geometry')
            return 1
        
    def _compute_cos(self):
        """Compute the common offset stack."""
        
        self._check_processing_ready()
        
        if self._cosst != None:
            self._st = self._cosst.copy()
            self._logger.info('Loaded existing common offset stack')
        else:
            # get min data key
            mk = np.fromiter(self._data.keys(), dtype=int).min()
            
            # get receiver geometry
            cmd = """SELECT receiver_index_number rin, polarity
                     FROM receivers r LEFT JOIN geometry g
                     ON r.station_id == g.station_id"""
            recs = self.slh.read_data(cmd)
            
            # get absolute offsets from database
            cmd = """SELECT absolute_offset
                     FROM applied_geometry"""
            ao = self.slh.read_data(cmd).to_numpy()
            
            # compute unique absolute offsets
            self._aoffs = np.unique(ao)

            self._cosst = Stream()
            
            for i, uao in enumerate(self._aoffs):
                # ~ print('aoffset', uao)
                # select traces with given absolute offset
                self.select('aoffset %d' % (uao), auto=True)
            
                tmp = self._st.copy()
                # ~ print(len(tmp))
                if len(tmp) > 1:
                    # stack traces
                    stacked_data = np.zeros_like(tmp[0].data, dtype=float)
                    for tr in tmp:
                        # ~ print('RIN', tr.stats.RIN)
                        state = int(recs[recs.rin == tr.stats.RIN].polarity)
                        # ~ print('state', state)
                        stacked_data += (tr.data * state)
                        # ~ print(stacked_data)
                    stacked_data /= len(tmp)
                else:
                    stacked_data = np.array([0] * self._npts)
                
                stacked_trace = self._data[mk][0]
                stacked_trace.data = stacked_data
                stacked_trace.stats.channel = '%d' % (uao)
                stacked_trace.stats.xtl = '%d' % (uao)
                stacked_trace.stats.RIN = recs[recs.polarity==1].iloc[0].rin
                
                # add stacked trace to constant offset stream
                self._cosst += stacked_trace.copy()
                # ~ print(np.max(stacked_data))
                # progress bar
                print_progressbar(i + 1, len(self._aoffs))
            print('')
            
            self._st = self._cosst.copy()
            self._logger.info('Computed the common offset stack')
        self._pst = self._st.copy()
        self._cos = True
        self._mode = 'common offset stack'.upper()
        self._selected = ''
        self._xlabel = 'AOFFSET'
        self._procmode = PROC_MODES.inactive
    
    def _compute_autopicks(self, st, numfb, offset = 0):
        """Compute the first break picks for the given data.
        
        Parameters
        ----------
        st: obspy Stream
            Stream object holding the Traces for which the first break picks
            are to be computed.
        
        numfb: int or float
            Total number of first break picks to compute.
        
        offset: int or float, default 0
            Offset for total number of first break picks in case multiple
            data sets are considered.
        """
        
        gate_length = .01 # (%/100)
        
        # normalize trace data in stream object
        st.normalize()
        
        # ~ # detrend trace data
        # ~ for i, tr in enumerate(st):
            # ~ st[i].data -= np.median(tr.data)
            
        # determine first break picks
        with self.slh.dbc:
            for i, tr in enumerate(st):
                tr.data -= np.median(tr.data)
                
                # compute envelope of trace data
                env = envelope(tr.data)
                
                # smooth the envelope with a boxcar function
                # ~ bc = boxcar(int((tr.stats.npts) * gate_length))
                # ~ # bc = boxcar(int((tr.stats.npts - 1) * .1))
                # ~ env = np.convolve(env, bc, mode='same')
                
                # compute energy ratio in overlapping windows
                er = []
                for vals in slidingwindow(tr, gate_length):
                    vals = np.array(vals)
                    
                    idx = int(np.floor(vals.shape[0] / 2))
                    er.append(np.sum(env[vals[idx:]]) / \
                              np.sum(env[vals[0:idx]]))
                
                # determine position of the first onset
                ifb = np.argmax(er)/tr.stats.sampling_rate
                cmd = """UPDATE fbpicks
                         SET traveltime = %.6f
                         WHERE pickset == \'%s\' AND
                               shot_index_number == %d AND
                               receiver_index_number == %d""" % \
                               (ifb, self._activeps, 
                                tr.stats.SIN, tr.stats.RIN)
                
                self.slh.dbc.execute(cmd)
                
                print_progressbar(i + 1 + offset, numfb)

        # ~ return fbs

    def _manage_autopicking(self, params):
        """Mange the autopicking process.
        
        Parameters
        ----------
        params: numpy array
            Array holding the parameters deciding which autopicking strategy
            to use.
        """
        
        self._check_processing_ready()
        
        cmd = """SELECT *
                 FROM shots"""
        shts = self.slh.read_data(cmd)
        
        # create pickset if not found in project file
        cmd = """SELECT DISTINCT pickset
                 FROM fbpicks"""
        ps = self.slh.read_data(cmd)
        if params[1] in ps.values:
            self._activate_pickset([params[1]])
        else:
            self._create_pickset([params[1]])
            
        if params[0] == 'cur':
            self._check_plotting_ready()
            self._compute_autopicks(self._st.copy(),
                                    numfb=len(self._st))
            print('')
            
        elif params[0] == 'all':
            for i, (k, self._st) in enumerate(self._data.items()):
                if self._filtered != '':
                    self.filter(self._filtered.lower(), auto=True)
                self._compute_autopicks(self._st,
                                        numfb=self._numpicks,
                                        offset=i*len(self._st))
            print('')
            self._st = None
        
        # ~ self._load_picks()

    def _parse_compute_params(self, options):
        """Check the basic validity of the compute parameters and return them
        as numpy array.
        
        Parameters
        ----------
        options: str, default ''
            Options for different computations.
        
        Returns
        -------
        params: numpy array
            Parameters deciding which plot method to use.
        """
        
        if not isinstance(options, str):
            self._logger.error('Parameter of datatype str expected')
            return 0
            
        if options == "":
            self._logger.critical('No compute command provided')
            return 0
        
        option = options.lower().rstrip()
        params = options.split(" ")    
        
        if params[0] not in COMPUTE_KEYWORDS:
            self._logger.error('Invalid keyword given')
            return 0
        
        if params[0] == "autopick":
            if params[1] not in ["all", "cur"]:
                self._logger.error('Invalid autopick option')
                return 0
            
            if len(params) != 3:
                self._logger.error('Insufficient number of parameters')
        
        return params

    def compute(self, options):
        """Perform various computations based on the available data.
        
        Parameters
        ----------
        options: str
            The string has to contain a valid keyword ('cos', 'autopick',
            'geometry') to call the respective methods.
        """
        
        self._logger.input('compute ' + options)
        
        params = self._parse_compute_params(options)
        if params == 0: return
        
        if params[0] == 'autopick':
            self._manage_autopicking(params[1::])
        elif params[0] == 'geometry':
            self._read_geometry()
            self._apply_geometry()
            if self._data == None: self._read_data()
        elif params[0] == 'cos':
            self._compute_cos()

    def _parse_filter_params(self, options):
        """Check the basic validity of the plot parameters and return them
        as numpy array.
        
        Parameters
        ----------
        options: str, default ''
            Options for plotting.
        
        Returns
        -------
        params: numpy array
            Parameter deciding which plot method to use.
        """
        
        if not isinstance(options, str):
            self._logger.error('Parameter of datatype str expected')
            return 0
            
        if options == "":
            self._logger.error('No filter criteria provided')
            return 0
        
        option = options.lower().rstrip()
        params = options.split(" ")
        
        if params[0] not in FILTER_KEYWORDS:
            self._logger.error('Invalid keyword given')
            return 0
            
        if params[0] in ["lp", "hp", "3pf", "hold"]:
            if len(params) != 2:
                self._logger.error('Insufficient number of filter parameters')
                return 0
            
        if params[0] in ["bp", "bs"]:
            if len(params) != 3:
                self._logger.error('Insufficient number of filter parameters')
                return 0
        
        if params[0] in FILTER_KEYWORDS[0:4]:
            for p in params[1::]:
                try:
                    float(p)
                except: 
                    self._logger.error('Filter parameters need to be numeric')
                    return 0
        
        return params

    def filter(self, options, auto=False):
        """Filter the currently selected trace data.
        
        Parameters
        ----------
        options: str
            Parameters used for filtering the trace data. The string has to
            contain a valid keyword ('lp', 'hp', 'bp', 'bs',
            '3pf', 'hold', 'remove') at the beginning followed by the 
            required parameters.
            TODO: describe parameters for each keyword.
        
        auto: bool
            True if filter method was invoked automatically.
            
        Examples
        --------
        >>> from useis import SeismicRefractionManager
        >>> srm = SeismicRefractionManager(<project_directory>)
        >>> srm.filter('filter lp 120')
        """
        
        if auto:
            self._logger.auto('filter ' + options)
        else:
            self._logger.input('filter ' + options)
        
        if not self._check_processing_ready():
            return
            
        if not self._check_plotting_ready():
            return
        
        params = self._parse_filter_params(options)
        if params == 0: return
        
        if params[0] in ['bp','bs']:
            # determine filter type
            ft = 'bandpass' if params[0] == 'bp' \
                else 'bandstop'
                
            # get min and max frequency
            freqs = np.array([float(params[1]), float(params[2])])
            freqmin = np.min(freqs)
            freqmax = np.max(freqs)
            
            # apply filter
            if self._pst != None:
                self._st = self._pst.copy()
            self._st.filter(ft,
                            freqmin = freqmin,
                            freqmax = freqmax,
                            corners = 2,
                            zerophase = True)
            self._filtered = ' '.join(params).upper()
            logmsg = 'Applied %s filter (%.1f to %.1f Hz)' % \
                (ft, freqmin, freqmax)
            
        elif params[0] in ['hp','lp']:
            # determine filter type
            ft = 'highpass' if params[0] == 'hp' \
                else 'lowpass'
                
            # get min and max frequency
            freq = float(params[1])
            
            # apply filter
            if self._pst != None:
                self._st = self._pst.copy()
            self._st.filter(ft,
                            freq = freq,
                            corners = 2,
                            zerophase = True)
            self._filtered = " ".join(params).upper()
            logmsg = 'Applied %.1f Hz %s filter' % (freq, ft)
            
        elif params[0] == 'hold':
            if self._filtered != "":
                if params[1] == 'on':
                    self._filterhold = True
                elif params[1] == 'off':
                    self._filterhold = False
                else:
                    self._logger.error('Unknown filter option \'%s\'' % 
                        (params[1]))
                    return
                logmsg = 'Set filter hold %s' % (params[1])
                
            else:
                self._logger.error('No filter parameters provided')
                return
            
        elif params[0] == 'remove':
            if self._filtered != '':
                self._st = self._pst.copy()
                self._filtered = ''
                self._filterhold = False
                logmsg = 'Filter settings removed'
            else:
                logmsg = 'No filter settings to remove'
                    
        elif params[0] == '3pf':
            if self._pst != None:
                self._st = self._pst.copy()
            for tr in self._st:
                tr.data = threepointfilter(tr.data, int(params[1]))
            self._filtered = " ".join(params).upper()
            logmsg = '3-point filter recursively applied' \
                     ' %d times' % (int(params[1]))
            
        if auto: 
            self._logger.autoinfo(logmsg)
        else:
            self._logger.info(logmsg)
            
        if self._ax is not None:
            self._ax.cla()
            self._draw_seismograms()

    def _load_picks(self, load_all=False):
        """Load picks for loaded picksets from database.
        
        Parameters
        ----------
        load_all: bool
            If True the entire fbpicks table is loaded.
        """
        
        if load_all:
            cmd = """SELECT *
                     FROM fbpicks"""
            self._picks = self.slh.read_data(cmd)
        else:
            if self._st != None:
                wps = ', '.join(['\'%s\'' % (e) for e in self._picksets])
                wselin = ', '.join(['\'%d %d\'' % (s, r) 
                    for s, r in zip(self._selin.sin, self._selin. rin)])
                cmd = """SELECT *
                         FROM fbpicks
                         WHERE pickset IN (%s) AND
                               shot_index_number || " " || 
                               receiver_index_number IN
                               (%s)""" % (wps, wselin)
                self._picks = self.slh.read_data(cmd)
            
    def _save_picks(self, saveall=False):
        """Save traveltimes of currently selected traces in the database."""
        
        for ps in self._picksets:
            picks = self._picks[self._picks.pickset==ps]
            for i, pck in picks.iterrows():
                cmd = """UPDATE fbpicks
                         SET traveltime = %.8f
                         WHERE shot_index_number == %d AND
                               receiver_index_number == %d AND
                               pickset == \'%s\'""" % \
                               (pck.traveltime,
                                pck.shot_index_number,
                                pck.receiver_index_number,
                                ps)
                with self.slh.dbc:
                    self.slh.dbc.execute(cmd)

    def _activate_pickset(self, params):
        """Make the given pickset the active one, i.e., picks will 
        be stored in this pickset.
        
        Parameters
        ----------
        params: numpy array
            Array of length 1 holding the name (str) of the pickset to be 
            made the activate one.
        """
        
        if not self._geomapp:
            self._logger.critical('Geometry not applied')
            return
            
        if len(params) != 1:
            self._logger.error('Insufficient number of arguments')
            return
        
        if params[0] not in self._picksets:
            self._load_pickset(params)
            
        self._activeps = params[0]
        # ~ self._picking = True
        self._procmode = PROC_MODES.pick
        self._logger.info('\'%s\' set as active pickset' % (params[0]))

    def _copy_pickset(self, params):
        """Create a copy of a give pickset.
        
        Parameters
        ----------
        params: numpy array
            Array holding the names (str) of the source and the destination
            pickset.
        """
        if len(params) != 2:
            self._logger.error('Insufficient number of parameters')

        cmd = """SELECT *
                 FROM fbpicks
                 WHERE pickset==\'%s\'""" % (params[0])
        src = self.slh.read_data(cmd)
        
        if src.empty:
            self._logger.error('Source pickset \'%s\' does not exist' % 
                               (params[0]))
            return
        
        self._create_pickset(params[1::])
        src['pickset'] = params[1]
        
        self._load_picks(load_all=True)
        self._picks = self._picks[self._picks.pickset!=params[1]]
        self._picks = pd.concat([self._picks, src], axis=0,
                                ignore_index=True)
        self.slh.write_data('fbpicks', self._picks, mode='replace')
        
        self._logger.info('Copied pickset \'%s\' to \'%s\'' % 
            (params[0], params[1]))
            
        # ~ self._load_picks()

    def _create_pickset(self, params):
        """Create a new pickset in the project file.
        
        Parameters
        ----------
        params: numpy array
            Array of length 1 holding the name (str) of the pickset to be 
            created.
            Create multiple picksets in future version?
        """
        
        if len(params) > 1:
            self._logger.error('Insufficient number of parameters')
            return
        
        if params[0] == 'picks':
            self._logger.critical('Cannot create pickset with the name of ' \
                'the default pickset (\'%s\')' % (params[0]))
            return
        
        cmd = """SELECT DISTINCT pickset
                 FROM fbpicks
                 WHERE pickset==\'%s\'""" % (params[0])
        df = self.slh.read_data(cmd)
        
        if not df.empty:
            self._logger.error('Pickset \'%s\' already exists' % 
                               (params[0]))
            return 0
        
        cmd = """SELECT shot_index_number, receiver_index_number
                 FROM applied_geometry"""
        fbpicks = self.slh.read_data(cmd)
        fbpicks['pickset'] = params[0]
        fbpicks['traveltime'] = -1
        
        self.slh.write_data('fbpicks', fbpicks)
            
        self._logger.info('Created new pickset \'%s\'' % (params[0]))
        # ~ self._load_pickset(params)
        self._activate_pickset(params)
        
        return 1

    def _delete_pickset(self, params):
        """Delete one or more picksets from the project file.
        
        Parameters
        ----------
        params: numpy array
            Array holding the names (str) of the picksets to be deleted from 
            the project file.
        """
        if len(params) > 1:
            self._logger.error('Insufficient number of parameters')
            return
                
        if params[0] == "picks":
            self._logger.critical('Cannot delete default pickset (\'picks\')')
            return
            
            
        if params[0] in self._picksets:
            self._logger.critical('Cannot delete a loaded pickset')
            return
        
        cmd = """DELETE FROM fbpicks
                 WHERE pickset==\'%s\'""" % (params[0])
        with self.slh.dbc:
            self.slh.dbc.execute(cmd)
        self._logger.info('Pickset \'%s\' deleted' % (params[0]))

    def _export_pickset(self, params):
        """Export picks from a given pickset to pck file.
        
        Parameters
        ----------
            params: numpy array
                Array of length 1 holding the name (str) of the pickset to be 
                exported.
                Export multiple picksets in future version?
        """
        
        if len(params) != 1:
            self._logger.error('Insufficient number of parameters')
            return
        
        # get picks from database
        cmd = """SELECT s.station_id ssid, r.station_id rsid, fbp.traveltime tt
                 FROM fbpicks fbp
                    LEFT JOIN shots s
                        ON fbp.shot_index_number==s.shot_index_number
                    LEFT JOIN receivers r
                        ON fbp.receiver_index_number==r.receiver_index_number
                 WHERE fbp.pickset==\'%s\'""" % (params[0])
        expcks = self.slh.read_data(cmd)
        
        if expcks.empty:
            self._logger.error('Pickset \'%s\' does not exist' % 
                               (params[0]))
            return
        
        # get station geometry
        cmd = """SELECT x, y, z
                 FROM geometry"""
        stations = self.slh.read_data(cmd)
            
        # create file name and path
        pickfile = params[0] + '_' \
            + dt.utcnow().strftime('%Y%m%dT%H%M%S') + '.pck'
                   
        # create subdirectory to store pck files (if not exists)
        pickpath = os.path.join(self._procdir, 'picks')
        if not os.path.exists(pickpath): os.makedirs(pickpath)
        
        with open(os.path.join(pickpath, pickfile), 'w') as f:
            # number of stations
            f.write('%d\n' % (stations.shape[0]))
            
            # header for sensor positions
            f.write('# x y z\n')
            
            # receiver positions
            for idx, stat in stations.iterrows():
                f.write('%.3f\t%.3f\t%.3f\n' % (stat[0], stat[1], stat[2]))

            validpcks = expcks[expcks.tt>=0]
            f.write('%d\n' % (validpcks.shape[0]))
            
            # header for first break picks
            f.write('# s g t valid\n')
            
            for idx, p in validpcks.iterrows():
                f.write('%03d\t%03d\t%.8f\t%d\n' % (p.ssid, p.rsid, p.tt, 1))
                  
            # EOF marker
            f.write('0\n')
        
        self._logger.info('pickset \'%s\' saved to %s' % 
            (params[0], pickfile))

    def _import_pickset(self, params):
        """Import picks from a pck file to a new pickset.
        
        Parameters
        ----------
        params: numpy array
            Array holding the name of the pck file and the pickset.
        """
        
        if len(params) != 2:
            self._logger.error('Insufficient number of parameters')
            return
            
        # check if pick file exists/is a file
        pckf = os.path.join(os.path.join(self._procdir, 'picks'), params[0])
        
        if not os.path.isfile(pckf):
            self._logger.critical('Pick file \'%s\' not found' % (params[0]))
            return
            
        # read pick data
        imp = Refraction(pckf)
        srtData = pg.DataContainer(pckf, 
                                   sensorTokens='g s',
                                   removeInvalid=False)
        data = np.vstack((np.array(srtData['s']) + 1,
                          np.array(srtData['g']) + 1,
                          np.array(srtData['t']))).T
        
        # check if number of picks is (theoretically) valid
        if data.shape[0] > self._numpicks:
            self._logger.error('Number of picks to import exceeds the ' \
                'geometry based number of picks')
            return
            
        if not self._create_pickset([params[1]]):
            return
        
        imppicks = pd.DataFrame(data, columns=['shot_index_number',
                                              'receiver_index_number',
                                              'traveltime'])
        imppicks['pickset'] = params[1]
        
        # get shot_index_number from table shots
        cmd = """SELECT shot_index_number sin, station_id sid
                 FROM shots"""
        shts = self.slh.read_data(cmd)
        smap = dict(zip(shts.sid, shts.sin))
        imppicks['shot_index_number'] = imppicks['shot_index_number'].map(smap)
        
        # get receiver_index_number from table receivers
        cmd = """SELECT receiver_index_number rin, station_id sid
                 FROM receivers"""
        recs = self.slh.read_data(cmd)
        rmap = dict(zip(recs.sid, recs.rin))
        imppicks['receiver_index_number'] = \
            imppicks['receiver_index_number'].map(rmap)
            
        cmd = """SELECT *
                 FROM fbpicks
                 WHERE pickset=\'%s\'""" % (params[1])
        dbpicks = self.slh.read_data(cmd)
        
        newpicks = dbpicks.merge(imppicks, 'outer')
        newpicks = newpicks.sort_values(
            'traveltime', ascending=False).drop_duplicates(
                ['shot_index_number', 'receiver_index_number'], 
                keep='first').sort_values(
                    ['shot_index_number', 
                     'receiver_index_number']).reset_index(drop=True)

        self._load_picks(load_all=True)
        self._picks = self._picks[self._picks.pickset!=params[1]]
        self._picks = pd.concat([self._picks, newpicks], axis=0,
                                ignore_index=True)
        
        self.slh.write_data('fbpicks', self._picks, mode='replace')
        
        self._logger.info('Imported \'%s\' to pickset \'%s\'' % 
                          (params[0], params[1]))

    def _load_pickset(self, params):
        """Load one or more picksets from the project file.
        
        Parameters
        ----------
        params: numpy array or list
            Names (str) of the picksets to be loaded from the project file.
        """
        
        db_picksets = self.slh.read_data("""SELECT DISTINCT pickset
                                            FROM fbpicks""")
        
        for p in params:
            if p in self._picksets:
                self._logger.info('Pickset \'%s\' already loaded' % (p))
                continue
            
            # check if pickset exists in database
            if (db_picksets.pickset==p).sum():
                self._picksets.append(p)
                self._logger.info('Pickset \'%s\' loaded' % (p))
            else:
                self._logger.warn('Pickset \'%s\' does not exist' % (p))
        
        # Seems to work but looks like a quite complicated solution...
        if self._selected != "":
            if not self._filterhold:
                self._filterhold = True
                self.select(self._selected.lower(), auto=True)
                self._filterhold = False
            else: self.select(self._selected.lower(), auto=True)

    def _print_picksets_info(self):
        """Print information about all picksets saved in the project file."""
        
        if not self._geomapp:
            self._logger.error('Geometry not applied')
            return
        
        db_picksets = self.slh.read_data("""SELECT DISTINCT pickset
                                            FROM fbpicks""")
        
        print('%s\t%s\t%s' % ('pickset   ','loaded    ','active    '))
        print('------------------------------------------')
        
        for ps in db_picksets.pickset:
            l = ps in self._picksets
            a = ps == self._activeps
            
            print('%s\t%s\t%s' % ('{txt: <{w}}'.format(txt=ps,
                                                       w=10),
                                  '{txt: <{w}}'.format(txt='Y' 
                                                            if l else 'N',
                                                       w=10),
                                  '{txt: <{w}}'.format(txt='Y'
                                                            if a else 'N',
                                                       w=10)))           

    def _rename_pickset(self, params):
        """Assign a new name to a given pickset.
        
        Parameters
        ----------
         params: numpy array
            Array holding the old and new name (str) of the pickset.       
        """
        
        if len(params) != 2:
            self._logger.error('Insufficient number of parameters')
            return
    
        if params[0] == 'picks':
            self._logger.error('Cannot rename default pickset (\'picks\')')
            return
        elif params[0] == 'picks':
            self._logger.critical('Cannot rename to name of default ' \
                                  'pickset (\'picks\')')
            return
        
        if params[0] in self._picksets:
            self._logger.error('Cannot rename a currently loaded pickset')
            return
            
        db_picksets = self.slh.read_data("""SELECT DISTINCT pickset
                                            FROM fbpicks""")
                                            
        if not (db_picksets.pickset==params[0]).sum():
            self._logger.error('pickset \'%s\' does not exist' % 
                               (params[0]))
            return
        
        if (db_picksets.pickset==params[1]).sum():
            self._logger.error('Pickset \'%s\' already exists' % 
                               (params[1]))
            return
            
        cmd = """UPDATE fbpicks
                 SET pickset = \'%s\'
                 WHERE pickset == \'%s\'
              """ % (params[1], params[0])
        with self.slh.dbc:
            self.slh.dbc.execute(cmd)
        
        self._logger.info('Info: pickset \'%s\' renamed to \'%s\'' % 
                          (params[0], params[1]))
                
    def _unload_pickset(self, params):
        """Unload one or more picksets from the MethodManager object.
        
        Parameters
        ----------
        params: numpy array
            Names (str) of the picksets to be unloaded from the object.        
        """
        
        for p in params:
            if p == self._activeps:
                self._logger.error('Cannot unload the active pickset')
                continue
                
            if p not in self._picksets:
                continue
                
            self._picksets.remove(p)
            
            # is this if condition really necessary?
            if p in self._picks: del self._picks[p]
            
            self._logger.info('Pickset \'%s\' removed from workflow' % (p))
        
        # check difference to _load_pickset()!    
        if self._selected != "":
            self.select(self._selected.lower(), auto=True)

    def _parse_pickset_params(self, options):
        """Check the basic validity of the picksets parameters and return them
        as numpy array.
        
        Parameters
        ----------
        options: str, default ''
            Options for accessing the picksets.
        
        Returns
        -------
        params: numpy array
            Parameters for further processing.
        """
            
        if not isinstance(options, str):
            self._logger.error('Parameter of datatype str expected')
            return []
            
        if options == "":
            return np.array([""])
        
        options = options.lower().rstrip()    
        params = options.split(" ")
        
        if params[0] not in PICKSET_KEYWORDS:
            self._logger.error('Invalid keyword given')
            return []
            
        return params

    def picksets(self, options=""):
        """Work with and manipulate picksets.
        
        Parameters:
        options: str
            Parameters used for working with picksets. The string has to
            contain a valid keyword ('create', 'load', 'unload', 'rename',
            'use', 'copy', 'rename','export', 'import') at the beginning
            followed by the required parameters.
            TODO: describe parameters for each keyword.
        """
        
        self._logger.input('picksets ' + options)
        
        if not self._check_processing_ready():
            return
            
        params = self._parse_pickset_params(options)
        if len(params) == 0: return 
        
        if params[0] == "":
            self._print_picksets_info()
        elif params[0] == 'load':
            self._load_pickset(params[1::])
        elif params[0] == 'create':
            self._create_pickset(params[1::])
        elif params[0] == 'unload':
            self._unload_pickset(params[1::])
        elif params[0] == 'delete':
            self._delete_pickset(params[1::])
        elif params[0] == 'copy':
            self._copy_pickset(params[1::])
        elif params[0] == 'rename':
            self._rename_pickset(params[1::])
        elif params[0] == 'use':
            self._activate_pickset(params[1::])
        elif params[0] == 'export':
            self._export_pickset(params[1::])
        elif params[0] == 'import':
            self._import_pickset(params[1::])

    def _plot_pickpercentage(self, **kwargs):
        """Plot the picking progress in the active pickset.
        
        Parameters
        ----------
        **kwargs: arbitrary keyword arguments
            Keyword arguments forwarded to the underlying plotting methods of
            the matplotlib, e.g., 'fontsize', 'color'         
        """
        cmap = kwargs.pop('cmap', kwargs.pop('cMap', 'Greens'))
        c = kwargs.pop('color', kwargs.pop('c', 'k'))
        
        cmd = """SELECT pickset,
                        shot_index_number sin,
                        receiver_index_number rin,
                        traveltime
                 FROM fbpicks
                 WHERE pickset == \'%s\'""" % (self._activeps)
        fbpicks = self.slh.read_data(cmd)
       
        # compute picking progress for each shot
        perc = np.empty((len(self._data), 2))
        
        for i, sin in enumerate(fbpicks.sin.unique()):
            tt = fbpicks[fbpicks.sin == sin].traveltime
            perc[i, 0] = sin
            perc[i, 1] = tt[tt>=0].count()/tt.count()
        
        perc = perc[np.argsort(perc[:, 0])]
        
        # create mapper for color scheme
        norm = mpl.colors.Normalize(vmin=0, vmax=1)
        mapper = mpl.cm.ScalarMappable(norm=norm, cmap=cmap) 
        
        # plot picking progress
        self._ppfig, self._ppax = plt.subplots(1, figsize=(8, 2))

        self._ppax.bar(x=perc[:, 0],
                       height=perc[:, 1],
                       width=.8,
                       color=mapper.to_rgba(perc[:, 1]),
                       edgecolor=".2")
        
        self._ppax.set_xlim((np.min(perc[:, 0] - .5),
                             np.max(perc[:, 0] + .5)))
        self._ppax.set_ylim((0, 1))
        
        self._ppax.xaxis.tick_bottom()
        self._ppax.set_yticks([0, 1])
        self._ppax.set_yticklabels([0, 100])

        self._ppax.spines['top'].set_visible(False)
        self._ppax.spines['right'].set_visible(False)

        self._ppax.set_xlabel('SIN')
        self._ppax.set_ylabel('(%)')
        # ~ self._ppax.set_title('Picking percentage', fontsize=self._fontsize)

        plt.tight_layout()
        plt.show(block=False)

    def _save_3d_pseudosection(self):
        """Internal function for saving screenshots of 3d pseudosections."""
        
        self.ps3d.ren_win.OffScreenRenderingOn()
        self.ps3d.enable_anti_aliasing()
        self.ps3d.screenshot('~/try.png', window_size=(10000, 10000))

    def _plot_pseudosection(self, **kwargs):
        """Plot the apparent velocites as pseudosection.
        
        Parameters
        ----------
        **kwargs: arbitrary keyword arguments
            Keyword arguments forwarded to the underlying plotting methods of
            the matplotlib, e.g., 'fontsize', 'color'    
        """
        
        self._check_processing_ready()
        
        ps = kwargs.pop('point_size', 13)
        
        # set colormap
        cmap = kwargs.pop('cmap', kwargs.pop('cMap', 'Spectral_r'))
        
        # get information from database
        cmd = """SELECT ROUND(absolute_offset/traveltime, 1) vapp,
                        midpoint_x mp_x, midpoint_y mp_y,
                        ROUND(absolute_offset/3, 1) pd
                 FROM fbpicks fbp, applied_geometry ag
                 ON fbp.shot_index_number == ag.shot_index_number AND
                    fbp.receiver_index_number == ag.receiver_index_number
                 WHERE traveltime > 0 AND
                       pd > 0 AND
                       pickset == \'%s\'""" % (self._activeps)
        psec = self.slh.read_data(cmd)
        
        cmd = """SELECT x, y, z
                 FROM geometry"""
        stations = self.slh.read_data(cmd)
        stations.z += psec.pd.max()*.05
        
        vmin = kwargs.pop('vmin', 
            kwargs.pop('vMin', psec.vapp.quantile(.05)))
            # ~ kwargs.pop('vMin', psec.vapp.min()))
        vmax = kwargs.pop('vmax', 
            kwargs.pop('vMax', psec.vapp.quantile(.95)))
            # ~ kwargs.pop('vMax', psec.vapp.max()))
        ymax = psec.pd.max()

        # find unique mp-pd pairs
        mult = psec[['mp_x', 'mp_y', 'pd']].to_numpy()
        u, idx, inv, cnt = np.unique(mult,
                                     return_counts=True,
                                     return_index=True,
                                     return_inverse=True,
                                     axis=0)
                                     
        if self._3d:
            import pyvista as pv
            import pyvistaqt as pvqt
            
            psec['pd'] *= -1
            
            # configure pyvista
            thm = pv.themes.DocumentTheme()
            thm.font.family = 'arial'
            thm.font.color = 'black'
            thm.font.title_size = self._fontsize + 4
            thm.font.label_size = self._fontsize + 4
            thm.cmap = cmap
            thm.colorbar_horizontal.width = 0.4
            thm.colorbar_horizontal.position_x = 0.55
            pv.global_theme.load_theme(thm)
            
            # colorbar format
            cbargs = dict(
                # ~ title=r'Apparent velocity (m/s)\n',
                title="""Apparent velocity (m/s)
                      """,
                title_font_size=self._fontsize + 4,
                label_font_size=self._fontsize + 2,
                fmt='%.1f'
                )
            
            # plot with pyvista
            self.ps3d = pvqt.BackgroundPlotter(**kwargs)
            self.ps3d.show_bounds(grid=False, location='back')
            for i in np.arange(np.max(cnt)) + 1:
                pts = pv.PolyData(
                        (psec.iloc[idx[cnt == i]][['mp_x', 
                                                   'mp_y',
                                                   'pd']]).values)
                pts['vapp'] = (psec.iloc[idx[cnt==i]]['vapp']).values
                pts.set_active_scalars('vapp')
                self.ps3d.add_mesh(
                   pts, 
                   point_size=ps, clim=[vmin, vmax],
                   show_edges=True, render_points_as_spheres=True,
                   scalar_bar_args=cbargs)
                if i > 1:
                    # get idices of multiple pairs
                    mi = [np.where(np.logical_and(
                                mult[:,0]==mult[j, 0], \
                                mult[:,1]==mult[j, 1]) == 1)[0][-1] \
                            for j in idx[cnt == i]]
                    pts = pv.PolyData(
                        (psec.iloc[mi][['mp_x', 'mp_y', 'pd']]).values)
                    pts['vapp'] = (psec.iloc[mi]['vapp']).values
                    pts.set_active_scalars('vapp')
                    self.ps3d.add_mesh(
                       pts, 
                       point_size=ps / i ** 3, clim=[vmin, vmax],
                       show_edges=True, render_points_as_spheres=True,
                       scalar_bar_args=cbargs) 
            
            stats = pv.PolyData(stations.values)
            self.ps3d.add_mesh(
                stats, point_size=ps-3, color='#373737', 
                render_points_as_spheres=True,
                show_scalar_bar=False)
            
            labels = dict(
                zlabel='Pseudodepth (m)', xlabel='X (m)', ylabel='Y (m)')
            self.ps3d.show_grid(**labels)
            self.ps3d.show()
        else:
            norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
            self._vafig, self._vaax = plt.subplots(1, figsize=(10, 4.5))
            
            # plot pseudosection
            for i in np.arange(np.max(cnt)) + 1:
                scat = self._vaax.scatter(psec.iloc[idx[cnt == i]]['mp_x'],
                                    psec.iloc[idx[cnt == i]]['pd'],
                                    c=psec.iloc[idx[cnt==i]]['vapp'], s=40,
                                    cmap=cmap,
                                    norm=norm,
                                    **kwargs)
                if i > 1:
                    # get idices of multiple pairs
                    mi = [np.where(np.logical_and(
                        mult[:,0]==mult[j, 0], \
                        mult[:,1]==mult[j, 1]) == 1)[0][-1]\
                            for j in idx[cnt == i]]
                            
                    self._vaax.scatter(psec.iloc[mi]['mp_x'], 
                                       psec.iloc[mi]['pd'],
                                 c=psec.iloc[mi]['vapp'], s=40 / i ** 3,
                                 cmap=cmap, 
                                 norm=norm,
                                 **kwargs)

            axcb = self._vaax.inset_axes(
                bounds=[.05, .02, .2, .03],
                transform=self._vaax.transAxes)
                                
            cb = plt.colorbar(scat, cax=axcb,
                              orientation='horizontal')
            cb.set_label('Apparent velocity (m/s)', fontweight='bold')
            axcb.xaxis.set_ticks_position('top')
            axcb.xaxis.set_label_position('top')

            self._vaax.invert_yaxis()
            self._vaax.xaxis.set_ticks_position('bottom')

            self._vaax.spines['top'].set_visible(False)
            self._vaax.spines['right'].set_visible(False)

            self._vaax.set_xlabel('x (m)')
            self._vaax.set_ylabel('Pseudodepth (m)')

            plt.tight_layout()
            plt.show(block=False)

    def _draw_estimated_velocities(self):
        """Draw the velocity estimation lines."""
        
        if self._vellines is not None:
            # ~ self._vellines.remove()
            for a in self._vellinetext:
                a.remove()
            self._vellinetext = []
            
            self._ax.add_collection(self._vellines)
            
            for seg in self._vellines.get_segments():
                
                # get geometry information from database
                if self._cos:
                    geom = np.vstack(([self._aoffs[int(round(seg[0, 0], 0))],
                                       0, 0],
                                      [self._aoffs[int(round(seg[1, 0], 0))],
                                       0, 0]))
                else:
                    if self._xlabel == 'SIN':
                        cmd = """SELECT x, y, z
                                 FROM shots s JOIN geometry g
                                 ON s.station_id == g.station_id
                                 WHERE shot_index_number in (%d, %d)""" % \
                                 (self._picks.iloc[int(round(seg[0, 0], 
                                    0))].shot_index_number, 
                                  self._picks.iloc[int(round(seg[1, 0], 
                                    0))].shot_index_number)
                    else:
                        cmd = """SELECT x, y, z
                                 FROM receivers r JOIN geometry g
                                 ON r.station_id == g.station_id
                                 WHERE receiver_index_number in (%d, %d)""" % \
                                 (self._picks.iloc[int(round(seg[0, 0], 
                                    0))].receiver_index_number, 
                                  self._picks.iloc[int(round(seg[1, 0], 
                                    0))].receiver_index_number)
                    geom = self.slh.read_data(cmd).to_numpy()
                    
                # compute velocity
                d = ((geom[0, 0] - geom[1, 0])**2 +
                     (geom[0, 1] - geom[1, 1])**2 +
                     (geom[0, 2] - geom[1, 2])**2)**(1/2)
                vel = d/(seg[0, 1] - seg[1, 1])
                avel = np.around(abs(vel) / 10) * 10
                
                # plot velocity as text
                if avel < 100000:
                    text = '%d $\mathregular{m\/s^{-1}}$' % (avel)
                else:
                    text = 'inf.'
                
                xt = (seg[1, 0] + seg[0, 0]) / 2
                yt = (seg[1, 1] + seg[0, 1]) / 2
                
                self._vellinetext.append(
                    self._ax.annotate(text,
                                      xy=(xt, yt - yt * .02),
                                      textcoords='data',
                                      color='k', fontweight='bold',
                                      bbox=dict(boxstyle="round", 
                                                fc=(1.0, 1.0, 1.0),
                                                ec=".7",
                                                alpha=.7),
                                      horizontalalignment= \
                                        'left' if vel > 0 else "right",
                                      verticalalignment='center'))
                
            self._fig.canvas.draw()

    def _draw_pickline(self):
        """Draw the line used for setting multiple picks at once."""
        
        if self._picklineplot is not None:
            self._picklineplot.remove()
        
        self._picklineplot, =  self._ax.plot(self._pickline[0],
                                             self._pickline[1],
                                             c='darkgray',
                                             lw=2.,
                                             solid_capstyle='round')
        self._fig.canvas.draw()

    def _draw_picks(self, **kwargs):
        """Draw the picks for the currently selected seismograms.
        
        Parameters
        ----------
        **kwargs: arbitrary keyword arguments
            Keyword arguments forwarded to the underlying plotting methods of
            the matplotlib, e.g., 'fontsize', 'color'
        """
        
        # ~ print(self._pickscatter)
        
        # ~ picks = self._picks[self._picks.pickset==self._activeps]
        # ~ if self._pickscatter is None:
            # ~ print('first time')
            # ~ self._pickscatter = self._ax.scatter(
                             # ~ np.arange(len(picks.traveltime)), 
                             # ~ picks.traveltime,
                             # ~ marker=self._marker,
                             # ~ color='g',
                             # ~ label=self._activeps,
                             # ~ picker=True,
                             # ~ animated=True,
                             # ~ **kwargs)
            # ~ self._ax.legend(loc='lower right')
        # ~ else:
            # ~ print('on pick')
            # ~ self._fig.canvas.restore_region(self._bg)
            # ~ self._pickscatter.set_offsets(np.vstack((
                                            # ~ np.arange(len(picks.traveltime)),
                                            # ~ picks.traveltime
                                            # ~ )).T)
        
        
            # ~ self._ax.draw_artist(self._pickscatter)
            # ~ self._fig.canvas.blit(self._fig.bbox)
            # ~ self._fig.canvas.flush_events()
        
        if self._pickscatter is not None:
            for e in self._pickscatter: 
                e.remove()
        
        self._pickscatter = []
        
        cmap = plt.cm.tab10(np.arange(10))
        
        lgndelems = []
        for i, ps in enumerate(self._picksets):
            c = 'g' if ps == self._activeps else cmap[i % 10]
            
            picks = self._picks[self._picks.pickset==ps]
            
            self._pickscatter = np.append(self._pickscatter, 
                self._ax.scatter(np.arange(len(picks.traveltime)), 
                                 picks.traveltime,
                                 marker=self._marker,
                                 color=c,
                                 label=ps,
                                 picker=True,
                                 # ~ animated=True,
                                 **kwargs))
        
        self._ax.legend(loc='lower right')
        self._fig.canvas.draw()

    def _draw_seismograms(self, **kwargs):
        """Do the actual plotting.
        
        Parameters
        ----------
        **kwargs: arbitrary keyword arguments
            Keyword arguments forwarded to the underlying plotting methods of
            the matplotlib, e.g., 'fontsize', 'color' 
        """
        
        st = self._st
        
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            st.normalize()
        
        self._title = self._mode \
           + ("SELECT: " + self._selected if self._selected != ""
                                          else self._selected) \
           + (" | FILTER: " + self._filtered if self._filtered != "" 
                                             else self._filtered)
        
        t = np.arange(st[0].stats.npts) / st[0].stats.sampling_rate
        
        if not self._pvd:
            cmd = """SELECT receiver_index_number rin, polarity
                     FROM receivers"""
            recs = self.slh.read_data(cmd)
        
        # plot seismograms
        X, Y = [], []
        fill = []
        plotset = np.empty((0, 3))
        vd = np.empty((len(t), len(st)))
        xtls = np.empty(len(st), dtype=int)
        for i, tr in enumerate(st):
			# prepare trace data for plotting
            trd = tr.data.copy()
            
            # ~ f = np.exp(1/1000 * t)
            # ~ trd *= (f.max() - f)
            # ~ trd /= (max(pg.abs(trd)))
            
            if not self._pvd:
                pol = int(recs[recs.rin==tr.stats.RIN].polarity)
                trd *= pol
            trd = trd * self._scaling + i
            
            # clip amplitudes
            trd[trd > .45 + i] = .45 + i; trd[trd < i-.45] = i-.45
            
            X.append(list(trd))
            Y.append(list(t))
            
            if self._plotmode == PLOT_MODES.seismograms:
                self._psg_fb.append((
                    self._ax.fill_betweenx(t, i, trd,
                            where=trd > i,
                            color="dodgerblue",
                            alpha=.5, lw=.5),
                    self._ax.fill_betweenx(t, i, trd,
                            where=trd < i,
                            color="tomato",
                            alpha=.5, lw=.5)))
            
            # add ticklabel
            if self._pvd:
                if self._syndata:
                    xtls[i] = int(tr.stats.channel)
                else:
                    xtls[i] = int(tr.stats.seg2['CHANNEL_NUMBER'])
            elif self._xlabel == 'SIN': xtls[i] = int(tr.stats.SIN)
            elif self._xlabel == 'RIN': xtls[i] = int(tr.stats.RIN)
            elif self._xlabel == 'AOFFSET': xtls[i] = int(tr.stats.xtl)
            else: xtls[i] = int(tr.stats.RIN)
            
            if not self._pvd:
                # add element to picking dictionary
                plotset = np.vstack((plotset,
                                    np.array([i,
                                              tr.stats.SIN,
                                              tr.stats.RIN])))
            
            # populate variable density array
            vd[:, i] = trd - i
        
        self._plotset = pd.DataFrame(plotset,
                                     columns=['idx', 'sin', 'rin'])

        # add seismograms to plot
        seismograms = LineCollection((list(zip(x, y)) for x, y in zip(X, Y)),
                                     colors='k', linestyles='-',
                                     linewidths=.5)
        seismograms.set_picker(True)
        self._ax.add_collection(seismograms)
        
        # ~ self._psg_sg = LineCollection((list(zip(x, y)) for x, y in zip(X, Y)),
                                     # ~ colors='k', linestyles='-',
                                     # ~ linewidths=.5)
        # ~ self._psg_sg.set_picker(True)
        # ~ self._ax.add_collection(self._psg_sg)
        
        # show variable density information
        if self._plotmode == PLOT_MODES.var_dens:
            self._psg_isvd = self._ax.imshow(vd,
                        extent=(-.5, i + .5, t[-1], t[0]),
                        cmap='seismic_r',
                        aspect='auto',
                        vmin=-1, vmax=1,
                        interpolation='bicubic')
            
        self._ax.set_xlim((-1, i + 1))
        self._ax.set_ylim(self._ylim)
        
        self._ax.set_xlabel(self._xlabel)
        self._ax.xaxis.set_label_position('top')
        self._ax.set_ylabel('Time (s)')
        
        self._ax.set_xticks(np.arange(0, len(st))[0:-1:5])
        self._ax.set_xticklabels(xtls[0:-1:5])
        
        if not self._pvd:
            self._draw_estimated_velocities()
            
        if not self._cos and not self._pvd:
            self._draw_picks(**kwargs)
        
        plt.title(self._title, fontsize=self._fontsize, loc='left')
        
        self._fig.canvas.draw()
        
    def _draw_velline(self):
        """Draw the velocity estimation line that is currently being
        created.
        """
        
        if self._vellineplot is not None:
            self._vellineplot.remove()
        
        if self._velline is not None:
            self._vellineplot, = self._ax.plot(self._velline[0],
                                               self._velline[1],
                                               c='darkgray',
                                               lw=2.,
                                               solid_capstyle='round')    
        
        self._fig.canvas.draw()

    def _onclose_seismograms(self, event):
        """Cleanup after closing the plot figure."""
        
        if not self._pvd:
            self._save_picks()
        self._title = ''
        self._ylim = (self._st[0].stats.npts / \
            self._st[0].stats.sampling_rate, 0)
        self._scaling = 1
        self._pickscatter = None
        self._picklineplot = None
        self._vellineplot = None
        self._vellines = None
        self._velline = None
        self._vellinetext = []
        # ~ self._procmode = PROC_MODES.inactive
        self._fig = None
        self._ax = None
        
    def _onkeyrelease_seismograms(self, event):
        """Process key input in plot figure."""
        
        if event.key == 'v' and not self._pvd:
            if self._procmode == PROC_MODES.vel_est:
                if self._cos:
                    self._procmode = PROC_MODES.inactive
                else:
                    self._procmode = PROC_MODES.pick
            else:
                self._procmode = PROC_MODES.vel_est
            
        elif event.key == 'a':
            if self._scrollmode == SCROLL_MODES.zoom:
                self._scrollmode = SCROLL_MODES.amplitude
            else:
                self._scrollmode = SCROLL_MODES.zoom
            
        elif event.key == 'r' and not (self._pvd or self._cos):
            if self._procmode == PROC_MODES.trace_reverse:
                self._procmode = PROC_MODES.pick
            else:
                self._procmode = PROC_MODES.trace_reverse
            
        elif event.key == 'm' and not (self._pvd or self._cos):
            if self._procmode == PROC_MODES.trace_mute:
                self._procmode = PROC_MODES.pick
            else:
                self._procmode = PROC_MODES.trace_mute
            
        elif event.key == 'up' or event.key == 'down':
            it = 1 if event.key == 'up' else -1
            self._plotmode += it
            if self._plotmode < PLOT_MODES.min():
                self._plotmode = PLOT_MODES.max()
            elif self._plotmode > PLOT_MODES.max():
                self._plotmode = PLOT_MODES.min()
            
            # update the plot window
            self._ax.cla()
            self._draw_seismograms()
            self._fig.canvas.draw()
            
        elif event.key == 'left' or event.key == 'right':
            inc = 1 if event.key == 'right' else -1

            if self._pvd:
                keys = list(self._data.keys())
                
                # determine next valid index
                newkey = self._pvk + inc
                if newkey < min(keys):
                    newkey = max(keys)
                elif newkey > max(keys):
                    newkey = min(keys)
                
                # select corresponding shot data
                self._pvk = newkey
                self._st = self._data[self._pvk]
                self._pst = self._st.copy()
                self._selected = self._fnames[self._pvk]
                
                # apply filter
                if self._filterhold and self._filtered != "":
                    self.filter(self._filtered.lower(), auto=True)
                else:
                    self.filter('remove', auto=True)
                    
            else:
                by = (self._selected.split(" ")[0]).lower()
                curin = float(self._selected.split(" ")[1])
                newin = curin + inc
                if by == 'sin':
                    cmd = """SELECT shot_index_number
                             FROM shots"""
                    sins = self.slh.read_data(cmd).to_numpy()
                    
                    if newin < sins.min(): newin = sins.max()
                    elif newin > sins.max(): newin = sins.min()
                    
                elif by == 'rin':
                    cmd = """SELECT receiver_index_number
                             FROM receivers"""
                    rins = self.slh.read_data(cmd).to_numpy()
                    
                    if newin < rins.min(): newin = rins.max()
                    elif newin > rins.min(): newin = rins.min()
                
                if not self._pvd:
                    self._save_picks()
                self.select(by + " " + str(int(newin)), auto=True)
                
                self._pickscatter = None
                self._picklineplot = None
                self._vellineplot = None
                self._vellines = None
                self._velline = None
                self._vellinetext = []
            
            # update the plot window
            self._ax.cla()
            self._draw_seismograms()
            self._fig.canvas.draw()
        
        # update status bar text
        try:
            self._fig.canvas.motion_notify_event(
                *self._ax.transData.transform_point([event.xdata,
                                                     event.ydata]))
        except TypeError as e: pass

    def _onmove_seismograms(self, event):
        """Process mouse move in plot figure."""
        
        if (self._fig.canvas.manager.toolbar.mode == "" and
            self._pressed):
            self._moved = True
                
            if (self._procmode == PROC_MODES.pick and
                event.button in [1, 3]):
                self._pickline[0][1] = event.xdata
                self._pickline[1][1] = event.ydata
                self._draw_pickline()
            elif (self._procmode == PROC_MODES.vel_est and
                  event.button == 1):
                self._velline[0][1] = event.xdata
                self._velline[1][1] = event.ydata
                self._draw_velline()

    def _onpick_seismograms(self, event):
        """Process picking in seismogram figure."""
        
        # ~ bg = self._fig.canvas.copy_from_bbox(self._fig.bbox)
        # ~ self._fig.canvas.blit(self._fig.bbox)
        
        if (self._fig.canvas.manager.toolbar.mode == "" and
            len(event.ind) == 1):
            ind = int(event.ind)
            me = event.mouseevent
            if (self._procmode == PROC_MODES.vel_est and
                me.button == 3):
                
                try: # workaround for occasionally occuring IndexError
                    segs = self._vellines.get_segments()
                    segs.pop(ind)
                    self._vellines.set_segments(segs)                
                    self._vellinetext[int(event.ind)].remove()
                    self._vellinetext.pop(int(event.ind))
                except IndexError: pass
            else:
                if self._pvd or self._cos:
                    return
                
                if (np.isclose(ind, me.xdata, atol=.1) and
                    not self._procmode == PROC_MODES.inactive):
                    if self._procmode == PROC_MODES.pick:
                        picks = self._picks[self._picks.pickset == \
                                            self._activeps]
                        idx = self._picks.index == picks.iloc[ind].name
                        
                        if me.button == 1:
                            self._picks.loc[idx, 'traveltime'] = me.ydata
                        elif me.button == 3:
                            self._picks.loc[idx, 'traveltime'] = -1
                            
                        self._draw_picks()
                        
                        # update traveltime plot
                        if self._ttfig != None:
                            self._save_picks()
                            self._plot_traveltimes()
                            
                    elif me.button == 1:
                        rin = self._picks.iloc[ind].receiver_index_number
                        
                        cmd = """SELECT polarity
                                 FROM receivers
                                 WHERE receiver_index_number==%d""" % (rin)
                        pol = int(self.slh.read_data(cmd).polarity)
                        
                        if self._procmode == PROC_MODES.trace_reverse:
                            if pol != 0: pol *= -1
                        elif self._procmode == PROC_MODES.trace_mute:
                            pol = 1 if pol==0 else 0
                        
                        cmd = """UPDATE receivers
                                 SET polarity = %d
                                 WHERE receiver_index_number==%d""" % (pol, 
                                                                       rin)
                        with self.slh.dbc:
                            self.slh.dbc.execute(cmd)

                        self._ax.cla()
                        self._draw_seismograms()
                        
            # ~ self._fig.canvas.blit(self._fig.bbox)
            # ~ self._fig.canvas.flush_events()
            self._fig.canvas.draw()

    def _onpress_seismograms(self, event):
        """Process mouse button pressed in seismogram figure."""
        
        if self._fig.canvas.manager.toolbar.mode == "":
            self._pressed = True
            if self._procmode == PROC_MODES.pick:
                self._pickline = [[event.xdata, event.xdata],
                                  [event.ydata, event.ydata]]
            elif (self._procmode == PROC_MODES.vel_est and
                  event.button == 1):
                self._velline = [[event.xdata, event.xdata],
                                 [event.ydata, event.ydata]]
                self._draw_velline()

    def _update_estimated_velocities(self):
        """Add the newly created velocity estimation line."""
        
        if self._vellines is None:
            if self._velline is not None:
                self._vellines = LineCollection(
                    [[[self._velline[0][0], self._velline[1][0]],
                      [self._velline[0][1], self._velline[1][1]]]],
                    colors='k', linewidths=2.,
                    linestyles='solid', capstyle='round')
                self._vellines.set_picker(True)
                
        else:
            segs = self._vellines.get_segments()
            segs.append([[self._velline[0][0],
                          self._velline[1][0]],
                         [self._velline[0][1],
                          self._velline[1][1]]])
            self._vellines.set_segments(segs)

    def _onrelease_seismograms(self, event):
        """Process mouse button release in seismogram figure."""
        
        if (self._fig.canvas.manager.toolbar.mode == "" and
            self._pressed):

            if (self._procmode == PROC_MODES.pick and
                self._picklineplot is not None):
                    
                spl = np.array([self._pickline[0][0],
                                self._pickline[1][0]])
                epl = np.array([self._pickline[0][1],
                                self._pickline[1][1]])
                                
                # check direction of picking line
                if spl[0] > epl[0]:
                    spl, epl = epl.copy(), spl.copy()
                
                for i in self._selin.index:
                    if i >= spl[0] and i <= epl[0]:
                        ssg = np.array([i, self._ylim[0]])
                        esg = np.array([i, self._ylim[1]])
                        
                        isec = compute_line_intersection(spl, epl,
                                                         ssg, esg)
                        
                        picks = self._picks[self._picks.pickset == \
                                            self._activeps]
                        idx = self._picks.index == picks.iloc[i].name
                        
                        if event.button == 1:
                            self._picks.loc[idx, 'traveltime'] = isec[1]
                        elif event.button == 3:
                            self._picks.loc[idx, 'traveltime'] = -1
                
                self._picklineplot.remove()
                self._draw_picks()
                self._picklineplot = None
                
                # update traveltime plot
                if self._ttfig != None:
                    self._save_picks()
                    self._plot_traveltimes()
                
            elif (self._procmode == PROC_MODES.vel_est and 
                  self._vellineplot is not None and
                  self._moved):
                self._vellineplot.remove()
                self._vellineplot = None
                self._update_estimated_velocities()
                self._draw_estimated_velocities()
                
            self._pressed = False
            self._moved = False

    def _onscroll_seismograms(self, event):
        """Process scrolling in the seismogram figure."""
        
        # ~ if self._amplitude:
        if self._scrollmode == SCROLL_MODES.amplitude:
            if event.button == 'up':
                self._scaling *= 2
            elif event.button == 'down':
                self._scaling /= 2
            
            self._ylim = event.inaxes.get_ylim()
            
            # ~ for line in self._psg_sg:
                # ~ print(line.get_ydata())
                # ~ line.set_ydata(line.get_ydata() * self._scaling)
            
            
            # ~ self._fig.canvas.restore_region(self._bg)
            # ~ self._fig.canvas.blit(self._fig.bbox)
            # ~ self._fig.canvas.flush_events()
            
            self._ax.cla()
            self._draw_seismograms()
        else:
            # get axis
            ax = event.inaxes
            
            # get current y-axis limits
            cur_ylim = ax.get_ylim()

            if event.button == 'up':
                self._ylim = (cur_ylim[1]
                              + (cur_ylim[0] - cur_ylim[1]) * (1 / 1.5),
                              cur_ylim[1])
            elif event.button == 'down':
                self._ylim = (cur_ylim[1]
                              + (cur_ylim[0]-cur_ylim[1]) * 1.5,
                              cur_ylim[1])
            ax.set_ylim(self._ylim)
            plt.draw()

    def _set_statusbartext(self, x, y):
        """Show relevant processing information in the status bar of the 
        plot window.

        """
        
        proc = ["Inactive",
                "Fb pick",
                "Trc mute",
                "Trc rev",
                "Vel est"][self._procmode]
        scroll = ["Zoom", "Amp scal"][self._scrollmode]
        time = 'Time (s) = %.3f' % (y)
        
        text = "Proc: %s | Scroll: %s | %s" % \
               (proc, scroll, time)
        
        return text

    def _plot_seismograms(self, **kwargs):
        """Create and setup figure for plotting the selected seismograms.
        
        Parameters
        ----------
        **kwargs: arbitrary keyword arguments
            Keyword arguments forwarded to the underlying plotting methods of
            the matplotlib, e.g., 'fontsize', 'color' 
        """
        
        self._check_plotting_ready()
                
        plt.rcParams['xtick.bottom'] = False
        plt.rcParams['xtick.labelbottom'] = False
        plt.rcParams['xtick.top'] = False
        plt.rcParams['xtick.labeltop'] = True
        
        self._fig, self._ax = plt.subplots(1, constrained_layout=True,
                                           figsize=(8, 6))
        # ~ self._procmode = PROC_MODES.pick
        self._draw_seismograms(**kwargs)
        
        # ~ if not self._pvd:
            # ~ self._draw_estimated_velocities()
            
        # ~ if not self._cos and not self._pvd:
            # ~ self._draw_picks(**kwargs)
        
        # register pick event handler
        self._fig.canvas.mpl_connect('pick_event',
                                       self._onpick_seismograms)
        
        # ~ # register onclick event handler
        self._fig.canvas.mpl_connect('button_release_event',
                                     self._onrelease_seismograms)
        
        # ~ # register onpress event handler
        self._fig.canvas.mpl_connect('button_press_event',
                                     self._onpress_seismograms)
        
        # register onmove event handler
        self._fig.canvas.mpl_connect('motion_notify_event',
                                     self._onmove_seismograms)
        
        # register onscroll event handler
        self._fig.canvas.mpl_connect('scroll_event',
                                     self._onscroll_seismograms)
        
        # register onclose event handler
        self._fig.canvas.mpl_connect('close_event',
                                     self._onclose_seismograms)
        
        # register keypress event handler
        self._fig.canvas.mpl_connect('key_release_event',
                                     self._onkeyrelease_seismograms)
        
        # ~ fig_manager = plt.get_current_fig_manager()
        #~ if hasattr(fig_manager, 'window'):
            #~ fig_manager.window.showMaximized()
        #~ self._fig.tight_layout(rect=[0, 0.03, 1, 0.95])
        
        self._ax.format_coord = self._set_statusbartext
        
        plt.show(block=False)
        
        # ~ plt.pause(0.1)
        
        # ~ self._bg = self._fig.canvas.copy_from_bbox(self._fig.bbox)
        # ~ if self._psg_sg != None: self._ax.draw_artist(self._psg_sg)
        # ~ if self._psg_fb != []: 
            # ~ for fb in self._psg_fb:
                # ~ self._ax.draw_artist(fb[0])
                # ~ self._ax.draw_artist(fb[1])
        # ~ if self._psg_isvd != None: self._ax.draw_artist(self._psg_isvd)
        
        # ~ self._fig.canvas.blit(self._fig.bbox)
        
        # ~ self._draw_seismograms(**kwargs)
        
        # ~ self._fig.canvas.blit(self._fig.bbox)

    def _plot_spectrum(self, **kwargs):
        """Plot the (stacked) frequency content of a currently selected traces.
        
        Parameters
        ----------
        **kwargs: arbitrary keyword arguments
            Keyword arguments forwarded to the underlying plotting methods of
            the matplotlib, e.g., 'fontsize', 'color'         
        """
        
        self._check_processing_ready()
        self._check_plotting_ready()
        
        c = kwargs.pop('c', kwargs.pop('color', 'g'))
        lw = kwargs.pop('lw', kwargs.pop('linewidth', 2.))
        a = kwargs.pop('a', kwargs.pop('alpha', .2))
        
        stacked_psd = np.array([])
        for tr in self._st:
            # compute fft
            tmp_fft = fftpack.fft(tr.data)
            
            # compute power spectral density
            tmp_psd = np.abs(tmp_fft)**2
            
            if len(stacked_psd) > 0: stacked_psd += tmp_psd
            else: stacked_psd = tmp_psd
            
        stacked_psd /= len(self._st)
        
        # get corresponding frequencies
        fftfreq = fftpack.fftfreq(len(stacked_psd),
                                  self._st[0].stats.delta)
        
        # omit negative frequencies
        idcs = fftfreq > 0
        
        # express in dB
        psd_db = 10 * np.log10(stacked_psd[idcs])
            
        # plot sprectrum
        self._spfig, self._spax = plt.subplots(1, figsize=(8, 4),
                                               constrained_layout=True)
                                               
        self._spax.plot(fftfreq[idcs], psd_db,
                 c=c,
                 lw=lw, **kwargs)
        self._spax.fill_between(fftfreq[idcs], psd_db, np.min(psd_db),
                         color=c,
                         alpha=a,
                         **kwargs)
        
        self._spax.set_xlim((np.min(fftfreq[idcs]),
                            np.max(fftfreq[idcs])))
        self._spax.set_ylim((np.min(psd_db),
                            np.max(psd_db)))
        
        self._spax.set_xlabel('Frequency (Hz)')
        self._spax.set_ylabel('PSD (dB)')
        
        self._spax.xaxis.set_ticks_position('bottom')
        
        self._spax.spines['top'].set_visible(False)
        self._spax.spines['right'].set_visible(False)
            
        plt.show(block=False)

    def _onclose_traveltimes(self,event):
        """Clean up after closing the traveltime figure."""
        
        self._ttfig = None
        self._ttax = None

    def _onpick_traveltimes(self, event):
        """Process picking in traveltime figure."""
        
        me = event.mouseevent
        rsid = round(me.xdata, 0)
        # ~ ind = event.ind
        sin = int(float(event.artist.properties()['label']))
        
        cmd = """SELECT first_geophone fg
                 FROM shots s LEFT JOIN geometry g
                 ON s.station_id == g.station_id
                 WHERE shot_index_number == %d""" % (sin)
        fg = int(self.slh.read_data(cmd).fg)
        
        cmd = """SELECT receiver_index_number rin
                 FROM receivers
                 WHERE station_id == %d""" % (rsid)
        rin = int(self.slh.read_data(cmd).rin)
        
        # ~ print('RSIN', rsid)
        # ~ print('RIN', rin)
        # ~ print('FG', fg)
        # ~ print('RIN', ind + fg)
        cmd = """UPDATE fbpicks
                 SET traveltime = -1
                 WHERE pickset == \'%s\' AND
                       shot_index_number == %d AND
                       receiver_index_number == %d AND
                       traveltime >= 0""" % (self._activeps, sin, rin)
                       # ~ traveltime >= 0""" % (self._activeps, sin, ind+fg)
        with self.slh.dbc:
            self.slh.dbc.execute(cmd)
        
        self._plot_traveltimes()
        if self._ax is not None:
            self._load_picks()
            self._draw_picks()

    def _plot_traveltimes(self, **kwargs):
        """Plot the traveltime curves for the currently active pickset.
        
        Parameters
        ----------
        **kwargs: arbitrary keyword arguments
            Keyword arguments forwarded to the underlying plotting methods of
            the matplotlib, e.g., 'fontsize', 'color' 
        """
        
        self._check_processing_ready()
        
        if self._ttfig == None:
            # create figure
            self._ttfig, self._ttax = plt.subplots(2, 1, figsize=(12, 9),
                                                   sharex=True, 
                                                   gridspec_kw={
                                                    'height_ratios': [.1, 1]})
        else:
            self._ttax[0].cla()
            self._ttax[1].cla()
        
        ls = kwargs.get('ls', kwargs.get('linestyle', '-'))
        lw = kwargs.get('lw', kwargs.get('linewidth', 1.5))
        
        cmd = """SELECT fbp.shot_index_number sin,
                        fbp.receiver_index_number rin,
                        fbp.traveltime tt,
                        gs.x sx,
                        gr.x rx,
                        gs.station_id ssid,
                        gr.station_id rsid
                 FROM fbpicks fbp
                 LEFT JOIN shots s
                     ON fbp.shot_index_number == s.shot_index_number
                 LEFT JOIN receivers r
                     ON fbp.receiver_index_number == r.receiver_index_number
                 LEFT JOIN geometry gs
                     ON s.station_id == gs.station_id
                 LEFT JOIN geometry gr
                     ON r.station_id == gr.station_id
                 WHERE fbp.pickset == \'%s\' AND
                       tt >= 0""" % (self._activeps)
        fbpicks = self.slh.read_data(cmd)
        
        # define colormap
        cmap = plt.cm.tab10(np.arange(10))

        # initialize station array
        # ~ stations = np.array([])
        cmd = """SELECT station_id
                 FROM geometry"""
        stations = self.slh.read_data(cmd)

        # plot receiver stations
        if len(fbpicks) > 0: ymax = fbpicks.tt.max()
        else: ymax = 0.
        
        cmd = """SELECT s.shot_index_number sin, 
                        g.x sx, 
                        g.station_id sid
                 FROM shots s LEFT JOIN geometry g
                 ON s.station_id == g.station_id"""
        shts = self.slh.read_data(cmd)
        
        cmd = """SELECT r.receiver_index_number sin, 
                        g.x rx,
                        g.station_id sid
                 FROM receivers r LEFT JOIN geometry g
                 ON r.station_id == g.station_id"""
        recs = self.slh.read_data(cmd)
        
        for i, r in recs.iterrows():
            # ~ if not r[1] in shts.iloc[:, 1].values:
            if not r[2] in shts.iloc[:, 2].values:
                # ~ self._ttax[0].plot(r[1], 0,
                self._ttax[0].plot(r[2], 0,
                                   marker='v', c='k',
                                   markersize=8,
                                   markeredgecolor='k')
                # ~ stations = np.append(stations, r[1])
                # ~ stations = np.append(stations, r[2])
                
            self._ttax[0].annotate(str(int(r[0])),
                                  # ~ (r[1], 0),
                                  (r[2], 0),
                                  xytext=(0, -20),
                                  textcoords='offset pixels',
                                  ha='center',
                                  fontsize=self._fontsize - 2)
        
        X, Y, colors = [], [], []
                    
        # plot traveltime diagram
        for i, s in shts.iterrows():
            # draw station symbol
            # ~ if s[1] in recs.iloc[:, 1].values:
            if s[2] in recs.iloc[:, 2].values:
                marker = 'o'
                # ~ stations = np.append(stations, s[1])
                # ~ stations = np.append(stations, s[2])
            else:
                marker = '*'
                # ~ stations = np.append(stations, s[2])
                
            # ~ self._ttax[0].plot(s[1], 0,
            self._ttax[0].plot(s[2], 0,
                               marker=marker,
                               c=cmap[int(s[0]) % 10],
                               markersize=10,
                               markeredgecolor='k')
            
            self._ttax[0].annotate(str(int(s[0])),
                                   # ~ (s[1], 0),
                                   (s[2], 0),
                                   xytext=(0, 10), textcoords='offset pixels',
                                   ha='center',
                                   fontsize=self._fontsize - 2)
            
            # draw traveltime curve
            tt = fbpicks[fbpicks.sin == s[0]].tt
            # ~ rx = fbpicks[fbpicks.sin == s[0]].rx
            rsid = fbpicks[fbpicks.sin == s[0]].rsid
            
            if len(tt) > 0:
                # ~ X.append(list(rx))
                X.append(list(rsid))
                Y.append(list(tt))
                colors.append(cmap[int(s[0]) % 10])
            
                # ~ self._ttax[1].scatter(rx, tt,
                self._ttax[1].scatter(rsid, tt,
                                      color=cmap[int(s[0]) % 10],
                                      marker=self._marker,
                                      picker=True,
                                      label=str(s[0]),
                                      **kwargs)
        
        lines = LineCollection((list(zip(x, y)) for x, y in zip(X, Y)),
                               colors=colors, linestyles='-')
        self._ttax[1].add_collection(lines)
        
        for stat in stations.to_numpy():
            self._ttax[1].axvline(stat,
                                  c='lightgray', lw=.5, zorder=0)
            self._ttax[0].plot([stat, stat],
                               [0, -.1],
                               c='lightgray',
                               lw=.5,
                               zorder=0)
        
        # format plot
        # ~ xmin = np.min([shts.sx.min(), recs.rx.min()])
        # ~ xmax = np.max([shts.sx.max(), recs.rx.max()])
        xmin = np.min([shts.sid.min(), recs.sid.min()])
        xmax = np.max([shts.sid.max(), recs.sid.max()])
        dx = xmax - xmin
             
        self._ttax[0].set_xlim((xmin - dx*.01, xmax + dx*.01))
        self._ttax[0].set_ylim((-.1, .1))
        self._ttax[0].set_yticks([-.0425, .035])
        self._ttax[0].set_yticklabels(['RIN', 'SIN'],
                                      fontsize=self._fontsize-2)
        self._ttax[0].yaxis.set_ticks_position('none')
        self._ttax[0].spines['top'].set_visible(False)
        self._ttax[0].spines['bottom'].set_visible(False)
        self._ttax[0].spines['left'].set_visible(False)
        self._ttax[0].spines['right'].set_visible(False)

        self._ttax[1].set_ylabel('Traveltime (s)')
        self._ttax[1].set_xticks([])
        self._ttax[1].spines['top'].set_visible(False)
        self._ttax[1].spines['bottom'].set_visible(False)
        self._ttax[1].spines['right'].set_visible(False)
        self._ttax[1].grid()
        self._ttax[1].invert_yaxis()

        self._ttfig.set_tight_layout(True)
        self._ttfig.subplots_adjust(hspace=0)

        self._ttfig.canvas.draw()
        
        # register pick event handler
        self._ttfig.canvas.mpl_connect('pick_event',
            self._onpick_traveltimes)
        
        # register onclose event handler
        self._ttfig.canvas.mpl_connect('close_event',
            self._onclose_traveltimes)
        
        plt.show(block=False)
        
        return self._ttax

    def _preview_data(self):
        """Plot the seismograms without geometry information directly from
        the raw data files.
        """
        
        # set preview flag
        self._logger.info('Starting in data preview mode')
        # ~ self._pvd = True
        self._procmode = PROC_MODES.inactive
        
        # read all raw data files
        self._read_data()
        
        # select first shot file
        self._pvk = min(list(self._data.keys()))
        self._selected = self._fnames[self._pvk]
        self._xlabel = 'RSN'
        self._st = self._data[self._pvk]
        self._pst = self._st.copy()
        
        self._plot_seismograms()

    def _parse_plot_params(self, options):
        """Check the basic validity of the plot parameters and return them
        as numpy array.
        
        Parameters
        ----------
        options: str, default ''
            Options for plotting.
        
        Returns
        -------
        params: numpy array
            Parameter deciding which plot method to use.
        """
        
        if not isinstance(options, str):
            self._looger.error('Parameter of datatype str expected')
            return 0
            
        if options == "":
            return ""
            
        option = options.lower().rstrip()
        params = options.split(" ")
        
        if len(params) > 1:
            self._logger.critical('More than one plot option provided')
            return 0
        
        if params[0] not in PLOT_KEYWORDS:
            self._logger.error('Invalid keyword given')
            return 0
        
        return params[0]

    def plot(self, options='', **kwargs):
        """Plot different visualizations of the data.
        
        Parameters
        ----------
        options: str
            Parameters defining which visualization to plot. An empty string
            ('') will plot the currently selected traces. Other valid keywords
            are 'traveltimes', 'pseudosection', 'pickperc', 'spectrum'.
            
        **kwargs: arbitrary keyword arguments
            Keyword arguments forwarded to the underlying plotting methods of
            the matplotlib, e.g., 'fontsize', 'color' etc.
        """
        
        self._logger.input('plot ' + options)

        param = self._parse_plot_params(options)
        if param == 0: return
        
        # set matplotlib keywords
        self._fontsize = kwargs.pop('fontsize', 10)
        self._marker = kwargs.pop('marker', 'x')
        
        font = {'family': kwargs.pop('fontfamily', 'sans-serif'),
                'weight': kwargs.pop('fontweight', 'normal'),
                'size': self._fontsize}
        mpl.rc('font', **font)
                
        if param == "pickperc":
            return self._plot_pickpercentage(**kwargs)
        elif param == "traveltimes":
            return self._plot_traveltimes(**kwargs)
        elif param == "pseudosection":
            return self._plot_pseudosection(**kwargs)
        elif not self._check_plotting_ready():
            return
                
        if param == "":
            return self._plot_seismograms(**kwargs)
        elif param == "spectrum":
            return self._plot_spectrum(**kwargs)

    def _select_aoffset(self, aoffset):
        """Select traces based on the absolute offset.
        
        Parameters
        ----------
        aoffset: float
            The absolute offset in m.
        """
        
        cmd = """SELECT shot_index_number sin,
                        receiver_index_number rin
                 FROM applied_geometry
                 WHERE absolute_offset==%.3f""" % (aoffset)
        self._selin = self.slh.read_data(cmd)
        
        for sin, rin in zip(self._selin.sin, self._selin.rin):
            # iterate all shots
            for k, st in self._data.items():
                if st[0].stats.SIN == sin:
                    for tr in st:
                        if tr.stats.RIN == rin:
                            self._st += tr.copy()
            
        self._st.sort(['SIN'])
        self._xlabel = 'SIN'

    def _select_midpoint(self, midpoint):
        """Select traces based on the midpoint.
        
        Parameters
        ----------
        midpoint: float
            The midpoint in m.
        """
        
        cmd = """SELECT shot_index_number sin,
                        receiver_index_number rin
                 FROM applied_geometry
                 WHERE midpoint_x==%.3f""" % (midpoint)
        self._selin = self.slh.read_data(cmd)
        
        for sin, rin in zip(self._selin.sin, self._selin.rin):
            # iterate all shots
            for k, st in self._data.items():
                if st[0].stats.SIN == sin:
                    for tr in st:
                        if tr.stats.RIN == rin:
                            self._st += tr.copy()
                            
        self._st.sort(['SIN'])
        self._xlabel = 'SIN'

    def _select_sin(self, sin):
        """Select traces based on the shot index number.
        
        Parameters
        ----------
        sin: float
            The shot index number.
        """

        cmd = """SELECT shot_index_number sin,
                        receiver_index_number rin
                 FROM applied_geometry
                 WHERE shot_index_number==%d""" % (sin)
        self._selin = self.slh.read_data(cmd)

        if self._selin.empty:
            self._logger.error('SIN %d not found' % (sin))
            return
        
        for k, st in self._data.items():
            if len(self._st) != 0:
                break
            
            if st[0].stats.SIN == sin:
                for tr in st:
                    tr.stats.channel = '%03d' % (tr.stats.RIN)
                    self._st += tr.copy()
                            
        self._st.sort(['channel'])
        self._xlabel = 'RIN'
            
    def _select_rin(self, rin):
        """Select traces based on the receiver index number.
        
        Parameters
        ----------
        rin: float
            The receiver index number.
        """

        cmd = """SELECT shot_index_number sin,
                        receiver_index_number rin
                 FROM applied_geometry
                 WHERE receiver_index_number==%d""" % (rin)
        self._selin = self.slh.read_data(cmd)

        if self._selin.empty:
            self._logger.error('RIN %d not found' % (rin))
            return
        
        for k, st in self._data.items():
            for tr in st:
                if tr.stats.RIN == rin:
                    tr.stats.channel = '%03d' % (tr.stats.SIN)
                    self._st += tr.copy()
                    
        self._st.sort(['channel'])
        self._xlabel = 'SIN'

    def _parse_select_condition(self, condition):
        """Check the basic validity of the select condition and return 
        the parameters as numpy array.
        
        Parameters
        ----------
        condition: str
            Condition for trace selection.
            
        Returns
        -------
        parameters: tuple (str, float)
            Keyword and corresponding number.
        """
        
        # check the data type
        if isinstance(condition, str):
            condition = condition.lower().rstrip()
            params = condition.split(" ")
        else:
            self._logger.error('Parameter of datatype str expected')
            return 0
            
        # check the number of arguments
        if len(params) != 2:
            self._logger.error('Insufficient number of parameters')
            return 0
        
        # check if valid keyword is provided
        if params[0] not in SELECT_KEYWORDS:
            self._logger.error('Invalid keyword given')
            return 0
        
        # check if second parameter is a numerical value
        try:
            float(params[1])
        except ValueError:
            self._logger.error('Second parameter needs to be a ' \
                'numerical value')
            return 0
        
        return params

    def select(self, condition, auto=False):
        """Select traces based and make the data accessible for 
        filtering and picking.
        
        Parameters
        ----------
        condition: str
            Condition for trace selection. The string has to contain
            'sin', 'rin' or 'aoffset' followed by a numerical value.
            'sin', 'rin', 'aoffset' refer to the shot index number
            (SIN), the receiver index number (RIN) and the absolute
            offset, respectively.
            
        auto: bool
            True if select method was invoked automatically.
            
        Examples
        --------
        >>> from useis import SeismicRefractionManager
        >>> srm = SeismicRefractionManager(<project_directory>)
        >>> srm.select('sin 3')
            
        """
        
        if auto:
            self._logger.auto('select ' + condition)
        else:
            self._logger.input('select ' + condition)
        if not self._check_processing_ready():
            return
        params = self._parse_select_condition(condition)
        if params == 0: return
        
        self._actssns = np.array([])
        self._actrsns = np.array([])
        self._st = Stream()

        if params[0] == 'aoffset':
            self._select_aoffset(float(params[1]))

        elif params[0] == 'sin':
            self._select_sin(int(params[1]))
            
        elif params[0] == 'rin':
            self._select_rin(int(params[1]))
            
        elif params[0] == 'midpoint':
            self._select_midpoint(float(params[1]))
        
        if auto: 
            self._logger.autoinfo('%d traces selected' % len(self._st))
        else:
            self._logger.info('%d traces selected' % len(self._st))
        
        # set processing and plotting related attributes
        self._title = ''
        self._mode = ''
        self._cos = False
        self._selected = ' '.join(params).upper()
        if self._procmode == PROC_MODES.inactive:
            self._procmode = PROC_MODES.pick
        
        # load picks
        self._load_picks()
        
        # apply/remove filter and plot the data (if autoplot is on)
        self._pst = self._st.copy()
        if len(self._st) > 0:
            if self._filterhold and self._filtered != '':
                self.filter(self._filtered.lower(), auto=True)
            else:
                self.filter('remove', auto=True)
