from obspy import UTCDateTime, Trace, Stream
from pygimli.physics.seismics import *
import yaml
from yaml.loader import SafeLoader

from . import DataModeler
from ..utils import *

class SeismicWaveformModeler(DataModeler):
    'Class allowing to generate seismic waveform data'
    
    def __init__(self, workdir):
        """Create an instance of the SeismicWaveformModeler class.
        
        Parameters
        ----------
        workdir: str
            Path of the project directory.
        """
        
        if not DataModeler.__init__(self, workdir):
            return

        self._logger = create_logger('useis.SeismicWaveformModeler',
                                     os.path.join(self._workdir, 
                                                  self._logfile))
        
        self._init_base_attributes()
        self._logger.info('Created instance of SeismicWaveformModeller')
    
    def _init_base_attributes(self):
        """Create and initialize the fundamental attributes."""
        
        self._cfg = None
        self._wavelet = None
        self.mesh = None
        self.velmod = None
        self._scheme = None
        
        self._skipshots = []
        self._brokengeo = []
        self._wrongpol = []       
    
    def _check_modeling_ready(self):
        """Check if the object is ready for data processing."""
        
        if self._cfg is None:
            self._logger.error('Configuration not loaded')
            return 0
        
        if self.mesh is None:
            self._logger.error('Mesh not loaded')
            return 0
            
        if self.velmod is None:
            self._logger.error('Velocity model not available')
            return 0
            
        if self._scheme is None:
            self._logger.error('Measurement scheme not loaded')
            return 0
            
        if self._wavelet is None:
            self.logger.error('Wavelet not generated')
            return 0
        
        return 1
    
    def _export_station_coordinates(self):
        """Write the station coordinates to a csv file."""
        
        # get station coordinates
        stats = np.array(self._scheme.sensorPositions())
        stats[:, 1], stats[:, 2] = stats[:, 2].copy(), stats[:, 1].copy()
        stats = np.hstack(((np.arange(stats.shape[0])+1).reshape(-1, 1),
                            stats))
        
        np.savetxt(os.path.join(self._dsdir, 'data/station_coords.csv'),
           stats,
           delimiter=',',
           fmt='%d,%.1f,%.1f,%.1f',
           header='station,x (m), y (m), z (m)')
    
    def _write_protocol(self, protocol, shtgeo):
        """Write the protocol for the given dataset.
        
        Parameters
        ----------
        protocol: numpy array
            Contains shot number and correspondig shot indices.
            
        shtgeo: numpy array
            Shot and geophone indices in first and second column, respectively.
        """
        
        dataset = self._dsdir.split('/')[-1]
            
        with open(os.path.join(self._dsdir, 'data/protocol.txt'), 'w+') as f:
            f.write(' #############################\n')
            f.write('   Line: SYN_%s             \n' % (dataset))
            f.write('   Sampling rate: %d Hz    \n' % 
                (self._cfg['wavelet']['sampling_rate']))
            f.write('   Recording length: %.3f s \n' % 
                (self._cfg['wavelet']['length']))
            f.write('   Number of geophones: %d   \n' % (self._ngeo))
            f.write('   Geophone spacing: %d m     \n' % (self._geospc))
            f.write(' #############################\n')
            
            f.write('\n     File number | Station\n')
            
            # get unique shot and geophone indices
            s = np.unique(np.array(self._scheme['s']))
            g = np.unique(np.array(self._scheme['g']))
            
            # iterate protocol entries
            for p in protocol:
                # check if shot was conducted at geophone position
                if p[1] in g:
                    station = 'G%03d' % (p[1] + 1)
                else: # shot was conducted between geophone positions
                    # determine indices of closest geophones                    
                    didx = np.abs(g - p[1])
                    clogeo = g[didx==didx.min()]
                    
                    station = 'G%03d | G%03d' % (clogeo[0], clogeo[1])
                
                f.write('         %04d    |   %s\n' % 
                        (p[0], station))
                        
        with open(os.path.join(self._dsdir, 'info.txt'), 'w+') as f:
            f.write('Number of geophones: %d\n' % (self._ngeo))
            f.write('Number of shots: %d\n' % (self._nsht))
            
            f.write('Recording length (s): %.3f\n' % 
                (self._cfg['wavelet']['length']))
            f.write('Sampling frequency (Hz): %d\n' % 
                (self._cfg['wavelet']['sampling_rate']))
            f.write('Wavelet type: %s\n' % ('Ricker'))
            f.write('Frequency of the wavelet (Hz): %d\n\n' % 
                (self._cfg['wavelet']['frequency']))
            
            f.write('Missing shot(s): ' + 
                     ', '.join([str(np.where(
                                        np.unique(shtgeo[:, 0])==e)[0][0]+1) 
                                    for e in self._skipshots]) + '\n')
            f.write('Broken geophone(s): ' + 
                     ', '.join([str(np.where(
                                        np.unique(shtgeo[:, 1])==e)[0][0]+1)
                                    for e in self._brokengeo]) + '\n')
            f.write('Wrong polarity geophone(s): ' + 
                     ', '.join([str(np.where(
                                        np.unique(shtgeo[:, 1])==e)[0][0]+1)
                                    for e in self._wrongpol]) + '\n')

    def _set_dataset_problems(self, shtgeo):
        """Determine settings and values regarding missing shots, broken 
           geophones, wrong polarity.
           
        Parameters
        ----------
        shtgeo: numpy array
            Shot and geophone indices in first and second column, respectively.
        """
        
        # determine maximum number of affected shots
        naffsht = round(self._nsht * .05)
        
        # determine shots to skip
        if self._cfg['dataset']['missing_shots']:
            nskipshots = int(np.random.randint(1, naffsht + 1, 1))
            self._skipshots = np.random.choice(shtgeo[:, 0],
                                size=nskipshots, replace=False)

        # determine maximum number of affected traces
        nafftrc = round(self._ngeo * .05)
                                
        # determine broken geophones
        if self._cfg['dataset']['broken_geophones']:
            nbrokengeo = int(np.random.randint(1, nafftrc + 1, 1))
            self._brokengeo = np.random.choice(shtgeo[:, 1],
                                size=nbrokengeo, replace=False)
                                
        # determine traces with wrong polarity for non-broken geophones
        if self._cfg['dataset']['wrong_polarity']:
            nwrongpol = int(np.random.randint(1, nafftrc + 1, 1))
            self._wrongpol = np.random.choice(
                np.setdiff1d(shtgeo[:, 1], self._brokengeo),
                size=nwrongpol, replace=False)

    def _create_data(self, sht, geos):
        """Create seismic waveform data for a single shot position.
        
        Parameter
        ---------
        sht: int
            Index of the shot.
            
        geos: numpy array
            Indices of the geophones.
            
        Returns
        -------
        st: obspy Stream
            Contains the synthetic waveform data for the given shot
            position.
        """
        
        stats = self._scheme.sensorPositions()
                
        u = np.array(
            solvePressureWave(self.mesh, self.velmod,
                              times=self._wavelet[:, 0],
                              sourcePos=stats[sht],
                              uSource=self._wavelet[:, 1],
                              verbose=True))
        
        st = Stream()
        for gi in geos:
            if gi not in self._brokengeo:
                data = u[:, self.mesh.findNearestNode(stats[gi])].copy()

                # add offset dependent noise
                if self._cfg['dataset']['noise']:
                    w = float(pg.utils.dist(
                        p=[stats[gi]], c=stats[sht]))
                    data += (self._cfg['dataset']['noise_level'] * \
                        np.random.randn(self._wavelet.shape[0]) * w)
                
                if gi in self._wrongpol:
                    data *= -1
            else:
                data = np.ones_like(self._wavelet[:, 0]) * \
                    np.random.randn(self._wavelet.shape[0])
            
            data = np.require(data, dtype=np.float32)
            
            header = {}
            header = {'station': str(sht + 1001), 
                      'channel': str(gi + 1),
                      'npts': len(data),
                      'sampling_rate':
                          self._cfg['wavelet']['sampling_rate'],
                      'mseed': {'dataquality': 'D'}}
            header['starttime'] = UTCDateTime('1986-07-14T00:00:00')
            
            st += Trace(data=data, header=header)
        
        return st
    
    def _check_dataset_parameters(self):
        """Check dataset parameters provided in the config file."""
        
        if 'dataset' not in self._cfg.keys():
            self._logger.error('No dataset parameters provided')
            return 0
            
        if 'names' in self._cfg['dataset'].keys():
            if not isinstance(self._cfg['dataset']['names'], list):
                self._logger.error("""Dataset name(s) need to be provided as 
                                      list""")
                return 0
                
            self._cfg['dataset']['number'] = \
                len(self._cfg['dataset']['names'])
            return 1
                
        elif 'number' in self._cfg['dataset'].keys():
            if not (isinstance(self._cfg['dataset']['number'], int) or
                    isinstance(self._cfg['dataset']['number'], float)):
                self._logger.error("Number of datasets needs to be " +
                                   "numeric")
                return 0
                        
            if self._cfg['dataset']['number'] <= 0:
                self._logger.error('Number of datasets needs to be positive')
                return 0
                
            self._set_dataset_names(
                ['dataset_%003d' % (i + 1) 
                 for i in np.arange(self._cfg['dataset']['number'])])
            return 1
        else:
            self._logger.error('Neither dataset names nor number provided')
            return 0
    
    def _create_waveforms(self):
        """Create seismic waveform dataset(s)."""
        
        if self._cfg is None:
            if not self._load_config():
                return
        
        if self.mesh is None:
            if not self._load_mesh():
                self._logger.info('Trying to create a mesh...')
                if not self._create_mesh():
                    return
            
        if self.velmod is None:
            if not self._create_velocity_model():
                return
            
        if self._scheme is None:
            if not self._load_scheme():
                return
            
        if self._wavelet is None:
            if not self._create_wavelet():
                return
        
        if not self._check_dataset_parameters():
            return
        
        sg = np.vstack((np.array(self._scheme['s']),
                        np.array(self._scheme['g']))).T
        
        for ds in self._cfg['dataset']['names']:
            self._dsdir = os.path.join(self._outdir, ds)
            if not os.path.exists(self._dsdir): os.makedirs(self._dsdir)
            if not os.path.exists(os.path.join(self._dsdir, 'data')):
                os.makedirs(os.path.join(self._dsdir, 'data'))

            self._set_dataset_problems(sg)
            
            prot = np.empty((0, 2))
            for i, si in enumerate(np.unique(sg[:, 0])):
                if si in self._skipshots:
                    continue
                
                st = self._create_data(si, sg[sg[:, 0] == si, 1])
            
                st.write(os.path.join(self._dsdir, 'data/Shot_%04d.syn' %
                         (i + 1001)), format='MSEED')
                prot = np.vstack((prot, np.array([i + 1001, si])))
                
            self._export_station_coordinates()
            self._write_protocol(prot, sg)
            
            self._logger.info('Dataset %s created' % (ds))

    def _create_velocity_model(self):
        """Create the velocity model based on the mesh and the velocity map
           provided in the config file."""

        if self._cfg is None:
            if not self._load_config():
                return
           
        if 'model' not in self._cfg.keys():
            self._logger.error('No model parameters provided in config file')
            return 0
        
        if 'velmap' not in self._cfg['model'].keys():
            self._logger.error('No velocity map provided in config file')
            return 0
            
        self.velmod = pg.solver.parseMapToCellArray(
            self._cfg['model']['velmap'], 
            self.mesh)
            
        self._logger.info('Velocity model created')
        return 1

    def _check_wavelet_parameters(self):
        """Check wavelet parameters provided in the config file."""
        
        if self._cfg is None:
            if not self._load_config():
                return 0
        
        if 'wavelet' not in self._cfg.keys():
            self._logger.error('No wavelet parameters provided')
            return 0
            
        if np.sum(np.in1d(['length', 'frequency', 'sampling_rate'],
                          list(self._cfg['wavelet'].keys()))) != 3:
            self._logger.error("Wavelet length, frequency and " + 
                               "sampling_rate need to be provided")
            return 0
            
        return 1
    
    def _create_wavelet(self):
        """Create the wavelet used for the modeling of the seismic waveforms.
           At this point only ricker wavelets are supported.
        """
        
        if not self._check_wavelet_parameters():
            return 0
        
        t = np.arange(0, self._cfg['wavelet']['length'], 
                      1. / self._cfg['wavelet']['sampling_rate'])
        
        # determine pretrigger
        pretrigger = self._cfg['wavelet'].get(
            'pretrigger', self._cfg['wavelet']['length']*.02)
        
        wl = ricker(self._cfg['wavelet']['frequency'], 
                    t, 
                    1./self._cfg['wavelet']['frequency'] + pretrigger)

        self._wavelet = np.stack((t, wl)).T
        
        self._logger.info('Wavelet created')
        
        return 1
    
    def _check_traveltime_parameters(self):
        """Check traveltime parameters provided in the config file."""
        
        if 'traveltimes' not in self._cfg.keys():
            self._logger.error('No traveltime parameters provided')
            return 0
            
        if 'names' in self._cfg['dataset'].keys():
            if not isinstance(self._cfg['dataset']['names'], list):
                self._logger.error("""Dataset name(s) need to be provided as 
                                      list""")
                return 0
                
            self._cfg['dataset']['number'] = \
                len(self._cfg['dataset']['names'])
            return 1
                
        elif 'number' in self._cfg['dataset'].keys():
            if not (isinstance(self._cfg['dataset']['number'], int) or
                    isinstance(self._cfg['dataset']['number'], float)):
                self._logger.error("Number of datasets needs to be " +
                                   "numeric")
                return 0
                        
            if self._cfg['dataset']['number'] <= 0:
                self._logger.error('Number of datasets needs to be positive')
                return 0
                
            self._set_dataset_names(
                ['dataset_%003d' % (i + 1) 
                 for i in np.arange(self._cfg['dataset']['number'])])
            return 1
        else:
            self._logger.error('Neither dataset names nor number provided')
            return 0
    
    def _create_traveltimes(self):
        """Create the first break traveltimes for a given model."""
        
        if self._cfg is None:
            if not self._load_config():
                return
        
        if self.mesh is None:
            if not self._load_mesh():
                self._logger.info('Trying to create a mesh...')
                if not self._create_mesh():
                    return
            
        if self.velmod is None:
            if not self._create_velocity_model():
                return
            
        if self._scheme is None:
            if not self._load_scheme():
                return
        
        for ds in self._cfg['dataset']['names']:
            self._dsdir = os.path.join(self._outdir, ds)
            if not os.path.exists(self._dsdir): os.makedirs(self._dsdir)
        
            # compute first break traveltimes
            ttmgr = TravelTimeManager(verbose=False)
            
            noise_rel, noise_abs = 0., 0.
            if 'traveltimes' in self._cfg.keys():
                noise_rel = self._cfg['traveltimes'].get('noise_relative', 0.)
                noise_abs = self._cfg['traveltimes'].get('noise_absolute', 0.)
            
            ttdata = ttmgr.simulate(
                self.mesh, 
                scheme=self._scheme, 
                slowness=1. / self.velmod,
                noiseLevel=noise_rel, noiseAbs=noise_abs,
                secNodes=self._cfg['model'].get('sec_nodes', 3),
                verbose=False)
            
            # add pre-trigger
            if self._check_wavelet_parameters():
                pretrigger = self._cfg['wavelet'].get(
                    'pretrigger', self._cfg['wavelet']['length'] * .02)
                t = np.array(ttdata['t']) + pretrigger
                ttdata.set('t', t)
                
            ttdata.set('valid', np.ones(ttdata.size()))
            
            ttdata.save(os.path.join(self._dsdir, '%s_tt.pck' % (ds)))
    
    def _parse_create_params(self, options):
        """Check the basic validity of the create parameters and return them
        as numpy array.
        
        Parameters
        ----------
        options: str, default ''
            Options for different create activities.
        
        Returns
        -------
        params: numpy array
            Parameters deciding which create method to use.
        """
        
        if not isinstance(options, str):
            self._logger.error('Parameter of datatype str expected')
            return 0
            
        if options == "":
            self._logger.critical('No compute command provided')
            return 0
        
        option = options.lower().rstrip()
        params = options.split(" ")    
        
        if params[0] not in CREATE_KEYWORDS:
            self._logger.error('Invalid keyword given')
            return 0
        
        return params
    
    def create(self, options):
        """Create/model different kind of objects/data.
        
        Parameters
        ----------
        options: str
            The string has to contain a valid keyword ('waveforms', 'mesh',
            'velmod', 'wavelet') to call the respective methods.
        """
        
        self._logger.input('create ' + options)
        
        params = self._parse_create_params(options)
        if params == 0: return
        
        if params[0] == 'waveforms':
            self._create_waveforms()
        elif params[0] == 'mesh':
            self._create_mesh()
        elif params[0] == 'velmod':
            self._create_velocity_model()
        elif params[0] == 'wavelet':
            self._create_wavelet()
        elif params[0] == 'traveltimes':
            self._create_traveltimes()
    
    def _load_config(self, filename=None):
        """Read yaml file containing the configuration for the forward 
        modeling.
        
        Parameters
        ----------
        filename: str, default None
           Name of the config file to load
        """
           
        if filename is None:           
            files = glob(os.path.join(self._inpdir, '*.yml'))
            
            if not len(files):
                self._logger.error('No configuration file found')
                return 0
                
            if len(files) > 1:
                self._logger.warning('Multiple configuration files found')
                self._logger.info('Provide specific configuration file name')
                return 0
                
            cfgfile = files[0]
        else:
            cfgfile = os.path.join(self._inpdir, filename)
            if not os.path.exists(cfgfile):
                self._logger.error('Config file %s not found' % (cfgfile))
                return 0
            
        with open(cfgfile, "r") as f:
            self._cfg = yaml.load(f, Loader=SafeLoader)
            
        self._logger.info('Configuration loaded')
        return 1
        
    def _load_mesh(self, filename=None):
        """Read the mesh to be used for the forward modeling.
        
        Parameters
        ----------
        filename: str, default None
           Name of the mesh file to load
        """
        
        if filename is None:
            files = glob(os.path.join(self._inpdir, '*.bms'))
            
            if not len(files):
                self._logger.warning('No mesh file found')
                return 0
            if len(files) > 1:
                self._logger.warning('Multiple mesh files found')
                self._logger.info('Provide specific mesh file name')
                return 0
                
            meshfile = files[0]
        else:
            meshfile = os.path.join(self._inpdir, filename)
            if not os.path.exists(meshfile):
                self._logger.error('Mesh file %s not found' % (meshfile))
                return 0
                
            if not meshfile.endswith('.bms'):
                self._logger.error('Mesh with unsupported file format')
                return 0
            
        self.mesh = pg.load(meshfile)
        self._logger.info('Mesh loaded')
        return 1 
    
    def _create_scheme_datacontainer(self, data):
        """Create a pygimli DataContainer from scheme csv file.
        
        Parameters
        ----------
        inp: numpy array
            Columns 1-3 contain geometry of the stations, columns 4 and
            5 defines whether this station is a geophone or shot 
            position.
            
        """
        
        self._scheme = pg.DataContainer()
        self._scheme.registerSensorIndex('s')
        self._scheme.registerSensorIndex('g')
        
        for r in data[['x', 'z', 'y']].to_numpy():
            self._scheme.createSensor(r)
            
        S, G = [], []
        for s in data.index[data['s']==1].to_list():
            for g in data.index[data['g']==1].to_list():
                S.append(s)
                G.append(g)
                
        self._scheme.resize(len(S))
        self._scheme.set('s', S)
        self._scheme.set('g', G)
        self._scheme.set('valid', np.abs(np.sign(
            self._scheme['g'] - self._scheme['s'])))
    
    def _load_scheme(self, filename=None):
        """Read the measurement scheme.
        
        Parameters
        ----------
        filename: str, default None
           Name of the measurement scheme file to load
        """
        
        if not self._check_directory_structure(): return 0
        
        if filename is None:
            files = glob(os.path.join(self._inpdir, 'scheme.*'))
            
            if not len(files):
                self._logger.error('No measurement scheme found')
                return
                
            if len(files) > 1:
                self._logger.error('Multiple measurement schemes found')
                return
                
            schemefile = files[0]
        else:
            schemefile = os.path.join(self._inpdir, filename)
            if not os.path.exists(os.path.join(self._inpdir, schemefile)):
                self._logger.error('scheme file %s not found' % (schemefile))
                return 0
        
        if schemefile.lower().endswith('.dat'):
            self._logger.info('Measurement scheme in unified data format ' \
                              'found')
            self._scheme = pg.load(schemefile)
            
        elif schemefile.lower().endswith('.csv'):
            inp = pd.read_csv(schemefile, 
                              sep=None, engine='python',
                              names=['x', 'y', 'z', 'g', 's'])
            
            self._create_scheme_datacontainer(inp)
        else:
            self._logger.error('Measurement scheme with unsupported ' \
                               'file format')
        
        # determine number of geophones
        self._ngeo = len(np.unique(np.array(self._scheme['g'])))
        
        # determine number of shots
        self._nsht = len(np.unique(np.array(self._scheme['s'])))
        
        # determine geophone spacing
        stats = np.array(self._scheme.sensorPositions())
        g = np.array(self._scheme['g'])
        self._geospc = np.median(np.diff(np.sqrt(
                (stats[g[0], 0] - stats[g, 0])**2 + 
                (stats[g[0], 1] - stats[g, 1])**2 +
                (stats[g[0], 2] - stats[g, 2])**2)))
        
        self._logger.info('Measurement scheme loaded')
        return 1
        
    def _parse_load_params(self, options):
        """Check the basic validity of the load parameters and return them
        as numpy array.
        
        Parameters
        ----------
        options: str, default ''
            Options for different load activities.
        
        Returns
        -------
        params: numpy array
            Parameters deciding which load method to use.
        """
        
        if not isinstance(options, str):
            self._logger.error('Parameter of datatype str expected')
            return 0
            
        if options == "":
            self._logger.critical('No load command provided')
            return 0
        
        option = options.lower().rstrip()
        params = options.split(" ")    
        
        if params[0] not in LOAD_KEYWORDS:
            self._logger.error('Invalid keyword given')
            return 0

        return params

    def load(self, options):
        """Load different objects, data, resources, information etc.
        
        Parameters
        ----------
        options: str
            The string has to contain a valid keyword ('config', 'mesh',
            'scheme') to call the respective methods.
        """
        
        self._logger.input('set ' + options)
        
        params = self._parse_load_params(options)
        if not params: return
        
        if params[0] == 'config':
            filename = None if len(params) != 2 else params[1]
            self._load_config(filename=filename)
        elif params[0] == 'mesh':
            filename = None if len(params) != 2 else params[1]
            self._load_mesh(filename=filename)
        elif params[0] == 'scheme':
            self._load_scheme()
