import pygimli.meshtools as mt
from ..utils import *

class DataModeler():
    'Class allowing to generate synthetic geophysical data'
    
    def __init__(self, workdir):
        """Create an instance of the DataModeler class.
        
        Parameters
        ----------
        workdir: str
            Path of the project directory.
        """
        
        if os.path.exists(workdir):
            self._workdir = workdir
            self._inpdir = os.path.join(self._workdir, 'in')
            self._logfile = 'history.log'
            
            self._cfg = None
            return 1
        else:
            print('\033[31mDIRECTORY NOT FOUND\033[0;0m')
            return 0
    
    def _check_directory_structure(self):
        """Check if the directory structure is valid."""
        
        if not os.path.exists(self._inpdir):
            self._logger.error('Directory structure invalid')
            return 0
    
        self._outdir = os.path.join(self._workdir, 'out')
        if not os.path.exists(self._outdir):
            os.makedirs(self._outdir)
            
        return 1
            
    def _create_mesh(self):
        """Create forward modeling mesh."""
        
        if self._cfg is None:
            if not self._load_config():
                return 0
                
        if self._scheme is None:
            if not self._load_scheme():
                return 0
        
        if 'model' not in self._cfg.keys():
            self._logger.error('No model parameters provided')
            return 0
        
        if 'layers' not in self._cfg['model'].keys():
            self._logger.error('No layer information provided')
            return 0
        
        topo = np.array(self._scheme.sensorPositions())[:, [0, 2]]
        
        dx = topo[:, 0].max() - topo[:, 0].min()
        
        topo = np.vstack((
            np.array([topo[0, 0] - 1*dx, topo[0, 1]]),
                np.vstack((topo,
                           np.array([topo[-1, 0] + 1*dx, topo[-1, 1]])))))
        
        model = []
        model.append(topo)
        
        geom = pg.Mesh(dim=2, isGeometry=True)
        for i, l in enumerate(self._cfg['model']['layers']):
            model.append(model[i].copy())

            thk = l[1]
            model[i+1][:, 1] -= thk
            poly = mt.createPolygon(np.vstack((model[i], model[i+1][::-1])),
                                    isClosed=True, marker=l[0],
                                    boundaryMarker=pg.core.MARKER_BOUND_MIXED)
            geom = mt.mergePLC([geom, poly])
            
        mesh = mt.createMesh(geom,
            quality=self._cfg['model'].get('quality', 32), 
            area=self._cfg['model'].get('area', 2), 
            smooth=self._cfg['model'].get('smooth', [1, 10]))
        
        self.mesh = mt.appendTriangleBoundary(mesh, marker=2,
                                              xbound=100*dx, ybound=100*dx,
                                              isSubSurface=True)

        self._logger.info('Created a 2D mesh')
        return 1
