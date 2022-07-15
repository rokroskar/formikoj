from formikoj.utils import *

class MethodManager():
    'Common base class for all Manager classes'
    
    def __init__(self, workdir):
        """Create an instance of the SeismicRefraction class.
        
        Parameters
        ----------
        workdir: str
            Path of the project directory.
        """
        
        if os.path.exists(workdir):
            self._workdir = workdir
            self._init_base_attributes()
            return 1
        else:
            print('\033[31mDIRECTORY NOT FOUND\033[0;0m')
            return 0

    def _init_base_attributes(self):
        """Create and initialize base attributes relevant for any 
        kind of Manager.
        """
        
        # path to database file
        self._dbfile = os.path.join(self._workdir, 'prj.db')
        
        # set directories
        self._datadir = os.path.join(self._workdir, '01_data/raw')
        self._geomdir = os.path.join(self._workdir, '02_geom')
        self._procdir = os.path.join(self._workdir, '03_proc')

    def _check_directory_structure(self):
        """Check if the directory structure is valid."""
        
        if not (os.path.exists(self._datadir) and \
                os.path.exists(self._geomdir) and \
                os.path.exists(self._procdir)):
            return 0
        
        return 1
