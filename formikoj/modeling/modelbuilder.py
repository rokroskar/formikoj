
class ModelBuilder():
    'Class allowing to build your own models'
    
    def __init__(self, workdir):
        """Create an instance of the ModelBuilder class.
        
        Parameters
        ----------
        workdir: str
            Path of the project directory.
        """
        
        if os.path.exists(workdir):
            self._workdir = workdir
