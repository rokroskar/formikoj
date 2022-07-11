import numpy as np
import pandas as pd
import sqlite3 as sl

class SQLiteHandler():
    'Class for handling communication with and manipulation of SQLite database'
    
    def __init__(self, database):
        """Create an instance of the SQLiteHandler class.
        
        Parameters
        ----------
        database: string
            Name of the database
        """
    
        self.database = database
        self.open_connection(self.database)

    def open_connection(self, database):
        """Open the connection to the database.
        
        Parameters
        ----------
        database: string
            Name of the database
        """
        self.dbc = sl.connect(database, timeout=10)
        with self.dbc:
            self.dbc.execute('PRAGMA foreign_keys = ON')

    def close_connection(self):
        """Close the connection to the database."""
        
        self.dbc.close()

    def create_table(self, table, cols, dtypes, pks, fks=[]):
        """Create a table in the database from a pandas dataframe.
        
        Parameters
        ----------
        table: pandas DataFrame
            Contains relevant information to create the table.
        """
        
        # table
        cmd = ["CREATE TABLE %s (" % (table)]
        
        # columns
        for c, d in zip(cols, dtypes):
            tmp = '%s %s, ' if c not in pks else '%s %s NOT NULL, '
            cmd.append(tmp % (c, d))
        
        # foreign key(s)
        for fk in fks:
            cmd.append('FOREIGN KEY (%s) REFERENCES %s (%s), ' % 
                       (fk[0], fk[1], fk[2]))
        
        # primary key(s)
        cmd.append('PRIMARY KEY (')
        for pk in pks[:-1]:
            cmd.append('%s, ' % (pk))
        cmd.append('%s)' % (pks[-1]))
        
        cmd.append(');')
        
        with self.dbc:
            self.dbc.execute("".join(cmd))

    def write_data(self, table, data, mode='append'):
        """Write data to database table.
        
        Parameters
        ----------
        table: string
            Table to write into
            
        data: pandas DataFrame
            Data to write into database
        """
        
        data.to_sql(table, self.dbc, if_exists=mode, index=False)
        
    def read_data(self, cmd):
        """Read data from the database.
        
        Returns
        -------
        df: pandas DataFrame
            Data read from the database
        """
        
        df = pd.read_sql(cmd, self.dbc)
        return df
