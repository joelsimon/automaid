# -*- coding: utf-8 -*-
#
# Part of automaid -- a Python package to process MERMAID files
# pymaid environment (Python v2.7)
#
# Author: Joel. D Simon
# Contact: jdsimon@alumni.princeton.edu | joeldsimon@gmail.com
# Last modified: 15-Dec-2020, Python 2.7.15, Darwin-18.7.0-x86_64-i386-64bit

import csv
import pickle

class GeoCSV:
    """Organize MERMAID metadata and write them to GeoCSV format

    Args:
        dives (list): List of dives.Dive instances (def: None)
        version (str): GeoCSV version (def: '2.0')
        delimiter (str): GeoCSV delimiter (def: ',')

    Attributes:
        dives (list): List of dives.Dive instances

    """

    def __init__(self, dives=None, version='2.0', delimiter=','):
        self.dives = dives
        self.version = version
        self.delimiter = ','

        # Define first two header lines
        self.dataset_header = '#dataset: GeoCSV ' + self.version
        self.delimiter_header = '#delimiter: ' + self.delimiter

        
    def write(self, filename='test.geocsv'):
        """Write GeoCSV file 
        
        Args:
            filename (str): GeoCSV filename

        """
        #filename.append('geocsv') if filename.split('.')[-1] != 'geocsv' else pass

        # Open as as 'wb' in Python 2 rather than 'w' with newline='' in Python 3
        # https://docs.python.org/2/library/csv.html#csv.writer
        with open(filename, 'wb') as csvfile:
            csvwriter = csv.writer(csvfile, delimiter=self.delimiter)
            csvwriter.writerow([self.dataset_header])
            csvwriter.writerow([self.delimiter_header])
                        

if __name__ == "__main__":
    dive_list = pickle.load( open('mdives.p', 'rb'))
    geocsv = GeoCSV(dive_list)