# Class to read an IDL database.
# Based on the IDL implementation from
# http://idlastro.gsfc.nasa.gov/ftp/pro/database/
#
# It reads a header file (.dbh), extracts information
# about the fields and reads all records from the
# database file (.dbf). Currently only the data types
# that I encountered are supported, only one database
# is supported (no links to others).
#
# Laurens Keek, Sept. 2007

import struct
import numpy as n
from numpy import rec
from astropy.table import Table
import logging

class IDLDatabase:
    """
    Read the contents of an IDL database.
    """
    
    def __init__(self, filename):
        """
        Start reading an IDLdatabase with given
        filename. Note that the filename should
        not include a suffix (eg. ".dbf").
        """
        # Read header file
        f = open(filename+'.dbh', 'rb')
        self.read_hdf_header(f)
        self.read_hdf_items(f)
        f.close()
        
        self.create_format()
        
        # Read database file
        f = open(filename+'.dbf', 'rb')
        self.read_dbf_header(f)
        self.read_dbf_items(f)
        f.close()
    
    def read_hdf_header(self, file):
        """
        Read the header of the header file
        """
        fmt = '>19s61s2HI2HI2BHIBI9x2B'
        result = struct.unpack(fmt, file.read(struct.calcsize(fmt)))
        self.header = result
        
        if result[-2]==1: # New format (4 byte record size)
            self.record_size = result[-3]
        else: # Old format (2 byte record size)
            self.record_size = result[3]
    
    def read_hdf_items(self, file):
        """
        Read the items following the header in the
        header file.
        """
        items = []
        for i in range(self.header[2]):
            fmt = '>20s4Hb69sHb19s6s45s4H21x'
            result = struct.unpack(fmt, file.read(200))
            items.append(result)
        
        self.header_items = items
        self.field_names = [i[0].strip().lower() for i in items]
        if not isinstance(self.field_names[0], str):
            self.field_names = [i.decode('ascii') for i in self.field_names] # Python 3 gives string as 'bytes': must cast to string for records

        labels = [i[6].strip(b' \x00').decode('ascii') for i in items]
        
        self.field_labels = dict(zip(self.field_names, labels))
        formats = [self.idl2py_format(i[10].strip()) for i in items]
        self.field_format = dict(zip(self.field_names, formats))
    
    def read_dbf_header(self, file):
        """
        Read the header of the database file.
        """
        fmt = '>2I'
        result = struct.unpack(fmt, file.read(struct.calcsize(fmt)))
        
        self.num_records = result[0]
    
    def read_dbf_items(self, file):
        """
        Read all items from the database file.
        """
        fmt = self.record_format
        size = self.record_size
        
        file.seek(size)
        # 'big' below corresponds to 'big-endian', character '>'
        records = rec.array(file, formats=fmt, shape=self.num_records,
                          names=self.field_names, byteorder='big')

        # For consistency with the text input, we convert to an astropy Table here
        # This step is primarily to allow the strings to be treated as strings (not
        # byte arrays, as they are for the base rec.array)

        self.records = Table(records)

        # Also strip out the trailing spaces in the various strings
        # This is designed for the MINBAR databases, for which these columns are
        # present; but we check here in case we're using for other databases

        for column in ['name','obsid','notes']:
            if column in self.field_names:
                self.records[column] = [x.strip() for x in self.records[column]]


    def create_format(self):
        """
        From the specification given in the header file,
        generate the format to unpack the records in the
        database file.
        """
        # Map IDL data types to convention of struct module.
        datatype = {2:'i2', 3:'u4', 4:'f4', 5:'f8', 7:'a'}
        fmt = ''
        total = 0
        
        for i in self.header_items:
            dtype = datatype[i[1]]
            if dtype == 'a':
                dtype = 'a%i'%i[4]
            if i[2]>1:
                dtype = '%i%s' % (i[2], dtype)
            fmt += dtype
            fmt += ','
            total += int(dtype[1:])
        
        # See if we need some padding at the end
        pad = (self.record_size - total)
        if pad<0:
            logging.warn('idldatabase:Sum of items is larger than record size?')
        elif pad>0:
            fmt += 'a%i'%pad
            fmt += ','
        
        self.record_format = fmt[:-1]
        self.num_fields = len(fmt.split(','))
    
    def idl2py_format(self, format):
        """
        Convert an IDL string format string to python's
        equivalent. Assumes only one 'field' is formatted
        and that each 'f' field contains a '.'.
        """
        if format[0] == 'a':
            return '%%%ss'%format[1:]
        elif format[0] == 'i':
            return '%%%si'%format[1:]
        elif format[0] == 'f':
            part = format[1:].split('.')
            return '%%%s.%sf'%(part[0], part[1])
