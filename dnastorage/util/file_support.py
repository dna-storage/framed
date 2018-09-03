#!/usr/bin/python
import os
import sys

"""
Write a file by receiving (index,value) packets of data of a specific size at a given index. Value
is of packetSize length.  index specifies where in the file to write the data.  Packets can be
received in any order. But, the file cannot be written until the file is complete, in other words,
it has received all size/packetSize packets.  The indices are assumed to be densely numbered from 0
to size/packetSize by 1.  Pythonically: range(0,size/packetSize+(size-size/packetSize*packetSize),1).
"""
class WritePacketizedFilestream:
    def __init__(self,fd,size,packetSize):
        self.__fd = fd
        self.size = size
        self.__set_packetSize(packetSize)
        self.__data = {}
        
    def has_key(self,key):
        return self.__data.has_key(key)
    def __setitem__(self,key,value):
        if key >= 0 and key < self.numberOfPackets:
            self.__data[key] = value
    def __getitem__(self,key):
        assert key >= 0 and key < self.numberOfPackets
        return self.__data[key]

    @property
    def numberOfPackets(self):
        if (self.size % self.packetSize) > 0:
            return self.size / self.packetSize + 1
        else:
            return self.size / self.packetSize
        
    # packet property
    def __set_packetSize(self,val):
        self.__packetSize = val
    def __get_packetSize(self):
        return self.__packetSize
    packetSize = property(__get_packetSize,__set_packetSize)

    @property
    def lastPacketSize(self):
        # this function should never return 0, not the same as modulo
        # should return a number >= 1 and <= self.packetSize
        return self.size - (self.numberOfPackets-1)*self.packetSize

    @property
    def complete(self):
        keys = self.__data.keys()
        keys.sort()
        if len(keys) < self.numberOfPackets:
            return False
        for i in range(self.numberOfPackets):
            if keys[i] != i:
                return False
        return True

    def getMissingKeys(self):
        keys = self.__data.keys()
        keys.sort()
        missing = []
        i = 0
        for k in keys:
            if i != k:
                while i < k:
                    missing.append(i)
                    i+=1
            i+=1
        return missing
        

    ## Warning: requires buffering the whole file!
    def write(self):
        items = self.__data.items()
        items.sort(cmp=lambda x,y: cmp(x[0],y[0]))
        i = 0
        emptyString = '\x00'*self.packetSize
        for key,value in items:
            if i < key:
                while i < key:                    
                    self.__fd.write(emptyString)
                    i+=1
            if i == self.numberOfPackets-1:
                self.__fd.write(value[0:self.lastPacketSize])
            else:
                self.__fd.write(value)
            i+=1
        if self.__fd != sys.stdout:
            self.__fd.close()
    
    def dummy_write(self):
        items = self.__data.items()
        items.sort(cmp=lambda x,y: cmp(x[0],y[0]))
        i = 0
        emptyString = '\x00'*self.packetSize
        output_data=""
        for key,value in items:
            if value is not type(str):
                value2="".join(value)
            if i < key:
                while i < key:                    
                    output_data+=emptyString
                    i+=1
            if i == self.numberOfPackets-1:
                output_data+=value2[0:self.lastPacketSize]
            else:
                output_data+=value2
            i+=1
        return output_data


            
class WritePacketizedFile(WritePacketizedFilestream):
    def __init__(self,filename,size,packetSize):
        WritePacketizedFilestream.__init__(self,open(filename,"wb"),size,packetSize)
        self.__filename = filename



"""
Read a file by breaking it up into packetSize pieces. If the last piece is smaller
than packetSize, it's padded with zeros. Clients are guaranteed all packets are uniform
size.

This is an iterable object. Packets can be read using the iterator or by requesting
a specific index between 0 and size/packetSize. 

This class expects a file descriptor. There is a derived class that accepts a filename.
"""
class ReadPacketizedFilestream:    
    def __init__(self,fd):
        self.__fd = fd
        self.__set_packetSize(120)
        self.__read_size = 0

    def read(self):
        b = self.__fd.read(self.packetSize)
        self.__read_size += len(b)
        if b and len(b) != self.packetSize:
            b = b.ljust(self.packetSize,'\x00')           
        return b

    # packet property
    def __set_packetSize(self,val):
        self.__packetSize = val
    def __get_packetSize(self):
        return self.__packetSize
    packetSize = property(__get_packetSize,__set_packetSize)

    # iterator
    def __iter__(self):
        # reset to begnning of file
        self.__fd.seek(0)
        return self

    def next(self):
        b = self.read()
        if b:
            return b
        else:
            raise StopIteration()        
    
    def __getitem__(self,key):
        self.__fd.seek(key*self.packetSize)
        return self.read()        
       
    @property
    def size(self):
        return os.fstat(self.__fd.fileno()).st_size
    @property
    def bytes_read(self):
        return self.__read_size
        
    @property
    def numberOfPackets(self):
        if (self.size % self.packetSize)==0:
            return self.size / self.packetSize
        else:
            return self.size / self.packetSize + 1
            

class ReadPacketizedFile(ReadPacketizedFilestream):    
    def __init__(self,filename):
        ReadPacketizedFilestream.__init__(self,open(filename,"rb"))
        self.__filename = filename
    @property
    def filename(self):
        return self.__filename


if __name__ == "__main__":
    import os
    import sys
    from random import randint
    packetFile = ReadPacketizedFile(sys.argv[1])
    
    out = WritePacketizedFile("output.d",packetFile.size,120)

    assert out.complete==False

    i = 0
    for p in packetFile:
        out[i] = p
        i += 1
    
    assert out.complete==True
    out.write()
