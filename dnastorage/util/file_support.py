#!/usr/bin/python
import os

class WritePacketizedFile:
    def __init__(self,filename,size,packetSize):
        self.__filename = filename
        self.size = size
        self.__set_packetSize(packetSize)
        self.__data = {}
        
    def has_key(self,key):
        return self.__data.has_key(key)
    def __setitem__(self,key,value):
        assert key >= 0 and key < self.numberOfPackets
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
        return self.size - (self.numberOfPackets)*self.packetSize

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
        f = open(self.__filename,"wb")
        items = self.__data.items()
        items.sort(cmp=lambda x,y: cmp(x[0],y[0]))
        i = 0
        emptyString = '\x00'*self.packetSize
        for key,value in items:
            if i < key:
                while i < key:                    
                    f.write(emptyString)
                    i+=1
            if i == self.numberOfPackets-1:
                f.write(value[0:self.lastPacketSize])
            else:
                f.write(value)
            i+=1
        f.close()
        

class ReadPacketizedFile:    
    def __init__(self,filename):
        self.__filename = filename
        self.__fd = open(self.__filename,"rb")
        self.__set_packetSize(120)

    @property
    def filename(self):
        return self.__filename

    def read(self):
        b = self.__fd.read(self.packetSize)
        print(self.packetSize)
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
    def numberOfPackets(self):
        if (self.size % self.packetSize)==0:
            return self.size / self.packetSize
        else:
            return self.size / self.packetSize + 1
            
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
