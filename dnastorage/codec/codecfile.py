#!/usr/bin/python
from dnastorage.util.file_support import *
from dnastorage.codec import base_conversion
from dnastorage.codec import dense

class EncodePacketizedFile:
    def __init__(self,packetizedFile,CodecObj=None):
        self._packetizedFile = packetizedFile
        self._iterating = False
        if CodecObj == None:
            self._Codec = dense.DenseCodec()
        else:
            self._Codec = CodecObj

    @property
    def bytes_encoded(self):
        return self._packetizedFile.bytes_read
         
    # Derived classes should override only this method, unless other
    # functionality also needs to be altered
    # returns key, value
    def _encode(self):
        return self._packetizedFile.next()

    # Ideally, this would not be overriden, but if a subclass wants to
    # fully control all steps of encoding, it must be overriden
    def encode(self):
        return self._Codec.encode(self._encode())

    # Sub-classes that alter key names, may need to override this
    # Right now, only using one file descriptor, so this will
    # cause iteration to break. Either use this interface or iteration,
    # not both!
    def __getitem__(self,key):
        if self._iterating:
            assert False and "Warning: should not use __getitem__ while iterating!"
        return self._packetizedFile[key]

    # iterator
    def __iter__(self):
        self._iterating = True
        self._packetizedFile.__iter__()
        # reset to begnning of file
        return self

    def next(self):
        packet = self.encode()
        if packet:
            return packet
        else:
            self._iterating = False
            raise StopIteration()        

class DecodePacketizedFile:
    def __init__(self,packetizedFile,CodecObj=None):
        self._packetizedFile = packetizedFile
        self._data = {}
        if CodecObj == None:
            self._Codec = dense.DenseCodec()
        else:
            self._Codec = CodecObj
        
    def has_key(self,key):
        return self._data.has_key(key)

    #def __setitem__(self,key,value):
    #    self._data[key] = value
    #def __getitem__(self,key):
    #    return self._data[key]
        
    def writeToFile(self,key,value):
        self._packetizedFile[key] = value

    # Ideally, derived classes will only override this implementation
    def _decode(self,key,value):
        self.writeToFile(key,value)

    def decode(self,strand):
        key,val = self._Codec.decode(strand)
        self._decode(key,val)

    @property
    def complete(self):
        return self._packetizedFile.complete
    
    def write(self):
        self._packetizedFile.write()

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
