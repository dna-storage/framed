class BaseCodec:
    def __init__(self,CodecObj=None):
        self._Obj = CodecObj
    def _encode(self, s):
        return s
    def encode(self,s):
        if self._Obj != None:
            return self._encode(self._Obj.encode(s))
        else:
            return self._encode(s)
    def _decode(self, s):
        return s
    def decode(self,s):
        s = self._decode(s)
        if self._Obj != None:
            return self._Obj.decode(s)
        else:
            return s
