from math import log,floor

class TableCodecBuilder:
    def __init__(self, filename, name, codewords):
        # None
        self._filename = filename
        self._name = name
        self._codewords = codewords
        self._cwLen = len(codewords[0])
        self._args = ", ".join([filename,name])
        self._nbytes = log(len(codewords),2)
        assert( self._nbytes - floor(self._nbytes) == 0)
        return

    def _classdecl(self,f):
        f.write("from dnastorage.codec.base import TableCodec\n\n")
        f.write("class "+self._name+"(TableCodec):                    \n")
        f.write("   def __init__(self,CodecObj=None,keyWidth={}):     \n".format(self._cwLen*2))

        args = "keyWidth,keyWidth/{},{},{}".format(self._cwLen,self._cwLen,int(log(len(self._codewords),2)/8))

        f.write("      TableCodec.__init__(self,CodecObj," + args + ")\n")
        f.write("      self._keyWidth = keyWidth\n");
        f.write("      self._etab = [");
        for i,c in enumerate(self._codewords):
            f.write("'{}'".format(c))
            if i<len(self._codewords)-1:
                f.write(',')
                if (i+1)%8==0:
                    f.write('\n                    ')
        f.write(']\n\n\n')
        f.write("      self._dtab = {");
        for i,c in enumerate(self._codewords):
            tabs = ""
            if i>0:
                tabs = "                    "
            f.write(tabs + "'{}' : {}".format(c,i))
            if i<len(self._codewords)-1:
                f.write(',\n')
        f.write('}\n\n\n')

    def _build_enctab(self,f):
        s =     "   def _enctab(self, val):    \n"
        s = s + "       return self._etab[val]  \n"
        s = s + " \n"
        f.write(s)

    def _build_dectab(self,f):
        s =     "   def _dectab(self, s):     \n"
        s = s + "       return self._dtab[s]  \n"
        s = s + " \n"
        f.write(s)

    def build(self):
        f = open(self._filename,"w")
        f.write("# Auto generated file by TableCodecBuilder     \n")
        f.write("# TableCodecBuilder("+self._args+", codewords) \n")
        f.write("#\n")
        self._classdecl(f)
        self._build_enctab(f)
        self._build_dectab(f)
        f.close()

if __name__ == "__main__":
    b = TableCodecBuilder("test.py","TestCodec",["AA","GG"])
    b.build()
