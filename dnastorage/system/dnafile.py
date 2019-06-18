from dnastorage.system.formats import *
from dnastorage.system.header import *
import io
import sys

class DNAFile:
    def __init__(self):
        return
    
    @classmethod
    def open(self, filename, op, primer5, primer3, format_name=""):
        # check if we are reading or writing
        if op=="r":            
            # 1. filename is the input file of strands
            # 2. primer_info optionally tells us the primers we're looking,
            #    if not specified, just assume first first 20 bases are
            #    are the primer.
            # 3. format is optional for reading and often ignored. Instead we look
            #    at the header to deduce what to do.
            #
            #b = io.BytesIO()
            #pf = WritePacketizedFilestream(b)
            #s = get_strands(name,primer_info)
            #h = decode_file_header(s)            
            #f = self(name,op,primer_info,h)
            return ReadDNAFile(input=filename, primer5=primer5, primer3=primer3)
        elif op=="w":
            return WriteDNAFile(output=filename,primer5=primer5,primer3=primer3, format_name=format_name)
        else:
            return None

    def flush(self):
        return

    def close(self):
        return

    def readable(self):
        assert False
    def writable(self):
        assert False

def get_strands(in_fd):
    strands = []
    while True:
        s = in_fd.readline()        
        if len(s) == 0:
            break
        s = s.strip()
        if s.startswith('%'):
            continue
        strands.append(s)
    return strands
    


class ReadDNAFile(DNAFile):
    # ReadDNAFile reads a set of strands from a file.  It finds the header,
    # determines compatibility and encoding type, and then decodes the file.
    #
    def __init__(self,**kwargs):     
        DNAFile.__init__(self)

        if kwargs.has_key('input'):
            self.input_filename = kwargs['input']
            self.in_fd = open(self.input_filename,"r")
        elif kwargs.has_key('in_fd'):
            self.in_fd = kwargs['in_fd']
            self.input_filename = ""

        assert kwargs.has_key('primer5') and kwargs.has_key('primer3')
        self.primer5 = kwargs['primer5']
        self.primer3 = kwargs['primer3']

        # get all the strands ( this will be bad for large files )
        strands = get_strands(self.in_fd)
        h = decode_file_header(strands,self.primer5,self.primer3)
        self.strands = pick_nonheader_strands(strands,self.primer5)
        
        assert h['version'][0] <= system_version['major']
        assert h['version'][1] <= system_version['minor']

        self.formatid = h['formatid']
        self.header = h 
        self.size = h['size']
        
        dec_func = file_system_decoder(self.formatid)
        
        self.mem_buffer = BytesIO()
        self.pf = WritePacketizedFilestream(self.mem_buffer,self.size,file_system_format_packetsize(self.formatid))
        self.dec = dec_func(self.pf,kwargs['primer5'],kwargs['primer3'])

        for s in self.strands:
            if s.startswith(self.primer5):
                self.dec.decode(s)

        self.dec.write()
        assert self.dec.complete

        self.mem_buffer.seek(0,0) # set read point at beginning of buffer
        return

    def read(self, n=1):        
        return self.mem_buffer.read(n)
    def readline(self, n=-1):
        return self.mem_buffer.readline(n)

    def readable(self):
        return True
    def writable(self):
        return False
        
class WriteDNAFile(DNAFile):
    # WriteDNAFile writes a set of strands.
    def __init__(self,**kwargs):     
        DNAFile.__init__(self)
        if kwargs.has_key('formatid'):            
            enc_func = file_system_encoder(kwargs['formatid'])
            self.formatid = kwargs['formatid']
        elif kwargs.has_key('format_name'):
            enc_func = file_system_encoder_by_abbrev(kwargs['format_name'])
            self.formatid = file_system_formatid_by_abbrev(kwargs['format_name'])

        self.mem_buffer = BytesIO()
        self.pf = ReadPacketizedFilestream(self.mem_buffer)
        self.enc = enc_func(self.pf,kwargs['primer5'],kwargs['primer3'])

        assert kwargs.has_key('primer5') and kwargs.has_key('primer3')
        self.primer5 = kwargs['primer5']
        self.primer3 = kwargs['primer3']

        if kwargs.has_key('output'):
            self.output_filename = kwargs['output']
            self.out_fd = open(self.output_filename,"w")
        elif kwargs.has_key('out_fd'):
            self.out_fd = kwargs['out_fd']
            self.output_filename = ""
            
        self.size = 0
        self.strands = []
        return

    def _encode_buffer(self):
        for e in self.enc:
            self.strands.append(e)

    def read(self, size):
        assert False and "Do not read at the same time as writing"
        
    def writable(self):
        return True
    def readable(self):
        return False
    
    def write(self, buff):
        self.size += len(buff)
        tell = self.mem_buffer.tell()
        self.mem_buffer.seek(0,2)
        self.mem_buffer.write(buff)
        self.mem_buffer.seek(tell,0)
        return
    
    def flush(self):
        self._encode_buffer()
        return

    def close(self):
        self.flush()
        hdr = encode_file_header(self.output_filename,self.formatid,self.size,\
                                 "",self.primer5,self.primer3)
        for i,h in enumerate(hdr):
            self.strands.insert(i,h)

        comment = encode_file_header_comments(self.output_filename,self.formatid,self.size,"",\
                                              self.primer5,self.primer3)
        self.out_fd.write(comment)
        for s in self.strands:
            self.out_fd.write("{}\n".format(s))

        if self.out_fd != sys.stdout and self.out_fd != sys.stderr:
            self.out_fd.close()
        return


class SegmentedWriteDNAFile(WriteDNAFile):
    # SegmentedWriteDNAFile writes a set of strands.
    def __init__(self,**kwargs):     
        WriteDNAFile.__init__(self,**kwargs)
        self.segments = []
        self.beginIndex = 0
        return

    def _record_segment(self):
        self.segments += [[ self.formatid, self.size, self.primer5, self.primer3, self.beginIndex ]]

    def new_segment(self, format_name, primer5, primer3):
        self._encode_buffer()  # write everything in the buffer to the file
        self._record_segment() # remember new segment
        self.primer5 = primer5
        self.primer3 = primer3
        # get current index + 1
        self.beginIndex = self.enc.index
        #print "beginIndex={}".format(self.enc.index)

        
        enc_func = file_system_encoder_by_abbrev(format_name)
        self.formatid = file_system_formatid_by_abbrev(format_name)
        # we consumed the prior buffer, so just make a new one to avoid
        # cursor positioning problems (JMT: not sure if this is the best way)
        self.mem_buffer = BytesIO()
        self.pf = ReadPacketizedFilestream(self.mem_buffer)
        self.enc = enc_func(self.pf,primer5,primer3,self.beginIndex)
        self.size=0

    def encode_segments_header(self,segments):
        oprimer5 = segments[0][2]
        oprimer3 = segments[0][3]
        assert len(segments) <= 256
        if len(segments[1:]) == 0:
            return []
        hdr = [ len(segments[1:]) ]
        for s in segments[1:]:
            hdr += convertIntToBytes(s[0],2)
            hdr += encode_size_and_value( s[1] )
            hdr += encode_size_and_value( s[4] )
            hdr += encode_primer_diff(oprimer5,s[2])
            hdr += encode_primer_diff(oprimer3,s[3])

        return hdr

    def close(self):
        self.flush()
        self._record_segment() # record last segment
        
        hdr_other = self.encode_segments_header(self.segments)

        formatid = self.segments[0][0]
        size = self.segments[0][1]
        primer5 = self.segments[0][2]
        primer3 = self.segments[0][3]
        
        hdr = encode_file_header(self.output_filename,formatid,size,hdr_other,primer5,primer3)
        for i,h in enumerate(hdr):
            self.strands.insert(i,h)

        comment = encode_file_header_comments(self.output_filename,formatid,size,hdr_other,primer5,primer3)
        self.out_fd.write(comment)
        for s in self.strands:
            self.out_fd.write("{}\n".format(s))

        if self.out_fd != sys.stdout and self.out_fd != sys.stderr:
            self.out_fd.close()
        return


class SegmentedReadDNAFile(ReadDNAFile):

    def decode_segments_header(self,other_data):
        val = other_data
        numSeg = val[0]
        if numSeg==0:
            return
        pos = 1
        allSegs = []
        for i in range(numSeg):
            # get format of this segment
            seg = [convertBytesToInt(val[pos:pos+2])]
            pos+=2

            # get size in bytes
            v,p = decode_size_and_value(val,pos)
            pos += p
            seg += [v]

            # get begin index
            v,p = decode_size_and_value(val,pos)
            pos += p
            seg += [v]

            primer5,p = decode_primer_diff(val[pos:], self.primer5)
            pos += p
            primer3,p = decode_primer_diff(val[pos:], self.primer3)
            pos += p
            seg.append(primer5)
            seg.append(primer3)
            allSegs.append(seg)            
        return allSegs



    # ReadDNAFile reads a set of strands from a file.  It finds the header,
    # determines compatibility and encoding type, and then decodes the file.
    #    
    def __init__(self,**kwargs):     
        ReadDNAFile.__init__(self,**kwargs)

        if len(self.header['other_data'])==0:
            return

        # restore cursor to end of buffer for writing
        self.mem_buffer.seek(0,2)
        
        segs = self.decode_segments_header(self.header['other_data'])
        self.segments = segs

        #print "segments=",segs

        for s in segs:
            formatid = s[0]
            size = s[1]
            bindex = s[2]
            primer5 = s[3]
            primer3 = s[4]
            if kwargs.has_key('use_single_primer'):
                # strands from sequencing should all have the same primer
                primer5 = self.primer5
                primer3 = self.primer3

            dec_func = file_system_decoder(formatid)
            #self.mem_buffer = BytesIO()
            self.pf = WritePacketizedFilestream(self.mem_buffer,size,\
                                                file_system_format_packetsize(formatid),\
                                                minKey=bindex)
            self.dec = dec_func(self.pf,primer5,primer3,bindex)

            for s in self.strands:
                if s.startswith(primer5):
                    #print s
                    self.dec.decode(s)

            write_anyway = kwargs.has_key('write_incomplete_file') and \
                kwargs['write_incomplete_file']==True
            if self.dec.complete or write_anyway:
                self.dec.write()            
            #assert self.dec.complete
            
        #print [x for x in self.mem_buffer.getvalue()]
        self.mem_buffer.seek(0,0) # set read point at beginning of buffer        
        return

    def read(self, n=1):        
        return self.mem_buffer.read(n)
    def readline(self, n=-1):
        return self.mem_buffer.readline(n)

    def readable(self):
        return True
    def writable(self):
        return False

    
if __name__ == "__main__":
    import io
    import sys
    
    wf = SegmentedWriteDNAFile(primer3='T'*19+'G',primer5='A'*19+'G',format_name='RS+CFC8',output="out.dna")
    
    #wf = DNAFile.open(primer3='TTTG',primer5='AAAG',format_name='FSMD',filename="out.dna2",op="w")
    
    for i in range(10):
        wf.write( "".join([chr(x) for x in convertIntToBytes(i,4)]) )

    wf.new_segment('RS+ROT','AT'+'A'*17+'G','TA'+'T'*17+'G')
        
    for i in range(10,30):
        wf.write( "".join([chr(x) for x in convertIntToBytes(i,4)]) )

    wf.close()

    rf = SegmentedReadDNAFile(primer3='T'*19+'G',primer5='A'*19+'G',input="out.dna")

    while True:
        s = rf.read(4)
        if len(s)==0:
            break
        n = convertBytesToInt([ord(x) for x in s])
        print n

    sys.exit(0)
    
        
    b = io.BytesIO()
    print b.readable(), b.writable()
    #sys.exit(0)
    readpos = 0
    for _ in range(4):
        print "-"*20
        b.seek(0,2)
        for i in range(ord('A'),ord('H')+1):
            b.write("".join([chr(i)]))
        b.seek(readpos)
        while True:
            s = b.read(1)
            if len(s)==0:
                break
            print s
        readpos = b.tell()
    print b.getvalue()
    sys.exit(0)
