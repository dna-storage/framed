from dnastorage.codec.base_conversion import convertIntToBytes,convertBytesToInt
from dnastorage.exceptions import *
from dnastorage.strand_representation import *
import math

#This is the base codec that allows for nesting of codecs, should be used for basic codecs that pass along a single strand
class BaseCodec(object):
    def __init__(self,CodecObj=None,Policy=None):
        self._Obj = CodecObj
        self.reverse=True 
        if Policy is None:
            self._Policy = AllowAll()
        else:
            self._Policy = Policy
    def set_object(self,obj):
        self._Obj=obj
    def _encode(self, s):
        return s
    def encode(self,s):
        if self.reverse: #this is the original direction
            if self._Obj != None:
                return self._encode(self._Obj.encode(s))
            else:
                return self._encode(s)
        else: #this is opposite direction
            x = self._encode(s)
            if self._Obj != None:
                return self._Obj.encode(x)
            return x
    def _decode(self, s):
        return s
    def decode(self,s):
        if self.reverse:
            s = self._decode(s)
            if self._Obj != None:
                return self._Obj.decode(s)
            else:
                return s
        else:
            
            if self._Obj != None:
                x =self._Obj.decode(s)
                return self._decode(x)
            else:
                return self._decode(s)

    def _encode_header(self):
        return []

    #Header data is byte-data that stores important properties that should survive between an encoding instantiation of a codec and the decoding instantiation
    def encode_header(self):
        if self._Obj!=None:
            return self._encode_header()+self._Obj.encode_header()
        else:
            return self._encode_header()

    def _decode_header(self,buff):
        return buff
    
    def decode_header(self,buff):
        b = self._decode_header(buff)
        if self._Obj!=None:
            return self._Obj.decode_header(b)
        else:
            return b
        
'''
Note: The BaseOuterCodec most easily works for systematic so that parity and data are easily separatd. Codes 
that modify data and/or have special indices, like non-systematic LT fountain codes, will need to be handled specially in the derived class.
A systematic code allows for easy iteration on outer codec cascades since strands uncovered from error correction are exactly strands stored in the storage system,
while non-systematic codes tranform the data from 1 domain to another, and strands before non-systematic codes will not match those after the non-systematic code.
This Base class assumes that parity added will be the same size of that of the data sub-packets. E.g sub packets of 10 strands will be assume 10 strand parity packets.
This is a reasonable assumption for things like XOR, RS, and systematic-LT type codes. I suppose it shoudld be easy enough to determine how much padding should be added to each parity
to ensure regularity at other steups. 
'''
#This is a base codec for outer encodings, has similiar structure to BaseCodec but needs to handle indexing needed for sub-packets
class BaseOuterCodec(BaseCodec):
    def __init__(self,packet_divisor,OuterCodecObj=None,Policy=None,level=None):
        BaseCodec.__init__(self,OuterCodecObj,Policy)
        self._packet_divisor = packet_divisor
        self._zero_range=() #determines the range of strands for this level that are zero-strands
        self._level=level
        self._num_data_sub_packets = 1 #determines the length of data packets before ECC is applied
        self._total_sub_packets  = 1 #determines how many sub-packets there are at this level, helps with determining indexing 
        self._index_bits=None

    def encode_header(self): #serialization for outer codecs
        buf=[]
        buf+=convertIntToBytes(self._index_bits,1)
        buf+=convertIntToBytes(self._level,1)
        buf+=convertIntToBytes(len(self._zero_range),1)
        if len(self._zero_range)>0:
            for i in self._zero_range:
                buf+=convertIntToBytes(i,4)
        buf+=convertIntToBytes(self._num_data_sub_packets,2)
        buf+=convertIntToBytes(self._total_sub_packets,2)
        return buf
        
    def decode_header(self,buff): #serialization for outer codecs
        pos=0
        self._index_bits=convertBytesToInt(buff[pos:pos+1])
        pos+=1
        self._level = convertBytesToInt(buff[pos:pos+1])
        pos+=1
        num_zero_values=convertBytesToInt(buff[pos:pos+1])
        pos+=1
        _=[]
        for i in range(num_zero_values):
            _.append(convertBytesToInt(buff[pos:pos+4]))
            pos+=4
        self._zero_range=tuple(_)
        self._num_data_sub_packets=convertBytesToInt(buff[pos:pos+2])
        pos+=2
        self._total_sub_packets=convertBytesToInt(buff[pos:pos+2])
        return buff[pos+2:]
        
    def decode(self,packet):
        #first thing to do is to group things together based on sub-packet indexes
        sub_packets={}
        for s in packet:
            sub_packet_index=s.index_ints[:self._level]
            sub_packets[tuple(sub_packet_index)]= sub_packets.get(tuple(sub_packet_index),[])+[s]
        #call down to the next level to try to patch up the packet
        if self._Obj != None:
            for s in sub_packets:
                sub_packets[s] = self._Obj.decode(sub_packets[s])
        packet=self._decode(sub_packets)
        final_packet=[]
        for s in packet:
            if not self._zero_range is () and tuple(s.index_ints)[self._level-1:]>=self._zero_range[0] and tuple(s.index_ints)[self._level-1:]<=self._zero_range[1]:
                continue
            else:
                final_packet.append(s)
        return final_packet

    
    def encode(self,packet):
        #Take strand packet, divide it by divisor, and then call the next level for each divided packet
        divided_packets=[]
        index_prefix=None
        if (len(packet)//self._packet_divisor)==0:
            raise DNAStorageError("Packet length {} too small for packet divisor {}".format(len(packet),self._packet_divisor))
        for i in range(0,len(packet),len(packet)//self._packet_divisor):
            divided_packets.append(packet[i:i+len(packet)//self._packet_divisor])
            for s in divided_packets[-1]:
                index_prefix=s.index_ints
                s.index_ints=s.index_ints+(len(divided_packets)-1,)
                if self._level==None:
                    self._level=len(s.index_ints)
                
        length_before_ecc = len(divided_packets)
        pad = None
        pad_length=0
        pad_start_index=None
        pad_end_boundary=None
        first_pad_strand=None
        #Note down some padding that needs to be done to fill out the outer encoding, code be useful to note down in metadata
        if not len(packet)%self._packet_divisor==0:
            pad_start_index=len(divided_packets[-1])
            while len(divided_packets[-1])<len(divided_packets[-2]):
                divided_packets[-1].append(BaseDNA(codewords=[0]*len(packet[0].codewords),index_ints=divided_packets[-1][0].index_ints))
                pad_length+=1
                divided_packets[-1][-1].is_zero=True #flag to help filter out zeros
            pad_end_boundary=pad_length+pad_start_index
            pad=divided_packets[-1][0].index_ints

        self._num_data_sub_packets=len(divided_packets)
        #put ecc on these sub packets
        divided_packets=self._encode(divided_packets)
        for index, x in enumerate(divided_packets[length_before_ecc::]):
            for s in x:
                s.index_ints=index_prefix+(index+length_before_ecc,)
                 
        '''
        Note: the way indexes are handled here may not be the most time efficient, O(Number of strands * Layers) number of operations
        '''
        if self._index_bits is None:
            self._index_bits = math.ceil(math.log(len(divided_packets),2))
            self._total_sub_packets=len(divided_packets)

        else:
            if not self._index_bits== math.ceil(math.log(len(divided_packets),2)): raise DNAStorageError("Indexes not consistent")
            if not self._total_sub_packets==len(divided_packets): raise DNAStorageError("Total Subpackets Inconsistent")
            
        protected_sub_packets=[]
        #call next level to protect the divided packets
        if self._Obj !=None:
            for dp in divided_packets:
                protected_sub_packets.append(self._Obj.encode(dp))

        else:
            protected_sub_packets=divided_packets

        #gets rid of 0 strands and notes their positions
        out_packet=[]
        for i, dp in enumerate(protected_sub_packets):
            p=[]
            zeros=[]
            for sub_packet_index, s in enumerate(dp):
                if hasattr(s,'is_zero') and s.is_zero and i==(length_before_ecc-1) and pad!=None and sub_packet_index>=pad_start_index and sub_packet_index<pad_end_boundary:
                    zeros.append(s)
                    continue
                else:
                    p.append(s) #only want to synthesize non zero strands
                #now write down the zero range
            if self._zero_range==() and len(zeros)>0:
                self._zero_range=(zeros[0].index_ints[self._level-1:],zeros[-1].index_ints[self._level-1:])
            out_packet=out_packet+p
        return out_packet
        
    def get_index_bits(self):
        #recursively get the number of bits needed for indexing
        if self._Obj !=None:
            return (self._index_bits,)+self._Obj.get_index_bits()
        else:
            return (self._index_bits,)


    def get_zero_pad(self):
        #return location of zero pad strands so it can be possibly noted by the file system for decoding.
        if self._Obj !=None:
            return (self._zero_range,)+self._Obj.get_zero_pad()
        else:
            return (self.zero_range,)


    def set_zero_pad(self,pad_list):
        self._zero_range=pad_list[0]
        if self._Obj !=None:
            return self._Obj.set_zero_pad(pad_list[1::])


    def filter_data(self, strands):
        #filters strands that fall in this levels data range
        filtered_strands=[]
        for s in strands:
            level_index = s.index_ints[self._level-1]
            if level_index<self._num_data_subpackets:
                filtered_strands.append(s)
        if self._Obj!=None:
            return self._Obj.filter_data(filtered_strands)
        else:
            return filtered_strands

    def get_next_index(self,index):
        listed_index=list(index)
        is_zero=False
        #given an index try to determine the next index, and tell whether if its a zero strand or not
        if self._Obj!=None:
            after_index,is_zero = self._Obj.get_next_index(index)
        else:
            after_index=listed_index
            after_index[self._level-1]+=1
        this_level_index = after_index[self._level-1]
        if this_level_index==self._total_sub_packets:
            after_index[self._level-1]=0
            after_index[self._level-2]+=1
        #check if this level index is actually in a zero range, allows for easy detection of zero strands 
        if not is_zero and self._zero_range is not () and (tuple(after_index)[self._level-1:]>=self._zero_range[0] and tuple(after_index)[self._level-1:]<=self._zero_range[1]):
            is_zero=True
        else:
            is_zero=False
        return after_index,is_zero
            
    def get_previous_index(self,index):
        #given an index try to determine the previous index, and tell whether if its a zero strand or not
        listed_index=list(index)
        is_zero=False
        #given an index try to determine the next index, and tell whether if its a zero strand or not
        if self._Obj!=None:
            after_index,is_zero = self._Obj.get_previous_index(index)
        else:
            after_index = listed_index
            after_index[self._level-1]-=1
        this_level_index = after_index[self._level-1]
        if this_level_index<0:
            after_index[self._level-1]=self._total_sub_packets
            after_index[self._level-2]-=1
        #check if this level index is actually in a zero range, allows for easy detection of zero strands 
        if not is_zero and self._zero_range is not () and (tuple(after_index)[self._level-1:]>=self._zero_range[0] and tuple(after_index)[self._level-1:]<=self._zero_range[1]):
            is_zero=True
        return after_index,is_zero

    def get_total_divisor(self):
        if self._Obj != None:
            return self._packet_divisor*self._Obj.get_total_divisor()
        else:
            return self._packet_divisor

    def valid(self,index):
        if self._Obj==None:
            return index[self._level-1]<self._total_sub_packets 
        else:
            return index[self._level-1]<self._total_sub_packets and self._Obj.valid(index)
    def is_zero(self,index):
        if self._zero_range is ():
            return False
        if self._Obj==None:
            return (tuple(after_index)[self._level-1:]>=self._zero_range[0] and tuple(after_index)[self._level-1:]<=self._zero_range[1])
        else:
            return (tuple(after_index)[self._level-1:]>=self._zero_range[0] and tuple(after_index)[self._level-1:]<=self._zero_range[1]) or self._Obj.is_zero(index)
    
class TableCodec(BaseCodec):
    def __init__(self,CodecObj=None,keyEncWidth=20,keyDecWidth=4,cwEncWidth=5,cwDecWidth=1,Policy=None):
        BaseCodec.__init__(self,CodecObj=CodecObj,Policy=Policy)
        self._keyEncWidth = keyEncWidth
        self._keyDecWidth = keyDecWidth
        self._cwEncWidth = cwEncWidth
        self._cwDecWidth = cwDecWidth
        assert keyEncWidth % cwEncWidth == 0
        assert keyEncWidth / cwEncWidth == keyDecWidth

    def _enctab(self, val):
        assert 0 and "not implemented"

    def _dectab(self, seq):
        assert 0 and "not implemented"

    def _encode(self, s):
        key = convertIntToBytes(s[0],self._keyDecWidth)
        payload = bytearray(s[1])
        #print ("{}:{}".format(key,payload))
        enc = []
        for i in range(0,len(key),self._cwDecWidth):
            enc.append( self._enctab( convertBytesToInt(key[i:i+self._cwDecWidth])) )
        for i in range(0,len(payload),self._cwDecWidth):
            enc.append( self._enctab( convertBytesToInt(payload[i:i+self._cwDecWidth])))
        return "".join(enc)

    def _decode(self, s):
        bytes = []
        for i in range(0,len(s),self._cwEncWidth):
            val = self._dectab(s[i:i+self._cwEncWidth])
            val = convertIntToBytes(val,self._cwDecWidth)
            bytes = bytes + val
        key = convertBytesToInt(bytes[0:self._keyDecWidth])
        #bytes = [ chr(b) for b in bytes[self._keyDecWidth:] ]
        return [key,bytes[self._keyDecWidth:]]
        #return [key,bytearray(bytes)]
