#!/usr/bin/python
from dnastorage.codec.base import *
from dnastorage.codec import base_conversion
from dnastorage.codec import rs
from random import randint
from dnastorage.codec.codecfile import *
from math import ceil

rs_initialized = False

def rs_encode(message,n,k):
    global rs_initialized
    if not rs_initialized:
        rs.init_tables(0x11d)
        rs_initialized = True
    assert k == len(message)
    return rs.rs_encode_msg(message,n-k);

def rs_decode(message,n,k):
    global rs_initialized
    if not rs_initialized:
        rs.init_tables(0x11d)
        rs_initialized = True
    assert k == len(message)
    erasures = [ i for i in range(len(message)) if message[i]==-1 ]
    # correct the message
    corrected_message, corrected_ecc = rs.rs_correct_msg(message,n-k, erase_pos=erasures)
    return corrected_message

class ReedSolomonInnerCodec(BaseCodec):
    """
    ReedSolomonInnerCodec takes a key,value pair as input to the _encode function and
    produces a Reed-Solomon encoded message as a byte array. Both key and value are
    protected. If the value is shorter than expected, it is padded to padWidth with
    random data.

    This is an "Inner" Codec because it only can correct errors within a strand.

    This class is hard coded to use GF(256).
    """
    def __init__(self,numberBytes,lengthMessage,CodecObj=None):
        BaseCodec.__init__(self,CodecObj)
        global rs_initialized
        if not rs_initialized:
            rs.init_tables(0x11d)
            rs_initialized = True
        self._numberBytes = numberBytes
        self._lengthMessage = lengthMessage

        assert(self._lengthMessage <= 255) # required by GF(256)

    """
    packet is a tuple with the key at position 0 and value at position 1. This should return
    the Reed-Solomon encoded byte array.
    """
    def _encode(self,packet):
        assert packet[0] < 256**3
        assert len(packet[1])<=self._numberBytes

        l = base_conversion.convertIntToBytes(packet[0],3)
        a = bytearray(packet[1])

        # construct message
        message = l + [x for x in a]

        # pad the message (wish we didn't have to do this!)
        if len(message) < self._numberBytes:
            message += [ randint(0,10) for _ in range(self._numberBytes - len(message)) ]

        k = self._numberBytes

        # encoded the message using the RS library
        mesecc = rs.rs_encode_msg(message, self._lengthMessage - len(message));
        # separate out key and value
        # convert mesecc into a string
        return packet[0],"".join([ chr(_) for _ in mesecc[3:]])

    """
    This function expects a list of unsigned integers in GF(256). For now, erasures are
    denoted with -1.
    """
    def _decode(self,packet):
        l = base_conversion.convertIntToBytes(packet[0],3)
        a = [x for x in bytearray(packet[1])]
        message = l + a
        # find the -1s in the list
        erasures = [ i for i in range(len(message)) if message[i]==-1 ]
        # correct the message
        corrected_message, corrected_ecc = rs.rs_correct_msg(message,self._lengthMessage-self._numberBytes, erase_pos=erasures)
        value = corrected_message[3:]
        return packet[0],value


class ReedSolomonInnerOuterEncoder(EncodePacketizedFile):
    """
    ReedSolomonInnerOuterCodec takes a key,value pair as input to the _encode
    function and produces a Reed-Solomon encoded message as a byte array. Both key
    and value are protected. If the value is shorter than expected, it is padded
    to padWidth with random data.

    This is an Inner-Outer Codec. It protects within strands and across strands.

    In this implementation, we specify the number of errors to correct accross
    strands, e_outer, and the number within strands, e_inner. We add 2*e_outer
    error correction symbols for each 255-2*e_outer strands and 2*e_inner error
    correction bytes per strand.

    This class is hard coded to use GF(256). This limits the error correction to
    255 length messages.

    Arguments:
    k_datastrand - # of data bytes in a strand
    e_inner  - # of error correction codes used to correct in a strand
    k_index  - # of data bytes to devote to address
    k_outer  - # of strands in one block protected by the RS outer code
    e_outer  - # of error correction codes used to protect an outer block

    The total length of a data strand will be k_index + k_datastrand + e_inner.
    The total size of a outer block will k_datastrand*k_outer protected by k_datastrand*e_outer codes.
    """
    def __init__(self,pf,CodecObj=None,k_datastrand=15,e_inner=2,k_index=3,k_outer=251,e_outer=4):
        EncodePacketizedFile.__init__(self,pf,CodecObj)
        global rs_initialized
        if not rs_initialized:
            rs.init_tables(0x11d)
            rs_initialized = True
        self._k_strand = k_datastrand
        self._k_index = k_index
        self._e_inner = e_inner
        self._k_outer = k_outer
        self._e_outer = e_outer
        self._packetizedFile.packetSize = k_outer*k_datastrand
        self._packetizedFile._RS=True;
        assert (k_datastrand + k_index + e_inner  <= 255) # required by GF(256)
        assert (k_outer + e_outer  <= 255) # required by GF(256)
        self.strand_length = k_datastrand + k_index + e_inner
        self.index = 0
        self.strands = []
        self._transpose = False
    """
    Overload function to prevent incorrect behavior
    """
    def _encode(self,packet):
        assert 0
        return None

    def _pop_strand(self):
        return self.strands.pop(0)

    def encode(self):
        if len(self.strands):
            return self._pop_strand()

        # self._packetizedFile will pad it to match the requested packetSize.
        raw = self._packetizedFile.next()
        raw = [ _ for _ in raw ]

        #pad raw out to make sure we have enough data for strands
        while(len(raw)%self._k_strand != 0):
            raw.append('\x00')

        matrix = []

        # distribute adjacent bytes in the file across strands to make loss of a strand
        # less severe
        if self._transpose:
            pass

        # compute inner code
        #Changed _k_strand*_k_outer to len(raw) !! May need to fix this
        for x in range(0,len(raw),self._k_strand):
            r = raw[x:x+self._k_strand]
            ind = base_conversion.convertIntToBytes(x/self._k_strand + self.index,self._k_index)
            message = ind + [x for x in bytearray(r)]
            mesecc = rs.rs_encode_msg(message, self._e_inner);
            matrix.append(mesecc)


        self.index += self._k_outer*self._k_strand/self._k_strand

        num_real_strands=len(matrix)
        #pad out the matrix with dummy strands
        dummy_index=num_real_strands
        while(len(matrix)<self._k_outer):
            dummy=['\x00']*(self._k_strand)
            ind=base_conversion.convertIntToBytes(dummy_index,self._k_index)
            dummy_message=ind+[_ for _ in bytearray(dummy)]
            dummy_mesecc=rs.rs_encode_msg(message,self._e_inner)
            matrix.append(dummy_mesecc)
            dummy_index+=1


        # number of error strands
        n_error_strands = int(ceil(self.strand_length * self._e_outer / float(self._k_strand)))


        # outer error codes
        error_codes = []

        # compute outer code over all bytes, including index
        for x in range(self.strand_length):
            message = [ matrix[_][x] for _ in range(len(matrix)) ]
            #print "{} - {}".format(len(message),message)
            mesecc = rs.rs_encode_msg(message,self._e_outer);
            error_codes += mesecc[-self._e_outer:]
            #print "{}".format(mesecc[-self._e_outer:])

        #print "encode: {} {} {}".format(len(error_codes),error_codes,n_error_strands)

        # pad with zeros
        if len(error_codes) % self._k_strand != 0:
            pad = [ 0 for _ in range(int(self._k_strand - len(error_codes) % self._k_strand)) ]
            error_codes += pad

        for i in range(n_error_strands):
            err = [ error_codes[j] for j in range(i,len(error_codes),n_error_strands) ]
            ind = base_conversion.convertIntToBytes(self.index,self._k_index)
            message = ind + err
            mesecc = rs.rs_encode_msg(message, self._e_inner);
            #print "{} mesecc = {}".format(len(mesecc),mesecc)
            matrix.append(mesecc)
            self.index += 1

        for i, m in enumerate(matrix):
            #only append strands if they are real and are error correction strands
            if i<num_real_strands or i>=self._k_outer:
                tup = (base_conversion.convertBytesToInt(m[:self._k_index]),m[self._k_index:])
                #codecs expect a (index,value) tuple
                self.strands.append(self._Codec.encode(tup))

        return self._pop_strand()



class ReedSolomonInnerOuterDecoder(DecodePacketizedFile):
    """
    ReedSolomonInnerOuterDecoder takes a key,value pair as input to the decode
    function and produces a Reed-Solomon decoded message as a byte array. Both key
    and value are protected.

    This is an Inner-Outer Codec. It protects within strands and across strands. The parameters
    to the class must match the ReedSolomonInnerOuterEncoder exactly for decoding to correct
    errors.

    In this implementation, we specify the number of errors to correct accross
    strands, e_outer, and the number within strands, e_inner. We add 2*e_outer
    error correction symbols for each 255-2*e_outer strands and 2*e_inner error
    correction bytes per strand.

    This class is hard coded to use GF(256). This limits the error correction to
    255 length messages.

    Arguments:
    k_datastrand - # of data bytes in a strand
    e_inner  - # of error correction codes used to correct in a strand
    k_index  - # of data bytes to devote to address
    k_outer  - # of strands in one block protected by the RS outer code
    e_outer  - # of error correction codes used to protect an outer block

    The total length of a data strand will be k_index + k_datastrand + e_inner.
    The total size of a outer block will k_datastrand*k_outer protected by k_datastrand*e_outer codes.
    """
    def __init__(self,pf,CodecObj=None,k_datastrand=15,e_inner=2,k_index=3,k_outer=251,e_outer=4):
        DecodePacketizedFile.__init__(self,pf,CodecObj)
        global rs_initialized
        if not rs_initialized:
            rs.init_tables(0x11d)
            rs_initialized = True
        self._file_size=self._packetizedFile.size
        self._file_size=self._packetizedFile.size
        self._k_strand = k_datastrand
        self._k_index = k_index
        self._e_inner = e_inner
        self._k_outer = k_outer
        self._e_outer = e_outer
        self._packetizedFile.packetSize = k_outer*k_datastrand
        assert (k_datastrand + k_index + e_inner  <= 255) # required by GF(256)
        assert (k_outer + e_outer  <= 255) # required by GF(256)
        self.strand_length = k_datastrand + k_index + e_inner
        self.index = 0
        self.strands = []
        self._transpose = False
        self._rsMap = {}
        self.n_error_strands = int(ceil(self.strand_length * self._e_outer / float(self._k_strand)))
        self.outer_block = k_outer+self.n_error_strands
        self.decodedMap = {}
        self._build_dummy_strands()



    #need to build and insert dummy strands into rsMap
    def _build_dummy_strands(self):
        block_size_B=self._packetizedFile.packetSize
        last_block_size_B=self._file_size%block_size_B
        num_real_data_strands_last_block=int(ceil(last_block_size_B/float(self._k_strand)))
        num_blocks=int(ceil(self._file_size/float(block_size_B)))
        strands_per_block_data=int(ceil(block_size_B/float(self._k_strand)))
        strands_per_block_error=self.n_error_strands
        #total_strands_per_block will be the number of indexes per block
        total_strands_per_block=strands_per_block_data+strands_per_block_error

        #this is the first index in the last block
        last_block_base_index=total_strands_per_block*(num_blocks-1)

        #this is the first index of dummy strands
        dummy_start_index=last_block_base_index+num_real_data_strands_last_block

        #upper_bound for the dummy strands
        dummy_upper_bound=(total_strands_per_block*num_blocks)-strands_per_block_error


        #insert dummy values into the rsMap table
        for index in range(dummy_start_index,dummy_upper_bound):
            dummy_index=base_conversion.convertIntToBytes(index,self._k_index)
            dummy_data=['\x00']*self._k_strand
            dummy_message=dummy_index+[_ for _ in bytearray(dummy_data)]
            dummy_mesecc=rs.rs_encode_msg(dummy_message,self._e_inner)
            self._rsMap[index]=dummy_mesecc



    def _get_base(self, index):
        base = index/self.outer_block*self.outer_block
        return base

    def _is_block_decoded(self,index):
        base = self._get_base(index)
        return self.decodedMap.has_key(base)

    def _is_block_ready_to_decode(self,index):
        base = self._get_base(index)
        for x in range(base,base+self.outer_block):
            if not self._rsMap.has_key(x):
                return False
        # change this later to allow for missing strands
        return True

    def is_error_free(self, message, nsym):
       # print message
        if None in message:
            return True
        return max(rs.rs_calc_syndromes(message,nsym)) == 0

    def check_inner_strand(self, message):
        # might need try-catch here
        erasures = [ i for i in range(len(message)) if message[i]==-1 or message[i]==None ]
        # correct the message

        try:
            corrected_message, corrected_ecc = rs.rs_correct_msg(message,self._e_inner, erase_pos=erasures)
        except Exception:
            corrected_message=message[0:self._k_strand+self._k_index]
            corrected_ecc=message[-self._e_inner:]


        return corrected_message+corrected_ecc

    def _collect_error_codes(self, error_strands):
        collect = []
        es = []
        for e in error_strands:
            collect += e[self._k_index:self._k_index+self._k_strand]
        assert len(collect) == self.n_error_strands * self._k_strand
        rearranged = []
        for i in range(self.strand_length):
            rearranged += [ collect[j] for j in range(i,len(collect),self._k_strand) ]
            es += [ collect[j] for j in range(i,len(collect),self._k_strand) ]

        rearranged = [ es[i:i+self._e_outer] for i in range(0,self.strand_length*self._e_outer,self._e_outer) ]

        #print collect
        #print rearranged
        return rearranged


    def _find_erasures(self, message):
        return [ i for i in range(len(message)) if message[i]==-1 or message[i]==None ]

    def _getBlockIndex(self, index):
        block = self._get_base(index) / self.outer_block
        return block

    def _decode_block(self,index):
        base = self._get_base(index)
        matrix = []
        for x in range(base,base+self.outer_block-self.n_error_strands):
            matrix.append(self._rsMap[x])

        # check for inner errors on all the strands
        #for m in matrix:
            # better way?
            #if not self.is_error_free(m,self._e_inner):
                #s = self.check_inner_strand(m)
                #for i in range(len(s)):
                    #m[i] = s[i]

        error_strands = []
        for x in range(base+self.outer_block-self.n_error_strands,base+self.outer_block):
            s = self._rsMap[x]
            # check for inner errors on the error_strands
            #if not self.is_error_free(s,self._e_inner):
             #   s = self.check_inner_strand(s)
            error_strands.append(s)

        # check for outer errors on the block
        #   - note: this could modify an index
        error_codes = self._collect_error_codes(error_strands)
        #print error_codes

        # compute outer code over all bytes, including index
        for x in range(self.strand_length):
            message = [ matrix[_][x] for _ in range(len(matrix)) ]
            message += error_codes[x]

            #print "before messsage = {}".format(message)

            erasures = self._find_erasures(message)
            # correct the message
            #print "messsage = {} {}".format(message,erasures)

            #if number of erasures is too big give up on correction process
            try:
                corrected_message, corrected_ecc = rs.rs_correct_msg(message,self._e_outer,erase_pos=erasures)
            except Exception:
                corrected_message=message[0:self._k_outer]
                corrected_ecc=message[-self._e_outer:]

            #print corrected_message
            for i in range(len(corrected_message)):
                # correct the matrix
                matrix[i][x] = corrected_message[i]


        # got here, everything worked out!
        # form large packet
        key = self._getBlockIndex(index)
        value=""
        for m in matrix:
            slice_data=m[self._k_index:self._k_index+self._k_strand]
            for data_element in slice_data:
                if data_element is None or  data_element==-1:
                    value+='\x00'
                else:
                    value+=chr(data_element)

        #value += "".join([ chr(_) for _ in m[self._k_index:self._k_index+self._k_strand] ])

        assert len(value) == self._packetizedFile.packetSize
        self.writeToFile(key,value)
        self.decodedMap[base] = True

    def _decode(self,key,value):
        # merge index bytes and value bytes
        value2 = base_conversion.convertIntToBytes(key,self._k_index)
        # expand string into bytes
        #value2 += [ ord(x) for x in value ]
        value2+=value

        s=value2
        #do inner checking before placement into RS table, hopefully fix index errors before placement
        if not self.is_error_free(value2,self._e_inner):
                s = self.check_inner_strand(value2)
        _key=base_conversion.convertBytesToInt(s[:self._k_index])
        self._rsMap[_key] = s
        #print "rsMap[{}] = {}".format(key,value2)
        if self._is_block_ready_to_decode(_key) and not self._is_block_decoded(_key):
            self._decode_block(_key)
        return

    """
    If a block is incomplete, it will not yet have been decoded. Call this
    function to detect when strands are missing and it will insert erasures.
    After that, decoding may be able to succeed using the outer codes.
    """
    def attempt_final_decoding(self):
        numBlocks = int(ceil(float(self._packetizedFile.size) / self._packetizedFile.packetSize))
        for block in range(0,numBlocks*self.outer_block,self.outer_block):
            if self._is_block_decoded(block):
                continue
            if not self._is_block_ready_to_decode(block):
                # figure out why
                base = block
                found = 0
                for x in range(base,base+self.outer_block):
                    if not self._rsMap.has_key(x):
                        # create erasure block
                        erased = [-1 for _ in range(self.strand_length) ]
                        self._rsMap[x] = erased
                        found += 1
                assert self._is_block_ready_to_decode(block)==True
            # try to decode block in either case
            self._decode_block(block)

    def write(self):
        self.attempt_final_decoding()
        DecodePacketizedFile.write(self)

    #dummy write allows us to get a data structure back for comparison
    def dummy_write(self):
        self.attempt_final_decoding()
        return DecodePacketizedFile.dummy_write(self)

if __name__ == "__main__":
    import sys
    from random import randint
    from dnastorage.codec.codecfile import *
    import rs

    pf = ReadPacketizedFile(sys.argv[1])
    wpf = WritePacketizedFile("out.d",pf.size,15)

    enc = ReedSolomonInnerOuterEncoder(pf)
    dec = ReedSolomonInnerOuterDecoder(wpf)

    for i,s in enumerate(enc):
        #print "{}. ({}) - {}".format(s[0],len(s[1]),[ord(x) for x in s[1]])
        t = [ord(x) for x in s[1]]
        if randint(0,700) < 400:
            print "Alter strand!"
            t[randint(0,16)] = randint(0,255)

        if randint(0,1000) < 10:
            print "Throw away strand!"
        else:
            t2 = "".join([chr(x) for x in t])
            dec.decode((s[0],t2))

    dec.attempt_final_decoding()
    dec.write()

    sys.exit(0)
    pf = ReadPacketizedFile(sys.argv[1])
    pf.packetSize = 15

    r = ReedSolomonInnerCodec(18,20)
    enc = EncodePacketizedFile(pf,r)

    wpf = WritePacketizedFile("out.d",pf.size,None)
    dpf = DecodePacketizedFile(wpf,r)

    for i,s in enumerate(enc):
        print "{}. ({}) - {}".format(i,s[0],[ord(x) for x in s[1]])
        dpf.decode(s)

    dpf.write()
