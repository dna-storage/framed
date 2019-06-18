from dnastorage.codec.base_conversion import convertIntToBytes,convertBytesToInt
from dnastorage.arch.builder import *
import editdistance as ed
#from dnastorage.primer.primer_util import edit_distance
from io import BytesIO
from dnastorage.util.file_support import *
import math
import struct
from dnastorage.system.formats import *

### Designed to fit on a single strand for most use cases
### 
### Every header strand begins with special sequence that can't be used at the beginning of indices: CCATCCAT
###
### 1. 'CCATCCAT' [1]
### 2. short index - usually 0, 0-255, at most 16 strands [1]
### 3. major version (0-255) [1]
### 4. minor version (0-255) [1]
### 5. num bytes for size [1]
### 6. size [x]
### 7  num bytes for original filename [2]
### 8. null terminated string
### 9. encoding style        [2]
### 10. length of remaining record (2 bytes)      [2]
### 11. remaining record byte encoded  [?]
### Pad to final width using arbitrary sequence

system_version = { 'major': 0, 'minor':1 }
magic_header = 'CCATCCAT'

def encode_primer_diff(o,n):
    hdr = []
    baseVal = { 'A': 0, 'C': 1, 'G':2, 'T':3 }                            
    assert len(o)==len(n)
    for i,(oo,nn) in enumerate(zip(o,n)):
        if oo != nn:                                
            hdr +=  [ (baseVal[nn] << 6) | (i&0x3F) ]
    if len(hdr) == 0:
        return [0]
    hdr = [len(hdr)] + hdr
    return hdr

def decode_primer_diff(data,oprimer):
    sz = data[0]
    if sz==0:
        return oprimer,1    
    baseVal = [ 'A', 'C', 'G', 'T' ]
    nprimer = [ _ for _ in oprimer ]
    for i in range(sz):
        val = data[1+i]
        base = baseVal[(val&0xC0)>>6]
        pos = val&0x3F
        nprimer[pos] = base
    return "".join(nprimer),sz+1

def encode_size_and_value(val):
    data = []
    if val==0:
        data += [1, 0]
    else:
        data += convertIntToBytes(int(math.ceil(math.log(val+1,2)/8.0)),1)
        data += convertIntToBytes(int(val),int(math.ceil(math.log(val+1,2)/8.0)))
    return data

def decode_size_and_value(data,pos):
    size_bytes = data[pos]
    val = convertBytesToInt(data[pos+1:pos+1+size_bytes])
    return val,size_bytes+1


def encode_file_header_comments(filename,format_id,size,other_data,primer5,primer3):
    comment = "% dnastorage version {}.{}\n".format(system_version['major'],system_version['minor'])
    comment += "% {} \n".format(size)
    if len(filename) > 0:
        comment += "% {} \n".format(filename)
    else:
        comment += "% No filename recorded.\n"

    comment += "% 5-{} 3-{}\n".format(primer5,primer3)
    comment += "% Id-{} Description-{} \n".format(format_id,file_system_format_description(format_id))
    comment += "% {} bytes of additional data \n".format(len(other_data))
    return comment

def encode_file_header(filename,format_id,size,other_data,primer5,primer3):
    data =  [ system_version['major'], system_version['minor'] ]
    data += convertIntToBytes(int(math.ceil(math.log(size,2)/8.0)),1)
    data += convertIntToBytes(size,int(math.ceil(math.log(size,2)/8.0)))
    data += convertIntToBytes(len(filename)+1,2)
    data += [ ord(_) for _ in filename ] + [0]        
    data += convertIntToBytes(format_id,2)
    data += convertIntToBytes(len(other_data),2)
    data += other_data
    data += [0]*(80-len(data))
    data = "".join([chr(x) for x in data])

    pf = ReadPacketizedFilestream(BytesIO(data))
    enc_func = file_system_encoder_by_abbrev('FSMD')
    enc = enc_func(pf,primer5+magic_header,primer3)
    strands = []
    for e in enc:
        strands.append(e)
    
    return strands

def pick_nonheader_strands(strands,primer5):
    others = []
    picks = []
    for s in strands:
        if s.startswith(primer5):            
            if s.startswith(primer5+magic_header):
                pass
            else:
                others.append(s)
        else:
            others.append(s)
    return others

def pick_header_strands(strands,primer5):
    picks = []
    others = []
    for s in strands:
        if s.startswith(primer5+magic_header):
            picks.append(s)
        elif s.startswith(primer5):
            plen= len(primer5)
            possible_hdr = s[plen:plen+len(magic_header)]
            if ed.eval(possible_hdr,magic_header) < 2:
                ss = s[:]
                ss[plen:plen+len(magic_header)] = magic_header
                picks.append(ss)
            else:
                others.append(s)
    return picks,others


def decode_file_header(strands,primer5,primer3):
    picks,_ = pick_header_strands(strands,primer5)

    b = BytesIO()
    pf = WritePacketizedFilestream(b,80,20)
    dec_func = file_system_decoder_by_abbrev('FSMD')
    dec = dec_func(pf,primer5+magic_header,primer3)
    
    for s in picks:
        dec.decode(s)

    dec.write()
    assert dec.complete    

    data = [ ord(x) for x in b.getvalue() ]
    
    assert data[0] == system_version['major']
    assert data[1] == system_version['minor']

    header = {}
    header['version'] = [ data[0], data[1] ]

    size_bytes = data[2]
    pos = 3
    header['size'] = convertBytesToInt(data[pos:pos+size_bytes])
    pos += size_bytes

    size_filename = convertBytesToInt(data[pos:pos+2])
    pos+=2
    header['filename'] = "".join([chr(x) for x in data[pos:pos+size_filename]])
    pos += size_filename

    header['formatid'] = convertBytesToInt(data[pos:pos+2])
    pos += 2
    
    size_other_data = convertBytesToInt(data[pos:pos+2])
    pos += 2
    header['other_data'] = [ x for x in data[pos:pos+size_other_data] ]
    
    return header

if __name__ == "__main__":
    strands = encode_file_header("",0xA,2,[1,2,3,4],"A"*19+"G","T"*19+"G")
    for s in strands:
        print "{}: strand={}".format(len(s), s)
    
    print decode_file_header(strands,"A"*19+"G","T"*19+"G")
