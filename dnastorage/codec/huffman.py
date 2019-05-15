#!/usr/bin/python
from dnastorage.codec.base import *
from dnastorage.codec import base_conversion
from dnastorage.codec.huffman_table import HuffmanTable

import unittest

class HuffmanEncTester(unittest.TestCase):
    """Compare table generated by HuffmanTable to the textual table here"""
    """They should match."""
    def test_tables_match(self):
        syms = [ x for x in range(256) ]
        w = [ 0.1 for x in range(256) ]
        for i in range(ord('A'),ord('z'),1):
            w[i] = 0.2
            w[0] = 0.2

        ht3 = HuffmanTable(3, ['0','1', '2'], syms, w)
        enc,dec = ht3.get_tables()

        for k in range(256):
            assert enc[k]==huffman_enc_table[k]

_huffman_enc_map = { '0': 'A' ,
           '1': 'C' ,
           '2': 'G' }

_huffman_dec_map = { 'A':'0' ,
           'C':'1' ,
           'G':'2' }

# fix me! make it A,C, G and not 0,1,2
huffman_enc_table = [
'02010', '110210', '110211', '110212', '110220', '110221', '110222', '111000', 
'111001', '111002', '111010', '111011', '111012', '111020', '111021', '111022', 
'111100', '111101', '111102', '111110', '111111', '111112', '111120', '111121', 
'111122', '111200', '111201', '111202', '111210', '111211', '111212', '111220', 
'111221', '111222', '112000', '112001', '112002', '112010', '112011', '112012', 
'112020', '112021', '112022', '112100', '112101', '112102', '112110', '112111', 
'112112', '112120', '112121', '112122', '112200', '112201', '112202', '112210', 
'112211', '112212', '112220', '112221', '112222', '120000', '120001', '120002', 
'120010', '02011', '02012', '02020', '02021', '02022', '02100', '02101', 
'02102', '02110', '02111', '02112', '02120', '02121', '02122', '02200', 
'02201', '02202', '02210', '02211', '02212', '02220', '02221', '02222', 
'10000', '10001', '10002', '10010', '10011', '10012', '10020', '10021', 
'10022', '10100', '10101', '10102', '10110', '10111', '10112', '10120', 
'10121', '10122', '10200', '10201', '10202', '10210', '10211', '10212', 
'10220', '10221', '10222', '11000', '11001', '11002', '11010', '11011', 
'11012', '11020', '120011', '120012', '120020', '120021', '120022', '120100', 
'120101', '120102', '120110', '120111', '120112', '120120', '120121', '120122', 
'120200', '120201', '120202', '120210', '120211', '120212', '120220', '120221', 
'120222', '121000', '121001', '121002', '121010', '121011', '121012', '121020', 
'121021', '121022', '121100', '121101', '121102', '121110', '121111', '121112', 
'121120', '121121', '121122', '121200', '121201', '121202', '121210', '121211', 
'121212', '121220', '121221', '121222', '122000', '122001', '122002', '122010', 
'122011', '122012', '122020', '122021', '122022', '122100', '122101', '122102', 
'122110', '122111', '122112', '122120', '122121', '122122', '122200', '122201', 
'122202', '122210', '122211', '122212', '122220', '122221', '122222', '00000', 
'00001', '00002', '00010', '00011', '00012', '00020', '00021', '00022', 
'00100', '00101', '00102', '00110', '00111', '00112', '00120', '00121', 
'00122', '00200', '00201', '00202', '00210', '00211', '00212', '00220', 
'00221', '00222', '01000', '01001', '01002', '01010', '01011', '01012', 
'01020', '01021', '01022', '01100', '01101', '01102', '01110', '01111', 
'01112', '01120', '01121', '01122', '01200', '01201', '01202', '01210', 
'01211', '01212', '01220', '01221', '01222', '02000', '02001', '02002'
]

huffman_dec_table = {
'02010' : 0, '110210' : 1, '110211' : 2, '110212' : 3, 
'110220' : 4, '110221' : 5, '110222' : 6, '111000' : 7, 
'111001' : 8, '111002' : 9, '111010' : 10, '111011' : 11, 
'111012' : 12, '111020' : 13, '111021' : 14, '111022' : 15, 
'111100' : 16, '111101' : 17, '111102' : 18, '111110' : 19, 
'111111' : 20, '111112' : 21, '111120' : 22, '111121' : 23, 
'111122' : 24, '111200' : 25, '111201' : 26, '111202' : 27, 
'111210' : 28, '111211' : 29, '111212' : 30, '111220' : 31, 
'111221' : 32, '111222' : 33, '112000' : 34, '112001' : 35, 
'112002' : 36, '112010' : 37, '112011' : 38, '112012' : 39, 
'112020' : 40, '112021' : 41, '112022' : 42, '112100' : 43, 
'112101' : 44, '112102' : 45, '112110' : 46, '112111' : 47, 
'112112' : 48, '112120' : 49, '112121' : 50, '112122' : 51, 
'112200' : 52, '112201' : 53, '112202' : 54, '112210' : 55, 
'112211' : 56, '112212' : 57, '112220' : 58, '112221' : 59, 
'112222' : 60, '120000' : 61, '120001' : 62, '120002' : 63, 
'120010' : 64, '02011' : 65, '02012' : 66, '02020' : 67, 
'02021' : 68, '02022' : 69, '02100' : 70, '02101' : 71, 
'02102' : 72, '02110' : 73, '02111' : 74, '02112' : 75, 
'02120' : 76, '02121' : 77, '02122' : 78, '02200' : 79, 
'02201' : 80, '02202' : 81, '02210' : 82, '02211' : 83, 
'02212' : 84, '02220' : 85, '02221' : 86, '02222' : 87, 
'10000' : 88, '10001' : 89, '10002' : 90, '10010' : 91, 
'10011' : 92, '10012' : 93, '10020' : 94, '10021' : 95, 
'10022' : 96, '10100' : 97, '10101' : 98, '10102' : 99, 
'10110' : 100, '10111' : 101, '10112' : 102, '10120' : 103, 
'10121' : 104, '10122' : 105, '10200' : 106, '10201' : 107, 
'10202' : 108, '10210' : 109, '10211' : 110, '10212' : 111, 
'10220' : 112, '10221' : 113, '10222' : 114, '11000' : 115, 
'11001' : 116, '11002' : 117, '11010' : 118, '11011' : 119, 
'11012' : 120, '11020' : 121, '120011' : 122, '120012' : 123, 
'120020' : 124, '120021' : 125, '120022' : 126, '120100' : 127, 
'120101' : 128, '120102' : 129, '120110' : 130, '120111' : 131, 
'120112' : 132, '120120' : 133, '120121' : 134, '120122' : 135, 
'120200' : 136, '120201' : 137, '120202' : 138, '120210' : 139, 
'120211' : 140, '120212' : 141, '120220' : 142, '120221' : 143, 
'120222' : 144, '121000' : 145, '121001' : 146, '121002' : 147, 
'121010' : 148, '121011' : 149, '121012' : 150, '121020' : 151, 
'121021' : 152, '121022' : 153, '121100' : 154, '121101' : 155, 
'121102' : 156, '121110' : 157, '121111' : 158, '121112' : 159, 
'121120' : 160, '121121' : 161, '121122' : 162, '121200' : 163, 
'121201' : 164, '121202' : 165, '121210' : 166, '121211' : 167, 
'121212' : 168, '121220' : 169, '121221' : 170, '121222' : 171, 
'122000' : 172, '122001' : 173, '122002' : 174, '122010' : 175, 
'122011' : 176, '122012' : 177, '122020' : 178, '122021' : 179, 
'122022' : 180, '122100' : 181, '122101' : 182, '122102' : 183, 
'122110' : 184, '122111' : 185, '122112' : 186, '122120' : 187, 
'122121' : 188, '122122' : 189, '122200' : 190, '122201' : 191, 
'122202' : 192, '122210' : 193, '122211' : 194, '122212' : 195, 
'122220' : 196, '122221' : 197, '122222' : 198, '00000' : 199, 
'00001' : 200, '00002' : 201, '00010' : 202, '00011' : 203, 
'00012' : 204, '00020' : 205, '00021' : 206, '00022' : 207, 
'00100' : 208, '00101' : 209, '00102' : 210, '00110' : 211, 
'00111' : 212, '00112' : 213, '00120' : 214, '00121' : 215, 
'00122' : 216, '00200' : 217, '00201' : 218, '00202' : 219, 
'00210' : 220, '00211' : 221, '00212' : 222, '00220' : 223, 
'00221' : 224, '00222' : 225, '01000' : 226, '01001' : 227, 
'01002' : 228, '01010' : 229, '01011' : 230, '01012' : 231, 
'01020' : 232, '01021' : 233, '01022' : 234, '01100' : 235, 
'01101' : 236, '01102' : 237, '01110' : 238, '01111' : 239, 
'01112' : 240, '01120' : 241, '01121' : 242, '01122' : 243, 
'01200' : 244, '01201' : 245, '01202' : 246, '01210' : 247, 
'01211' : 248, '01212' : 249, '01220' : 250, '01221' : 251, 
'01222' : 252, '02000' : 253, '02001' : 254, '02002' : 255
}

def huffman_encode_byte(byte):        
    e = huffman_enc_table[byte]
    e = [ _huffman_enc_map[x] for x in e ]
    return "".join(e)

def huffman_decode_check(s):
    d = [ _huffman_dec_map[x] for x in s ]
    d = "".join(d)
    return huffman_dec_table.has_key(d)

def huffman_decode_byte(s):
    d = [ _huffman_dec_map[x] for x in s ]
    d = "".join(d)
    assert huffman_dec_table.has_key(d)==True        
    return huffman_dec_table[d]

def huffman_encode(packet):
    array = bytearray(packet)
    strand = []
    for b in array:
        strand.append( huffman_encode_byte(b) )    
    return "".join(strand)

def huffman_decode(strand):
    i = 0
    array = []
    while i < len(strand):
        s = strand[i:i+5]
        if huffman_decode_check(s):        
            byte = huffman_decode_byte(s)
            array.append(chr(byte))
            i+=5
        else:
            s = strand[i:i+6]
            if huffman_decode_check(s):        
                byte = huffman_decode_byte(s)
                array.append(chr(byte))
                i+=6
            else:
                array.append(None)
            
    packet = bytearray(array)
    return packet

def huffman_decode_limited(strand,numBytes):
    i = 0
    count = 0
    array = []
    while i < len(strand) and count < numBytes:
        s = strand[i:i+5]
        if huffman_decode_check(s):        
            byte = huffman_decode_byte(s)
            count += 1
            array.append(chr(byte))
            i+=5
        else:
            s = strand[i:i+6]
            if huffman_decode_check(s):        
                byte = huffman_decode_byte(s)
                array.append(chr(byte))
                count += 1
                i+=6
            else:
                #append None in case I cant figure out what the data is due to an error in the strand
                array.append(chr(0))
                i+=6
                count+=1
    while count < numBytes:
        array.append(chr(0))
        count+=1
        
    packet = bytearray(array)
    return packet

rotate_map = {
    'A' : 0,
    'C' : 1,
    'G' : 2,
    'T' : 3
}
rotate_table = { 
    'A' : [ 'C', 'G', 'T'],
    'C' : [ 'G', 'T', 'A'],
    'G' : [ 'T', 'A', 'C'],
    'T' : [ 'A', 'C', 'G']
}
unrotate_table = { 
    'A' : ['X','A','C','G'],
    'C' : ['G','X','A','C'],
    'G' : ['C','G','X','A'],
    'T' : ['A','C','G','X']
}

def rotate_encode(strand,prev='A'):
    r = []
    for s in strand:
        if s == 'T':
            print strand
        assert s != 'T'
        n = rotate_table[prev][rotate_map[s]]
        r.append( n )
        prev = n
    return "".join(r)

def rotate_decode(strand,prev='A'):
    r = []
    for s in strand:
        n = unrotate_table[prev][rotate_map[s]]
        #add filler if X happens, means there was an error
        if n =='X':
            n='A'
        r.append( n )
        prev = s
    return "".join(r)

class HuffmanCodec(BaseCodec):
    def __init__(self,numberBytes,CodecObj=None,keyWidth=20,valueWidth=140):
        BaseCodec.__init__(self,CodecObj)
        self._keyWidth=keyWidth
        self._valueWidth=valueWidth
        self._numberBytes = numberBytes

    def _encode(self,packet):
        assert len(packet[1])<=self._numberBytes
        key = base_conversion.convertBase(3,packet[0],self._keyWidth)
        value = huffman_encode(packet[1])
        assert len(value) <= self._valueWidth
        if len(value) < self._valueWidth:
            value = value + base_conversion.randomTernary(self._valueWidth-len(value))
        assert len(value)==self._valueWidth
        return key+value

    def _decode(self,s):
        key = base_conversion.convertFromBase(3,s[0:self._keyWidth])
        value = huffman_decode_limited(s[self._keyWidth:],self._numberBytes)
        return key,value

class RotateCodec(BaseCodec):
    def __init__(self,CodecObj=None,prev='A'):
        BaseCodec.__init__(self,CodecObj)
        self._prev = prev

    def _encode(self,s):
        assert isinstance(s,str) or "Didn't get expected type."
        return rotate_encode(s,self._prev)

    def _decode(self,s):
        assert isinstance(s,str) or "Didn't get expected type."
        return rotate_decode(s,self._prev)

if __name__ == "__main__":
    print huffman_encode_byte(255)
    s = "12345678"
    print huffman_encode(s)            
    assert s == huffman_decode(huffman_encode(s))

    s = 'AGCAGCAGC'
    print s
    print rotate_encode(s)
    print rotate_decode(rotate_encode(s))
    assert s == rotate_decode(rotate_encode(s))
