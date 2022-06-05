 #!/usr/bin/python
from copy import copy
from random import *

bases = ['A', 'C', 'G', 'T']

values = { 'A' : 0 ,
           'C' : 1 ,
           'G' : 2 ,
           'T' : 3   }

def randomTernary(length):
    return "".join([ bases[x%3] for x in range(length) ])

def convertTernaryHelper(dec,s):
    m = dec % 3
    q = dec // 3
    s = s + bases[m]
    if q > 0:
        return convertTernaryHelper(q,s)
    else:
        return s

def convertTernary(dec,length):
    s = convertTernaryHelper(dec,'')
    s = s.ljust(length,bases[0])
    return s

def convertQuarnaryHelper(dec,s):
    m = dec % 4
    q = dec // 4
    s = s + bases[m]
    if q > 0:
        return convertQuarnaryHelper(q,s)
    else:
        return s

def convertQuarnary(dec,length):
    s = convertQuarnaryHelper(dec,'')
    s = s.ljust(length,bases[0])
    return s

def convertBaseHelper(base : int, dec : int, s : str, symbols = bases):
    m = int(dec % base)
    q = int(dec / base)
    #print(s)
    s = s + symbols[m]
    #print(s)
    if q > 0:
        return convertBaseHelper(base,q,s,symbols)
    else:
        return s
        
def convertBytetoBinary(x,length):
    s=convertBytetoBinaryHelper(x,'')
    s=s.ljust(length,"0")
    #print s
    return s

def convertBytetoBinaryHelper(x,s): #convert byte x to a binary string
    binary=['0','1']
    m=int(x%2)
    q=x//2
    s=s+binary[m]
    if q >0:
        return convertBytetoBinaryHelper(q,s)
    else:
        return s


def convertBase(base, dec, length, symbols=bases):
    assert base <= len(symbols)
    s = convertBaseHelper(base,dec,'', bases)
    s = s.ljust(length,bases[0])
    return s

def convertToAnyBase(base : int, dec : int, length: int, symbols=bases):
    assert base == len(symbols)
    s = convertBaseHelper(base,dec,'', symbols)
    while len(s) < length:
        s += symbols[0]
    assert length == len(s)
    return s

def convertFromBase(base,s):
    val = 0
    power = 1
    for i in range(len(s)):
        val += values[s[i]]*power
        power *= base
    return val

def convertIntToBytes(val,num_bytes):
    val = int(val)
    if val == None:
        return [-1 for _ in range(num_bytes)]
    else:
        #print ("HERE!",val,num_bytes,range(num_bytes))
        l = [(val & (0xff << pos*8)) >> pos*8 for pos in range(num_bytes)]
        return l

def convertBytesToInt(l):
    _sum = 0
    for i,val in enumerate(l):
        _sum += val * (256**i)
    return _sum

ibases = ['A', 'C', 'T']

ivalues = { 'A' : 0 ,
           'C' : 1 ,
           'G' : 3 ,
           'T' : 2   }

def encodeWithExclusionHelper(dec,s,excluded):
    if len(s) < len(excluded):
        b = len(excluded[len(s)])
        m = dec % b
        q = dec / b
        s += excluded[len(s)][m]
    else:
        m = dec % 3
        q = dec / 3
        s += ibases[m]

    if q > 0 or len(s) < len(excluded):
        return encodeWithExclusionHelper(q,s,excluded)
    else:
        return s

def encodeWithExclusion(dec,length,primer):
    excluded = []
    plist = [b for b in primer]
    for i in range(len(primer)):
        b = copy(ibases)
        if plist[i] != 'G':
            b.remove(plist[i])
        excluded.append(b)
    #print excluded
    s = encodeWithExclusionHelper(dec,"",excluded)
    if len(s) < len(primer):
        for i in range(len(s),min(length,len(primer))):
            s += excluded[i][0]
    s = s.ljust(length,ibases[0])
    return s


def decodeWithExclusion(s,primer):
    val = 0
    power = 1
    excluded = []
    plist = [b for b in primer]
    for i in range(len(primer)):
        b = copy(ibases)
        if plist[i]!='G':
            b.remove(plist[i])
        excluded.append( { b[i] : i for i in range(len(b)) }  )
    #print excluded
    base = 3
    for i in range(len(s)):
        if i < len(excluded):
            base = len(excluded[i])
            val += excluded[i][s[i]] * power
        else:
            base = 3
            val += ivalues[s[i]] * power
        power *= base
    return val


#packs index bits into an array of byte values suitable for encoding
def pack_bits_to_bytes(index_set,index_bit_sizes):
    byte_array=[]
    bits_in_remaining_byte=8
    current_byte=0
    for t, index in enumerate(index_set):
        index_size=index_bit_sizes[t]
        index_value=index

        bits_to_pack=index_size

        while bits_to_pack>0:
            if bits_to_pack > bits_in_remaining_byte: 
                bits_from_current_iter = bits_in_remaining_byte
            else:
                bits_from_current_iter = bits_to_pack
            current_byte = (current_byte << bits_from_current_iter)|(index_value>>(bits_to_pack-bits_from_current_iter))
            if bits_from_current_iter == bits_in_remaining_byte:
                byte_array.append(current_byte)
                current_byte=0
                bits_in_remaining_byte = 8
            else:
                bits_in_remaining_byte -= bits_from_current_iter
            mask= 2**(bits_to_pack-bits_from_current_iter)-1
            index_value = index_value & mask
            bits_to_pack-=bits_from_current_iter
    
        assert index_value==0 and bits_to_pack==0

    #put the last byte into the array if last byte is not completely filled, push all bits to the top of the byte
    if bits_in_remaining_byte >0 and bits_in_remaining_byte!=8:
        current_byte=current_byte<<bits_in_remaining_byte
        byte_array.append(current_byte)
    
    return byte_array
        
#unpacks bytes into the correct index values
def unpack_bytes_to_indexes(index_bytes,index_bit_sizes):
    index_array=[]
    byte_index=0
    current_byte = index_bytes[byte_index]
    bits_remaining_in_byte = 8 #track how many bits we've taken from the current byte
    for index_size in index_bit_sizes:
        bits_to_get = index_size
        current_index_value=0
        while bits_to_get>0:
            if bits_to_get > bits_remaining_in_byte:
                current_iter_bits = bits_remaining_in_byte
            else:
                current_iter_bits = bits_to_get
            current_index_value = (current_index_value<<current_iter_bits)|(current_byte>>(8-current_iter_bits))
            current_byte = (current_byte<<current_iter_bits)&((2**8)-1)
            bits_remaining_in_byte -=current_iter_bits
            bits_to_get -=current_iter_bits
            if bits_remaining_in_byte==0 and bits_to_get>0:
                byte_index+=1
                current_byte=index_bytes[byte_index]
                bits_remaining_in_byte=8
        index_array.append(current_index_value)

    return index_array


if __name__ == "__main__":
    import math
    primer = "GTCTCGTGGGCTCGG"
    #for i in range(3**5):
    #    s = encodeWithExclusion(i,10,primer)
    #    print "{} = {}  ({})".format(i,s,decodeWithExclusion(s,primer))
    #    assert i == decodeWithExclusion(s,primer)

    for i in range(2**8):
        print (convertBase(2,i,8))

    for i in range (len(primer),60):
        val = (2**len(primer)) * (3 **(i-len(primer)))
        print ("Max({}) = {}".format(i,math.log(val,2)))

    for i in [2**64, 2**74, 2**80]:
        s = encodeWithExclusion(i,100,primer)
        print (s)


    #test out the packing/unpacking of indexing
    
    index_array = [ random.randint(0,255) for i in range(0,20)] #generate a list of random indices

    index_size_array = [ random.randint(math.ceil(math.log(i,2)),8) for i in index_array] #choose random sizes that will at least satisfy each index

    byte_array = pack_bits_to_bytes(index_array,index_size_array)

    index_array2 = unpack_bytes_to_indexes(byte_array,index_size_array)

    print(index_array)
    print(index_array2)
    
    assert index_array==index_array2

    

