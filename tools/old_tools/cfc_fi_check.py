from dnastorage.codec.commafreecodec import *
from dnastorage.primer.primer_util import *
from dnastorage.primer.design import *
from random import randint,shuffle
import sys
import generate


def substitution(sims=10000,rate=0.01):    
    codec = CommaFreeCodewords(10)

    match = 0
    mismatch = 0

    fails = 0
    should_fail = 0
    for _ in range(sims):
        arr = [ randint(0,255) for _ in range(10) ]
        cf_arr = codec.encode(arr)
        cf_arr = "".join(cf_arr)

        err = [ generate.rand() < rate for _ in range(len(cf_arr)) ]

        list_cf = [ x for x in cf_arr ]
        err_ins =0
        for i,e in enumerate(err):
            if e:
                err_ins += 1
                prior = list_cf[i]
                while list_cf[i] == prior:
                    list_cf[i] = random.choice('AGCT')

        err_cf = "".join(list_cf)

        assert len(err_cf) == len(cf_arr)

        print "+"*50
        print cf_arr

        s = ""
        for a,b in zip(cf_arr,err_cf):
            if a!=b:
                s += "|"
            else:
                s += " "
        print s
        print err_cf

        new_arr = codec.decode(err_cf)

        print arr
        print new_arr
        print "-"*50

        err = 0
        for a,b in zip(arr,new_arr):
            if a!=b:
                err += 1
                mismatch += 1
            else:
                match += 1

        if err_ins > 0:
            should_fail += 1
        if err > 0:
            fails += 1

    print "substition failed     : {:.2f}".format( fails / float(sims) * 100.0 )
    print "substition should fail:{:.2f}".format( should_fail / float(sims) * 100.0 )
    print "substition match rate : {:.2f}".format( match / (match + mismatch + 0.0) * 100.0 )

def deletion(sims=10000,rate=0.01):    
    codec = CommaFreeCodewords(10)
    
    match = 0
    mismatch = 0
    fails = 0
    should_fail = 0
    for _ in range(sims):
        arr = [ randint(0,255) for _ in range(10) ]
        cf_arr = codec.encode(arr)
        cf_arr = "".join(cf_arr)

        err = [ generate.rand() < rate for _ in range(len(cf_arr)) ]

        list_cf = [ x for x in cf_arr ]
        err_ins =0
        indices = []
        for i,e in enumerate(err):
            if e:
                err_ins += 1
                indices += [i]
        indices.reverse()

        for i in indices:
            err_ins += 1
            list_cf = list_cf[:i]+list_cf[i+1:]
                
        err_cf = "".join(list_cf)

        #assert len(err_cf) == len(cf_arr)

        print "+"*50
        print cf_arr

        s = ""
        for a,b in zip(cf_arr,err_cf):
            if a!=b:
                s += "|"
            else:
                s += " "
        print s
        print err_cf

        new_arr = codec.decode(err_cf)

        print arr
        print new_arr
        print "-"*50


        assert len(arr)==len(new_arr)
        
        err = 0
        for a,b in zip(arr,new_arr):
            if a!=b:
                mismatch += 1
                err += 1
            else:
                match += 1

        if err_ins > 0:
            should_fail += 1
        if err > 0:
            fails += 1

    print "deletion failed     : {:.2f}".format( fails / float(sims) * 100.0 )
    print "deletion should fail:{:.2f}".format( should_fail / float(sims) * 100.0 )

    print "deletion match rate : {:.2f}".format( match / (match + mismatch + 0.0) * 100.0 )
    
    #print "{:.2f}".format( fails / 10000.0 * 100.0 )
    #print "{:.2f}".format( should_fail / 10000.0 * 100.0 )

def insertion(sims=10000,rate=0.01):    
    codec = CommaFreeCodewords(10)
    
    match = 0
    mismatch = 0
    fails = 0
    should_fail = 0
    for _ in range(sims):
        arr = [ randint(0,255) for _ in range(10) ]
        cf_arr = codec.encode(arr)
        cf_arr = "".join(cf_arr)

        err = [ generate.rand() < rate for _ in range(len(cf_arr)) ]

        list_cf = [ x for x in cf_arr ]
        err_ins =0
        indices = []
        for i,e in enumerate(err):
            if e:
                err_ins += 1
                indices += [i]
        indices.reverse()

        for i in indices:
            err_ins += 1
            list_cf = list_cf[:i]+[random.choice('AGCT')]+list_cf[i:]
                
        err_cf = "".join(list_cf)

        #assert len(err_cf) == len(cf_arr)

        if err_ins > 0:
            print "+"*50
            print cf_arr

            s = ""
            for a,b in zip(cf_arr,err_cf):
                if a!=b:
                    s += "|"
                else:
                    s += " "
            print s
            print err_cf

            new_arr = codec.decode(err_cf)

            print arr
            print new_arr
            print "-"*50


        assert len(arr)==len(new_arr)
        
        err = 0
        for a,b in zip(arr,new_arr):
            if a!=b:
                mismatch += 1
                err += 1
            else:
                match += 1

        if err_ins > 0:
            should_fail += 1
        if err > 0:
            fails += 1

    print "deletion failed     : {:.2f}".format( fails / float(sims) * 100.0 )
    print "deletion should fail:{:.2f}".format( should_fail / float(sims) * 100.0 )
    print "deletion match rate : {:.2f}".format( match / (match + mismatch + 0.0) * 100.0 )
    
    #print "{:.2f}".format( fails / 10000.0 * 100.0 )
    #print "{:.2f}".format( should_fail / 10000.0 * 100.0 )


    

#substitution(1000)
#deletion(1000)
insertion(1000)
