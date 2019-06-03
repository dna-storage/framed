from PIL import Image,ImageFilter,ImageFont,ImageDraw
from numpy import mean
from scipy.fftpack import *
from numpy import *

from dnastorage.codec.huffman_table import *

import sys

Q = array([ [16, 11, 10, 16, 24, 40, 51, 61],
      [12, 12, 14, 19, 26, 58, 60, 55],
      [14, 13, 16, 24, 40, 57, 69, 56],
      [14, 17, 22, 29, 51, 87, 80, 62],
      [18, 22, 37, 56, 68, 109, 103, 77],
      [24, 35, 55, 64, 81, 104, 113, 92],
      [49, 64, 78, 87, 103, 121, 120, 101],
      [72, 92, 95, 98, 112, 100, 103, 99] ])

Q = ones( (8,8) )


order =   [ (0,0) ] \
        + [ (i, 1-i) for i in range(2) ] \
        + [ (2-i, i) for i in range(3) ] \
        + [ (i, 3-i) for i in range(4) ] \
        + [ (4-i, i) for i in range(5) ] \
        + [ (i, 5-i) for i in range(6) ] \
        + [ (6-i, i) for i in range(7) ] \
        + [ (i, 7-i) for i in range(8) ] \
        + [ (1+i, 7-i) for i in range(7) ]\
        + [ (7-i, i+2) for i in range(6) ]\
        + [ (i+3, 7-i) for i in range(5) ]\
        + [ (7-i, i+4) for i in range(4) ]\
        + [ (i+5, 7-i) for i in range(3) ]\
        + [ (7-i, i+6) for i in range(2) ]\
        + [ (i+7, 7-i) for i in range(1) ]


# assume all images are same size
def make_collage(imlist,labels=None,buff=50):
    x = max([ len(i) for i in imlist ])
    y = len(imlist)
    all_size = ((imlist[0][0].size[0]),(imlist[0][0].size[1]))
    size = ((imlist[0][0].size[0]+buff)*x,(imlist[0][0].size[1]+buff)*y)
    imstr = "".join( [ chr(255) for x in range(size[0]*size[1]*3) ] )
    image = Image.frombuffer('RGB',size,imstr,"raw",'RGB',0,1) #.convert('YCbCr')
    draw = ImageDraw.Draw(image)
    #font = ImageFont.truetype("sans-serif.ttf", 12)
    i = 0
    j = 0
    if labels==None:
        for il in imlist:
            for im in il:
                image.paste(im,(i,j,i+im.size[0],j+im.size[1]))
                draw.text( (i+im.size[0]/2,j+im.size[1]+15), "Test", (0,0,0))
                i += all_size[0]+buff
            j += all_size[1]+buff
            i = 0
    else:
        for il,lab in zip(imlist,labels):
            for im,la in zip(il,lab):
                image.paste(im,(i,j,i+im.size[0],j+im.size[1]))
                draw.text( (i+im.size[0]/2,j+im.size[1]+15),la, (0,0,0))
                i += all_size[0]+buff
            j += all_size[1]+buff
            i = 0

    return image


def getpixellist(im):
    out = []
    data = im.getdata()
    for y in range(0,im.size[1]):
        for x in range(0,im.size[0]):
            p = data.getpixel( (x,y) )
            out.append(p[0])
            out.append(p[1])
            out.append(p[2])                
    return out

def get2drange(im,tl,br):
    res = []
    data = im.getdata()
    for y in range(tl[1],br[1]):
        d = []
        for x in range(tl[0],br[0]):
            d.append(data.getpixel(x,y))
        res.append(d)
    return res


def getYpixels(image,x,yz):
    y = []
    #print x, yz, x+8, yz+8, image.size
    assert x % 8 == 0
    assert yz % 8 == 0
    assert x >= 0 and x+8 <= image.size[0] 
    assert yz >= 0 and yz+8 <= image.size[1]
    for j in range(yz,min(image.size[1],yz+8)):
        yy = []
        cbcb = []
        crcr = []
        for i in range(x,min(image.size[0],x+8)):
            p = image.getdata().getpixel( (i,j) )
            yy.append(p[0])
        while len(yy) < 8:
            yy.append( yy[-1] )
        y.append(yy)

    while len(y) < 8:
        y.append( y[-1] )
    return array(y)

def getCpixels(image,b,x,yz):
    y = []
    #print x, yz, x+8, yz+8, image.size
    assert b==1 or b==2
    #assert x % 8 == 0
    #assert yz % 8 == 0
    #assert x >= 0 and x+16 <= image.size[0] 
    #assert yz >= 0 and yz+16 <= image.size[1]
    for j in range(yz,min(image.size[1],yz+16),2):
        yy = []
        cbcb = []
        crcr = []
        for i in range(x,min(image.size[0],x+16),2):
            p = image.getdata().getpixel( (i,j) )
            yy.append(p[b])
        while len(yy) < 8:
            yy.append( yy[-1] )
        y.append(yy)
    while len(y) < 8:
        y.append( y[-1] )
    return array(y)

def setCpixels(image,b,x,yz,Y):
    rate = 8
    assert x % 8 == 0
    assert yz % 8 == 0
    if not( x >= 0 and x+8 <= image.size[0] ):
        return
    if not( yz >= 0 and yz+8 <= image.size[1] ):
        return
    for j in range(yz,min(image.size[1],yz+rate)):
        for i in range(x,min(image.size[0],x+rate)):
            p = image.getdata().getpixel( (i,j) )
            if b == 1:
                p = ( p[0] , int(Y[j/2][i/2]) , p[2] )
            if b == 2:
                p = ( p[0] , p[1], int(Y[j/2][i/2]) )
            image.putpixel( (i,j), p)
    return

def setYpixels(image,x,yz,Y):
    rate = 8
    assert x % 8 == 0
    assert yz % 8 == 0
    if not( x >= 0 and x+8 <= image.size[0] ):
        return
    if not( yz >= 0 and yz+8 <= image.size[1] ):
        return
    for j in range(yz,min(image.size[1],yz+rate)):
        for i in range(x,min(image.size[0],x+rate)):
            p = image.getdata().getpixel( (i,j) )
            p = ( int(Y[j][i]) , p[1] , p[2] )
            #p = ( int(Y[j][i]) , int(Y[j][i]) , int(Y[j][i]) )
            #print p
            image.putpixel( (i,j), p)
    return


def setYpixels2(image,x,yz,Y):
    rate = 8
    assert x % 8 == 0
    assert yz % 8 == 0
    if not( x >= 0 and x+8 <= image.size[0] ):
        return
    if not( yz >= 0 and yz+8 <= image.size[1] ):
        return
    for j in range(yz,min(image.size[1],yz+rate)):
        for i in range(x,min(image.size[0],x+rate)):
            p = image.getdata().getpixel( (i,j) )
            #p = ( int(Y[j][i]) , p[1] , p[2] )
            p = ( int(Y[j][i]), int(Y[j][i]), int(Y[j][i]) )
            #print p
            image.putpixel( (i,j), p)
    return


def dct2(y):
    global Q
    R = []
    for yy in y:
        # normalize vector to +-128
        a = array(yy) - 128
        # perform dct and divide by 16 (8*2)
        r = dct(a) / 16
        R.append(r)
    #print transpose(array(R))
    R = transpose(array(R))
    newR = []
    for r in R:
        rr = dct(r) / 16
        newR.append(rr)
    R = transpose(array(newR))
    R = rint( R / Q )
    R = R.astype(int)
    return R

def idct2(y):
    global Q
    R = []
    y = y * Q
    y = transpose(array(y))
    for yy in y:
        # normalize vector to +-128
        a = array(yy)
        # perform dct and divide by 16 (8*2)
        r = idct(a) 
        R.append(r)
    #print transpose(array(R))
    R = transpose(array(R))
    newR = []
    for r in R:
        rr = idct(r) 
        newR.append(rr)
    #R = transpose(array(newR))
    R = rint( array(newR) )
    R = array(R)+128
    return R


def hack_box_Y(image,gray,x,y,keep):
    box = image.crop((x,y,x+8,y+8))
    Y = getYpixels(box,0,0)
    Y_dct = dct2(Y)        
    z = 0
    sz = 0
    for i,o in enumerate(order):
        p = Y_dct[ o[0] ][ o[1] ]
        if sz+1 > keep:
            Y_dct[ o[0] ][ o[1] ] = 0
        if p==0:
            if z==15:
                sz += 1
                z = 0
            z += 1
        else:            
            sz += 1.5
            z = 0
    Y_idct = idct2(Y_dct)
    calc_size[0] += sz
    calc_size[1] += min(sz,keep)
    if gray!=None:
        gbox = gray.crop((x,y,x+8,y+8))
        setYpixels2(gbox,0,0,Y_idct)
        gray.paste(gbox,(x,y,x+8,y+8))
    else:
        setYpixels(box,0,0,Y_idct)
        image.paste(box,(x,y,x+8,y+8))

def hack_box_C(image,band,x,y,keep,gray=False):
    box = image.crop((x,y,x+16,y+16))
    Y = getCpixels(box,band,0,0)
    Y_dct = dct2(Y)        
    z = 0
    sz = 1
    for i,o in enumerate(order):
        p = Y_dct[ o[0] ][ o[1] ]
        if sz > keep:
            Y_dct[ o[0] ][ o[1] ] = 0
        if p==0:
            if z==15:
                sz += 1
                z = 0
            z += 1
        else:            
            sz += 1.5
            z = 0
    Y_idct = idct2(Y_dct)

    global calc_size
    calc_size[0] += sz
    if gray==False:
        calc_size[1] += min(sz,keep)
    #print (keep / float(sz) * 100)
    setCpixels(box,band,0,0,Y_idct)
    image.paste(box,(x,y,x+16,y+16))


def hack_image(filename,size,Y, C, gray=False):
    global calc_size
    calc_size = [0,0]
    im = Image.open(filename).convert('YCbCr')
    im.thumbnail( size )
    if gray==True:
        grayim = Image.open(filename)
        grayim.thumbnail(size)
    else:
        grayim = None
    C_strand = C
    Y_strand = Y
    #size = (im.size[0]/8*8,im.size[1]/8*8)
    for z in range(0,size[1],8):
        for x in range(0,size[0],8):
            hack_box_Y(im,grayim,x,z,Y_strand)
    for z in range(0,size[1],16):
        for x in range(0,size[0],16):
            hack_box_C(im,1,x,z,C_strand,gray)
            hack_box_C(im,2,x,z,C_strand,gray)

    if gray==False:
        return im
    else:
        #grayim.show()
        return grayim
    

if len(sys.argv)==1:
    filename = 'obama.jpg'
    f = 10
elif len(sys.argv) > 1:
    filename = sys.argv[1]
    f = 5
    if len(sys.argv) > 2:
        f = int(sys.argv[2])


imt = Image.open(filename).convert('YCbCr')
imt.thumbnail( (imt.size[0]/f, imt.size[1]/f) )
#imt.show()


calc_size = [0,0]

#tmp = [ imt ]
images = []
labels = []

for C in [100,5,1]:
    tmp = [ ]
    labe = []
    for Y in [100,20,5,1]:
        tmp.append( hack_image(filename,imt.size,Y,C) )
        labe.append( "{:.1f}KB,{:.0f}%,C={},Y={}".format(calc_size[0]/1000.,calc_size[1]/float(calc_size[0])*100,C,Y) )
    images.append(tmp)
    labels.append(labe)

tmp = []
labe = []
for Y in [100,20,5,1]:
    tmp.append( hack_image(filename,imt.size,Y,C,True) )
    labe.append( "{:.1f}KB,{:.0f}%,Gray,Y={}".format(calc_size[0]/1000.,calc_size[1]/float(calc_size[0])*100,Y) )
images.append(tmp)
labels.append(labe)

im = make_collage(images,labels)
im.show()
im.save("grid.jpg")

sys.exit(0)


im = Image.open(filename).convert('YCbCr')
im.thumbnail( (im.size[0]/f, im.size[1]/f) )

grayimage = Image.open(filename)
grayimage.thumbnail( (grayimage.size[0]/f, grayimage.size[1]/f) )

data = im.getdata()
size = (im.size[0]/8*8,im.size[1]/8*8)

rate = 8

C_strand = 2
Y_strand = 10

for z in range(0,size[1],rate):
    for x in range(0,size[0],rate):
        hack_box_Y(im,grayimage,x,z,Y_strand)
        #pass

for z in range(0,size[1],16):
    for x in range(0,size[0],16):
        #pass
        hack_box_C(im,1,x,z,C_strand)
        hack_box_C(im,2,x,z,C_strand)


#im.show()
#grayimage.show()

c = make_collage( [ [ imt, im, grayimage ] ] )
c.show()

sys.exit(0)

y = getYpixels(im,128,64)
y_dct = dct2(y)
y_prime = idct2(y_dct)

diff = array(y) - array(y_prime)

print y
print y_dct
print y_prime
print diff

sys.exit(0)

Y = []
Cb = []
Cr = []

for z in range(0,size[1],rate):
    for x in range(0,size[0],rate):
        y = []
        cb = []
        cr = []        
        for j in range(z,min(size[1],z+rate)):
            yy = []
            cbcb = []
            crcr = []
            for i in range(x,min(size[0],x+rate)):
                p = data.getpixel( (i,j) )
                yy.append(p[0])
                if i%2 == 0:
                    cbcb.append(p[1])
                    crcr.append(p[2])                
            y.append(yy)
            if j%2==0:
                cb.append(cbcb)
                cr.append(crcr)
        Y.append(y)
        Cb.append(cb)
        Cr.append(cr)


        
print order

f = {}

Im = []

for i,y in enumerate(Y):
    R = []
    for yy in y:
        a = array(yy) - 128
        r = dct(a)/16
        R.append(r)
    #print transpose(array(R))
    R = transpose(array(R))
    newR = []
    for r in R:
        rr = dct(r)/16
        newR.append(rr)
    R = transpose(array(newR))
    R = rint( R / Q )
    R = R.astype(int)
    #print i, R
    l = []
    for p,o in enumerate(order):
        r = R[ o[0] ][ o[1] ]
        l.append(r)
        if r != 0:
            f[r] = f.get(r,0) + 1        
    #print i, l

    Im.append( [i,l] )

print len(f.keys())

w = []
sym = []
for k,v in f.items():
    sym.append(k)
    w.append(v)

w = list(array(w)/float(sum(w)))

ht = HuffmanTable(2,['0', '1'],sym, w)
enc,dec = ht.get_tables()

Im2 = []
for i in Im:
    index = i[0]
    l = i[1]
    newl = []
    z = 0
    Len = 0
    for j,ll in enumerate(l):
        if ll != 0:
            newl.append( [z, len(enc[ll]), enc[ll]] )
            Len += 8 + len(enc[ll])
            z = 0
        else:
            if z==15:
                newl.append( [15,0,0] )
                z = 0
                Len += 8
            z += 1
    Im2.append(newl)
    print "{} - {}".format(ceil(Len/8.0),newl)




print ht.average_length()




sys.exit(0)

rate = 10

size = (im.size[0]/rate*rate,im.size[1]/rate*rate)

out = []

for y in range(0,size[1],rate):
    for x in range(0,size[0],rate):
        r = []
        g = []
        b = []
        for i in range(x,min(size[0],x+rate)):
            for j in range(y,min(size[1],y+rate)):
                p = data.getpixel( (i,j) )
                r.append(p[0])
                g.append(p[1])
                b.append(p[2])                
        res = [ int(mean(r)), int(mean(g)), int(mean(b)) ]
        out += res


imstr = "".join( [ chr(x) for x in out ] )

msize = (size[0]/rate,size[1]/rate)

im2 = Image.frombuffer('RGB',msize,imstr,"raw",'RGB',0,1)
im2.show()


imt = Image.open('obama.jpg')
imt.thumbnail( msize )

imt.show()
dt = imt.getdata()

dtout = getpixellist(imt)

print im2.size, imt.size

print len(out), len(dtout)
