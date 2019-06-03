import numpy as np
from PIL import Image
from scipy.fftpack import dct,idct
from skimage import measure
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import math

Q_Factor=1


#luminance quantization matrix from JPEG standard document
Q_luminance=[[16,11,10,16,24,40,51,61],[12,12,14,19,26,58,60,55],[14,13,16,24,40,57,69,56],
             [14,17,22,29,51,87,80,62],[18,22,37,56,68,109,103,77],[24,35,55,64,81,104,113,92],
             [49,64,78,87,103,121,120,101],[72,92,95,98,112,100,103,99]]


#AC lum nance table taken from the JPEG specification
AC_lum_table={
    "0:1":3,"0:2":4,"0:3":6,"0:4":8,"0:5":10,"0:6":12,"0:7":14,"0:8":18,"0:9":25,"0:10":26,
    "1:1":4,"1:2":5,"1:3":7,"1:4":9,"1:5":11,"1:6":16,"1:7":16,"1:8":16,"1:9":16,"1:10":16,
    "2:1":5,"2:2":8,"2:3":10,"2:4":12,"2:5":16,"2:6":16,"2:7":16,"2:8":16,"2:9":16,"2:10":16,
    "3:1":6,"3:2":9,"3:3":12,"3:4":16,"3:5":16,"3:6":16,"3:7":16,"3:8":16,"3:9":16,"3:10":16,
    "4:1":6,"4:2":10,"4:3":16,"4:4":16,"4:5":16,"4:6":16,"4:7":16,"4:8":16,"4:9":16,"4:10":16,
    "5:1":7,"5:2":11,"5:3":16,"5:4":16,"5:5":16,"5:6":16,"5:7":16,"5:8":16,"5:9":16,"5:10":16,
    "6:1":7,"6:2":12,"6:3":16,"6:4":16,"6:5":16,"6:6":16,"6:7":16,"6:8":16,"6:9":16,"6:10":16,
    "7:1":8,"7:2":12,"7:3":16,"7:4":16,"7:5":16,"7:6":16,"7:7":16,"7:8":16,"7:9":16,"7:10":16,
    "8:1":9,"8:2":15,"8:3":16,"8:4":16,"8:5":16,"8:6":16,"8:7":16,"8:8":16,"8:9":16,"8:10":16,
    "9:1":9,"9:2":16,"9:3":16,"9:4":16,"9:5":16,"9:6":16,"9:7":16,"9:8":16,"9:9":16,"9:10":16,
    "10:1":9,"10:2":16,"10:3":16,"10:4":16,"10:5":16,"10:6":16,"10:7":16,"10:8":16,"10:9":16,"10:10":16,
    "11:1":10,"11:2":16,"11:3":16,"11:4":16,"11:5":16,"11:6":16,"11:7":16,"11:8":16,"11:9":16,"11:10":16,
    "12:1":10,"12:2":16,"12:3":16,"12:4":16,"12:5":16,"12:6":16,"12:7":16,"12:8":16,"12:9":16,"12:10":16,
    "13:1":11,"13:2":16,"13:3":16,"13:4":16,"13:5":16,"13:6":16,"13:7":16,"13:8":16,"13:9":16,"13:10":16,
    "14:1":16,"14:2":16,"14:3":16,"14:4":16,"14:5":16,"14:6":16,"14:7":16,"14:8":16,"14:9":16,"14:10":16,
    "15:1":16,"15:2":16,"15:3":16,"15:4":16,"15:5":16,"15:6":16,"15:7":16,"15:8":16,"15:9":16,"15:10":16
}

ZRL=10 #number of bits for a zero run (16 zeros in a row) symbol

    
#2D dct and idct
def dct2(a):
    return dct(dct(a,axis=0,norm='ortho'),axis=1,norm='ortho')

def idct2(a):
    return idct(idct(a,axis=0,norm='ortho'),axis=1,norm='ortho')

def quant(a):
    for i in range(0,8):
        for j in range(0,8):
            a[i,j]=round(a[i,j]/(float(Q_luminance[i][j])/Q_Factor))

def iquant(a):
    for i in range(0,8):
        for j in range(0,8):
            a[i,j]=a[i,j]*(float(Q_luminance[i][j])/Q_Factor)



#block_statistics calculates statistics on blocks of data from the image
def block_statistics(image_statistics_array, image_name):
    fig=plt.figure()
    ax1=fig.add_subplot(111,projection='3d')
    block_size_data=np.zeros(len(image_statistics_array),dtype=float) #total size in bits for each block in the image
    block_num_coef=np.zeros(len(image_statistics_array),dtype=float) #number of non-zero coeficients in each block
    total_image_size=0 #total size of image in bits
    total_image_coef=0 #total number of coefficients in the image
    block_data_distribution=np.zeros((8,8),dtype=float) #8x8 block that accumulates total number of bits for each (i,j) position 
    for block_idx, block_data in enumerate(image_statistics_array):
        block_size=0
        block_coef=0
        for data_idx, data_coef in enumerate(block_data):
            total_image_size+=data_coef["coef_bits"]
            total_image_coef+=1
            block_coef+=1
            block_size+=data_coef["coef_bits"]
            block_data_distribution[data_coef["i"]][data_coef["j"]]+=data_coef["coef_bits"]
        block_size_data[block_idx]+=block_size
        block_num_coef[block_idx]+=block_coef
    avg_coef_per_block=np.mean(block_num_coef)
    avg_bits_per_block=np.mean(block_size_data)
    avg_bits_per_coef=float(total_image_size)/float(total_image_coef)
    print "avg block size in bits: {}".format(avg_bits_per_block)
    print "avg number coef per block: {}".format(avg_coef_per_block)
    print "total image size in bits: {}".format(total_image_size)
    print "total coef in image: {}".format(total_image_coef)
    print "avg bits per coef: {}".format(avg_bits_per_coef)
    xpos=[]
    ypos=[]
    dx=[]
    dy=[]
    dz=[]
    cdf_data={}
    zpos=np.zeros(64)
    for i in range(0,8):
        for j in range(0,8):
            xpos.append(i)
            ypos.append(j)
            dx.append(0.5)
            dy.append(0.5)
            dz.append(block_data_distribution[i][j]/float(total_image_size))
            distance=math.sqrt(float(i)**2+float(j)**2)
            if distance not in cdf_data:
                cdf_data[distance]=block_data_distribution[i][j]/float(total_image_size)
            else:
                cdf_data[distance]=(block_data_distribution[i][j]/float(total_image_size))+cdf_data[distance]
    distance_array=[]
    for distance in cdf_data:
        distance_array.append(distance)
    distance_array.sort()
    cdf_array=np.zeros(len(distance_array))
    for idx, distance in enumerate(distance_array):
        if idx==0:
            cdf_array[idx]=cdf_data[distance]
        else:
            cdf_array[idx]=cdf_array[idx-1]+cdf_data[distance]

    ax1.bar3d(xpos,ypos,zpos,dx,dy,dz)
    plt.xlabel('i position')
    plt.ylabel('j position')
    plt.title('% of data at (i,j)')
    plt.savefig("spatial_amount_"+image_name+".jpg")
    fig2=plt.figure()
    ax2=fig2.add_subplot(111)
    ax2.plot(distance_array,cdf_array)
    plt.xlabel('distance from (i=0,j=0)')
    plt.ylabel('Cumulative % of data')
    
    plt.savefig("cdf_data_"+image_name+".jpg")
    print "distance_array {}".format(distance_array)
    print "cdf_array {}".format(cdf_array)
    


#Process the data obtained for each block by the process_pixel_block function, calculates the total number of bits required of each coefficient in the block
def huffman_block_data(image_statistics_array,category_table):
    for block_ID, block_data in enumerate(image_statistics_array):
        #huffman encode each pieve of data for the block and figure out the number of bits required for each element
        for data_idx, data_coef in enumerate(block_data):
            #print data_coef
            if data_idx==0 and len(data_coef)>0:
                if block_ID==0:
                    DC_difference=data_coef["val"]
                else:
                    #print "block_ID {} image_statistics_array[block_ID] {}".format(block_ID,data_coef)
                    #print "image_statistics_array[block_ID-1] {}".format(image_statistics_array[block_ID-1])
                    if len(image_statistics_array[block_ID-1])>0 and len(data_coef)>0:
                        DC_difference=data_coef["val"]-image_statistics_array[block_ID-1][0]["val"]
                    else:
                        if len(data_coef)>0:
                            DC_difference=data_coef["val"]
                #print "DC diff {}".format(abs(DC_difference))
                data_coef["category"]=category_table[abs(DC_difference)]["category"]
                #record the number of bits required to encode the DC component
                data_coef["coef_bits"]=category_table[abs(DC_difference)]["DC_huffman_length"]+category_table[abs(DC_difference)]["2c_length"]
            else:
                AC_val=data_coef["val"]
                #have an AC component, need to use the AC tables
                data_coef["category"]=category_table[abs(AC_val)]["category"]
                zero_run=data_coef["num_zeros"]
                run_key=zero_run%16
                AC_table_key=str(run_key)+":"+str(data_coef["category"])
                #print "AC table key {}".format(AC_table_key)
                #AC bits are the sum of the 2's complement value encoding, the AC_lum_table entry value, and the ZRL multiplied by the number of 16 runs before the coefficient
                data_coef["coef_bits"]=category_table[abs(AC_val)]["2c_length"]+AC_lum_table[AC_table_key]+(zero_run/16)*ZRL
            
        
# Process a block of pixels
# B: block descriptor
# image_matrix: matrix format of the image being processed
# returns an array with a dictionary for each non-zero piece of the quantized matrix
def process_pixel_block(B,image_matrix):
    h_l=B["height"][0]
    h_u=B["height"][1]
    w_l=B["width"][0]
    w_u=B["width"][1]
    block_slice=image_matrix[h_l:h_u,w_l:w_u]
    slice_matrix=np.zeros(shape=(8,8))
    block_statistics_array=[] #each element holds a non-zero component of the quantized matrix
    w=block_slice.shape[1]
    h=block_slice.shape[0]
    for i in range(0,h):
        for j in range(0,w):
            slice_matrix[i,j]=block_slice[i,j]-127 #down shift
    #print "start"
    #print slice_matrix
    slice_matrix=dct2(slice_matrix)
    #print slice_matrix
    quant(slice_matrix)
    #print slice_matrix
    #up to this point we have a quantized matrix, next step is to gather information about non zero elements
    i_1=0
    j_1=0
    i_2=0
    j_2=0
    zero_counter=0
    #ierate over the matrix diagonally
    for k in range(0,15):
        if k%2==0:
            ii=i_2
            jj=j_2
        else:
            ii=i_1
            jj=j_1
        _break=0
        #print "i_1 {} j_1 {} i_2 {} j_2 {}".format(i_1,j_1,i_2,j_2)
        while 1:
            #print "k {}".format(k)
            if slice_matrix[ii][jj]==0:
                zero_counter+=1
            else:
                block_statistics_array.append({})
                #collect information on the non-zero array value, needed for huffman encoding
                block_statistics_array[-1]["num_zeros"]=zero_counter
                block_statistics_array[-1]["i"]=ii
                block_statistics_array[-1]["j"]=jj
                block_statistics_array[-1]["val"]=int(slice_matrix[ii][jj])
                zero_counter=0 #clear the zero counter to start a new run
            if k%2==0:
                if ii>0:
                    ii-=1
                if jj<7:
                    jj+=1
                if _break==1 or (i_1==0 and j_1==0) or (i_1==7 and j_1==7):
                    break
                if (ii==i_1 and jj==j_1):
                    _break=1
            else:
                if ii<7:
                    ii+=1
                if jj>0:
                    jj-=1
                if _break==1:
                    break
                if (ii==i_2 and jj==j_2):
                    _break=1
        #need to fix up i_1,j_1,i_2,j_2
        if j_1==7:
            i_1+=1
        else:
            j_1+=1
        if i_2==7:
            j_2+=1
        else:
            i_2+=1
    #assert(len(block_statistics_array)>0)
    #print "finished zig zag"
    return block_statistics_array
    
def filter_image(B,image_statistics_array,filter_length,new_image_data):
    h_l=B["height"][0]
    h_u=B["height"][1]
    w_l=B["width"][0]
    w_u=B["width"][1]
    block_slice=new_image_data[h_l:h_u,w_l:w_u]
    w=block_slice.shape[1]
    h=block_slice.shape[0]
    slice_matrix=np.zeros(shape=(8,8))

    #filter out coefficients for the particular block
    for coef in image_statistics_array:
        #print coef
        distance=math.sqrt(float(coef["i"])**2+float(coef["j"])**2)
        if distance<filter_length:
            slice_matrix[coef["i"],coef["j"]]=float(coef["val"])
    iquant(slice_matrix)
    slice_matrix=idct2(slice_matrix)

    for i in range(0,h):
        for j in range(0,w):
            if slice_matrix[i,j]<128:
                slice_matrix[i,j]+=127
            else:
                slice_matrix[i,j]=255
            if slice_matrix[i,j]<0:
                slice_matrix[i,j]=0

    #print slice_matrix
    slice_matrix=slice_matrix.astype(np.uint8,copy=False)
    new_image_data[h_l:h_u,w_l:w_u]=slice_matrix[0:h,0:w]


#generate lookup table to determine a coefficient's category
def gen_cat_table():
    i=0
    cat_length=[2,3,3,3,3,3,4,5,6,7,8,9]
    cat_table={}
    cat_table[0]={}
    cat_table[1]={}
    cat_table[0]["category"]=0
    cat_table[0]["2c_length"]=0
    cat_table[0]["DC_huffman_length"]=2
    cat_table[1]["category"]=1
    cat_table[1]["2c_length"]=2
    cat_table[1]["DC_huffman_length"]=3
    i=2
    for cat_counter in range(2,12):
        ii=i
        while ii<i*2:
            cat_table[ii]={}
            cat_table[ii]["category"]=cat_counter
            cat_table[ii]["2c_length"]=cat_counter+1
            cat_table[ii]["DC_huffman_length"]=cat_length[cat_counter]
            ii+=1
        i=i*2
    return cat_table


if __name__=="__main__":
    category_table=gen_cat_table()
    image_name='wuflab_grayscale'
    filter_length=5
    image_start=Image.open('clean_files/'+image_name+'.png')
    image_start_array=list(image_start.getdata())
    width, height = image_start.size
    image_matrix=np.zeros(shape=(height,width))
    new_image_data=np.zeros(shape=(height,width),dtype=np.uint8)
    list_counter=0
    for i in range(0,height):
        for j in range(0,width):
            #print "i: {} j: {} list_counter: {} image_start_array[list_counter]: {}".format(i,j,list_counter,image_start_array[list_counter])
            if(type(image_start_array[list_counter])==tuple):
                image_matrix[i,j]=image_start_array[list_counter][0]
            else:
                image_matrix[i,j]=image_start_array[list_counter]
            list_counter+=1

    block_dict={}
    block_counter=0
    i=0
    j=0
    #set up data collection and figure out pixel blocks
    while i<height:
        while j<width:
            #get ranges for each block
            block_dict[block_counter]={}
            block_dict[block_counter]["width"]=[j, min(width+1,j+8)]
            block_dict[block_counter]["height"]=[i, min(height+1,i+8)]
            block_dict[block_counter]["bits_jpeg_traditional"]=0 #this field will hold data that tells us how many bits are required to encode this block based on traditional jpeg 
            block_dict[block_counter]["bits_jpeg_DNA"]=0 #this field holds information that tells us the bits required for
            block_dict[block_counter]["jpeg_elements"]=[]
            block_counter+=1
            j+=8
        i+=8
        j=0


    image_statistics_array=[]
    for block_iter in range(0,block_counter):
        image_statistics_array.append(process_pixel_block(block_dict[block_iter],image_matrix))
        filter_image(block_dict[block_iter],image_statistics_array[block_iter],filter_length,new_image_data)

    #collect information on the amount of data required for each coefficient in a block
    huffman_block_data(image_statistics_array,category_table)
    

    #analze statistics of blocks in the image
    block_statistics(image_statistics_array,image_name)



    new_image=Image.fromarray(new_image_data)
    new_image.save('filtered_'+image_name+str(filter_length)+'.png')
    print "ssim comparison measure: {}".format(measure.compare_ssim(image_matrix.astype(np.uint8,copy=False),new_image_data))
