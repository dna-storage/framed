import argparse
import os
from dnastorage.primer.primer_util import read_resumed_runs, hamming_distance
import random

def submit(jobs, timeout, notime):
  bash1 = "bsub -W 6:00 python ./primer_design.py --ranged " 
  bash2 = " --use-nupack --o output"
  bash3 = " --notime" 
  bash4 = " --timeout "
  r = 0 
  j = int((4**20)/jobs)+1

  for i in range(jobs): 
    command= bash1 + str(r) + ' ' + str(j) + bash2 + str(i) + ".txt" 
    if(notime):
        command += bash3 
    else: 
        command += bash4 + str(timeout)
    os.system(command)
    r+=j

def resume(timeout):
  bash1 = "bsub -W 6:00 python ./primer_design.py --ranged " 
  bash2 = " --use-nupack --o "
  bash3 = " --timeout "
  bash4 = " --primers "

  for files in os.listdir('./Primers'): 
    filename = os.path.join("./Primers",files)
    f = open(filename,"r")
    s,r = f.readline().split() 
    f.close()
    if(s == "finished"): 
      continue
    start = int(s)
    finish = int(r)
    finish += start
    command = bash1 + s + " " + str(finish) + bash2 + files + bash3 +str(timeout)+ bash4 + files
    os.system(command)


def collect(outputfile):
  #condense primers from all files into one list
  all_primers = []
  files_searched = 0
  primers = 0
  for files in os.listdir('./Primers'):
    filename = os.path.join("./Primers",files)
    portion=read_resumed_run(filename)
    for primer in portion:
        all_primers.append(primer)
        primers+=1
   # os.remove(filename)
    files_searched+=1

  #Write all primers into one file
  f= open(outputfile, "w")
  for primer in all_primers: 
     f.write("%s\n" % primer)
  f.close() 


  f = open("runInfo.txt", 'a')
  f.write("Collection:\n")
  f.write('total number of primers: '+str(primers)+'\n')
  f.write('number of files: '+str(files_searched)+'\n')
  f.close()

def pool(dist,masterlist, filename):
  
  i = 0
  dup = 0
  all_primers = read_primers(masterlist)
  random.shuffle(all_primers)

  #check hamming distances
  while(i<len(all_primers)):
    j=i
    while(j<len(all_primers)):
        if(hamming_distance(all_primers[i],all_primers[j])==0 and i!=j):
            del all_primers[j] 
            j-=1
            dup+=1
        if(hamming_distance(all_primers[i],all_primers[j])<=dist and i!=j):
            del all_primers[j] 
            j-=1
        j+=1       
    i+=1   

  #write output of one possible batch
  output=open(filename,'w')
  for primer in all_primers: 
    output.write("%s\n" % primer)
  output.close()

  #write run information
  f = open("runInfo.txt", 'a')
  f.write("Pool test:\n")
  f.write('passed hamming distance check: '+str(len(all_primers))+'\n')
  f.write("Found "+ str(dup)+ " duplicates \n")
  f.close()


def value(filename): 
  L=read_primers(filename)

  for primer in L: 
   value=0
   i=19
   for letter in primer: 
       if(letter=='T'):
          value+=3*4**i
       elif(letter=='C'):
          value+=2*4**i
       elif(letter=='G'):
          value+=4**i
       i-=1
   if(value>=4**20):
      print("Error")
   else: 
      print value

def delete(extension):
   direct = os.listdir("./") 
   count =0

   for files in direct: 
       if files.endswith(extension): 
           os.remove(files)
           count+=1

   f = open("runInfo.txt", 'a')
   f.write("removed " + str(count) + " files with "+ extension)
   f.close()	
 	
parser = argparse.ArgumentParser(description="Select primer generator utility to run")

parser.add_argument('--submit',dest="submit", help="Submit a number of jobs to hpc", type = int, action='store')

parser.add_argument('--timeout',type=int,dest="timeout",action="store",default=60, help="(Submit Option)If no primers are produced after timeout tries, give up.")

parser.add_argument('--notime',action="store_true", help="(Submit Option)Get rid of time limitations for run")

parser.add_argument('--resume',action="store_true", help="Resume a previous run with files stored in /Primers")

parser.add_argument('--collect',help="Compile all primer files into one master list, first argument is output filename second is bool to delete all read files", type=str,action='store')

parser.add_argument('--pool',dest="pool",nargs=3, help="Generate a pool from Master list of certain hamming distance give hamming distance first then masterlist filename then output-file name", action='store')

parser.add_argument('--value',dest="value", help="Output the numerical values of primers in a file",type=str, action='store')

#parser.add_argument('--delete',dest="delete", help="deletes files of certain file extention",type=str, action='store')

args = parser.parse_args()

if(args.submit != None):
  submit(args.submit, args.timeout, args.notime)

if(args.resume): 
  resume(args.timeout)

if(args.collect!=None):
  collect(args.collect)

if(args.pool != None):
  pool(args.pool[0],args.pool[1], args.pool[2])

#if(args.delete != None): disabled for now as use is limited
 # delete(args.delete)

if(args.value != None):
  value(args.value)
