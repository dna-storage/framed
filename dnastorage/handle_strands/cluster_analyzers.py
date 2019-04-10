from dnastorage.fi import fault_injector
import Levenshtein as lv
import numpy as np
import matplotlib.pyplot as plt
import operator
#this cluster analysis technique tries to come up with a way to extract the correct strand using the context of all the other strands in the cluster. The idea is to try to understand what edits need to be done to transform a strand to the original strand.
def cluster_analyze_ed(strand_cluster,rate,expected_length,iterations):

    thresh=10#len(strand_cluster)*(1-rate**2*expected_length)
    _strand_cluster=strand_cluster
    for i in range(0,iterations):
        _strand_cluster_new=[]
        for start_strand in _strand_cluster:
            strand_ops={}
            strand_array=list(start_strand)
            for index, strand in enumerate(_strand_cluster):
                eds=lv.editops(start_strand,strand)
                for edit in eds:
                    edit_pos1=edit[1]
                    edit_pos2=edit[2]
                    edit_type=edit[0]
                    if edit_type=='delete':
                        edit_key=(edit_type,edit_pos1,'X')
                    else:
                        edit_char2=strand[edit_pos2]
                        edit_key=(edit_type,edit_pos1,edit_char2)
                    if edit_key not in strand_ops:
                        strand_ops[edit_key]=0
                        strand_ops[edit_key]+=1
                    else:
                        strand_ops[edit_key]+=1
            #print "Applyin operations\n"
            #go through the edit operations found and apply
            for op in sorted(strand_ops, key=lambda x:x[1],reverse=True):
                if strand_ops[op]>=thresh:
                    #print "opreation {} count {}\n".format(op,strand_ops[op])
                    #apply the edit
                    if op[0]=='replace':
                        #print "replace position: {} \n".format(op[1])
                        strand_array[op[1]]=op[2]
                    elif op[0]=='delete':
                        #print "delete"
                        left_side=strand_array[0:op[1]]
                        right_side=strand_array[op[1]+1:]
                        strand_array=left_side+right_side
                    elif op[0]=='insert':
                        #print "insert"
                        left_side=strand_array[0:op[1]]
                        right_side=strand_array[op[1]:]
                        strand_array=left_side+[op[2]]+right_side
            _strand_cluster_new.append("".join(strand_array))
            if i+1==iterations:
                return "".join(strand_array)
        _strand_cluster=list(_strand_cluster_new)


def cluster_analyze_majority(strand_cluster,strand_length):
    final_strand=[]
    pointer_list=[0]*len(strand_cluster)
    majority_vote={"A":0,"G":0,"C":0,"T":0}
    fix_error={"del":0,"inser":0,"sub":0}
    
    for final_strand_index in range(0,strand_length):
        majority_vote["A"]=0
        majority_vote["G"]=0
        majority_vote["C"]=0
        majority_vote["T"]=0
        for cluster_strand_index, cluster_strand in enumerate(strand_cluster):
            majority_vote[cluster_strand[pointer_list[cluster_strand_index]]]+=1 #sum up across all the current pointers

        #apply the majority vote
        majority_base=max(majority_vote.iteritems(), key=operator.itemgetter(1))[0];
        final_strand.append(majority_base)
        for cluster_strand_index, cluster_strand in enumerate(strand_cluster):
            if cluster_strand[pointer_list[cluster_strand_index]]!=majority_base:
                fix_error["del"]=0
                fix_error["inser"]=0
                fix_error["sub"]=0
                for cluster_strand_index_fix, cluster_strand_fix in enumerate(strand_cluster):
                    sub_index=min(pointer_list[cluster_strand_index]+1,len(cluster_strand)-1)
                    sub_index_fix=min(pointer_list[cluster_strand_index_fix]+1,len(cluster_strand_fix)-1)
                    del_index=pointer_list[cluster_strand_index]
                    del_index_fix=min(pointer_list[cluster_strand_index_fix]+1,len(cluster_strand_fix)-1)
                    insert_index=min(pointer_list[cluster_strand_index]+1,len(cluster_strand)-1)
                    insert_index_fix=pointer_list[cluster_strand_index_fix]
                    if cluster_strand_fix[pointer_list[cluster_strand_index_fix]]==majority_base:
                        if cluster_strand[sub_index]==cluster_strand_fix[sub_index_fix]:
                            fix_error["sub"]+=1
                        elif cluster_strand[del_index]==cluster_strand_fix[del_index_fix]:
                            fix_error["del"]+=1
                        elif cluster_strand[insert_index]==cluster_strand_fix[insert_index_fix]:
                            fix_error["inser"]+=1
                majority_fix=max(fix_error.iteritems(), key=operator.itemgetter(1))[0]
                if majority_fix=="inser":
                    pointer_list[cluster_strand_index]=min(len(cluster_strand)-1,pointer_list[cluster_strand_index]+2)
                elif majority_fix=="sub":
                    pointer_list[cluster_strand_index]=min(len(cluster_strand)-1,pointer_list[cluster_strand_index]+1)
            else:
                #just increment if it is in the majority
                pointer_list[cluster_strand_index]=min(len(cluster_strand)-1,pointer_list[cluster_strand_index]+1)
    return "".join(final_strand)


        
if __name__=="__main__":
    base_strand="ATCCTACTCCAGCGGGATCTTTTATCTAAAGACGATGAGAGGAGTATTCGTCAGGCCACATAGCTTTCATGTTCTGATCGGAACGATCGTTGGCGCCCGACCCTCGGATTCCGTAGTGAGTTCTTTGTCCGAGCCATTGTATGCGAGATCGGTAGACTGATAGGGGATGCAGTATATCCCTGGATGCAATAGACGGACAG"

    import argparse
    parser = argparse.ArgumentParser(description="Analyze different cluster analysis methods")
    parser.add_argument('--analysis_type', dest="cluster_type",action="store", default=None, help='type of analysis to perform')


    args = parser.parse_args()
    
    injector_args=fault_injector.fault_injector_arguments()
    injector_args.o="outfile.log"
    injector_args.input_file=None
    injector_args.fault_file=None
    injector_args.fault_model="fixed_rate"
    injector_args.run=False
    injector_args.faulty=0
    injector_args.fails=0
    injector_args.p1=0
    injector_args.p2=0
    fault_model=eval('fault_injector.'+"fixed_rate"+"(injector_args)")
    data_set={}
    
    for strand_count_step in range(60,80,5):
        if strand_count_step not in data_set:
            data_set[str(strand_count_step)+"_strands"]=[]
        base_set=[base_strand for i in range(0,strand_count_step)]
        for rate_step in np.arange(0.01,0.07,0.01):
            correct_count=0
            hist_fig, hist_axes=plt.subplots(nrows=1, ncols=1)
            distance_after_array=[]
            distance_before_array=[]
            for trial_number in range(0,100):
                fault_model.set_fault_rate(rate_step)
                fault_model.set_library(base_set)
                error_set=fault_model.Run()
                distance_before=lv.distance(error_set[0],base_strand)
                distance_before_array.append(distance_before)
                if args.cluster_type=="ed":
                    return_strand=cluster_analyze_ed(error_set,rate_step,200,1)
                elif args.cluster_type=="BMA":
                    return_strand=cluster_analyze_majority(error_set,200)
                distance_after=lv.distance(return_strand,base_strand)
                if distance_after==0:
                    correct_count+=1
                distance_after_array.append(distance_after)
                print "trial: {} complete\n".format(trial_number)
            #plot histograms for edit distance before and after
            hist_axes.clear()
            hist_axes.hist(distance_before_array,bins=10)
            hist_fig.savefig("hist_before_rate_"+str(rate_step)+".png")
            hist_axes.clear()
            hist_axes.hist(distance_after_array,bins=10)
            hist_fig.savefig("hist_after_rate_"+str(rate_step)+".png")

            
            percent_correct=100.0*(float(correct_count)/100.0)
            data_set[str(strand_count_step)+"_strands"].append(percent_correct)
            
    rate_array=np.arange(0.01,0.07,0.01)
    rate_fig, rate_ax=plt.subplots(nrows=1,ncols=1)
    for strand_step in data_set:
        rate_ax.plot(rate_array,data_set[strand_step], label=strand_step)
    rate_ax.legend()
    rate_fig.savefig("edit_distance_correction_50_70.png")



