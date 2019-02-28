#include "system.h"
#include "decoder.h"
#include "prep.h"
#include "generator.h"
#include "scheduler.h"
#include "storage_system.h"
#include "sequencer.h"
#include "utlist.h"
#include "parameters.h"
#include "stats.h"
#include <stdlib.h>


system_storage_t:: system_storage_t(storage_params_t storage_params){
  this->number_pools=storage_params.number_pools;
  this->pools=(pool_model_t*)malloc(sizeof(pool_model_t)*this->number_pools);
  this->pool_wait_time=storage_params.pool_wait_time;
  this->pool_copies=storage_params.pool_copies;
  this->number_reads=storage_params.number_reads;
  this->pool_write_time=storage_params.pool_write_time;
  this->average_pool_capacity=storage_params.average_pool_capacity;
  for(int i=0; i<this->number_pools;i++){
    this->pools[i].copies=(pool_char_t*)malloc(sizeof(pool_char_t)*this->pool_copies);
    //need to initialize the array of pool copies
    for(int j=0; j<this->pool_copies; j++){
      this->pools[i].copies[j].timer=0;
      this->pools[i].copies[j].remaining_reads=this->number_reads;
      this->pools[i].copies[j].in_use=0;
    }
  }
}

//free up space allocated for pools/files
system_storage_t::~system_storage_t(){
  for(int i=0; i<this->number_pools;i++) free(this->pools[i].copies);
  free(this->pools);
}





//looks to see if there is an availabel copy for the pool requested, returns the copy identifier of that pool if so, else return -1
int system_storage_t::storage_poolavailable(unsigned long pool_ID){
  pool_char_t* _pool_copies=this->pools[pool_ID].copies; //array of copies for the pool_ID

  for(int i=0; i<this->pool_copies; i++){
    if(_pool_copies[i].in_use!=1) return i;
  }
  return -1;

}

//decrement the read counter for the pool indicated by <pool_ID,pool_copy_ID>
//if the read counter hit 0 change in_use=1 and set the timer to the write time
//if the read counter != 0 change in_use=1 and set the timer for the non-write wait time
void system_storage_t::storage_readmanage(unsigned long pool_ID, unsigned long pool_copy_ID){
  unsigned long remaining_reads;
  this->pools[pool_ID].copies[pool_copy_ID].remaining_reads--;
  remaining_reads=this->pools[pool_ID].copies[pool_copy_ID].remaining_reads;
  if(remaining_reads==0){
    //used up all the reads, pool needs to be recovered in some way
    this->pools[pool_ID].copies[pool_copy_ID].in_use=1;
    this->pools[pool_ID].copies[pool_copy_ID].timer=this->pool_write_time-1;
  }
  else{
    this->pools[pool_ID].copies[pool_copy_ID].in_use=1;
    this->pools[pool_ID].copies[pool_copy_ID].timer=this->pool_wait_time-1;
  }
}


//restores pool copies that are marked as in use and have restoration timers at 0
// when pools are found to have this condition, if the number of reads left is equal to 0 that means that the pool required a re-write,
// thus refresh the number of reads to the base number of reads. If the number of reads is not equal to 0, just simply mark its in_use flag as 0 
void system_storage_t::storage_poolrestore(void){
  //iterate over all of the pools in the system
  for(int i=0; i<this->number_pools; i++){
    //iterate over all of the copies for each pool
    for(int j=0;j<this->pool_copies;j++){
      //restore timers if it is expired and it was in_use
      if(this->pools[i].copies[j].timer==0 && this->pools[i].copies[j].in_use==1){
	this->pools[i].copies[j].in_use=0;
	if(this->pools[i].copies[j].remaining_reads==0){
	  this->pools[i].copies[j].remaining_reads=this->number_reads; //restore the number of reads to the original value
	}
      }
    }
  }
}
  

