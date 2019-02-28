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
#include <stdio.h>
#include <math.h>


scheduler_t::scheduler_t(scheduler_params_t scheduler_params){

  this->_prep=scheduler_params._prep;
  this->_system=scheduler_params._system;
  this->_storage=scheduler_params._storage;
  this->system_queue=scheduler_params.system_queue;
  this->bytes_per_strand=scheduler_params.bytes_per_strand;
  this->sequencing_depth=scheduler_params.sequencing_depth;
  this->efficiency=scheduler_params.efficiency;
  this->reorder=&scheduler_t::reorder_none;
  this->scheduler=&scheduler_t::scheduler_anypool;
  this->strand_calculator=&scheduler_t::calc_singlefile;
  this->batch_size=scheduler_params.batch_size; 
}

//top level function called by the system sim
void scheduler_t::schedule_stage(void){
  (this->*reorder)(); //this function reorders I/O in the system_queue
  (this->*scheduler)(); //this function injects transactions into prep stations 
}



//empty function for no re-ordering default setting
void scheduler_t::reorder_none(void){
}

//this scheduling policy indicates that any pool can be put with one another 
void scheduler_t::scheduler_anypool(void){

  int prep_ID;
  unsigned long transaction_ID;
  system_sim_t* _system=this->_system;
  transaction_t* _window=_system->window;
  system_storage_t* _storage=this->_storage;
  trace_t* head;
  unsigned long desired_strands_sequenced;
  unsigned long undesired_strands_sequenced;
  unsigned long strands_to_sequence;
  prep_t* _prep=this->_prep;
  while((prep_ID=_prep->prep_stationavail()>0) && *(this->system_queue)!=NULL && !_system->window_full()){
    head=*(this->system_queue);
    //keep going while there are things on the system queue and there are prep stations available
    if(!_storage->storage_poolavailable(head->pool_ID)) break;

    //get a spot on the window
    transaction_ID=_system->window_add();
    
    //set up the transaction_t structure at window_head
    _system->window_init(transaction_ID);

    //loop until we can no longer batch anymore transactions
    //each transaction needs to calculate out its strands, and then be accumulated with the top level transaction indicated by transaction_ID and then added to the component linked list
    for(int i=0; i<this->batch_size; i++){
      if(*(this->system_queue)==NULL) break; //break if nothing left on the queue
      (this->*strand_calculator)(desired_strands_sequenced,undesired_strands_sequenced,
			      head->file_size, head->pool_ID, this->efficiency,
			      this->bytes_per_strand, this->sequencing_depth);
      _system->window_componentadd(transaction_ID, undesired_strands_sequenced,
				   desired_strands_sequenced, head->pool_ID,
				   head->file_size);
      //need to take the request off the top of the queue and remove the space for it
      _system->queue_pop();
    }
    //submit the transaction to the previously found prep station
    _prep->prep_stationsubmit(prep_ID,transaction_ID);
  }
}


/*
This strand calculator uses the file_size, efficiency, bytes_per_strand,
and sequencing_depth to calculate the amount of desired_strands_sequenced
and the undesired_strands_sequenced. This calculator should be used when 
modeling the effect of accessing single individual files from a pool. Another
calculator should be used to model sequencing the full pool.
*/
void scheduler_t::calc_singlefile(unsigned long& desired_strands_sequenced,
				  unsigned long& undesired_strands_sequenced,
				  unsigned long file_size,
				  unsigned long pool_ID,
				  float efficiency,
				  float bytes_per_strand,
				  unsigned long sequencing_depth)
{
  
  float desired_strands;
  float undesired_strands;
  float _file_size;
  float _sequencing_depth;

  _sequencing_depth=(float)sequencing_depth;
  _file_size=(float)file_size*(float)FILE_UNIT;
  desired_strands=ceil(_file_size/bytes_per_strand)*_sequencing_depth;
  undesired_strands=ceil(desired_strands/efficiency)-desired_strands;
  desired_strands_sequenced=(unsigned long)desired_strands;
  undesired_strands_sequenced=(unsigned long)undesired_strands;
}
