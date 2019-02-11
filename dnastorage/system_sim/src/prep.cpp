#include "dna_storage_attributes.h"
#include "decoder.h"
#include "prep.h"
#include "sequencer.h"
#include "utlist.h"
#include <stdlib.h>
prep_unit_t::prep_unit_t(unsigned long timer, unsigned long num_channels) : system_unit_t(num_channels){}

//constructor for the prep component of the system
prep_t::prep_t(unsigned long timer, unsigned long num_channels,
	       unsigned long buffer_size,unsigned long num_preps,list_entry_t* prep_seq_buffer,
	       system_sim_t* _system,system_storage_t* dna_storage){
  
  this->num_preps=num_preps;
  this->_system=_system;
  this->prep_seq_buffer=prep_seq_buffer;
  this->buffer_size=buffer_size;
  this->dna_storage=dna_storage;
  this->base_timer=timer;
  //make the list of prep units
  this->prep_set=(prep_unit_t**)malloc(sizeof(prep_unit_t*)*this->num_preps);
  for(int i=0; i<this->num_preps; i++) this->prep_set[i]=new prep_unit_t(num_channels);

  //assign the prep policy here
  prep_poolpolicy=&prep_t::single_pool;
}


prep_t::~prep_t(){
  //free up the prep list
  for(int i=0; i<this->num_preps; i++) delete this->prep_set[i];
  free(this->prep_set);
}


//look through and find done prep stations, move them to the buffer station
//also decrement the timer of done prep stations
void prep_t::prep_backend(void){
  prep_unit_t* _prep;
  //iterate through the prep stations
  for(unsigned long i=0; i<this->num_preps; i++){
    _prep=this->prep_set[i];
    unsigned long _timer=_prep->timer;
    int active=_prep->unit_active;
    if(active && _timer==0) this->prep_complete(i);
    else if(active && _timer>0) this->prep_timestep(i);
  }
}

//front end for the prep stations
/*
  Steps:
  1. Need to put transactions that are ready from the trace into an available prep station
     If a transaction from the trace buffer is used, need to replace with the next trace in the trace file
  2. When bringin a new transaction into the pipeline need to initialize the transaction_t structure that will represent it in the window[] array
  3. After using a transaction, need to make changes to the dna storage unit.
     Need to decrement the reads left on the pool used, if it reaches 0 need to start the write timer for that pool. If the number of reads is not exhausted, need to start the next-use-timer to model the time it takes to get the pool back for use.
  4. After we can no longer bring transactions into the pipeline, need to fix up the pools by decrementing timers and .

 */
void prep_t::prep_frontend(void){
  system_sim_t* _system=this->_system;
  system_storage_t* _dna_storage=this->dna_storage;
  prep_unit_t* _prep;
  trace_t* trace_trans;
  int trace_count=0;
  int pool_copy;
  int prep_ID;
  int pool_in_standby;
  unsigned long transaction_ID;
  LL_FOREACH(_system->trace_list_head,trace_trans){
    trace_count++;//count the number of traces, if we see later that it is 0 then the trace has been fully gone through
    pool_in_standby=this->prep_poolstandby(trace_trans->pool_ID);
    if(pool_in_standby<0){
      //look for a prep station that is available
      pool_copy=_dna_storage->prep_poolavailable(trace_trans->pool_ID);
      //check to see if the trace's pool is available
      if(pool_copy<0) continue;
      //check to see if there is a prep station that can be used
      prep_ID=this->prep_stationavail(trace_trans->pool_ID);
      if(prep_ID<0) continue;
      //made it to this point, use the prep station and pool copy found previously
      //at this point need to make the appropriate changes to the prep station to be used and the copy of the pool
      _prep=this->prep_set[prep_ID];
      if(_prep->next_open==0) _prep->standby_timer=this->base_standby_timer;
      _dna_storage->storage_readmanage(trace_trans->pool_ID,pool_copy);
    }
    else{
      //found a prep station that already has the pool on standby, so just use that without affecting the dna storage
      _prep=this->prep_set[pool_in_standby]; 
    }
    //need to allocate a transaction structure in the window
    transaction_ID=this->prep_windowallocate(trace_trans->pool_ID);
    //add the instruction to the found prep station
    _prep->transaction_slots[next_open]=transaction_ID;
    _prep->next_open++; 
  }
  if(trace_count==0)_system->trace_complete=1; //nothing left on the list   
  //need to got through the dna storage unit and fix pool counters and kickoff prep stations
  this->prep_kickoff();
  this->prep_standbystep();
  _dna_storage->storage_poolrestore(); //call to restore pools
}

//iterate through prep stations and find standby_timers==0. If so, kick them off by setting the run-time timer and the active flag
void prep_t::prep_kickoff(void){
  prep_unit_t* _prep;
  //iterate through the prep stations 
  for(int i=0; i<this->num_preps;i++){
    _prep=this->prep_set[i];
    //look for stations that have expired staby timers or have already used up all of their channels
    if((_prep->standby_timer==0 && _prep->next_open!=0) || _prep->next_open==_prep->num_channels){
      _prep->unit_active=1;
      _prep->timer=this->base_timer-1;
    }
  }
}
void prep_t::pre_standbystep(void){
  prep_unit_t* _prep;
  //iterate through the prep stations and decrement the standby_timer
  for(int i=0; i<this->num_preps;i++){
    _prep=this->prep_set[i];
    if(_prep->standby_timer>0 && !_prep->unit_active && _prep->next_open!=0) _prep->standby_timer--;
  }
}

//initialize a new transaction at the window tail using pool_ID
unsigned long prep_t::windowallocate(unsigned long pool_ID){
  transaction_t* _window=this->_system->window;
  unsigned long transaction_ID;

  transaction_ID=this->_system->window_tail;
  this->_system->window_tail++;
  if(this->_system->window_tail==WINDOW_SIZE) this->_system->window_tail=0;

  //initialize fields
  _window[transaction_ID].cracked_count=1; //cracked count = 1 signifies a transaction that has not been split yet
  _window[transaction_ID].pool_ID=pool_ID;
  
  //use the calculate policy in order to determine the strand count fields of the transaction
  calc_policy(transaction_ID); //this is a function pointer that must be assigned in the constructur for the prep_t class  
}

//look to see if there is a pool that we can use already in standby 
int prep_t::prep_poolstandby(unsigned long pool_ID){
  prep_unit_t* _prep;
  transaction_t* _window=this->_system->window;
  for(int i=0; i<this->num_preps; i++){
    _prep=this->prep_set[i];
    if(!_prep->unit_active && _prep->next_open!=_prep->num_channels){
      for(int j=0; j<_prep->next_open;j++){
	if(_window[_prep->transaction_slots[j]].pool_ID==pool_ID) return i; //return the prep station ID that has the pool on standby
      }
    }
  }
  return -1; //did not find a prep station with the pool on standby
}


//looks to see if there is a prep station that will accept the pool indicated by pool_ID
int prep_t::prep_stationavail(unsigned long pool_ID){
  prep_unit_t _prep;

  //iterate through the prep units
  for(int i=0; i<this->num_preps;i++){
    _prep=this->prep_set[i];
    if(!_prep->unit_active && _prep->next_open<_prep->num_channels){
      
      if(_prep->next_open==0) return i; //automatically know we can use this prep station
      else{
	//need to make sure the inquired pool fits the policy
	if(this->prep_poolpolicy(pool_ID,i)) return i;
      }
    }
  }
}

//move transactions from the prep unit to the prep_seq_buffer
void prep_t::prep_complete(unsigned long prep_ID){
  unsigned long number_transactions=this->prep_set[prep_ID]->next_open;
  prep_unit_t* _prep=this->prep_set[prep_ID];
  list_entry_t* _p_s=this->prep_seq_buffer;
  transaction_t* _window=this->_system->window;
  for(_prep->transaction_pointer;
      _prep->transaction_pointer < number_transactions;
      _prep->transaction_pointer++){

    unsigned long i = _prep->transaction_pointer;
    unsigned long transaction_number=_prep->transaction_slots[i];
    int prep_seq_index=this->get_prepseq(transaction_number);
    if(prep_seq_index==-1) return; //no more buffer spots
    if(prep_seq_index>=0){
      _p_s[prep_seq_index].used=1;
      _p_s[prep_seq_index].transaction_index=transaction_number;
    }
  }
  //relinquish the prep unit
  _prep->unit_active=0;
  _prep->next_open=0;
  _prep->transaction_pointer=0;

}

//find an open spot in the prep_seq_buffer
void prep_t::get_prepseq(unsigned long transaction_ID){
  list_entry_t* _p_s=this->prep_seq_buffer;

  //find an open prep_seq_buffer location
  for(int i=0; i<QUEUE_SIZE; i++){
    if(_p_s[i].used==0) return i; //return an unused list entry
  }
  return -1; //did not find an unused entry

}


//timestep the given prep unit
void prep_t::prep_timestep(unsigned long prep_ID){
  prep_unit_t* _prep=this->prep_set[prep_ID];
  _prep->timer--;

}

//return 1 if the pool does not exist in the prep station yet
//return 0 if the pool does exist in the prep station
int prep_t::single_pool(unsigned long pool_ID, unsigned long prep_ID){
  prep_unit_t* _prep=this->prep_set[prep_ID];
  unsigned long prep_trans;
  transaction_t* _window=this->_system-window;
  //iterate through all the transactions on a given prep station 
  for(int i=0; i<_prep->next_open;i++){
    prep_trans=_prep->transaction_slots[i];
    if(_window[prep_trans].pool_ID==pool_ID) return 0;
    
  }
  return 1;
}
