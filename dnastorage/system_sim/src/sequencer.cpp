#include "system.h"
#include "decoder.h"
#include "prep.h"
#include "generator.h"
#include "scheduler.h"
#include "storage_system.h"
#include "sequencer.h"
#include "utlist.h"
#include "buffer.h"
#include "parameters.h"
#include "stats.h"
#include <stdlib.h>
sequencer_unit_t::sequencer_unit_t(unsigned long max_sequencing,unsigned long num_channels):system_unit_t(num_channels){
  this->max_strands=max_sequencing;
  this->used_strands=0;
  this->wasted_strands=0;
  this->timeout=0;
}


sequencer_t:: sequencer_t(sequencer_params_t sequencer_params){
  this->num_sequencers=sequencer_params.num_sequencers;
  this->_system=sequencer_params._system;
  this->seq_dec_buffer=sequencer_params.seq_dec_buffer;
  this->prep_seq_buffer=sequencer_params.prep_seq_buffer;
  this->base_timer=sequencer_params.timer;
  this->base_timeout=sequencer_params.base_timeout;
  this->stats=sequencer_params.stats;
  //allocate sequencer units
  this->sequencer_set=(sequencer_unit_t**)malloc(sizeof(sequencer_unit_t*)*this->num_sequencers);
  for(int i=0; i<this->num_sequencers; i++) this->sequencer_set[i]=new sequencer_unit_t(sequencer_params.max_strands,1);
}

sequencer_t::~sequencer_t(){
  //free up the sequencer list
  for(int i=0; i<this->num_sequencers; i++) delete this->sequencer_set[i];
  free(this->sequencer_set);
}


//wrapper for the sequencer stage
void sequencer_t::sequencer_stage(void){
  this->sequencer_backend();
  this->sequencer_frontend();
}


 
//look through the sequencers and find finished ones
//implements the backend completion of sequenced transactions
void sequencer_t::sequencer_backend(void){
  sequencer_unit_t* _sequencer;
  //iterate through the sequencers
  for(unsigned long i=0;i<(this->num_sequencers);i++){
    _sequencer=this->sequencer_set[i];
    unsigned long _timer=_sequencer->timer;
    int active=_sequencer->unit_active;
    if(active && _timer==0){
      //found a done sequencer, need to move each transaction in the list to the seq_dec_buffer
      this->sequencer_complete(i); //this function handles all of the steps needed for opening a decoder
    }
    else if(active && _timer>0){
      //active sequencer, decrement the timer
      this->sequencer_timestep(i);
    }
  }
}

//this function implements the frontend of the sequencer
void sequencer_t::sequencer_frontend(void){
  sequencer_unit_t* _sequencer;
  unsigned long transaction_strands; //number of strands that are still needed to be sequenced
  transaction_t* _window=this->_system->window;
  unsigned long transaction_ID;
  buffer_t* _p_s=this->prep_seq_buffer;
  list_entry_t* prepped;
  //iterate over all of the transactions in the prep_seq_buffer
  for(_p_s->iter_start(); _p_s->iter_get()!=NULL; _p_s->iter_next()){
    prepped=_p_s->iter_get();
    if(prepped->used){
      transaction_ID=prepped->transaction_index;
      int sequencer_ID=this->sequencer_avail(transaction_ID);
      if(sequencer_ID>=0){
	this->sequencer_submit(transaction_ID,sequencer_ID,this->prep_seq_buffer->get_iterator()); //submit the transaction to the sequencer
      }
    }
  }
  //kickoff sequencers
  this->sequencer_kickoff();
  this->sequencer_timeoutstep();
}

//find sequencers that are setup to be kicked off
void sequencer_t::sequencer_kickoff(void){
  sequencer_unit_t* _sequencer;
  //iterate through the sequencers and find sequencers ready for kickoff
  for(int i=0; i<this->num_sequencers; i++){
    _sequencer=this->sequencer_set[i];
    if(!_sequencer->unit_active && ((_sequencer->timeout==0 && _sequencer->next_open!=0) || _sequencer->used_strands==_sequencer->max_strands || _sequencer->next_open==_sequencer->num_channels)){
      //printf("kickoff sequencer %i\n",i);
      _sequencer->unit_active=1; //activate the unit
      _sequencer->timer=this->base_timer-1; //initialize the timer for the sequencer
    }
  }
}

void sequencer_t::sequencer_timeoutstep(void){
  sequencer_unit_t* _sequencer;
  //decrement timeout timers for sequencers
  for(int i=0; i<this->num_sequencers;i++){
    _sequencer=this->sequencer_set[i];
    if(_sequencer->timeout>0 && !_sequencer->unit_active && _sequencer->next_open!=0) _sequencer->timeout--;
    //   if(_sequencer->timeout>0 && !_sequencer->unit_active && _sequencer->next_open!=0) //printf("sequencer_ID %i, sequencer_timeout %i\n",i,_sequencer->timeout);
  }
  
}


//function to place a transaction into a sequencer
void sequencer_t::sequencer_submit(unsigned long transaction_ID, unsigned long sequencer_ID, unsigned long prep_seq_buffer_index){
  sequencer_unit_t* _sequencer;
  unsigned long transaction_strands;
  transaction_t* _window=this->_system->window;
  
  //found a sequencer that works with the transaction, need to place the transaction in that sequencer
  _sequencer=this->sequencer_set[sequencer_ID];
  transaction_strands=_window[transaction_ID].strands_to_sequence;
  if(transaction_strands>(_sequencer->max_strands-_sequencer->used_strands)){
    //have more transaction strands than whats available in the sequencer
    _sequencer->used_strands=_sequencer->max_strands;
    _window[transaction_ID].strands_to_sequence-=(_sequencer->max_strands-_sequencer->used_strands);
    
    //put the undesired strands in the sequencer 
    if(_window[transaction_ID].undesired_strands_sequenced>0){
      if(_window[transaction_ID].undesired_strands_sequenced>_sequencer->max_strands-_sequencer->used_strands){
	_window[transaction_ID].undesired_strands_sequenced-=(_sequencer->max_strands-_sequencer->used_strands);
	_sequencer->wasted_strands+=(_sequencer->max_strands-_sequencer->used_strands);
      }
      else{
	_sequencer->wasted_strands+=_window[transaction_ID].undesired_strands_sequenced;
	_window[transaction_ID].undesired_strands_sequenced=0;
      }
    }
    _window[transaction_ID].cracked_count++; //split transaction up and take note of it
  }
  else{
    //transaction will fit in the remaining sequencer
    _sequencer->used_strands+=_window[transaction_ID].strands_to_sequence;
    _sequencer->wasted_strands+=_window[transaction_ID].undesired_strands_sequenced;
    //free up the _p_s entry
    this->prep_seq_buffer->free_entry(prep_seq_buffer_index);
  }
  //make some changes to the and add the transactions to the sequencer
  _sequencer->transaction_slots[_sequencer->next_open]=transaction_ID;
  if(_sequencer->next_open==0) _sequencer->timeout=this->base_timeout; //set the timeout counter for the sequencer that just got its first transaction 
  _sequencer->next_open++;
}



//this function tries to find a sequencer that will work with transaction requested
//if there is no possible sequencer, then return -1
int sequencer_t::sequencer_avail(unsigned long transaction_ID){
  sequencer_unit_t* _sequencer;  
  //iterate over all of the sequencers
  for(unsigned long i=0;i<(this->num_sequencers);i++){
    _sequencer=this->sequencer_set[i];
    //found an inactive sequencer
    if(!_sequencer->unit_active && _sequencer->used_strands < _sequencer->max_strands
       && _sequencer->next_open<_sequencer->num_channels){
      if(_sequencer->next_open==0){
	//know that this sequencer is unused and unscheduled by any other transaction
	return i;
      }
    }
  }
}





//complete the decoder by placing the transactions into the seq_dec_buffer
void sequencer_t::sequencer_complete(unsigned long sequencer_ID){
  unsigned long number_transactions=(this->sequencer_set[sequencer_ID])->next_open;
  sequencer_unit_t* _sequencer=this->sequencer_set[sequencer_ID];
  transaction_t* _window=this->_system->window;
  int test_result;
  unsigned long buffer_index;
  for(_sequencer->transaction_pointer;
      _sequencer->transaction_pointer < number_transactions;
      _sequencer->transaction_pointer++){
    
    unsigned long i=_sequencer->transaction_pointer;
    unsigned long transaction_number=_sequencer->transaction_slots[i];
    test_result= this->seq_dec_buffer->test_transaction(transaction_number,buffer_index);
    if(test_result==0 && buffer_index==-1) return; //no more buffer spots 
    
    if(test_result==0 && buffer_index!=-1){ //need a new buffer spot
      //put the transactions that can fit into the seq_dec_buffer, and mark it used
      this->seq_dec_buffer->init_entry(transaction_number,buffer_index);
      _window[transaction_number].cracked_count--;
    }
    else if(test_result==1){
      //transaction already in the buffer, just decrement the cracked count for the transaction
      _window[transaction_number].cracked_count--;
    }
    
  }
  //if we made it here, we moved all of the transactions from the sequencer to the buffer
  _sequencer->unit_active=0; //kill the activity of the unit
  _sequencer->used_strands=0; //reset the used strands
  _sequencer->next_open=0; //reset the next open pointer
  _sequencer->transaction_pointer=0;//reset the transaction pointer
}



//decreases the timer for a sequencer that is active
void sequencer_t::sequencer_timestep(unsigned long sequencer_ID){
  sequencer_unit_t* _sequencer=this->sequencer_set[sequencer_ID];
  _sequencer->timer--;
  //printf("sequencer_ID %i, sequencer timer %i\n",sequencer_ID,_sequencer->timer);
}
