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
prep_unit_t::prep_unit_t(unsigned long num_channels) : system_unit_t(num_channels){}

//constructor for the prep component of the system
prep_t::prep_t(prep_params_t prep_params){
  
  this->num_preps=prep_params.num_preps;
  this->_system=prep_params._system;
  this->prep_seq_buffer=prep_params.prep_seq_buffer;
  this->dna_storage=prep_params.dna_storage;
  this->base_timer=prep_params.timer;
  this->stats=prep_params.stats;
  //make the list of prep units
  this->prep_set=(prep_unit_t**)malloc(sizeof(prep_unit_t*)*this->num_preps);
  for(int i=0; i<this->num_preps; i++) this->prep_set[i]=new prep_unit_t(prep_params.num_channels);
}


//top level function wrapper for the prep stage 
void prep_t::prep_stage(void){
  this->prep_backend();
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


//looks to see if there is a prep station that will accept the pool indicated by pool_ID
int prep_t::prep_stationavail(void){
  prep_unit_t* _prep;
  //iterate through the prep units
  for(int i=0; i<this->num_preps;i++){
    _prep=this->prep_set[i];
    if(!_prep->unit_active && _prep->next_open<_prep->num_channels){
      
      if(_prep->next_open==0) return i; //open pool, return ID
    }
  }
  return -1; //inidcate no prep station available 
}

//submit a transaction to the specifiec station
void prep_t::prep_stationsubmit(unsigned long prep_ID, unsigned long transaction_ID){
  prep_unit_t* _prep=this->prep_set[prep_ID];
  _prep->timer=this->base_timer-1;
  _prep->unit_active=1;
  _prep->transaction_slots[_prep->next_open]=transaction_ID;
  _prep->next_open++;
}

//move transactions from the prep unit to the prep_seq_buffer
void prep_t::prep_complete(unsigned long prep_ID){
  unsigned long number_transactions=this->prep_set[prep_ID]->next_open;
  prep_unit_t* _prep=this->prep_set[prep_ID];
  transaction_t* _window=this->_system->window;
  for(_prep->transaction_pointer;
      _prep->transaction_pointer < number_transactions;
      _prep->transaction_pointer++){

    unsigned long i = _prep->transaction_pointer;
    unsigned long transaction_number=_prep->transaction_slots[i];
    int prep_seq_index=this->prep_seq_buffer->get_free();
    if(prep_seq_index==-1) return; //no more buffer spots
    if(prep_seq_index>=0){
      this->prep_seq_buffer->init_entry(transaction_number,prep_seq_index);
    }
  }
  //relinquish the prep unit
  _prep->unit_active=0;
  _prep->next_open=0;
  _prep->transaction_pointer=0;
}



//timestep the given prep unit
void prep_t::prep_timestep(unsigned long prep_ID){
  prep_unit_t* _prep=this->prep_set[prep_ID];
  _prep->timer--;

}
