#include "dna_storage_attributes.h"
#include "decoder.h"
#include "prep.h"
#include "sequencer.h"
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
