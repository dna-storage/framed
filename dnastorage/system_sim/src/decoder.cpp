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
#include <stdio.h>
#include <stdlib.h>
decoder_unit_t::decoder_unit_t(unsigned long num_channels) : system_unit_t(num_channels){}

decoder_t::decoder_t(decoder_params_t decoder_params)
{
  
  this->num_decoders=decoder_params.num_decoders;
  this->_system=decoder_params._system;
  this->seq_dec_buffer=decoder_params.seq_dec_buffer;
  this->base_timer=decoder_params.timer;
  this->stats=decoder_params.stats;
  //allocate the decoder set
  this->decoder_set=(decoder_unit_t**)malloc(sizeof(decoder_unit_t*)*this->num_decoders);
  for(int i=0; i<this->num_decoders; i++) this->decoder_set[i]=new decoder_unit_t(1);

}

//destructor for the decoder system component
decoder_t::~decoder_t(){
  //deallocate items
  for(int i=0; i<this->num_decoders; i++) delete this->decoder_set[i];
  free(this->decoder_set);
}


//wrapper function for the decoder stage
void decoder_t::decoder_stage(void){
  this->decoder_backend();
  this->decoder_frontend();
}



//function that looks through the decoders and see which transactions are finished
//increments the transactions complete counter at the end of checking the decoder units
void decoder_t::decoder_backend(void){
  unsigned long done_count=0;
  system_sim_t* _system_sim=this->_system;
  decoder_unit_t* _decoder;
  //iterate through the decoders and look for complete jobs
  for(unsigned long i=0; i<this->num_decoders;i++){
    _decoder=this->decoder_set[i];
    int active=_decoder->unit_active;
    unsigned long timer=_decoder->timer;
    if(timer==0 && active){
      this->decoder_complete(i); //complete the transaction within the decoder
      done_count++;
    }
    else if(active && timer>0) this->decoder_timestep(i); //decrease the time on the decoder
  }
}

//need to fill up non-active decoders with transactions in the seq_dec_buffer
void decoder_t::decoder_frontend(void){
  buffer_t* _s_d=this->seq_dec_buffer;
  transaction_t* _window=this->_system->window;
  list_entry_t* seq_buffer;
  
  for(_s_d->iter_start();_s_d->iter_get()!=NULL;_s_d->iter_next()){
    seq_buffer=_s_d->iter_get();
    if(seq_buffer->used && _window[seq_buffer->transaction_index].cracked_count==0){
      //found a location in the buffer that is ready
      unsigned long decoder_ID=this->decoder_avail();
      unsigned long transaction_ID=seq_buffer->transaction_index;
      if(decoder_ID==-1) break; //no decoder avilable, stop looking through the seq_dec_buffer
      this->init_decoder(decoder_ID,transaction_ID,_s_d->get_iterator()); //initialize the decoder found thats available
    }
  }
}



//reach into the decoder set and deactivate the decoder and complete the transaction
void decoder_t::decoder_complete(unsigned long decoder_ID){
  decoder_unit_t* _decoder=this->decoder_set[decoder_ID];
  transaction_t* temp;
  transaction_t* comp_head;
  transaction_t* _window=this->_system->window;
  int batch_complete=1;
  _decoder->unit_active=0;
  comp_head=_window[_decoder->transaction_slots[0]].components;
  _window[_decoder->transaction_slots[0]].components[_decoder->component_ID].transaction_finished=1;
  //find out when the component finished
  _window[_decoder->transaction_slots[0]].components[_decoder->component_ID].time_stamp_end=counter(time_step);
  //printf("decoder_ID %i complete\n",decoder_ID);
  //check to see if all components are done
  LL_FOREACH(comp_head,temp){
    if(temp->transaction_finished!=1) batch_complete=0;
  }
  if(batch_complete) _window[_decoder->transaction_slots[0]].transaction_finished=1; //all pieces done in the original batch
}

void decoder_t::decoder_timestep(unsigned long decoder_ID){
  (this->decoder_set[decoder_ID])->timer--;
  //printf("decoder %i timer %i\n",decoder_ID,(this->decoder_set[decoder_ID])->timer);
}

//search for deactivated decoders
unsigned long decoder_t::decoder_avail(void){
  decoder_unit_t* _decoder;
  for(unsigned long i=0; i<(this->num_decoders);i++){
    _decoder=this->decoder_set[i];
    if (_decoder->unit_active==0){
      return i;
    }
  }
  return -1;
}

//initialize the chosen decoder unit with timer, active signals etc
void decoder_t::init_decoder(unsigned long decoder_ID, unsigned long transaction_ID, unsigned long seq_dec_index){
  transaction_t* _window=this->_system->window;
  transaction_t* comp_head=_window[transaction_ID].components;
  transaction_t* comp_temp;
  decoder_unit_t* _decoder=this->decoder_set[decoder_ID];
  int decoding_left=0;
  int comp_ID=0;
  _decoder->next_open=0;
  _decoder->transaction_slots[_decoder->next_open]=transaction_ID;
  _decoder->unit_active=1;
  _decoder->timer=(this->base_timer);
  //need to make sure the sequencer gets the component ID
  LL_FOREACH(comp_head,comp_temp){
    if(comp_temp->component_decoded==0){
      _decoder->component_ID=comp_ID;
      comp_temp->component_decoded=1;
      break;
    }
    comp_ID++;
  }

  //take care of the buffer entry if it is no longer needed
  LL_FOREACH(comp_head,comp_temp){
    if(comp_temp->component_decoded==0){
      decoding_left=1;
    }
  }

  //free up buffer
  if(!decoding_left) this->seq_dec_buffer->free_entry(seq_dec_index);
}
