#include "dna_storage_attributes.h"
#include "decoder.h"
#include "prep.h"
#include "sequencer.h"
#include <stdlib.h>
decoder_unit_t::decoder_unit_t(unsigned long timer, unsigned long num_channels) : system_unit_t(num_channels){}

decoder_t::decoder_t(unsigned long timer, unsigned long num_channels,
		     unsigned long buffer_size, unsigned long num_decoders,
		     list_entry_t* seq_dec_buffer,system_sim_t* _system )
{
  
  this->num_decoders=num_decoders;
  this->_system=_system;
  this->seq_dec_buffer=seq_dec_buffer;
  this->buffer_size=buffer_size;
  this->base_timer=timer;
  
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

//function that looks through the decoders and see which transactions are finished
unsigned long decoder_t::decoder_backend(void){
  unsigned long done_count=0;
  //iterate through the decoders and look for complete jobs
  for(unsigned long i=0; i<this->num_decoders;i++){
    int active=(this->decoder_set[i])->unit_active;
    unsigned long timer=(this->decoder_set[i])->timer;
    if(timer==0 && active){
      this->decoder_complete(i); //complete the transaction within the decoder
      done_count++;
    }
    else this->decoder_timestep(i); //decrease the time on the decoder
  }
  
  return done_count;
}

//need to fill up non-active decoders with transactions in the seq_dec_buffer
void decoder_t::decoder_frontend(void){
  for(unsigned long i=0; i<(this->buffer_size);i++){
    if(this->seq_dec_buffer[i].used && (this->_system)->window[this->seq_dec_buffer[i].transaction_index].cracked_count==0){
      //found a location in the buffer that is ready
      unsigned long decoder_ID=this->decoder_avail();
      unsigned long transaction_ID=this->seq_dec_buffer[i].transaction_index;
      if(decoder_ID==-1) break; //no decoder avilable, stop looking through the seq_dec_buffer
      this->init_decoder(decoder_ID,transaction_ID); //initialize the decoder found thats available
    }
  }
}



//reach into the decoder set and deactivate the decoder and complete the transaction
void decoder_t::decoder_complete(unsigned long decoder_ID){
  (this->decoder_set[decoder_ID])->unit_active=0;
  for(unsigned long i=0; i<next_open;i++) (this->_system)->window[(this->decoder_set[decoder_ID])->transaction_slots[i]].transaction_finished=1;
}

void decoder_t::decoder_timestep(unsigned long decoder_ID){
  (this->decoder_set[decoder_ID])->timer--;
}

//search for deactivated decoders
unsigned long decoder_t::decoder_avail(void){
  for(unsigned long i=0; i<(this->num_decoders);i++){
    if ((this->decoder_set[i])->unit_active==0){
      return i;
    }
  }
  return -1;
}

//initialize the chosen decoder unit with timer, active signals etc
void decoder_t::init_decoder(unsigned long decoder_ID, unsigned long transaction_ID){
  (this->decoder_set[decoder_ID])->next_open=0;
  (this->decoder_set[decoder_ID])->transaction_slots[next_open]=transaction_ID;
  (this->decoder_set[decoder_ID])->unit_active=1;
  (this->decoder_set[decoder_ID])->timer=(this->base_timer);
}
