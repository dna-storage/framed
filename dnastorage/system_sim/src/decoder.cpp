#include "dna_storage_attributes.h"
#include "decoder.h"
#include "prep.h"
#include "sequencer.h"
#include <stdlib.h>
decoder_unit_t::decoder_unit_t(unsigned long timer, unsigned long num_channels) : system_unit_t(timer,num_channels){}

decoder_t::decoder_t(unsigned long timer, unsigned long num_channels, unsigned long buffer_size, unsigned long num_decoders,list_entry_t* seq_dec_buffer,system_sim_t* _system ){
  this->num_decoders=num_decoders;
  this->_system=_system;
  this->seq_dec_buffer=seq_dec_buffer;
  this->buffer_size=buffer_size;

  //allocate the decoder set
  this->decoder_set=(decoder_unit_t**)malloc(sizeof(decoder_unit_t*)*this->num_decoders);
  for(int i=0; i<this->num_decoders; i++) this->decoder_set[i]=new decoder_unit_t(timer,1);

}

//destructor for the decoder system component
decoder_t::~decoder_t(){
  //deallocate items
  for(int i=0; i<this->num_decoders; i++) delete this->decoder_set[i];
  free(this->decoder_set);
}
