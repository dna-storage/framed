#include "dna_storage_attributes.h"
#include "decoder.h"
#include "prep.h"
#include "sequencer.h"
#include <stdlib.h>
sequencer_unit_t::sequencer_unit_t(unsigned long timer,unsigned long max_sequencing,
				   unsigned long num_channels):system_unit_t(num_channels){
  this->max_strands=max_sequencing;
  this->used_strands=0;
  this->wasted_strands=0;
  this->utilization=0;
}


sequencer_t:: sequencer_t(unsigned long timer,unsigned long max_strands,
			  unsigned long num_channels,unsigned long buffer_size,
			  unsigned long num_sequencers,list_entry_t* seq_dec_buffer,
			  list_entry_t* prep_seq_buffer, system_sim_t* _system){

  
  this->num_sequencers=num_sequencers;
  this->_system=_system;
  this->seq_dec_buffer=seq_dec_buffer;
  this->prep_seq_buffer=prep_seq_buffer;
  this->buffer_size=buffer_size;
  this->base_timer=timer;
  //allocate sequencer units
  this->sequencer_set=(sequencer_unit_t**)malloc(sizeof(sequencer_unit_t*)*this->num_sequencers);
  for(int i=0; i<this->num_sequencers; i++) this->sequencer_set[i]=new sequencer_unit_t(max_strands,0);
 
}


sequencer_t::~sequencer_t(){
  for(int i=0; i<this->num_sequencers; i++) delete this->sequencer_set[i];
  free(this->sequencer_set);
}
