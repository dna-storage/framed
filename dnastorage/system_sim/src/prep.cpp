#include "dna_storage_attributes.h"
#include "decoder.h"
#include "prep.h"
#include "sequencer.h"
prep_unit_t::prep_unit_t(unsigned long timer, unsigned long num_channels, system_sim* system_descriptor):system_unit(timer,num_channels,system_descriptor){}

//constructor for the prep component of the system
prep_t::prep_t(unsigned_long timer, unsigned long num_channels, unsigned long buffer_size,unsigned long num_preps,list_entry_t* prep_seq_buffer,system_sim_t* _system,system_storage_t* dna_storage){
  
  this->num_preps=num_preps;
  this->_system=_system;
  this->prep_seq_buffer=prep_seq_buffer;
  this->buffer_size=buffer_size;
  this->dna_storage=dna_storage;

  //make the list of prep units
  this->prep_set=(prep_unit_t**)malloc(sizeof(prep_unit_*)*this->preps);
  for(int i=0; i<this->num_preps; i++) this->prep_set[i]=new prep_unit_t(timer,num_channels);

  //assign the prep policy here

}


prep_t::~prep_t(){
  //free up the prep list
  for(int i=0; i<this->_numpreps; i++) delete this->prep_set[i];
  free(this->prep_set);


}
