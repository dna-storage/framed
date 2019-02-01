#include <stdlib.h>
#include <stdio.h>
#include "dna_storage_attributes.h"

//constructor for the whole system
system::system(FILE* trace_pointer,unsigned long num_preps, unsigned long  prep_channels, unsigned long  num_sequencers, unsigned long  max_strands_sequencer,unsigned long num_decoders,float seq_efficiency, unsigned long prep_time, unsigned long seq_time, unsigned long dec_time, float seq_eff, unsigned long files_per_pool, unsigned long number_pools, unsigned long average_file_size, unsigned long bytes_per_strand){
  //initialize the system class members
  this->timer_tick=0;
  this->current_data_buffer=0;
  this->transactions_completed=0;
  this->window_head=0;
  this->window_tail=0;
  this->sequencing_efficiency=seq_eff;
  this->trace_pointer=trace_pointer;
  //initialize lists
  for(int i=0; i<QUEUE_SIZE; i++){
    this->prep_seq_list[i].used=0;
    this->seq_dec_list[i].used=0;
  }

  //will need to set the function pointers to different policies 
  this->sequencing_policy=&system::seq_no_overlap;
  this->decoder_policy=&system::decoder_no_split;
  this->&system::prep_no_pipeline;
  //////////////////////////////////////////////////////////////
  this->sequencers=num_sequencers;
  this->preps=num_preps;
  this->decoders=num_decoders;

  //instantiate the system units
  this->sequencer_set=(sequencer**)malloc(sizeof(sequencer*)*this->sequencers);
  this->prep_set=(prep**)malloc(sizeof(prep*)*this->preps);
  this->decoder_set=(decoder**)malloc(sizeof(decoder*)*this->decoders);
  for(int i=0; i<this->sequencers; i++) this->sequencer_set[i]=new sequencer(seq_time,max_strands_sequencer,0,this);
  for(int i=0; i<this->preps; i++) this->prep_set[i]=new prep(prep_time,prep_channels,this);
  for(int i=0; i<this->decoders; i++) this->decoder_set[i]=new decoder(dec_time,1,this);


  //instantiate the storage model
  dna_storage=new system_storage(seq_eff,files_per_pool,number_pools,average_file_size,bytes_per_strand);

  
}

//destructor for the system class
system::~system(){
  //delete the previously allocated system units
  for(int i=0; i<this->sequencers; i++) delete this->sequencer_set[i];
  for(int i=0; i<this->preps; i++) delete this->prep_set[i];
  for(int i=0; i<this->decoders; i++) delete this->decoder_set[i];
  //now free the array of pointers previously allocated
  free(this->sequencer_set);
  free(this->prep_set);
  free(this->decoder_set);

  //delete the dna_storage member
  delete dna_storage;

}



//constructor for the system_unit parent class
system_unit::system_unit(unsigned long timer, unsigned long num_channels, system* system_descriptor){
  this->timer=timer;
  this->num_channels=num_channels;
  this->unit_active=0;
  this->system_descriptor=system_descriptor;
}

system_storage:: system_storage(float sequencing_efficiency, unsigned long files_per_pool, unsigned long number_pools, unsigned long average_file_size, unsigned long bytes_per_strand){
  this->sequencing_efficiency=sequencing_efficiency;
  this->number_pools=number_pools;
  this->bytes_per_strand=bytes_per_strand;
  this->pools=(pool_model*)malloc(sizeof(pool_model)*this->number_pools);
  for(int i=0; i<this->number_pools,i++){
    this->pools[i].in_use=0;
    this->pools[i].number_files=files_per_pool;
    this->pools[i].files=(int*)malloc(sizeof(int)*files_per_pool);
  }

}

//free up space allocated for pools/files
system_storage::~system(storage){
    for(int i=0; i<this->number_pools,i++) free(this->pools[i].files);
    free(this->pools);
}


system::simulate(){

}
