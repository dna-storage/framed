#include <stdlib.h>
#include <stdio.h>
#include "dna_storage_attributes.h"
#include "decoder.h"
#include "prep.h"
#include "sequencer.h"
#include "stdio.h"

//constructor for the whole system
system_sim_t::system_sim_t(FILE* trace_pointer,unsigned long num_preps,
			   unsigned long  prep_channels, unsigned long  num_sequencers,
			   unsigned long  max_strands_sequencer,unsigned long num_decoders,
			   float seq_efficiency, unsigned long prep_time, unsigned long seq_time,
			   unsigned long dec_time, float seq_eff, unsigned long files_per_pool,
			   unsigned long number_pools, unsigned long average_file_size,
			   unsigned long bytes_per_strand, unsigned long pool_write_time,
			   unsigned long pool_wait_time, unsigned long pool_copies,
			   unsigned long number_reads){
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

  //////////////////////////////////////////////////////////////
  this->sequencers=num_sequencers;
  this->preps=num_preps;
  this->decoders=num_decoders;

  //instantiate the storage model
  this->dna_storage=new system_storage_t(seq_eff,files_per_pool,number_pools,average_file_size,bytes_per_strand,pool_copies,number_reads,pool_write_time,pool_wait_time);
  //instantiate system components
  this->decoder=new decoder_t(dec_time,1,QUEUE_SIZE,this->decoders,seq_dec_list,this);
  this->prep=new prep_t(prep_time,prep_channels,QUEUE_SIZE,this->preps,prep_seq_list,this,this->dna_storage);
  this->sequencer=new sequencer_t(seq_time,max_strands_sequencer,0,QUEUE_SIZE,this->sequencers,seq_dec_list,this->prep_seq_list,this);
}

//destructor for the system class
system_sim_t::~system_sim_t(){
  delete this->decoder;
  delete this->prep;
  delete this->sequencer;
  //delete the dna_storage member
  delete this->dna_storage;
  
}



//constructor for the system_unit parent class
system_unit_t::system_unit_t(unsigned long timer, unsigned long num_channels){
  this->timer=timer;
  this->num_channels=num_channels;
  this->unit_active=0;
  this->next_open=0;
  //allocate an array for each system unit
  this->transaction_slots=(unsigned long*)malloc(sizeof(unsigned long)*this->num_channels);
}

system_unit_t::~system_unit_t(){
  free(transaction_slots);
  printf("Destructor of parent units\n");
}

system_storage_t:: system_storage_t(float sequencing_efficiency, unsigned long files_per_pool,
				    unsigned long number_pools, unsigned long average_file_size,
				    unsigned long bytes_per_strand,unsigned long pool_copies,
				    unsigned long number_reads,unsigned long pool_write_time,
				    unsigned long pool_wait_time){
  this->sequencing_efficiency=sequencing_efficiency;
  this->number_pools=number_pools;
  this->bytes_per_strand=bytes_per_strand;
  this->pools=(pool_model_t*)malloc(sizeof(pool_model_t)*this->number_pools);
  this->pool_wait_time=pool_wait_time;
  this->pool_copies=pool_copies;
  this->number_reads=number_reads;
  this->pool_write_time=pool_write_time;
  for(int i=0; i<this->number_pools;i++){
    this->pools[i].in_use=0;
    this->pools[i].number_files=files_per_pool;
    this->pools[i].files=(unsigned long*)malloc(sizeof(unsigned long)*files_per_pool);
  }

}

//free up space allocated for pools/files
system_storage_t::~system_storage_t(){
  for(int i=0; i<this->number_pools;i++) free(this->pools[i].files);
  free(this->pools);
}



//this is the top level simulator, and calls the different system unit functions
void system_sim_t::simulate(){

  
  //clean up the active transaction list
  this->cleanup_active_list();
  //add on the number decoders done, decoder_backend will relinquish decoders provide a
  this->transactions_completed+=(this->decoder)->decoder_backend();
  (this->decoder)->decoder_frontend();
  (this->sequencer)->sequencer_backend();
  (this->sequencer)->sequencer_frontend();
  (this->prep)->prep_backend();
  //FIX ME: Need to be able to create a list of the access trace in order to allow the prep_frontend function to traverse it and find possible jobs to run
  (this->prep)->prep_frontend();
  //FIX ME:add in some support for the dna storage system structure here
  
  //FIX ME:need to increment time_step
  


  
}

 
//look through the transactions in window[], and clean it up by advancing the head pointer
void system_sim_t::cleanup_active_list(void){
  unsigned long start_point=this->window_head;
  while(start_point!=this->window_tail){
    if(this->window[start_point].transaction_finished){
      this->window_head++;
      start_point++;
    }
    else break;
  }
}







