#include <stdlib.h>
#include <stdio.h>
#include "dna_storage_attributes.h"
#include "decoder.h"
#include "prep.h"
#include "sequencer.h"
#include "stdio.h"
#include "utlist.h"

//constructor for the whole system
system_sim_t::system_sim_t(FILE* trace_file,unsigned long num_preps,
			   unsigned long  prep_channels, unsigned long  num_sequencers,
			   unsigned long  max_strands_sequencer,unsigned long num_decoders,
			   float seq_efficiency, unsigned long prep_time, unsigned long seq_time,
			   unsigned long dec_time, float seq_eff,
			   unsigned long average_pool_capacity,
			   unsigned long number_pools, 
			   unsigned long bytes_per_strand, unsigned long pool_write_time,
			   unsigned long pool_wait_time, unsigned long pool_copies,
			   unsigned long number_reads){
  //initialize the system class members
  this->timer_tick=0;
  this->current_data_buffer=0;
  this->transactions_completed=0;
  this->window_head=0;
  this->window_tail=0;
  this->window_empty=1;
  this->sequencing_efficiency=seq_eff;
  this->trace_file=trace_file;
  this->queue_head=NULL;
  this->trace_exhausted=0;
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
  this->dna_storage=new system_storage_t(seq_eff,average_pool_capacity,number_pools,bytes_per_strand,pool_copies,number_reads,pool_write_time,pool_wait_time);
  //instantiate system components
  this->decoder=new decoder_t(dec_time,1,QUEUE_SIZE,this->decoders,seq_dec_list,this);
  this->prep=new prep_t(prep_time,prep_channels,QUEUE_SIZE,this->preps,prep_seq_list,this,this->dna_storage);
  this->sequencer=new sequencer_t(seq_time,max_strands_sequencer,0,QUEUE_SIZE,this->sequencers,seq_dec_list,this->prep_seq_list,this);

  this->trace_exhausted=this->add_traces(TRACE_BUFFER_SIZE); //start up the trace list

  
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
  this->transaction_pointer=0;
}

system_unit_t::~system_unit_t(){
  free(transaction_slots);
  printf("Destructor of parent units\n");
}

system_storage_t:: system_storage_t(float sequencing_efficiency,
				    unsigned long average_pool_capacity,
				    unsigned long number_pools,
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
  this->average_pool_capacity=average_pool_capacity;
  for(int i=0; i<this->number_pools;i++){
    this->pools[i].copies=(pool_char_t*)malloc(sizeof(pool_char_t)*this->pool_copies);
    //need to initialize the array of pool copies
    for(int j=0; j<this->pool_copies; j++){
      this->pools[i].copies[j].timer=0;
      this->pools[i].copies[j].remaining_reads=this->number_reads;
      this->pools[i].copies[j].in_use=0;
    }
  }
}

//free up space allocated for pools/files
system_storage_t::~system_storage_t(){
  for(int i=0; i<this->number_pools;i++) free(this->pools[i].copies);
  free(this->pools);
}


//function that manages the head pointer when adding items
unsigned long system_sim_t::window_add(void){
  unsigned long out= this->window_tail;
  this->window_tail++;
  if(this->window_tail>=WINDOW_SIZE) this->window_tail=0;
  if(this->window_head==this->window_tail) this->window_empty=0;
  return out;
}

//function that creates and inserts a component into a window element
void system_sim_t::window_componentadd(unsigned long transaction_ID,
				       unsigned long undesired_strands_sequenced,
				       unsigned long desired_strands_sequenced,
				       unsigned long pool_ID,
				       unsigned long digital_data_size){
  transaction_t* new_component;
  new_component=(transaction_t*)malloc(sizeof(transaction_t)); //allocate the space
  new_component->next=NULL;
  new_component->undesired_strands_sequenced=undesired_strands_sequenced;
  new_component->desired_strands_sequenced=desired_strands_sequenced;
  new_component->strands_to_sequence=desired_strands_sequenced+undesired_strands_sequenced;
  new_component->pool_ID=pool_ID;
  new_component->component_decoded=0;
  LL_APPEND(this->window[transaction_ID].components,new_component); //append the component to the linked list
  //accumulate component values into the overall transaction
  this->window[transaction_ID].strands_to_sequence+=new_component->strands_to_sequence;
  this->window[transaction_ID].digital_data_size+=digital_data_size;
  this->window[transaction_ID].undesired_strands_sequenced+=undesired_strands_sequenced;
  this->window[transaction_ID].desired_strands_sequenced+=desired_strands_sequenced;
}


//initialize a entry in the window
void system_sim_t::window_init(unsigned long transaction_ID){
  this->window[transaction_ID].cracked_count=1;
  this->window[transaction_ID].strands_to_sequencer=0;
  this->window[transaction_ID].undesired_strands_sequenced=0;
  this->window[transaction_ID].desired_strands_sequenced=0;
  this->window[transaction_ID]. component_decoded=0;
  this->window[transaction_ID].transaction_finished=0;
  this->window[transaction_ID].components=NULL;
  this->window[transaction_ID].next=NULL;
  this->window[transaction_ID].digital_data_size=0;
}


//pop a request off of the queue head
void system_sim_t::queue_pop(void){
  trace_t* temp;
  temp=this->queue_head;
  LL_DELETE(this->queue_head,this->queue_head);
  free(temp);
}

//add new_trans to the end of the queue
void system_sim_t::queue_append(trace_t* new_trans){
  LL_APPEND(this->queue_head,new_trans);
}


//this is the top level simulator, and calls the different system unit functions
void system_sim_t::simulate(){

  //run the simulator for the duration of the simulation time
  while(this->timer_tick<=this->sim_time){
    //clean up the active transaction list
    this->cleanup_active_list();
    this->decoder->decoder_stage();
    this->sequencer->sequencer_stage();
    this->prep->prep_stage();
    this->scheduler->schedule_stage();
    this->generator->generator_stage();
    this->timer_tick++;
  }
}

 
//look through the transactions in window[], and clean it up by advancing the head pointer
void system_sim_t::cleanup_active_list(void){
  unsigned long start_point=this->window_head;
  transaction_t* component_head;
  transaction_t* component_temp1;
  transaction_t* component_temp2;
  while(start_point!=this->window_tail){
    if(this->window[start_point].transaction_finished){
      //deallocate space for the components of the transactions
      LL_FOREACH_SAFE(component_head,component_temp1,component_temp2){
	LL_DELETE(component_head,component_head);
	free(component_head);
      }
      this->window_head++;
      start_point++;
    }
    else break;
  }
}


