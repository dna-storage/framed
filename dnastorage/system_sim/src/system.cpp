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
#include <stdio.h>


//constructor for the whole system
system_sim_t::system_sim_t(system_sim_params_t system_sim_params){
  prep_params_t prep_params;
  sequencer_params_t sequencer_params;
  storage_params_t storage_params;
  scheduler_params_t scheduler_params;
  generator_params_t generator_params;
  decoder_params_t decoder_params;
 
  
  //initialize the system class members
  this->timer_tick=0;
  this->window_size=system_sim_params.window_size;
  this->window_head=0;
  this->window_tail=0;
  this->window_length=0;
  this->system_queue=NULL;
 
  //create the window
  this->window_size=system_sim_params.window_size;  
  this->window=(transaction_t*)malloc((this->window_size+1)*sizeof(transaction_t));
  
  this->buffers[0]=new buffer_t(system_sim_params.prep_seq_buffer_size);
  this->buffers[1]=new buffer_t(system_sim_params.seq_dec_buffer_size);
  this->stats=new stats_t(system_sim_params.stats_log,system_sim_params.phase_log);
  
  //setup storage parameters
  storage_params.sequencing_efficiency=system_sim_params.seq_efficiency;
  storage_params.average_pool_capacity=system_sim_params.average_pool_capacity;
  storage_params.number_pools=system_sim_params.number_pools;
  storage_params.bytes_per_strand=system_sim_params.bytes_per_strand;
  storage_params.pool_copies=system_sim_params.pool_copies;
  storage_params.pool_write_time=system_sim_params.pool_write_time;
  storage_params.pool_wait_time=system_sim_params.pool_wait_time;
  storage_params.stats=this->stats;

  //instantiate the storage model
  this->dna_storage=new system_storage_t(storage_params);
  
  //setup prep parameters
  prep_params.timer=system_sim_params.prep_time;
  prep_params.num_channels=system_sim_params.prep_channels;
  prep_params.prep_seq_buffer=this->buffers[0];
  prep_params._system=this;
  prep_params.dna_storage=this->dna_storage;
  prep_params.num_preps=system_sim_params.num_preps;
  prep_params.stats=this->stats;
  
  //setup decoder parameters
  decoder_params.timer=system_sim_params.dec_time;
  decoder_params.num_channels=1;
  decoder_params.seq_dec_buffer=this->buffers[1];
  decoder_params.num_decoders=system_sim_params.num_decoders;
  decoder_params.stats=this->stats;
  
  //setup sequencer parameters
  sequencer_params.timer=system_sim_params.seq_time;
  sequencer_params.max_strands=system_sim_params.max_strands_sequencer;
  sequencer_params.num_channels=system_sim_params.seq_channels;
  sequencer_params.base_timeout=system_sim_params.sequencer_timeout;
  sequencer_params.num_sequencers=system_sim_params.num_sequencers;
  sequencer_params.seq_dec_buffer=this->buffers[1];
  sequencer_params.prep_seq_buffer=this->buffers[0];;
  sequencer_params._system=this;
  sequencer_params.stats=this->stats;

  //instantiate the decoder, prep and sequencer
  this->decoder = new decoder_t(decoder_params);
  this->prep = new prep_t(prep_params);
  this->sequencer = new sequencer_t(sequencer_params);


  //setup scheduler parameters
  scheduler_params._storage=this->dna_storage;
  scheduler_params.system_queue=&(this->system_queue);
  scheduler_params._prep=this->prep;
  scheduler_params._system=this;
  scheduler_params.batch_size=system_sim_params.batch_size;
  scheduler_params.bytes_per_strand=system_sim_params.bytes_per_strand;
  scheduler_params.sequencing_depth=system_sim_params.sequencing_depth;
  scheduler_params.efficiency=system_sim_params.seq_efficiency;
  scheduler_params.stats=this->stats;

  //setup the generator parameters
  generator_params.max_file_size=system_sim_params.max_file_size;
  generator_params.min_file_size=system_sim_params.min_file_size;
  generator_params.unique_pools=system_sim_params.number_pools;
  generator_params.random_seed=system_sim_params.seed;
  generator_params._system=this;
  generator_params.rate=system_sim_params.rate;
  generator_params.stats=this->stats;
  
  this->scheduler = new scheduler_t(scheduler_params);
  this->generator = new generator_t(generator_params);
}


//destructor for the system class
system_sim_t::~system_sim_t(){
  delete this->decoder;
  delete this->prep;
  delete this->sequencer;
  //delete the dna_storage member
  delete this->dna_storage;
  delete this->generator;
  delete this->scheduler;
  free(this->window);
  for(int i=0;i<2;i++){
    delete this->buffers[i]; //free the created buffers
  }
}



//constructor for the system_unit parent class
system_unit_t::system_unit_t(unsigned long num_channels){
  this->timer=0;
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



//pop a request off of the queue head
void system_sim_t::queue_pop(void){
  trace_t* temp;
  temp=this->system_queue;
  LL_DELETE(this->system_queue,this->system_queue);
  free(temp);
}

//add new_trans to the end of the queue
void system_sim_t::queue_append(trace_t* new_trans){
  LL_APPEND(this->system_queue,new_trans);
}


//this is the top level simulator, and calls the different system unit functions
void system_sim_t::simulate(){

  //run the simulator for the duration of the simulation time, and until both the system_queue and window are empty
  while(this->timer_tick<=this->sim_time && this->system_queue!=NULL && !this->window_empty()){
    //clean up the active transaction list
    this->cleanup_active_list();
    this->decoder->decoder_stage();
    this->sequencer->sequencer_stage();
    this->prep->prep_stage();
    this->scheduler->schedule_stage();
    this->generator->generator_stage();
    this->timer_tick++;
    inc_counter(time_step);
  }
}

 
//look through the transactions in window[], and clean it up by advancing the head pointer
void system_sim_t::cleanup_active_list(void){
  unsigned long start_point=this->window_head;
  transaction_t* component_head;
  transaction_t* component_temp1;
  transaction_t* component_temp2;
  while(start_point!=this->window_length){
    if(this->window[start_point].transaction_finished){
      //deallocate space for the components of the transactions
      LL_FOREACH_SAFE(component_head,component_temp1,component_temp2){
	add_counter(total_latency,(component_head->time_stamp_end-component_head->time_stamp_start));
	LL_DELETE(component_head,component_head);
	free(component_head);
	inc_counter(finished_requests); //increment the number of requests finished
      }
      add_counter(data_decoded,this->window[start_point].digital_data_size); //add to the total data decoded
      this->window_pop();
      start_point++;
    }
    else break;
  }
}



//check to see if the window is full
int system_sim_t::window_full(void){
  unsigned long _tail;
  _tail=(this->window_tail+1)%(this->window_size+1);
  return _tail==this->window_head; //check whether the head is the next position if so it is full
}

//remove element from the window
void system_sim_t::window_pop(void){
  this->window_head=(this->window_head+1)%(this->window_size+1);
  this->window_length--;
}

//function that manages the tail pointer when adding items
unsigned long system_sim_t::window_add(void){
  unsigned long out= this->window_tail;
  this->window_tail=(this->window_tail+1)%(this->window_size+1);
  this->window_length++;
  return out;
}

//return 1 if window is empty
int system_sim_t::window_empty(void){
  return this->window_head==this->window_tail;
}


//function that creates and inserts a component into a window element
void system_sim_t::window_componentadd(unsigned long transaction_ID,
				       unsigned long undesired_strands_sequenced,
				       unsigned long desired_strands_sequenced,
				       unsigned long pool_ID,
				       unsigned long digital_data_size,
				       unsigned long time_stamp){
  transaction_t* new_component;
  new_component=(transaction_t*)malloc(sizeof(transaction_t)); //allocate the space
  new_component->next=NULL;
  new_component->undesired_strands_sequenced=undesired_strands_sequenced;
  new_component->desired_strands_sequenced=desired_strands_sequenced;
  new_component->strands_to_sequence=desired_strands_sequenced+undesired_strands_sequenced;
  new_component->pool_ID=pool_ID;
  new_component->component_decoded=0;
  new_component->time_stamp_start=time_stamp;
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
  this->window[transaction_ID].strands_to_sequence=0;
  this->window[transaction_ID].undesired_strands_sequenced=0;
  this->window[transaction_ID].desired_strands_sequenced=0;
  this->window[transaction_ID].transaction_finished=0;
  this->window[transaction_ID].components=NULL;
  this->window[transaction_ID].next=NULL;
  this->window[transaction_ID].digital_data_size=0;
}
