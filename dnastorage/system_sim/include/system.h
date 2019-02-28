#ifndef DNA_STORE_SYS
#define DNA_STORE_SYS
#include <stdlib.h>
#include <stdio.h>


//forward declarations for classes 
class sequencer_t;
class prep_t;
class decoder_t;
class system_storage_t;
class scheduler_t;
class generator_t;
class buffer_t;


typedef struct transaction_t{
  unsigned long pool_ID; //pool_ID of the original transaction 
  unsigned long cracked_count; //keeps track of how many times the transaction was split up due to sequecer space limitations
  unsigned long strands_to_sequence;// total strands to be sequenced: read_depth*(file_size/bytes_per_strand)+undesired_strands
  unsigned long undesired_strands_sequenced; //strands that are junk: (1/eff.-1)*desired
  unsigned long desired_strands_sequenced; //strands that we want: read_depth*(file_size/bytes_per_strand)
  int transaction_finished; //flag used to indicate if the transaction in the list is finished or not
  unsigned long digital_data_size; //digital data size
  unsigned long time_stamp; //indicates the time a transaction or component became first available to the system
  int component_decoded; //inidicates if the component has been moved to decoding
  struct transaction_t* components; // components of a single large transaction created at the scheduling point, this member is the head of a linked list of components for a transaction
  struct transaction_t* next;
} transaction_t; //represents a transaction in the pipeline





typedef struct trace_t{
  unsigned long pool_ID; //ID of the pool to be accessed
  unsigned long time_stamp; //time at which the transaction became available to the system
  unsigned long file_size; //size of the file that is being accessed (in KB)
  struct trace_t* next;//this structure is going to be used in a linked list
} trace_t; //structure that defines the transactions in the system queue


typedef struct{
  unsigned long num_preps;
  unsigned long prep_channels;
  unsigned long num_sequencers;
  unsigned long max_strands_sequencer;
  unsigned long seq_channels;
  unsigned long num_decoders;
  float seq_efficiency;
  unsigned long prep_time;
  unsigned long seq_time;
  unsigned long dec_time;
  unsigned long average_pool_capacity;
  unsigned long number_pools;
  float bytes_per_strand;
  unsigned long pool_write_time;
  unsigned long pool_wait_time;
  unsigned long pool_copies;
  unsigned long number_reads;
  unsigned long sequencer_timeout;
  FILE* trace_file;
  unsigned long window_size;
  unsigned long prep_seq_buffer_size;
  unsigned long seq_dec_buffer_size;
  unsigned long batch_size;
  unsigned long max_file_size;
  unsigned long min_file_size;
  float rate;
  int seed;
  unsigned long sequencing_depth;
  
} system_sim_params_t; //bundled system parameters


//variables that describe the overall system
class system_sim_t{
 public:
  transaction_t* window; //window of transactions allowed into the system
  
   system_sim_t(system_sim_params_t system_sim_params);
  ~system_sim_t();
  void queue_pop(void); //take a request from the head of the system queue
  void queue_append(trace_t* new_trans); //append a request to the system queue
  void simulate(void); //top level call for main
  void cleanup_active_list(void); //cleanup the window
  int window_full(void); //check to see if the window is full
  void window_pop(void); //pop a transaction from the window
  unsigned long window_add(void); //add a transaction to the window
  void window_componentadd(unsigned long transaciton_ID,
			   unsigned long undesired_strands_sequenced,
			   unsigned long desired_strands_sequenced,
			   unsigned long pool_ID,
			   unsigned long digital_data_size); //add a component to the transaction_ID
  void window_init(unsigned long transaction_ID); // initialize the window entry at transaction_ID

  unsigned long timer_tick; //current time of the simulator

 private:
  unsigned long sim_time; //denotes the maximum amount of time the simulator should run for
  trace_t* system_queue; //system_queue used at the front end of the system for scheduling
  unsigned long window_size; //max size of the window
  unsigned long window_length; //number of active items in the window
  unsigned long window_head; //head pointer for the window
  unsigned long window_tail; //tail pointer for the window
  sequencer_t* sequencer; //pointer to the sequencer stage
  prep_t* prep; //pointer to the prep stage
  decoder_t* decoder; //pointer to the decode stage
  system_storage_t* dna_storage; //pointer to the dna storage unit
  scheduler_t* scheduler; //pointer to the scheduler
  generator_t* generator; //pointer to the generator 
 
  buffer_t* buffers[2]; //array of buffers we create
 


};




class system_unit_t{ //decoder,sequencer,prep unit types are all derived from this class
  public:
  unsigned long timer; //timer that indicates how much time is left for the unit
  unsigned long num_channels; //number of separate transactions on the unit
  unsigned long next_open; //next spot open in transaction_slots 
  unsigned long* transaction_slots; //array of transactions on the unit
  unsigned long timeout; //timeout to delay the start of a unit after transactions have been put into it 
  unsigned long transaction_pointer;//used to empty transaction_slots
  int unit_active; //1-> unit is busy, 0 -> unit is not busy
  system_unit_t(unsigned long num_channels);
  ~system_unit_t();
};





#endif
