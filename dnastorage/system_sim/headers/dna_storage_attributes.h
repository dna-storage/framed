#ifndef DNA_STORE_ATT
#define DNA_STORE_ATT

#include<stdio.h>





#define WINDOW_SIZE 512
#define QUEUE_SIZE 512

//forward declarations for child classes 
class sequencer_t;
class prep_t;
class decoder_t;
class system_storage_t;


//transaction struct for transactions entered to the system 
typedef struct{
  unsigned long pool_ID;
  unsigned long file_ID;
  unsigned long transaction_ID;
  //conuter that will increment to keep track of the divisions of a request
  unsigned long cracked_count;
  unsigned long decoded;

  //number of strands that will need to be sequenced for this file, will be a function of desired read depth
  unsigned long strands_to_sequence;
  
} transaction_t;

//list entry type used to manage transactions between stages
typedef struct{
  unsigned long transaction_index;
  int used;
} list_entry_t;

//structure used to model a pool containing files
typedef struct{
  unsigned long* files;
  int in_use;
  unsigned long number_files;
} pool_model_t;
 



//variables that describe the overall system
class system_sim_t{
 public:
  //overall system variables
  unsigned long current_data_buffer;
  unsigned long timer_tick;
  unsigned long transactions_completed;
  float sequencing_efficiency;
  //trace pointer is used to access the file that contains the trace
  FILE* trace_pointer;
  //list used for finding and retiring cracked operations
  transaction_t* window[WINDOW_SIZE];
  unsigned long window_head;
  unsigned long window_tail;
  //lists for each component in the system
  sequencer_t* sequencer;
  prep_t* prep;
  decoder_t* decoder;
  system_storage_t* dna_storage;

  //these lists are used to manage the transactions between stages
  list_entry_t prep_seq_list[QUEUE_SIZE];
  list_entry_t seq_dec_list[QUEUE_SIZE];
  
  //number of each unit
  int sequencers;
  int preps;
  int decoders;
  
  system_sim_t(FILE* trace_pointer,unsigned long num_preps, unsigned long  prep_channels, unsigned long  num_sequencers, unsigned long  max_strands_sequencer,unsigned long num_decoders,float seq_efficiency, unsigned long prep_time, unsigned long seq_time, unsigned long dec_time, float seq_eff, unsigned long files_per_pool, unsigned long number_pools, unsigned long average_file_size, unsigned long bytes_per_strand);
  ~system_sim_t();
  
  //this function will be the main function that processes each step in the pipeline
  void simulate(void);
    
};



//parent class for all units in the system
class system_unit_t{
  public:
  unsigned long timer;
  //number of channels that can fit inside the unit
  unsigned long num_channels;
  //numbers of transactions in the current run of a system units
  unsigned long current_num;
  //trasactions mapped to the unit
  unsigned long* transaction_ID;

  // indicates the unit is busy
  int unit_active;

  //constructor
  system_unit_t(unsigned long timer, unsigned long num_channels);
};





//class that models the system storage
class system_storage_t{
 public:
  //parameters of the system
  unsigned int bytes_per_strand;
  unsigned long number_pools;
  float sequencing_efficiency;
  //pool object that holds files
  pool_model_t* pools;
  system_storage_t(float sequencing_efficiency, unsigned long files_per_pool, unsigned long number_pools, unsigned long average_file_size, unsigned long bytes_per_strand);
  ~system_storage_t();
};



#endif
