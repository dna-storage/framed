#ifndef DNA_STORE_ATT
#define DNA_STORE_ATT

#include<stdio.h>





#define WINDOW_SIZE 512
#define QUEUE_SIZE 512

//forward declarations for child classes 
class sequencer;
class prep;
class decoder;
class system_storage;


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
  
} transaction;

//list entry type used to manage transactions between stages
typedef struct{
  unsigned long transaction_index;
  int used;
} list_entry;

//structure used to model a pool containing files
typedef struct{
  unsigned long* files;
  int in_use;
  unsigned long number_files;
} pool_model;
 



//variables that describe the overall system
class system{
 public:
  //overall system variables
  unsigned long current_data_buffer;
  unsigned long timer_tick;
  unsigned long transactions_completed;
  float sequencing_efficiency;
  //trace pointer is used to access the file that contains the trace
  FILE* trace_pointer;
  //list used for finding and retiring cracked operations
  transaction* window[WINDOW_SIZE];
  unsigned long window_head;
  unsigned long window_tail;
  //lists for each component in the system
  sequencer** sequencer_set;
  prep** prep_set;
  decoder** decoder_set;
  system_storage* dna_storage;

  //these lists are used to manage the transactions between stages
  list_entry prep_seq_list[QUEUE_SIZE];
  list_entry seq_dec_list[QUEUE_SIZE];
  
  //number of each unit
  int sequencers;
  int preps;
  int decoders;
  
  system(unsigned long num_preps, unsigned long  prep_channels, unsigned long  num_sequencers, unsigned long  max_strands_sequencer,unsigned long num_decoders,float seq_efficiency, unsigned long prep_time, unsigned long seq_time, unsigned long dec_time);
  ~system();
  
  //this function will be the main function that processes each step in the pipeline
  void simulate(void);
  
 private:
  //function pointers used to make switching policies easier
  typedef void (system::*policy)(void);
  policy sequencing_policy;
  policy prep_selection_policy;
  policy decoder_policy;
  //policies that can be used
  void prep_no_pipeline(void);
  void decoder_no_split(void);
  void seq_no_overlap(void);
  
}



//parent class for all units in the system
class system_unit{
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
  
  //pointer to the system class
  system* system_descriptor;

  //constructor
  system_unit(unsigned long timer, unsigned long num_channels, system* system_descriptor);
}





//class that models the system storage
class system_storage{
 public:
  //parameters of the system
  unsigned int bytes_per_strand;
  unsigned long number_pools;
  float sequencing_efficiency;
  //pool object that holds files
  pool_model* pools;
  system_storage(float sequencing_efficiency, unsigned long files_per_pool, unsigned long number_pools, unsigned long average_file_size, unsigned long bytes_per_strand);
  ~system_storage();
}



#endif
