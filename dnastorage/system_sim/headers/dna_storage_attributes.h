#ifndef DNA_STORE_ATT
#define DNA_STORE_ATT

#define WINDOW_SIZE 512
#define QUEUE_SIZE 512
#define TRACE_BUFFER_SIZE 512

//forward declarations for classes 
class sequencer_t;
class prep_t;
class decoder_t;
class system_storage_t;


//transaction struct for transactions entered to the system 
typedef struct{
  unsigned long pool_ID; //pool_ID of the original transaction
  //conuter that will increment to keep track of the divisions of a request, should be 0 after all pieces are sequenced 
  unsigned long cracked_count;
  //number of strands that will need to be sequenced for this file, will be a function of desired read depth
  unsigned long strands_to_sequence;// total strands to be sequenced: read_depth*(file_size/bytes_per_strand)+undesired_strands
  unsigned long undesired_strands_sequenced; //strands that are junk: (1/eff.-1)*desired
  unsigned long desired_strands_sequenced; //strands that we want: read_depth*(file_size/bytes_per_strand)

  int transaction_finished; //flag used to indicate if the transaction in the list is finished or not
  
} transaction_t;

//list entry type used to manage transactions between stages
typedef struct{
  unsigned long transaction_index;
  int used;
} list_entry_t;


typedef struct{
  unsigned long timer; //timer keeping track of whether the pool is usable
  unsigned long remaining_reads; //keeps track of the remaining reads of this pool
  unsigned long in_use; //flag that indicates if pool can be used, make sure back to back uses are blocked
} pool_char_t;

//structure used to model a pool containing files
typedef struct{
  pool_char_t* copies;
} pool_model_t;


typedef struct trace_t{
  unsigned long pool_ID; //ID of the pool to be accessed
  unsigned long time_stamp; //time at which the transaction will be available to the system
  unsigned long file_size; //size of the file that is being accessed (in MB)
  struct trace_t* next;//this structure is going to be used in a linked list
} trace_t;


//variables that describe the overall system
class system_sim_t{
 public:
  //overall system variables
  unsigned long current_data_buffer;
  unsigned long timer_tick;
  unsigned long transactions_completed;
  float sequencing_efficiency;
  //trace pointer is used to access the file that contains the trace
  FILE* trace_file;
  trace_t* trace_list_head;
  //list used for finding and retiring cracked operations
  transaction_t window[WINDOW_SIZE];
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
  int trace_exhuasted;//indicates that the trace file is exhausted
  int trace_complete; //indicates that all of the trace transactions have been put into the pipeline
  
  system_sim_t(FILE* trace_file,unsigned long num_preps, unsigned long  prep_channels,
	       unsigned long num_sequencers, unsigned long  max_strands_sequencer,
	       unsigned long num_decoders,float seq_efficiency, unsigned long prep_time,
	       unsigned long seq_time, unsigned long dec_time, float seq_eff,
	       unsigned long average_pool_capacity, unsigned long number_pools,
	       unsigned long average_file_size, unsigned long bytes_per_strand,
	       unsigned long pool_write_time,unsigned long pool_wait_time,
	       unsigned long pool_copies, unsigned long number_reads);
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
  //next open slot in the channels list, serves also as a way to determine the number of transactions in the certain unit
  unsigned long next_open; //used to fill up transaction_slots
  //trasactions mapped to the unit
  unsigned long* transaction_slots;

  unsigned long transaction_pointer;//used to empty transaction_slots
  // indicates the unit is busy
  int unit_active;

  unsigned long standby_timer;
  
  //constructor
  system_unit_t(unsigned long num_channels);
  ~system_unit_t();
};





//class that models the system storage
class system_storage_t{
 public:
  //parameters of the system
  unsigned long bytes_per_strand; //used to model density
  unsigned long number_pools;
  unsigned long average_pool_capacity; //average capacity of the pool (MB), need to model whole pool sequencing 
  float sequencing_efficiency; //models how well a certain 
  //parameters of pool re-use
  unsigned long pool_copies;//number of copies the pool has
  unsigned long number_reads;//number of reads that a pool can withstand 
  unsigned long pool_write_time;//time it takes to replenish pool after using all of its reads
  unsigned long pool_wait_time;//time it takes for a pool to be available after using it
  //pool object that holds files
  pool_model_t* pools;
  system_storage_t(float sequencing_efficiency, unsigned long average_pool_capacity,
		   unsigned long number_pools,
		   unsigned long bytes_per_strand,unsigned long pool_copies,
		   unsigned long number_reads, unsigned long pool_write_time,
		   unsigned long pool_wait_time);
  ~system_storage_t();
  int storage_poolavailable(unsigned long pool_ID);
  void storage_readmanage(unsigned long pool_ID, unsigned long copy_ID);
  void storage_poolrestore(void);
};



#endif
