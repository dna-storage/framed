#ifndef SCHEDULER
#define SCHEDULER

class system_storage_t;
class prep_t;
class system_sim_t;
struct transaction_t;
struct trace_t;

typedef struct{
  system_storage_t* _storage;
  struct trace_t** system_queue;
  prep_t* _prep;
  system_sim_t* _system;
  unsigned long batch_size;
  float bytes_per_strand;
  unsigned long sequencing_depth;
  float efficiency;
}scheduler_params_t;


//class that envelopes functionality for the scheduler
class scheduler_t{
 public:
  scheduler_t(scheduler_params_t scheduler_params);
  void schedule_stage(void); //wrapper function for the system simulator


 private:
  prep_t* _prep; //pointer to the prep stage
  system_sim_t* _system; //pointer to the system
  system_storage_t* _storage; //pointer to the storage unit
  struct trace_t** system_queue; //reference to the head node of the system_queue
  unsigned long batch_size;//max number of components for a transaction
  unsigned long sequencing_depth; //depth of sequencing required
  float bytes_per_strand; //density of the encoding
  float efficiency; //denotes how many garbage strands will occur for a given file
  typedef void(scheduler_t::*reorder_policy)(void);//the idea of the reorder function is to reorder I/O operations in the system_queue
  typedef void(scheduler_t::*schedule_policy)(void); //the scheduler policy dictates how to group transactions together
  typedef void(scheduler_t::*calcstrands)(unsigned long&,unsigned long&, unsigned long,
					  unsigned long,float, float, unsigned long); //calculate desired/undesired strands for a certain file request

  reorder_policy reorder;
  schedule_policy scheduler;
  calcstrands strand_calculator;
  
  void reorder_none(void); //no reordering policy
  void scheduler_anypool(void); //allow any pool to be batched together
  void calc_singlefile(unsigned long& desired_strands_sequenced,
		       unsigned long& undesired_strands_sequenced,
		       unsigned long file_size,
		       unsigned long pool_ID,
		       float efficiency,
		       float bytes_per_strand,
		       unsigned long sequencing_depth); //calculate the strands for a certain file size

  
};

#endif
