#ifndef SCHEDULER
#define SCHEDULER

class system_storage_t;
class prep_t;
class system_sim_t;
struct transaction_t;



//class that envelopes functionality for the scheduler
class scheduler_t{
  prep_t* _prep;
  system_sim_t* _system;
  system_storage_t* _storage; 
  struct transaction_t* system_queue;
  unsigned long number_components;//max number of components for a transaction
  scheduler_t(system_storage_t* _storage, transaction_t* system_queue, prep_t*  _prep,
	      system_sim_t* _system, unsigned long number_components); 
  typedef void(*scheduler_t::reorder_policy)(void);//the idea of the reorder function is to reorder I/O operations in the system_queue
  typedef void(*scheduler_t::schedule_policy)(void); //the scheduler policy dictates how to send instructions into the pipeline and how to group transactions together
  typedef void(*scheduler_t::calcstrands)(unsigned long&,unsigned long&, unsigned long,
					  unsigned long,float, float)// determine the desired strands and undesired strands fields of each request such that we know how many strands are used for the file that we actually want
  reorder_policy reorder;//will be used as the pointer to a scheduling policy
  schedule_policy scheduler; //scheduler function pointer
  calcstrands strand_calculator;
  void schedule_stage(void); //wrapper function for the system simulator
}

#endif
