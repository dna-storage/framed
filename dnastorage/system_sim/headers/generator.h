#ifndef GENERATOR
#define GENERATOR
#include <random>

class system_sim_t;
struct transaction_t;

typedef struct{
  float rate;
  unsigned long max_file_size;
  unsigned long min_file_size;
  unsigned long unique_pools;
  int random_seed;
  system_sim_t* _system;
} generator_params_t;




class generator_t{//class that implements the generator, places transactions into the syste queue
  float rate; //rate at which transactions will be generated
  int random_seed; //seed for the request generator
  unsigned long max_file_size; //max size a file can be
  unsigned long min_file_size; //min size a file can be
  unsigned long unique_pools; //number of unique pools in the system
  struct transaction_t* system_queue; //system queue to add transactions to
  system_sim_t* _system;

  std::default_random_engine* rand_pois;
  std::poisson_distribution<int>* poisson_transactions;
  std::default_random_engine* rand_file;
  std::default_random_engine* rand_pool;
  FILE* trace_file; //trace file pointer in case the generator is reading from a trace
  generator_t(generator_params_t generator_params);
  ~generator_t();
  typedef void(*generator_t::generator_source)(void);
  generator_source gen;
  void generator_stage(void);
  
}




#endif
