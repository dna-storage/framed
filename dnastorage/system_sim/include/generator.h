#ifndef GENERATOR
#define GENERATOR
#include <random>

class system_sim_t;
class stats_t;

struct transaction_t;

typedef struct{
  float rate;
  unsigned long max_file_size;
  unsigned long min_file_size;
  unsigned long unique_pools;
  int random_seed;
  system_sim_t* _system;
  stats_t* stats;
} generator_params_t; //bundle generator parameters




class generator_t{//class that implements the generator, places transactions into the syste queue
 public:

  generator_t(generator_params_t generator_params);
  ~generator_t();
  void generator_stage(void); //interface with the top level simulator
  void generator_stop(void); //signal to the generator to stop
 private:
  std::default_random_engine* rand_tool; //used to help the poisson distrubution function
  std::poisson_distribution<unsigned long>* poisson_transactions; //used to generate poisson random numbers
  std::uniform_int_distribution<unsigned long>* rand_pool;//uniform distributions for pools
  std::uniform_int_distribution<unsigned long>* rand_file;//uniform distribtion for file sizes
  FILE* trace_file; //trace file pointer in case the generator is reading from a trace
  float rate; //rate at which transactions will be generated
  int random_seed; //seed for the request generator
  unsigned long max_file_size; //max size a file can be
  unsigned long min_file_size; //min size a file can be
  unsigned long unique_pools; //number of unique pools in the system
  int stop; //signals whether the generator should keep pushing out requests
  struct transaction_t* system_queue; //system queue to add transactions to
  system_sim_t* _system;
  stats_t* stats; //pointer to stat handling object
  typedef void(generator_t::*generator_source)(void);
  generator_source gen;
  void gen_poisson(void); //generator that uses the poisson distribution
  

  
};




#endif