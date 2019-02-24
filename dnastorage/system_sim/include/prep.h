#ifndef PREP
#define PREP

class system_unit_t;
class system_sim_t;
class system_storage_t;

struct list_entry_t;


typedef struct{
  unsigned long timer;
  unsigned long num_channels;
  unsigned long buffer_size;
  unsigned long num_preps;
  struct list_entry_t* prep_seq_buffer;
  system_sim_t* _system;
  system_storage_t* dna_storage;
} prep_params_t; //bundled prep parameters

class prep_unit_t: public system_unit_t{
 public:
  prep_unit_t(unsigned long num_channels);
 
};


//overall class that holds the prep unit array and supporting member functions to operate on the prep units
class prep_t{
 public:
  prep_t(prep_params_t prep_params);
  ~prep_t();
  void prep_stage(void); //use this function to interface with the top level simulator
  void prep_stationsubmit(unsigned long prep_ID, unsigned long transaction_ID); //submit a transaction to the specified prep station
    int prep_stationavail(void); //return an available prep station, if none available return -1



 private: 
  prep_unit_t** prep_set;
  unsigned long num_preps; //number of prep stations in the system
  system_sim_t* _system; //pointer to the top level simulator
  struct list_entry_t* prep_seq_buffer; //pointer to the prep_seq_buffer 
  unsigned long buffer_size; //size of the prep_seq_buffer
  unsigned long base_timer; //initial value of the timer for each prep station 
  system_storage_t* dna_storage;
  void prep_backend(void); //function that will remove jobs from prep stations and place them into the prep_seq_buffer
  void prep_complete(unsigned long prep_ID); //service a completed prep station
  int get_prepseq(unsigned long transaction_ID); //find an open spot in the prepseq buffer
  void prep_timestep(unsigned long prep_ID); //decrement the timer for the specified prep unit
    
};







#endif
