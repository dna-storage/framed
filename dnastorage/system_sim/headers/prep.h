#ifndef PREP
#define PREP

class system_unit_t;
class system_sim_t;
class system_storage_t;


//typedef list_entry_t;

class prep_unit_t: public system_unit_t{
 public:
  prep_unit_t(unsigned long num_channels);
 
};


//overall class that holds the prep unit array and supporting member functions to operate on the prep units
class prep_t{
 public:
  prep_unit_t** prep_set;
  unsigned long num_preps;
  system_sim_t* _system;
  list_entry_t* prep_seq_buffer;
  unsigned long buffer_size;
  system_storage_t* dna_storage;
  unsigned long base_timer;
  
  
  prep_t(unsigned long timer, unsigned long num_channels, unsigned long buffer_size,
	 unsigned long num_preps,list_entry_t* prep_seq_buffer,
	 system_sim_t* _system,system_storage_t* dna_storage);
  ~prep_t();

  
  //functions that are called from the top level simulator
  void prep_backend(void); //function that will remove jobs from prep stations and place them into the prep_seq_buffer
  void prep_frontend(void); //function that will try to service jobs that are waiting 

};







#endif
