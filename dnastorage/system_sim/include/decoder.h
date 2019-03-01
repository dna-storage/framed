#ifndef DEC
#define DEC

//forward declarations
class system_unit_t;
class system_sim_t;
class buffer_t;
class stats_t;

typedef struct{
  unsigned long timer;
  unsigned long num_channels;
  unsigned long buffer_size;
  unsigned long num_decoders;
  buffer_t* seq_dec_buffer;
  system_sim_t* _system;
  stats_t* stats;
} decoder_params_t; //bundled decoder parameters


class decoder_unit_t:public system_unit_t{
 public:
  unsigned long component_ID; //identifies what component of a transaction a decoder is working on
  decoder_unit_t(unsigned long num_channels);
 
};

//overall class that holds the decoder unit array and supporting member functions to operate on the decoder units
class decoder_t{
 public:
  decoder_t(decoder_params_t decoder_params);
  ~decoder_t();
  void decoder_stage(void); //top level function called by the system simulator
 private:
  decoder_unit_t** decoder_set; //array of decoder_unit_t pointers
  unsigned long num_decoders; //number of decoders in the system
  system_sim_t* _system; //pointer to the system simulator
  buffer_t* seq_dec_buffer; //buffer between the sequencer and the decoders
  unsigned long base_timer; //initial value for decoder timers, used to refresh exhausted timer counters
  
  void decoder_backend(void); //remove finished transactions from decoders
  void decoder_frontend(void); //move transactions from the seq_dec_buffer to an open decoder
  void decoder_complete(unsigned long decoder_ID); //service a completed decoder
  void decoder_timestep(unsigned long decoder_ID); //decrement the timer counter for a decoder
  unsigned long decoder_avail(void); //check to see if there is a decoder free
  void init_decoder(unsigned long decoder_ID, unsigned long transaction_ID, unsigned long seq_dec_idnex);//setup a newly selected decoder
  stats_t* stats;
};


#endif
