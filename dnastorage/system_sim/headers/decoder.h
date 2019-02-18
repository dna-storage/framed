#ifndef DEC
#define DEC

class system_unit_t;
class system_sim_t;

struct list_entry_t;

typedef struct{
  unsigned long timer;
  unsigned long num_channels;
  unsigned long buffer_size;
  unsigned long num_decoders;
  struct list_entry_t* seq_dec_buffer;
  system_sim_t* _system;
} decoder_params_t;


class decoder_unit_t:public system_unit_t{
 public:
  unsigned long component_ID; //identifies what component of a transaction a decoder is working on
  decoder_unit_t(unsigned long num_channels);
 
};

//overall class that holds the prep unit array and supporting member functions to operate on the decoder units
class decoder_t{
 public:
  decoder_unit_t** decoder_set;
  unsigned long num_decoders;
  system_sim_t* _system;
  struct list_entry_t* seq_dec_buffer;
  unsigned long buffer_size;
  unsigned long base_timer;
  decoder_t(decoder_params_t decoder_params);
  ~decoder_t();

  typedef void (decoder_t::*policy)(void);
  policy decoder_policy; //function pointer ensures that we dont change the interface for other routines upon the desire of a different policy

  //functions to be called from the top level simulator
  unsigned long decoder_backend(void); //decoder backend will remove jobs from each decoder,and return the number of jobs completed
  void decoder_frontend(void);//search the seq_dec_buffer to find jobs that can be placed into decoders
  //function that can be used to deactivate a decoder
  void decoder_complete(unsigned long decoder_ID);
  void init_decoder(unsigned long decoder_ID,unsigned long transaction_ID);
  unsigned long decoder_avail(void);
};


#endif
