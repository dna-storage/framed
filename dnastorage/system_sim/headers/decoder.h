#ifndef DEC
#define DEC

class system_unit_t;
class system_sim_t;

typedef list_entry_t;


class decoder_unit_t: public system_unit{
 public:
  decoder_unit(unsigned long timer, unsigned long num_channels);
 
};


//overall class that holds the prep unit array and supporting member functions to operate on the decoder units
class decoder_t{
 public:
  decoder_unit_t** decoder_set;
  unsigned long num_decoders;
  system_sim_t* _system;
  list_entry_t* seq_dec_buffer;
  unsigned long buffer_size;
  decoder_t(unsigned long timer, unsigned long num_channels, unsigned long buffer_size, unsigned long num_decoders,list_entry_t* seq_dec_buffer,system_sim_t* _system );
  ~decoder_t();

  typedef void (decoder_t::*policy)(void);
  policy decoder_policy;
};


#endif
