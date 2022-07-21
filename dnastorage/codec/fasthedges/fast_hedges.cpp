#include <cmath>
#include <cstdlib>
#include <cstdint>
#include <cassert>
#include <vector>
#include <string>
#include <iostream>
#include <queue>
  
#include "fast_hedges.hpp"

using namespace hedges;

using hr = hedge_rate;

namespace hedges {

uint64_t ranhash(uint64_t u)
{
  /* Logic adapted from Press et al. */
  uint64_t v = u * 3935559000370003845ul + 2691343689449507681ul;
  v ^= v >> 21;
  v ^= v << 37;
  v ^= v >> 4;
  v *= 4768777513237032717;
  v ^= v << 20;
  v ^= v >> 41;
  v ^= v << 5;
  return  v;
}

uint64_t digest(uint64_t prev, uint64_t prev_bits,
		uint64_t index, uint64_t index_bits,
		uint64_t salt, uint64_t salt_bits,
		uint64_t mod)
{
  uint64_t prev_mask = (1ULL << prev_bits) - 1;
  uint64_t index_mask = (1ULL << index_bits) - 1;
  uint64_t salt_mask = (1ULL << salt_bits) - 1;

  //std::cout << "digest: " << std::hex << prev_mask << " " << index_mask << " " << salt_mask << " " << salt_bits << " " << (1ULL << salt_bits) << std::endl;
  
  uint64_t t =  ((((index&index_mask) << prev_bits) | (prev & prev_mask)) << salt_bits) | (salt&salt_mask);
  t = ranhash(t) % mod;
  return t;
}

  
  hedge::hedge(double rate,
	       int seq_bytes,
	       int message_bytes,
	       int pad_bits,
	       int prev_bits,
	       int salt_bits,
	       int codeword_sync_period,
	       int parity_period
	       )
{
  this->raw_rate = rate;
  this->codeword_sync_period=codeword_sync_period;
  this->rate = pick_rate(rate);
  this->seq_bytes = seq_bytes;
  this->message_bytes = message_bytes;
  this->pad_bits = pad_bits;
  this->prev_bits = prev_bits;
  this->salt_bits = salt_bits;
  this->parity_period = parity_period;

  //Take into account parity when determining adjusted message bits so that the index can be correctly calcualted
  uint32_t total_parameter_data_bits = message_bytes*8+pad_bits;
  if(this->parity_period>0){
    assert(this->parity_period>1 && "Parity can't be set to a period of just 1");
    //parity is interleaved into pad_bits speced by user, so that needs to be taken into account in the total bits
    total_parameter_data_bits = (this->parity_period+1) + (total_parameter_data_bits-this->parity_period)/(this->parity_period-1) + total_parameter_data_bits;   
  }
  
  this->adj_seq_bits = seq_bytes*8 + check_if_padding_needed(seq_bytes*8);
  this->adj_pad_bits = this->pad_bits + check_if_padding_needed(total_parameter_data_bits);
  this->adj_message_bits = total_parameter_data_bits + (this->adj_pad_bits-this->pad_bits);
  max_index = get_index_range(adj_message_bits + adj_seq_bits);
  
}



void Constraint::next(char c) {
  if ( last[index] == c )
    run++;
  else
    run = 1;

  if ( c == 'A' || c == 'T')
    AT++;
  else
    GC++;

  index = (index+1)%12;  
  if (last[index] == 'A' || last[index] == 'T')
    AT--;
  // we need the check due to initialization
  else if (last[index] == 'C' || last[index] == 'G')
    GC--;
  
  last[index] = c;
  
  assert (run < 4);
}

std::vector<char> Constraint::get(int mod) const
{
  auto choose = std::vector<char>({'A','C','G','T'});    
  if (mod>2) {
    return choose;
  } else {
    if ( run > 2 ) {
      switch(last[index]) {
      case 'A': choose = std::vector<char>({'C','G','T'}); break;
      case 'C': choose = std::vector<char>({'A','G','T'}); break;
      case 'G': choose = std::vector<char>({'A','C','T'}); break;
      case 'T': choose = std::vector<char>({'A','C','G'}); break;
      default: assert(false);
      }
    }      
    if ( choose.size() == 4 && std::abs(AT-GC)>3 ) {
      if ( AT > GC )
	{
	  //std::cout << "here-GC!" << std::endl;
	  choose = std::vector<char>({'C','G'});
	}
      else 
	{
	  //std::cout << "here-AT!" << std::endl;		  
	  choose = std::vector<char>({'A','T'});
	}
    }      
    return choose;
  }      
}


void hedge::print(bool extra)
{
  std::cout << "-------- hedges --------" << std::endl;
  std::cout << "rate:               " << this->raw_rate << std::endl;
  std::cout << "rate id:            " << int(this->rate) << std::endl;
  std::cout << "seq_bytes:          " << this->seq_bytes << std::endl;
  std::cout << "message_bytes:      " << this->message_bytes << std::endl;
  std::cout << "salt bits:          " << this->salt_bits << std::endl;
  std::cout << "prev bits:          " << this->prev_bits << std::endl;
  std::cout << "index range:        " << 0 << " to " << max_index << std::endl;

  if (extra)
    {
      std::cout << "Extra details: " << std::endl;
      std::cout << "adj message bits:    " << this->adj_message_bits << std::endl;
      std::cout << "adj seq bits:        " << this->adj_seq_bits << std::endl;
      std::cout << "adj pad bits:        " << this->adj_pad_bits << std::endl;
    }
}



std::string hedge::encode(std::vector<uint8_t> seqId, std::vector<uint8_t> message)
{
  bitwrapper seq(seqId);    
  bitwrapper mess(message);
  std::string buff;
  
  std::vector<int> pat = patterns[ (int) rate ];
  int len = pattern_length[ (int) rate ];
  
  uint32_t bit = 0;
  uint32_t index = 0;
  
  context<Constraint> state(prev_bits, salt_bits);
  // initial static salt; not really required
  // could be good to base this on the primer
  // or file id
  state.salt = 0xA5A5A5A5;

  std::cout << std::hex;
  
  // Encode seqId
  while(index < get_index_range(adj_seq_bits)) {
    int nbits = pat[index%len];
    int val = seq.get_bits(bit, bit+nbits);

    //std::cout << std::dec << index << ": " << std::hex << state << std::endl;
    
    char c = state.nextSymbolWithUpdate(nbits,val,0);
    buff.push_back(c);
    
    index ++;
    bit += nbits;
  }

  //std::cout << "after seq: " << std::dec << index << ": " << std::hex << state << std::endl;  
    
  // Copy prev over to be the salt
  state.salt = state.prev;
  // ? Should we set prev to 0 or not? Since we have the salt now, we don't need the prev.
  //   Probably doesn't make much of a difference since it's just redundant information.
  state.prev = 0;
  bit = 0;

  // Encode message
  while(index < max_index) {
    int nbits = pat[index%len];
    int val = mess.get_bits(bit, bit+nbits);

    //std::cout << std::dec << index << ": " << std::hex << state << std::endl; 
    char c = state.nextSymbolWithUpdate(nbits,val,0);
    buff.push_back(c);
    
    index++;
    bit += nbits;
  }

  std::cout << std::dec;
  
  return buff;
}



} // end namespace
  
int main()
{
  int message_size = 15;
  
  hedge h(1/4.0,4,message_size,4,8,8,0,0);

  std::vector<uint8_t> arr = {100,102,103,106};
  std::string buff;

  std::vector<uint8_t> data;
  for(int i=0; i<message_size; i++) {
    uint8_t b = (uint8_t)(ranhash(i)&0xFF);
    //data[i] = b;
    data.push_back(b);
    std::cout << int(data[i]) << ", ";
  }
  std::cout << std::endl;
    
  for (int i=0; i<1; i++)
    buff = h.encode(arr,data);

  std::cout << "correct: " << buff << std::endl;

  std::vector<uint8_t> mess(10), seq(4);

  std::string observed = buff;

  std::cout << "length of buff = " << buff.size() << std::endl;

  observed[1] = 'A';
  observed[30] = 'C';
  observed[35] = 'T';

  bool t;
  for (int i=0; i<1; i++)
    {
      t = h.decode(observed, seq, mess, 100000);
      if ((i+1)%100==0)
	std::cout << "." ;
    }

  if (t) std::cout << "correct: " << buff << std::endl;  
  if (t) std::cout << "decoding succeeded." << std::endl;

  for(auto s: seq)
    std::cout << int(s) << ", ";
  std::cout << std::endl;

  for(auto s: mess)
    std::cout << int(s) << ", ";
  std::cout << std::endl;
  
  return 0;
}
