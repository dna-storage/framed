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
  uint64_t prev_mask = (1 << prev_bits) - 1;
  uint64_t index_mask = (1 << index_bits) - 1;
  uint64_t salt_mask = (1 << salt_bits) - 1;
  uint64_t t =  ((((index&index_mask) << prev_bits) | (prev & prev_mask)) << salt_bits) | (salt&salt_mask);
  t = ranhash(t) % mod;
  return t;
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

template<typename DNAConstraint, typename Reward>
search_tree<DNAConstraint,Reward>::search_tree(hedge *_h,
				  search_tree *parent,
				  const context<DNAConstraint> &_c,
				  float _score,
				  uint32_t _offset,
				  std::string *_observed,
				  int nbits,
				  int bit,
				  char base,
				  bool insertion,
				  char _kind )
  :h(_h),c(_c),score(_score),offset(_offset), observed(_observed), kind(_kind), guess(base)
{
  if (!insertion) {    
    char tmp = c.nextSymbolWithUpdate(nbits, bit);
    assert (tmp == base);
    if (c.index == h->get_index_range(h->adj_seq_bits))
      {
	c.salt = c.prev;
	c.prev = 0;
      }
    if (nbits==2) {
      bits = std::make_shared<bit_tree<uint8_t>>(parent->bits,bit&1);
      bits = std::make_shared<bit_tree<uint8_t>>(bits,(bit>>1)&1);
    } else if (nbits==1) {
      bits = std::make_shared<bit_tree<uint8_t>>(parent->bits,bit&1);
    } else {
      bits = parent->bits;
    }
    bases = std::make_shared<bit_tree<char>>(parent->bases,base);
  } else {
    // no update to corrected strand or message bits
    bases = parent->bases;
    bits = parent->bits;
    kind = 'i';
  }  
}

template<typename DNAConstraint, typename Reward>
search_tree<DNAConstraint,Reward>::search_tree(hedge *_h,
				  const context<DNAConstraint> &_c,
				  float _score,
				  uint32_t _offset,
				  std::string *_observed)
  :h(_h),c(_c),score(_score),offset(_offset), observed(_observed), bits(nullptr), bases(nullptr)
{
}

template<typename DNAConstraint, typename Reward>
search_tree<DNAConstraint,Reward>
search_tree<DNAConstraint,Reward>::addMatch(char base, uint8_t nbits, uint8_t val)
{
  Reward r(h->rate);
  //std::cout << "Got match?" << std::endl;
  
  return search_tree<DNAConstraint,Reward>(h,    // use same decoder settings
					   this, // this node is parent node, needed for
					         // bit_trees
					   c,
					   score+r.getReward(), // reward for matching
					   offset+1, // move forward in string
					   observed, // unchanged
					   nbits,
					   val,
					   base,
					   false // not an insertion
					   );
}

template<typename DNAConstraint, typename Reward>
search_tree<DNAConstraint,Reward>
search_tree<DNAConstraint,Reward>::addSubst(char base, uint8_t nbits, uint8_t val)
{
  Reward r(h->rate);

  //std::cout << "Got sub?" << std::endl;
    
  return search_tree<DNAConstraint,Reward>(h,    // use same decoder settings
					   this, // this node is parent node, needed for
					         // bit_trees
					   c,
					   score+r.getSubPenalty(), // reward for matching
					   offset+1, // move forward in string
					   observed, // unchanged
					   nbits,
					   val,
					   base,
					   false, // not an insertion
					   's'
					   );
}

template<typename DNAConstraint, typename Reward>
search_tree<DNAConstraint,Reward>
search_tree<DNAConstraint,Reward>::addDel(char base, uint8_t nbits, uint8_t val)
{
  
  Reward r(h->rate);

  //std::cout << "Got del?" << std::endl;
  
  return search_tree<DNAConstraint,Reward>(h,    // use same decoder settings
					   this, // this node is parent node, needed for
					         // bit_trees
					   c,
					   score+r.getDelPenalty(), // penalty
					   offset,   // deletion, so offset stays fixed
					   observed, // unchanged
					   nbits,
					   val,
					   base,
					   false, // not an insertion
					   'd'
					   );
}

template<typename DNAConstraint, typename Reward>
search_tree<DNAConstraint,Reward>
search_tree<DNAConstraint,Reward>::addIns()
{
  Reward r(h->rate);
  //std::cout << "Got ins?" << std::endl;
  
  return search_tree<DNAConstraint,Reward>(h,    // use same decoder settings
					   this, // this node is parent node, needed for
					         // bit_trees
					   c,
					   score+r.getInsPenalty(), // penalty
					   offset+1, // insertion, move down 1 in array
					   observed, // unchanged
					   0,
					   0,
					   '_',
					   true, // not an insertion
					   'i'
					   );
}

template<typename DNAConstraint, typename Reward>
search_tree<DNAConstraint,Reward>
search_tree<DNAConstraint,Reward>::addIns2(char base, uint8_t nbits, uint8_t val, double penalty)
{
  Reward r(h->rate);
  //std::cout << "Got ins2?" << std::endl;

  return search_tree<DNAConstraint,Reward>(h,    // use same decoder settings
					   this, // this node is parent node, needed for
					         // bit_trees
					   c,
					   
					   // penalty + reward for matching			
					   score+penalty,
					   
					   offset+2,   // deletion, so offset stays fixed
					   observed,   // unchanged
					   nbits,
					   val,
					   base,
					   false, // not an insertion
					   'I'
					   );
}


template<typename DNAConstraint, typename Reward>
std::vector<search_tree<DNAConstraint,Reward>> search_tree<DNAConstraint,Reward>::make0bitGuesses()
{
  char c = this->c.getNextSymbol(0, 0);

  std::vector<search_tree<DNAConstraint,Reward>> ret;
  
  if ( c == observedAt(offset) )
    {
      // guess match
      //std::cout << "guess 0 bit match! " << this->c.index << " " << c << std::endl;
      ret.push_back( addMatch(c,0,0) );
    }
  else
    {
      // guess substitution
      ret.push_back( addSubst(c,0,0) );
    }
  
  // guess deletion
  ret.push_back( addDel(c,0,0) );      
  

  Reward r;
  if ( c == observedAt(offset+1) )
    {
      // guess specific insertion and skip 2 ahead
      ret.push_back( addIns2(c,0,0, r.getReward()+r.getInsPenalty()) );
    }
  else
    {
      // guess unknown insertion
      ret.push_back( addIns2(c,0,0, r.getSubPenalty()+r.getInsPenalty()) );      
    }
  
  return ret;  
}

template<typename DNAConstraint, typename Reward>
void search_tree<DNAConstraint,Reward>::guessHelper(std::vector<search_tree<DNAConstraint,Reward>> &ret,
					       char c,
					       int nbits,
					       int val)
{
  Reward r;
  
  if ( c == observedAt(offset) )
    {
      // guess match
      //std::cout << "guess " << nbits << " bit match! " << this->c.index << " " << c << std::endl;
      
      ret.push_back( addMatch(c,nbits,val) );
    }
  else
    {
      // guess substitution
      ret.push_back( addSubst(c,nbits,val) );
    }

  // guess deletion
  ret.push_back( addDel(c,nbits,val) );

  if ( c == observedAt(offset+1) )
    {
      // guess specific insertion and skip 2 ahead
      ret.push_back( addIns2(c,nbits,val, r.getReward()+r.getInsPenalty()) );
    }
  else
    {
      // guess unknown insertion
      ret.push_back( addIns2(c,nbits,val, r.getSubPenalty()+r.getInsPenalty()) );      
    }
  
}


template<typename DNAConstraint, typename Reward>
std::vector<search_tree<DNAConstraint,Reward>> search_tree<DNAConstraint,Reward>::make1bitGuesses()
{
  char c0 = c.getNextSymbol(1, 0);
  char c1 = c.getNextSymbol(1, 1);

  std::vector<search_tree<DNAConstraint,Reward>> ret;

  bool inserted = false;

  guessHelper(ret, c0, 1, 0);
  guessHelper(ret, c1, 1, 1);

  // For 1 bit case we can be smarter with insertion guesses and save a little time by
  // jumping ahead if we know the next position matches:
  // if ( c0 != observedAt(offset) && c0 == observedAt(offset+1) )
  //   {
  //     // guess specific insertion and skip 2 ahead
  //     ret.push_back( addIns2(c0,1,0) );
  //     inserted = true;
  //   }

  // if ( c1 != observedAt(offset) && c1 == observedAt(offset+1) )
  //   {
  //     // guess specific insertion and skip 2 ahead
  //     ret.push_back( addIns2(c1,1,1) );
  //     inserted = true;
  //   }      
 
  // if (!inserted)
  //   {
  //     // guess unknown insertion
  //     ret.push_back( addIns() );      
  //   }

  return ret;  
}

template<typename DNAConstraint, typename Reward>
std::vector<search_tree<DNAConstraint,Reward>> search_tree<DNAConstraint,Reward>::make2bitGuesses()
{
  char c0 = c.getNextSymbol(2, 0);
  char c1 = c.getNextSymbol(2, 1);
  char c2 = c.getNextSymbol(2, 2);
  char c3 = c.getNextSymbol(2, 3);

  std::vector<search_tree<DNAConstraint,Reward>> ret; 

  guessHelper(ret, c0, 2, 0);
  guessHelper(ret, c1, 2, 1);
  guessHelper(ret, c2, 2, 2);
  guessHelper(ret, c3, 2, 3);

  // guess unknown insertion
  //ret.push_back( addIns() );      

  return ret;  
}

template<typename DNAConstraint, typename Reward>
std::vector<search_tree<DNAConstraint,Reward>> search_tree<DNAConstraint,Reward>::makeGuesses()
{
  uint8_t nbits = h->get_n_bits(c.index);
  //std::cout << "index=" << c.index << " bits=" << int(nbits) << std::endl;
  if ( nbits == 0 ) {
    return make0bitGuesses();    
  } else if ( nbits == 1 ) {
    return make1bitGuesses();
  } else if ( nbits == 2 ) {
    return make2bitGuesses();
  } else {
    assert (false && "invalid case");
  }
  return {};
}

template<typename DNAConstraint, typename Reward>
bool search_tree<DNAConstraint,Reward>::checkPad()
{
  return bits.get()->check_for_val(h->adj_pad_bits,0);
}

hedge::hedge(double rate,
	     int seq_bytes,
	     int message_bytes,
	     int pad_bits,
	     int prev_bits,
	     int salt_bits)
{
  this->raw_rate = rate;
  this->rate = pick_rate(rate);
  this->seq_bytes = seq_bytes;
  this->message_bytes = message_bytes;
  this->pad_bits = pad_bits;
  this->prev_bits = prev_bits;
  this->salt_bits = salt_bits;
  
  this->adj_seq_bits = seq_bytes*8 + check_if_padding_needed(seq_bytes*8);
  this->adj_pad_bits = pad_bits + check_if_padding_needed(message_bytes*8 + pad_bits);   
  this->adj_message_bits = message_bytes*8 + this->adj_pad_bits;
  
  max_index = get_index_range(adj_message_bits + adj_seq_bits);
  
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


std::ostream& operator << (std::ostream& os, context<Constraint> &c)
{
  os << std::hex;
  os << c.prev << " " << c.salt << " " << c.index ;
  return os;
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
    
    char c = state.nextSymbolWithUpdate(nbits,val);
    buff.push_back(c);
    
    index ++;
    bit += nbits;
  }

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
    char c = state.nextSymbolWithUpdate(nbits,val);
    buff.push_back(c);
    
    index++;
    bit += nbits;
  }

  std::cout << std::dec;
  
  return buff;
}

std::ostream & operator << (std::ostream &o, search_tree<Constraint, Reward> &node)
{
  o << "i:" << std::dec << node.c.index << " " << "score:" << node.score;

  if (node.offset < node.observed->size())
    o << " obs[" << node.offset << "]=" << (*node.observed)[node.offset] ;
  else
    o << " [X]";

  o  << " " << node.kind << ":" << node.guess ;

  o << node.c << " ";
  
  o << " bits: ";
  
  std::list<uint8_t> res2 = node.bits.get()->get_n(16);
  for(auto r : res2)
    std::cout << int(r) ;

  o << " bases: ";
  std::list<char> res3 = node.bases.get()->get_all();
  for(auto r : res3)
    std::cout << r ;

  std::cout << std::endl;
  
  return o;
}

template <typename Constraint, typename Reward>
uint32_t hedge::decode(std::string &observed,
		   std::vector<uint8_t> &seqId,
		   std::vector<uint8_t> &message,
		   int max_guesses)
{
  using node = search_tree<Constraint, Reward>;
  
  std::priority_queue<node,std::vector<node>,BestScore<node> > heap;

  while( observed.size() < max_index )
    observed.append("A");
  
  context<Constraint> state(prev_bits,salt_bits);
  state.salt = 0xA5A5A5A5;    
  heap.push(node(this,state,1000,0,&observed));
    
  while( !heap.empty() && max_guesses > 0 )
    {
      node tmp = heap.top();
     
      //std::cout << "score=" << tmp.score << " index=" << tmp.c.index
      //	<< " offset=" << tmp.offset << std::endl;
      
      if ( tmp.isWinner() )
	{
	  break;	   
	}

      // pop so that new guesses don't rise to the top of the queue
      heap.pop();

      //std::cout << tmp;
      
      if ( tmp.isIncomplete() )
	{
	  std::vector<node> guesses = tmp.makeGuesses();
	  //std::cout << "make " << guesses.size() << " new guesses." << std::endl;
	  for(auto n : guesses) {
	    //std::cout << "score=" << n.score << " index=" << n.c.index << std::endl;
	    //std::cout << "    " << n;
	    heap.push(n);
	  }
	}

      // if (true) {
      // 	std::list<uint8_t> res2 = tmp.bits.get()->get_all();
      // 	std::cout << res2.size() << std::endl;
      // 	for(auto r : res2)
      // 	  std::cout << int(r) ;
      // 	std::cout << std::endl;
      // 	std::list<char> res3 = tmp.bases.get()->get_all();
      // 	std::cout << res3.size() << std::endl;
      // 	for(auto r : res3)
      // 	  std::cout << r ;
      // 	std::cout << std::endl;
      // }
      
      //std::cout << heap.size() << " " << max_guesses << std::endl;
      max_guesses--;
    }

  //std::cout << max_guesses << std::endl;

  // may either find a winner or run out of guesses.
  // in either case, collect the result
  node best = heap.top();
  bool result = best.isWinner();

  uint32_t ret;
  
  // collect seqId and message
  //std::list<char> res = best.bases.get()->get_all();
  //std::string corrected(res.begin(), res.end());
  
  std::list<uint8_t> res2 = best.bits.get()->get_all();
  bitwrapper seq(seqId);
  bitwrapper mess(message);

  if ( res2.size() == (message_bytes + seq_bytes)*8 + pad_bits )
    ret = message_bytes + seq_bytes;
  else if (res2.size() > (message_bytes + seq_bytes)*8)
    ret = message_bytes + seq_bytes;
  else
    // Since we only take whole bytes, this result needs to be truncated to
    // the nearest byte. The last partial byte would be invalid.
    ret = res2.size() / 8;
  
  // Extract message and seqId
  int k=0;
  std::list<uint8_t>::iterator b = res2.begin();
  while(b != res2.end()) {
    if (k<seq_bytes*8)
      seq.set_bit(k,*b);
    k++;
    b++;
    if (k==adj_seq_bits)
      break;
  }
  
  k = 0;
  while(b != res2.end()) {    
    if (k<message_bytes*8)
      mess.set_bit(k,*b);
    else
      break;
    k++;
    b++;
  }
  
  return ret;
}

} // end namespace
  
int main()
{
  int message_size = 15;
  
  hedge h(1/4.0,4,message_size,4,8,8);

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
