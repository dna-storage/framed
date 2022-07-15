
  template<typename DNAConstraint, typename Reward,template <typename> class Context>
  bool search_tree<DNAConstraint,Reward,Context>::isIncomplete()
{
  return c.index < h->max_index;
}


template<typename DNAConstraint, typename Reward, template <typename> class Context>
  search_tree<DNAConstraint,Reward,Context>::search_tree(hedge *_h,
				  search_tree *parent,
				  const Context<DNAConstraint> &_c,
				  float _score,
				  uint32_t _offset,
				  std::string *_observed,
				  int nbits,
				  uint32_t bit,
				  char base,
				  bool insertion,
				  char _kind )
  :h(_h),c(_c),score(_score),offset(_offset), observed(_observed), kind(_kind), guess(base)
{
  if (!insertion) {    
    char tmp = c.nextSymbolWithUpdate(nbits, bit, base);
    assert (tmp == base);
    if (c.index == h->get_index_range(h->adj_seq_bits))
      {
	c.salt = c.prev;
	c.prev = 0;
      }
    if(nbits>0){
      std::shared_ptr<bit_tree<uint8_t>> p = parent->bits;
      for(int i=0; i<nbits; i++){
	bits=std::make_shared<bit_tree<uint8_t>>(p,(bit>>i)&0x01);
	p=bits;
      }
    }
    else{
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

  
template<typename DNAConstraint, typename Reward, template <typename> class Context>
  search_tree<DNAConstraint,Reward, Context>::search_tree(hedge *_h,
				  const Context<DNAConstraint> &_c,
				  float _score, 
				  uint32_t _offset,
				  std::string *_observed)
  :h(_h),c(_c),score(_score),offset(_offset), observed(_observed), bits(nullptr), bases(nullptr)
{
}



  template<typename DNAConstraint, typename Reward,  template <typename> class Context>
search_tree<DNAConstraint,Reward,Context>
search_tree<DNAConstraint,Reward,Context>::addMatch(char base, uint8_t nbits, uint32_t val)
{
  Reward r(h->rate);
  //std::cout << "Got match?" << std::endl;
  
  return search_tree<DNAConstraint,Reward,Context>(h,    // use same decoder settings
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


template<typename DNAConstraint, typename Reward, template <typename> class Context>
search_tree<DNAConstraint,Reward,Context>
search_tree<DNAConstraint,Reward,Context>::addSubst(char base, uint8_t nbits, uint32_t val)
{
  Reward r(h->rate);

  //std::cout << "Got sub?" << std::endl;
    
  return search_tree<DNAConstraint,Reward,Context>(h,    // use same decoder settings
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

template<typename DNAConstraint, typename Reward,  template <typename> class Context>
search_tree<DNAConstraint,Reward,Context>
search_tree<DNAConstraint,Reward,Context>::addDel(char base, uint8_t nbits, uint32_t val)
{
  
  Reward r(h->rate);

  //std::cout << "Got del?" << std::endl;
  
  return search_tree<DNAConstraint,Reward,Context>(h,    // use same decoder settings
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

template<typename DNAConstraint, typename Reward,  template <typename> class Context>
search_tree<DNAConstraint,Reward,Context>
search_tree<DNAConstraint,Reward,Context>::addIns()
{
  Reward r(h->rate);
  //std::cout << "Got ins?" << std::endl;
  
  return search_tree<DNAConstraint,Reward,Context>(h,    // use same decoder settings
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

template<typename DNAConstraint, typename Reward, template <typename> class Context>
search_tree<DNAConstraint,Reward,Context>
search_tree<DNAConstraint,Reward,Context>::addIns2(char base, uint8_t nbits, uint32_t val, double penalty)
{
  Reward r(h->rate);
  //std::cout << "Got ins2?" << std::endl;

  return search_tree<DNAConstraint,Reward,Context>(h,    // use same decoder settings
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


template<typename DNAConstraint, typename Reward,  template <typename> class Context>
std::vector<search_tree<DNAConstraint,Reward,Context>> search_tree<DNAConstraint,Reward,Context>::make0bitGuesses()
{
  char c = this->c.getNextSymbol(0, 0);

  std::vector<search_tree<DNAConstraint,Reward,Context>> ret;
  
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

//TODO:KV: In the future we may need to create a different helper that is more based towards codewords, e.g. checking codeword values for a sliding verification
template<typename DNAConstraint, typename Reward,  template <typename> class Context>
void search_tree<DNAConstraint,Reward,Context>::guessHelper(std::vector<search_tree<DNAConstraint,Reward,Context>> &ret,
					       char c,
					       int nbits,
					       uint32_t val)
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


template<typename DNAConstraint, typename Reward,  template <typename> class Context>
std::vector<search_tree<DNAConstraint,Reward,Context>> search_tree<DNAConstraint,Reward,Context>::make1bitGuesses()
{
  char c0 = c.getNextSymbol(1, 0);
  char c1 = c.getNextSymbol(1, 1);

  std::vector<search_tree<DNAConstraint,Reward,Context>> ret;

  bool inserted = false;

  guessHelper(ret, c0, 1, 0);
  guessHelper(ret, c1, 1, 1);

  return ret;  
}

  template<typename DNAConstraint, typename Reward, template <typename> class Context>
std::vector<search_tree<DNAConstraint,Reward,Context>> search_tree<DNAConstraint,Reward,Context>::makeCWGuesses(){
  std::vector<search_tree<DNAConstraint,Reward,Context>> ret;
  char next_guess;
  uint32_t val = NULL_VALUE;
  //this bits value identifies how much the current codeword counts towards the number of stored bits 
  uint32_t bits_next_state = c.at_codeword_end() ? h->get_n_bits(c.index+1) : h->get_n_bits(c.index); //we need to check if we are transitioning from a complete codeword
  uint32_t bits_this_state = h->get_n_bits(c.index);
  while((next_guess = c.getNextSymbol(bits_next_state,val))!=0){
    uint32_t bits = bits_this_state;
    //do something with this guess
    if(val == NULL_VALUE) bits = 0; //nullify the number of bits, we haven't reached a codeword yet
    guessHelper(ret,next_guess,bits,val);
  }
  return ret;
}


  
  template<typename DNAConstraint, typename Reward, template <typename> class Context>
std::vector<search_tree<DNAConstraint,Reward,Context>> search_tree<DNAConstraint,Reward,Context>::makeGuesses()
{
  uint8_t nbits = h->get_n_bits(c.index);
  hedges::hedge_rate hedge_type = h->get_rate();
  //std::cout << "index=" << c.index << " bits=" << int(nbits) << std::endl;

  if(hedge_type==hedges::hedge_rate::codewords){
    //base guesses for codewords
    return makeCWGuesses();
  }
  else if ( nbits == 0 ) {
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

template<typename DNAConstraint, typename Reward,template <typename> class Context>
bool search_tree<DNAConstraint,Reward,Context>::checkPad()
{
  return bits.get()->check_for_val(h->adj_pad_bits,0);
}




template<template <typename> class Context = context>
std::ostream& operator << (std::ostream& os, Context<Constraint> &c)
{
  os << std::hex;
  os << c.prev << " " << c.salt << " " << c.index ;
  return os;
}

template<typename Constraint = Constraint, typename Reward = Reward, template <typename> class Context=context>
std::ostream & operator << (std::ostream &o, search_tree<Constraint, Reward, Context> &node)
{
  o << "i:" << std::dec << node.c.index << " " << "score:" << node.score;

  if (node.offset < node.observed->size())
    o << " obs[" << node.offset << "]=" << (*node.observed)[node.offset] ;
  else
    o << " [X]";

  o  << " " << node.kind << ":" << node.guess ;

  o << " context: " << node.c << " ";
  
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


template<typename DNAConstraint, typename Reward, template <typename> class Context>
std::vector<search_tree<DNAConstraint,Reward,Context>> search_tree<DNAConstraint,Reward,Context>::make2bitGuesses()
{
  char c0 = c.getNextSymbol(2, 0);
  char c1 = c.getNextSymbol(2, 1);
  char c2 = c.getNextSymbol(2, 2);
  char c3 = c.getNextSymbol(2, 3);

  std::vector<search_tree<DNAConstraint,Reward,Context>> ret; 

  guessHelper(ret, c0, 2, 0);
  guessHelper(ret, c1, 2, 1);
  guessHelper(ret, c2, 2, 2);
  guessHelper(ret, c3, 2, 3);

  // guess unknown insertion
  //ret.push_back( addIns() );      

  return ret;  
}

template <typename Constraint, typename Reward, template <typename> class Context>
uint32_t hedge::decode(std::string &observed,
		   std::vector<uint8_t> &seqId,
		   std::vector<uint8_t> &message,
		   int max_guesses)
{
  using node = search_tree<Constraint, Reward, Context>;
  
  std::priority_queue<node,std::vector<node>,BestScore<node> > heap;

  while( observed.size() < max_index )
    observed.append("A");
  
  Context<Constraint> state(prev_bits,salt_bits);
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
  while(b != res2.end() && adj_seq_bits>0){
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



