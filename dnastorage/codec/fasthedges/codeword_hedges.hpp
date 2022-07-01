"""
Filename: codeword_hedges.hpp
Description: Support for hedges to allow for codeword based implementations. Support for codeword+assembly will also be included at some point.
Author: Kevin Volkel
"""

#include <iostream>
#include <vector>
#include <cstdint>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include <string>
#include "fast_hedges.hpp"


namespace codeword_hedges{
  
uint8_t convert_base_to_bits(const char& base){ //binary representation for bases
  switch(base){
  case "A":
    return 0x00;
  case "G":
    return 0x01;
  case "C":
    return 0x02;
  case "T":
    return 0x03;
  }
  assert(0); //shouldn't reach here
}

char convert_bits_to_base(const uint8_t& base){ //binary representation for bases
  switch(base){
  case 0x00:
    return "A";
  case 0x01:
    return "G";
  case 0x02:
    return "C";
  case 0x03:
    return "T";
  }
  assert(0); //shouldn't reach here
}

//convert DNA string to byte array that will be compatible with bitwrapper object
std::vector<uin8_t> convert_DNA_to_bytes(const std::string& DNA){
  uint32_t pad_bits = 8 - (DNA.size()*2)%8;
  if (pad_bits==8){
    pad_bits=0; //pad_bits can only be 8 if % is 0
  }
  std::vector<uint8_t> DNA_bytes((DNA.size()*2+pad_bits)/8); //make a vector of bytes large enough to hold binary representation of DNA
  uint8_t shift=0;
  uint32_t byte_index=0;
  for(int i=0; i<DNA.size(); i++){
    assert(shift<8); //make sure shift doesn't fly off the rails
    uint8_t base_bits=convert_base_to_bits(DNA[i]);
    base_bits = (base_bits&0x01)<<1 | (base_bits>>1); //flipping the endianess here so when we use bitwrapper it comes out the right way
    DNA_bytes[byte_index]|=(base_bits<<shift);
    shift+=2;
    if(shift%8==0){
      shift=0;
      byte_index++;
    }
  }
  return DNA_bytes;
}

#define MAX_BASES 4
#define NULL_VALUE 0xffffffff
  
  //DNAtrie holds a tree of DNA strings as a compressed bit representation, necessary for tracking guesses along the tree for decoding 
  class DNAtrie{
  public:
    DNAtrie(std::vector<uint8_t> dna,uint8_t dna_length):_dna(dna),_prefix_length(dna_length){ //initialize node
      for(int i=0; i<MAX_BASES; i++) children[i] = nullptr;
    }
    DNAtrie():_prefix_length(0){ //initialize empty node
      for(int i=0; i<MAX_BASES; i++) children[i] = nullptr;
    }
    
    bool insert(const std::string& DNA,uint32_t value){
      assert(this->is_root()); //searches should start from the root
      std::vector<uint8_t> DNA_as_bytes = convert_DNA_to_bytes(DNA);
      hedges::bitwrapper DNA_bits(DNA_as_bytes);
      bool found=true;
      uint32_t input_length = DNA.size();
      uint32_t total_indexes_covered=0;
      DNAtrie* current_node=this;
      while(total_indexes_covered<input_length){
	current_node = current_node->children[DNA_bits.get_bits(total_indexes_covered*2,total_indexes_covered*2+2)].get();
	//check if node needs to be split
	hedges::bitwrapper node_bits(_dna);
	for(int i=0; i<_prefix_length; i++){
	  uint32_t node_base = node_bits.get_bits(i*2,i*2+2);
	  if(total_indexes_covered>=input_length) break;
	  uint32_t DNA_base = DNA_bits.get_bits(i*2,i*2+2);
	  if((DNA_base^node_base)!=0){

	    uint32_t bases_remaining_DNA = input_length-total_indexes_covered;
	    uint32_t bases_remaining_node = current_node->_prefix_length-i;

	    uint32_t node_pad_bits = (8 - (bases_remaining_node*2)%8)==8 ? 0 : (8 - (bases_remaining_node*2)%8);
	    uint32_t DNA_pad_bits  = (8 - (bases_remaining_DNA*2)%8)==8 ? 0 :  (8 - (bases_remaining_DNA*2)%8);

	    std::vector<uint8_t> new_node_strand((bases_remaining_node*2+node_pad_bits)/8);
	    std::vector<uint8_t> new_DNA_strand((bases_remaining_DNA*2+DNA_pad_bits)/8);

	    hedges::bitwrapper new_node_bits(new_node_strand);
	    hedges::bitwrapper new_DNA_bits(new_DNA_strand);
	    //we need to split this node at position i
	    for(int j=i; j<current_node->_prefix_length; j++){ //first make new DNA string for the new node
	      uint32_t new_node_base = node_bits.get_bits(j*2,j*2+2);
	      new_node_bits.set_bit((j-i)*2,(new_node_base>>1&0x01)); //most sig bit should be first bitin new node
	      new_node_bits.set_bit((j-i)*2+1,new_node_base&0x01);
	    }
	    //now make the new string for the input string
	    for(int j=total_indexes_covered; j<DNA.size(); j++){
	      uint32_t new_node_base = DNA_bits.get_bits(j*2,j*2+2);
	      new_node_bits.set_bit((j-total_indexes_covered)*2,(new_node_base>>1&0x01));
	      new_node_bits.set_bit((j-total_indexes_covered)*2+1,new_node_base&0x01);
	    }
	    std::shared_ptr<DNAtrie> DNA_node = std::make_shared<DNAtrie>(new_DNA_strand,bases_remaining_DNA);
	    std::shared_ptr<DNAtrie> new_node = std::make_shared<DNAtrie>(new_node_strand,bases_remaining_node);
	    DNA_node->_value = value;

	    //now connect to the new-internal node everything that was on this current node
	    for(int j=0; j<MAX_BASES; j++){
	      new_node->set_child(current_node->children[j],j);
	      current_node->children[j]=nullptr; //null out old connections
	    }
	    current_node->children[node_base] = new_node; //finally insert the new node as the split
	    current_node->children[DNA_base] = DNA_node; //insert the node that arrisses from the input string
	    if(current_node->_value != NULL_VALUE){
	      new_node->_value = current_node->_value; //transfer value ownership down the tree
	      current_node -> _value = NULL_VALUE; //this current node shouldn't have ownership of a value anymore 
	    }
	    found=false; //the addition of a node means that the string was indeed not found
	    current_node = current_node->children[DNA_base].get(); //go to the next node which was made for the input strand so algorithm can finish out
	    current_node->_prefix_length= current_node->prefix_length-bases_remaining_node;
	    break;
	  }
	  total_indexes_covered++; //tracking where at in the input strand we have traversed
	}
      }
      return found;
    }
   

    bool is_leaf(){
      return _value!=NULL_VALUE; //leaf nodes should only exist where non-NULL values exist
    }

    bool is_root(){
      return _prefix_length==0; //root should be only node empty
    }

    void set_child(std::shared_ptr<DNAtrie> child_node,uint8_t child_base){children[child_base]=child_node;}

    char get_index(uint8_t index){ //gets the base at the index for this node
      if(index>=_prefix_length) return 0x00; //return nothing if there is no more string to this node
      hedges::bitwrapper dna_bits(_dna);
      return convert_bits_to_base(dna_bits.get_bits(index*2,index*2+2)); //return the base
    }


    std::vector<char> get_transition_bases(){//get the transition bases for this node
      std::vector<char> return_bases;
      for(int i=0; i<MAX_BASES; i++){
	if(children[i]!=nullptr) return_bases.push_back(convert_bits_to_base(i&0x0f));
      }
      return return_bases;
    }

    uint32_t get_value(){ return _value;} //return the value of 
    
  protected:
    std::vector<uint8_t> _dna; //local dna string for this node
    std::shared_ptr<DNAtrie> children[MAX_BASES];//pointer to children
    uint8_t _prefix_length; //length of the string at this node
    uint32_t _value=NULL_VALUE; //value stored for the DNA string
  };


  class basic_codeword_hedges_context : public hedges::context{
  public:
    //this class is going to provide the main mechanism for traversing the tree while we make guesses in a strand
    codeword_hedges_context(int _prev_bits, int _salt_bits)
      :hedges::context(_prev_bits,_salt_bits)
    {
      //TODO: Need to grab pointers from the hedges module to the necessary DNATrie object
    }
    
  private:
    DNAtrie* _codeword_trie; //pointer to the codeword tree we will use to make guesses

    
  };

    
  
} //end namespace codeword_hedges
