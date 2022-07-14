#ifndef SHARED_HPP
#define SHARED_HPP
#include <Python.h>
#include "fast_hedges.hpp"
#include "codeword_hedges.hpp"

long getLong(PyObject *Object, const char *attr)
{
  if (PyObject_HasAttrString(Object, attr))
    {
      PyObject *a = PyObject_GetAttrString(Object,attr);
      if (PyLong_Check(a))
	{
	  return PyLong_AsLong(a);
	}
    }
  return -1;
}

double getDouble(PyObject *Object, const char *attr)
{
  if (PyObject_HasAttrString(Object, attr))
    {
      PyObject *a = PyObject_GetAttrString(Object,attr);
      if (PyFloat_Check(a))
	{
	  return PyFloat_AsDouble(a);
	}
    }
  return -1.0;
}

static hedges::hedge make_hedge_from_pyobject(PyObject *object)
{
  double rate = getDouble(object,"rate");
  int seq_bytes = getLong(object, "seq_bytes");
  int message_bytes = getLong(object, "message_bytes");
  int pad_bits = getLong(object, "pad_bits");
  int prev_bits = getLong(object, "prev_bits");
  int salt_bits = getLong(object, "salt_bits");

  if (pad_bits == -1) {
    if (rate > 0.33)
      pad_bits = 8;
    else if (rate > 0.125)
      pad_bits = 4;
    else
      pad_bits = 1;
  }
  if (prev_bits == -1)
    prev_bits = 8;
  if (salt_bits == -1)
    salt_bits = 8;
  
  hedges::hedge h(rate, seq_bytes, message_bytes, pad_bits, prev_bits, salt_bits);
  return h;
}



template<typename Constraint = hedges::Constraint, typename Reward = hedges::Reward,  template <typename> class Context = hedges::context>
static PyObject *
shared_decode(PyObject *self, PyObject *args)
{
    const char *strand;
    PyObject *hObj;
    int guesses = 100000;
    
    if (PyArg_ParseTuple(args, "sO|i", &strand, &hObj, &guesses)) {

      hedges::hedge h = make_hedge_from_pyobject(hObj);
      
      std::vector<uint8_t> mess(h.message_bytes), seq(h.seq_bytes);

      std::string sstrand(strand);
      uint32_t t = h. template decode<Constraint,Reward,Context>(sstrand,seq,mess,guesses);

      // std::cout << "[";
      // for(auto i : seq)
      // 	std::cout << int(i) << ", ";
      // for(auto i : mess)
      // 	std::cout << int(i) << ", ";
      // std::cout << "]" << std::endl; 

      int sz = seq.size() + mess.size();
      PyObject *list = PyList_New(sz);
      for(auto i=0; i<sz; i++)
	{	  
	  if (i<seq.size()) {
	    if (i >= t) {
	      PyObject *item = Py_BuildValue("s",NULL);
	      Py_INCREF(item);
	      PyList_SetItem(list,i,item);

	    } else {
	      PyObject *item = Py_BuildValue("i",seq[i]);
	      Py_INCREF(item);
	      PyList_SetItem(list,i,item);
	    }
	  } else {
	    if (i >= t) {
	      PyObject *item = Py_BuildValue("s",NULL);
	      Py_INCREF(item);
	      PyList_SetItem(list,i,item);
	    } else {
	      PyObject *item = Py_BuildValue("i",mess[i-seq.size()]);
	      Py_INCREF(item);
	      PyList_SetItem(list,i,item);
	    }
	  }
	  
	}
      return list;
    }

    return Py_BuildValue("s",NULL);
}




static PyObject *
create_codebook(PyObject*codebook, const char* name, PyObject* exception)
{
  codeword_hedges::DNAtrie* trie_root = new codeword_hedges::DNAtrie();
  PyObject* PyDNA;
  PyObject* PyKey;
  Py_ssize_t pos = 0;
  while(PyDict_Next(codebook,&pos,&PyKey,&PyDNA)){
    if(!PyUnicode_Check(PyDNA)){
      PyErr_SetString(exception,"Dictionary values are not unicode");
      return NULL;
    }
    const char* strand = PyUnicode_AsUTF8(PyDNA);
    std::string DNA(strand);
    uint32_t value = (uint32_t)PyLong_AsLong(PyKey);
    trie_root->insert(DNA,value);
  }
  //add the trie to the module
  std::string codebook_name(name);

  if(codeword_hedges::CodebookMap.find(codebook_name)!=codeword_hedges::CodebookMap.end()){
    PyErr_WarnEx(PyExc_RuntimeWarning,"Codebook already exists, deleting old codebook to avoid memory leak",1);
    delete codeword_hedges::CodebookMap[codebook_name];
    codeword_hedges::CodebookMap.erase(codebook_name);
  }
  
  codeword_hedges::CodebookMap[codebook_name] = trie_root;

#ifdef DEBUG
  trie_root->print();
#endif
  
  return Py_BuildValue("s",NULL);
}

static PyObject* destroy_codebook(const char* name, PyObject* exception){ //removes a codebook associated with the module
  std::string codebook_name(name);
  if(codeword_hedges::CodebookMap.find(codebook_name) == codeword_hedges::CodebookMap.end()){
    PyErr_WarnEx(PyExc_RuntimeWarning,"Codebook does not exist on delete, returning to avoid double deletion",1);
    return NULL;
  }
  if(codeword_hedges::CodebookMap[codebook_name]==nullptr){
    PyErr_SetString(exception,"Nullptr being freed on codebook destroy operation");
    return NULL;
  }
  delete codeword_hedges::CodebookMap[codebook_name];
  codeword_hedges::CodebookMap.erase(codebook_name);			       
  return Py_BuildValue("s",NULL);
}


#endif
