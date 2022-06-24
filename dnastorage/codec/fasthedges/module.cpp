#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <limits>
#include "fast_hedges.hpp"

using namespace hedges;

static PyObject *FasthedgesError;

// static PyObject *
// fasthedges_system(PyObject *self, PyObject *args)
// {
//     const char *command;
//     int sts;

//     if (!PyArg_ParseTuple(args, "s", &command))
//         return NULL;
//     sts = system(command);
//     if (sts < 0) {
//         PyErr_SetString(FasthedgesError, "System command failed");
//         return NULL;
//     }
//     return PyLong_FromLong(sts);
// }

// static PyObject *
// fasthedges_parser(PyObject *self, PyObject *args)
// {
//     const char *command;
//     PyObject * num;
//     int sts=0;

//     if (PyArg_ParseTuple(args, "O", &num)) {
      
//       printf("Got string = %s.\n",Py_TYPE(num)->tp_name);
//       //return PyLong_FromLong(0);
//     }

//     //if (PyArg_ParseTuple(args, "s", &command)) {
      
//     //printf("Got string = %s.\n",command);
//       //return PyLong_FromLong(0);
//     //}



//     return Py_BuildValue("s",NULL);
// }

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

static hedge make_hedge_from_pyobject(PyObject *object)
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
  
  hedge h(rate, seq_bytes, message_bytes, pad_bits, prev_bits, salt_bits);
  return h;
}

static PyObject *
fasthedges_encode2(PyObject *self, PyObject *args)
{
    const char *command;
    Py_buffer *buff = (Py_buffer *) malloc(sizeof(*buff));
    double rate;
    unsigned long num_index_bytes;
    unsigned long num_message_bytes;
    PyObject *hObj;

    if (PyArg_ParseTuple(args, "y*O", buff, &hObj)) {

      hedge h = make_hedge_from_pyobject(hObj);
      
      std::vector<uint8_t> index;
      std::vector<uint8_t> message;
      
      uint8_t *buf = (uint8_t*) buff->buf;
      
      for(int i=0; i<buff->len && i<h.seq_bytes; i++)
	{
	  //printf("%d,",buf[i]);
	  index.push_back(buf[i]);
	}
      for(int i=h.seq_bytes; i<buff->len && i<h.message_bytes+h.seq_bytes; i++)
	{
	  //printf("%d,",buf[i]);
	  message.push_back(buf[i]);
	}

      std::string buff = h.encode(index, message);

      PyObject *o = Py_BuildValue("s",buff.c_str());
      return o;      
    }

    return Py_BuildValue("s",NULL);
}

static PyObject *
fasthedges_decode2(PyObject *self, PyObject *args)
{
    const char *strand;
    PyObject *hObj;
    int guesses = 100000;
    
    if (PyArg_ParseTuple(args, "sO|i", &strand, &hObj, &guesses)) {

      hedge h = make_hedge_from_pyobject(hObj);
      
      std::vector<uint8_t> mess(h.message_bytes), seq(h.seq_bytes);

      std::string sstrand(strand);
      uint32_t t = h.decode(sstrand,seq,mess,guesses);

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

static PyObject * decode_item(const char *strand, hedge &h, int guesses=1000000)
{
  std::string sstrand(strand);

  std::vector<uint8_t> mess(h.message_bytes), seq(h.seq_bytes);
  bool t = h.decode(sstrand,seq,mess,guesses);

  // std::cout << "[";
  // for(auto i : seq)
  //   std::cout << int(i) << ", ";
  // for(auto i : mess)
  //   std::cout << int(i) << ", ";
  // std::cout << "]" << std::endl;

  uint32_t size = h.seq_bytes + h.message_bytes;
  PyObject *list = PyList_New(size);
  for(auto i=0; i<size; i++)
    {
      if (i<h.seq_bytes) {
	PyObject *item = Py_BuildValue("i",seq[i]);
	Py_INCREF(item);
	PyList_SetItem(list,i,item);
      } else {
	PyObject *item = Py_BuildValue("i",mess[i-seq.size()]);
	Py_INCREF(item);
	PyList_SetItem(list,i,item);
      }
      
    }
  return list;
}

 
static PyObject *
fasthedges_bulk_decode(PyObject *self, PyObject *args)
{
    const char *strand;
    PyObject *l;
    PyObject *hObj;
    int guesses = 1000000;
    
    if (PyArg_ParseTuple(args, "OO|i", &l, &hObj, &guesses)) {

      hedge h = make_hedge_from_pyobject(hObj);
      
      if (l!=NULL && !PyList_Check(l)) {
         PyErr_SetString(FasthedgesError, "Expected first argument to be a list of DNA strands.");
         return NULL;	
      }

      PyObject *bulk_list = PyList_New(PyList_Size(l));

      for(Py_ssize_t j=0; j<PyList_Size(l); j++)
	{
	  PyObject *item = PyList_GetItem(l, j);
	  
	  if (PyUnicode_Check(item))
	    {
	      const char *strand = PyUnicode_AsUTF8(item);
	      PyObject *decoded_list = decode_item(strand, h);
	      Py_INCREF(decoded_list);
	      PyList_SetItem(bulk_list,j,decoded_list);	      
	    }
	}

      return bulk_list;
    }

    //if (PyArg_ParseTuple(args, "s", &command)) {
      
    //printf("Got string = %s.\n",command);
      //return PyLong_FromLong(0);
    //}

    return Py_BuildValue("s",NULL);
}





static PyObject *
fasthedges_echo(PyObject *self, PyObject *args)
{
    PyObject *hObj;
    
    if (PyArg_ParseTuple(args, "O", &hObj))
      {
	hedge h = make_hedge_from_pyobject(hObj);

	h.print(true);
      }
    return Py_BuildValue("s",NULL);
}

static PyMethodDef FasthedgesMethods[] = {
    {"encode",  fasthedges_encode2, METH_VARARGS, "Encode into DNA."},
    {"decode",  fasthedges_decode2, METH_VARARGS, "Decode from DNA back into bytes."},
    {"bulk_decode",  fasthedges_bulk_decode, METH_VARARGS, "Bulk decode a list of strands from DNA back into bytes."},
    {"echo",  fasthedges_echo, METH_VARARGS, "Echo hedges configuration."},
    
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

static struct PyModuleDef fasthedgesmodule = {
    PyModuleDef_HEAD_INIT,
    "fasthedges",   /* name of module */
    NULL, /* module documentation, may be NULL */
    -1,       /* size of per-interpreter state of the module,
                 or -1 if the module keeps state in global variables. */
    FasthedgesMethods
};

PyMODINIT_FUNC PyInit_fasthedges(void)
{
    PyObject *m;

    m = PyModule_Create(&fasthedgesmodule);
    if (m == NULL)
        return NULL;

    FasthedgesError = PyErr_NewException("fasthedges.error", NULL, NULL);
    Py_XINCREF(FasthedgesError);
    if (PyModule_AddObject(m, "error", FasthedgesError) < 0) {
        Py_XDECREF(FasthedgesError);
        Py_CLEAR(FasthedgesError);
        Py_DECREF(m);
        return NULL;
    }

    return m;
}
