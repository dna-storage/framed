#include <Python.h>
#include <time.h>

//Python module seed C random number generator
static PyObject* generate_seed(PyObject *self){
  srand(time(NULL));
  Py_RETURN_NONE;
}


//Python module to help with random integers
static PyObject* generate_rand_in_range(PyObject *self, PyObject *args){
  int lower,upper;
  if(!PyArg_ParseTuple(args,"ii",&lower,&upper)){
    return NULL;
  }
  
  return Py_BuildValue("i",(rand()%(upper+1))+lower);
}


static PyMethodDef generate_methods[] = {
   { "rand_in_range", (PyCFunction)generate_rand_in_range, METH_VARARGS, NULL },
   {"seed",(PyCFunction)generate_seed,METH_NOARGS,NULL},
   { NULL, NULL, 0, NULL }
};

PyMODINIT_FUNC initgenerate() {
   Py_InitModule3("generate", generate_methods, "Faster random number generation");
}
