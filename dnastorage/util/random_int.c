#include <Python.h>
#include <time.h>

//Python module seed C random number generator
static PyObject* generate_seed(PyObject *self){
  srand(time(NULL));
  Py_RETURN_NONE;
}


//Python module to help with random integers, generates a random integer in a specified range
static PyObject* generate_rand_in_range(PyObject *self, PyObject *args){
  int lower,upper;
  if(!PyArg_ParseTuple(args,"ii",&lower,&upper)){
    return NULL;
  }
  
  return Py_BuildValue("i",(rand()%(upper+1))+lower);
}

//Python module function to generate a random number between [0,1] (double)
static PyObject* generate_rand(PyObject *self){

  return Py_BuildValue("d",rand()/(double)RAND_MAX);
  
}


static PyMethodDef generate_methods[] = {
   { "rand_in_range", (PyCFunction)generate_rand_in_range, METH_VARARGS, NULL },
   {"seed",(PyCFunction)generate_seed,METH_NOARGS,NULL},
   {"rand",(PyCFunction)generate_rand,METH_NOARGS,NULL},
   { NULL, NULL, 0, NULL}
};


static struct PyModuleDef generate={
				    PyModuleDef_HEAD_INIT,
				    "generate",
				    NULL,
				    -1,
				    generate_methods
};

PyMODINIT_FUNC PyInit_generate(void){
  return PyModule_Create(&generate);
}
