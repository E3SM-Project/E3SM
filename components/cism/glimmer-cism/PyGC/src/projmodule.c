/* 

   Copyright (C) 2004, 2009, 2010
   Glimmer-CISM contributors - see AUTHORS file for list of contributors

   This file is part of Glimmer-CISM.

   Glimmer-CISM is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 2 of the License, or (at
   your option) any later version.

   Glimmer-CISM is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with Glimmer-CISM.  If not, see <http://www.gnu.org/licenses/>.

   Glimmer-CISM is hosted on BerliOS.de:
   https://developer.berlios.de/projects/glimmer-cism/
*/

#include <Python.h>
#include <numpy/arrayobject.h>
#include <projects.h>
#include <stdio.h>

typedef struct {
    PyObject_HEAD
    PJ *projection;
}Proj;

static void Proj_dealloc(Proj* self)
{
  pj_free(self->projection);
  self->ob_type->tp_free((PyObject*)self);
}

static PyObject * Proj_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  Proj *self;
  
  self = (Proj *)type->tp_alloc(type, 0);
  if (self != NULL) 
    {
      self->projection = NULL;
    }
  
  return (PyObject *)self;
}

static int Proj_init(Proj *self, PyObject *args, PyObject *kwds)
{
  PyListObject  *proj_params;
  PyObject *o;
  static char *kwlist[] = {"params"};

  char **params;
  int num_params;
  int i;

  /* parsing arguments */
  if (!PyArg_ParseTuple(args, "O!",&PyList_Type, &proj_params))
    return -1;
  
  /* setting up projection */
  num_params = PyList_Size((PyObject *) proj_params);
  params = (char **) malloc(num_params*sizeof(char *));
  for (i=0;i<num_params;i++)
    {
      o = PyList_GetItem((PyObject *) proj_params, i);
      if (PyString_Check(o)) 
	*(params+i) = PyString_AsString(o);
      else
	{
	  PyErr_SetString(PyExc_ValueError,"Projection parameter list contains non-string types");
	  free(params);
	  return -1;
	}
    }

  /* setup projection */
  if (! (self->projection=pj_init(num_params,params)))
    {
      PyErr_SetString(PyExc_RuntimeError,pj_strerrno(pj_errno));
      free(params);
      return -1;
    }

  return 0;
}

static PyObject *Proj_fwd(Proj* self, PyObject * args)
{
  PyListObject  *input;
  PyObject *o;
  int n;
  projUV data;
  
  /* parsing arguments */
  if (!PyArg_ParseTuple(args, "O!", &PyList_Type, &input))
    return NULL;
  n = PyList_Size((PyObject *) input);
  if (n!=2)
    {
      PyErr_SetString(PyExc_ValueError,"Input must be list of two values");
      return NULL;
    }
  if ((o = PyNumber_Float(PyList_GetItem((PyObject *) input, 0)))==NULL)
    {
      PyErr_SetString(PyExc_ValueError,"Must be a number.");
      return NULL;
    }
  data.u = PyFloat_AsDouble((PyObject *) o) * DEG_TO_RAD;

  if ((o = PyNumber_Float(PyList_GetItem((PyObject *) input, 1)))==NULL)
    {
      PyErr_SetString(PyExc_ValueError,"Must be a number.");
      return NULL;
    }
  data.v = PyFloat_AsDouble((PyObject *) o) * DEG_TO_RAD;
  
  data = pj_fwd(data,self->projection);
  if (data.u==HUGE_VAL && data.v==HUGE_VAL)
    {
      PyErr_SetString(PyExc_RuntimeError,pj_strerrno(pj_errno));
      return NULL;
    }
  
  return Py_BuildValue("[dd]",data.u,data.v);
}

static PyObject *Proj_inv(Proj* self, PyObject * args)
{
  PyListObject  *input;
  PyObject *o;
  int n;
  projUV data;
  
  /* parsing arguments */
  if (!PyArg_ParseTuple(args, "O!", &PyList_Type, &input))
    return NULL;
  n = PyList_Size((PyObject *) input);
  if (n!=2)
    {
      PyErr_SetString(PyExc_ValueError,"Input must be list of two values");
      return NULL;
    }
  if ((o = PyNumber_Float(PyList_GetItem((PyObject *) input, 0)))==NULL)
    {
      PyErr_SetString(PyExc_ValueError,"Must be a number.");
      return NULL;
    }
  data.u = PyFloat_AsDouble((PyObject *) o);

  if ((o = PyNumber_Float(PyList_GetItem((PyObject *) input, 1)))==NULL)
    {
      PyErr_SetString(PyExc_ValueError,"Must be a number.");
      return NULL;
    }
  data.v = PyFloat_AsDouble((PyObject *) o);
  
  data = pj_inv(data,self->projection);

  if (data.u==HUGE_VAL && data.v==HUGE_VAL)
    {
      PyErr_SetString(PyExc_RuntimeError,pj_strerrno(pj_errno));
      return NULL;
    }
  
  return Py_BuildValue("[dd]",data.u*RAD_TO_DEG,data.v*RAD_TO_DEG);
}

static PyObject *Proj_gridfwd(Proj* self, PyObject * args)
{
  PyTupleObject  *input;
  PyArrayObject *PyX, *PyY;  
  PyObject *o;  
  int n,nx,ny,i;
  projUV data;

  /* parsing arguments */
  if (!PyArg_ParseTuple(args, "O!", &PyTuple_Type, &input))
    return NULL;
  /* get tupe size */
  n = PyTuple_Size((PyObject *) input);
  if (n!=2)
    {
      PyErr_SetString(PyExc_ValueError,"Input must be tuple with two values");
      return NULL;
    }
  /* convert tuple to two arrays */
  o = PyTuple_GetItem((PyObject *) input, 0);
  if ((PyX = (PyArrayObject *) PyArray_CopyFromObject(o, NPY_DOUBLE, 1, 0)) == NULL)
    return NULL;
  nx = PyArray_Size((PyObject *) PyX);
  o = PyTuple_GetItem((PyObject *) input, 1);
  if ((PyY = (PyArrayObject *) PyArray_CopyFromObject(o, NPY_DOUBLE, 1, 0)) == NULL)
    return NULL;
  ny = PyArray_Size((PyObject *) PyY);
  /* checking they are the same size */
  if (nx!=ny)
    {
      PyErr_SetString(PyExc_ValueError,"Input arrays must have the same size");
      return NULL;
    }
  
  for (i=0;i<nx;i++)
    {
      data.u = (*(double *) (PyX->data+i*PyX->strides[0])) * DEG_TO_RAD;
      data.v = (*(double *) (PyY->data+i*PyY->strides[0])) * DEG_TO_RAD;
      
      data = pj_fwd(data,self->projection);
      
      (*(double *) (PyX->data+i*PyX->strides[0])) = data.u;
      (*(double *) (PyY->data+i*PyY->strides[0])) = data.v;
    }
  
  return Py_BuildValue("(OO)",PyArray_Return(PyX),PyArray_Return(PyY));
}

static PyObject *Proj_gridinv(Proj* self, PyObject * args)
{
  PyTupleObject  *input;
  PyArrayObject *PyX, *PyY;  
  PyObject *o;  
  int n,nx,ny,i;
  projUV data;

  /* parsing arguments */
  if (!PyArg_ParseTuple(args, "O!", &PyTuple_Type, &input))
    return NULL;
  /* get tupe size */
  n = PyTuple_Size((PyObject *) input);
  if (n!=2)
    {
      PyErr_SetString(PyExc_ValueError,"Input must be tuple with two values");
      return NULL;
    }
  /* convert tuple to two arrays */
  o = PyTuple_GetItem((PyObject *) input, 0);
  if ((PyX = (PyArrayObject *) PyArray_CopyFromObject(o, NPY_DOUBLE, 1, 0)) == NULL)
    return NULL;
  nx = PyArray_Size((PyObject *) PyX);
  o = PyTuple_GetItem((PyObject *) input, 1);
  if ((PyY = (PyArrayObject *) PyArray_CopyFromObject(o, NPY_DOUBLE, 1, 0)) == NULL)
    return NULL;
  ny = PyArray_Size((PyObject *) PyY);
  /* checking they are the same size */
  if (nx!=ny)
    {
      PyErr_SetString(PyExc_ValueError,"Input arrays must have the same size");
      return NULL;
    }
  
  for (i=0;i<nx;i++)
    {
      data.u = (*(double *) (PyX->data+i*PyX->strides[0]));
      data.v = (*(double *) (PyY->data+i*PyY->strides[0]));
      
      data = pj_inv(data,self->projection);
      
      (*(double *) (PyX->data+i*PyX->strides[0])) = data.u * RAD_TO_DEG;
      (*(double *) (PyY->data+i*PyY->strides[0])) = data.v * RAD_TO_DEG;
    }
  
  return Py_BuildValue("(OO)",PyArray_Return(PyX),PyArray_Return(PyY));
}

static int  Proj_print(Proj *self, FILE *file, int flags)
{
  char *info;
  
  if ((info=pj_get_def(self->projection,0)) == NULL) 
    {
      PyErr_SetString(PyExc_RuntimeError,pj_strerrno(pj_errno));
      return -1;
    }
  fprintf(file,"%s",info);
  return 0;
}


static PyMethodDef Proj_methods[] = {
    {"fwd", (PyCFunction)Proj_fwd, METH_VARARGS, "Do the forward projection."},
    {"inv", (PyCFunction)Proj_inv, METH_VARARGS, "Do the inverse projection."},
    {"gridfwd", (PyCFunction)Proj_gridfwd, METH_VARARGS, "Do the forward projection of an entire array."},
    {"gridinv", (PyCFunction)Proj_gridinv, METH_VARARGS, "Do the forward projection of an entire array."},
    {NULL}  /* Sentinel */
};

static PyTypeObject ProjType = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "_proj4.Proj",             /*tp_name*/
    sizeof(Proj),              /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)Proj_dealloc,  /*tp_dealloc*/
    (printfunc)Proj_print,     /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    0,                         /*tp_repr*/
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash */
    0,                         /*tp_call*/
    0,                         /*tp_str*/
    0,                         /*tp_getattro*/
    0,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,        /*tp_flags*/
    "python type wrapping libproj",             /* tp_doc */
    0,		               /* tp_traverse */
    0,		               /* tp_clear */
    0,		               /* tp_richcompare */
    0,		               /* tp_weaklistoffset */
    0,		               /* tp_iter */
    0,		               /* tp_iternext */
    Proj_methods,              /* tp_methods */
    0,                         /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)Proj_init,       /* tp_init */
    0,                         /* tp_alloc */
    Proj_new,                  /* tp_new */
};

static PyMethodDef proj_methods[] = {
    {NULL}  /* Sentinel */
};

#ifndef PyMODINIT_FUNC	/* declarations for DLL import/export */
#define PyMODINIT_FUNC void
#endif
PyMODINIT_FUNC
init_proj4(void) 
{
  PyObject* m;
  
  ProjType.tp_new = PyType_GenericNew;
  if (PyType_Ready(&ProjType) < 0)
    return;
  
  m = Py_InitModule3("_proj4", proj_methods,"Proj4 Python module");
  
  Py_INCREF(&ProjType);
  PyModule_AddObject(m, "Proj", (PyObject *)&ProjType);
  import_array();
}

