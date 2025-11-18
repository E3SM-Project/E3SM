#ifndef EAMXX_FIELD_DLPACK_HPP
#define EAMXX_FIELD_DLPACK_HPP

#include <Python.h>
#include <dlpack/dlpack.h>

#include "share/field/field.hpp"
// #include "share/core/eamxx_pysession.hpp"

#include <ekat_assert.hpp>

namespace scream {

// Define the C++ structure that represents our Python object
typedef struct {
  PyObject_HEAD
  // DLManagedTensorVersioned dlmtv; // Embedded DLManagedTensor
  DLManagedTensor dlmtv; // Embedded DLManagedTensor
} PyDLTensor;

// The __dlpack__ implementation
static PyObject* PyDLTensor_dlpack(PyObject* self, PyObject* noargs) {
  PyDLTensor* wrapper = (PyDLTensor*)self;

  // Check if already consumed by checking the deleter function
  if (wrapper->dlmtv.deleter == nullptr && wrapper->dlmtv.manager_ctx != nullptr) {
    PyErr_SetString(PyExc_RuntimeError, "DLManagedTensor has already been consumed.");
    return NULL;
  }
  // Create the PyCapsule pointing to the DLManagedTensor
  // The capsule takes ownership of calling the deleter.
  PyObject* capsule = PyCapsule_New(&wrapper->dlmtv, "dltensor", NULL);
  
  // Crucially, we MUST clear the context/deleter pointers AFTER creating the capsule
  // to mark it as consumed and prevent double-deletion (though here it's "double no-op")
  wrapper->dlmtv.deleter = nullptr; 
  wrapper->dlmtv.manager_ctx = nullptr; // Also clear context

  return capsule;
}

// Method definitions for our custom Python type
static PyMethodDef PyDLTensor_methods[] = {
  {"__dlpack__", PyDLTensor_dlpack, METH_NOARGS, 
   "Returns a PyCapsule containing the DLManagedTensor."},
  {NULL, NULL, 0, NULL}
};

// Deallocator for the Python wrapper object itself
static void PyDLTensor_dealloc(PyDLTensor* self) {
  // No action needed here regarding the external buffer.
  // We only deallocate the PyObject wrapper structure itself.
  Py_TYPE(self)->tp_free((PyObject*)self);
}

static PyTypeObject PyDLTensorType = {
  PyVarObject_HEAD_INIT(NULL, 0) /* Head initialization macro */
  "DLTensor",           /* tp_name: For example, 'my_module.DLTensor' */
  sizeof(PyDLTensor),   /* tp_basicsize: Size of the structure */
  0,                         /* tp_itemsize: 0 for fixed-size objects */
  
  /* Method pointers: */
  (destructor)PyDLTensor_dealloc, /* tp_dealloc: Object destructor */
  0,                         /* tp_vectorcall_offset: Offset of tp_vectorcall, or 0 */
  0,                         /* tp_getattr: Old-style attribute lookup, usually 0 */
  0,                         /* tp_setattr: Old-style attribute setting, usually 0 */
  0,                         /* tp_as_async: Pointer to async/await functions (PyAsyncMethods*) */
  0,                         /* tp_repr: Object representation (repr()) */
  
  /* Number protocols */
  0,                         /* tp_as_number: Pointer to number methods (PyNumberMethods*) */
  0,                         /* tp_as_sequence: Pointer to sequence methods (PySequenceMethods*) */
  0,                         /* tp_as_mapping: Pointer to mapping methods (PyMappingMethods*) */
  
  /* Hash, Call, Str */
  0,                         /* tp_hash: Object hash (hash()) */
  0,                         /* tp_call: Call object as function (obj()) */
  0,                         /* tp_str: String representation (str()) */
  
  /* Get/Set attribute (modern style) */
  PyObject_GenericGetAttr,   /* tp_getattro: Generic get attribute */
  PyObject_GenericSetAttr,   /* tp_setattro: Generic set attribute */
  
  0,                         /* tp_as_buffer: Pointer to buffer methods (PyBufferProcs*) */
  
  /* Flags, Docstring, Order */
  Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_GC, /* tp_flags: Standard flags */
  "DLTensor object to enable zero-copy data transfer.", /* tp_doc: Docstring */
  
  /* Garbage collection (GC) */
  0,                         /* tp_traverse: GC traversal function */
  0,                         /* tp_clear: GC clear function */
  
  /* Comparisons */
  0,                         /* tp_richcompare: Rich comparison function */
  
  /* Weak references */
  0,                         /* tp_weaklistoffset: Offset of weak reference list, or 0 */
  
  /* Iteration */
  0,                         /* tp_iter: Get iterator (iter()) */
  0,                         /* tp_iternext: Next item in iterator (next()) */
  
  /* Descriptors/Methods */
  PyDLTensor_methods,   /* tp_methods: Object methods */
  0,                         /* tp_members: Object members (PyMemberDef*) */
  0,                         /* tp_getset: Object getters/setters (PyGetSetDef*) */
  0,                         /* tp_base: Base class (e.g., &PyBaseObject_Type) */
  0,                         /* tp_dict: Class dictionary */
  0,                         /* tp_descr_get: Descriptor get function */
  0,                         /* tp_descr_set: Descriptor set function */
  0,                         /* tp_dictoffset: Offset to instance dictionary */
  
  /* Initialization and Allocation */
  0,                         /* tp_init: Object initializer (__init__) */
  0,                         /* tp_alloc: Object allocator */
  PyType_GenericNew,         /* tp_new: Object constructor (__new__) */
   /* Type-specific memory management (Python 3.8+ specific fields might be here depending on Python version) */
  0,                         /* tp_free */
  0,                         /* tp_is_gc */
  0,                         /* tp_bases */
  0,                         /* tp_mro */
  0,                         /* tp_cache */
  0,                         /* tp_subclasses */
  0,                         /* tp_weaklist */
  0,                         /* tp_del */
  0,                         /* tp_version_tag */
  0,                         /* tp_finalize */
  // Note: Python versions > 3.4 added 'tp_finalize'. The list above matches ~Python 3.8 structure.
};

// ------------------- Create PyDLTensor from a Field ------------------ //

inline std::vector<std::int64_t> get_strides (const FieldHeader& fh)
{
  const auto& fid = fh.get_identifier();
  const auto& fl  = fid.get_layout();
  const int rank = fl.rank();

  std::vector<std::int64_t> strides(rank);
  if (rank==0) {
    return strides;
  }

  const auto& fap = fh.get_alloc_properties ();
  auto p = fh.get_parent();
  if (p) {
    strides = get_strides(*p);
    const auto& sv_info = fap.get_subview_info();
    strides.erase(strides.begin()+sv_info.dim_idx);
  } else {
    auto dims = fl.dims();
    dims.back() = fap.get_last_extent();
    strides.back() = 1;
    for (int i=rank-1; i>=1; --i) {
      strides[i-1] = strides[i]*dims[i];
    }
  }

  return strides;
}

// Get offset of field actual data from the start of the internal data pointer
inline std::uint64_t get_offset (const FieldHeader& fh)
{
  auto p = fh.get_parent();
  if (not p) {
    return 0;
  }

  std::uint64_t p_offset = get_offset(*p);

  const auto& fap = fh.get_alloc_properties ();
  const auto& sv_info = fap.get_subview_info();
  if (sv_info.dim_idx>0) {
    return p_offset;
  }

  return p_offset*get_strides(fh)[0]*sv_info.slice_idx;
}

PyObject* create_dl_tensor (Field& f)
{
  if (PyType_Ready(&PyDLTensorType) < 0) {
    // Handle error (e.g., return false, throw exception)
    return nullptr;
  }
  // Create an instance of the Python wrapper type
  // PyDLTensor* wrapper = (PyDLTensor*)PyDLTensorType.tp_new(&PyDLTensorType, NULL, NULL);
  // auto wrapper_obj = PyDLTensorType.tp_new(&PyDLTensorType, NULL, NULL);
  PyDLTensor* wrapper = (PyDLTensor*)PyDLTensorType.tp_new(&PyDLTensorType, NULL, NULL);
  if (!wrapper) return NULL;

  DLTensor& dl_tensor = wrapper->dlmtv.dl_tensor;
  dl_tensor = f.get_header().get_dltensor();
  dl_tensor.data = f.get_internal_view_data();
  // wrapper->dlmtv.version.major = DLPACK_MAJOR_VERSION;
  // wrapper->dlmtv.version.minor = DLPACK_MINOR_VERSION;
  // if (f.is_read_only()) {
  //   wrapper->dlmtv.flags |= DLPACK_FLAG_BITMASK_READ_ONLY;
  // } else {
  //   wrapper->dlmtv.flags &= ~DLPACK_FLAG_BITMASK_READ_ONLY;
  // }

  printf("for f=%s, device=%d, device_id=%d\n",f.name().c_str(),dl_tensor.device.device_type,dl_tensor.device.device_id);

  return (PyObject*) wrapper;
}

} // namespace scream

#endif // EAMXX_FIELD_DLPACK_HPP
