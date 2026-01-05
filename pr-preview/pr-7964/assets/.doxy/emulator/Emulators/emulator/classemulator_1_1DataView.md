

# Class emulator::DataView

**template &lt;typename T&gt;**



[**ClassList**](annotated.md) **>** [**emulator**](namespaceemulator.md) **>** [**DataView**](classemulator_1_1DataView.md)



_Non-owning view over contiguous data._ [More...](#detailed-description)

* `#include <data_view.hpp>`





































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**DataView**](#function-dataview-15) () <br>_Default constructor: empty view._  |
|   | [**DataView**](#function-dataview-25) (T \* data, std::size\_t size, [**DataLayout**](namespaceemulator.md#enum-datalayout) layout=DataLayout::ROW\_MAJOR) <br>_Construct from raw pointer and size._  |
|   | [**DataView**](#function-dataview-35) (const T \* data, std::size\_t size, [**DataLayout**](namespaceemulator.md#enum-datalayout) layout=DataLayout::ROW\_MAJOR) <br>_Construct from raw pointer and size (const version)_  |
|   | [**DataView**](#function-dataview-45) (std::vector&lt; T &gt; & vec, [**DataLayout**](namespaceemulator.md#enum-datalayout) layout=DataLayout::ROW\_MAJOR) <br>_Construct from std::vector (non-owning)_  |
|   | [**DataView**](#function-dataview-55) (const std::vector&lt; T &gt; & vec, [**DataLayout**](namespaceemulator.md#enum-datalayout) layout=DataLayout::ROW\_MAJOR) <br>_Construct from const std::vector (non-owning)_  |
|  T \* | [**begin**](#function-begin-12) () <br>_Iterator support._  |
|  const T \* | [**begin**](#function-begin-22) () const<br> |
|  T \* | [**data**](#function-data-12) () <br>_Get raw pointer to data._  |
|  const T \* | [**data**](#function-data-22) () const<br> |
|  bool | [**empty**](#function-empty) () const<br>_Check if view is empty._  |
|  T \* | [**end**](#function-end-12) () <br> |
|  const T \* | [**end**](#function-end-22) () const<br> |
|  [**DataLayout**](namespaceemulator.md#enum-datalayout) | [**layout**](#function-layout) () const<br>_Get memory layout._  |
|  T & | [**operator[]**](#function-operator) (std::size\_t idx) <br>_Element access (bounds checking in debug builds)_  |
|  const T & | [**operator[]**](#function-operator_1) (std::size\_t idx) const<br> |
|  std::size\_t | [**size**](#function-size) () const<br>_Get number of elements._  |




























## Detailed Description


## Public Functions Documentation




### function DataView [1/5]

_Default constructor: empty view._ 
```C++
inline emulator::DataView::DataView () 
```




<hr>



### function DataView [2/5]

_Construct from raw pointer and size._ 
```C++
inline emulator::DataView::DataView (
    T * data,
    std::size_t size,
    DataLayout layout=DataLayout::ROW_MAJOR
) 
```




<hr>



### function DataView [3/5]

_Construct from raw pointer and size (const version)_ 
```C++
inline emulator::DataView::DataView (
    const T * data,
    std::size_t size,
    DataLayout layout=DataLayout::ROW_MAJOR
) 
```




<hr>



### function DataView [4/5]

_Construct from std::vector (non-owning)_ 
```C++
inline explicit emulator::DataView::DataView (
    std::vector< T > & vec,
    DataLayout layout=DataLayout::ROW_MAJOR
) 
```




<hr>



### function DataView [5/5]

_Construct from const std::vector (non-owning)_ 
```C++
inline explicit emulator::DataView::DataView (
    const std::vector< T > & vec,
    DataLayout layout=DataLayout::ROW_MAJOR
) 
```




<hr>



### function begin [1/2]

_Iterator support._ 
```C++
inline T * emulator::DataView::begin () 
```




<hr>



### function begin [2/2]

```C++
inline const T * emulator::DataView::begin () const
```




<hr>



### function data [1/2]

_Get raw pointer to data._ 
```C++
inline T * emulator::DataView::data () 
```




<hr>



### function data [2/2]

```C++
inline const T * emulator::DataView::data () const
```




<hr>



### function empty 

_Check if view is empty._ 
```C++
inline bool emulator::DataView::empty () const
```




<hr>



### function end [1/2]

```C++
inline T * emulator::DataView::end () 
```




<hr>



### function end [2/2]

```C++
inline const T * emulator::DataView::end () const
```




<hr>



### function layout 

_Get memory layout._ 
```C++
inline DataLayout emulator::DataView::layout () const
```




<hr>



### function operator[] 

_Element access (bounds checking in debug builds)_ 
```C++
inline T & emulator::DataView::operator[] (
    std::size_t idx
) 
```




<hr>



### function operator[] 

```C++
inline const T & emulator::DataView::operator[] (
    std::size_t idx
) const
```




<hr>



### function size 

_Get number of elements._ 
```C++
inline std::size_t emulator::DataView::size () const
```




<hr>

------------------------------
The documentation for this class was generated from the following file `components/emulator_comps/common/src/data_view.hpp`

