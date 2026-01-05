

# Class emulator::EmulatorContext



[**ClassList**](annotated.md) **>** [**emulator**](namespaceemulator.md) **>** [**EmulatorContext**](classemulator_1_1EmulatorContext.md)



_Singleton context for managing emulator component instances._ [More...](#detailed-description)

* `#include <emulator_context.hpp>`





































## Public Functions

| Type | Name |
| ---: | :--- |
|  void | [**clean\_up**](#function-clean_up) () <br>_Remove all objects from the context._  |
|  T & | [**create**](#function-create) (Args &&... args) <br>_Create and register a new component instance._  |
|  const T & | [**get**](#function-get) () const<br>_Get a const reference to an existing component._  |
|  T & | [**getNonConst**](#function-getnonconst) () <br>_Get a non-const reference to an existing component._  |
|  bool | [**has**](#function-has) () const<br>_Check if a component of type T exists._  |


## Public Static Functions

| Type | Name |
| ---: | :--- |
|  [**EmulatorContext**](classemulator_1_1EmulatorContext.md) & | [**singleton**](#function-singleton) () <br>_Get the singleton instance._  |


























## Detailed Description


## Public Functions Documentation




### function clean\_up 

_Remove all objects from the context._ 
```C++
inline void emulator::EmulatorContext::clean_up () 
```



Should be called during shutdown to release resources. 


        

<hr>



### function create 

_Create and register a new component instance._ 
```C++
template<typename T, typename... Args>
inline T & emulator::EmulatorContext::create (
    Args &&... args
) 
```



Creates an instance of type T with the given constructor arguments and stores it in the context.




**Template parameters:**


* `T` Component type to create 
* `Args` Constructor argument types 



**Parameters:**


* `args` Arguments to pass to T's constructor 



**Returns:**

Reference to the newly created instance 




**Exception:**


* `std::runtime_error` if an instance of T already exists 




        

<hr>



### function get 

_Get a const reference to an existing component._ 
```C++
template<typename T>
inline const T & emulator::EmulatorContext::get () const
```





**Template parameters:**


* `T` Component type to retrieve 



**Returns:**

Const reference to the component 




**Exception:**


* `std::runtime_error` if no instance of T exists 




        

<hr>



### function getNonConst 

_Get a non-const reference to an existing component._ 
```C++
template<typename T>
inline T & emulator::EmulatorContext::getNonConst () 
```





**Template parameters:**


* `T` Component type to retrieve 



**Returns:**

Reference to the component 




**Exception:**


* `std::runtime_error` if no instance of T exists 




        

<hr>



### function has 

_Check if a component of type T exists._ 
```C++
template<typename T>
inline bool emulator::EmulatorContext::has () const
```





**Template parameters:**


* `T` Component type to check 



**Returns:**

true if an instance of T exists in the context 





        

<hr>
## Public Static Functions Documentation




### function singleton 

_Get the singleton instance._ 
```C++
static inline EmulatorContext & emulator::EmulatorContext::singleton () 
```





**Returns:**

Reference to the global [**EmulatorContext**](classemulator_1_1EmulatorContext.md) 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `components/emulator_comps/common/src/emulator_context.hpp`

