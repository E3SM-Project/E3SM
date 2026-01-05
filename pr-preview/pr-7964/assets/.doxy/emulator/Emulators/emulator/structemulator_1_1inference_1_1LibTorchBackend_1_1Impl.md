

# Struct emulator::inference::LibTorchBackend::Impl



[**ClassList**](annotated.md) **>** [**emulator**](namespaceemulator.md) **>** [**inference**](namespaceemulator_1_1inference.md) **>** [**LibTorchBackend**](classemulator_1_1inference_1_1LibTorchBackend.md) **>** [**Impl**](structemulator_1_1inference_1_1LibTorchBackend_1_1Impl.md)



_Private implementation details for_ [_**LibTorchBackend**_](classemulator_1_1inference_1_1LibTorchBackend.md) _._[More...](#detailed-description)






















## Public Attributes

| Type | Name |
| ---: | :--- |
|  torch::Device | [**device**](#variable-device)   = `torch::kCPU`<br>_Execution device (CPU or CUDA)_  |
|  torch::ScalarType | [**dtype**](#variable-dtype)   = `torch::kFloat32`<br>_Model precision._  |
|  torch::jit::script::Module | [**model**](#variable-model)  <br>_Loaded TorchScript model._  |
|  bool | [**model\_loaded**](#variable-model_loaded)   = `false`<br>_Whether model was successfully loaded._  |












































## Detailed Description


Uses PIMPL idiom to hide LibTorch types from the header, avoiding the need to include torch headers in dependent code. 


    
## Public Attributes Documentation




### variable device 

_Execution device (CPU or CUDA)_ 
```C++
torch::Device emulator::inference::LibTorchBackend::Impl::device;
```




<hr>



### variable dtype 

_Model precision._ 
```C++
torch::ScalarType emulator::inference::LibTorchBackend::Impl::dtype;
```




<hr>



### variable model 

_Loaded TorchScript model._ 
```C++
torch::jit::script::Module emulator::inference::LibTorchBackend::Impl::model;
```




<hr>



### variable model\_loaded 

_Whether model was successfully loaded._ 
```C++
bool emulator::inference::LibTorchBackend::Impl::model_loaded;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `components/emulator_comps/common/src/inference/libtorch_backend.cpp`

