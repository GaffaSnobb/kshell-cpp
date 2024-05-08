#include <vector>
#include <hip/hip_runtime.h>

using std::vector;

namespace hip_wrappers
{
template <typename T>
void hipMalloc(T** devPtr, size_t size)
{
    hipError_t result = ::hipMalloc(devPtr, size);
    if (result != hipSuccess)
    {
        throw std::runtime_error(hipGetErrorString(result));
    }
}

template <typename T1>
void hipMemcpyToSymbol(T1& dev_const_mem_ptr, vector<uint16_t> src)
{
    const size_t offset = 0;
    hipError_t result = ::hipMemcpyToSymbol(
        HIP_SYMBOL(dev_const_mem_ptr),
        src.data(),
        sizeof(uint16_t)*src.size(),
        offset,
        hipMemcpyHostToDevice
    );
    if (result != hipSuccess)
    {
        throw std::runtime_error(hipGetErrorString(result));
    }
}

template <typename T1>
void hipMemcpyToSymbol(T1& dev_const_mem_ptr, vector<int16_t> src)
{
    /*
    I coud have templated int16_t instead of overloading, but this way
    I'm sure that T1 and the data type in the vector is the same. I cant
    use T1 in vector<T1> because `dev_const_mem_ptr` is a constant size
    T1 array of type T1 [<size>], not T1*.
    */
    const size_t offset = 0;
    hipError_t result = ::hipMemcpyToSymbol(
        HIP_SYMBOL(dev_const_mem_ptr),
        src.data(),
        sizeof(int16_t)*src.size(),
        offset,
        hipMemcpyHostToDevice
    );
    if (result != hipSuccess)
    {
        throw std::runtime_error(hipGetErrorString(result));
    }
}

template <typename T1>
void hipMemcpyToSymbol(T1& dev_const_mem_ptr, const uint16_t* src, const size_t size)
{
    /*
    Have to use a template for the device constant pointer because the
    type of the array is not just a pointer of the data type it
    contains, the type is specific to the length of the array.
    */
    const size_t offset = 0;
    hipError_t result = ::hipMemcpyToSymbol(
        HIP_SYMBOL(dev_const_mem_ptr),
        src,
        size,
        offset,
        hipMemcpyHostToDevice
    );
    if (result != hipSuccess)
    {
        throw std::runtime_error(hipGetErrorString(result));
    }
}
} // namespace hip_wrappers