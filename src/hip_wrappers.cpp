#include <hip/hip_runtime.h>

/*
Wrapper functions which deals with error handling so that the error
handling does not need to happen inside the actual code, and so that I
can get the pesky unused result warnings to shut up.
*/

namespace hip_wrappers
{
void hipMemcpy(void* dst, const void* src, size_t sizeBytes, hipMemcpyKind kind)
{
    hipError_t result = ::hipMemcpy(dst, src, sizeBytes, kind);
    if (result != hipSuccess)
    {
        throw std::runtime_error(hipGetErrorString(result));
    }
}

void hipDeviceSynchronize()
{
    hipError_t result = ::hipDeviceSynchronize();
    if (result != hipSuccess)
    {
        throw std::runtime_error(hipGetErrorString(result));
    }
}

void hipFree(void* ptr)
{
    hipError_t result = ::hipFree(ptr);
    if (result != hipSuccess)
    {
        throw std::runtime_error(hipGetErrorString(result));
    }
}

void hipMemset(void *dst, int value, size_t sizeBytes)
{
    hipError_t result = ::hipMemset(dst, value, sizeBytes);
    if (result != hipSuccess)
    {
        throw std::runtime_error(hipGetErrorString(result));
    }
}

} // namespace hip_wrappers