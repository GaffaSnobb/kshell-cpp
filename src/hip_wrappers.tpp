#include <hip/hip_runtime.h>

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
} // namespace hip_wrappers