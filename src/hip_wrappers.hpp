#ifndef HIP_WRAPPERS_HPP
#define HIP_WRAPPERS_HPP

#include <hip/hip_runtime.h>

namespace hip_wrappers
{
template <typename T>
void hipMalloc(T** devPtr, size_t size);
void hipMemcpy(void* dst, const void* src, size_t sizeBytes, hipMemcpyKind kind);
void hipDeviceSynchronize();
void hipFree(void* ptr);
} // namespace hip_wrappers

#include "hip_wrappers.tpp"
#endif