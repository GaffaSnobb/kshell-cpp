#ifndef HIP_WRAPPERS_HPP
#define HIP_WRAPPERS_HPP

#include <hip/hip_runtime.h>

using std::vector;

namespace hip_wrappers
{
template <typename T>
void hipMalloc(T** devPtr, size_t size);
void hipMemcpy(void* dst, const void* src, size_t sizeBytes, hipMemcpyKind kind);
void hipDeviceSynchronize();
void hipFree(void* ptr);

template <typename T1>
void hipMemcpyToSymbol(T1& to_become_symbol, vector<T1> src);
template <typename T1>
void hipMemcpyToSymbol(T1& to_become_symbol, vector<int16_t> src);
} // namespace hip_wrappers

#include "hip_wrappers.tpp"
#endif