#ifndef HIP_WRAPPERS_HPP
#define HIP_WRAPPERS_HPP

#include <hip/hip_runtime.h>

using std::vector;

namespace hip_wrappers
{
void hipMemcpy(void* dst, const void* src, size_t sizeBytes, hipMemcpyKind kind);
void hipDeviceSynchronize();
void hipFree(void* ptr);
void hipMemset(void *dst, int value, size_t sizeBytes);
void hipGetLastError();

template <typename T>
void hipMalloc(T** devPtr, size_t size);

template <typename T1>
void hipMemcpyToSymbol(T1& to_become_symbol, vector<T1> src);

template <typename T1>
void hipMemcpyToSymbol(T1& to_become_symbol, vector<int16_t> src);

template <typename T1>
void hipMemcpyToSymbol(T1& dev_const_mem_ptr, vector<double> src);

template <typename T1>
void hipMemcpyToSymbol(T1& dev_const_mem_ptr, const uint16_t* src, const size_t size);
} // namespace hip_wrappers

#include "hip_wrappers.tpp"
#endif