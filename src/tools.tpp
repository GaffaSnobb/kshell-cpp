#include <iostream>
#include <iomanip>
#include <vector>
#include <type_traits>
#include <numeric>
#include <omp.h>

using std::cout;
using std::endl;

template <typename T1, typename T2>
void print_flattened_2d_array(T1* arr, T2 size)
{
    for (int row_idx = 0; row_idx < size; row_idx++)
    {
        for (int col_idx = 0; col_idx < size; col_idx++)
        {
            cout << arr[row_idx*size + col_idx] << ", ";
        }
        cout << endl;
    }
    cout << endl;
}

template <typename T1, typename T2>
bool compare_arrays(T1* arr1, T1* arr2, T2 size)
{
    for (int i = 0; i < size; i++) if (std::abs(arr1[i] - arr2[i]) > 1e-13) return false;
    return true;
}

template <typename T>
void print_bit_representation(const T& value)
{
    static_assert(std::is_integral<T>::value, "T must be an integral type.");
    cout << "Value: " << value << endl;
    cout << "Number of bits: " << sizeof(T) * 8 << endl;
    cout << "Bit representation: " << std::bitset<sizeof(T) * 8>(value) << endl;
}

template <typename T>
void print_vector(const std::vector<T>& vec)
{   
    bool first_value = true;
    cout << "[";
    for (const T& val : vec)
    {
        if (not first_value) cout << ", ";
        cout << val;
        first_value = false;
    }
    cout << "]" << endl;
}

template <typename T>
void print_vector(const std::vector<std::vector<T>>& nested_vector)
{
    for (const auto& inner_vector : nested_vector)
    {
        for (const auto& element : inner_vector)
        {
            cout << element << " ";
        }
        cout << endl;
    }
    // cout << "(";
    // for (const auto& inner_vector : nested_vector)
    // {
    //     cout << "(";
    //     for (const auto& element : inner_vector)
    //     {
    //         cout << element << ", ";
    //     }
    //     cout << "),";
    // }
    // cout << ")" << endl;
}

template <typename T>
void print_vector(std::string name, const std::vector<T>& vec)
{   
    cout << name << " = ";
    for (const T& val : vec)
    {
        cout << val << " ";
    }
    cout << endl;
}

template <typename T>
void print(std::string name, T value)
{
    cout << name << " = " << value << endl;
    return;
}

template <typename T0, typename T1, typename T2, typename T3>
std::vector<T0> range(T1 start, T2 stop, T3 step)
{
    static_assert(
        std::is_integral<T1>::value and std::is_integral<T2>::value and std::is_integral<T3>::value,
        "Start, stop, and step must be integral types."
    );
    std::vector<T0> range;

    if (step == 0)
    {
        return range;
    }
    if (step > 0)
    {
        for (T0 i = start; i < stop; i += step)
        {
            range.push_back(i);
        }
    }
    else
    {
        for (T0 i = start; i > stop; i += step)
        {
            range.push_back(i);
        }
    }
    return range;
}

template <typename T>
double mean(std::vector<T> vec)
{
    return std::accumulate(vec.begin(), vec.end(), 0.0)/vec.size();
}

template <typename T0, typename T1, typename T2>
void print_loop_timer(std::vector<T0>& loop_timings, T1 idx, T2 n_iterations)
{
    double mean_time = mean(loop_timings)/1000;
    int num_threads = omp_get_num_threads();
    
    cout << "\r[" << idx << " of ≈ " << (double)n_iterations/num_threads << "]";
    cout << " [loop time: ";
    cout << std::setfill(' ') << std::setw(5) << loop_timings.back()/1000. << " s";
    cout << " - mean loop time: ";
    cout << std::setfill(' ') << std::setw(10) << mean_time << " s]";
    cout << " [est. time left: ";
    cout << std::setfill(' ') << std::setw(10) << (n_iterations - (idx + 1))*mean_time << " s]" << std::flush;
}