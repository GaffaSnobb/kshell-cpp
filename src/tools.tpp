#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <type_traits>
#include <numeric>
#include <omp.h>
#include <stdint.h>

using std::cout;
using std::endl;
using std::vector;
using std::string;

template <typename T1, typename T2>
void print_flattened_2d_array(const T1* arr, const T2 size)
{
    for (size_t row_idx = 0; row_idx < size; row_idx++)
    {
        for (size_t col_idx = 0; col_idx < size; col_idx++)
        {
            cout << arr[row_idx*size + col_idx] << ", ";
        }
        cout << endl;
    }
    cout << endl;
}

template <typename T1, typename T2>
void write_flattened_2d_array_to_file(T1* arr, const T2 size, const std::string& filename)
{   
    std::ofstream file(filename);

    if (!file.is_open())
    {
        std::cerr << "Error opening file for writing." << std::endl;
        return;
    }

    for (size_t row_idx = 0; row_idx < size; row_idx++)
    {
        /*
        Copies the values from the upper triangle to the lower triangle.
        */
        for (size_t col_idx = row_idx + 1; col_idx < size; col_idx++)
        {
            arr[col_idx*size + row_idx] = arr[row_idx*size + col_idx];
        }
    }

    size_t max_width = 15;

    for (size_t row_idx = 0; row_idx < size; row_idx++)
    {
        for (size_t col_idx = 0; col_idx < size; col_idx++)
        {
            file << std::left << std::setw(max_width) << arr[row_idx*size + col_idx];
        }
        file << '\n';
    }
    file << '\n';

    file.close();
}


template <typename T1, typename T2, typename T3>
bool compare_arrays(T1* arr1, T2* arr2, T3 size)
{
    for (size_t i = 0; i < size; i++) if (std::abs(arr1[i] - arr2[i]) > 1e-13) return false;
    return true;
}

template <typename T1, typename T2, typename T3>
bool compare_arrays_upper_triangle(T1* arr1, vector<T2> arr2, T3 size)
{
    for (size_t row_idx = 0; row_idx < size; row_idx++)
    {
        for (size_t col_idx = row_idx; col_idx < size; col_idx++)
        {
            size_t flat_idx = row_idx*size + col_idx;
            if (std::abs(arr1[flat_idx] - arr2[flat_idx]) > 1e-3)
            {   
                // cout << "flat_idx: " << flat_idx << ", arr1[flat_idx]: " << arr1[flat_idx] << ", arr2[flat_idx]: " << arr2[flat_idx] << ", " << std::abs(arr1[flat_idx] - arr2[flat_idx]) << endl;
                // cout << "LOL" << endl;
                // cout << (arr1[flat_idx] - arr2[flat_idx]) << endl;
                // exit(0);
                return false;
            }
        }
    }
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
void print_vector(const vector<T>& vec)
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
void print_vector(const vector<vector<T>>& nested_vector)
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
void print_vector(string name, const vector<T>& vec)
{   
    cout << name << " = ";
    for (const T& val : vec)
    {
        cout << val << " ";
    }
    cout << endl;
}

template <typename T>
void print(string name, T value)
{
    cout << name << " = " << value << endl;
    return;
}

template <typename T0, typename T1, typename T2, typename T3>
vector<T0> range(T1 start, T2 stop, T3 step)
{
    static_assert(
        std::is_integral<T1>::value and std::is_integral<T2>::value and std::is_integral<T3>::value,
        "Start, stop, and step must be integral types."
    );
    vector<T0> range;

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
double mean(vector<T> vec)
{
    return std::accumulate(vec.begin(), vec.end(), 0.0)/vec.size();
}

template <typename T0, typename T1, typename T2>
void print_loop_timer(vector<T0>& loop_timings, T1 idx, T2 n_iterations)
{
    double mean_time = mean(loop_timings)/1000;
    int32_t num_threads = omp_get_num_threads();
    
    cout << "\r[" << idx << " of ≈ " << (double)n_iterations/num_threads << "]";
    cout << " [loop time: ";
    cout << std::setfill(' ') << std::setw(5) << loop_timings.back()/1000. << " s";
    cout << " - mean loop time: ";
    cout << std::setfill(' ') << std::setw(10) << mean_time << " s]";
    cout << " [est. time left: ";
    cout << std::setfill(' ') << std::setw(10) << (n_iterations - (idx + 1))*mean_time << " s]" << std::flush;
}