#include <iostream>
#include <vector>
#include <type_traits>

using std::cout;
using std::endl;

template <typename T>
void print_vector(const std::vector<T>& vec) {
    for (const T& val : vec) {
        cout << val << " ";
    }
    cout << endl;
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