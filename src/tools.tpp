#include <iostream>
#include <vector>

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

template <typename T>
std::vector<T> range(int start, int stop, int step) {
    std::vector<T> range;

    if (step == 0) {
        return range;
    }

    if (step > 0) {
        for (T i = start; i < stop; i += step) {
            range.push_back(i);
        }
    } else {
        for (T i = start; i > stop; i += step) {
            range.push_back(i);
        }
    }

    return range;
}