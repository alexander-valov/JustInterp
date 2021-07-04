# JustInterp

JustInterp is the simple lightweight header-only C++17 interpolation library.

- [Integration](#integration)
- [Linear Interpolation](#linear-interpolation)
  - [Input data format](#input-data-format)
  - [Example of usage](#example-of-usage)
  - [How to set data](#how-to-set-data)
  - [Extrapolation](#extrapolation)
- [Bilinear Interpolation](#bilinear-interpolation)
  - [Input data format](#input-data-format-1)
  - [Example of usage](#example-of-usage-1)
  - [How to set data](#how-to-set-data-1)
  - [Extrapolation](#extrapolation-1)
- [Table Interpolation](#table-interpolation)
  - [Input data format](#input-data-format-2)
  - [Example of usage](#example-of-usage-2)
  - [How to set data](#how-to-set-data-2)
  - [Extrapolation](#extrapolation-2)

## Integration

JustInterp is the single-header library, hence, you just need to copy `JustInterp.hpp` from `include_single/JustInterp` or release tab to your project.

```c++
#include "JustInterp/JustInterp.hpp"
```

Set the necessary options to enable C++17, for example:

- Compiler flag: `-std=c++17` for GCC and Clang, `/std:c++17` for MSVC
- Using CMake:
    - Specify compile features for specific target: `target_compile_features(<target> <PRIVATE|PUBLIC|INTERFACE> cxx_std_17)`
    - Set global property: `set(CMAKE_CXX_STANDARD 17 CACHE STRING "The C++ standard to use")`

## Linear Interpolation

Piece-wise linear interpolation (https://en.wikipedia.org/wiki/Linear_interpolation).

### Input data format

`RealContainer` is real-type iterable containers i.e. should support `begin(), end(), data(), size()` methods, for example, `std::vector` or `std::array`.

- `RealContainer x`: data points
- `RealContainer y`: function values

### Example of usage

```c++
#include "JustInterp/JustInterp.hpp"

std::vector<double> x, y;

// ...
// set points x and values y
// ...

// Define linear interpolator and set data
JustInterp::LinearInterpolator<double> interpolator(x, y);

// Interpolate at the specific point
double result = interpolator(3.5);

// Interpolate at the array of points
std::vector<double> point_array{1.5, 2.5, 3.5};
std::vector<double> results = interpolator(point_array);
```

### How to set data

1. Using constructor:
    ```c++
    JustInterp::LinearInterpolator<double> interpolator(x, y);
    ```
    `x.size()` and `y.size()` must be equal.
2. Using `SetData` method:
    ```c++
    JustInterp::LinearInterpolator<double> interpolator;
    interpolator.SetData(x, y);
    ```
    `x.size()` and `y.size()` must be equal.
3. Using `SetData` method for data pointers:
    ```c++
    JustInterp::LinearInterpolator<double> interpolator;
    interpolator.SetData(x.size(), x.data(), y.data());
    ```
    It is assumed that arrays `x.data()` and `y.data()` of the same size of `x.size()`.

### Extrapolation

There are two avaliable types of [extrapolation](https://en.wikipedia.org/wiki/Extrapolation):
1. (Default) `JustInterp::ConstantExtrapolation`
    ```c++
    JustInterp::LinearInterpolator<double, JustInterp::ConstantExtrapolation> interpolator;
    /* or */
    JustInterp::LinearInterpolator<double> interpolator;
    ```
2. `JustInterp::LinearExtrapolation`
    ```c++
    JustInterp::LinearInterpolator<double, JustInterp::LinearExtrapolation> interpolator;
    ```

![Avaliable extrapolation types](./doc/img/linear_interpolator_extrapolation.png)

## Bilinear Interpolation

Bilinear interpolation (https://en.wikipedia.org/wiki/Bilinear_interpolation).

### Input data format

`RealContainer` is real-type iterable containers i.e. should support `begin(), end(), data(), size()` methods, for example, `std::vector` or `std::array`.

- `RealContainer x_1d` - grid coordinates along x-axis
- `RealContainer y_1d` - grid coordinates along y-axis
- `RealContainer z_all` - values at all grid points. 1D array of size `x_1d.size() * y_1d.size()`.

There are two avaliable storage orders for `z_all`:
- `YMajor` - default storage order, see left figure.
    ```c++
    JustInterp::BilinearInterpolator<double, JustInterp::YMajor> interpolator(x_1d, y_1d, z_ymajor);
    ```
    or equivalent
    ```c++
    JustInterp::BilinearInterpolator<double> interpolator(x_1d, y_1d, z_ymajor);
    ```
- `XMajor` - see right figure
    ```c++
    JustInterp::BilinearInterpolator<double, JustInterp::XMajor> interpolator(x_1d, y_1d, z_xmajor);
    ```

![Avaliable 2D storage orders](./doc/img/ymajor_vs_xmajor.png)

### Example of usage

```c++
#include "JustInterp/JustInterp.hpp"

std::vector<double> x_1d, y_1d;
std::vector<double> z_ymajor;

// ...
// set grid points x_1d, y_1d along corresponding direction
// ...

// ...
// set z_ymajor values according to YMajor order
// ...

// Define bilinear interpolator and set data
JustInterp::BilinearInterpolator<double> interpolator(x_1d, y_1d, z_major);

// Interpolate at the specific point
double result = interpolator(3.5, 2.0);

// Interpolate at the array of points
std::vector<double> points_x{ 0.3,  0.6, 1.5};
std::vector<double> points_y{-1.3, -0.4, 3.2};
// result is the array of size points_x.size() = points_y.size()
std::vector<double> results = interpolator(points_x, points_y);

// Interpolate to destination grid
std::vector<double> x_1d_dest, y_1d_dest;
// ...
// fill destination grid
// ...
// results_dest is the array of size x_1d_dest.size() * y_1d_dest.size()
// storage order is the same to interpolator's storage order
std::vector<double> results_dest = interpolator.GridInterpolation(x_1d_dest, y_1d_dest);
```

### How to set data

1. Using constructor:
    ```c++
    JustInterp::BilinearInterpolator<double, JustInterp::YMajor> interpolator(x_1d, y_1d, z_ymajor);
    ```
    or
    ```c++
    JustInterp::BilinearInterpolator<double, JustInterp::XMajor> interpolator(x_1d, y_1d, z_xmajor);
    ```
    Size of `z_ymajor` or `z_xmajor` must be equal to `x_1d.size() * y_1d.size()`.
2. Using `SetData` method:
    ```c++
    JustInterp::BilinearInterpolator<double, JustInterp::YMajor> interpolator;
    interpolator.SetData(x_1d, y_1d, z_ymajor);
    ```
    or
    ```c++
    JustInterp::BilinearInterpolator<double, JustInterp::XMajor> interpolator;
    interpolator.SetData(x_1d, y_1d, z_xmajor);
    ```
3. Using `SetData` method for data pointers:
    ```c++
    JustInterp::BilinearInterpolator<double, JustInterp::YMajor> interpolator;
    interpolator.SetData(x_1d.size(), y_1d.size(), x_1d.data(), y_1d.data(), z_ymajor.data());
    ```
    or
    ```c++
    JustInterp::BilinearInterpolator<double, JustInterp::XMajor> interpolator;
    interpolator.SetData(x_1d.size(), y_1d.size(), x_1d.data(), y_1d.data(), z_xmajor.data());
    ```

### Extrapolation
Only bilinear extrapolation based on data from the nearest grid cell is implemented.

## Table Interpolation

The table interpolation is 2D bilinear interpolation for data where the rows do not have the same number of columns, i.e., when the data for interpolation does not constitute a two-dimensional rectangular grid (see figure). If the data is a rectangular grid, then the proper solution is to use `BilinearInterpolator`, even though in this case `TableInterpolator` and `BilinearInterpolator` give the same result.

The interpolation algorithm looks as follow:
1. For given point `(x*, y*)` find neighboring `x_1d[i]` and `x_1d[i + 1]`.
2. For given `y*` perform linear interpolation along y-axis for `i` and `i + 1` to yield values `f1` and `f2`.
3. Use data `(x_1d[i], f1)` and `(x_1d[i + 1], f2)` for linear interpolation along x-axis.

![Table interpolation algorithm](./doc/img/table_interpolator.png)

### Input data format

`RealContainer` is real-type iterable containers i.e. should support `begin(), end(), data(), size()` methods, for example, `std::vector` or `std::array`.

Input data should be stored as follows:
- `RealContainer x_1d`: 1D array of data points along x-axis.
- `std::vector<RealContainer> y_1d_arrays`: vector of 1D arrays of data points along y-axis for each x-value from `x_1d` array. `y_1d_arrays.size()` must be equal `x_1d.size()`.
- `std::vector<RealContainer> z_1d_arrays`: vector of 1D arrays of function values along y-axis for each x-value from `x_1d` array. `z_1d_arrays.size()` must be equal `x_1d.size()` and `z_1d_arrays[i].size()` must be equal `y_1d_arrays[i].size()`.

### Example of usage

```c++
#include "JustInterp/JustInterp.hpp"

std::vector<double> x_1d;
std::vector<std::vector<double>> y_1d_arrays;
std::vector<std::vector<double>> z_1d_arrays;

// ...
// set grid points x_1d, y_1d_arrays
// ...

// ...
// set z_1d_arrays
// ...

// Define table interpolator and set data
JustInterp::TableInterpolator<double> interpolator(x_1d, y_1d_arrays, z_1d_arrays);

// Interpolate at the specific point
double result = interpolator(3.5, 2.0);

// Interpolate at the array of points
std::vector<double> points_x{ 0.3,  0.6, 1.5};
std::vector<double> points_y{-1.3, -0.4, 3.2};
// result is the array of size points_x.size() = points_y.size()
std::vector<double> results = interpolator(points_x, points_y);
```

### How to set data

1. Using constructor:
    ```c++
    JustInterp::TableInterpolator<double> interpolator(x_1d, y_1d_arrays, z_1d_arrays);
    ```
2. Using `SetData` method:
    ```c++
    JustInterp::TableInterpolator<double> interpolator;
    interpolator.SetData(x_1d, y_1d_arrays, z_1d_arrays);
    ```

### Extrapolation

There are two avaliable types of [extrapolation](https://en.wikipedia.org/wiki/Extrapolation):
1. (Default) `JustInterp::ConstantExtrapolation`
    ```c++
    JustInterp::TableInterpolator<double, JustInterp::ConstantExtrapolation> interpolator;
    /* or */
    JustInterp::TableInterpolator<double> interpolator;
    ```
2. `JustInterp::LinearExtrapolation`
    ```c++
    JustInterp::TableInterpolator<double, JustInterp::LinearExtrapolation> interpolator;
    ```

First, extrapolation is performed along the y-axis, then extrapolation is performed along the x-axis.