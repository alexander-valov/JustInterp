#include <vector>
#include <cmath>
#include <array>

#include "doctest.h"
#include "JustInterp/JustInterp.hpp"

TEST_SUITE("TableInterpolatorTest") {

TEST_CASE("TableInterpolatorTest.check_node_points") {
    std::vector<double> x_data{-3.784, -2.456, -0.134, 0, 3.5643, 28.4588};
    std::vector<std::vector<double>> y_arrays{{
        {-973.234, -64.234, -32.98374, -0.000234, 0, 3.546, 5.8745, 8576.8575},
        {-10000, 0, 10000},
        {228.282, 387.8363, 847453.847437},
        {-0.0002663, -0.000013456, 0, 0.0000056456},
        {-97456.9345, -3423, -0.00001},
        {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20},
    }};

    /* define test function */
    auto test_func = [](double x, double y) { return std::log(1 + std::abs(x + y)) * std::exp(std::cos(x) * std::sin(y)); };

    std::vector<std::vector<double>> z_arrays(x_data.size());
    for (std::size_t i = 0; i < x_data.size(); i++) {
        std::vector<double> z_array;
        for (const auto& y : y_arrays[i]) {
            z_array.push_back(test_func(x_data[i], y));
        }
        z_arrays[i] = z_array;
    }

    JustInterp::TableInterpolator<double> interpolator(x_data, y_arrays, z_arrays);

    for (std::size_t ix = 0; ix < x_data.size(); ix++) {
        for (std::size_t iy = 0; iy < y_arrays[ix].size(); iy++) {
            CHECK(interpolator(x_data[ix], y_arrays[ix][iy]) == doctest::Approx(test_func(x_data[ix], y_arrays[ix][iy])));
        }
    }
}

TEST_CASE("TableInterpolatorTest.bilinear_function") {
    /* define bilinear function */
    double a_coeff = 17.837;
    double b_coeff = 25.524;
    double c_coeff = -69.228;
    auto linear_func = [&](double x, double y) { return a_coeff * x + b_coeff * y + c_coeff; };

    /* mesh data and values */
    std::vector<double> x_data{-3.784, -2.456, -0.134, 0, 3.5643, 28.4588};
    std::vector<std::vector<double>> y_arrays{{
        {-973.234, -64.234, -32.98374, -0.000234, 0, 3.546, 5.8745, 8576.8575},
        {-10000, 0, 10000},
        {228.282, 387.8363, 847453.847437},
        {-0.0002663, -0.000013456, 0, 0.0000056456},
        {-97456.9345, -3423, -0.00001},
        {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20},
    }};
    std::vector<std::vector<double>> z_arrays(x_data.size());
    for (std::size_t i = 0; i < x_data.size(); i++) {
        std::vector<double> z_array;
        for (const auto& y : y_arrays[i]) {
            z_array.push_back(linear_func(x_data[i], y));
        }
        z_arrays[i] = z_array;
    }

    JustInterp::TableInterpolator<double> interpolator(x_data, y_arrays, z_arrays);

    /* check mid-points along y-axis */
    for (std::size_t ix = 0; ix < x_data.size(); ix++) {
        for (std::size_t iy = 0; iy < y_arrays[ix].size() - 1; iy++) {
            double y_mid = 0.5 * (y_arrays[ix][iy] + y_arrays[ix][iy + 1]);
            CHECK(interpolator(x_data[ix], y_mid) == doctest::Approx(linear_func(x_data[ix], y_mid)));
        }
    }

    /* check mid-points along x-axis and intermediate points along y-axis */
    for (std::size_t ix = 0; ix < x_data.size() - 1; ix++) {
        std::array<double, 2> interval_a{y_arrays[ix].front(), y_arrays[ix].back()};
        std::array<double, 2> interval_b{y_arrays[ix + 1].front(), y_arrays[ix + 1].back()};
        /* intersection of neighboring intervals along y-axis */
        if (interval_b[0] <= interval_a[1] && interval_a[0] <= interval_b[1]) {
            double ymin = std::max(interval_a[0], interval_b[0]);
            double ymax = std::min(interval_a[1], interval_b[1]);

            const int n_points = 5;
            double x = 0.5 * (x_data[ix] + x_data[ix + 1]);
            for (int iy = 0; iy < n_points; iy++) {
                double y = ymin + iy * (ymax - ymin) / (n_points - 1);
                CHECK(interpolator(x, y) == doctest::Approx(linear_func(x, y)));
            }
        }
    }
}

TEST_CASE("TableInterpolatorTest.extrapolation_test") {
    /* define bilinear function */
    double a_coeff = 17.837;
    double b_coeff = 25.524;
    double c_coeff = -69.228;
    auto linear_func = [&](double x, double y) { return a_coeff * x + b_coeff * y + c_coeff; };

    /* mesh data and values */
    std::vector<double> x_data{-1, 0, 1};
    std::vector<std::vector<double>> y_arrays{{
        {-2, 0, 2},
        {-1, 1},
        {-2, 0, 2}
    }};
    std::vector<std::vector<double>> z_arrays(x_data.size());
    for (std::size_t i = 0; i < x_data.size(); i++) {
        std::vector<double> z_array;
        for (const auto& y : y_arrays[i]) {
            z_array.push_back(linear_func(x_data[i], y));
        }
        z_arrays[i] = z_array;
    }

    /* define control points */
    std::vector<std::array<double, 2>> control_points{{
        {0, -2},
        {0, 2},
        {0, -3},
        {0, 3},
        {-2, -3},
        {2, -3},
        {-2, 0},
        {2, 0},
        {-2, 3},
        {2, 3}
    }};

    /* Linear extrapolation */
    JustInterp::TableInterpolator<double, JustInterp::LinearExtrapolation> interpolator_linear(x_data, y_arrays, z_arrays);
    for (const auto& point : control_points) {
        CHECK(interpolator_linear(point[0], point[1]) == doctest::Approx(linear_func(point[0], point[1])));
    }

    /* Constant extrapolation */
    std::vector<std::array<double, 2>> const_points{{
        {0, -1},
        {0, 1},
        {0, -1},
        {0, 1},
        {-1, -2},
        {1, -2},
        {-1, 0},
        {1, 0},
        {-1, 2},
        {1, 2}
    }};
    JustInterp::TableInterpolator<double, JustInterp::ConstantExtrapolation> interpolator_const(x_data, y_arrays, z_arrays);
    for (std::size_t i = 0; i < control_points.size(); i++) {
        CHECK(interpolator_const(control_points[i][0], control_points[i][1]) == doctest::Approx(linear_func(const_points[i][0], const_points[i][1])));
    }
}

TEST_CASE("TableInterpolatorTest.table_bilinear_comparison") {
    /* define bilinear function ceofficients */
    double a_coeff = 17.837;
    double b_coeff = 25.524;
    double c_coeff = -69.228;

    /* define domain */
    double x_min = -74.736;
    double x_max = 13.836;
    double y_min = -0.926;
    double y_max = 17.746;

    /* evenly spaced nodes */
    int nx = 27;
    std::vector<double> x_1d(nx);
    double dx = (x_max - x_min) / (nx - 1);
    for (int i = 0; i < nx; i++) {
        x_1d[i] = x_min + i * dx;
    }

    int ny = 5;
    std::vector<double> y_1d(ny);
    double dy = (y_max - y_min) / (ny - 1);
    for (int i = 0; i < ny; i++) {
        y_1d[i] = y_min + i * dy;
    }

    std::vector<std::vector<double>> y_arrays, z_arrays;
    for (int i = 0; i < nx; i++) {
        y_arrays.push_back(y_1d);
    }

    auto linear_func = [&](double x, double y) { return a_coeff * x + b_coeff * y + c_coeff; };

    /* YMajor: calculate node values */
    std::vector<double> z_ymajor(nx * ny);
    for (int ix = 0; ix < nx; ix++) {
        std::vector<double> z_data;
        for (int iy = 0; iy < ny; iy++) {
            z_ymajor[ix * ny + iy] = linear_func(x_1d[ix], y_1d[iy]);
            z_data.push_back(linear_func(x_1d[ix], y_1d[iy]));
        }
        z_arrays.push_back(z_data);
    }

    JustInterp::BilinearInterpolator<double, JustInterp::YMajor> interp_bilinear(x_1d, y_1d, z_ymajor);
    JustInterp::TableInterpolator<double> interp_table(x_1d, y_arrays, z_arrays);

    for (int ix = 0; ix < nx - 1; ix++) {
        double x = 0.5 * (x_1d[ix] + x_1d[ix + 1]);
        for (int iy = 0; iy < ny - 1; iy++) {
            double y = 0.5 * (y_1d[iy] + y_1d[iy + 1]);
            CHECK(interp_table(x, y) == doctest::Approx(interp_bilinear(x, y)));
        }
    }
}

}