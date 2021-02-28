#include <vector>
#include <array>
#include <cmath>

#include "doctest.h"
#include "JustInterp/BilinearInterpolator.hpp"

TEST_SUITE("BilinearInterpolatorTest") {

TEST_CASE("check_node_points") {
    /* define domain */
    double x_min = -98.983;
    double x_max = 273.635;
    double y_min = -83.726;
    double y_max = 13.259;

    /* evenly spaced nodes */
    int nx = 73;
    std::vector<double> x_1d(nx);
    double dx = (x_max - x_min) / (nx - 1);
    for (int i = 0; i < nx; i++) {
        x_1d[i] = x_min + i * dx;
    }

    int ny = 69;
    std::vector<double> y_1d(ny);
    double dy = (y_max - y_min) / (ny - 1);
    for (int i = 0; i < ny; i++) {
        y_1d[i] = y_min + i * dy;
    }

    /* define test function */
    auto test_func = [](double x, double y) { return std::log(1 + std::abs(x + y)) * std::exp(std::cos(x) * std::sin(y)); };

    /* YMajor: calculate node values */
    std::vector<double> z_ymajor(nx * ny);
    for (int ix = 0; ix < nx; ix++) {
        for (int iy = 0; iy < ny; iy++) {
            z_ymajor[ix * ny + iy] = test_func(x_1d[ix], y_1d[iy]);
        }
    }

    /* XMajor: calculate node values */
    std::vector<double> z_xmajor(nx * ny);
    for (int iy = 0; iy < ny; iy++) {
        for (int ix = 0; ix < nx; ix++) {
            z_xmajor[iy * nx + ix] = test_func(x_1d[ix], y_1d[iy]);
        }
    }

    auto interpYMajor = JustInterp::BilinearInterpolator<double, JustInterp::YMajor>(x_1d, y_1d, z_ymajor);
    auto interpXMajor = JustInterp::BilinearInterpolator<double, JustInterp::XMajor>(x_1d, y_1d, z_xmajor);

    /* check YMajor nodes values */
    for (int ix = 0; ix < nx; ix++) {
        for (int iy = 0; iy < ny; iy++) {
            CHECK(interpYMajor(x_1d[ix], y_1d[iy]) == doctest::Approx(test_func(x_1d[ix], y_1d[iy])));
        }
    }

    /* check XMajor nodes values */
    for (int ix = 0; ix < nx; ix++) {
        for (int iy = 0; iy < ny; iy++) {
            CHECK(interpXMajor(x_1d[ix], y_1d[iy]) == doctest::Approx(test_func(x_1d[ix], y_1d[iy])));
        }
    }

    /* YMajor: check GridInterpolation */
    auto result_ymajor = interpYMajor.GridInterpolation(x_1d, y_1d);
    for (int i = 0; i < nx * ny; i++) {
        CHECK(result_ymajor[i] == doctest::Approx(z_ymajor[i]));
    }

    /* XMajor: check GridInterpolation */
    auto result_xmajor = interpXMajor.GridInterpolation(x_1d, y_1d);
    for (int i = 0; i < nx * ny; i++) {
        CHECK(result_xmajor[i] == doctest::Approx(z_xmajor[i]));
    }
}

TEST_CASE("bilinear_function") {
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

    auto linear_func = [&](double x, double y) { return a_coeff * x + b_coeff * y + c_coeff; };

    /* YMajor: calculate node values */
    std::vector<double> z_ymajor(nx * ny);
    for (int ix = 0; ix < nx; ix++) {
        for (int iy = 0; iy < ny; iy++) {
            z_ymajor[ix * ny + iy] = linear_func(x_1d[ix], y_1d[iy]);
        }
    }

    /* XMajor: calculate node values */
    std::vector<double> z_xmajor(nx * ny);
    for (int iy = 0; iy < ny; iy++) {
        for (int ix = 0; ix < nx; ix++) {
            z_xmajor[iy * nx + ix] = linear_func(x_1d[ix], y_1d[iy]);
        }
    }

    auto interpYMajor = JustInterp::BilinearInterpolator<double, JustInterp::YMajor>(x_1d, y_1d, z_ymajor);
    auto interpXMajor = JustInterp::BilinearInterpolator<double, JustInterp::XMajor>(x_1d, y_1d, z_xmajor);

    /* test grid */
    int nx_test = 83;
    std::vector<double> x_1d_test(nx_test);
    double dx_test = (x_max - x_min) / (nx_test - 1);
    for (int i = 0; i < nx_test; i++) {
        x_1d_test[i] = x_min + i * dx_test;
    }

    int ny_test = 57;
    std::vector<double> y_1d_test(ny_test);
    double dy_test = (y_max - y_min) / (ny_test - 1);
    for (int i = 0; i < ny_test; i++) {
        y_1d_test[i] = y_min + i * dy_test;
    }

    /* check YMajor test values */
    for (int ix = 0; ix < nx_test; ix++) {
        for (int iy = 0; iy < ny_test; iy++) {
            CHECK(interpYMajor(x_1d_test[ix], y_1d_test[iy]) == doctest::Approx(linear_func(x_1d_test[ix], y_1d_test[iy])));
        }
    }

    /* check XMajor test values */
    for (int ix = 0; ix < nx_test; ix++) {
        for (int iy = 0; iy < ny_test; iy++) {
            CHECK(interpXMajor(x_1d_test[ix], y_1d_test[iy]) == doctest::Approx(linear_func(x_1d_test[ix], y_1d_test[iy])));
        }
    }
}

TEST_CASE("extrapolation_test") {
    /* define bilinear function ceofficients */
    double a_coeff = 17.837;
    double b_coeff = 25.524;
    double c_coeff = -69.228;

    /* define domain */
    double x_min = -2;
    double x_max = 3;
    double y_min = -1.5;
    double y_max = 2.5;

    /* test points */
    std::array<double, 2> point_1{-3, 0};
    std::array<double, 2> point_2{-3, -3};
    std::array<double, 2> point_3{0, -3};
    std::array<double, 2> point_4{5, -3};
    std::array<double, 2> point_5{5, 0};
    std::array<double, 2> point_6{5, 3};
    std::array<double, 2> point_7{0, 3};
    std::array<double, 2> point_8{-3, 3};

    /* evenly spaced nodes */
    int nx = 17;
    std::vector<double> x_1d(nx);
    double dx = (x_max - x_min) / (nx - 1);
    for (int i = 0; i < nx; i++) {
        x_1d[i] = x_min + i * dx;
    }

    int ny = 32;
    std::vector<double> y_1d(ny);
    double dy = (y_max - y_min) / (ny - 1);
    for (int i = 0; i < ny; i++) {
        y_1d[i] = y_min + i * dy;
    }

    auto linear_func = [&](double x, double y) { return a_coeff * x + b_coeff * y + c_coeff; };

    /* YMajor: calculate node values */
    std::vector<double> z_ymajor(nx * ny);
    for (int ix = 0; ix < nx; ix++) {
        for (int iy = 0; iy < ny; iy++) {
            z_ymajor[ix * ny + iy] = linear_func(x_1d[ix], y_1d[iy]);
        }
    }

    /* XMajor: calculate node values */
    std::vector<double> z_xmajor(nx * ny);
    for (int iy = 0; iy < ny; iy++) {
        for (int ix = 0; ix < nx; ix++) {
            z_xmajor[iy * nx + ix] = linear_func(x_1d[ix], y_1d[iy]);
        }
    }

    auto interpYMajor = JustInterp::BilinearInterpolator<double, JustInterp::YMajor>(x_1d, y_1d, z_ymajor);
    auto interpXMajor = JustInterp::BilinearInterpolator<double, JustInterp::XMajor>(x_1d, y_1d, z_xmajor);

    /* Check outside points */
    CHECK(interpYMajor(point_1[0], point_1[1]) == doctest::Approx(linear_func(point_1[0], point_1[1])));
    CHECK(interpYMajor(point_2[0], point_2[1]) == doctest::Approx(linear_func(point_2[0], point_2[1])));
    CHECK(interpYMajor(point_3[0], point_3[1]) == doctest::Approx(linear_func(point_3[0], point_3[1])));
    CHECK(interpYMajor(point_4[0], point_4[1]) == doctest::Approx(linear_func(point_4[0], point_4[1])));
    CHECK(interpYMajor(point_5[0], point_5[1]) == doctest::Approx(linear_func(point_5[0], point_5[1])));
    CHECK(interpYMajor(point_6[0], point_6[1]) == doctest::Approx(linear_func(point_6[0], point_6[1])));
    CHECK(interpYMajor(point_7[0], point_7[1]) == doctest::Approx(linear_func(point_7[0], point_7[1])));
    CHECK(interpYMajor(point_8[0], point_8[1]) == doctest::Approx(linear_func(point_8[0], point_8[1])));

    CHECK(interpXMajor(point_1[0], point_1[1]) == doctest::Approx(linear_func(point_1[0], point_1[1])));
    CHECK(interpXMajor(point_2[0], point_2[1]) == doctest::Approx(linear_func(point_2[0], point_2[1])));
    CHECK(interpXMajor(point_3[0], point_3[1]) == doctest::Approx(linear_func(point_3[0], point_3[1])));
    CHECK(interpXMajor(point_4[0], point_4[1]) == doctest::Approx(linear_func(point_4[0], point_4[1])));
    CHECK(interpXMajor(point_5[0], point_5[1]) == doctest::Approx(linear_func(point_5[0], point_5[1])));
    CHECK(interpXMajor(point_6[0], point_6[1]) == doctest::Approx(linear_func(point_6[0], point_6[1])));
    CHECK(interpXMajor(point_7[0], point_7[1]) == doctest::Approx(linear_func(point_7[0], point_7[1])));
    CHECK(interpXMajor(point_8[0], point_8[1]) == doctest::Approx(linear_func(point_8[0], point_8[1])));

    /* check x-constant-y-linear data */
    std::vector<double> x_const_lin{0};
    std::vector<double> y_const_lin{-1, 0, 1, 2, 3};
    std::vector<double> z_all_const_lin{-2, 0, 2, 4, 6};
    auto interpYMajor_const_lin = JustInterp::BilinearInterpolator<double, JustInterp::YMajor>(x_const_lin, y_const_lin, z_all_const_lin);
    auto interpXMajor_const_lin = JustInterp::BilinearInterpolator<double, JustInterp::XMajor>(x_const_lin, y_const_lin, z_all_const_lin);
    std::vector<double> y_const_lin_test_set{-2, 0.5, 2.5, 4};
    for (auto y0 : y_const_lin_test_set) {
        CHECK(interpYMajor_const_lin(-1, y0) == doctest::Approx(2 * y0));
        CHECK(interpXMajor_const_lin(-1, y0) == doctest::Approx(2 * y0));
    }

    /* check x-linear-y-const data */
    std::vector<double> x_lin_const{-7, -2, 0.5, 9, 13};
    std::vector<double> y_lin_const{17};
    std::vector<double> z_all_lin_const{-3.5, -1, 0.25, 4.5, 6.5};
    auto interpYMajor_lin_const = JustInterp::BilinearInterpolator<double, JustInterp::YMajor>(x_lin_const, y_lin_const, z_all_lin_const);
    auto interpXMajor_lin_const = JustInterp::BilinearInterpolator<double, JustInterp::XMajor>(x_lin_const, y_lin_const, z_all_lin_const);
    std::vector<double> x_lin_const_test_set{-10, -5, 10, 20};
    for (auto x0 : x_lin_const_test_set) {
        CHECK(interpYMajor_lin_const(x0, 0) == doctest::Approx(0.5 * x0));
        CHECK(interpXMajor_lin_const(x0, 0) == doctest::Approx(0.5 * x0));
    }

    /* check x-constant-y-constant data */
    std::vector<double> x_const_const{0};
    std::vector<double> y_const_const{0};
    std::vector<double> z_all_const_const{-3};
    auto interpYMajor_const_const = JustInterp::BilinearInterpolator<double, JustInterp::YMajor>(x_const_const, y_const_const, z_all_const_const);
    auto interpXMajor_const_const = JustInterp::BilinearInterpolator<double, JustInterp::XMajor>(x_const_const, y_const_const, z_all_const_const);
    CHECK(interpYMajor_const_const(-1, -1) == doctest::Approx(z_all_const_const.front()));
    CHECK(interpXMajor_const_const(-1, -1) == doctest::Approx(z_all_const_const.front()));
}

}