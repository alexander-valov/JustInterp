#include <vector>
#include <array>
#include <cmath>

#include "gtest/gtest.h"
#include "BilinearInterpolator.hpp"

TEST(BiinearInterpolatorTest, check_node_points) {
    const double TOL = 1.0e-14;

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
            EXPECT_NEAR(test_func(x_1d[ix], y_1d[iy]), interpYMajor(x_1d[ix], y_1d[iy]), TOL);
        }
    }

    /* check XMajor nodes values */
    for (int ix = 0; ix < nx; ix++) {
        for (int iy = 0; iy < ny; iy++) {
            EXPECT_NEAR(test_func(x_1d[ix], y_1d[iy]), interpXMajor(x_1d[ix], y_1d[iy]), TOL);
        }
    }

    /* YMajor: check GridInterpolation */
    auto result_ymajor = interpYMajor.GridInterpolation(x_1d, y_1d);
    for (int i = 0; i < nx * ny; i++) {
        EXPECT_NEAR(z_ymajor[i], result_ymajor[i], TOL);
    }

    /* XMajor: check GridInterpolation */
    auto result_xmajor = interpXMajor.GridInterpolation(x_1d, y_1d);
    for (int i = 0; i < nx * ny; i++) {
        EXPECT_NEAR(z_xmajor[i], result_xmajor[i], TOL);
    }
}

TEST(BiinearInterpolatorTest, bilinear_function) {
    const double TOL = 1.0e-12;

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
            EXPECT_NEAR(linear_func(x_1d_test[ix], y_1d_test[iy]), interpYMajor(x_1d_test[ix], y_1d_test[iy]), TOL);
        }
    }

    /* check XMajor test values */
    for (int ix = 0; ix < nx_test; ix++) {
        for (int iy = 0; iy < ny_test; iy++) {
            EXPECT_NEAR(linear_func(x_1d_test[ix], y_1d_test[iy]), interpXMajor(x_1d_test[ix], y_1d_test[iy]), TOL);
        }
    }
}

TEST(BiinearInterpolatorTest, extrapolation_test) {
    const double TOL = 1.0e-12;

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
    EXPECT_NEAR(linear_func(point_1[0], point_1[1]), interpYMajor(point_1[0], point_1[1]), TOL);
    EXPECT_NEAR(linear_func(point_2[0], point_2[1]), interpYMajor(point_2[0], point_2[1]), TOL);
    EXPECT_NEAR(linear_func(point_3[0], point_3[1]), interpYMajor(point_3[0], point_3[1]), TOL);
    EXPECT_NEAR(linear_func(point_4[0], point_4[1]), interpYMajor(point_4[0], point_4[1]), TOL);
    EXPECT_NEAR(linear_func(point_5[0], point_5[1]), interpYMajor(point_5[0], point_5[1]), TOL);
    EXPECT_NEAR(linear_func(point_6[0], point_6[1]), interpYMajor(point_6[0], point_6[1]), TOL);
    EXPECT_NEAR(linear_func(point_7[0], point_7[1]), interpYMajor(point_7[0], point_7[1]), TOL);
    EXPECT_NEAR(linear_func(point_8[0], point_8[1]), interpYMajor(point_8[0], point_8[1]), TOL);

    EXPECT_NEAR(linear_func(point_1[0], point_1[1]), interpXMajor(point_1[0], point_1[1]), TOL);
    EXPECT_NEAR(linear_func(point_2[0], point_2[1]), interpXMajor(point_2[0], point_2[1]), TOL);
    EXPECT_NEAR(linear_func(point_3[0], point_3[1]), interpXMajor(point_3[0], point_3[1]), TOL);
    EXPECT_NEAR(linear_func(point_4[0], point_4[1]), interpXMajor(point_4[0], point_4[1]), TOL);
    EXPECT_NEAR(linear_func(point_5[0], point_5[1]), interpXMajor(point_5[0], point_5[1]), TOL);
    EXPECT_NEAR(linear_func(point_6[0], point_6[1]), interpXMajor(point_6[0], point_6[1]), TOL);
    EXPECT_NEAR(linear_func(point_7[0], point_7[1]), interpXMajor(point_7[0], point_7[1]), TOL);
    EXPECT_NEAR(linear_func(point_8[0], point_8[1]), interpXMajor(point_8[0], point_8[1]), TOL);

    /* check x-constant-y-linear data */
    std::vector<double> x_const_lin{0};
    std::vector<double> y_const_lin{-1, 0, 1, 2, 3};
    std::vector<double> z_all_const_lin{-2, 0, 2, 4, 6};
    auto interpYMajor_const_lin = JustInterp::BilinearInterpolator<double, JustInterp::YMajor>(x_const_lin, y_const_lin, z_all_const_lin);
    auto interpXMajor_const_lin = JustInterp::BilinearInterpolator<double, JustInterp::XMajor>(x_const_lin, y_const_lin, z_all_const_lin);
    std::vector<double> y_const_lin_test_set{-2, 0.5, 2.5, 4};
    for (auto y0 : y_const_lin_test_set) {
        EXPECT_NEAR(2 * y0, interpYMajor_const_lin(-1, y0), TOL);
        EXPECT_NEAR(2 * y0, interpXMajor_const_lin(-1, y0), TOL);
    }

    /* check x-linear-y-const data */
    std::vector<double> x_lin_const{-7, -2, 0.5, 9, 13};
    std::vector<double> y_lin_const{17};
    std::vector<double> z_all_lin_const{-3.5, -1, 0.25, 4.5, 6.5};
    auto interpYMajor_lin_const = JustInterp::BilinearInterpolator<double, JustInterp::YMajor>(x_lin_const, y_lin_const, z_all_lin_const);
    auto interpXMajor_lin_const = JustInterp::BilinearInterpolator<double, JustInterp::XMajor>(x_lin_const, y_lin_const, z_all_lin_const);
    std::vector<double> x_lin_const_test_set{-10, -5, 10, 20};
    for (auto x0 : x_lin_const_test_set) {
        EXPECT_NEAR(0.5 * x0, interpYMajor_lin_const(x0, 0), TOL);
        EXPECT_NEAR(0.5 * x0, interpXMajor_lin_const(x0, 0), TOL);
    }

    /* check x-constant-y-constant data */
    std::vector<double> x_const_const{0};
    std::vector<double> y_const_const{0};
    std::vector<double> z_all_const_const{-3};
    auto interpYMajor_const_const = JustInterp::BilinearInterpolator<double, JustInterp::YMajor>(x_const_const, y_const_const, z_all_const_const);
    auto interpXMajor_const_const = JustInterp::BilinearInterpolator<double, JustInterp::XMajor>(x_const_const, y_const_const, z_all_const_const);
    EXPECT_NEAR(z_all_const_const.front(), interpYMajor_const_const(-1, -1), TOL);
    EXPECT_NEAR(z_all_const_const.front(), interpXMajor_const_const(-1, -1), TOL);
}