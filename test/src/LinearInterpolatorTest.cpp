#include <vector>
#include <cmath>

#include "doctest.h"
#include "JustInterp/JustInterp.hpp"

TEST_SUITE("LinearInterpolatorTest") {

TEST_CASE("LinearInterpolatorTest.check_node_points") {
    /* define domain */
    double x_min = -98.983;
    double x_max = 273.635;

    /* evenly spaced nodes */
    int n_nodes = 37;
    std::vector<double> x_uniform(n_nodes);
    double dx_uniform = (x_max - x_min) / (n_nodes - 1);
    for (int i = 0; i < n_nodes; i++) {
        x_uniform[i] = x_min + i * dx_uniform;
    }

    /* non-uniform nodes distribution */
    int n_test_nodes = 4 * (n_nodes - 1) + 1;
    std::vector<double> x_test(n_test_nodes);
    for (int i = 0; i < n_nodes - 1; i++) {
        x_test[4 * i    ] = x_uniform[i];
        x_test[4 * i + 1] = x_uniform[i] + 0.5   * dx_uniform;
        x_test[4 * i + 2] = x_uniform[i] + 0.75  * dx_uniform;
        x_test[4 * i + 3] = x_uniform[i] + 0.875 * dx_uniform;
    }
    x_test.back() = x_uniform.back();

    /* define test function */
    auto test_func = [](double x) { return std::cos(x) * std::log(3.765 * std::abs(x) + 1.0); };

    /* calculate node values */
    std::vector<double> y_test(n_test_nodes);
    for (int i = 0; i < n_test_nodes; i++) {
        y_test[i] = test_func(x_test[i]);
    }

    auto interpolator = JustInterp::LinearInterpolator<double>(x_test, y_test);

    for (int i = 0; i < n_test_nodes; i++) {
        CHECK(interpolator(x_test[i]) == doctest::Approx(y_test[i]));
    }
}

TEST_CASE("LinearInterpolatorTest.linear_function") {
    /* define linear function ceofficients */
    double k = 8.9876;
    double b = 3.875;

    /* define domain */
    double x_min = 15.84;
    double x_max = 87.38;

    /* define number of sample points */
    int n_known_points = 30;
    int n_test_points = 64;

    auto linear_func = [&](double x) { return k * x + b; };

    /* generate known points */
    std::vector<double> x_known(n_known_points);
    std::vector<double> y_known(n_known_points);
    for (int i = 0; i < n_known_points; i++) {
        double x = x_min + i * (x_max - x_min) / (n_known_points - 1);
        double y = linear_func(x);
        x_known[i] = x;
        y_known[i] = y;
    }

    /* generate sample points */
    std::vector<double> x_test(n_test_points);
    std::vector<double> y_test(n_test_points);
    for (int i = 0; i < n_test_points; i++) {
        double x = x_min + i * (x_max - x_min) / (n_test_points - 1);
        double y = linear_func(x);
        x_test[i] = x;
        y_test[i] = y;
    }

    /* linear interpolator */
    auto interpolator = JustInterp::LinearInterpolator<double>(x_known, y_known);

    /* check test values */
    for (int i = 0; i < n_test_points; i++) {
        CHECK(interpolator(x_test[i]) == doctest::Approx(y_test[i]));
    }
}

TEST_CASE("LinearInterpolatorTest.extrapolation_test") {
    /* define domain */
    double x_min = -2.0;
    double x_max = 1.0;

    /* sample points */
    double x1 = -2.5;
    double x2 = 2.3;

    /* evenly spaced nodes */
    int n_nodes = 27;
    std::vector<double> x_uniform(n_nodes);
    double dx_uniform = (x_max - x_min) / (n_nodes - 1);
    for (int i = 0; i < n_nodes; i++) {
        x_uniform[i] = x_min + i * dx_uniform;
    }

    /* check extrapolation for constant function */
    double const_func = 3.14;
    std::vector<double> y_const_func(n_nodes);
    for (int i = 0; i < n_nodes; i++) {
        y_const_func[i] = const_func;
    }
    auto interp_const = JustInterp::LinearInterpolator<double, JustInterp::LinearExtrapolation>(x_uniform, y_const_func);
    CHECK(interp_const(x1) == doctest::Approx(const_func));
    CHECK(interp_const(x2) == doctest::Approx(const_func));

    /* check extrapolation for linear function */
    double k = 17.6252;
    double b = 0.762524;
    auto linear_func = [&](double x) { return k * x + b; };
    std::vector<double> y_linear_func(n_nodes);
    for (int i = 0; i < n_nodes; i++) {
        y_linear_func[i] = linear_func(x_uniform[i]);
    }
    auto interp_linear = JustInterp::LinearInterpolator<double, JustInterp::LinearExtrapolation>(x_uniform, y_linear_func);
    CHECK(interp_linear(x1) == doctest::Approx(linear_func(x1)));
    CHECK(interp_linear(x2) == doctest::Approx(linear_func(x2)));

    /* check extrapolation for x^2 */
    std::vector<double> y_square_func(n_nodes);
    for (int i = 0; i < n_nodes; i++) {
        y_square_func[i] = x_uniform[i] * x_uniform[i];
    }
    auto interp_square = JustInterp::LinearInterpolator<double, JustInterp::LinearExtrapolation>(x_uniform, y_square_func);
    auto linear_left = [&](double x) {
        double m = (y_square_func[1] - y_square_func[0]) / (x_uniform[1] - x_uniform[0]);
        return y_square_func[0] + m * (x - x_uniform[0]);
    };
    auto linear_right = [&](double x) {
        double m = (y_square_func[n_nodes - 1] - y_square_func[n_nodes - 2]) / (x_uniform[n_nodes - 1] - x_uniform[n_nodes - 2]);
        return y_square_func[n_nodes - 2] + m * (x - x_uniform[n_nodes - 2]);
    };
    CHECK(interp_square(x1) == doctest::Approx(linear_left(x1)));
    CHECK(interp_square(x2) == doctest::Approx(linear_right(x2)));
}

}