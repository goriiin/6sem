#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <numeric> // For std::iota
#include <stdexcept> // For error handling
#include <algorithm> // For std::lower_bound

//----------------------------------------------------------------------------
// Part 1: Lagrange Interpolation
//----------------------------------------------------------------------------

// Calculates the Lagrange Polynomial value at point xp
double lagrangeInterpolation(const std::vector<double>& x_nodes, const std::vector<double>& y_nodes, double xp) {
    if (x_nodes.size() != y_nodes.size() || x_nodes.empty()) {
        throw std::invalid_argument("Invalid input nodes for Lagrange interpolation.");
    }
    size_t n = x_nodes.size();
    double interpolated_value = 0.0;

    for (size_t i = 0; i < n; ++i) {
        double basis_polynomial = 1.0;
        for (size_t j = 0; j < n; ++j) {
            if (i != j) {
                if (std::abs(x_nodes[i] - x_nodes[j]) < 1e-9) {
                     throw std::runtime_error("Duplicate x_nodes found, Lagrange interpolation requires distinct nodes.");
                }

                basis_polynomial *= (xp - x_nodes[j]) / (x_nodes[i] - x_nodes[j]);
            }
        }
        interpolated_value += y_nodes[i] * basis_polynomial;
    }
    return interpolated_value;
}

void runLagrangeExample() {
    std::cout << "============================================" << std::endl;
    std::cout << "      Lagrange Interpolation Example        " << std::endl;
    std::cout << "============================================" << std::endl;

    // Data from Variant 12, Table 1
    std::vector<double> x_nodes = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
    std::vector<double> y_nodes = {9.0, 8.2, 7.4, 6.6, 5.8, 4.9, 4.1, 3.3, 2.5, 1.8};

    std::cout << "Input Data Points:\n";
    std::cout << std::fixed << std::setprecision(1);
    std::cout << "xi\t| yi\n";
    std::cout << "--------|-------\n";
    for (size_t i = 0; i < x_nodes.size(); ++i) {
        std::cout << x_nodes[i] << "\t " << y_nodes[i] << std::endl;
    }
    std::cout << std::endl;

    // Points to evaluate the Lagrange polynomial
    std::vector<double> x_eval = {0.1, 0.2, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 8.6, 8.8};

    std::cout << "Lagrange Polynomial Evaluation:\n";
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "xp\t| P(xp)\n";
    std::cout << "--------|-------------\n";
    try {
        for (double xp : x_eval) {
            double yp = lagrangeInterpolation(x_nodes, y_nodes, xp);
            std::cout << xp << "\t " << yp << std::endl;
        }
    } catch (const std::exception& e) {
        std::cerr << "Error during Lagrange interpolation: " << e.what() << std::endl;
    }
     std::cout << std::endl;
}


//----------------------------------------------------------------------------
// Part 2: Cubic Spline Interpolation
//----------------------------------------------------------------------------

// The function to be interpolated
double func(double x) {
    // f(x) = sin(x^2) * e^(-x^2)
    // Note: exp(-x*x) can be small, sin(x*x) oscillates.
    return sin(x * x) * exp(-x * x);
}

// Structure to hold spline data
struct SplineData {
    std::vector<double> x; // Nodes xi
    std::vector<double> y; // Function values yi = f(xi)
    std::vector<double> M; // Second derivatives Mi at nodes
    double h;              // Step size (assuming uniform grid)
};

// Solves a tridiagonal system using Thomas Algorithm (TDMA)
// a: sub-diagonal (a[0] unused)
// b: main diagonal
// c: super-diagonal (c[n-1] unused)
// d: right-hand side
// x: solution vector (output)
// n: size of the system (number of unknowns)
void solveTridiagonal(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c, const std::vector<double>& d, std::vector<double>& x, int n) {
    if (n <= 0) return;

    std::vector<double> c_prime(n);
    std::vector<double> d_prime(n);

    // Forward elimination
    c_prime[0] = c[0] / b[0];
    d_prime[0] = d[0] / b[0];
    for (int i = 1; i < n; ++i) {
        double m = 1.0 / (b[i] - a[i] * c_prime[i - 1]);
        c_prime[i] = (i < n - 1) ? (c[i] * m) : 0.0; // Handle last element of c
        d_prime[i] = (d[i] - a[i] * d_prime[i - 1]) * m;
    }

    // Back substitution
    x[n - 1] = d_prime[n - 1];
    for (int i = n - 2; i >= 0; --i) {
        x[i] = d_prime[i] - c_prime[i] * x[i + 1];
    }
}


// Builds the natural cubic spline
SplineData buildNaturalCubicSpline(double a, double b, int n) {
    if (n < 3) {
        throw std::invalid_argument("Need at least 3 points for cubic spline.");
    }
    SplineData spline;
    spline.x.resize(n);
    spline.y.resize(n);
    spline.M.resize(n); // Second derivatives
    spline.h = (b - a) / (n - 1);

    // Generate nodes and function values
    for (int i = 0; i < n; ++i) {
        spline.x[i] = a + i * spline.h;
        spline.y[i] = func(spline.x[i]);
    }

    // Set up the tridiagonal system for M_1, ..., M_{n-2}
    // Natural spline boundary conditions: M_0 = 0, M_{n-1} = 0
    int system_size = n - 2; // We solve for M_1 to M_{n-2}
    if (system_size <= 0) { // Handle n=3 case separately if needed, or rely on TDMA returning empty
         spline.M[0] = 0.0;
         if (n > 1) spline.M[n - 1] = 0.0;
         if (n == 3) { // Need M1 for n=3
             double d1 = (6.0 / spline.h) * ((spline.y[2] - spline.y[1]) / spline.h - (spline.y[1] - spline.y[0]) / spline.h);
             double b1 = 4.0 * spline.h; // = 2*(h+h)
             spline.M[1] = d1 / b1;
         }
         return spline;
    }


    std::vector<double> sub_diag(system_size);   // a_i for equation i (corresponds to M_{i+1})
    std::vector<double> main_diag(system_size);  // b_i
    std::vector<double> super_diag(system_size); // c_i
    std::vector<double> rhs(system_size);        // d_i

    // Equation for M_i: h*M_{i-1} + 4h*M_i + h*M_{i+1} = 6/h * (y_{i+1} - 2y_i + y_{i-1})
    // System indices k = 0 to system_size-1 correspond to M indices i = 1 to n-2
    for (int k = 0; k < system_size; ++k) {
        int i = k + 1; // M index
        main_diag[k] = 4.0 * spline.h; // Coefficient of M_i
        rhs[k] = (6.0 / spline.h) * (spline.y[i + 1] - 2.0 * spline.y[i] + spline.y[i - 1]);

        if (k > 0) {
            sub_diag[k] = spline.h; // Coefficient of M_{i-1} = M_{k}
        }
        if (k < system_size - 1) {
            super_diag[k] = spline.h; // Coefficient of M_{i+1} = M_{k+2}
        }
    }

    // Adjust RHS for known boundary conditions M_0 = 0 and M_{n-1} = 0
    // rhs[0] -= sub_diag[0]*M_0 = 0 (no change needed as M_0=0)
    // rhs[system_size-1] -= super_diag[system_size-1]*M_{n-1} = 0 (no change needed M_{n-1}=0)
    // However, TDMA requires sub_diag[0] and super_diag[n-1] to be passed,
    // but they are effectively zero in the equations due to boundary M's
    // Adjusting indices for TDMA function:
    // a[i] coeff of x[i-1], b[i] coeff of x[i], c[i] coeff of x[i+1]
    std::vector<double> tdma_a(system_size), tdma_b(system_size), tdma_c(system_size), tdma_d(system_size);
    std::vector<double> solution_M(system_size); // Holds M_1 to M_{n-2}

    for(int k = 0; k < system_size; ++k) {
        tdma_b[k] = main_diag[k];
        tdma_d[k] = rhs[k];
        if (k > 0) tdma_a[k] = sub_diag[k]; else tdma_a[k] = 0; // a[0] is unused by TDMA func but maybe accessed
        if (k < system_size - 1) tdma_c[k] = super_diag[k]; else tdma_c[k] = 0; // c[n-1] is unused
    }


    solveTridiagonal(tdma_a, tdma_b, tdma_c, tdma_d, solution_M, system_size);

    // Store results back into spline.M
    spline.M[0] = 0.0; // Natural boundary condition
    for (int k = 0; k < system_size; ++k) {
        spline.M[k + 1] = solution_M[k];
    }
    spline.M[n - 1] = 0.0; // Natural boundary condition

    return spline;
}

// Evaluates the cubic spline S(xp) at a given point xp
double evaluateSpline(const SplineData& spline, double xp) {
    size_t n = spline.x.size();
    if (n < 2) {
         throw std::runtime_error("Spline not built or has too few points.");
    }

    // Find the interval [x_i, x_{i+1}] containing xp
    // Use lower_bound for efficiency (requires sorted x)
    auto it = std::lower_bound(spline.x.begin(), spline.x.end(), xp);
    int i = std::distance(spline.x.begin(), it);

    // Handle edge cases: xp exactly at a node or outside the range
    if (i == n) {
        i = n - 1; // Use the last interval if xp is >= x_{n-1}
    }
    if (i > 0 && (xp < spline.x[i] || std::abs(xp - spline.x[i]) < 1e-9)) {
        i--; // Use the interval [x_{i-1}, x_i]
    }
     // Now xp is in [x_i, x_{i+1}] (or at x_i)

    // Clamp i to be a valid index for the interval [x_i, x_{i+1}]
    i = std::max(0, std::min((int)n - 2, i));

    double xi = spline.x[i];
    double xi1 = spline.x[i+1];
    double yi = spline.y[i];
    double yi1 = spline.y[i+1];
    double Mi = spline.M[i];
    double Mi1 = spline.M[i+1];
    double h = spline.h; // xi1 - xi; Assumes uniform grid

    if (std::abs(h) < 1e-9) {
         throw std::runtime_error("Spline interval width is zero.");
    }

    // Cubic spline formula for S(x) in [xi, xi+1]
    double term1 = Mi * pow(xi1 - xp, 3) / (6.0 * h);
    double term2 = Mi1 * pow(xp - xi, 3) / (6.0 * h);
    double term3 = (yi / h - Mi * h / 6.0) * (xi1 - xp);
    double term4 = (yi1 / h - Mi1 * h / 6.0) * (xp - xi);

    return term1 + term2 + term3 + term4;
}


void runCubicSplineExample() {
    std::cout << "============================================" << std::endl;
    std::cout << "     Cubic Spline Interpolation Example     " << std::endl;
    std::cout << "============================================" << std::endl;

    double a = 0.0;
    double b = 3.0;
    int n = 20;
    std::cout << "Function: f(x) = sin(x^2) * exp(-x^2)" << std::endl;
    std::cout << "Interval: [" << a << ", " << b << "]" << std::endl;
    std::cout << "Number of nodes (n): " << n << std::endl << std::endl;

    SplineData spline;
    try {
        spline = buildNaturalCubicSpline(a, b, n);
    } catch (const std::exception& e) {
        std::cerr << "Error building spline: " << e.what() << std::endl;
        return;
    }

    std::cout << "Interpolation Nodes (xi, yi=f(xi)) and Second Derivatives (Mi):\n";
    std::cout << std::fixed << std::setprecision(8);
    std::cout << " xi\t\t| yi=f(xi)\t| Mi\n";
    std::cout << "--------|---------------|---------------|---------------\n";
    for (int i = 0; i < n; ++i) {
        std::cout << spline.x[i] << "\t " << spline.y[i] << "\t " << spline.M[i] << std::endl;
    }
    std::cout << std::endl;


    std::cout << "Comparison at Midpoints:\n";
    std::cout << " x_mid\t\t| S(x_mid)\t| f(x_mid)\t| |S-f|\n";
    std::cout << "---------------|---------------|---------------|---------------\n";

    for (int i = 0; i < n - 1; ++i) {
        double x_mid = (spline.x[i] + spline.x[i+1]) / 2.0;
        double spline_val = 0.0;
        double exact_val = func(x_mid);
        try {
            spline_val = evaluateSpline(spline, x_mid);
        } catch (const std::exception& e) {
             std::cerr << "Error evaluating spline at x=" << x_mid << ": " << e.what() << std::endl;
             spline_val = NAN; // Not a number
        }
        double error = std::abs(spline_val - exact_val);

        std::cout  << x_mid << "\t " << spline_val << "\t " << exact_val << "\t " << error << std::endl;
    }
     std::cout << std::endl;
}


//----------------------------------------------------------------------------
// Main Function
//----------------------------------------------------------------------------

int main() {
    // Run Lagrange Interpolation Example
    runLagrangeExample();

    std::cout << "\n\n"; // Separator

    // Run Cubic Spline Interpolation Example
    runCubicSplineExample();

    return 0;
}