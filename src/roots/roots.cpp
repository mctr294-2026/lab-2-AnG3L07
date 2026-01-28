#include "roots.hpp"
#include <cmath>

constexpr double TOL = 1e-6;
constexpr int MAX_ITERS = 1000000;

/*BISECTION*/
bool bisection(std::function<double(double)> f,
               double a, double b,
               double* root)
{
    double fa = f(a);
    double fb = f(b);

    if (fa * fb >= 0)
        return false;

    for (int i = 0; i < MAX_ITERS; ++i)
    {
        double m = 0.5 * (a + b);
        double fm = f(m);

        if (std::abs(fm) < TOL || (b - a) / 2 < TOL)
        {
            *root = m;
            return true;
        }

        if (fa * fm < 0)
        {
            b = m;
            fb = fm;
        }
        else
        {
            a = m;
            fa = fm;
        }
    }

    return false;
}

/*REGULA FALSI*/
bool regula_falsi(std::function<double(double)> f,
                  double a, double b,
                  double* root)
{
    double fa = f(a);
    double fb = f(b);

    if (fa * fb >= 0)
        return false;

    for (int i = 0; i < MAX_ITERS; ++i)
    {
        double x = (a * fb - b * fa) / (fb - fa);
        double fx = f(x);

        if (std::abs(fx) < TOL)
        {
            *root = x;
            return true;
        }

        if (fa * fx < 0)
        {
            b = x;
            fb = fx;
        }
        else
        {
            a = x;
            fa = fx;
        }
    }

    return false;
}

/*NEWTON-RAPHSON*/
bool newton_raphson(std::function<double(double)> f,
                    std::function<double(double)> g,
                    double a, double b, double c,
                    double* root)
{
    double x = c;

    for (int i = 0; i < MAX_ITERS; ++i)
    {
        double fx = f(x);
        double gx = g(x);

        if (std::abs(gx) < 1e-12)
            return false;

        double x_next = x - fx / gx;

        if (x_next < a || x_next > b)
            return false;

        if (std::abs(x_next - x) < TOL)
        {
            *root = x_next;
            return true;
        }

        x = x_next;
    }

    return false;
}

/*SECANT*/
bool secant(std::function<double(double)> f,
            double a, double b, double c,
            double* root)
{
    double x0 = c;
    double x1 = b;

    for (int i = 0; i < MAX_ITERS; ++i)
    {
        double f0 = f(x0);
        double f1 = f(x1);

        if (std::abs(f1 - f0) < 1e-12)
            return false;

        double x2 = x1 - f1 * (x1 - x0) / (f1 - f0);

        if (x2 < a || x2 > b)
            return false;

        if (std::abs(x2 - x1) < TOL)
        {
            *root = x2;
            return true;
        }

        x0 = x1;
        x1 = x2;
    }

    return false;
}
