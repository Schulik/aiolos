#ifndef _BRENT_H_
#define _BRENT_H_

#include <cmath>
#include <limits>
#include <stdexcept>

/* class Brent
 *   Solve a 1D root finding problem. Follows Numerical Recipes.
 */
class Brent {
   public:
    Brent(){};
    Brent(double tol) : _tol(tol){};
    Brent(double tol, int max_iter=128) : _tol(tol), _max_iter(max_iter){};

    template <class System>
    double solve(double a, double b, System& sys) {
        using std::abs;
        using std::min;
        double fa = sys(a);
        double fb = sys(b);

        if (not(fa * fb <= 0))
            throw std::invalid_argument("Root must be bracketed");

        double EPS = std::numeric_limits<double>::epsilon();
        double c = b, fc = fb;
        double d=0, e=0;
        for (int i = 0; i < _max_iter; i++) {
            // Order interval
            if (fb * fc > 0) {
                c = a;
                fc = fa;

                e = d = b - a;
            }
            if (abs(fc) < abs(fb)) {
                a = b;
                b = c;
                c = a;
                fa = fb;
                fb = fc;
                fc = fa;
            }

            double tol1 = 2 * EPS * abs(b) + 0.5 * _tol;
            double xm = 0.5 * (c - b);
            if (abs(xm) <= tol1) return b;

            // Interpolate to new root
            if (abs(e) >= tol1 && abs(fa) > abs(fb)) {
                double s = fb / fa;
                double p, q;
                if (a == c) {  // Use linear
                    p = 2 * xm * s;
                    q = 1 - s;
                } else {  // Use inverse quadratic
                    double r;
                    q = fa / fc;
                    r = fb / fc;

                    p = s * (2 * xm * q * (q - r) - (b - a) * (r - 1));
                    q = (q - 1) * (r - 1) * (s - 1);
                }
                if (p > 0)
                    q = -q;
                else
                    p = -p;

                // Check whether the interpolation succeeded
                if (2 * p < min(3 * xm * q - abs(tol1 * q), abs(e * q))) {
                    e = d;  // Accept interpolation
                    d = p / q;
                } else {
                    d = xm;  // Use bisection
                    e = d;
                }
            }  // Slow convergence
            else {
                d = xm;
                e = d;
            }
            // Store previous best guess
            a = b;
            fa = fb;
            // Evaluate new root
            if (abs(d) > tol1)
                b += d;
            else {
                if (xm > 0)
                    b += tol1;
                else
                    b -= tol1;
            }
            fb = sys(b);
        }
        throw std::runtime_error("Brent's method did not converge");
    }

   private:
    double _tol = 1e-12;
    int _max_iter = 128;
};

#endif  //_BRENT_H_