/*
 * Copyright (C) ????-2015 Authors
 *
 * This file is part of dgeom.
 * 
 * dgeom is free software: you can redistribute 
 * it and/or modify it under the terms of the GNU General Public 
 * License as published by the Free Software Foundation, either 
 * version 3 of the License, or (at your option) any later version.
 * 
 * dgeom is distributed in the hope that it will 
 * be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with dgeom.  If not, see <http://www.gnu.org/licenses/>.
 */

module geom.solve_bezier;

import std.math;
import geom.bezier;
import geom.interval;

private const size_t MAX_DEPTH = 22;

/* Find the zeros of a Bezier.  The code subdivides until it is happy with the linearity of the
 * function.  This requires an O(degree^2) subdivision for each step, even when there is only one
 * solution.
 * 
 * We try fairly hard to correctly handle multiple roots. */

Coord[] find_bezier_roots(in Bezier b, double left_t, double right_t)
{
    Bezier bz = Bezier(b);
    Coord[] solutions;

    // a constant bezier, even if identically zero, has no roots
    if (bz.isConstant()) {
        return solutions;
    }

    while (bz[0] == 0) {
        bz = bz.deflate();
        solutions ~= 0;
    }

    if (bz.degree() == 1) {
        if (SGN(bz[0]) != SGN(bz[1])) {
            Coord d = bz[0] - bz[1];
            if (d != 0) {
                Coord r = bz[0] / d;
                if (0 <= r && r <= 1)
                    solutions ~= r;
            }
        }
        return solutions;
    }

    find_bernstein_roots(solutions, bz, 0, left_t, right_t);
    return solutions;
}

void find_bernstein_roots(ref Coord[] solutions, Bezier bz, size_t depth, Coord left_t, Coord right_t)
{
    size_t n_crossings = 0;

    int old_sign = SGN(bz[0]);
    for (size_t i = 1; i < bz.size(); ++i) {
        int sign = SGN(bz[i]);
        if (sign != 0) {
            if (sign != old_sign && old_sign != 0) {
                ++n_crossings;
            }
            old_sign = sign;
        }
    }

    if (n_crossings == 0) return;
    if (n_crossings == 1) {
        // Unique solution
        // Stop recursion when the tree is deep enough
        // if deep enough, return 1 solution at midpoint
        if (depth > MAX_DEPTH) {
            const Coord Ax = right_t - left_t;
            const Coord Ay = bz.at1() - bz.at0();

            solutions ~= left_t - Ax*bz.at0() / Ay;
            return;
        }

        Coord r = secant(bz); //XXX
        solutions ~= r*right_t + (1 - r)*left_t;
        return;
    }

    /* Otherwise, solve recursively after subdividing control polygon  */
    Bezier Left = Bezier.withOrder(bz.order()), Right = bz;
    Coord split_t = (left_t + right_t) * 0.5;

    // If subdivision is working poorly, split around the leftmost root of the derivative
    if (depth > 2) {
        Bezier dbz = bz.derivative();
        Coord[] dsolutions = dbz.roots(Interval(left_t, right_t));
        Coord dsplit_t = 0.5;

        if (dsolutions.length > 0) {
            dsplit_t = dsolutions[0];
            split_t = left_t + (right_t - left_t)*dsplit_t;
        
        }

        Bezier[2] LR = bz.subdivide(dsplit_t);
        Left = LR[X];
        Right = LR[Y];
    } else {
        // split at midpoint, because it is cheap
        Left[0] = Right[0];
        for (size_t i = 1; i < bz.size(); ++i) {
            for (size_t j = 0; j < bz.size()-i; ++j) {
                Right[j] = (Right[j] + Right[j+1]) * 0.5;
            }
            Left[i] = Right[0];
        }
    }

    Left = reverse(Left);
    while (Right.order() > 0 && fabs(Right[0]) <= 1e-10) {
        Right = Right.deflate();
        Left = Left.deflate();
        solutions ~= split_t;
    }

    Left = reverse(Left);
    if (Right.order() > 0) {
        find_bernstein_roots(solutions, Left, depth+1, left_t, split_t);
        find_bernstein_roots(solutions, Right, depth+1, split_t, right_t);
    }
}

// FIXME document me
Coord secant(Bezier bz)
{
    Coord s = 0, t = 1;
    Coord e = 1e-14;
    int side = 0;
    Coord r, fs = bz.at0(), ft = bz.at1();

    for (size_t n = 0; n < 100; ++n) {
        r = (fs*t - ft*s) / (fs - ft);
        if (fabs(t-s) < e * fabs(t+s)) {
            return r;
        }

        Coord fr = horner(bz, r);

        if (fr * ft > 0) {
            t = r; ft = fr;
            if (side == -1) fs /= 2;
            side = -1;
        }
        else if (fs * fr > 0) {
            s = r;  fs = fr;
            if (side == +1) ft /= 2;
            side = +1;
        }
        else break;
    }
    return r;
}

// suggested by Sederberg.
Coord horner(Bezier bz, double t)
{
    import geom.choose;

    Coord u, tn, tmp;
    u = 1.0 - t;
    tn = 1.0;
    tmp = bz.at0() * u;
    for (size_t i = 1; i < bz.degree(); ++i) {
        tn *= t;
        tmp = (tmp + tn * choose!Coord(bz.order(), i) * bz[i]) * u;
    }
    return (tmp + tn*t*bz.at1());
}

private int SGN(T)(T x) { return (x > 0 ? 1 : (x < 0 ? -1 : 0)); }

/*
  Local Variables:
  mode:d
  c-file-style:"stroustrup"
  c-file-offsets:((innamespace . 0)(inline-open . 0)(case-label . +))
  indent-tabs-mode:nil
  fill-column:99
  End:
*/
// vim: filetype=d:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8 :
