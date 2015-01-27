/**
 * Linear fragment function class
 *
 *  Authors:
 *    Nathan Hurst <njh@mail.csse.monash.edu.au>
 *    Michael Sloan <mgsloan@gmail.com>
 *    Liam P. White
 *
 * Copyright (C) 2006-2015 authors
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

module geom.linear;

public import geom.coord;

import math = std.math;
import geom.interval;
import geom.rect;

/* This class feels _way_ too much like D2. */
struct Linear
{
    /* Let's have a vagueness competition! */
    this(Coord aa, Coord b) { a = [aa, b]; }
    this(Coord aa) { a = [aa, aa]; }
    
    ref inout(Coord) opIndex(size_t i) inout
    { return a[i]; }

    bool isZero(Coord eps=  EPSILON) const { return are_near(a[0], 0., eps) && are_near(a[1], 0., eps); }
    bool isConstant(Coord eps = EPSILON) const { return are_near(a[0], a[1], eps); }
    bool isFinite() const { return math.isFinite(a[0]) && math.isFinite(a[1]); }

    Coord at0() const { return a[0]; }
    Coord at1() const { return a[1]; }

    Coord valueAt(Coord t) const { return lerp(t, a[0], a[1]); }
    Coord operator()(Coord t) const { return valueAt(t); }

    Interval bounds_exact() const { return Interval(a[0], a[1]); }
    Interval bounds_fast() const { return bounds_exact(); }
    Interval bounds_local(Coord u, Coord v) const { return Interval(valueAt(u), valueAt(v)); }

    Coord tri() const 
    { return a[1] - a[0]; }

    Coord hat() const
    { return (a[1] + a[0])/2; }
    
    bool opEquals(in Linear b) const
    { return a[0] == b[0] && a[1] == b[1]; }
    
    Linear opUnary(string s)() const if (s == "-") { return Linear(-a[0], -a[1]); } /* negation */
    
    Linear opBinary(string op, T)(T rhs) const
    {
        static if (op == "+") { return Linear(a[0] + rhs[0], a[1] + rhs[1]); }
        else static if (op == "-") { return Linear(a[0] - rhs[0], a[1] - rhs[1]); }
        else static assert(false, "Linear operator "~op~" not implemented");
    }
    
    Linear opBinary(string op, T : Coord)(T rhs) const
    {
        static if (op == "+") { return Linear(a[0] + rhs, a[1] + rhs); }
        else static if (op == "-") { return Linear(a[0] - rhs, a[1] - rhs); }
        else static if (op == "*") { return Linear(a[0] * rhs, a[1] * rhs); }
        else static if (op == "/") { return Linear(a[0] / rhs, a[1] / rhs); }
        else static assert(false, "Linear operator "~op~" not implemented");
    }

    void opOpAssign(string op, T)(T rhs)
    { mixin("this = this "~op~" rhs;"); }

    private Coord[2] a = [0, 0];
}

/** Linear interpolation between two values. */
Coord lerp(Coord t, Coord a, Coord b)
{ return a*(1-t) + b*t; }

/** Gives a reversed line. */
Linear reverse(in Linear a)
{ return Linear(a[1], a[0]); }

import geom.sbasis;

/** Compute the sine of a to k terms
 * @param b linear function
 * @return sbasis sin(a)
 * It is recommended to use the piecewise version unless you have good reason.
 */
SBasis sin(Linear b, int k)
{
    SBasis s = SBasis(k+2, Linear());
    s[0] = Linear(math.sin(b[0]), math.sin(b[1]));
    double tr = s[0].tri();
    double t2 = b.tri();
    s[1] = Linear(math.cos(b[0])*t2 - tr, -math.cos(b[1])*t2 + tr);

    t2 *= t2;
    for (int i = 0; i < k; i++) {
        auto bo = Linear(4*(i+1)*s[i+1][0] - 2*s[i+1][1],
                  -2*s[i+1][0] + 4*(i+1)*s[i+1][1]);
        bo -= s[i]*(t2/(i+1));


        s[i+2] = bo / cast(Coord)i+2;
    }

    return s;
}

/** Compute the cosine of a
 * @param b linear function
 * @return sbasis cos(a)
 * It is recommended to use the piecewise version unless you have good reason.
 */
SBasis cos(Linear bo, int k)
{ return sin(Linear(bo[0] + math.PI/2, bo[1] + math.PI/2), k); }

/* Really lame unittest */
unittest
{
    Linear l = Linear(0, 0);
    Linear j = Linear(1, 1);
    
    Linear k = j + l;
    assert(k[0] == 1);
    assert(k[1] == 1);
}

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
