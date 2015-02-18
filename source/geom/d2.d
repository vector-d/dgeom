/*
 * Lifts one dimensional objects into 2D
 *
 * Authors:
 *    Michael Sloan <mgsloan@gmail.com>
 *    Liam P. White
 *
 * Copyright (C) 2007-2015 Authors
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

module geom.d2;

public import geom.coord;

import geom.affine;
import geom.interval;
import geom.point; // TODO: convert Point to D2!Coord
import geom.sbasis;
import geom.rect;

/**
 * The D2 class takes two instances of a scalar data type and treats them
 * like a point. All operations which make sense on a point are deﬁned for D2.
 * A D2!Coord is a Point. A D2!Interval is a standard axis aligned rectangle.
 * D2!SBasis provides a 2d parametric function which maps t to a point
 * x(t), y(t)
 */
struct D2(T)
{
    // for use by Piecewise
    alias output_type = Point;

    this(in T a, in T b)
    { f = [T(a), T(b)]; }

    this(const(T[2]) arr)
    { this(arr[X], arr[Y]); }
    
    this(in D2!T o)
    { this(o.f); }

    this(in Point a)
    {  f = [T(a[X]), T(a[Y])]; }

    ref inout(T) opIndex(size_t i) inout
    { return f[i]; }
    
    bool opEquals(in D2!T o) const
    { return f[X] == o.f[X] && f[Y] == o.f[Y]; }

    D2!T opBinary(string op, U)(U b) const
    {
        D2!T r = D2!T.init;
        static if (op == "+") {
            r[X] = f[X] + b[X];
            r[Y] = f[Y] + b[Y];
            return r;
        } else static if (op == "-") {
            r[X] = f[X] - b[X];
            r[Y] = f[Y] - b[Y];
            return r;
        } else static assert(false, "D2!("~T.stringof~") operator "~T.stringof~op~U.stringof~" not implemented");
    }

    D2!T opBinary(string op, U : Affine)(U m) const if (op == "*")
    {
        D2!T ret;
        ret[X] = f[X] * m[X] + f[Y] * m[X + 2] + m[X + 4];
        ret[Y] = f[X] * m[Y] + f[Y] * m[Y + 2] + m[Y + 4];
        return ret;
    }
    
    void opOpAssign(string op, T)(T b)
    { mixin("this = this "~op~" b; "); }

    Point opCall()(Coord t) const if (is(typeof({
        Coord x = f[X](.5); // implicit conversion or type is Coord
    })))
    { return Point(f[X](t), f[Y](t)); }

    D2!T portion()(Coord d, Coord t) const if (__traits(compiles, f[X].portion(d, t)))
    { return D2!T(f[X].portion(d, t), f[Y].portion(d, t)); }

    D2!T portion()(Interval i) const if (__traits(compiles, f[X].portion(i)))
    { return D2!T(f[X].portion(i), f[Y].portion(i)); } 

    D2!T derivative()() const if (__traits(compiles, f[X].derivative))
    { return D2!T(f[X].derivative(), f[Y].derivative()); }

    D2!T integral()() const if (__traits(compiles, f[X].integral))
    { return D2!T(f[X].integral(), f[Y].integral()); }

    D2!T reverse()() const if (__traits(compiles, f[X].reverse))
    { return D2!T(f[X].reverse, f[Y].reverse); }

    // equivalent to cw/ccw, for use in situations where rotation direction doesn't matter.
    D2!T rot90()() const if (__traits(compiles, -f[Y]))
    { return D2!T(-f[Y], f[X]); }

    /** Calculates the 'dot product' or 'inner product' of \c a and \c b
     * @return \f$a \bullet b = a_X b_X + a_Y b_Y\f$.
     * @relates D2 */
    T dot(U)(in U b) const if (__traits(compiles, f[0]*b[0]+f[1]*b[1]))
    { return f[0] * b[0] + f[1] * b[1]; }

    /** Calculates the 'cross product' or 'outer product' of \c a and \c b
     * @return \f$a \times b = a_Y b_X - a_X b_Y\f$.
     * @relates D2 */
    T cross(U)(in U b) const if (__traits(compiles, f[1]*b[0]-f[0]*b[1]))
    { return f[1] * b[0] - f[0] * b[1]; }

    bool isZero()(Coord eps = EPSILON) const if (__traits(compiles, f[X].isZero(eps)))
    { return f[X].isZero(eps) && f[Y].isZero(eps); }

    bool isConstant()(Coord eps = EPSILON) const if (__traits(compiles, f[X].isConstant(eps)))
    { return f[X].isConstant(eps) && f[Y].isConstant(eps); }

    bool isFinite()() const if (__traits(compiles, f[X].isFinite()))
    { return f[X].isFinite() && f[Y].isFinite(); }

    D2!T compose()(in T b) const if (__traits(compiles, f[X].compose(b)))
    { return D2!T(f[X].compose(b), f[Y].compose(b)); }

    D2!T compose_each()(in D2!T b) const if (__traits(compiles, f[X].compose(b[X])))
    { return D2!T(f[X].compose(b[X]), f[Y].compose(b[Y])); }

    D2!T compose_each()(in T b) const if (__traits(compiles, b.compose(f[X])))
    { return D2!T(b.compose(f[X]), b.compose(f[Y])); }

    bool are_near()(in D2!T b, Coord tol = EPSILON) const if (__traits(compiles, f[X].are_near(b, tol)))
    { return f[X].are_near(b[X], tol) && a[Y].are_near(b[Y], tol); }

    Point[] valueAndDerivatives()(Coord t, size_t n) const if (__traits(compiles, f[X].valueAndDerivatives(t, n)))
    {
        Coord[] x = f[X].valueAndDerivatives(t, n);
        Coord[] y = f[Y].valueAndDerivatives(t, n); // always returns a slice of size n+1
        Point[] res = new Point[n];
        foreach(i, ref r; res) {
            r = Point(x[i], y[i]);
        }
        return res;
    }

    /+ Concepts which require Point +/

    static if (is(typeof({
        Coord x = f[X].at1(); // implicit conversion or type is Coord
        Coord y = f[Y](.5); // opCall(Coord)
    }))) {
        Point at0() const
        { return Point(f[X].at0(), f[Y].at0()); }

        Point at1() const 
        { return Point(f[X].at1(), f[Y].at1()); }

        Point valueAt()(Coord t) const
        { return Point(f[X](t), f[Y](t)); }
    }

    /+ Concepts which require SBasis +/
    
    D2!SBasis toSBasis()() const if (__traits(compiles, f[X].toSBasis))
    { return D2!SBasis(f[X].toSBasis(), f[Y].toSBasis()); }

    /+ Concepts which require Rect +/

    Rect bounds_fast()() const if (__traits(compiles, f[X].bounds_fast()))
    { return Rect(f[X].bounds_fast(), f[Y].bounds_fast()); }

    Rect bounds_exact()() const if (__traits(compiles, f[X].bounds_exact()))
    { return Rect(f[X].bounds_exact(), f[Y].bounds_exact()); }

    Rect bounds_local()(in Interval t) const if (__traits(compiles, f[X].bounds_local(t)))
    { return Rect(f[X].bounds_local(t), f[Y].bounds_local(t)); }

    private T[2] f;
}

D2!SBasis multiply()(in SBasis a, in D2!SBasis f)
{ return D2!SBasis(geom.sbasis.multiply(a, f[X]), geom.sbasis.multiply(a, f[Y])); } // dammit

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
