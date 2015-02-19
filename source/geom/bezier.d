/*
 * Bezier polynomial
 *
 * Authors:
 *   MenTaLguY <mental@rydia.net>
 *   Michael Sloan <mgsloan@gmail.com>
 *   Nathan Hurst <njh@njhurst.com>
 *   Liam P. White
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

module geom.bezier;

import geom.coord;
import geom.interval;
import geom.point;

struct Bezier
{
    // for use by Piecewise
    alias output_type = Coord;

    // Copy constructor
    this(in Bezier b) { c_ = b.c_.dup; }

    // Construct an n-order bezier
    this(A...)(A a) { foreach (t; a) c_ ~= t; }
    this(Coord[] a...) { c_ = a.dup; }

    static Bezier withOrder(size_t order)
    { Bezier b; b.resize(order + 1); return b; }

    bool isZero(double eps = EPSILON) const
    {
        foreach (c; c_)
            if (!are_near(c, 0., eps)) return false;
        return true;
    }
    
    bool isConstant(double eps = EPSILON) const
    {
        foreach (c; c_[1 .. $])
            if (!are_near(c, c_[0], eps)) return false;
        return true;
    }
    
    bool isFinite() const
    {
        import std.math : isInfinity;
        foreach (c; c_)
            if (isInfinity(c)) return false;
        return true;
    }

    /** The size of the returned slice equals n_derivs+1. */
    Coord[] valueAndDerivatives(Coord t, size_t n_derivs) const
    {
        /* This is inelegant, as it uses several extra stores.  I think there might be a way to
         * evaluate roughly in situ. */

         // initialize return array with zeroes, such that we only need to replace the non-zero derivs
        Coord[] val_n_der;
        val_n_der.length = n_derivs + 1;
        foreach (ref i; val_n_der) i = 0;

        // initialize temp storage variables
        Coord[] d_ = c_.dup; 
        d_.length = order() + 1;

        size_t nn = n_derivs + 1;
        if (n_derivs > order()) {
            nn = order() + 1; // only calculate the non zero derivs
        }

        for (size_t di = 0; di < nn; di++) {
            val_n_der[di] = bernsteinValueAt(t, d_, order() - di);
            for (size_t i = 0; i < order() - di; i++) {
                d_[i] = (order()-di)*(d_[i+1] - d_[i]);
            }
        }

        return val_n_der;
    }

    Bezier[2] subdivide(Coord t) const
    {
        Coord[][2] s = casteljau_subdivision(t, c_);
        Bezier left = from_array(s[X]);
        Bezier right = from_array(s[Y]);
        return [left, right];
    }

    Coord[] roots() const
    {
        import geom.solve_bezier;
        Coord[] solutions = find_bezier_roots(this, 0, 1);
        return solutions;
    }

    Coord[] roots(in Interval ivl) const
    {
        import geom.solve_bezier;
        Coord[] solutions;
        find_bernstein_roots(solutions, Bezier(this), 0, ivl.min(), ivl.max());
        return solutions;
    }

    Bezier forward_difference(uint k)
    {
        import geom.choose;
        Bezier fd = withOrder(order() - k);
        size_t n = fd.size();
        
        for (uint i = 0; i < n; i++) {
            fd[i] = 0;
            for (uint j = i; j < n; j++) {
                fd[i] += (((j)&1)?-c_[j]:c_[j]) * choose!Coord(n, j-i);
            }
        }
        return fd;
    }
    
    Coord valueAt(Coord t) const
    { return bernsteinValueAt(t, c_, order()); }

    Bezier elevate_degree() const
    {
        Bezier ed = withOrder(order() + 1);
        size_t n = size();
        ed[0] = c_[0];
        ed[n] = c_[n-1];
        for (uint i = 1; i < n; i++) {
            ed[i] = (i*c_[i-1] + (n - i)*c_[i])/(n);
        }
        return ed;
    }

    Bezier elevate_to_degree(size_t newDegree) const
    {
        Bezier ed = Bezier(this);
        foreach (i; degree() .. newDegree)
            ed = ed.elevate_degree();
        return ed;
    }

    Bezier reduce_degree() const
    {
        if (order() == 0) return Bezier(this);
        Bezier ed = withOrder(order() - 1);
        size_t n = size();
        ed[0] = c_[0];
        ed[n-1] = c_[n]; // ensure exact endpoints
        size_t middle = n / 2;

        for (size_t i = 1; i < middle; ++i) {
            ed[i] = (n*c_[i] - i*ed[i-1])/(n-i);
        }

        for (size_t i = n - 1; i >= middle; --i) {
            ed[i] = (n*c_[i] - i*ed[n-i])/(i);
        }

        return ed;
    }

    Bezier deflate() const
    {
        if (order() == 0) return Bezier(this);
        size_t n = order();
        Bezier b = withOrder(n - 1);
        for (size_t i = 0; i < n; i++) {
            b[i] = (n * c_[i+1])/(i+1); // not confusing at all!
        }
        return b;
    }

    Coord opCall(Coord t) const { return valueAt(t); }
    
    Bezier opBinary(string op, T)(T b) const
    {
        static if (op == "+") {
            Bezier result = Bezier(this);
            foreach (ref t; result.c_) t += b;
            return result;
        } else static if (op == "-") {
            Bezier result = Bezier(this);
            foreach (ref t; result.c_) t += b;
            return result;
        }
    }
    
    // These are the only direct mutators
    ref inout(Coord) opIndex(size_t ix) inout { return c_[ix]; }
    void setCoeff(size_t ix, Coord val) { c_[ix] = val; }

    Coord at0() const { return c_[0]; }
    Coord at1() const { return c_[$-1]; }

    size_t degree() const { return order(); }
    size_t order() const { return c_.length - 1; }
    size_t size() const { return c_.length; }
    
    void clear() { c_ = []; }
    void resize(size_t n, Coord v = 0)
    {
        Coord[] x; x.length = n;
        foreach (ref t; x) t = v;
        c_ = x;
    }

    protected static Bezier from_array(in Coord[] c) { Bezier b; b.c_ = c.dup; return b; }

    /** coefficients of the polynomial */
    private Coord[] c_ = [];
}

/** Perform De Casteljau subdivision of a Bezier polynomial.
 * Given an array of coefficients and a time value, computes two new Bernstein-Bezier basis
 * polynomials corresponding to the \f$[0, t]\f$ and \f$[t, 1]\f$ intervals of the original one.
 * @param t Time value
 * @param v Array of input coordinates
 *    The order of the input polynomial, equal to one less the number of coefficients, is
 *    inferred from the length of the input array.
 * @return Pair of polynomials
 *    [0] Polynomial corresponding to \f$[0, t]\f$
 *    [1] Polynomial corresponding to \f$[t, 1]\f$
 */
private T[][2] casteljau_subdivision(T)(Coord t, in T[] v)
{
    import geom.linear : lerp;

    size_t order = v.length - 1;
    T[] left = new T[order + 1];
    T[] right = v.dup;

    for (size_t i = 1; i <= order; ++i) {
        left[i-1] = right[0];
        for (size_t j = i; j > 0; --j) {
            right[j-1] = lerp(t, right[j-1], right[j]);
        }
    }

    // The Horner-like scheme gives very slightly different results, but we need
    // the result of subdivision to match exactly with Bezier's valueAt function.
    T val = bernsteinValueAt(t, v, order);
    left[order] = val;
    right[0] = val;

    return [left, right];
}

private T bernsteinValueAt(T)(Coord t, in T[] c_, size_t n)
{
    Coord u = 1.0 - t;
    Coord bc = 1;
    Coord tn = 1;
    T tmp = c_[0] * u;
    for (size_t i = 1; i < n; i++) {
        tn = tn*t;
        bc = bc*(n-i+1)/i;
        tmp = (tmp + tn*bc*c_[i])*u;
    }
    return (tmp + tn*t*c_[n]);
}

Bezier multiply(in Bezier f, in Bezier g)
{
    import geom.choose;
    size_t m = f.order();
    size_t n = g.order();
    Bezier h = Bezier.withOrder(m + n);
    // h_k = sum_(i+j=k) (m i)f_i (n j)g_j / (m+n k)
    
    for (size_t i = 0; i <= m; i++) {
        const Coord fi = choose!Coord(m,i)*f[i];
        for (size_t j = 0; j <= n; j++) {
            h[i+j] += fi * choose!Coord(n,j)*g[j];
        }
    }
    for (size_t k = 0; k <= m+n; k++) {
        h[k] /= choose!Coord(m+n, k);
    }
    return h;
}

Bezier reverse(in Bezier a)
{
    Bezier result = Bezier.withOrder(a.order());
    for (size_t i = 0; i <= a.order(); i++)
        result[i] = a[a.order() - i];
    return result;
}

Bezier portion(in Bezier a, Coord from, Coord to)
{
    import std.algorithm : swap, reverse;

    Bezier ret = Bezier(a);

    bool reverse_result = false;
    if (from > to) {
        swap(from, to);
        reverse_result = true;
    }

    Coord[] empty_slice = [];

    do {
        if (from == 0) {
            if (to == 1) {
                break;
            }
            ret.c_ = casteljau_subdivision(to, ret.c_)[X];
            break;
        }
        ret.c_ = casteljau_subdivision(from, ret.c_)[Y];
        if (to == 1) break;
        ret.c_ = casteljau_subdivision((to - from) / (1 - from), ret.c_)[X];
        // to protect against numerical inaccuracy in the above expression, we manually set
        // the last coefficient to a value evaluated directly from the original polynomial
        ret.c_[ret.order()] = a.valueAt(to);
    } while(0);

    if (reverse_result) {
        reverse(ret.c_);
    }
    return ret;
}

import geom.d2;

// Todo: how to handle differing orders
Point[] bezier_points(in D2!Bezier a)
{
    Point[] result;
    for (size_t i = 0; i <= a[0].order(); i++) {
        Point p = [ a[0][i], a[1][i] ];
        result ~= p;
    }
    return result;
}

Bezier derivative(in Bezier a)
{
    //if (a.order() == 1) return Bezier(0.0);
    if (a.order() == 1) return Bezier(a.c_[1]-a.c_[0]);
    Bezier der = Bezier.withOrder(a.order()-1);

    for (size_t i = 0; i < a.order(); i++) {
        der.c_[i] = a.order()*(a.c_[i+1] - a.c_[i]);
    }
    return der;
}

Bezier integral(in Bezier a)
{
    Bezier inte = Bezier.withOrder(a.order() + 1);

    inte[0] = 0;
    for (size_t i = 0; i < inte.order(); i++) {
        inte[i+1] = inte[i] + a[i]/(inte.order());
    }
    return inte;
}

import geom.sbasis;

SBasis toSBasis(in Bezier b)
{
    import geom.sbasis_to_bezier;
    SBasis sb;
    bezier_to_sbasis(sb, b);
    return sb;
}

import geom.interval;

Interval bounds_fast(in Bezier b)
{ return Interval.from_array(b.c_); }

//TODO: better bounds exact
Interval bounds_exact(in Bezier b)
{
    return geom.sbasis.bounds_exact(b.toSBasis());
}

Interval bounds_local(in Bezier b, Interval i)
{
    //return bounds_local(b.toSBasis(), i);
    return ( i.isEmpty() ? Interval.empty() : bounds_fast(portion(b, i.min(), i.max())) );
}

unittest
{
    bool are_equal(in Bezier A, in Bezier B)
    {
        const maxSize = A.size() > B.size() ? A.size() : B.size();
        Coord t = 0., dt = 1./maxSize;
        foreach (i; 0 .. maxSize) {
            assert(are_near(A.valueAt(t), B.valueAt(t)));
            t += dt;
        }
        return true;
    }

    const zero = Bezier(0,0);
    const unit = Bezier(0,1);
    const hump = Bezier(0,1,0);
    const wiggle = Bezier(0,1,-2,3);

    const(Bezier[4]) fragments = [zero, unit, hump, wiggle];

    /+ Basics +/
    assert(Bezier(0,0,0,0).isZero());
    assert(Bezier(0,1,2,3).isFinite());
    assert(3u == Bezier(0,2,4,5).order());
    assert(2u == hump.degree());
    assert(3u == hump.size());

    /+ valueAt +/
    assert(0.0 == wiggle.at0());
    assert(3.0 == wiggle.at1());
    assert(0.0 == wiggle.valueAt(0.5));
    assert(0.0 == wiggle(0.5));

    /+ Mutation +/
    Bezier bigun = Bezier.withOrder(30);
    bigun.setCoeff(5,10.0);
    foreach (i; 0 .. bigun.size()) {
        assert(((i == 5) ? 10 : 0) == bigun[i]);
    }

    bigun[5] = -3;
    foreach (i; 0 .. bigun.size()) {
        assert(((i == 5) ? -3 : 0) == bigun[i]);
    }

    /+ multiderivative +/
    Coord[] vnd = wiggle.valueAndDerivatives(0.5, 5);
    assert(vnd == [0,0,12,72,0,0]);

    /+ degree elevation +/
    assert(are_equal(wiggle, wiggle));
    Bezier Q = Bezier(wiggle);
    Bezier P = Q.elevate_degree();
    assert(P.size() == Q.size()+1);
    assert(are_equal(Q, P));
    Q = Bezier(wiggle);
    P = Q.elevate_to_degree(10);
    assert(10u == P.order());
    assert(are_equal(Q, P));

    Bezier linear_root(Coord t) { return Bezier(0-t, 1-t); }
    Bezier array_roots(Coord[] x)
    {
        auto b = Bezier(1);
        foreach (i; x) {
            b = multiply(b, linear_root(i));
        }
        return b;
    }

    /+ deflation +/
    Bezier b = array_roots([0,0.25,0.5]);
    assert(are_near(0, b.at0()));
    b = b.deflate();
    assert(are_near(0, b.valueAt(0.25)));
    b = b.subdivide(0.25)[Y];
    assert(are_near(0, b.at0()));
    b = b.deflate();
    const rootposition = (0.5-0.25) / (1-0.25);
    assert(are_near(0, b.valueAt(rootposition)));
    b = b.subdivide(rootposition)[Y];
    assert(are_near(0, b.at0()));

    bool array_eq(in Coord[] a, in Coord[] b)
    {
        if (a.length != b.length) return false;
        foreach (i; 0 .. a.length)
            if (!are_near(a[i], b[i])) return false;
        return true;
    }

    /+ roots +/
    Coord[][] tests = [
        [0.],
        [0.,0],
        [.5],
        [.5,.5],
        [.1,.1],
        // [.1,.1,.1], // our results are close but not close enough (1e-2)
        [.25,.75],
        [.5,.5],
        // [0,.2,.6,.6,1],
        // [.1,.2,.3,.4,.5,.6],
        [.25,.25,.25,.75,.75,.75],
    ];

    foreach (i, test; tests) {
        b = array_roots(test);
        assert(array_eq(test, b.roots));
    }

    /+ operators +/
    Bezier reverse_wiggle = reverse(wiggle);
    assert(reverse_wiggle.at0() == wiggle.at1());
    assert(reverse_wiggle.at1() == wiggle.at0());
    assert(are_equal(reverse(reverse_wiggle), wiggle));
    assert(are_equal(derivative(integral(wiggle)), wiggle));
    assert(array_eq([0.5], derivative(hump).roots));

    /+ bounds +/
    assert(bounds_fast(hump).contains(Interval(0,hump.valueAt(0.5))));
    assert(Interval(0,hump.valueAt(0.5)) == bounds_exact(hump));

    T min(T)(T o, T t) { return o < t ? o : t; } // irritating

    auto tight_local_bounds = Interval(min(hump.valueAt(0.3),hump.valueAt(0.6)), hump.valueAt(0.5));
    assert(bounds_local(hump, Interval(0.3, 0.6)).contains(tight_local_bounds));
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
