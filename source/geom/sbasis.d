/**
 * Defines S-power basis function class
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

module geom.sbasis;

public import geom.coord;
import math = std.math;
import geom.linear;
import geom.interval;

/**
* S-power basis function class
*
* An empty SBasis is identically 0.
*/
struct SBasis
{
    // for use by Piecewise
    alias output_type = Coord;

    /* Construct an SBasis from a single value. */
    this(Coord a) { push_back(Linear(a)); }

    /* Construct an SBasis from a linear fragment. */
    this(Coord a, Coord b) { push_back(Linear(a, b)); }

    /* Construct an SBasis from a linear fragment. */
    this(in Linear bo) { push_back(bo); }

    /* Construct from another SBasis. */
    this(in SBasis a) { d = a.d.dup; }

    /* Construct from an array of linear fragments. */
    this(const(Linear[]) ls) { d = ls.dup; }

    this(size_t n, in Linear l)
    {
        d.length = n;
        foreach (ref i; d) i = l;
    }

    /+ Get information about the SBasis +/

    ref inout(Coord) at0() inout { return d[0][0]; }
    ref inout(Coord) at1() inout { return d[0][1]; }

    int degreesOfFreedom() const { return cast(int)size() * 2; }
    
    Coord valueAt(Coord t) const
    {
        import std.range : retro;
        Coord s = t * (1 - t);
        Coord p0 = 0, p1 = 0;
        foreach (lin; retro(d)) {
            p0 = p0 * s + lin[0];
            p1 = p1 * s + lin[1];
        }
        return (1 - t)*p0 + t*p1;
    }

    /* Test to see if any of the pieces are zero. */
    bool isZero(Coord eps = EPSILON) const
    {
        if (empty()) return true;
        foreach (i; d)
            if (!i.isZero(eps)) return false;
        return true;
    }

    bool isConstant(Coord eps = EPSILON) const
    {
        if (empty()) return true;
        if (!d[0].isConstant(eps)) return false;
        foreach (i; d[1..$])
            if (!i.isZero(eps)) return false;
        return true;
    }

    bool isFinite() const
    {
        foreach (i; d)
            if(!i.isFinite()) return false;
        return true;
    }

    /** Returns a function which reverses the domain
     * useful for reversing a parameteric curve.
     */
    SBasis reverse() const
    {
        SBasis result = SBasis(this);
        foreach(k, m; d)
            result.d[k] = m.reverse();
        return result;
    }

    /** Compute the value and the first n derivatives
     * @param t position to evaluate
     * @param n number of derivatives (not counting value)
     * @return an array with the value and the n derivative evaluations
     *
     * There is an elegant way to compute the value and n derivatives for a polynomial using
     * a variant of horner's rule.  Someone will someday work out how for sbasis.
     */
    Coord[] valueAndDerivatives(Coord t, uint n) const
    {
        Coord[] ret;
        ret.length = n+1;
        ret[0] = valueAt(t);
        SBasis tmp = SBasis(this);
        for (uint i = 1; i < n + 1; ++i) {
            tmp.derive();
            ret[i] = tmp.valueAt(t);
        }
        return ret;
    }

    SBasis toSBasis() const { return SBasis(this); }

    bool opEquals(in SBasis o) const { return d == o.d; }

    /** bound the error from term truncation
     * @param tail first term to chop
     * @return the largest possible error this truncation could give
     */
    Coord tailError(uint tail) const
    {
        Interval bs = bounds_fast(this, tail);
        return math.fmax(math.fabs(bs.min()), math.fabs(bs.max()));
    }

    // compute f(g)
    SBasis opCall(in SBasis g) const
    { return compose(this, g); }

    Coord opCall(Coord t) const { return valueAt(t); }

    // remove extra zeros
    void normalize()
    {
        while (!empty() && 0 == back()[0] && 0 == back()[1])
            pop_back();
    }

    void truncate(uint k) { if (k < size()) d = d[0..k]; }
    
    SBasis opBinary(string op, T)(T b) const
    {
        static if (op == "+") {
            const uint out_size = cast(uint)math.fmax(size(), b.size());
            const uint min_size = cast(uint)math.fmin(size(), b.size());
            SBasis result;
            result.resize(out_size);

            for (uint i = 0; i < min_size; i++)
                result[i] = this[i] + b[i];
            for (uint i = min_size; i < size(); i++)
                result[i] = this[i];
            for (uint i = min_size; i < b.size(); i++)
                result[i] = b[i];

            assert(result.size() == out_size);
            return result;
        } else static if (op == "-") {
            const uint out_size = cast(uint)math.fmax(size(), b.size());
            const uint min_size = cast(uint)math.fmin(size(), b.size());
            SBasis result;
            result.resize(out_size);

            for (uint i = 0; i < min_size; i++)
                result[i] = this[i] + b[i];
            for (uint i = min_size; i < size(); i++)
                result[i] = this[i];
            for (uint i = min_size; i < b.size(); i++)
                result[i] = b[i];

            assert(result.size() == out_size);
            return result;
        } else static if (op == "*") {
            SBasis c = SBasis(size(), Linear());
            for (size_t i = 0; i < size(); i++)
                c[i] = d[i] * b;
            return c;
        } else static if (op == "/") {
            return opBinary!"*"(1./b);
        } else static assert(false, "SBasis operator "~op~" not implemented");
    }
    
    SBasis opBinary(string op, T : SBasis)(in T b) const if (op == "*")
    { return multiply(this, b); }
    
    SBasis opBinary(string op, T : Coord)(T k) const
    {
        static if (op == "+") {
            if (isZero()) return SBasis(Linear(k, k));
            SBasis a = SBasis(this);
            a[0] += k;
            return a;
        } else static if (op == "-") {
            if (isZero()) return SBasis(Linear(-k, -k));
            SBasis a = SBasis(this);
            a[0] -= k;
            return a;
        } else static if (op == "*") {
            SBasis c;
            c.resize(size());
            for (uint i = 0; i < size(); i++)
                c[i] = this[i] * k;
            return c;
        } else static if (op == "/") { return this * (1./k); }
        else static assert(false, "SBasis operator "~op~" not implemented");
    }

    SBasis opUnary(string op)() const if (op == "-")
    {
        if (isZero()) return SBasis.init;
        SBasis result = SBasis(size());
            
        foreach (i; 0 .. size()) {
            result[i] = -this[i];
        }
        return result;
    }

    void opOpAssign(string op, T)(T rhs)
    { mixin("this = this "~op~" rhs;"); }

    /+ Array-like operations +/

    ref inout(Linear) opIndex(size_t i) inout { return d[i]; }

    ref inout(Linear) back() inout { return d[$-1]; }
    ref inout(Linear) at(size_t i) inout { return d[i]; }

    bool empty() const { return d.length == 0; }
    void pop_back() { if (size() > 0) d.length--; }
    void clear() { d = []; }
    void resize(in size_t new_sz) { d.length = new_sz; }
    size_t size() const { return d.length; }

    /++++++++++++++++++++++++/

    private Linear[] d;
    private void push_back(in Linear l) { d ~= l; }

    // in place version
    private void derive()
    { 
        if (isZero()) return;
        for (uint k = 0; k < size() - 1; k++) {
            double d = (2*k+1)*(this[k][1] - this[k][0]);
            
            this[k][0] = d + (k+1)*this[k+1][0];
            this[k][1] = d - (k+1)*this[k+1][1];
        }

        int k = cast(int) size() - 1;
        double d = (2*k+1)*(this[k][1] - this[k][0]);
        if (d == 0)
            pop_back();
        else {
            this[k][0] = d;
            this[k][1] = d;
        }
    }
}

// a(b(t))

/** Compute a composed with b
 * @param a,b sbasis functions
 * @return sbasis a(b(t))
 *
 * return a0 + s(a1 + s(a2 +...  where s = (1-u)u; ak =(1 - u)a^0_k + ua^1_k
 */
SBasis compose(in SBasis a, in SBasis b)
{
    import std.range : iota, retro;	

    SBasis ctmp = SBasis(Linear(1, 1)) - b;
    SBasis s = multiply(ctmp, b);
    SBasis r;

    foreach (i; iota(0, a.size).retro) {
        r = multiply_add(r, s, 
            SBasis(Linear(a[i][0]))
             - b*a[i][0]
             + b*a[i][1]);
    }

    return r;
}

/** Compute a composed with b to k terms
 * @param a,b sbasis functions
 * @return sbasis a(b(t))
 *
 * return a0 + s(a1 + s(a2 +...  where s = (1-u)u; ak =(1 - u)a^0_k + ua^1_k
 */
SBasis compose(in SBasis a, in SBasis b, uint k)
{
    SBasis r = compose(a, b);
    r.truncate(k);
    return r;
}

SBasis portion(in SBasis t, Coord from, Coord to)
{
    Coord fv = t.valueAt(from);
    Coord tv = t.valueAt(to);
    SBasis ret = compose(t, SBasis(Linear(from, to)));
    ret.at0() = fv;
    ret.at1() = tv;
    return ret;
}

//SBasis portion(in SBasis t, Interval ivl) { return compose(t, SBasis(Linear(ivl.min(), ivl.max()))); }

/* Inversion algorithm. The notation is certainly very misleading. The
pseudocode should say:

c(v) := 0
r(u) := r_0(u) := u
for i:=0 to k do
  c_i(v) := H_0(r_i(u)/(t_1)^i; u)
  c(v) := c(v) + c_i(v)*t^i
  r(u) := r(u) ? c_i(u)*(t(u))^i
endfor
*/

/** find the function a^-1 such that a^-1 composed with a to k terms is the identity function
 * @param a sbasis function
 * @return sbasis a^-1 s.t. a^-1(a(t)) = 1
 *
 * The function must have 'unit range'("a00 = 0 and a01 = 1") and be monotonic.
 */
SBasis inverse(SBasis a, int k)
{
    debug(debug_inversion) import std.stdio;

    assert(a.size() > 0);
    Coord a0 = a[0][0];
    if (a0 != 0) {
        a -= a0;
    }
    Coord a1 = a[0][1];
    assert(a1 != 0); // not invertable.

    if(a1 != 1) {
        a /= a1;
    }
    SBasis c;
    c.resize(k);                           // c(v) := 0
    if (a.size() >= 2 && k == 2) {
        c[0] = Linear(0,1);
        Linear t1 = Linear(1+a[1][0], 1-a[1][1]);    // t_1
        c[1] = Linear(-a[1][0]/t1[0], -a[1][1]/t1[1]);
    } else if (a.size() >= 2) {                      // non linear
        SBasis r = Linear(0,1);             // r(u) := r_0(u) := u
        Linear t1 = Linear(1./(1+a[1][0]), 1./(1-a[1][1]));    // 1./t_1
        Linear one = Linear(1,1);
        Linear t1i = one;                   // t_1^0
        SBasis one_minus_a = SBasis(one) - a;
        SBasis t = multiply(one_minus_a, a); // t(u)
        SBasis ti = SBasis(one);                     // t(u)^0

        debug(debug_inversion) {
            writeln("a=", a);
            writeln("1-a=", one_minus_a);
            writeln("t1=", t1);
        }

        //c.resize(k+1, Linear(0,0));
        for (uint i = 0; i < cast(uint)k; i++) {   // for i:=0 to k do
            debug(debug_inversion) {
                writeln("-------", i, ": ---------");
                writeln("r=", r);
                writeln("c=", c);
                writeln("ti=", ti);
            }

            if(r.size() <= i)                // ensure enough space in the remainder, probably not needed
                r.resize(i+1);
            Linear ci = Linear(r[i][0]*t1i[0], r[i][1]*t1i[1]); // c_i(v) := H_0(r_i(u)/(t_1)^i; u)

            debug(debug_inversion) {
                writeln("t1i=", t1i);
                writeln("ci=", ci);
            }

            for(int dim = 0; dim < 2; dim++) // t1^-i *= 1./t1
                t1i[dim] *= t1[dim];
            c[i] = ci; // c(v) := c(v) + c_i(v)*t^i
            // change from v to u parameterisation
            SBasis civ = one_minus_a*ci[0] + a*ci[1];
            // r(u) := r(u) - c_i(u)*(t(u))^i
            // We can truncate this to the number of final terms, as no following terms can
            // contribute to the result.
            r -= multiply(civ,ti);
            r.truncate(k);
            if(r.tailError(i) == 0)
                break; // yay!
            ti = multiply(ti,t);
        }
        debug(debug_inversion) writeln("##########################");
    } else {
        c = SBasis(Linear(0,1)); // linear
    }

    c -= a0; // invert the offset
    c /= a1; // invert the slope
    return c;
}

/** Compute the pointwise product of a and b (Exact)
 * @param a,b sbasis functions
 * @returns sbasis equal to a*b
 */
SBasis multiply(in SBasis a, in SBasis b)
{
    SBasis c;
    c.resize(a.size() + b.size());
    if (a.isZero() || b.isZero())
        return c;
    return multiply_add(a, b, c);
}

/** Compute the pointwise product of a and b adding c (Exact)
 * @param a,b,c sbasis functions
 * @return sbasis equal to a*b+c
 *
 * The added term is almost free
 */
static SBasis multiply_add(in SBasis a, in SBasis b, SBasis c)
{
    if (a.isZero() || b.isZero())
        return c;
    c.resize(a.size() + b.size());

    for (uint j = 0; j < b.size(); j++) {
        for (uint i = j; i < a.size() + j; i++) {
            Coord tri = b[j].tri() * a[i-j].tri();
            c[i+1/*shift*/] += Linear(-tri);
        }
    }
    for (uint j = 0; j < b.size(); j++) {
        for (uint i = j; i < a.size() + j; i++) {
            c[i][X] += b[j][X]*a[i-j][X];
            c[i][Y] += b[j][Y]*a[i-j][Y];
        }
    }
    c.normalize();
    //assert(!(0 == c.back()[0] && 0 == c.back()[1]));
    return c;
}

Interval bounds_fast(in SBasis sb, int order = 0)
{
    Interval res = Interval(0,0); // an empty sbasis is 0.

    for (int j = cast(int)sb.size() - 1; j >= order; j--) {
        Coord a = sb[j][0];
        Coord b = sb[j][1];

        Coord v, t = 0;
        v = res[0];
        if (v < 0) t = ((b-a)/v+1)*0.5;
        if (v >= 0 || t<0 || t>1) {
            res.setMin(math.fmin(a, b));
        } else {
            res.setMin(lerp(t, a+v*t, b));
        }

        v = res[1];
        if (v > 0) t = ((b-a)/v+1)*0.5;
        if (v <= 0 || t < 0 || t > 1) {
            res.setMax(math.fmax(a,b));
        } else {
            res.setMax(lerp(t, a+v*t, b));
        }
    }

    if (order > 0) res *= math.pow(.25, order);
    return res;
}

Interval bounds_exact(in SBasis a)
{
    import geom.sbasis_roots;
    Interval result = Interval(a.at0(), a.at1());
    SBasis df = a.derivative();
    Coord[] extrema = roots(df);
    foreach (i; extrema) {
        result.expandTo(a(i));
    }
    return result;
}


/** Find a small interval that bounds a(t) for t in i to order order
 * @param sb sbasis function
 * @param i domain interval
 * @param order number of terms
 * @return interval
 */
Interval bounds_local(in SBasis sb, in Interval i, int order = 0)
{
    Coord t0=i.min(), t1=i.max(), lo=0., hi=0.;
    for (int j = cast(int)sb.size()-1; j >= order; j--) {
        Coord a=sb[j][0];
        Coord b=sb[j][1];

        Coord t = 0;
        if (lo<0) t = ((b-a)/lo+1)*0.5;
        if (lo>=0 || t<t0 || t>t1) {
            lo = math.fmin(a*(1-t0)+b*t0+lo*t0*(1-t0),a*(1-t1)+b*t1+lo*t1*(1-t1));
        } else {
            lo = lerp(t, a+lo*t, b);
        }

        if (hi>0) t = ((b-a)/hi+1)*0.5;
        if (hi<=0 || t<t0 || t>t1) {
            hi = math.fmax(a*(1-t0)+b*t0+hi*t0*(1-t0),a*(1-t1)+b*t1+hi*t1*(1-t1));
        } else {
            hi = lerp(t, a+hi*t, b);
        }
    }
    Interval res = Interval(lo,hi);
    if (order > 0) res*= math.pow(.25,order);
    return res;
}

/** Compute the integral of c (Exact)
 * @param c sbasis functions
 * @return sbasis integral(c)
 */
SBasis integral(in SBasis c)
{
    SBasis a;
    a.resize(c.size() + 1);
    a[0] = Linear(0,0);

    foreach (k; 1 .. c.size() + 1) {
        Coord ahat = -c[k-1].tri()/(2*k);
        a[k][0] = a[k][1] = ahat;
    }

    Coord aTri = 0;
    for (long k = c.size()-1; k >= 0; k--) {
        aTri = (c[k].hat() + (k+1)*aTri/2)/(2*k+1);
        a[k][0] -= aTri/2;
        a[k][1] += aTri/2;
    }
    a.normalize();
    return a;
}

SBasis derivative(in SBasis x)
{ SBasis o = SBasis(x); o.derive(); return o; }

unittest
{
    import geom.bezier;
    import geom.sbasis_roots;

    auto zero = SBasis(Bezier(0.0).toSBasis());
    auto unit = SBasis(Bezier(0.0,1.0).toSBasis());
    auto hump = SBasis(Bezier(0,1,0).toSBasis());
    auto wiggle = SBasis(Bezier(0,1,-2,3).toSBasis());

    assert(Bezier(0,0,0,0).toSBasis().isZero());
    assert(Bezier(0,1,2,3).toSBasis().isFinite());

    // note: "size" of sbasis equals half the number of coefficients
    assert(2u == Bezier(0,2,4,5).toSBasis().size());
    assert(2u == hump.size());

    assert(0.0 == wiggle.at0());
    assert(3.0 == wiggle.at1());
    assert(0.0 == wiggle.valueAt(0.5));
    assert(0.0 == wiggle(0.5));

    Coord[] vnd = wiggle.valueAndDerivatives(0.5, 5);
    assert(vnd == [0, 0, 12, 72, 0, 0]);

    assert(roots(wiggle) == [0, 0.5, 0.5]);
    
    SBasis linear_root(double t)
    { return SBasis(Linear(0-t, 1-t)); }

    SBasis array_roots(Coord[] x)
    {
        SBasis b = SBasis(1);
        for (size_t i = 0; i < x.length; i++) {
            b = geom.sbasis.multiply(b, linear_root(x[i]));
        }
        return b;
    }
    
    bool are_near(in Coord[] a, in Coord[] b, Coord eps = EPSILON)
    {
        if (a.length != b.length) return false;
        foreach(i, x; a)
            if (!geom.coord.are_near(x, b[i], eps)) return false;
        return true;
    }
    
    // The results of our rootfinding are at the moment fairly inaccurate.
    Coord eps = 5e-4;

    Coord[][] tests = [
         [0.],
         [0.5],
         [0.25, 0.75],
         [0.5, 0.5], ];
         //[0, 0.2, 0.6, 0.6, 1], ];
         //[.1,.2,.3,.4,.5,.6], // there's a 1 on the end (?)
         //[0.25,0.25,0.25,0.75,0.75,0.75] ];

    foreach (test_i; tests) {
        SBasis b = array_roots(test_i);
        assert(are_near(test_i, roots(b), eps));
    }

    SBasis reverse_wiggle = wiggle.reverse();
    assert(reverse_wiggle.at0() == wiggle.at1());
    assert(reverse_wiggle.at1() == wiggle.at0());
    assert(reverse_wiggle.valueAt(0.5) == wiggle.valueAt(0.5));
    assert(reverse_wiggle.valueAt(0.25) == wiggle.valueAt(0.75));
    assert(reverse_wiggle.reverse() == wiggle);

    for (size_t i = 1; i < 10000; ++i) {
        Coord t = math.fabs(i / 10000. - 1.);
        SBasis input = wiggle;
        SBasis[2] result;
        result[X] = geom.sbasis.portion(input, 0, t);
        result[Y] = geom.sbasis.portion(input, t, 1);
        
        // the endpoints must correspond exactly
        assert(result[X].at0() == input.at0());
        assert(result[X].at1() == result[Y].at0());
        assert(result[Y].at1() == input.at1());

        // ditto for valueAt
        assert(result[X].valueAt(0) == input.valueAt(0));
        assert(result[X].valueAt(1) == result[Y].valueAt(0));
        assert(result[Y].valueAt(1) == input.valueAt(1));
    }
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
