/**
 * Bezier curve
 *
 * Authors:
 *   MenTaLguY <mental@rydia.net>
 *   Marco Cecchetti <mrcekets at gmail.com>
 *   Krzysztof Kosiński <tweenk.pl@gmail.com>
 *   Liam P. White
 * 
 * Copyright 2007-2015 Authors
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

module geom.bezier_curve;

public import geom.coord;

import geom.affine;
import geom.bezier;
import geom.curve;
import geom.d2;
import geom.interval;
import geom.point;
import geom.rect;
import geom.sbasis;
import geom.transforms;

class BezierCurve : Curve
{
    /+ Access and modify control points +/

    /** Get the order of the Bezier curve.
     * A Bezier curve has order() + 1 control points. */
    size_t order() const { return inner[X].order(); }

    /** Get the control points.
     * @return Vector with order() + 1 control points. */
    Point[] points() const { return bezier_points(inner); }

    /** Modify a control point.
     * @param ix The zero-based index of the point to modify. Note that the caller is responsible for checking that this value is <= order().
     * @param v The new value of the point */
    void setPoint(size_t ix, Point v)
    {
        inner[X].setPoint(ix, v[X]);
        inner[Y].setPoint(ix, v[Y]);
    }
    /** Set new control points.
     * @param ps Vector which must contain order() + 1 points.
     *           Note that the caller is responsible for checking the size of this vector.
     * @throws LogicalError Thrown when the size of the vector does not match the order. */
    void setPoints(in Point[] ps)
    {
        // HLineSegment will need to redefine this
        if (ps.length != order() + 1)
            throw new Exception("BezierCurve::setPoints: incorrect number of points in vector");
        for (size_t i = 0; i <= order(); i++) {
            setPoint(i, ps[i]);
        }
    }
    /** Access control points of the curve.
     * @param ix The (zero-based) index of the control point. Note that the caller is responsible for checking that this value is <= order().
     * @return The control point. No-reference return, use setPoint() to modify control points. */
    Point opIndex(size_t ix) const { return Point(inner[X][ix], inner[Y][ix]); }

    /+ Curve interface, documented in curve.d +/

    Point initialPoint() const { return inner.at0(); }
    Point finalPoint() const { return inner.at1(); }
    bool isDegenerate() const { return inner.isConstant(); }
    void setInitial(in Point v) { setPoint(0, v); }
    void setFinal(in Point v) { setPoint(order(), v); }
    Rect boundsFast() const { return Rect(inner[X].bounds_fast(), inner[Y].bounds_fast()); }
    Rect boundsExact() const { return Rect(inner[X].bounds_exact(), inner[Y].bounds_exact()); }
    Rect boundsLocal(in Interval i, uint deg) const
    {
        if (i.isEmpty()) return Rect.empty();
        if (i.min() == 0 && i.max() == 1) return boundsFast();
        if (deg == 0) return Rect(inner[X].bounds_local(i), inner[Y].bounds_local(i));
        // TODO: UUUUUUGGGLLY
        if (deg == 1 && order() > 1) return Rect(bounds_local(inner[X].derivative(), i),
                                                 bounds_local(inner[Y].derivative(), i));
        return Rect.empty();
    }

    BezierCurve duplicate() const
    { return new BezierCurve(this); }

    // FIXME overload resolution in D2 sucks, fix it one day

    BezierCurve portion(Coord f, Coord t) const
    { return new BezierCurve(inner[X].portion(f, t), inner[Y].portion(f, t)); }

    BezierCurve reverse() const
    { return new BezierCurve(inner[X].reverse(), inner[Y].reverse()); }

    void transform(in Affine m)
    {
        Point[] ps = points();
        foreach (ref i; ps) {
            i = i * m;
        }
        setPoints(ps);
    }

    BezierCurve opBinary(string op)(in Translate m) const if (op == "*")
    {
        BezierCurve ret = new BezierCurve(this);
        ret.inner += m.vector();
        return ret;
    }

    BezierCurve derivative() const
    { return new BezierCurve(inner[X].derivative(), inner[Y].derivative()); }

    int degreesOfFreedom() const
    { return cast(int)(2 * (order() + 1)); }

    Coord[] roots(Coord v, size_t d) const { return (inner[d] - v).roots();  }

    Coord length(Coord tolerance = EPSILON) const
    {
        switch (order()) {
        case 0:
            return 0;
        case 1:
            return distance(initialPoint(), finalPoint());
        default:
            return bezier_length(points(), tolerance);
        }
    }

    Point pointAt(Coord t) const { return inner.valueAt(t); }
    Point[] pointAndDerivatives(Coord t, uint n) const { return inner.valueAndDerivatives(t, n); }
    Coord valueAt(Coord t, size_t d) const { return inner[d].valueAt(t); }
    D2!SBasis toSBasis() const { return D2!SBasis(inner[X].toSBasis(), inner[Y].toSBasis()); }

protected:
    D2!Bezier inner;
    this() { Bezier b; inner = [b, b]; }
    this(in BezierCurve o) { inner = D2!Bezier(o.inner); }
    this(in D2!Bezier b) { inner = D2!Bezier(b); }
    this(in Bezier x, in Bezier y) { inner = D2!Bezier(Bezier(x), Bezier(y)); }
    this(in Point[] pts);
}

class BezierCurveN(size_t N) : BezierCurve
{
    /+ Construct Bezier curves +/

    /** Construct a Bezier curve of the specified order with all points zero. */
    this() { super(Bezier.withOrder(degree), Bezier.withOrder(degree)); }

    /** Construct from 2D Bezier polynomial. */
    this(in D2!Bezier x) { super(x); }

    /** Construct from two 1D Bezier polynomials of the same order. */
    this(Bezier x, Bezier y) { super(x, y); }

    /** Construct a Bezier curve from a vector of its control points. */
    this(in Point[] points...)
    {
        size_t ord = points.length - 1;
        if (ord != degree) throw new Exception("BezierCurve!N does not match number of points");
        for (size_t d = 0; d < 2; ++d) {
            inner[d] = Bezier.withOrder(ord);
            for (size_t i = 0; i <= ord; i++)
                inner[d][i] = points[i][d];
        }
    }
    
    /** Construct an n-order Bezier curve with compiler-determined degree */
    static auto fromPoints(A...)(A a)
    { return BezierCurveN!(a.length)(a); }

    /** Divide a Bezier curve into two curves
     * @param t Time value
     * @return Pair of Bezier curves \f$(\mathbf{D}, \mathbf{E})\f$ such that
     *         \f$\mathbf{D}[ [0,1] ] = \mathbf{C}[ [0,t] ]\f$ and
     *         \f$\mathbf{E}[ [0,1] ] = \mathbf{C}[ [t,1] ]\f$ */
    BezierCurveN!N[2] subdivide(Coord t) const
    {
        Bezier[2] sx = inner[X].subdivide(t);
        Bezier[2] sy = inner[Y].subdivide(t);
        return [BezierCurveN!N(sx[0], sy[0]), BezierCurveN!N(sx[1], sy[1])];
    }
    
    /+ Curve interface +/

    BezierCurveN!N duplicate() const
    { return new BezierCurveN!N(this); }

    BezierCurveN!N portion(Coord f, Coord t) const
    {
        static if (N == 1) {
            return new BezierCurveN!1(pointAt(f), pointAt(t));
        } else {
            return new BezierCurveN(inner[X].portion(f, t), inner[Y].portion(f, t));
        }
    }

    BezierCurveN!N reverse() const
    {
        static if (degree == 1) {
            return new BezierCurveN!1(finalPoint(), initialPoint()); 
        } else {
            return new BezierCurveN!N(inner[X].reverse(), inner[Y].reverse());
        }
    }

    BezierCurveN!N transformed(in Affine m) const
    {
        static if (degree == 1) {
            return new BezierCurveN!1(initialPoint() * m, finalPoint() * m);
        } else {
            BezierCurveN ret = new BezierCurveN();
            Point[] ps = points();
            for (size_t i = 0;  i <= degree; i++) {
                ps[i] = ps[i] * m;
            }
            ret.setPoints(ps);
            return ret;
        }
    }

    BezierCurve!N opBinary(string op)(in Translate m) const if (op == "*")
    {
        BezierCurve!N ret = new BezierCurve!N(this);
        ret.inner += m.vector();
        return ret;
    }

    BezierCurve!N derivative() const
    { return new BezierCurveN!(N-1)(inner[X].derivative(), inner[Y].derivative()); }
}

/** Compute the length of a bezier curve given by an array of its control points */
Coord bezier_length(in Point[] points, Coord tolerance)
{
    if (points.length < 2) return 0;
    Point[] v1 = points.dup;
    return bezier_length_internal(v1, tolerance);
}

private Coord bezier_length_internal(ref Point[] v1, Coord tolerance)
{
    /* The Bezier length algorithm used in 2Geom utilizes a simple fact:
     * the Bezier curve is longer than the distance between its endpoints
     * but shorter than the length of the polyline formed by its control
     * points. When the difference between the two values is smaller than the
     * error tolerance, we can be sure that the true value is no further than
     * 2*tolerance from their arithmetic mean. When it's larger, we recursively
     * subdivide the Bezier curve into two parts and add their lengths.
     */
    Coord lower = distance(v1[0], v1[$]);
    Coord upper = 0.0;
    for (size_t i = 0; i < v1.length - 1; ++i) {
        upper += distance(v1[i], v1[i+1]);
    }

    // termination point for recursion
    if (upper - lower < 2*tolerance) {
        return (lower + upper) / 2;
    }
        

    Point[] v2 = v1.dup;

    /* Compute the right subdivision directly in v1 and the left one in v2.
     * Explanation of the algorithm used:
     * We have to compute the left and right edges of this triangle in which
     * the top row are the control points of the Bezier curve, and each cell
     * is equal to the arithmetic mean of the cells directly above it
     * to the right and left. This corresponds to subdividing the Bezier curve
     * at time value 0.5: the left edge has the control points of the first
     * portion of the Bezier curve and the right edge - the second one.
     * In the example we subdivide a curve with 5 control points (order 4).
     *
     * Start:
     * 0 1 2 3 4
     *  ? ? ? ?
     *   ? ? ?
     *    ? ?
     *     ?
     * # means we have overwritten the value, ? means we don't know
     * the value yet. Numbers mean the value is at i-th position in the vector.
     *
     * After loop with i==1
     * # 1 2 3 4
     *  0 ? ? ? -> write 0 to v2[1]
     *   ? ? ?
     *    ? ?
     *     ?
     *
     * After loop with i==2
     * # # 2 3 4
     *  # 1 ? ?
     *   0 ? ? -> write 0 to v2[2]
     *    ? ?
     *     ?
     *
     * After loop with i==3
     * # # # 3 4
     *  # # 2 ?
     *   # 1 ?
     *    0 ? -> write 0 to v2[3]
     *     ?
     *
     * After loop with i==4, we have the right edge of the triangle in v1,
     * and we write the last value needed for the left edge in v2[4].
     */

    for (size_t i = 1; i < v1.length; ++i) {
        for (size_t j = i; j > 0; --j) {
            v1[j-1] = (v1[j-1] + v1[j]) * 0.5;
        }
        v2[i] = v1[0];
    }

    return bezier_length_internal(v1, 0.5*tolerance) + bezier_length_internal(v2, 0.5*tolerance);
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
