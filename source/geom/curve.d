/*
 * Curve interface.
 *
 *
 * Authors:
 *   MenTaLguY <mental@rydia.net>
 *   Marco Cecchetti <mrcekets at gmail.com>
 *   Krzysztof Kosiński <tweenk.pl@gmail.com>
 * 
 * Copyright 2007-2009 Authors
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

module geom.curve;

import geom.coord;
import geom.affine;
import geom.d2;
import geom.interval;
import geom.point;
import geom.rect;
import geom.sbasis;
import geom.transforms;

/**
 * Abstract continuous curve on a plane defined on [0,1].
 *
 * Formally, a curve in 2Geom is defined as a function
 * \f$\mathbf{C}: [0, 1] \to \mathbb{R}^2\f$
 * (a function that maps the unit interval to points on a 2D plane). Its image (the set of points
 * the curve passes through) will be denoted \f$\mathcal{C} = \mathbf{C}[ [0, 1] ]\f$.
 * All curve types available in 2Geom are continuous and differentiable on their
 * interior, e.g. \f$(0, 1)\f$. Sometimes the curve's image (value set) is referred to as the curve
 * itself for simplicity, but keep in mind that it's not strictly correct.
 * 
 * It is common to think of the parameter as time. The curve can then be interpreted as
 * describing the position of some moving object from time \f$t=0\f$ to \f$t=1\f$.
 * Because of this, the parameter is frequently called the time value.
 *
 * Some methods return pointers to newly allocated curves. They are expected to be freed
 * by the caller when no longer used. Default implementations are provided for some methods.
 *
 * @ingroup Curves
 */
interface Curve
{
    /+ Evaluate the curve +/

    /** Retrieve the start of the curve.
    * @return The point corresponding to \f$\mathbf{C}(0)\f$. */
    Point initialPoint() const;

    /** Retrieve the end of the curve.
     * @return The point corresponding to \f$\mathbf{C}(1)\f$. */
    Point finalPoint() const;

    /** Check whether the curve has exactly zero length.
     * @return True if the curve's initial point is exactly the same as its final point, and it contains
     *         no other points (its value set contains only one element).
     */
    bool isDegenerate() const;

    /** Evaluate the curve at a specified time value.
     * @param t Time value
     * @return \f$\mathbf{C}(t)\f$ */
    Point pointAt(Coord t) const;
    /+ { return pointAndDerivatives(t, 0)[0]; } +/

    /** Evaluate one of the coordinates at the specified time value.
     * @param t Time value
     * @param d The dimension to evaluate
     * @return The specified coordinate of \f$\mathbf{C}(t)\f$ */
    Coord valueAt(Coord t, size_t d) const;
    /+  { return pointAt(t)[d]; } +/

    /** Evaluate the function at the specified time value. Allows curves to be used
     * as functors. */
    final Point opCall(Coord t) const { return pointAt(t); }

    /** Evaluate the curve and its derivatives.
     * This will return a vector that contains the value of the curve and the specified number
     * of derivatives. However, the returned vector might contain less elements than specified
     * if some derivatives do not exist.
     * @param t Time value
     * @param n The number of derivatives to compute
     * @return Vector of at most \f$n+1\f$ elements of the form \f$[\mathbf{C}(t),
        \mathbf{C}'(t), \mathbf{C}''(t), \ldots]\f$ */
    Point[] pointAndDerivatives(Coord t, uint n) const;

    /+ Change the curve's endpoints +/

    /** Change the starting point of the curve.
     * After calling this method, it is guaranteed that \f$\mathbf{C}(0) = \mathbf{p}\f$,
     * and the curve is still continuous. The precise new shape of the curve varies with curve
     * type.
     * @param p New starting point of the curve */
    void setInitial(in Point v);

    /** Change the ending point of the curve.
     * After calling this method, it is guaranteed that \f$\mathbf{C}(0) = \mathbf{p}\f$,
     * and the curve is still continuous. The precise new shape of the curve varies
     * with curve type.
     * @param p New ending point of the curve */
    void setFinal(in Point v);

    /+ Compute the bounding box +/

    /** Quickly compute the curve's approximate bounding box.
     * The resulting rectangle is guaranteed to contain all points belonging to the curve,
     * but it might not be the smallest such rectangle. This method is usually fast.
     * @return A rectangle that contains all points belonging to the curve. */

    Rect boundsFast() const;
    /** Compute the curve's exact bounding box.
     * This method can be dramatically slower than boundsExact() depending on the curve type.
     * @return The smallest possible rectangle containing all of the curve's points. */
    Rect boundsExact() const;

    // I have no idea what the 'deg' parameter is for, so this is undocumented for now.
    Rect boundsLocal(in Interval i, uint deg) const;

    /** Compute the bounding box of a part of the curve.
     * Since this method returns the smallest possible bounding rectangle of the specified portion,
     * it can also be rather slow.
     * @param a An interval specifying a part of the curve, or nothing.
     *          If \f$[0, 1] \subseteq a\f$, then the bounding box for the entire curve
     *          is calculated.
     * @return The smallest possible rectangle containing all points in \f$\mathbf{C}[a]\f$,
     *         or nothing if the supplied interval is empty. */
    final Rect boundsLocal(in Interval a) const { return boundsLocal(a, 0); }


    /+ Create new curves based on this one
     + IMPLEMENTORS' NOTE: please use covariant return types! +/

    /** Create an exact copy of this curve.
     * @return A newly allocated curve, identical to the original */
    Curve duplicate() const;

    /** Create a curve transformed by an affine transformation.
     * This method returns a new curve instead modifying the existing one, because some curve
     * types are not closed under affine transformations. The returned curve may be of different
     * underlying type (as is the case for horizontal and vertical line segments).
     * @param m Affine describing the affine transformation
     * @return Pointer to a new, transformed curve */
    void transform(in Affine m);

    /** Translate the curve (i.e. displace by Point)
     * This method modifies the curve; all curve types are closed under
     * translations (the result can be expressed in its own curve type).
     * This function yields the same result as transformed(m).
     * @param p Point by which to translate the curve
     * @return reference to self */
    Curve opBinary(string op)(Translate rhs) const if (op == "*");

    /** Create a curve that corresponds to a part of this curve.
     * For \f$a > b\f$, the returned portion will be reversed with respect to the original.
     * The returned curve will always be of the same type.
     * @param a Beginning of the interval specifying the portion of the curve
     * @param b End of the interval
     * @return New curve \f$\mathbf{D}\f$ such that:
     * - \f$\mathbf{D}(0) = \mathbf{C}(a)\f$
     * - \f$\mathbf{D}(1) = \mathbf{C}(b)\f$
     * - \f$\mathbf{D}[ [0, 1] ] = \mathbf{C}[ [a?b] ]\f$,
     *   where \f$[a?b] = [\min(a, b), \max(a, b)]\f$ */
    Curve portion(Coord a, Coord b) const;

    /** A version of that accepts an Interval. */
    final Curve portion(in Interval i) const { return portion(i.min(), i.max()); }

    /** Create a reversed version of this curve.
     * The result corresponds to <code>portion(1, 0)</code>, but this method might be faster.
     * @return Pointer to a new curve \f$\mathbf{D}\f$ such that
     *         \f$\forall_{x \in [0, 1]} \mathbf{D}(x) = \mathbf{C}(1-x)\f$ */
    final Curve reverse() const { return portion(1, 0); }

    /** Create a derivative of this curve.
     * It's best to think of the derivative in physical terms: if the curve describes
     * the position of some object on the plane from time \f$t=0\f$ to \f$t=1\f$ as said in the
     * introduction, then the curve's derivative describes that object's speed at the same times.
     * The second derivative refers to its acceleration, the third to jerk, etc.
     * @return New curve \f$\mathbf{D} = \mathbf{C}'\f$. */
    Curve derivative() const;

    /+ Advanced operations +/

    /** Compute a time value at which the curve comes closest to a specified point.
     * The first value with the smallest distance is returned if there are multiple such points.
     * @param p Query point
     * @param a Minimum time value to consider
     * @param b Maximum time value to consider; \f$a < b\f$
     * @return \f$q \in [a, b]: ||\mathbf{C}(q) - \mathbf{p}|| = 
               \inf(\{r \in \mathbb{R} : ||\mathbf{C}(r) - \mathbf{p}||\})\f$ */
    final Coord nearestPoint(in Point p, Coord a = 0, Coord b = 1) const
    {
        import geom.nearest_time;
        return nearest_time(p, toSBasis(), a, b);
    }

    /** A version that takes an Interval. */
    final Coord nearestPoint(in Point p, in Interval i) const
    { return nearestPoint(p, i.min(), i.max()); }

    /** Compute time values at which the curve comes closest to a specified point.
     * @param p Query point
     * @param a Minimum time value to consider
     * @param b Maximum time value to consider; \f$a < b\f$
     * @return Vector of points closest and equally far away from the query point */
    final Coord[] allNearestPoints(in Point p, Coord from = 0, Coord to = 1) const
    {
        import geom.nearest_time;
        return all_nearest_times(p, toSBasis(), from, to);
    }

    /** A version that takes an Interval. */
    final Coord[] allNearestPoints(in Point p, in Interval i)
    { return allNearestPoints(p, i.min(), i.max()); }

    /** Compute the arc length of this curve.
     * For a curve \f$\mathbf{C}(t) = (C_x(t), C_y(t))\f$, arc length is defined for 2D curves as
     * \f[ \ell = \int_{0}^{1} \sqrt { [C_x'(t)]^2 + [C_y'(t)]^2 }\, \text{d}t \f]
     * In other words, we divide the curve into infinitely small linear segments
     * and add together their lengths. Of course we can't subdivide the curve into
     * infinitely many segments on a computer, so this method returns an approximation.
     * Note that there is usually no closed form solution to such integrals, so this
     * method might be slow.
     * @param tolerance Maximum allowed error
     * @return Total distance the curve's value travels on the plane when going from 0 to 1 */
    Coord length(Coord tolerance = 0.01) const;

    /** Computes time values at which the curve intersects an axis-aligned line.
     * @param v The coordinate of the line
     * @param d Which axis the coordinate is on. X means a vertical line, Y a horizontal line. */
    Coord[] roots(Coord v, size_t d) const;

    /** Compute the winding number of the curve at the specified point.
     * Does this makes sense for curves at all? */
    final int winding(in Point p) const
    {
        import std.algorithm : sort;

        Coord[] ts = roots(p[Y], Y);
        if (ts.length == 0) return 0;
        sort(ts);

        // skip endpoint roots when they are local maxima on the Y axis
        // this follows the convention used in other winding routines,
        // i.e. that the bottommost coordinate is not part of the shape
        bool ignore_0 = unitTangentAt(0)[Y] <= 0;
        bool ignore_1 = unitTangentAt(1)[Y] >= 0;

        int wind = 0;
        foreach (t; ts) {
            if ((t == 0 && ignore_0) || (t == 1) && ignore_1) continue;
            if (valueAt(t, X) > p[X]) { // root is ray intersection
                Point tangent = unitTangentAt(t);
                if (tangent[Y] > p[X]) {
                    // at the point of intersection, curve goes in +Y direction,
                    // so it winds in the direction of positive angles
                    ++wind;
                } else if (tangent[Y] < 0) {
                    --wind;
                }
            }
        }
        return wind;
    }

    /** Compute a vector tangent to the curve.
     * This will return an unit vector (a Point with length() equal to 1) that denotes a vector
     * tangent to the curve. This vector is defined as
     * \f$ \mathbf{v}(t) = \frac{\mathbf{C}'(t)}{||\mathbf{C}'(t)||} \f$. It is pointed
     * in the direction of increasing \f$t\f$, at the specfied time value. The method uses
     * l'Hopital's rule when the derivative is zero. A zero vector is returned if no non-zero
     * derivative could be found.
     * @param t Time value
     * @param n The maximum order of derivative to consider
     * @return Unit tangent vector \f$\mathbf{v}(t)\f$
     * @bug This method might currently break for the case of t being exactly 1. A workaround
     *      is to reverse the curve and use the negated unit tangent at 0 like this:
     * @code
        Curve c_reverse = c.reverse();
        Point tangent = - c_reverse.unitTangentAt(0); @endcode */
    final Point unitTangentAt(Coord t, uint n = 3) const
    {
        Point[] derivs = pointAndDerivatives(t, n);
        foreach (deriv; derivs[1 .. $]) {
            Coord length = deriv.length();
            if (!are_near(length, 0)) {
                return deriv / length;
            }
        }
        return Point(0, 0);
    }

    /** Convert the curve to a symmetric power basis polynomial.
     * Symmetric power basis polynomials (S-basis for short) are numerical representations
     * of curves with excellent numerical properties. Most high level operations provided by 2Geom
     * are implemented in terms of S-basis operations, so every curve has to provide a method
     * to convert it to an S-basis polynomial on two variables. See SBasis class reference
     * for more information. */
    D2!SBasis toSBasis() const;
    
    /+  Miscellaneous +/

    /** Return the number of independent parameters required to represent all variations
     * of this curve. For example, for Bezier curves it returns the curve's order
     * multiplied by 2. */
    int degreesOfFreedom() const;
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
