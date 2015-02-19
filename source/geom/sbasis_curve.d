/**
 *  Symmetric power basis curve
 *
 * Authors:
 *   MenTaLguY <mental@rydia.net>
 *   Marco Cecchetti <mrcekets at gmail.com>
 *   Krzysztof Kosiński <tweenk.pl@gmail.com>
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

module geom.sbasis_curve;

import geom.curve;
import geom.coord;
import geom.affine;
import geom.d2;
import geom.interval;
import geom.point;
import geom.rect;
import geom.sbasis;
import geom.transforms;
import geom.nearest_time;
import geom.sbasis_roots;

/** Symmetric power basis curve.
 *
 * Symmetric power basis (S-basis for short) polynomials are a versatile numeric
 * representation of arbitrary continuous curves. They are the main representation of curves
 * in 2Geom.
 *
 * S-basis is defined for odd degrees and composed of the following polynomials:
 * \f{align*}{
         P_k^0(t) &= t^k (1-t)^{k+1} \\
         P_k^1(t) &= t^{k+1} (1-t)^k \f}
 * This can be understood more easily with the help of the chart below. Each square
 * represents a product of a specific number of \f$t\f$ and \f$(1-t)\f$ terms. Red dots
 * are the canonical (monomial) basis, the green dots are the Bezier basis, and the blue
 * dots are the S-basis, all of them of degree 7.
 *
 * The S-Basis has several important properties:
 * - S-basis polynomials are closed under multiplication.
 * - Evaluation is fast, using a modified Horner scheme.
 * - Degree change is as trivial as in the monomial basis. To elevate, just add extra
 *   zero coefficients. To reduce the degree, truncate the terms in the highest powers.
 *   Compare this with Bezier curves, where degree change is complicated.
 * - Conversion between S-basis and Bezier basis is numerically stable.
 *
 * More in-depth information can be found in the following paper:
 * J Sanchez-Reyes, "The symmetric analogue of the polynomial power basis".
 * ACM Transactions on Graphics, Vol. 16, No. 3, July 1997, pages 319--357.
 * http://portal.acm.org/citation.cfm?id=256162
 *
 */

class SBasisCurve : Curve
{
    this(in D2!SBasis sb) { inner = D2!SBasis(sb); }
    this(in SBasisCurve m) { inner = D2!SBasis(m.inner); }
    this(in Curve other) { inner = other.toSBasis(); }

    /+ Curve interface +/
    SBasisCurve duplicate() const { return new SBasisCurve(this); }
    Point initialPoint() const    { return inner.at0(); }
    Point finalPoint() const      { return inner.at1(); }
    bool isDegenerate() const     { return inner.isConstant(); }
    Point pointAt(Coord t) const  { return inner.valueAt(t); }
    Point[] pointAndDerivatives(Coord t, uint n) const
    { return inner.valueAndDerivatives(t, n); }

    Coord valueAt(Coord t, size_t d) const { return inner[d].valueAt(t); }

    void setInitial(in Point v) { for (uint d = 0; d < 2; d++) { inner[d][0][0] = v[d]; } }
    void setFinal(in Point v) { for (uint d = 0; d < 2; d++) { inner[d][0][1] = v[d]; } }

    Rect boundsFast() const  { return Rect(inner[X].bounds_fast(), inner[Y].bounds_fast()); }
    Rect boundsExact() const { return Rect(inner[X].bounds_exact(),inner[Y].bounds_exact()); }
    Rect boundsLocal(in Interval i, uint deg) const { return Rect(inner[X].bounds_local(i, deg), inner[Y].bounds_local(i, deg)); }

    Coord[] roots(Coord v, size_t d) const { return geom.sbasis_roots.roots(inner[d] - v); }

    // XXX implement the integral solution for this
    Coord length(Coord tolerance = 0.01) const
    {
        import geom.bezier_curve, geom.sbasis_to_bezier;
        Point[] pts;
        sbasis_to_bezier(pts, inner);
        return bezier_length(pts, tolerance);
    }
    SBasisCurve portion(Coord f, Coord t) const { return new SBasisCurve(inner.portion(f, t)); }

    void transform(in Affine m) { inner *= m; }
    
    SBasisCurve opBinary(string op)(in Translate m) const if (op == "*")
    {
        SBasisCurve ret = new SBasisCurve(this);
        ret.inner += m.vector();
        return ret;
    }

    SBasisCurve derivative() const { return new SBasisCurve(D2!SBasis(inner[X].derivative(), inner[Y].derivative())); }
    D2!SBasis toSBasis() const { return D2!SBasis(inner); }

    int degreesOfFreedom() const { return inner[0].degreesOfFreedom() + inner[1].degreesOfFreedom(); }

    private D2!SBasis inner;
};

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
