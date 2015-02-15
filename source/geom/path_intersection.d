/*
 * Path intersection
 *
 * Authors:
 *   ? <?@?.?>
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

module geom.path_intersection;

import geom.affine;
import geom.coord;
import geom.curve;
import geom.crossing;
import geom.interval;
import geom.path;
import geom.path_sequence;
import geom.point;
import geom.rect;

bool contains(in Path p, in Point i, bool evenodd = true)
{
    return (evenodd ? p.winding(i) % 2 : p.winding(i)) != 0;
}

Crossings curve_sweep(T)(in Path a, in Path b)
{
    import std.container : make;

    T t = make!T();
    Crossings ret;
    Rect[] bounds_a = bounds(a), bounds_b = bounds(b);
    size_t[][] ixs = sweep_bounds(bounds_a, bounds_b);
    foreach (i; 0 .. a.size) {
        foreach (jp; ixs[i]) {
            Crossings cc = t.crossings(a[i], b[jp]);
            offset_crossings(cc, i, jp);
            ret ~= cc;
        }
    }
    return ret;
}

class SimpleCrosser : Crosser!Path
{
    Crossings crossings(in Curve a, in Curve b) const
    {
        Crossings ret;
        pair_intersect(a, 0, 1, b, 0, 1, ret);
        return ret;
    }
    override Crossings crossings(in Path a, in Path b) const { return curve_sweep!SimpleCrosser(a, b); }
    CrossingSet crossings(in PathSequence a, in PathSequence b) const { return Crosser!(Path).crossings(a._data, b._data); }
}

alias DefaultCrosser = SimpleCrosser;

Crossings crossings(in Curve a, in Curve b)
{
    DefaultCrosser c = new DefaultCrosser();
    return c.crossings(a, b);
}

Crossings crossings(in Path a, in Path b)
{
    DefaultCrosser c = new DefaultCrosser();
    return c.crossings(a, b);
}

CrossingSet crossings(in PathSequence a, in PathSequence b)
{
    DefaultCrosser c = new DefaultCrosser();
    return c.crossings(a, b);
}

/** Finds the intersection between the lines defined by A0 & A1, and B0 & B1.
 * Returns through the last 3 parameters, returning the t-values on the lines
 * and the cross-product of the deltas (a useful byproduct).  The return value
 * indicates if the time values are within their proper range on the line segments.
 */
bool linear_intersect(in Point A0, in Point A1, in Point B0, in Point B1, ref Coord tA, ref Coord tB, ref Coord det)
{
    import std.math : fabs;

    bool both_lines_non_zero = (!are_near(A0, A1)) && (!are_near(B0, B1));

    // Cramer's rule as cross products
    Point Ad = A1 - A0,
          Bd = B1 - B0,
           d = B0 - A0;
    det = cross(Ad, Bd);

    Coord det_rel = det; // Calculate the determinant of the normalized vectors
    if (both_lines_non_zero) {
        det_rel /= Ad.length();
        det_rel /= Bd.length();
    }

    if (fabs(det_rel) < 1e-12) { // If the cross product is NEARLY zero,
        // Then one of the linesegments might have length zero
        if (both_lines_non_zero) {
            // If that's not the case, then we must have either:
            // - parallel lines, with no intersections, or
            // - coincident lines, with an infinite number of intersections
            // Either is quite useless, so we'll just bail out
            return false;
        } // Else, one of the linesegments is zero, and we might still be able to calculate a single intersection point
    } // Else we haven't bailed out, and we'll try to calculate the intersections

    Coord detinv = 1.0 / det;
    tA = cross(d, Bd) * detinv;
    tB = cross(d, Ad) * detinv;
    return (tA >= 0.) && (tA <= 1.) && (tB >= 0.) && (tB <= 1.);
}

/**
 * This uses the local bounds functions of curves to generically intersect two.
 * It passes in the curves, time intervals, and keeps track of depth, while
 * returning the results through the Crossings parameter.
 */
void pair_intersect(in Curve A, Coord Al, Coord Ah, in Curve B, Coord Bl, Coord Bh, ref Crossings ret, uint depth = 0)
{
    Rect Ar = A.boundsLocal(Interval(Al, Ah));
    if (Ar.isEmpty) return;

    Rect Br = B.boundsLocal(Interval(Bl, Bh));
    if (Br.isEmpty) return;
    
    if (!Ar.intersects(Br)) return;
    
    // Checks the general linearity of the function
    if (depth > 12) { // || (A.boundsLocal(Interval(Al, Ah), 1).maxExtent() < 0.1 
                    //&&  B.boundsLocal(Interval(Bl, Bh), 1).maxExtent() < 0.1)) {
        Coord tA, tB, c;
        if (linear_intersect(A.pointAt(Al), A.pointAt(Ah), B.pointAt(Bl), B.pointAt(Bh), tA, tB, c)) {
            tA = tA * (Ah - Al) + Al;
            tB = tB * (Bh - Bl) + Bl;
            intersect_polish_root(A, tA, B, tB);
            if (depth % 2)
                ret ~= new Crossing(tB, tA, c < 0);
            else
                ret ~= new Crossing(tA, tB, c > 0);
            return;
        }
    }
    if (depth > 12) return;
    Coord mid = (Bl + Bh)/2;
    pair_intersect(B, Bl, mid, A, Al, Ah, ret, depth+1);
    pair_intersect(B, mid, Bh, A, Al, Ah, ret, depth+1);
}

Crossings pair_intersect(in Curve A, in Interval Ad, in Curve B, in Interval Bd)
{
    Crossings ret;
    pair_intersect(A, Ad.min(), Ad.max(), B, Bd.min(), Bd.max(), ret);
    return ret;
}

private void intersect_polish_root(in Curve A, ref Coord s, in Curve B, ref Coord t)
{
    Point[] as, bs;
    as = A.pointAndDerivatives(s, 2);
    bs = B.pointAndDerivatives(t, 2);
    Point F = as[0] - bs[0];
    double best = dot(F, F);
    
    foreach (i; 0 .. 4) {
        
        /**
           we want to solve
           J*(x1 - x0) = f(x0)
           
           |dA(s)[0]  -dB(t)[0]|  (X1 - X0) = A(s) - B(t)
           |dA(s)[1]  -dB(t)[1]| 
        **/

        // We're using the standard transformation matricies, which is numerically rather poor.  Much better to solve the equation using elimination.

        Affine jack = Affine(as[1][0], as[1][1], -bs[1][0], -bs[1][1], 0, 0);
        Point soln = F*jack.inverse();
        Coord ns = s - soln[0];
        Coord nt = t - soln[1];

        if (ns < 0) ns = 0;
        else if (ns > 1) ns = 1;
        if (nt < 0) nt = 0;
        else if (nt > 1) nt = 1;
        
        as = A.pointAndDerivatives(ns, 2);
        bs = B.pointAndDerivatives(nt, 2);
        F = as[0] - bs[0];
        Coord trial = dot(F, F);
        if (trial > best*0.1) // we have standards, you know
            // At this point we could do a line search
            break;
        best = trial;
        s = ns;
        t = nt;
    }
}

Crossings self_crossings(in Path p)
{
    Crossings ret;
    size_t[][] cull = sweep_bounds(bounds(p));
    foreach (i; 0 .. cull.length) {
        Crossings res = curve_self_crossings(p[i]);
        offset_crossings(res, i, i);
        ret ~= res;
        foreach (jx; 0 .. cull[i].length) {
            size_t j = cull[i][jx];
            res = [];
            pair_intersect(p[i], 0, 1, p[j], 0, 1, res);
            
            //if(fabs(int(i)-j) == 1 || fabs(int(i)-j) == p.size()-1) {
                Crossings res2;
                foreach (k; 0 .. res.length) {
                    if (res[k].ta != 0 && res[k].ta != 1 && res[k].tb != 0 && res[k].tb != 1) {
                        res2 ~= res[k];
                    }
                }
                res = res2;
            //}
            offset_crossings(res, i, j);
            ret ~= res;
        }
    }
    return ret;
}

Crossings curve_self_crossings(in Curve a)
{
    Crossings res;
    Coord[] spl = 0 ~ curve_mono_splits(a) ~ 1;
    foreach (i; 1 .. spl.length)
        foreach (j; i+1 .. spl.length)
            pair_intersect(a, spl[i-1], spl[i], a, spl[j-1], spl[j], res);
    return res;
}

/** This returns the times when the x or y derivative is 0 in the curve. */
Coord[] curve_mono_splits(in Curve d)
{
    import std.algorithm : sort;

    Curve deriv = d.derivative();
    Coord[] rs = deriv.roots(0, X) ~ deriv.roots(0, Y);
    rs.sort();
    return rs;
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
