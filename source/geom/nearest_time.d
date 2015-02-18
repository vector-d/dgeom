/**
 *  Nearest time routines for D2<SBasis> and Piecewise<D2<SBasis>>
 *
 * Authors:
 *   Marco Cecchetti <mrcekets at gmail.com>
 *   Liam P. White
 *
 * Copyright 2007-2015  authors
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

module geom.nearest_time;

import std.algorithm : swap, uniq;

public import geom.coord;
import geom.d2;
import geom.piecewise;
import geom.point;
import geom.sbasis;
import geom.sbasis_roots;

////////////////////////////////////////////////////////////////////////////////
// D2<SBasis> versions

/*
 * Given a line L specified by a point A and direction vector v,
 * return the point on L nearest to p. Note that the returned value
 * is with respect to the _normalized_ direction of v!
 */
Coord nearest_time(in Point p, in Point A, in Point v)
{
    Point d = Point(p - A);
    return d[0] * v[0] + d[1] * v[1];
}

/*
 * Returns the parameter t of a nearest point on the portion of the curve "c",
 * related to the interval [from, to], to the point "p".
 * The needed curve derivative "deriv" is passed as a parameter.
 * The function returns the first nearest point to "p" that is found.
 */
Coord nearest_time(in Point p, in D2!SBasis c, in D2!SBasis dc, Coord from = 0, Coord to = 1)
{
    if (from > to) swap(from, to);
    if (from < 0 || to > 1) {
        throw new Exception("[from,to] interval out of bounds");
    }

    if (c.isConstant()) return from;

    SBasis dd = (c - p).dot(dc);
    Coord[] zeros = roots(dd);

    Coord closest = from;
    Coord min_dist_sq = L2sq(c(from) - p);

    foreach (i; zeros) {
        Coord distsq = L2sq(c(i) - p);
        if (min_dist_sq > distsq) {
            closest = i;
            min_dist_sq = distsq;
        }
    }
    if (min_dist_sq > L2sq(c(to) - p))
        closest = to;

    return closest;
}


Coord nearest_time(in Point p, in D2!SBasis c, Coord from = 0, Coord to = 1)
{ return nearest_time(p, c, D2!SBasis(c[X].derivative(), c[Y].derivative()), from, to); }

/*
 * Return the parameters t of all the nearest times on the portion of
 * the curve "c", related to the interval [from, to], to the point "p".
 * The needed curve derivative "dc" is passed as a parameter.
 */
Coord[] all_nearest_times(in Point p, in D2!SBasis c, in D2!SBasis dc, Coord from = 0, Coord to = 1)
{
    import object : reserve;

    if (from > to) swap(from, to);
    if (from < 0 || to > 1) {
        throw new Exception("[from,to] interval out of bounds");
    }

    Coord[] result;
    if (c.isConstant()) {
        result ~= from;
        return result;
    }

    SBasis dd = (c - p).dot(dc);

    Coord[] zeros = roots(dd);
    Coord[] candidates;
    candidates ~= from ~ zeros ~ to;

    Coord[] distsq;
    distsq.reserve(candidates.length);

    foreach (i; candidates) {
        distsq ~= L2sq(c(i) - p);
    }

    size_t closest = 0;
    Coord dsq = distsq[0];
    for (size_t i = 1; i < candidates.length; ++i) {
        if (dsq > distsq[i]) {
            closest = i;
            dsq = distsq[i];
        }
    }

    for (size_t i = 0; i < candidates.length; ++i) {
        if (distsq[closest] == distsq[i]) {
            result ~= candidates[i];
        }
    }

    return result;
}

Coord[] all_nearest_times(in Point p, in D2!SBasis c, Coord from = 0, Coord to = 1)
{ return all_nearest_times(p, c, D2!SBasis(c[X].derivative(), c[Y].derivative()), from, to); }

////////////////////////////////////////////////////////////////////////////////
// Piecewise< D2<SBasis> > versions

Coord nearest_time(in Point p, in Piecewise!(D2!SBasis) c, Coord from, Coord to)
{
    import geom.rect;

    if (from > to) swap(from, to);
    if (from < c.cuts[0] || to > c.cuts[c.size()]) {
        throw new Exception("[from,to] interval out of bounds");
    }

    size_t si = c.segN(from);
    size_t ei = c.segN(to);
    if (si == ei) {
        Coord nearest = nearest_time(p, c[si], c.segT(from, si), c.segT(to, si));
        return c.mapToDomain(nearest, si);
    }

    Coord t;
    Coord nearest = nearest_time(p, c[si], c.segT(from, si));
    size_t ni = si;
    Coord dsq;
    Coord mindistsq = geom.point.distanceSq(p, c[si](nearest));

    Rect bb = Rect.empty();
    for (size_t i = si + 1; i < ei; ++i) {
        bb = Rect(c[i][X].bounds_fast(), c[i][Y].bounds_fast());
        dsq = distanceSq(p, bb);
        if (mindistsq <= dsq) continue;

        t = nearest_time(p, c[i]);
        dsq = geom.point.distanceSq(p, c[i](t));
        if (mindistsq > dsq) {
            nearest = t;
            ni = i;
            mindistsq = dsq;
        }
    }

    bb = Rect(c[ei][X].bounds_fast(), c[ei][Y].bounds_fast());
    dsq = distanceSq(p, bb);
    if (mindistsq > dsq) {
        t = nearest_time(p, c[ei], 0, c.segT(to, ei));
        dsq = geom.point.distanceSq(p, c[ei](t));
        if (mindistsq > dsq) {
            nearest = t;
            ni = ei;
        }
    }
    return c.mapToDomain(nearest, ni);
}

Coord nearest_time(in Point p, in Piecewise!(D2!SBasis) c)
{
    return nearest_time(p, c, c.cuts[0], c.cuts[$]);
}


Coord[] all_nearest_times(in Point p, in Piecewise!(D2!SBasis) c, Coord from, Coord to)
{
    import std.array : array;
    import geom.rect;

    if (from > to) swap(from, to);
    if (from < c.cuts[0] || to > c.cuts[c.size()]) {
        throw new Exception("[from,to] interval out of bounds");
    }

    size_t si = c.segN(from);
    size_t ei = c.segN(to);
    if (si == ei) {
        Coord[] all_nearest = all_nearest_times(p, c[si], c.segT(from, si), c.segT(to, si));
        foreach(ref i; all_nearest) {
            i = c.mapToDomain(i, si);
        }
        return all_nearest;
    }

    Coord[] all_t;
    Coord[][] all_np;
    all_np ~= all_nearest_times(p, c[si], c.segT(from, si));
    size_t[] ni;
    ni ~= si;
    Coord dsq;
    Coord mindistsq = geom.point.distanceSq(p, c[si](all_np[0][0]));
    Rect bb = Rect.empty();

    for (size_t i = si + 1; i < ei; ++i) {
        bb = Rect(c[i][X].bounds_fast(), c[i][Y].bounds_fast());
        dsq = distanceSq(p, bb);
        if (mindistsq < dsq) continue;
        all_t = all_nearest_times(p, c[i]);
        dsq = geom.point.distanceSq( p, c[i](all_t[0]));

        if (mindistsq > dsq) {
            all_np = [];
            all_np ~= all_t;
            ni = [];
            ni ~= i;
            mindistsq = dsq;
        } else if (mindistsq == dsq) {
            all_np ~= all_t;
            ni ~= i;
        }
    }

    bb = Rect(c[ei][X].bounds_fast(), c[ei][Y].bounds_fast());
    dsq = distanceSq(p, bb);
    if (mindistsq >= dsq) {
        all_t = all_nearest_times(p, c[ei], 0, c.segT(to, ei));
        dsq = geom.point.distanceSq( p, c[ei](all_t[0]) );
        if (mindistsq > dsq) {
            foreach (ref i; all_t) {
                i = c.mapToDomain(i, ei);
            }
            return all_t;
        } else if (mindistsq == dsq) {
            all_np ~= all_t;
            ni ~= ei;
        }
    }

    Coord[] all_nearest;
    foreach (i, x; all_np) {
        foreach (j; x) {
            all_nearest ~= c.mapToDomain(j, ni[i]);
        }
    }

    all_nearest = uniq(all_nearest).array;
    return all_nearest;
}

Coord[] all_nearest_times(in Point p, in Piecewise!(D2!SBasis) c)
{ return all_nearest_times(p, c, c.cuts[0], c.cuts[$]); }

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
