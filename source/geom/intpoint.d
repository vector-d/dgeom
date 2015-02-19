/*
 * Cartesian point / 2D vector with integer coordinates
 *
 * Authors:
 *   Krzysztof Kosiński <tweenk.pl@gmail.com>
 *   Liam P. White
 *
 * Copyright (C) 2011-2015 Authors
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

module geom.intpoint;

import geom.coord;

/**
 * Two-dimensional point with integer coordinates.
 *
 * This class is an exact equivalent of Point, except it stores integer coordinates.
 * Integer points are useful in contexts related to rasterized graphics, for example
 * for bounding boxes when rendering SVG.
 *
 * See: Point
 * In group: Primitives */
struct IntPoint
{
    /* Creating integer points */
    this(IntCoord x, IntCoord y)
    { _pt = [x, y]; }
    
    this(const(IntCoord[2]) arr)
    { _pt = arr; }


    /* Access the coordinates of a point */
    ref inout(IntCoord) opIndex(size_t i) inout
    { return _pt[i]; }

    ref inout(IntCoord) x() inout { return _pt[X]; }
    ref inout(IntCoord) y() inout { return _pt[Y]; }

    /* Various utilities */
    IntPoint opBinary(string op)(IntPoint rhs) const
    {
        static if (op == "+") { return IntPoint(_pt[X] + rhs[X], _pt[Y] + rhs[Y]); }
        else static if (op == "-") { return IntPoint(_pt[X] - rhs[X], _pt[Y] - rhs[Y]); }
        else static assert(0, "IntPoint operator "~op~" not implemented");
    }
    
    void opOpAssign(string op)(IntPoint rhs)
    { mixin("this = this "~op~" rhs;"); }

    /** Lexicographical ordering for points.
     * Y coordinate is regarded as more significant. When sorting according to this
     * ordering, the points will be sorted according to the Y coordinate, and within
     * points with the same Y coordinate according to the X coordinate. */
    int opCmp(ref const(IntPoint) rhs) const
    { return ((_pt[Y] < rhs[Y]) || (( _pt[Y] == rhs[Y] ) && ( _pt[X] < rhs[X] ))) ? -1 : 1; }

    private IntCoord[2] _pt;
}

unittest
{
    IntPoint p = IntPoint(0, 0);
    IntPoint q = [1, 1];
    IntPoint r = IntPoint(2, 2);
    
    assert(p[X] == 0);
    assert(q[Y] == 1);
    
    IntPoint s = q + r;
    assert(s.x() == 3);
    
    s += p;
    assert(s.y() == 3);
    
    assert(p < s);
    assert(s > p);
    
    IntPoint t = [0, 0];
    assert(t == p);
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
