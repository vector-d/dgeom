/**
 * Cartesian point / 2D vector with integer coordinates
 *
 * Copyright 2011 Krzysztof Kosiński <tweenk.pl@gmail.com>
 *
 * This library is free software; you can redistribute it and/or
 * modify it either under the terms of the GNU Lesser General Public
 * License version 2.1 as published by the Free Software Foundation
 * (the "LGPL") or, at your option, under the terms of the Mozilla
 * Public License Version 1.1 (the "MPL"). If you do not alter this
 * notice, a recipient may use your version of this file under either
 * the MPL or the LGPL.
 *
 * You should have received a copy of the LGPL along with this library
 * in the file COPYING-LGPL-2.1; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 * You should have received a copy of the MPL along with this library
 * in the file COPYING-MPL-1.1
 *
 * The contents of this file are subject to the Mozilla Public License
 * Version 1.1 (the "License"); you may not use this file except in
 * compliance with the License. You may obtain a copy of the License at
 * http://www.mozilla.org/MPL/
 *
 * This software is distributed on an "AS IS" basis, WITHOUT WARRANTY
 * OF ANY KIND, either express or implied. See the LGPL or the MPL for
 * the specific language governing rights and limitations.
 */

module geom.intpoint;

public import geom.coord;

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
