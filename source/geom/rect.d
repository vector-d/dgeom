/*
 * Axis-aligned rectangle
 *
 * Authors:
 *   Michael Sloan <mgsloan@gmail.com>
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
 *
 * Authors of original rect class:
 *   Lauris Kaplinski <lauris@kaplinski.com>
 *   Nathan Hurst <njh@mail.csse.monash.edu.au>
 *   bulia byak <buliabyak@users.sf.net>
 *   MenTaLguY <mental@rydia.net>
 */

module geom.rect;

public import geom.coord;

import math = std.math;
import std.traits;

import geom.affine;
import geom.intpoint;
import geom.point;
import geom.interval;

struct GenericRect(C)
{
    alias CInterval = GenericInterval!C;
    static if (C.stringof == IntCoord.stringof) alias CPoint = IntPoint;
    static if (C.stringof == Coord.stringof) alias CPoint = Point;

    /** Create a rectangle that contains only the point at (0,0). */
    @disable this();

    /** Create a rectangle from X and Y intervals. */
    this(CInterval a, CInterval b) { f[X] = a;  f[Y] = b; }

    /** Create a rectangle from two points. */
    this(in CPoint a, in CPoint b)
    {
        f[X] = CInterval(a[X], b[X]);
        f[Y] = CInterval(a[Y], b[Y]);
    }

    /** Create rectangle from coordinates of two points. */
    this(C x0, C y0, C x1, C y1)
    {
        f[X] = CInterval(x0, x1);
        f[Y] = CInterval(y0, y1);
    }

    /** Create rectangle from origin and dimensions. */
    static GenericRect!C from_xywh(C x, C y, C w, C h)
    {
        CPoint xy = [x, y];
        CPoint wh = [w, h];
        return GenericRect!C(xy, xy + wh);
    }

    /** Create rectangle from origin and dimensions. */
    static GenericRect!C from_xywh(in CPoint xy, in CPoint wh)
    { return GenericRect!C(xy, xy + wh); }

    /** Create infinite rectangle. */
    static GenericRect!C infinite()() if (!isIntegral!C)
    { return GenericRect!C(-C.infinity, C.infinity, -C.infinity, C.infinity); }

    /** Specialization of infinity() for integers */
    static GenericRect!C infinite()() if (isIntegral!C)
    { return GenericRect!C(-C.max, C.max, -C.max, C.max); }

    /** Create an empty rectangle.
     * This rectangle is special in that intersecting with it will always return
     * the empty rectangle, and unioning will do nothing. */
    static GenericRect!C empty()
    { return GenericRect(CInterval.empty(), CInterval.empty()); }

    /+ Inspect dimensions +/

    /** Access an interval by its index. */
    ref inout(GenericInterval!C) opIndex(size_t i) inout
    { return f[i]; }

    /** Get the corner of the rectangle with smallest coordinate values.
     * In 2Geom standard coordinate system, this means upper left. */
    CPoint min() const { return CPoint(f[X].min(), f[Y].min()); }
    /** Get the corner of the rectangle with largest coordinate values.
     * In 2Geom standard coordinate system, this means lower right. */
    CPoint max() const { return CPoint(f[X].max(), f[Y].max()); }
    /** Return the n-th corner of the rectangle.
     * Returns corners in the direction of growing angles, starting from
     * the one given by min(). For the standard coordinate system used
     * in 2Geom (+Y downwards), this means clockwise starting from
     * the upper left. */
    CPoint corner(uint i) const
    {
        switch(i % 4) {
            case 0:  return CPoint(f[X].min(), f[Y].min());
            case 1:  return CPoint(f[X].max(), f[Y].min());
            case 2:  return CPoint(f[X].max(), f[Y].max());
            default: return CPoint(f[X].min(), f[Y].max());
        }
    }

    // We should probably remove these - they're coord sys gnostic

    /** Return top coordinate of the rectangle (+Y is downwards). */
    C top() const { return f[Y].min(); }
    /** Return bottom coordinate of the rectangle (+Y is downwards). */
    C bottom() const { return f[Y].max(); }
    /** Return leftmost coordinate of the rectangle (+X is to the right). */
    C left() const { return f[X].min(); }
    /** Return rightmost coordinate of the rectangle (+X is to the right). */
    C right() const { return f[X].max(); }

    /** Get the horizontal extent of the rectangle. */
    C width() const { return f[X].extent(); }
    /** Get the vertical extent of the rectangle. */
    C height() const { return f[Y].extent(); }
    /** Get the ratio of width to height of the rectangle. */
    Coord aspectRatio() const { return (cast(Coord)width()) / (cast(Coord)height()); }

    /** Get rectangle's width and height as a point.
     * @return Point with X coordinate corresponding to the width and the Y coordinate
     *         corresponding to the height of the rectangle. */
    CPoint dimensions() const { return CPoint(f[X].extent(), f[Y].extent()); }
    /** Get the point in the geometric center of the rectangle. */
    CPoint midpoint() const { return CPoint(f[X].middle(), f[Y].middle()); }

    /** Compute rectangle's area. */
    C area() const { return f[X].extent() * f[Y].extent(); }
    /** Check whether the rectangle has zero area. */
    bool hasZeroArea() const { return (area() == 0); }
    /** Check whether the rectangle is effectively empty */
    bool isEmpty() const { return f[X].isEmpty() || f[Y].isEmpty(); }

    /** Get the larger extent (width or height) of the rectangle. */
    C maxExtent() const { return cast(C)math.fmax(f[X].extent(), f[Y].extent()); }
    /** Get the smaller extent (width or height) of the rectangle. */
    C minExtent() const { return cast(C)math.fmin(f[X].extent(), f[Y].extent()); }

    /+ Test other rectangles and points for inclusion. +/

    /** Check whether the rectangles have any common points.
     * Empty rectangles will not intersect with any other rectangle. */
    bool intersects(in GenericRect!C r) const
    { return f[X].intersects(r[X]) && f[Y].intersects(r[Y]); }

    /** Check whether the rectangle includes all points in the given rectangle.
     * Empty rectangles will be contained in any non-empty rectangle. */
    bool contains(in GenericRect!C r) const
    { return f[X].contains(r[X]) && f[Y].contains(r[Y]); }

    /** Check whether the given point is within the rectangle. */
    bool contains(in CPoint p) const
    { return f[X].contains(p[X]) && f[Y].contains(p[Y]); }
    
    /+ Modify the rectangle +/

    /** Set the minimum X coordinate of the rectangle. */
    void setLeft(C val)
    { f[X].setMin(val); }

    /** Set the maximum X coordinate of the rectangle. */
    void setRight(C val)
    { f[X].setMax(val); }

    /** Set the minimum Y coordinate of the rectangle. */
    void setTop(C val)
    { f[Y].setMin(val); }

    /** Set the maximum Y coordinate of the rectangle. */
    void setBottom(C val)
    { f[Y].setMax(val); }

    /** Set the upper left point of the rectangle. */
    void setMin(in CPoint p)
    {
        f[X].setMin(p[X]);
        f[Y].setMin(p[Y]);
    }

    /** Set the lower right point of the rectangle. */
    void setMax(in CPoint p)
    {
        f[X].setMax(p[X]);
        f[Y].setMax(p[Y]);
    }

    /** Enlarge the rectangle to contain the given point. */
    void expandTo(in CPoint p)
    { 
        f[X].expandTo(p[X]);
        f[Y].expandTo(p[Y]);
    }

    /** Enlarge the rectangle to contain the argument.
     * Unioning with an empty rectangle results in no changes. */
    void unionWith(in GenericRect!C b)
    {
        f[X].unionWith(b[X]);
        f[Y].unionWith(b[Y]);
    }
    
    static GenericRect!C intersect(in GenericRect!C lhs, in GenericRect!C rhs)
    {
        CInterval x = lhs.f[X] & rhs.f[X], y = lhs.f[Y] & rhs.f[Y];
        return GenericRect!C(x, y);
    }

    /** Expand the rectangle in both directions by the specified amount.
     * Note that this is different from scaling. Negative values wil shrink the
     * rectangle. If <code>-amount</code> is larger than
     * half of the width, the X interval will contain only the X coordinate
     * of the midpoint; same for height. */
    void expandBy(C amount)
    { expandBy(amount, amount); }

    /** Expand the rectangle in both directions.
     * Note that this is different from scaling. Negative values wil shrink the
     * rectangle. If <code>-x</code> is larger than
     * half of the width, the X interval will contain only the X coordinate
     * of the midpoint; same for height. */
    void expandBy(C x, C y)
    { 
        f[X].expandBy(x);
        f[Y].expandBy(y);
    }

    /** Expand the rectangle by the coordinates of the given point.
     * This will expand the width by the X coordinate of the point in both directions
     * and the height by Y coordinate of the point. Negative coordinate values will
     * shrink the rectangle. If <code>-p[X]</code> is larger than half of the width,
     * the X interval will contain only the X coordinate of the midpoint;
     * same for height. */
    void expandBy(in CPoint p)
    { expandBy(p[X], p[Y]); }

    /+ Operators +/

    GenericRect!C opBinary(string op, T)(T rhs) const
    {
        /* Offset the rectangle by a vector. */
        static if (op == "+") { return GenericRect!C(this.f[X] + rhs[X], this.f[Y] + rhs[Y]); }
        else static if (op == "-") { return GenericRect!C(this.f[X] - rhs[X], this.f[Y] - rhs[Y]); }
        else static if (op == "|") {
            GenericRect!C r = this;
            r.unionWith(rhs);
            return r;
        }
        else static if (op == "&") { return intersect(this, rhs); }
        else static assert(false, "GenericRect!"~C.stringof~" operator "~op~" not implemented");
    }
    
    GenericRect!C opBinary(string op, T : Affine)(T rhs) const if (op == "*")
    {
        Point[4] pts;
        foreach (ref p; pts) p = this.corner(i) * rhs;
        Coord minx = math.min(math.min(pts[0][X], pts[1][X]), math.min(pts[2][X], pts[3][X]));
        Coord miny = math.min(math.min(pts[0][Y], pts[1][Y]), math.min(pts[2][Y], pts[3][Y]));
        Coord maxx = math.max(math.max(pts[0][X], pts[1][X]), math.max(pts[2][X], pts[3][X]));
        Coord maxy = math.max(math.max(pts[0][Y], pts[1][Y]), math.max(pts[2][Y], pts[3][Y]));
        f[X].setMin(minx); f[X].setMax(maxx);
        f[Y].setMin(miny); f[Y].setMax(maxy);
    }

    void opOpAssign(string op, T)(T rhs)
    { mixin("this = this "~op~" rhs;"); }

    protected CInterval[2] f;
}

Coord distanceSq(in Point p, in Rect rect)
{
    Coord dx = 0, dy = 0;
    if ( p[X] < rect.left() ) {
        dx = p[X] - rect.left();
    } else if ( p[X] > rect.right() ) {
        dx = rect.right() - p[X];
    }
    if (p[Y] < rect.top() ) {
        dy = rect.top() - p[Y];
    } else if (  p[Y] > rect.bottom() ) {
        dy = p[Y] - rect.bottom();
    }
    return dx*dx+dy*dy;
}

Coord distance(in Point p, in Rect rect)
{
    // copy of distanceSq, because we need to use hypot()
    Coord dx = 0, dy = 0;
    if ( p[X] < rect.left() ) {
        dx = p[X] - rect.left();
    } else if ( p[X] > rect.right() ) {
        dx = rect.right() - p[X];
    }
    if (p[Y] < rect.top() ) {
        dy = rect.top() - p[Y];
    } else if (  p[Y] > rect.bottom() ) {
        dy = p[Y] - rect.bottom();
    }
    return math.hypot(dx, dy);
}

/** Return the smallest integer rectangle which contains this one. */
IntRect roundOutwards(in Rect r)
{ return IntRect(geom.interval.roundOutwards(r[X]), geom.interval.roundOutwards(r[Y])); }
/** Return the largest integer rectangle which is contained in this one. */
IntRect roundInwards(in Rect r)
{ return IntRect(geom.interval.roundInwards(r[X]), geom.interval.roundInwards(r[Y])); }

alias Rect = GenericRect!Coord;
alias IntRect = GenericRect!IntCoord;

unittest
{
    IntRect a = IntRect(0, 0, 10, 10), a2 = a, b = IntRect(-5, -5, 5, 5);
    IntRect empty = IntRect.empty();

    assert(a == a);
    assert(a == a2);
    assert(empty == empty);
    assert(a != empty);
    assert(a != b);

    auto c = IntRect(-10, -10, -1, -1);
    auto d = IntRect(1, 1, 9, 9);

    assert(a.intersects(a));
    assert(a.intersects(b));
    assert(b.intersects(a));
    assert(b.intersects(c));
    assert(c.intersects(b));
    assert(a.intersects(d));
    assert(d.intersects(a));
    assert(!a.intersects(c));
    assert(!c.intersects(a));
    assert(!c.intersects(d));
    assert(!empty.intersects(empty));
    assert(!empty.intersects(a));

    /* JonCruz failure: (10, 20)-(55,30) and (45,20)-(100,30) should intersect. */
    a = IntRect(10, 20, 55, 30);
    b = IntRect(45, 20, 100,30);

    assert(a.intersects(a));
    assert(a.intersects(b));
    assert(b.intersects(a));

    /* Intersection operator */
    a = IntRect(0, 0, 10, 10);
    b = IntRect(-5, -5, 5, 5);
    c = IntRect(-10, -10, -1, -1);
    d = IntRect(1, 1, 9, 9);
    
    auto int_ab = IntRect(0, 0, 5, 5);
    auto int_bc = IntRect(-5, -5, -1, -1);

    assert((a & a) == a);
    assert((a & b) == int_ab);
    assert((b & c) == int_bc);
    assert((a & c) == a.empty());
    assert((a & d) == d);
    assert((a & empty) == a.empty());
    assert((empty & empty) == empty);

    a &= b;
    assert(a == int_ab);
    a &= empty;
    assert(a == a.empty());
    
    /* Union operator */
    a = IntRect(0, 0, 10, 10);
    b = IntRect(-5, -5, 5, 5);

    IntRect old_a = a;
    IntRect uni_ab = IntRect(-5, -5, 10, 10), uni_bc = IntRect(-10, -10, 5, 5);

    assert((a | b) == uni_ab);
    assert((b | c) == uni_bc);
    assert((a | a) == a);
    assert((a | d) == a);
    assert((a | int_ab) == a);
    assert((b | int_ab) == b);
    assert((uni_ab | a) == uni_ab);
    assert((uni_bc | c) == uni_bc);
    assert((a | empty) == a);
    assert((empty | empty) == empty);

    a |= b;
    assert(a == uni_ab);
    a = old_a;
    a |= empty;
    assert(a == old_a);
    
    /* Area */
    auto zero = IntRect(0,0,0,0);

    assert(a.area() == 100);
    assert(a.area() == a.width() * a.height());
    assert(b.area() == 100);
    assert(c.area() == 81);
    assert(d.area() == 64);
    assert(!a.hasZeroArea());
    assert(zero.hasZeroArea());

    /* Dimensions */
    a = IntRect(-10, -20, 10, 20);
    b = IntRect(-15, 30, 45, 90);

    assert(a.width() == 20);
    assert(a.height() == 40);
    assert(a.left() == -10);
    assert(a.top() == -20);
    assert(a.right() == 10);
    assert(a.bottom() == 20);
    assert(a.min() == a.CPoint(-10, -20));
    assert(a.max() == a.CPoint(10, 20));
    assert(a.minExtent() == a.width());
    assert(a.maxExtent() == a.height());
    assert(a.dimensions() == a.CPoint(20, 40));
    assert(a.midpoint() == a.CPoint(0, 0));

    assert(b.width() == 60);
    assert(b.height() == 60);
    assert(b.left() == -15);
    assert(b.top() == 30);
    assert(b.right() == 45);
    assert(b.bottom() == 90);
    assert(b.min() == b.CPoint(-15, 30));
    assert(b.max() == b.CPoint(45, 90));
    assert(b.minExtent() == b.maxExtent());
    assert(b.dimensions() == b.CPoint(60, 60));
    assert(b.midpoint() == b.CPoint(15, 60));
    
    /* Modification */
    a = IntRect(-1, -1, 1, 1);
    a.expandBy(9);
    assert(a == IntRect(-10, -10, 10, 10));
    a.setMin(a.CPoint(0, 0));
    assert(a == IntRect(0, 0, 10, 10));
    a.setMax(a.CPoint(20, 30));
    assert(a == IntRect(0, 0, 20, 30));
    a.setMax(a.CPoint(-5, -5));
    assert(a == IntRect(-5, -5, -5, -5));
    a.expandTo(a.CPoint(5, 5));
    assert(a == IntRect(-5, -5, 5, 5));
    a.expandTo(a.CPoint(0, 0));
    assert(a == IntRect(-5, -5, 5, 5));
    a.expandTo(a.CPoint(0, 15));
    assert(a == IntRect(-5, -5, 5, 15));
    a.expandBy(-10);
    assert(a == IntRect(0, 5, 0, 5));
    assert(a.midpoint() == a.CPoint(0, 5));
    a.unionWith(IntRect(-20, 0, -10, 20));
    assert(a == IntRect(-20, 0, 0, 20));
    
    /* Offset */
    a = IntRect(0, 0, 5, 5);
    old_a = a;
    auto app1 = IntRect(-5, 0, 0, 5);
    auto amp1 = IntRect(5, 0, 10, 5);
    auto app2 = IntRect(5, -10, 10, -5);
    auto amp2 = IntRect(-5, 10, 0, 15);
    auto p1 = a.CPoint(-5, 0), p2 = a.CPoint(5, -10);

    assert(a + p1 == app1);
    assert(a + p2 == app2);
    assert(a - p1 == amp1);
    assert(a - p2 == amp2);

    a += p1;
    assert(a == app1);
    a = old_a;
    a += p2;
    assert(a == app2);
    a = old_a;
    a -= p1;
    assert(a == amp1);
    a = old_a;
    a -= p2;
    assert(a == amp2);
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
