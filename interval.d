/**
 * Closed interval
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

module geom.interval;
public import geom.coord;

import math = std.math;

/**
 * A range of numbers which is never empty.
 * @ingroup Primitives
 */
struct GenericInterval(C)
{
    /+ Create intervals. +/

    @disable this();

    /** Create an interval from another interval. */
    this(in GenericInterval!C o) { _b = o._b; }
    /** Create an interval that contains a single point. */
    this(C u) { _b = [u, u]; }
    /** Create an interval that contains all points between @c u and @c v. */
    this(C u, C v)
    {
        if (u <= v) {
            _b[0] = u; _b[1] = v;
        } else {
            _b[0] = v; _b[1] = u;
        }
    }

    /** Create an empty interval. */
    static GenericInterval!C empty() { return GenericInterval!C(cast(C)0, cast(C)0); }

    /+ Inspect endpoints +/

    C min() const { return _b[0]; }
    C max() const { return _b[1]; }
    C extent() const { return max() - min(); }
    C middle() const { return (max() + min()) / 2; }
    bool isSingular() const { return min() == max(); }
    bool isEmpty() const { return this == empty(); }

    /+ Test coordinates and other intervals for inclusion. +/

    /** Check whether the interval includes this number. */
    bool contains(C val) const
    { return min() <= val && val <= max(); }
    /** Check whether the interval includes the given interval. */
    bool contains(in GenericInterval!C val) const
    { return min() <= val.min() && val.max() <= max(); }
    /** Check whether the intervals have any common elements. */
    bool intersects(in GenericInterval!C val) const
    { return contains(val.min()) || contains(val.max()) || val.contains(this); }

    /** Check whether the interior of the interval includes this number.
     * Interior means all numbers in the interval except its ends. */
    bool interiorContains(C val) const { return min() < val && val < max(); }
    /** @brief Check whether the interior of the interval includes the given interval.
     * Interior means all numbers in the interval except its ends. */
    bool interiorContains(in GenericInterval!C val) const { return min() < val.min() && val.max() < max(); }
    /** @brief Check whether the interiors of the intervals have any common elements.  A single point in common is not considered an intersection. */
    bool interiorIntersects(in GenericInterval!C val) const { return math.fmax(min(), val.min()) < math.fmin(max(), val.max()); }
    
    
    /+ Modify the interval +/

    // TODO: NaN handleage for the next two?
    /** Set the lower boundary of the interval.
     * When the given number is larger than the interval's largest element,
     * it will be reduced to the single number @c val. */
    void setMin(C val)
    {
        if(val > _b[1]) {
            _b[0] = _b[1] = val;
        } else {
            _b[0] = val;
        }
    }
    /** Set the upper boundary of the interval.
     * When the given number is smaller than the interval's smallest element,
     * it will be reduced to the single number @c val. */
    void setMax(C val)
    {
        if(val < _b[0]) {
            _b[1] = _b[0] = val;
        } else {
            _b[1] = val;
        }
    }
    /** Extend the interval to include the given number. */
    void expandTo(C val)
    {
       if(val < _b[0]) _b[0] = val;
       if(val > _b[1]) _b[1] = val;  // no else, as we want to handle NaN
    }

    /** Expand or shrink the interval in both directions by the given amount.
     * After this method, the interval's length (extent) will be increased by
     * <code>amount * 2</code>. Negative values can be given; they will shrink the interval.
     * Shrinking by a value larger than half the interval's length will create a degenerate
     * interval containing only the midpoint of the original. */
    void expandBy(C amount)
    {
        _b[0] -= amount;
        _b[1] += amount;
        if (_b[0] > _b[1]) {
            C halfway = (_b[0]+_b[1])/2;
            _b[0] = _b[1] = halfway;
        }
    }

    /** Union the interval with another one.
     * The resulting interval will contain all points of both intervals.
     * It might also contain some points which didn't belong to either - this happens
     * when the intervals did not have any common elements. */
    void unionWith(in GenericInterval!C a)
    {
        if(a._b[0] < _b[0]) _b[0] = a._b[0];
        if(a._b[1] > _b[1]) _b[1] = a._b[1];
    }

    static GenericInterval!C intersect(in GenericInterval!C t, in GenericInterval!C o)
    {
        C u = cast(C)math.fmax(t.min(), o.min());
        C v = cast(C)math.fmin(t.max(), o.max());
        if (u <= v) {
            return GenericInterval!C(u, v);
        }
        return empty();
    }

    GenericInterval!C opBinary(string op, T)(T rhs) const
    {
        static if (op == "+") { return GenericInterval!C(_b[0] + rhs, _b[1] + rhs); }
        else static if (op == "-") { return GenericInterval!C(_b[0] - rhs, _b[1] - rhs); }
        else static if (op == "*") { return GenericInterval!C(_b[0] * rhs, _b[1] * rhs); }
        else static if (op == "/") { return GenericInterval!C(_b[0] / rhs, _b[1] / rhs); } // TODO division by zero?
        else static if (op == "|") { 
            auto r = GenericInterval!C(this);
            r.unionWith(rhs);
            return r;
        }
        else static if (op == "&") { return intersect(this, rhs); }
        else static assert(false, "GenericInterval!"~C.stringof~" operator "~op~" not implemented");
    }

    void opOpAssign(string op, T)(T rhs)
    { mixin("this = this "~op~" rhs;"); }

    protected C[2] _b;
}

alias IntInterval = GenericInterval!IntCoord;
alias Interval = GenericInterval!Coord;

/** Return the smallest integer interval which contains this one. */
IntInterval roundOutwards(in Interval i)
{ return IntInterval(cast(IntCoord)math.floor(i.min()), cast(IntCoord)math.ceil(i.max())); }

/** Return the largest integer interval which is contained in this one. */
IntInterval roundInwards(in Interval i)
{
    IntCoord u = cast(IntCoord)math.ceil(i.min()), v = cast(IntCoord)math.floor(i.max());
    if (u > v) { return IntInterval.empty(); }
    return IntInterval(u, v);
}

unittest
{
    Interval i = Interval.empty();
    assert(i.min() == 0);
    assert(i.max() == 0);
    assert(i.extent() == 0);
    assert(i.middle() == 0);
    assert(i.isSingular());
    
    i = Interval(1, 4);
    Interval j = Interval(2, 3);
    
    assert(i.intersects(j));
    assert(j.intersects(i));
    
    Interval k = i & j;
    assert(k.min() == 2);
    assert(k.max() == 3);
    
    Interval l = i | j;
    assert(l.min() == 1);
    assert(l.max() == 4);

    k = Interval(1, 2);
    l = Interval(3, 4);
    j = k & l;
    assert(j.isEmpty());
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
