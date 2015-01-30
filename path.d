/*
 * Path - a sequence of contiguous curves
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

module geom.path;

public import geom.coord;

import geom.affine;
import geom.bezier_curve;
import geom.curve;
import geom.d2;
import geom.interval;
import geom.piecewise;
import geom.point;
import geom.rect;
import geom.sbasis;

/** PathPosition (generalized time value) in the path.
 *
 * This class exists because mapping the range of multiple curves onto the same interval
 * as the curve index, we lose some precision. For instance, a path with 16 curves will
 * have 4 bits less precision than a path with 1 curve. If you need high precision results
 * in long paths, either use this class and related methods instead of the standard methods
 * pointAt(), nearestTime() and so on, or use curveAt() to first obtain the curve, then
 * call the method again to obtain a high precision result.
 */
struct PathPosition
{
    Coord t; /// Time value in the curve
    size_t curve_index; /// Index of the curve in the path
    this(size_t idx, Coord tval) { t = tval; curve_index = idx; }
};

enum Stitching {
    NO_STITCHING = 0,
    STITCH_DISCONTINUOUS
};

/** Sequence of contiguous curves, aka spline.
 *
 * Path represents a sequence of contiguous curves, also known as a spline.
 * It corresponds to a "subpath" in SVG terminology. It can represent both
 * open and closed subpaths. The final point of each curve is exactly
 * equal to the initial point of the next curve.
 *
 * The path always contains a linear closing segment that connects
 * the final point of the last "real" curve to the initial point of the
 * first curve. This way the curves form a closed loop even for open paths.
 * If the closing segment has nonzero length and the path is closed, it is
 * considered a normal part of the path data.
 *
 * Since the methods for inserting, erasing and replacing curves can cause a path
 * to become non-contiguous, they have an additional parameter, called @a stitching,
 * which determines whether non-contiguous segments are joined with additional
 * linear segments that fill the created gaps.
 *
 * Note that this class cannot represent arbitrary shapes, which may contain holes.
 * To do that, use PathVector, which is more generic.
 *
 * It's not very convenient to create a Path directly. To construct paths more easily,
 * use PathBuilder.
 */
class Path
{
    alias ClosingSegment = BezierCurveN!1;
    alias StitchSegment = BezierCurveN!1;

    /** Construct an empty path starting at the specified point. */
    this(Point p = Point())
    {
        _closing_seg = new ClosingSegment(p, p);
        _closed = false;
        _curves ~= _closing_seg;
    }

    this(in Path o)
    {
        import object : reserve;

        _closed = o._closed;
        _curves.reserve(o._curves.length);

        foreach (t; o._curves[0 .. $-1])
            _curves ~= t.duplicate();

        _closing_seg = o._closing_seg.duplicate();
        _curves ~= _closing_seg;
    }

    /** Swap contents with another path */
    nothrow void swap(ref Path other)
    {
        import std.algorithm : swap;
        swap(other._curves, _curves);
        swap(other._closing_seg, _closing_seg);
        swap(other._closed, _closed);
    }

    /** Swap contents of two paths. */
    static nothrow void swap(ref Path a, ref Path b) { a.swap(b); }

    /** Access a curve by index */
    ref const(Curve) opIndex(size_t i) const { return _curves[i]; }
    /** Access a curve by index */
    ref const(Curve) at(size_t i) const { return _curves[i]; }

    /** Access the first curve in the path.
     * Since the curve always contains at least a degenerate closing segment,
     * it is always safe to use this method. */
    ref const(Curve) front() const { return _curves[0]; }

    /** Access the last curve in the path.
     * Since the curve always contains at least a degenerate closing segment,
     * it is always safe to use this method. */
    ref const(Curve) back() const { return back_default(); }
    ref const(Curve) back_open() const { return _curves[$-2]; }
    ref const(Curve) back_closed() const { return _curves[$-1]; }
    ref const(Curve) back_default() const { return (_closed ? back_closed() : back_open()); }

    size_t size_open() const { return _curves.length-1; }
    size_t size_closed() const { return _curves.length; }
    size_t size_default() const { return _includesClosingSegment() ? size_closed() : size_open(); }
    size_t size() const { return size_default(); }

    /** Check whether path is empty.
     * The path is empty if it contains only the closing segment, which according
     * to the continuity invariant must be degenerate. Note that unlike standard
     * containers, two empty paths are not necessarily identical, because the
     * degenerate closing segment may be at a different point, affecting the operation
     * of methods such as appendNew(). */
    bool empty() const { return (_curves.length == 1); }
    /** Check whether the path is closed. */
    bool closed() const { return _closed; }
    /** Set whether the path is closed. */
    void close(bool closed = true) { _closed = closed; }

    /** Remove all curves from the path.
     * The initial and final points of the closing segment are set to (0,0). */
    void clear()
    {
        _closing_seg.setInitial(Point(0, 0));
        _closing_seg.setFinal(Point(0, 0));
        _curves = [_closing_seg];
    }

    /** Get the approximate bounding box.
     * The rectangle returned by this method will contain all the curves, but it's not
     * guaranteed to be the smallest possible one */
    Rect boundsFast() const
    {
        Rect bounds = Rect.empty();
        if (empty())
            return bounds;
        bounds = _curves[0].boundsFast();

        // the closing path segment can be ignored, because it will always lie within the bbox of the rest of the path
        // note that if we are !empty() there is always more than one element in _curves
        foreach (i; _curves[1 .. $-1]) {
            bounds.unionWith(i.boundsFast());
        }

        return bounds;
    }

    /** Get a tight-fitting bounding box.
     * This will return the smallest possible axis-aligned rectangle containing
     * all the curves in the path. */
    Rect boundsExact() const
    {
        Rect bounds = Rect.empty();
        if (empty())
            return bounds;
        bounds = _curves[0].boundsExact();

        // the closing path segment can be ignored, because it will always lie within the bbox of the rest of the path
        // note that if we are !empty() there is always more than one element in _curves
        foreach (i; _curves[1 .. $-1]) {
            bounds.unionWith(i.boundsFast());
        }
        return bounds;
    }

    Piecewise!(D2!SBasis) toPwSb() const
    {
        Piecewise!(D2!SBasis) ret;
        ret.push_cut(0);
        uint i = 1;
        bool degenerate = true;
        // pw!(d2!()) is always open. so if path is closed, add closing segment as well to pwd2.
        foreach (it; _curves[0 .. (_closed ? $-1 : $)]) {
            if (!it.isDegenerate()) {
                ret.push(it.toSBasis(), i++);
                degenerate = false;
            }
        }
        if (degenerate) {
            // if path only contains degenerate curves, no second cut is added
            // so we need to create at least one segment manually
            ret = Piecewise!(D2!SBasis)(initialPoint());
        }
        return ret;
    }

    Path opBinary(string op, T : Affine)(in T m) const if (op == "*")
    {
        Path ret = new Path(this);
        foreach (ref c; ret._curves) c.transform(m);
        checkContinuity(); // Affines cannot break spline continuity, why this?
        return ret;
    }
    
    void opOpAssign(string op, T : Affine)(in T m) if (op == "*")
    {
        foreach (ref c; _curves) c.transform(m);
        checkContinuity();
    }

    /** Get values for which pointAt() and valueAt() yield valid results. */
    Interval timeRange() const { return Interval(0, size_default()); }

    /** Get the curve at the specified time value. */
    ref const(Curve) curveAt(Coord t) const { return at(_getPosition(t).curve_index); }
    ref const(Curve) curveAt(in PathPosition pos) const { return at(pos.curve_index); }

    /** Get the point at the specified time value.
     * Note that this method has reduced precision with respect to calling pointAt()
     * directly on the curve. If you want high precision results, use the version
     * that takes a PathPosition parameter.
     * 
     * Allowed time values range from zero to the number of curves; you can retrieve
     * the allowed range of values with timeRange(). */
    Point pointAt(Coord t) const { return pointAt(_getPosition(t)); }
    Point pointAt(in PathPosition pos) const { return at(pos.curve_index).pointAt(pos.t); }

    /** Get one coordinate (X or Y) at the specified time value. */
    Coord valueAt(Coord t, size_t d) const { return valueAt(_getPosition(t), d); }
    Coord valueAt(in PathPosition pos, size_t d) const { return at(pos.curve_index).valueAt(pos.t, d); }

    Point opCall(Coord t) const { return pointAt(t); }

    Coord[] roots(Coord v, size_t d) const
    {
        Coord[] res;
        for (size_t i = 0; i <= size(); i++) {
            Coord[] temp = _curves[i].roots(v, d);
            foreach (j; temp) {
                res ~= j + i; // ?
            }
        }
        return res;
    }

    /** Determine the winding number at the specified point.
     * 
     * The winding number is the number of full turns made by a ray that connects the passed
     * point and the path's value (i.e. the result of the pointAt() method) as the time increases
     * from 0 to the maximum valid value. Positive numbers indicate turns in the direction
     * of increasing angles.
     *
     * Winding numbers are often used as the definition of what is considered "inside"
     * the shape. Typically points with either nonzero winding or odd winding are
     * considered to be inside the path. */
    int winding(in Point p) const
    {
        int wind = 0;

        /* To handle all the edge cases, we consider the maximum Y edge of the bounding box
         * as not included in box. This way paths that contain linear horizontal
         * segments will be treated correctly. */
        foreach (i; _curves) {
            Rect bounds = i.boundsFast();

            if (bounds.height() == 0) continue;
            if (p[X] > bounds.right() || !bounds[Y].lowerContains(p[Y])) {
                // Ray doesn't intersect bbox, so we ignore this segment
                continue;
            }

            if (p[X] < bounds.left()) {
                /* Ray intersects the curve's bbox, but the point is outside it.
                 * The winding contribution is exactly the same as that
                 * of a linear segment with the same initial and final points. */
                Point ip = i.initialPoint();
                Point fp = i.finalPoint();
                Rect eqbox = Rect(ip, fp);

                if (eqbox[Y].lowerContains(p[Y])) {
                    /* The ray intersects the equivalent linear segment.
                     * Determine winding contribution based on its derivative. */
                    if (ip[Y] < fp[Y]) {
                        ++wind;
                    } else if (ip[Y] > fp[Y]) {
                        --wind;
                    } else {
                        // should never happen, because bounds.height() was not zero
                        assert(false);
                    }
                }
            } else {
                // point is inside bbox
                wind += i.winding(p);
            }
        }
        return wind;
    }

    Coord nearestTime(in Point p, Coord *dist = null) const
    {
        PathPosition pos = nearestPosition(p, dist);
        return pos.curve_index + pos.t;
    }

    /** Returns the nearest time for each curve in this path. */
    Coord[] nearestTimePerCurve(in Point p) const
    {
        Coord[] np;
        for (size_t i = 0; i < size_default(); ++i) {
            np ~= at(i).nearestPoint(p);
        }
        return np;
    }

    PathPosition nearestPosition(in Point p, Coord *dist = null) const
    {
        Coord mindist = Coord.max;
        PathPosition ret;

        if (_curves.length == 1) {
            // naked moveto
            ret.curve_index = 0;
            ret.t = 0;
            if (dist) *dist = distance(_closing_seg.initialPoint(), p);
            return ret;
        }

        for (size_t i = 0; i < size_default(); ++i) {
            const(Curve) c = at(i);
            if (distance(p, c.boundsFast()) >= mindist) continue;

            Coord t = c.nearestPoint(p);
            Coord d = distance(c.pointAt(t), p);
            if (d < mindist) {
                mindist = d;
                ret.curve_index = i;
                ret.t = t;
            }
        }
        if (dist) *dist = mindist;

        return ret;
    }

    void appendPortionTo(ref Path ret, Coord from, Coord to) const
    {
        import std.math : modf;

        if (!(from >= 0 && to >= 0))
            throw new Exception("from and to must be >=0 in Path.appendPortionTo");
        if (to == 0)
            to = size() + 0.999999;
        if (from == to)
            return;

        real fi, ti;
        Coord ff = modf(from, fi), tf = modf(to, ti);
        if (tf == 0) {
            ti--;
            tf = 1;
        }

        size_t fromi = cast(size_t)fi;
        size_t toi = cast(size_t)ti;
        if (fi == ti && from < to) {
            // we've narrowed it down to a portion of one curve
            const Curve v = _curves[fromi].portion(ff, tf);
            ret.append(v, Stitching.STITCH_DISCONTINUOUS);
            return;
        }

        if (ff != 1.) {
            const Curve fromv = _curves[fromi].portion(ff, 1.);
            ret.append(fromv, Stitching.STITCH_DISCONTINUOUS);
        }
        if (from >= to) {
            size_t ender = _curves.length-1;
            if (_curves[$-1].isDegenerate())
                ender--;

            ret.insert(ret.size(), _curves[++fromi .. ender], Stitching.STITCH_DISCONTINUOUS);
            ret.insert(ret.size(), _curves[0 .. toi], Stitching.STITCH_DISCONTINUOUS);
        } else {
            ret.insert(ret.size(), _curves[++fromi .. toi], Stitching.STITCH_DISCONTINUOUS);
        }

        Curve tov = _curves[toi].portion(0., tf);
        ret.append(tov, Stitching.STITCH_DISCONTINUOUS);
    }

    /** Get a subset of the current path.
     * Note that @a f can be smaller than @a t, in which case the returned part of the path
     * will go in the opposite direction.
     * @param f Time value specifying the initial point of the returned path
     * @param t Time value specifying the final point of the returned path
     * @return Portion of the path */
    Path portion(Coord f, Coord t) const
    {
        Path ret;
        ret.close(false);
        appendPortionTo(ret, f, t);
        return ret;
    }

    /** Get a subset of the current path.
     * This version takes an Interval. */
    Path portion(Interval i) const { return portion(i.min(), i.max()); }

    /** Obtain a reversed version of the current path.
     * The final point of the current path will become the initial point
     * of the reversed path, unless it is closed and has a non-degenerate
     * closing segment. In that case, the new initial point will be the final point
     * of the last "real" segment. */
    Path reversed() const
    {
        import std.range : retro;
        Path ret;
        ret._curves = []; // clear the array
        foreach (c; retro(_curves)) {
            ret._curves ~= c.reverse();
        }
        ret._closing_seg = ret._closing_seg.reverse();
        ret._curves ~= ret._closing_seg;
        return ret;
    }

    void insert(size_t pos, in Curve curve, Stitching stitching = Stitching.NO_STITCHING) { insert(pos, [curve], stitching); }

    void insert(size_t pos, in Curve[] curve, Stitching stitching = Stitching.NO_STITCHING)
    {
        import object : reserve;
        import core.exception;
        if (pos >= _curves.length) throw new RangeError;

        Curve[] source;
        source.reserve(curve.length);

        // Manually copy the input array
        foreach (t; curve) {
            source ~= t.duplicate();
        }

        if (stitching)
            stitch(pos, pos, source);

        // Slice the array in half at the insertion point.
        Curve[] before = _curves[0 .. pos];
        Curve[] after = _curves[pos .. $];
        _curves = before ~ source ~ after;

        assert(_curves[$-1] == _closing_seg);
        do_update();
    }

    void erase(size_t pos, Stitching stitching = Stitching.NO_STITCHING) { erase(pos, pos, stitching); }

    void erase(size_t first, size_t last, Stitching stitching = Stitching.NO_STITCHING)
    {
        import core.exception;
        import std.stdio : writeln;
        if (last < first || last >= _curves.length) throw new RangeError;

        // Slice the array in half at the erasure point.
        Curve[] before = _curves[0 .. first];
        Curve[] after = _curves[last .. $];

        // Optionally stitch at the point of erasure to prevent discontinuity.
        if (stitching) {
            Curve[] fixup;
            stitch(first, last, fixup);
            _curves = before ~ fixup ~ after;
        } else {
            _curves = before ~ after;
        }

        assert(_curves[$-1] == _closing_seg);
        do_update();
    }

    /** erase last segment of path */
    void erase_last() { erase(size() - 1); }

    /** Get the first point in the path. */
    Point initialPoint() const { return _closing_seg[1]; }
    /** Get the last point in the path.
     * If the path is closed, this is always the same as the initial point. */
    Point finalPoint() const { return _closing_seg[_closed ? 1 : 0]; }

    void append(in Curve curve, Stitching stitching = Stitching.NO_STITCHING)
    {
        if (stitching)
            stitchTo(curve.initialPoint());
        do_append(curve.duplicate());
    }

    /** Append a stitching segment ending at the specified point. */
    void stitchTo(in Point p)
    {
        if (!empty() && finalPoint() != p) {
            do_append(new StitchSegment(finalPoint(), p));
        }
    }

    /** Append a new curve to the path.
     *
     * This method will automaticaly use the current final point of the path
     * as the first argument of the new curve's constructor. To call this method,
     * you'll need to write e.g.:
     * @code
       path.appendNew!CubicBezier(control1, control2, end_point);
       @endcode
     * It is important to note that the coordinates passed to appendNew should be finite!
     * If one of the coordinates is infinite, 2geom will throw a ContinuityError exception.
     */
    void appendNew(CurveType, A...)(A a) { do_append(new CurveType(finalPoint(), a)); }

    /** Verify the continuity invariant.
     * If the path is not contiguous, this will throw a CountinuityError. 
     *
     * If this function throws, the path is necessarily considered to be in an invalid state.
     * Attempting to operate further on paths in invalid states will cause undefined behavior.
     */
    void checkContinuity() const
    {
        if (_curves[0].initialPoint() != _curves[$-1].finalPoint()) {
            throw new ContinuityError;
        }

        for (size_t i = 0; i < _curves.length-1; ++i) {
            if (_curves[i].finalPoint() != _curves[i+1].initialPoint())
                throw new ContinuityError;
        }
    }

private:

    bool _includesClosingSegment() const { return _closed && !_closing_seg.isDegenerate(); }

    // n.b. takes ownership of curve object
    void do_append(Curve c)
    {
        if (_curves[0] == _closing_seg) {
            _closing_seg.setFinal(c.initialPoint());
        } else {
            if (c.initialPoint() != finalPoint()) {
                throw new ContinuityError;
            }
        }
        _curves = _curves[0 .. $-1] ~ c ~ _closing_seg;
        _closing_seg.setInitial(c.finalPoint());
    }

    PathPosition _getPosition(Coord t) const
    {
        import std.math : modf;

        size_t sz = size_default();
        if (t < 0 || t > sz)
            throw new ContinuityError;

        PathPosition ret;
        real k;
        ret.t = modf(t, k);
        ret.curve_index = cast(size_t)k; // real <-!-> ulong
        if (ret.curve_index == sz) {
            --ret.curve_index;
            ret.t = 1;
        }
        return ret;
    }

    void stitch(size_t first_replaced, size_t last_replaced, ref Curve[] source)
    {
        const Curve first = _curves[first_replaced];
        const Curve last = _curves[last_replaced];

        // if the array contains a path...
        if (source.length != 0) {

            // don't bother stitching inserts at the start
            if (first_replaced != 0) {
                if (first.finalPoint() != source[0].initialPoint()) {
                    Curve stitch = new StitchSegment(first.initialPoint(), source[0].initialPoint());
                    source = stitch ~ source;
                }
            }

            // don't bother stitching inserts at the end
            if (last_replaced != _curves.length - 1) {
                if (last.finalPoint() != source[$-1].finalPoint()) {
                    Curve stitch = new StitchSegment(source[$-1].finalPoint(), last.initialPoint());
                    source ~= stitch;
                }
            }
        } else if (first != last && first_replaced != 0 && last_replaced != _curves.length - 1) {
            // the input array contains no path, fill it with a segment stitching first to last
            if (first.initialPoint() != _curves[last_replaced - 1].finalPoint()) {
                Curve stitch = new StitchSegment(_curves[last_replaced - 1].finalPoint(), first.initialPoint());
                source = stitch ~ source;
            }
        }
    }
    
    void do_update()
    {
        if (_curves[0] != _closing_seg) {
            _closing_seg.setPoint(0, back().finalPoint());
            _closing_seg.setPoint(1, front().initialPoint());
        }

        checkContinuity();
    }

    Curve[] _curves;
    ClosingSegment _closing_seg;
    bool _closed;
}
unittest
{
    import geom.transforms;

    // you better believe we're gonna unittest this sucker

    auto a = bezierFromPoints(Point(4,0), Point(1,0), Point(2, 2));
    auto b = bezierFromPoints(Point(6,0), Point(3,-1), Point(0, 7));
    auto c = bezierFromPoints(Point(5,2), Point(9,0), Point(8,3));

    Path p = new Path; // default path, open, contains closing segment
    p.insert(0, a);
    
    // test stitching
    p.insert(0, b, Stitching.STITCH_DISCONTINUOUS);
    p.insert(1, c, Stitching.STITCH_DISCONTINUOUS);

    auto d = bezierFromPoints(Point(10,0), Point(11,0), Point(12,0));
    p.insert(p.size(), d, Stitching.STITCH_DISCONTINUOUS);

    // test continuous inserts
    auto e = bezierFromPoints(Point(12,0), Point(13,1), Point(14,0));
    auto f = bezierFromPoints(Point(14,0), Point(15,-1), Point(16,0));
    p.insert(p.size(), e);
    p.insert(p.size(), f);

    // appending
    p.appendNew!(BezierCurveN!3)(Point(17,0), Point(18,1), Point(19,-3));

    // erasing
    p.erase(0);
    p.erase_last();
    p.erase(p.size()-2, Stitching.STITCH_DISCONTINUOUS);

    // transformations
    Affine m = Translate(1, 0);
    Point ip = p.initialPoint();
    Point fp = p.finalPoint();
    
    p *= m;

    assert(p.initialPoint() == ip * m);
    assert(p.finalPoint() == fp * m);

    // make sure that path knows to stay continuous when not stitching
    try {
        p.insert(2, b);
        assert(false, "Failed to throw non-continuity exception");
    } catch (ContinuityError e) {}
}


class ContinuityError : Exception
{
    @safe pure nothrow this(Throwable next = null)
    { super("non-contiguous path", next); }
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
