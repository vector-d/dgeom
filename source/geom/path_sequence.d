/**
 * PathSequence - a sequence of subpaths
 *
 * Authors:
 *   Johan Engelen <j.b.c.engelen@alumnus.utwente.nl>
 *   Krzysztof Kosiński <tweenk.pl@gmail.com>
 *   Liam P. White
 * 
 * Copyright 2008-2015 authors
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

module geom.path_sequence;

public import geom.coord;
import geom.interval;
import geom.path;
import geom.point;
import geom.rect;

/** Position (generalized time value) in the path sequence.
 *
 * This class exists because mapping the range of multiple curves onto the same interval
 * as the curve index, we lose some precision. For instance, a path with 16 curves will
 * have 4 bits less precision than a path with 1 curve. If you need high precision results
 * in long paths, use this class and related methods instead of the standard methods
 * pointAt(), nearestTime() and so on.
 */
struct PathSequencePosition
{
    size_t path_index = 0; /// Index of the path in the sequence
    size_t curve_index = 0; /// Index of the curve in the path
    Coord t = 0; /// Time value in the curve

    this(size_t _i, size_t _c, Coord _t) { curve_index = _c; t = _t; path_index = _i; }
}

/** Sequence of subpaths.
 *
 * This class corresponds to the SVG notion of a path:
 * a sequence of any number of open or closed contiguous subpaths.
 * Unlike Path, this class is closed under boolean operations.
 *
 * If you want to represent an arbitrary shape, this is the best class to use.
 * Shapes with a boundary that is composed of only a single contiguous
 * component can be represented with Path instead.
 */
class PathSequence
{
    this() {}
    this(in Path[] i) { foreach (p; i) _data ~= new Path(p); }
    this(in PathSequence o) { foreach (p; o._data) _data ~= new Path(p); }

    /// Check whether the sequence contains any paths.
    bool empty() const { return size() == 0; }

    /// Get the number of paths in the sequence.
    size_t size() const { return _data.length; }

    /// Get the total number of curves in the sequence.
    size_t curveCount() const
    {
        size_t n = 0;
        foreach (it; _data)
            n += it.size();
        return n;
    }

    ref inout(Path) opIndex(size_t i) inout { return _data[i]; }
    ref inout(Path) at(size_t i) inout { return _data[i]; }

    /// Append a path at the end.
    void push_back(in Path path) { _data ~= new Path(path); }

    /// Remove the last path.
    void pop_back() { _data = _data[0 .. $-1]; }

    /// Remove all paths from the sequence.
    void clear() { _data = []; }

    /** Change the number of paths.
     * If the array size increases, it is passed with paths that contain only
     * a degenerate closing segment at (0,0). */
    void resize(size_t n) { _data.length = n; }

    /** Reverse the direction of paths in the sequence.
     * @param reverse_paths If this is true, the order of paths is reversed as well;
     *                      otherwise each path is reversed, but their order in the
     *                      PathSequence stays the same */
    void reverse(bool reverse_paths = true)
    {
        if (reverse_paths) {
            import std.array : array;
            import std.range : retro;
            _data = retro(_data).array;
        }
        foreach (p; _data) p = p.reversed();
    }

    /** Get a new sequence with reversed direction of paths.
     * @param reverse_paths If this is true, the order of paths is reversed as well;
     *                      otherwise each path is reversed, but their order in the
     *                      PathSequence stays the same */
    PathSequence reversed(bool reverse_paths = true)
    {
        PathSequence ret = new PathSequence(this);
        ret.reverse(reverse_paths);
        return ret;
    }

    /// Get the range of allowed time values.
    Interval timeRange() const { return Interval(0, curveCount()); }

    /// Get the first point in the first path of the sequence. 
    Point initialPoint() const { return _data[0].initialPoint(); }

    /// Get the last point in the last path of the sequence.
    Point finalPoint() const { return _data[$-1].finalPoint(); }

    /** Determine the winding number at the specified point.
     * This is simply the sum of winding numbers for constituent paths. */
    int winding(in Point p) const
    {
        int wind = 0;
        foreach (i; _data)
            wind += i.winding(p);
        return wind;
    }

    /*Coord nearestTime(in Point p) const
    {
        XXX: this was never implemented in 2geom, design our own?
    }*/

    Rect boundsFast() const
    {
        Rect bound = Rect.empty();
        if (empty()) return bound;

        bound = _data[0].boundsFast();
        foreach (it; _data[1 .. $])
            bound.unionWith(it.boundsFast());
        return bound;
    }

    Rect boundsExact() const
    {
        Rect bound = Rect.empty();
        if (empty()) return bound;

        bound = _data[0].boundsExact();
        foreach (it; _data[1 .. $])
            bound.unionWith(it.boundsExact());
        return bound;
    }

    PathSequence opOpAssign(string op, T)(T b)
    {
        static if (op == "*") {
            if (!empty())
                foreach (ref it; _data)
                    it *= b;
            return this;
        } else static if (op == "~") { 
            _data ~= b;
            return this;
        }
        else static assert(false, "PathSequence operator "~op~" not implemented");
    }

    PathSequence opBinary(string op, T)(T b) const
    {
        auto ret = new PathSequence(this);
        mixin("ret "~op~"= b");
        return ret;
    }

private:

    PathSequencePosition _getPosition(Coord t) const
    {
        import std.math : modf;

        PathSequencePosition ret;
        real rest = 0;
        ret.t = modf(t, rest);
        ret.curve_index = cast(size_t)rest;
        for (; ret.path_index < size(); ++ret.path_index) {
            size_t s = _data[ret.path_index].size_default();
            if (s > ret.curve_index) break;
            // special case for the last point
            if (s == ret.curve_index && ret.path_index + 1 == size()) {
                --ret.curve_index;
                ret.t = 1;
                break;
            }
            ret.curve_index -= s;
        }
        return ret;
    }
    Path[] _data;
}

unittest
{
    import geom.bezier_curve;
    import core.exception;

    PathSequence p = new PathSequence;
    assert(p.empty() == true);
    assert(p.size() == 0);
    assert(p.curveCount() == 0);
    
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
