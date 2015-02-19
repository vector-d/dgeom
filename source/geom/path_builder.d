/*
 * callback interface for SVG path data
 *
 * Authors:
 *      MenTaLguY <mental@rydia.net>
 *      Liam P. White
 *
 * Copyright (C) 2007-2015 authors
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

module geom.path_builder;

import geom.coord;
import geom.curve;
import geom.path;
import geom.path_sequence;
import geom.point;
import geom.rect;

/** Callback interface for processing path data.
 *
 * PathSink provides an interface that allows one to easily write
 * code which processes path data, for instance when converting
 * between path formats used by different graphics libraries.
 * It is also useful for writing algorithms which must do something
 * for each curve in the path.
 *
 * To store a path in a new format, implement the methods
 * for segments in a derived class and call feed().
 */
abstract class PathSink {

    /** Move to a different point without creating a segment.
     * Usually starts a new subpath. */
    void moveTo(in Point p);

    /// Output a line segment.
    void lineTo(in Point p);

    /// Output a quadratic Bezier segment.
    void curveTo(in Point c0, in Point c1, in Point p);

    /// Output a cubic Bezier segment.
    void quadTo(in Point c, in Point p);

    /** Output an elliptical arc segment.
     * See the EllipticalArc class for the documentation of parameters. */
    void arcTo(Coord rx, Coord ry, Coord angle, bool large_arc, bool sweep, in Point p);

    /// Close the current path with a line segment.
    void closePath();

    /** Flush any internal state of the generator.
     * This call should implicitly finish the current subpath.
     * Calling this method should be idempotent, because the default
     * implementations of path() and pathvector() will call it
     * multiple times in a row. */
    void flush();

    // Get the current point, e.g. where the initial point of the next segment will be.
    abstract Point currentPoint() const;

    /** Undo the last segment.
     * This method is optional.
     * @return true true if a segment was erased, false otherwise. */
    bool backspace() { return false; }

}

class PathBuilder : PathSink
{
    import geom.bezier_curve;
    import geom.svg_elliptical_arc;

    this()
    {
        _path = new Path;
        _pathset = new PathSequence;
    }

    override void moveTo(in Point p)
    {
        flush();
        _path.start(p);
        _start_p = p;
        _in_path = true;
    }
//TODO: what if _in_path = false?

    override void lineTo(in Point p)
    {
        // check for implicit moveto, like in: "M 1,1 L 2,2 z l 2,2 z"
        if (!_in_path)
            moveTo(_start_p);
        _path.appendNew!LineSegment(p);
    }

    override void quadTo(in Point c, in Point p)
    {
        // check for implicit moveto, like in: "M 1,1 L 2,2 z l 2,2 z"
        if (!_in_path)
            moveTo(_start_p);
        _path.appendNew!QuadraticBezier(c, p);
    }

    override void curveTo(in Point c0, in Point c1, in Point p)
    {
        // check for implicit moveto, like in: "M 1,1 L 2,2 z l 2,2 z"
        if (!_in_path)
            moveTo(_start_p);
        _path.appendNew!CubicBezier(c0, c1, p);
    }

    override void arcTo(Coord rx, Coord ry, Coord angle, bool large_arc, bool sweep, in Point p)
    {
        // check for implicit moveto, like in: "M 1,1 L 2,2 z l 2,2 z"
        if (!_in_path)
            moveTo(_start_p);
        _path.appendNew!SVGEllipticalArc(rx, ry, angle, large_arc, sweep, p);
    }

    override bool backspace()
    {
        if (_in_path && _path.size() > 0) {
            _path.erase_last();
            return true;
        }
        return false;
    }

    void append(in Path other, Stitching stitching = Stitching.NO_STITCHING)
    {
        if (!_in_path)
            moveTo(other.initialPoint());
        _path.append(other, stitching);
    }

    override void closePath()
    {
        _path.close();
        flush();
    }

    override void flush()
    {
        if (_in_path) {
            _in_path = false;
            _pathset ~= new Path(_path);
            _path.clear();
        }
    }

    override Point currentPoint() const { return _path.finalPoint(); }

    /// Retrieve the path
    const(PathSequence) peek() const { return _pathset; }

    /// Clear the stored path sequence
    void clear()
    {
        _in_path = false;
        _path.clear();
        _pathset.clear();
    }

protected:
    bool _in_path = false;
    PathSequence _pathset;
    Path _path;
    Point _start_p;
}

unittest
{
    PathBuilder pb = new PathBuilder;
    pb.moveTo(Point());
    pb.lineTo(Point(1,1));
    pb.quadTo(Point(2,2), Point(3,3));
    pb.curveTo(Point(4,4), Point(5,5), Point(6,6));
    pb.arcTo(7, 7, 0, false, true, Point(20,20));
    pb.backspace();
    Point p = pb.currentPoint();
    assert(p == Point(6, 6));
    pb.flush();
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
