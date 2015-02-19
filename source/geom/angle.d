/*
 *  Various trigoniometric helper functions
 *
 *  Authors:
 *   Johan Engelen <goejendaagh@zonnet.nl>
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

module geom.angle;

import geom.coord;
import geom.point;
import std.math;

alias M_PI = std.math.PI;

/** Wrapper for angular values.
 *
 * This class is a convenience wrapper that implements the behavior generally expected of angles,
 * like addition modulo \f$2\pi\f$. The value returned from the default conversion
 * to <tt>double</tt> is in the range \f$[-\pi, \pi)\f$ - the convention used by C's
 * math library.
 *
 * @ingroup Primitives
 */
struct Angle
{
    this(Coord v) { _angle = v; _normalize(); } // this can be called implicitly
    this(Point p) { _angle = atan2(p); _normalize(); }
    this(Point a, Point b) { _angle = angle_between(a, b); _normalize(); }
    alias radians this;

    /** Get the angle as radians.
     * @return Number in range \f$[-\pi, \pi)\f$. */
    Coord radians() const
    { return _angle >= M_PI ? _angle - 2*M_PI : _angle; }

    /** Get the angle as positive radians.
     * @return Number in range \f$[0, 2\pi)\f$. */
    Coord radians0() const { return _angle; }

    /** Get the angle as degrees in math convention.
     * @return Number in range [-180, 180) obtained by scaling the result of radians()
     *         by \f$180/\pi\f$. */
    Coord degrees() const { return radians() * (180.0 / M_PI); }

    /** Get the angle as degrees in clock convention.
     * This method converts the angle to the "clock convention": angles start from the +Y axis
     * and grow clockwise. This means that 0 corresponds to \f$\pi/2\f$ radians,
     * 90 to 0 radians, 180 to \f$-\pi/2\f$ radians, and 270 to \f$\pi\f$ radians.
     * @return A number in the range [0, 360).
     */
    Coord degreesClock() const
    {
        Coord ret = 90.0 - _angle * (180.0 / M_PI);
        if (ret < 0) ret += 360;
        return ret;
    }

    /** Create an angle from its measure in radians. */
    static Angle from_radians(Coord d) { return Angle(d); }

    /** Create an angle from its measure in degrees. */
    static Angle from_degrees(Coord d) { return Angle(d * (M_PI / 180.0)); }

    /** Create an angle from its measure in degrees in clock convention.
     * See: Angle.degreesClock() */
    static Angle from_degrees_clock(Coord d)
    {
        // first make sure d is in [0, 360)
        d = fmod(d, 360.0);
        if (d < 0) d += 360.0;
        Coord rad = M_PI/2 - d * (M_PI / 180.0);
        if (rad < 0) rad += 2*M_PI;
        Angle a;
        a._angle = rad;
        return a;
    }

    Angle opBinary(string op, T : Angle)(T b) const
    {
        Angle x = Angle(this);
        static if (op == "+") {
            x._angle += b._angle;
            x._normalize();
            return x;
        } else static if (op == "-") {
            x._angle -= b._angle;
            x._normalize();
            return x;
        } else static assert(false, "Angle operator "~op~" not implemented");
    }

    void opOpAssign(string op, T)(T b)
    { mixin("this = this "~op~" Angle(b);"); }

private:
    void _normalize()
    { _angle -= floor(_angle * (1.0/(2*M_PI))) * 2*M_PI; }

    Coord _angle = 0; // this is always in [0, 2pi)
}

/** Directed angular interval.
 *
 * Wrapper for directed angles with defined start and end values. Useful e.g. for representing
 * the portion of an ellipse in an elliptical arc. Both extreme angles are contained
 * in the interval (it is a closed interval). Angular intervals can also be interptered
 * as functions \f$f: [0, 1] \to [-\pi, \pi)\f$, which return the start angle for 0,
 * the end angle for 1, and interpolate linearly for other values. Note that such functions
 * are not continuous if the interval contains the zero angle.
 *
 * This class is immutable - you cannot change the values of start and end angles
 * without creating a new instance of this class.
 */

/* this needs to be a class so that EllipticalArc can inherit from it */
class AngleInterval
{
    //@disable this();

    this(Angle s, Angle e, bool cw = false)
    { _start_angle = s; _end_angle = e; _sweep = cw; }

    this(Coord s, Coord e, bool cw = false)
    { _start_angle = s; _end_angle = e; _sweep = cw; }

    this(in AngleInterval o)
    { _start_angle = o._start_angle; _end_angle = o._end_angle; _sweep = o._sweep; }

    /** Get the angular coordinate of the interval's initial point
     * @return Angle in range \f$[0,2\pi)\f$ corresponding to the start of arc */
    Angle initialAngle() const { return _start_angle; }

    /** Get the angular coordinate of the interval's final point
     * @return Angle in range \f$[0,2\pi)\f$ corresponding to the end of arc */
    Angle finalAngle() const { return _end_angle; }
    bool is_degenerate() const { return initialAngle() == finalAngle(); }

    /** Get an angle corresponding to the specified time value. */
    Angle angleAt(Coord t) const
    {
        Coord span = extent();
        Angle ret = _start_angle.radians0() + span * (_sweep ? t : -t);
        return ret;
    }

    Angle opCall(Coord t) const { return angleAt(t); }

    /** Check whether the interval includes the given angle. */
    bool contains(Angle a) const
    {
        Coord s = _start_angle.radians0();
        Coord e = _end_angle.radians0();
        Coord x = a.radians0();
        if (_sweep) {
            if (s < e) return x >= s && x <= e;
            return x >= s || x <= e;
        } else {
            if (s > e) return x <= s && x >= e;
            return x <= s || x >= e;
        }
    }

    /** Extent of the angle interval.
     * @return Extent in range \f$[0, 2\pi)\f$ */
    Coord extent() const
    {
        Coord d = _end_angle - _start_angle;
        if (!_sweep) d = -d;
        if (d < 0) d += 2*M_PI;
        return d;
    }

protected:
    Angle _start_angle;
    Angle _end_angle;
    bool _sweep;

    this() {}
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
