/*
 * Elliptical arc curve
 *
 * Authors:
 *    MenTaLguY <mental@rydia.net>
 *    Marco Cecchetti <mrcekets at gmail.com>
 *    Krzysztof Kosiński <tweenk.pl@gmail.com>
 *    Liam P. White
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

module geom.elliptical_arc;

import std.math;

public import geom.coord;
import geom.angle;
import geom.affine;
import geom.bezier_curve;
import geom.curve;
import geom.d2;
import geom.interval;
import geom.point;
import geom.rect;
import geom.sbasis;
import geom.transforms;

class EllipticalArc : AngleInterval, Curve
{
    /** Creates an arc with all variables set to zero, and both flags to true. */
    this()
    {
        super(0, 0, true);
        _initial_point = Point(0,0);
        _final_point = Point(0,0);
        _rays = Point(0,0);
        _center = Point(0,0);
        _rot_angle = 0;
        _large_arc = true;
    }

    /** Create a new elliptical arc.
     * @param ip Initial point of the arc
     * @param rx First ray of the ellipse
     * @param ry Second ray of the ellipse
     * @param rot Angle of rotation of the X axis of the ellipse in radians
     * @param large If true, the large arc is chosen (always >= 180 degrees), otherwise
     *              the smaller arc is chosen
     * @param sweep If true, the clockwise arc is chosen, otherwise the counter-clockwise
     *              arc is chosen
     * @param fp Final point of the arc */
    this(Point ip, Coord rx, Coord ry, Coord rot_angle, bool large_arc, bool sweep, Point fp)
    {
        super(0, 0, sweep);
        _initial_point = ip;
        _final_point = fp;
        _rays = Point(rx, ry);
        _rot_angle = rot_angle;
        _large_arc = large_arc;
        _updateCenterAndAngles(false);
    }

    // methods new to EllipticalArc go here

    /+ Retrieve and modify parameters +/

    /** Get the interval of angles the arc contains
     * @return The interval between the final and initial angles of the arc */
    Interval angleInterval() const { return Interval(initialAngle(), finalAngle()); }

    /** Get the defining ellipse's rotation
     * @return Angle between the +X ray of the ellipse and the +X axis */
    Angle rotationAngle() const { return _rot_angle; }

    /** Get one of the ellipse's rays
     * @param d Dimension to retrieve
     * @return The selected ray of the ellipse */
    Coord ray(size_t d) const { return _rays[d]; }

    /** Get both rays as a point
     * @return Point with X equal to the X ray and Y to Y ray */
    Point rays() const { return _rays; }

    /** Whether the arc is larger than half an ellipse.
     * @return True if the arc is larger than \f$\pi\f$, false otherwise */
    bool largeArc() const { return _large_arc; }

    /** Whether the arc turns clockwise
     * @return True if the arc makes a clockwise turn when going from initial to final
     *         point, false otherwise */
    bool sweep() const { return _sweep; }

    /** Get the line segment connecting the arc's endpoints.
     * @return A linear segment with initial and final point correspoding to those of the arc. */
    LineSegment chord() const { return bezierFromPoints(_initial_point, _final_point); }

    /// Check whether both rays are nonzero
    bool isChord() const { return _rays[0] == 0 || _rays[Y] == 0; }

    /** Change the arc's parameters. */ 
    void set(in Point ip, Coord rx, Coord ry, Coord rot_angle, bool large_arc, bool sweep, in Point fp)
    {
        _initial_point = ip;
        _final_point = fp;
        _rays[X] = rx;
        _rays[Y] = ry;
        _rot_angle = Angle(rot_angle);
        _large_arc = large_arc;
        _sweep = sweep;
        _updateCenterAndAngles(isSVGCompliant());
    }

    void setRays(Coord rx, Coord ry)
    {
        _rays[X] = rx;
        _rays[Y] = ry;
        _updateCenterAndAngles(isSVGCompliant());
    }

    /** Change the initial and final point in one operation.
     * This method exists because modifying any of the endpoints causes rather costly
     * recalculations of the center and extreme angles.
     * @param ip New initial point
     * @param fp New final point */
    void setExtremes(in Point ip, in Point fp)
    {
        _initial_point = ip;
        _final_point = fp;
        _updateCenterAndAngles(isSVGCompliant());
    }

    /+ Access computed parameters of the arc +/

    Coord center(size_t d) const { return _center[d]; }

    /** Get the arc's center
     * @return The arc's center, situated on the intersection of the ellipse's rays */
    Point center() const { return _center; }

    /** Get the extent of the arc
     * @return The angle between the initial and final point, in arc's angular coordinates */
    Coord sweepAngle() const { return extent(); }

    /+ Angular evaluation +/

    /** Check whether the arc contains the given angle
     * @param t The angle to check
     * @return True if the arc contains the angle, false otherwise */
    bool containsAngle(Coord angle) const { return contains(Angle(angle)); }

    /** Evaluate the arc at the specified angular coordinate
     * @param t Angle
     * @return Point corresponding to the given angle */
    Point pointAtAngle(Coord t) const { return Point.polar(t) * unitCircleTransform(); }

    /** Evaluate one of the arc's coordinates at the specified angle
     * @param t Angle
     * @param d The dimension to retrieve
     * @return Selected coordinate of the arc at the specified angle */
    Coord valueAtAngle(Coord t, size_t d) const
    {
        Coord sinrot = sin(_rot_angle), cosrot = cos(_rot_angle);
        Coord cost = cos(t), sint = sin(t);

        if (d == X) {
            return    ray(X) * cosrot * cost
                    - ray(Y) * sinrot * sint
                    + center(X);
        } else {
            return    ray(X) * sinrot * cost
                    + ray(Y) * cosrot * sint
                    + center(Y);
        }
    }

    /** Retrieve the unit circle transform.
     * Each ellipse can be interpreted as a translated, scaled and rotate unit circle.
     * This function returns the transform that maps the unit circle to the arc's ellipse.
     * @return Transform from unit circle to the arc's ellipse */
    Affine unitCircleTransform() const
    {
        Affine ret = Scale(ray(X), ray(Y)) * Rotate(_rot_angle);
        ret.setTranslation(center());
        return ret;
    }

    /** Check whether the arc adheres to SVG 1.1 implementation guidelines */
    bool isSVGCompliant() const { return false; }

    EllipticalArc[2] subdivide(Coord t) const
    {
        EllipticalArc arc1 = portion(0, t);
        EllipticalArc arc2 = portion(t, 1);
        assert(arc1 !is null && arc2 !is null);
        return [arc1, arc2];
    }
    
    /+ Curve interface +/

    Point initialPoint() const { return _initial_point; }
    Point finalPoint() const { return _final_point; }
    EllipticalArc duplicate() const { return new EllipticalArc(this); }

    void setInitial(in Point p)
    {
        _initial_point = p;
        _updateCenterAndAngles(isSVGCompliant());
    }
    void setFinal(in Point p)
    {
        _final_point = p;
        _updateCenterAndAngles(isSVGCompliant());
    }

    bool isDegenerate() const { return ( are_near(ray(X), 0) || are_near(ray(Y), 0) ); }
    Rect boundsFast() const { return boundsExact(); }
    Rect boundsExact() const
    {
        import std.algorithm : swap;

        Coord extremes[4];
        Coord sinrot = sin(_rot_angle), cosrot = cos(_rot_angle);

        extremes[0] = atan2(-ray(Y) * sinrot, ray(X) * cosrot);
        extremes[1] = extremes[0] + M_PI;
        if (extremes[0] < 0) extremes[0] += 2*M_PI;
        extremes[2] = atan2(ray(Y) * cosrot, ray(X) * sinrot);
        extremes[3] = extremes[2] + M_PI;
        if (extremes[2] < 0) extremes[2] += 2*M_PI;

        Coord arc_extremes[4];
        arc_extremes[0] = initialPoint()[X];
        arc_extremes[1] = finalPoint()[X];
        if (arc_extremes[0] < arc_extremes[1])
            swap(arc_extremes[0], arc_extremes[1]);
        arc_extremes[2] = initialPoint()[Y];
        arc_extremes[3] = finalPoint()[Y];
        if (arc_extremes[2] < arc_extremes[3])
            swap(arc_extremes[2], arc_extremes[3]);

        if (!are_near(initialPoint(), finalPoint())) {
            for (uint i = 0; i < 4; ++i) {
                if (containsAngle(extremes[i])) {
                    arc_extremes[i] = valueAtAngle(extremes[i], (i >> 1) ? Y : X); // wtf, 1>>1
                }
            }
        }

        return Rect( Point(arc_extremes[1], arc_extremes[3]) ,
                     Point(arc_extremes[0], arc_extremes[2]) );
    }

    // TODO: native implementation of the following methods
    Rect boundsLocal(in Interval i, uint deg) const
    {
        import geom.sbasis_curve;
        auto x = new SBasisCurve(toSBasis());
        return x.boundsLocal(i, deg);
    }

    Coord[] roots(Coord v, size_t d) const
    {
        import geom.sbasis_curve;
        auto x = new SBasisCurve(toSBasis());
        return x.roots(v, d);
    }

    int degreesOfFreedom() const { return 7; }

    EllipticalArc derivative() const
    {
        // D(E(t,C),t) = E(t+PI/2,O), where C is the ellipse center
        // the derivative doesn't rotate the ellipse but there is a translation
        // of the parameter t by an angle of PI/2 so the ellipse points are shifted
        // of such an angle in the cw direction
        EllipticalArc result = duplicate();
        result._center[X] = result._center[Y] = 0;
        result._start_angle += M_PI/2;
        if (result._start_angle >= 2*M_PI)
            result._start_angle -= 2*M_PI;

        result._end_angle += M_PI/2;
        if (result._end_angle >= 2*M_PI)
            result._end_angle -= 2*M_PI;

        result._initial_point = result.pointAtAngle(result.initialAngle());
        result._final_point = result.pointAtAngle(result.finalAngle());
        return result;
    }
    
    Coord length(Coord tolerance = 0.01) const
    {
        import geom.sbasis_curve;
        auto x = new SBasisCurve(toSBasis());
        return x.length(tolerance);
    }

    EllipticalArc opBinary(string op)(in Translate m) const if (op == "*")
    {
        EllipticalArc d = new EllipticalArc(this);
        d._initial_point += m.vector();
        d._final_point += m.vector();
        d._center += m.vector();
        return d;
    }

    /**
    *  The size of the returned array equals n+1.
    */
    Point[] pointAndDerivatives(Coord t, uint n) const
    {
        import object : reserve;
        uint nn = n+1; // nn represents the size of the result vector.
        Point[] result;
        result.reserve(nn);
        Coord angle = map_unit_interval_on_circular_arc(t, initialAngle(),
                                                         finalAngle(), _sweep);
        EllipticalArc ea = duplicate();
        ea._center = Point(0,0);
        uint m = cast(uint)fmin(nn, 4u);

        for (uint i = 0; i < m; ++i) {
            result ~= ea.pointAtAngle(angle);
            angle += (_sweep ? M_PI/2 : -M_PI/2);
            if (angle >= 2*M_PI) angle -= 2*M_PI;
        }
        m = nn / 4;
        for (uint i = 1; i < m; ++i) {
            for (uint j = 0; j < 4; ++j)
                result ~= result[j];
        }
        m = nn - 4 * m;
        for (uint i = 0; i < m; ++i) {
            result ~= result[i];
        }
        if (result.length != 0) // nn != 0
            result[0] = pointAtAngle(angle);
        return result;
    }

    D2!SBasis toSBasis() const
    {
        import geom.linear;
        D2!SBasis arc;
        // the interval of parametrization has to be [0,1]
        Coord et = initialAngle().radians() + ( _sweep ? sweepAngle() : -sweepAngle() );
        Linear param = Linear(initialAngle(), et);
        Coord cos_rot_angle = std.math.cos(_rot_angle), sin_rot_angle = std.math.sin(_rot_angle);

        // order = 4 seems to be enough to get a perfect looking elliptical arc
        SBasis arc_x = cos(param,4) * ray(X);
        SBasis arc_y = sin(param,4) * ray(Y);
        arc[0] = arc_x * cos_rot_angle - arc_y * sin_rot_angle + SBasis(Linear(center(X),center(X)));
        arc[1] = arc_x * sin_rot_angle + arc_y * cos_rot_angle + SBasis(Linear(center(Y),center(Y)));

        // ensure that endpoints remain exact
        for (uint d = 0 ; d < 2; d++) {
            arc[d][0][0] = initialPoint()[d];
            arc[d][0][1] = finalPoint()[d];
        }

        return arc;
    }
    Coord valueAt(Coord t, size_t d) const { return valueAtAngle(angleAt(t), d); }
    Point pointAt(Coord t) const { return pointAtAngle(angleAt(t)); }

    EllipticalArc portion(Coord f, Coord t) const
    {
        // fix input arguments
        if (f < 0) f = 0;
        if (f > 1) f = 1;
        if (t < 0) t = 0;
        if (t > 1) t = 1;

        EllipticalArc arc = duplicate();

        // TODO: kill this, we're not THAT imprecise...
        if (are_near(f, t)) {
            arc._center = arc._initial_point = arc._final_point = pointAt(f);
            arc._start_angle = arc._end_angle = _start_angle;
            arc._rot_angle = _rot_angle;
            arc._sweep = _sweep;
            arc._large_arc = _large_arc;
            return arc;
        }

        arc._initial_point = pointAt(f);
        arc._final_point = pointAt(t);
        if (f > t) arc._sweep = !_sweep;
        if (_large_arc && fabs(sweepAngle() * (t-f)) < M_PI)
            arc._large_arc = false;
        arc._updateCenterAndAngles(arc.isSVGCompliant()); // TODO: be more clever
        return arc;
    }

    EllipticalArc reverse() const
    {
        // the arc is the same but traversed in the opposite direction
        EllipticalArc rarc = duplicate();
        rarc._sweep = !_sweep;
        rarc._initial_point = _final_point;
        rarc._final_point = _initial_point;
        rarc._start_angle = _end_angle;
        rarc._end_angle = _start_angle;
        rarc._updateCenterAndAngles(rarc.isSVGCompliant());
        return rarc;
    }
    
    void transform(in Affine m)
    {
        import geom.ellipse;
        auto e = Ellipse(center(X), center(Y), ray(X), ray(Y), _rot_angle);
        auto et = e.transformed(m);
        Point inner_point = pointAt(0.5);
        auto x = et.arc(initialPoint() * m, inner_point * m, finalPoint() * m, isSVGCompliant());
        _initial_point = x._initial_point; _final_point = x._final_point;
        _rays = x._rays;
        _rot_angle = x._rot_angle; _large_arc = x._large_arc;
        _sweep = x._sweep;
        _updateCenterAndAngles(isSVGCompliant());
    }

protected:

    /* NOTE: this implementation follows Standard SVG 1.1 implementation guidelines
     * for elliptical arc curves. See Appendix F.6.
     */
    void _updateCenterAndAngles(bool svg)
    {
        Point d = initialPoint() - finalPoint();
        if (svg) {
            if (initialPoint() == finalPoint()) {
                _rot_angle = _start_angle = _end_angle = Angle(0);
                _center = initialPoint();
                _rays = Point(0,0);
                _large_arc = _sweep = false;
                return;
            }

            _rays[X] = fabs(_rays[X]);
            _rays[Y] = fabs(_rays[Y]);

            if (are_near(ray(X), 0) || are_near(ray(Y), 0)) {
                _rays[X] = L2(d) / 2;
                _rays[Y] = 0;
                _rot_angle = Angle(atan2(d[Y], d[X]));
                _start_angle = Angle(0);
                _end_angle = Angle(M_PI);
                _center = middle_point(initialPoint(), finalPoint());
                _large_arc = false;
                _sweep = false;
                return;
            }
        } else {
            if (are_near(initialPoint(), finalPoint())) {
                if (are_near(ray(X), 0) && are_near(ray(Y), 0)) {
                    _start_angle = _end_angle = Angle(0);
                    _center = initialPoint();
                    return;
                } else {
                    throw new Exception("initial and final point are the same");
                }
            }
            if (are_near(ray(X), 0) && are_near(ray(Y), 0)) { // but initialPoint != finalPoint
                throw new Exception(
                    "there is no ellipse that satisfies the given constraints: "
                    "ray(X) == 0 && ray(Y) == 0 but initialPoint != finalPoint");
            }
            if (are_near(ray(Y), 0)) {
                Point v = initialPoint() - finalPoint();
                if (are_near(L2sq(v), 4 * ray(X) * ray(X))) {
                    auto angle = Angle(v);
                    if (are_near( angle, _rot_angle )) {
                        _start_angle = Angle(0);
                        _end_angle = Angle(M_PI);
                        _center = v/2 + finalPoint();
                        return;
                    }
                    angle -= M_PI;
                    if (are_near( angle, _rot_angle )) {
                        _start_angle = Angle(M_PI);
                        _end_angle = Angle(0);
                        _center = v/2 + finalPoint();
                        return;
                    }
                    throw new Exception(
                        "there is no ellipse that satisfies the given constraints: "
                        "ray(Y) == 0 "
                        "and slope(initialPoint - finalPoint) != rotation_angle "
                        "and != rotation_angle + PI");
                }
                bool gt = L2sq(v) > 4*ray(X)*ray(X);
                throw new Exception(
                    "there " ~ (gt ? "is no ellipse" : "are infinite ellipses") ~
                    " that satisfies the given constraints: "
                    "ray(Y) == 0 and distance(initialPoint, finalPoint) " ~ (gt ? ">" : "<") ~" 2*ray(X)");
            }

            if (are_near(ray(X), 0)) {
                Point v = initialPoint() - finalPoint();
                if (are_near(L2sq(v), 4*ray(Y)*ray(Y))) {
                    Coord angle = atan2(v[Y], v[X]);
                    if (angle < 0) angle += 2*M_PI;
                    Coord rot_angle = _rot_angle.radians() + M_PI/2;
                    if (rot_angle >= 2*M_PI) rot_angle -= 2*M_PI;
                    if (are_near(angle, rot_angle)) {
                        _start_angle = Angle(M_PI/2);
                        _end_angle = Angle(3*M_PI/2);
                        _center = v/2 + finalPoint();
                        return;
                    }
                    angle -= M_PI;
                    if (angle < 0) angle += 2*M_PI;
                    if (are_near( angle, rot_angle )) {
                        _start_angle = Angle(3*M_PI/2);
                        _end_angle = Angle(M_PI/2);
                        _center = v/2 + finalPoint();
                        return;
                    }
                    throw new Exception(
                        "there is no ellipse that satisfies the given constraints: "
                        "ray(X) == 0 "
                        "and slope(initialPoint - finalPoint) != rotation_angle + PI/2 "
                        "and != rotation_angle + (3/2)*PI");
                }
                bool gt = L2sq(v) > 4*ray(Y)*ray(Y);
                throw new Exception(
                    "there " ~ (gt ? "is no ellipse" : "are infinite ellipses") ~
                    " that satisfies the given constraints: "
                    "ray(X) == 0 and distance(initialPoint, finalPoint) " ~ (gt ? ">" : "<") ~" 2*ray(Y)");
            }
        }

        Rotate rm = Rotate(_rot_angle);
        Affine m = rm;
        m[1] = -m[1];
        m[2] = -m[2];

        Point p = (d / 2) * m;
        Coord rx2 = _rays[X] * _rays[X];
        Coord ry2 = _rays[Y] * _rays[Y];
        Coord rxpy = _rays[X] * p[Y];
        Coord rypx = _rays[Y] * p[X];
        Coord rx2py2 = rxpy * rxpy;
        Coord ry2px2 = rypx * rypx;
        Coord num = rx2 * ry2;
        Coord den = rx2py2 + ry2px2;
        assert(den != 0);
        Coord rad = num / den;
        Point c = Point(0,0);

        if (rad > 1) {
            rad -= 1;
            rad = sqrt(rad);

            if (_large_arc == _sweep) rad = -rad;
            c = Point(rxpy / ray(Y), -rypx / ray(X)) * rad;
            _center = c * rm + middle_point(initialPoint(), finalPoint());
        } else if (rad == 1 || svg) {
            Coord lamda = sqrt(1 / rad);
            _rays[X] *= lamda;
            _rays[Y] *= lamda;
            _center = middle_point(initialPoint(), finalPoint());
        } else {
            throw new Exception("there is no ellipse that satisfies the given constraints");
        }

        Point sp = Point((p[X] - c[X]) / ray(X), (p[Y] - c[Y]) / ray(Y));
        Point ep = Point((-p[X] - c[X]) / ray(X), (-p[Y] - c[Y]) / ray(Y));
        Point v = Point(1, 0);
        _start_angle = Angle(angle_between(v, sp));
        Coord sweep_angle = angle_between(sp, ep);
        if (!_sweep && sweep_angle > 0) sweep_angle -= 2*M_PI;
        if (_sweep && sweep_angle < 0) sweep_angle += 2*M_PI;

        _end_angle = _start_angle;
        _end_angle += sweep_angle;
    }

    Point _initial_point, _final_point;
    Point _rays, _center;
    Angle _rot_angle;
    bool _large_arc;
    
    this(in EllipticalArc o)
    {
        super(o._start_angle, o._end_angle, o._sweep);
        _initial_point = o._initial_point;
        _final_point = o._final_point;
        _rays = o._rays;
        _center = o._center;
        _rot_angle = o._rot_angle;
        _large_arc = o._large_arc;
    }

    Coord map_to_01(Coord angle) const
    { return map_circular_arc_on_unit_interval(angle, initialAngle(), finalAngle(), _sweep); }
}

/* start_angle and angle must belong to [0, 2PI[
 * and angle must belong to the cirsular arc defined by
 * start_angle, end_angle and with rotation direction cw
 */
private Coord map_circular_arc_on_unit_interval(Coord angle, Coord start_angle, Coord end_angle, bool cw = true)
{
    Coord d = end_angle - start_angle;
    Coord t = angle - start_angle;
    if (!cw) {
    	d = -d;
    	t = -t;
    }
    d = fmod(d, 2*M_PI);
    t = fmod(t, 2*M_PI);
    if (d < 0) d += 2*M_PI;
    if (t < 0) t += 2*M_PI;
    return t / d;
}

private Coord map_unit_interval_on_circular_arc(Coord t, Coord start_angle, Coord end_angle, bool cw = true)
{
    Coord sweep_angle = end_angle - start_angle;
    if (!cw) sweep_angle = -sweep_angle;
    sweep_angle = fmod(sweep_angle, 2*M_PI);
    if (sweep_angle < 0) sweep_angle += 2*M_PI;

	Coord angle = start_angle;
    if (cw)
        angle += sweep_angle * t;
    else
        angle -= sweep_angle * t;
    angle = fmod(angle, 2*M_PI);
    if (angle < 0) angle += 2*M_PI;
    return angle;
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
