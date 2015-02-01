/**
 * SVG 1.1-compliant elliptical arc curve
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

module geom.svg_elliptical_arc;

public import geom.coord;
import geom.angle;
import geom.d2;
import geom.curve;
import geom.elliptical_arc;
import geom.interval;
import geom.point;
import geom.rect;
import geom.sbasis;

class SVGEllipticalArc : EllipticalArc
{
    this() {}
    this(Point ip, Coord rx, Coord ry, Coord rot_angle, bool large_arc, bool sweep, Point fp)
    {
        _initial_point = ip;
        _final_point = fp;
        _rays[X] = rx; _rays[Y] = ry;
        _rot_angle = Angle(rot_angle);
        _large_arc = large_arc;
        _sweep = sweep;
        _updateCenterAndAngles(true);
    }

    override SVGEllipticalArc duplicate() const { return new SVGEllipticalArc(this); }

    override Coord valueAt(Coord t, size_t d) const
    {
        if (isChord()) return chord().valueAt(t, d);
        return super.valueAt(t, d);
    }

    override Point pointAt(Coord t) const
    {
        if (isChord()) return chord().pointAt(t);
        return super.pointAt(t);
    }

    override Point[] pointAndDerivatives(Coord t, uint n) const
    {
        if (isChord()) return chord().pointAndDerivatives(t, n);
        return super.pointAndDerivatives(t, n);
    }

    override Rect boundsExact() const
    {
        if (isChord()) return chord().boundsExact();
        return super.boundsExact();
    }

    override Rect boundsLocal(in Interval i, uint deg) const
    {
        if (isChord()) return chord().boundsLocal(i, deg);
        return super.boundsLocal(i, deg);
    }

    override Coord[] roots(Coord v, size_t d) const
    {
        if (isChord()) return chord().roots(v, d);
        return super.roots(v, d);
    }

    override D2!SBasis toSBasis() const
    {
        if (isChord()) return chord().toSBasis();
        return super.toSBasis();
    }

    override bool isSVGCompliant() const { return true; }

protected:
    // TODO move SVG-specific behavior here.
    //void _updateCenterAndAngles();

    this(in SVGEllipticalArc o)
    {
        _initial_point = o._initial_point;
        _final_point = o._final_point;
        _rays = o._rays;
        _rot_angle = o._rot_angle;
        _large_arc = o._large_arc;
        _sweep = o._sweep;
        _updateCenterAndAngles(true);
    }
} // end class SVGEllipticalArc

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
