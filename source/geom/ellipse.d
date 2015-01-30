/**
 * Ellipse Curve
 *
 * Authors:
 *      Marco Cecchetti <mrcekets at gmail.com>
 *      Liam P. White
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

module geom.ellipse;

import std.math;

public import geom.coord;
import geom.affine;
import geom.elliptical_arc;
import geom.point;

struct Ellipse
{
    Point center() const { return m_center; }
    Coord center(size_t d) const { return m_center[d]; }
    Coord ray(size_t d) const { return m_ray[d]; }
    Coord rot_angle() const { return m_angle; }

    this(Coord cx, Coord cy, Coord rx, Coord ry, Coord a) { set(cx, cy, rx, ry, a); }

    void set(Coord cx, Coord cy, Coord rx, Coord ry, Coord a)
    {
        m_center = Point(cx, cy);
        m_ray = Point(rx, ry);
        m_angle = a;
    }

    // build an ellipse by its implicit equation:
    // Ax^2 + Bxy + Cy^2 + Dx + Ey + F = 0
    this(Coord A, Coord B, Coord C, Coord D, Coord E, Coord F) { set(A, B, C, D, E, F); }

    void set(Coord A, Coord B, Coord C, Coord D, Coord E, Coord F)
    {
        import std.algorithm : swap;

        Coord den = 4*A*C - B*B;
        if (den == 0)
            throw new Exception("den == 0 while computing ellipse center");

        m_center[X] = (B*E - 2*C*D) / den;
        m_center[Y] = (B*D - 2*A*E) / den;

        // evaluate the a coefficient of the ellipse equation in normal form
        // E(x,y) = a*(x-cx)^2 + b*(x-cx)*(y-cy) + c*(y-cy)^2 = 1
        // where b = a*B , c = a*C, (cx,cy) == centre
        Coord num =    A * (m_center[X] * m_center[X])
                     + B * m_center[X] * m_center[Y]
                     + C * (m_center[Y] * m_center[Y])
                     - F;

        // evaluate ellipse rotation angle
        Coord rot = atan2( -B, -(A - C) )/2;
        bool swap_axes = false;
        if (are_near(rot, 0)) rot = 0;
        if (are_near(rot, std.math.PI/2) || rot < 0) {
            swap_axes = true;
        }

        // evaluate the length of the ellipse rays
        Coord cosrot = cos(rot);
        Coord sinrot = sin(rot);
        Coord cos2 = cosrot * cosrot;
        Coord sin2 = sinrot * sinrot;
        Coord cossin = cosrot * sinrot;

        den = A * cos2 + B * cossin + C * sin2;
        if (den == 0)
            throw new Exception("den == 0, while computing 'rx' coefficient");
        Coord rx2 =  num/den;
        if (rx2 < 0)
            throw new Exception("rx2 < 0, while computing 'rx' coefficient");

        Coord rx = sqrt(rx2);

        den = C * cos2 - B * cossin + A * sin2;
        if (den == 0)
            throw new Exception("den == 0, while computing 'ry' coefficient");

        Coord ry2 =  num/den;
        if (ry2 < 0)
            throw new Exception("ry2 < 0, while computing 'rx' coefficient");

        Coord ry = sqrt(ry2);

        // the solution is not unique so we choose always the ellipse
        // with a rotation angle between 0 and PI/2
        if (swap_axes) swap(rx, ry);
        if (    are_near(rot,  std.math.PI/2)
             || are_near(rot, -std.math.PI/2)
             || are_near(rx, ry)) {
            rot = 0;
        } else if (rot < 0) {
            rot += std.math.PI/2;
        }

        m_ray[X] = rx;
        m_ray[Y] = ry;
        m_angle = rot;
    }

    EllipticalArc arc(in Point initial, in Point inner, in Point fina, bool svg_compliant = true)
    {
        Point sp_cp = initial - center();
        Point ep_cp = fina    - center();
        Point ip_cp = inner   - center();

        Coord angle1 = angle_between(sp_cp, ep_cp);
        Coord angle2 = angle_between(sp_cp, ip_cp);
        Coord angle3 = angle_between(ip_cp, ep_cp);

        bool large_arc_flag = true;
        bool sweep_flag = true;

        if (angle1 > 0) {
            if (angle2 > 0 && angle3 > 0) {
                large_arc_flag = false;
                sweep_flag = true;
            } else {
                large_arc_flag = true;
                sweep_flag = false;
            }
        } else {
            if (angle2 < 0 && angle3 < 0) {
                large_arc_flag = false;
                sweep_flag = false;
            } else {
                large_arc_flag = true;
                sweep_flag = true;
            }
        }

        return new EllipticalArc(initial, ray(X), ray(Y), rot_angle(), large_arc_flag, sweep_flag, fina);
    }
    
    Ellipse transformed(in Affine m) const
    {
        import std.algorithm : swap;

        Coord cosrot = cos(rot_angle());
        Coord sinrot = sin(rot_angle());
        auto A = Affine(  ray(X) * cosrot, ray(X) * sinrot,
                  -ray(Y) * sinrot, ray(Y) * cosrot,
                   0,               0                );
        Point new_center = center() * m;
        Affine M = m.withoutTranslation();
        Affine AM = A * M;
        if (are_near(sqrt(fabs(AM.det())), 0)) {
            Coord angle;
            if (AM[0] != 0) {
                angle = atan2(AM[2], AM[0]);
            } else if (AM[1] != 0) {
                angle = atan2(AM[3], AM[1]);
            } else {
                angle = std.math.PI/2;
            }
            Point V = Point(cos(angle), sin(angle));
            V *= AM;
            Coord rx = L2(V);
            angle = atan2(V);
            return Ellipse(new_center[X], new_center[Y], rx, 0, angle);
        }

        Coord[] coeff = implicit_form_coefficients();
        auto Q = Affine( coeff[0],   coeff[1]/2,
                  coeff[1]/2, coeff[2],
                  0,          0   );

        Affine invm = M.inverse();
        Q = invm * Q;
        swap(invm[1], invm[2]);
        Q = Q * invm;
        Ellipse e = Ellipse(Q[0], 2*Q[1], Q[3], 0, 0, -1);
        e.m_center = new_center;
        return e;
    }

    Coord[6] implicit_form_coefficients() const
    {
        if (ray(X) == 0 || ray(Y) == 0)
            throw new Exception("a degenerate ellipse doesn't own an implicit form");

        Coord[6] coeff;
        Coord cosrot = cos(rot_angle());
        Coord sinrot = sin(rot_angle());
        Coord cos2 = cosrot * cosrot;
        Coord sin2 = sinrot * sinrot;
        Coord cossin = cosrot * sinrot;
        Coord invrx2 = 1 / (ray(X) * ray(X));
        Coord invry2 = 1 / (ray(Y) * ray(Y));

        coeff[0] = invrx2 * cos2 + invry2 * sin2;
        coeff[1] = 2 * (invrx2 - invry2) * cossin;
        coeff[2] = invrx2 * sin2 + invry2 * cos2;
        coeff[3] = -(2 * coeff[0] * center(X) + coeff[1] * center(Y));
        coeff[4] = -(2 * coeff[2] * center(Y) + coeff[1] * center(X));
        coeff[5] = coeff[0] * center(X) * center(X)
                 + coeff[1] * center(X) * center(Y)
                 + coeff[2] * center(Y) * center(Y)
                 - 1;
        return coeff;
    }

private:
    Point m_center, m_ray;
    Coord m_angle = 0;
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
