/**
 * Path outlining
 *
 * Authors:
 *   fred
 *   Liam P. White
 * 
 * Copyright 2003-2015 Authors
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

module geom.path_outline;

public import geom.coord;
import geom.curve;
import geom.d2;
import geom.bezier_curve;
import geom.ellipse;
import geom.elliptical_arc;
import geom.path;
import geom.path_builder;
import geom.path_sequence;
import geom.point;
import geom.sbasis;
import geom.sbasis_to_bezier;
import geom.svg_elliptical_arc;

    import std.stdio : writeln;

enum JoinType {
    JOIN_BEVEL,
    JOIN_ROUND,
    JOIN_MITER,
    JOIN_EXTRAPOLATE, // not in SVG 1
}

enum ButtType {
    BUTT_FLAT,
    BUTT_ROUND,
    BUTT_SQUARE,
    BUTT_PEAK, // why do we have this again?
}

PathSequence outline(in Path input, Coord width, Coord miter, JoinType join, ButtType butt)
{
    if (input.size == 0) return new PathSequence;

    PathBuilder res = new PathBuilder;
    Path with_dir = input.half_outline(width/2., miter, join);
    Path against_dir = input.reversed.half_outline(width/2., miter, join);

    res.moveTo(with_dir[0].initialPoint);
    res.append(with_dir);

    // glue caps
    if (!input.closed) {
        switch (butt) {
            case ButtType.BUTT_ROUND:
                res.arcTo((-width) / 2., (-width) / 2., 0., true, true, against_dir.initialPoint);
                break;
            case ButtType.BUTT_SQUARE:
            {
                Point end_deriv = -input[input.size-1].reverse.unitTangentAt(0.);
                Coord radius = 0.5 * with_dir.finalPoint.distance(against_dir.initialPoint);
                res.lineTo(with_dir.finalPoint + end_deriv*radius);
                res.lineTo(against_dir.initialPoint + end_deriv*radius);
                res.lineTo(against_dir.initialPoint);
                break;
            }
            case ButtType.BUTT_PEAK:
            {
                Point end_deriv = -input[input.size-1].reverse.unitTangentAt(0.);
                Coord radius = 0.5* with_dir.finalPoint.distance(against_dir.initialPoint);
                Point midpoint = ((with_dir.finalPoint + against_dir.initialPoint) * 0.5) + end_deriv*radius;
                res.lineTo(midpoint);
                res.lineTo(against_dir.initialPoint);
                break;
            }
            case ButtType.BUTT_FLAT:
            default:
                res.lineTo(against_dir.initialPoint);
                break;
        }
    } else {
        res.moveTo(against_dir.initialPoint);
    }

    res.append(against_dir);

    if (!input.closed) {
        switch(butt) {
            case ButtType.BUTT_ROUND:
                res.arcTo((-width) / 2., (-width) / 2., 0., true, true, with_dir.initialPoint);
                break;
            case ButtType.BUTT_SQUARE:
            {
                Point end_deriv = -input[0].unitTangentAt(0.);
                Coord radius = 0.5 * against_dir.finalPoint.distance(with_dir.initialPoint);
                res.lineTo(against_dir.finalPoint + end_deriv*radius);
                res.lineTo(with_dir.initialPoint + end_deriv*radius);
                res.lineTo(with_dir.initialPoint);
                break;
            }
            case ButtType.BUTT_PEAK:
            {
                Point end_deriv = -input[0].unitTangentAt(0.);
                Coord radius = 0.5 * against_dir.finalPoint.distance(with_dir.initialPoint);
                Point midpoint = ((against_dir.finalPoint + with_dir.initialPoint) * 0.5) + end_deriv*radius;
                res.lineTo(midpoint);
                res.lineTo(with_dir.initialPoint);
                break;
            }
            case ButtType.BUTT_FLAT:
            default:
                res.lineTo(with_dir.initialPoint);
        }
        res.closePath();
    }

    res.flush();
    const(PathSequence) ps = res.peek();
    return new PathSequence(res.peek());
}

private alias join_func = void function(ref Path res, in Curve outgoing, Coord miter, Coord width);

void bevel_join(ref Path res, in Curve outgoing, Coord /*miter*/, Coord /*width*/)
{
    res.appendNew!LineSegment(outgoing.initialPoint);
}

void round_join(ref Path res, in Curve outgoing, Coord /*miter*/, Coord width)
{
    res.appendNew!SVGEllipticalArc(width, width, 0, false, width > 0, outgoing.initialPoint());
}

void miter_join(ref Path res, in Curve outgoing, Coord miter, Coord width)
{
    const(Curve) incoming = res.back();
    Point tang1 = incoming.reverse.unitTangentAt(0);
    Point tang2 = outgoing.unitTangentAt(0);
    Point p = intersection_point(incoming.finalPoint, tang1, outgoing.initialPoint, tang2);
    if (p.isFinite) {
        // check size of miter
        Point point_on_path = incoming.finalPoint - tang1.rot90*width;
        Coord len = p.distance(point_on_path);
        if (len <= miter) {
            // miter OK, check to see if we can do a relocation
            if (auto line = cast(const(LineSegment))incoming) {
                Curve copy = line.duplicate;
                copy.setFinal(p);
                res.erase_last();
                res.append(copy);
            } else {
                res.appendNew!LineSegment(p);
            }
        }
    }
    res.appendNew!LineSegment(outgoing.initialPoint());
}

void join_inside(ref Path res, in Curve outgoing)
{
    res.appendNew!LineSegment(outgoing.initialPoint());
}

void outline_helper(ref Path res, in Path to_add, Coord width, Coord miter, JoinType join)
{
    Point tang1 = -res.back.reverse.unitTangentAt(0);
    Point tang2 = to_add[0].unitTangentAt(0);
    Point discontinuity_vec = to_add.initialPoint - res.finalPoint;
    bool on_outside = (tang1.dot(discontinuity_vec) >= 0);

    if (on_outside) {
        join_func jf;
        switch (join) {
            case JoinType.JOIN_BEVEL:
                jf = &bevel_join;
                break;
            case JoinType.JOIN_ROUND:
                jf = &round_join;
                break;
            default:
                jf = &miter_join;
        }
        jf(res, to_add[0], width, miter);
    } else {
        join_inside(res, to_add[0]);
    }

    res.append(to_add);
}

/**
 * Offset the input path by @a width.
 * Joins may behave oddly if the width is negative.
 *
 * @param input
 * @param width Amount to offset.
 * @param miter Miter limit. Only used with JOIN_EXTRAPOLATE and JOIN_MITER.
 * @param join
 * @param butt
 */
Path half_outline(in Path input, Coord width, Coord miter, JoinType join = JoinType.JOIN_BEVEL)
{
    Path res = new Path;
    if (input.size == 0) return res;

    Point tang1 = input[0].unitTangentAt(0);
    Point start = input.initialPoint + tang1 * width;

    res.start(start);

    // Do two curves at a time for efficiency, since the join function needs to know the outgoing curve as well
    const size_t k = input.size();
    for (size_t u = 0; u < k; u += 2) {
        Path temp = new Path;

        temp.offset_curve(input[u], width, miter);

        // on the first run through, there isn't a join
        if (u == 0) {
            res.append(temp);
        } else {
            outline_helper(res, temp, width, miter, join);
        }

        // odd number of paths
        if (u < k - 1) {
            temp = new Path;
            temp.offset_curve(input[u+1], width, miter);
            outline_helper(res, temp, width, miter, join);
        }
    }

    if (input.closed) {
        // handling these is so much fun...
        Curve c1 = null;
        Curve c2 = null;

        // TODO
    }

    return res;
}

void offset_curve(ref Path res, in Curve current, Coord width, Coord miter)
{
    const(Coord) tolerance = 0.0025;
    const(size_t) levels = 8;

    // not really a clean way to do this
    if (auto line = cast(const(LineSegment))current) {
        res.append(offset_line(line, width));
    } else if (auto quad = cast(const(QuadraticBezier))current) {
        res.offset_quadratic(quad, width, tolerance, levels);
    } else if (auto cub = cast(const(CubicBezier))current) {
        res.offset_cubic(cub, width, tolerance, levels);
    } else if (auto ell = cast(const(EllipticalArc))current) {
        EllipticalArc e = ell.duplicate();
        Point rays = e.rays();
        e.setRays(rays[X]+width, rays[Y]+width);
    } else {
        // cheat and use the SBasis representation
        D2!SBasis sb = current.toSBasis();
        Point[] temp;
        sbasis_to_bezier(temp, sb, 4);
        auto b = new CubicBezier(temp);
        res.offset_cubic(b, width, tolerance, levels);
    }
}

/// Offsetting a line segment is mathematically stable
LineSegment offset_line(in LineSegment l, Coord width)
{
    Point tang1 = l.unitTangentAt(0).rot90;
    Point tang2 = l.reverse.unitTangentAt(0).rot90;

    Point start = l.initialPoint + tang1 * width;
    Point end = l.finalPoint - tang2 * width;
    
    return new LineSegment(start, end);
}

/// Offsetting a cubic bezier is not
void offset_cubic(ref Path p, in CubicBezier bez, Coord width, Coord tol, size_t levels)
{
    import std.math : fabs;

    Point start_pos = bez.initialPoint();
    Point end_pos = bez.finalPoint();

    Point start_normal = bez.unitTangentAt(0).rot90;
    Point end_normal = -bez.reverse.unitTangentAt(0).rot90;

    // offset the start and end control points out by the width
    Point start_new = start_pos + start_normal*width;
    Point end_new = end_pos + end_normal*width;

    // --------
    Coord start_rad, end_rad;
    Coord start_len, end_len; // tangent lengths
    bez.get_cubic_data(0, start_len, start_rad);
    bez.get_cubic_data(1, end_len, end_rad);

    double start_off = 1, end_off = 1;
    // correction of the lengths of the tangent to the offset
    if (!are_near(start_rad, 0))
	    start_off += width / start_rad;
    if (!are_near(end_rad, 0))
	    end_off += width / end_rad;
    start_off *= start_len;
    end_off *= end_len;
    // --------

    Point mid1_new = start_normal.ccw*start_off;
    mid1_new = Point(start_new[X] + mid1_new[X]/3., start_new[Y] + mid1_new[Y]/3.);
    Point mid2_new = end_normal.ccw*end_off;
    mid2_new = Point(end_new[X] - mid2_new[X]/3., end_new[Y] - mid2_new[Y]/3.);

    // create the estimate curve
    auto c = new CubicBezier(start_new, mid1_new, mid2_new, end_new);

    // reached maximum recursive depth
    // don't bother with any more correction
    if (levels == 0) {
        p.append(c);
        return;
    }

    // check the tolerance for our estimate to be a parallel curve
    Point chk = c.pointAt(.5);
    Point req = bez.pointAt(.5) + bez.unitTangentAt(.5).rot90*width; // required accuracy

    const(Point) diff = req - chk;
    const(Coord) err = diff.dot(diff);

    if (err < tol) {
        // we're good, curve is accurate enough
        p.append(c);
        return;
    } else {
        // split the curve in two
        CubicBezier[2] s = bez.subdivide(.5);
        offset_cubic(p, s[X], width, tol, levels - 1);
        offset_cubic(p, s[Y], width, tol, levels - 1);
    }
}

void get_cubic_data(in CubicBezier bez, double time, ref double len, ref double rad)
{
    // get derivatives
    Point[] derivs = bez.pointAndDerivatives(time, 3);

    Point der1 = derivs[1]; // first deriv (tangent vector)
    Point der2 = derivs[2]; // second deriv (tangent's tangent)
    Coord l = der1.L2; // length

    len = rad = 0;

    if (are_near(l, 0, 1e-4)) {
        l = der2.L2;
        Point der3 = derivs[3]; // third deriv
        if (are_near(l, 0, 1e-4)) {
            l = der3.L2;
            if (are_near(l, 0)) {
                return; // this isn't a segment...
            }
	    rad = 1e8;
        } else {
            rad = -l * (der2.dot(der2) / der3.cross(der2));
        }
    } else {
        rad = -l * (der1.dot(der1) / der2.cross(der1));
    }
    len = l;
}

void offset_quadratic(ref Path p, in QuadraticBezier bez, Coord width, Coord tol, size_t levels)
{
    // cheat
    // it's extra code I saved...
    CubicBezier cub = bez.elevate_degree();
    p.offset_cubic(cub, width, tol, levels);
}

unittest
{

    Path p = new Path;
    p.start(Point(10,30));
    p.appendNew!LineSegment(Point(60,90));
    p.appendNew!LineSegment(Point(130,54));
    p.appendNew!CubicBezier(Point(150,55), Point(160,71), Point(150,92));
    Path out_p = outline(p, 1, 5, JoinType.JOIN_BEVEL, ButtType.BUTT_FLAT)[0];
    writeln("--------");
    foreach (i; 0 .. out_p.size) {
        const(BezierCurve) b = cast(const(BezierCurve))out_p[i];
        Point[] pts = b.points;
        foreach (pt; pts)
            writeln(pt);
        writeln("--");
    }
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
