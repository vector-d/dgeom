/**
 * Bezier Fitting
 *
 * Authors:
 *    Parcly Taxel (Jeremy Tan)
 *    Mihail K.
 *
 * Copyright (C) 2015 Authors
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

module geom.bezier_fitting;

import std.algorithm;
import std.exception : enforce;
import std.math;
import std.range;
import std.typecons;

import geom.bezier_curve;
import geom.coord;
import geom.point;
import geom.transforms;

/++
 + Helper function for bezierFit.
 + Creates a list of 4 points, representing the current curve.
 +
 + @param window Points to use in the curve calculation.
 + @return A list of points for the current curve.
 +/
Point[4] createCurve(Point[] window)
{
    return [
        window[0],
        window[0].lerp(window[1], 1.0 / 3.0),
        window[0].lerp(window[1], 2.0 / 3.0),
        window[1]
    ];
}

/++
 + Fit an array of points to Bezier(s).
 +
 + @pre At least two points must be given.
 + @param points Input slice of points to fit to a bezier
 + @param closed If true, smooth for a closed set of points
 +/
Point[] bezierFit(Point[] points, bool closed = false)
{
    // Curves need at least 2 points
    if(points.length < 2) {
        enforce("Curve needs at least two points!");
    }

    Point[] res;
    Point[] window = points[0 .. 2];
    int prevTrail = 0, trail = 0, lead = 1;

    // Current curve
    Point[4] curCurve = window.createCurve;

    bool maybeOver = false;
    while(lead + 1 < points.length) {
        lead++;
        window = points[trail .. lead + 1]; // Extend the window one more node
        Point v = window[$ - 3] - window[$ - 2];
        Point w = window[$ - 1] - window[$ - 2];

        Point[] newCurve;
        // 60 degrees or less
        if (v.dot(w) / v.L2 / w.L2 >= 0.5) {
            if (maybeOver) {
                newCurve = points[prevTrail .. lead].stress()[0];
                res[$ - 3 .. $] = newCurve[1 .. $];
                trail = lead - 1;
                maybeOver = false;
            } else {
                if(res.length == 0)
                    res ~= curCurve[0 .. 1];
                res ~= curCurve[1 .. $];
                prevTrail = trail;
            }
            trail = lead - 1;
            window = points[trail .. lead + 1];
            curCurve = window.createCurve;
        } else {
            bool over = false;
            if (window.length == 3) {
                Coord t = window.chords()[1];
                Point[] qcurve = [
                    window[0],
                    (window[1] - window[0] * (1 - t) ^^ 2 -
                            window[2] * t ^^ 2) / (2 * t * (1 - t)),
                    window[2]
                ];
                newCurve = [
                    qcurve[0],
                    qcurve[0].lerp(qcurve[1], 2.0 / 3.0),
                    qcurve[1].lerp(qcurve[2], 1.0 / 3.0),
                    qcurve[2]
                ];
            } else if(window.length == 4) {
                Point[4] tmp = window;
                newCurve = tmp.cubicFrom4;
            } else {
                auto product = window.stress;
                Coord shortSeg = iota(0, window.length - 1)
                        .map!(i => window[i].distance(window[i + 1])).reduce!min;

                // Stop condition: maximum error > 1 / 3 * minimum segment length
                if(product[1].reduce!max > 0.33 * shortSeg) over = true;
                else newCurve = product[0];
            }

            if(over) {
                maybeOver = true;
                if(res.length == 0)
                    res ~= curCurve[0 .. 1];
                res ~= curCurve[1 .. $];
                prevTrail = trail;
                trail = lead - 1;
                window = points[trail .. lead + 1];
                curCurve = window.createCurve;
            } else {
                // Update current state
                curCurve = newCurve;
                maybeOver = false;
            }
        }
    }

    if(maybeOver) {
        auto newCurve = points[prevTrail .. lead + 1].stress[0];
        res[$ - 3 .. $] = newCurve[1 .. $];
    } else {
        if(res.length == 0)
            res ~= curCurve[0 .. 1];
        res ~= curCurve[1 .. $];
    }

    // Remove the last element
    Point ouro = res.back;
    res.popBack;

    foreach (t; iota(0, res.length, 3)) {
        if (t != 0 || closed) {
            Point v = res[t - 1] - res[t];
            Point w = res[t + 1] - res[t];
            real angle = v.dot(w) / v.L2 / w.L2;

            // Cos 160
            if(angle <= -0.94) {
                real theta = (PI - acos(angle)) / 2;
                real sign = (abs(v.dirc - w.dirc) >= PI) ^
                        (v.dirc > w.dirc) ? 1 : -1;

                res[t - 1] = res[t] + v * Rotate(sign * theta);
                res[t + 1] = res[t] + v * Rotate(-sign * theta);
            }
        }
    }

    // Restore the last element
    return res ~= ouro;
}

real dirc(Point point0, Point point1 = Point())
{
    return atan2(
        point1.y - point0.y,
        point1.x - point0.x
    );
}

Point bezierPointAtT(Point[4] cubic, Coord t)
{
    auto b = bezierFromPoints(cubic);
    return b.pointAt(t);
}

/++
 + This function takes in a list of points and returns
 + a list of numbers between 0 and 1 corresponding to the relative positions
 + of said points (assuming consecutive points are linked by straight lines).
 + The first item is always 0.0 and the last one 1.0.
 +
 + @param points The list of points used in the chords calculation.
 + @return A list of numbers corresponding to relative positions of the points.
 +/
Coord[] chords(Point points[])
{
    if(points.length < 1) {
        return [0.0, 1.0];
    }

    Coord[] lengths, ratios;
    foreach(i; 0 .. points.length - 1)
        lengths ~= points[i + 1].distance(points[i]);

    ratios ~= 0.0;

    foreach(i; 0 .. lengths.length)
        ratios ~= lengths[0 .. i + 1].reduce!"a + b" / lengths.reduce!"a + b";
    ratios ~= 1.0;

    return ratios;
}

/++
 + Takes a list of four points and generates a curve passing
 + through all based on chords().
 +
 + @param points The points used in the calculation.
 + @param p Value for lambda.
 + @param q Value for mu.
 + @return A list of points that make up a curve.
 +/
Point[] cubicFrom4(Point points[4], Coord p, Coord q)
{
    Point l, m;
    Coord a, b, c, d;

    a = 3 * p * (1 - p) ^^ 2;
    b = 3 * (1 - p) * p ^^ 2;

    c = 3 * q * (1 - q) ^^ 2;
    d = 3 * (1 - q) * q ^^ 2;

    Point x = points[1] - points[0] * (1 - p) ^^ 3 - points[3] * p ^^ 3;
    Point y = points[2] - points[0] * (1 - q) ^^ 3 - points[3] * q ^^ 3;

    Coord det = a * d - b * c;
    l = (x * d - y * b) / det;
    m = (y * a - x * c) / det;

    return [points[0], l, m, points[3]];
}

/++
 + Convenience function for cubicFrom4.
 + The values of p and q are calculated automatically.
 +
 + @param points The points used in the calculation.
 + @return A list of points that make up a curve.
 +/
Point[] cubicFrom4(Point[4] points)
{
    Coord[] tmp = points.chords;
    return cubicFrom4(points, tmp[1], tmp[2]);
}

Tuple!(Point[4], Coord[]) stress(Point[] points)
{
    int middle = cast(int)(points.length / 2);
    Coord[] callipers = points.chords;

    Point[][] seeds;
    foreach(i; 1 .. middle) {
        seeds ~= cubicFrom4(
            [points[0], points[i], points[$ - i - 1], points[$ - 1]],
            callipers[i], callipers[$ - i - 1]
        );
    }

    Point a, b;
    foreach(seed; seeds) {
        a += seed[1];
        b += seed[2];
    }

    Coord[] errors;
    Point[4] curve = [
        points[0], a / seeds.length,
        b / seeds.length, points[$ - 1]
    ];

    // Initialize errors array
    errors.length = curve.length;

    // Refine projection
    foreach(i; 0 .. 5) {
        foreach(j; iota(0, middle - 1).retro) {
            Point delta1 = curve.project(points[j]);
            Point delta2 = curve.project(points[$ - j - 1]);

            curve[1] += delta1 * 2.5;
            curve[2] += delta2 * 2.5;

            errors[j] = delta1.L2;
            errors[$ - j - 1] = delta2.L2;
        }
    }

    if(points.length % 2)
        errors[middle] = curve.project(points[middle]).L2;
    return tuple(curve, errors);
}

/++
 + Projects a node onto a cubic curve.
 +
 + @param curve A list of 4 points that make up the curve.
 + @param point The point being projected.
 + @return The result of the projection.
 +/
Point project(Point[4] curve, Point point)
{
    Coord[] lookup;
    real samples = 200;
    real width = 1 / samples;

    foreach(i; 0 .. samples + 1) {
        lookup ~= curve.bezierPointAtT(i / samples).distance(point);
    }

    Coord mindist = lookup.reduce!(min);
    Coord t = lookup.countUntil(mindist) / samples;

    // TODO : Is this really the most effective
    // way of going about this?
    while(width > 1.0e-5) {
        Coord left  = curve.bezierPointAtT(max(t - width, 0)).distance(point);
        Coord right = curve.bezierPointAtT(min(t + width, 1)).distance(point);

        if(t == 0.0) left  = mindist + 1;
        if(t == 1.0) right = mindist + 1;

        if(left < mindist || right < mindist) {
            mindist = min(left, right);
            t = left < right ? max(t - width, 0) : min(t + width, 1);
        } else {
            width /= 2;
        }
    }

    // Apply the projection
    return point - curve.bezierPointAtT(t);
}

unittest
{
    Point[] input = [ Point(20.,40), Point(50.,100), Point(120.,110), Point(160.,60) ];
    Point[] expected = [ Point(20., 40), Point(11.8225095793, 120.582109841), Point(148.774429751, 145.987164121), Point(160., 60) ];
    Point[] output = bezierFit(input);
    foreach (i; 0 .. expected.length)
        assert(are_near(expected[i], output[i]));

    input = [ Point(20,180), Point(20,220), Point(50,240), Point(90,250) ];
    expected = [ Point(20,180), Point(-6.39036470691,230.093020882), Point(58.9142504083,245.099863871), Point(90,250) ];
    output = bezierFit(input);
    foreach (i; 0 .. expected.length)
        assert(are_near(expected[i], output[i]));
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
