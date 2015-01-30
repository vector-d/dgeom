/*
 * Affine transformation classes
 *
 * Authors:
 *   ? <?@?.?>
 *   Krzysztof Kosiński <tweenk.pl@gmail.com>
 *   Johan Engelen
 *   Liam P. White
 * 
 * Copyright ?-2015 Authors
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

module geom.transforms;

import math = std.math;
import geom.affine;
import geom.point;

/** Translation by a vector.
 * @ingroup Transforms */
struct Translate
{
    @disable this();

    /* Get a translation that doesn't do anything. */
    static Translate identity() { return Translate(0,0); }

    /* Construct a translation from its vector. */
    this(in Point p) { vec = p; }
    /* Construct a translation from its coordinates. */
    this(in Coord x, in Coord y) { vec = [x, y]; }

    ref inout(Coord) opIndex(size_t i) inout { return vec[i]; }

    Point vector() const { return vec; }
    /* Get the inverse translation. */
    Translate inverse() const { return Translate(-vec); }

    Translate opBinary(string op)(Translate rhs) const if (op == "*")
    { return Translate(this.vec + rhs.vec); }
    mixin self_assign;

    @property Affine toAffine() const
    { return Affine(1, 0, 0, 1, vec[X], vec[Y]); }
    
    alias toAffine this;

    private Point vec;
}

bool are_near(in Translate a, in Translate b, Coord eps = EPSILON)
{ return geom.coord.are_near(a[X], b[X], eps) && geom.coord.are_near(a[Y], b[Y], eps); }

/** Scaling from the origin.
 * During scaling, the point (0,0) will not move. To obtain a scale  with a different
 * invariant point, combine with translation to the origin and back.
 * @ingroup Transforms */
struct Scale
{
    @disable this();

    /* Get a scaling that doesn't do anything. */
    static Scale identity() { return Scale(1, 1); }

    /* Create a scaling from two scaling factors given as coordinates of a point. */
    this(in Point p) { vec = p; }
    /* Create a scaling from two scaling factors. */
    this(in Coord x, in Coord y) { vec = [x, y]; }
    /* Create an uniform scaling from a single scaling factor. */
    this(in Coord s) { vec = [s, s]; }

    // TODO: should we keep the mutators?
    ref inout(Coord) opIndex(size_t i) inout { return vec[i]; }

    Point vector() const { return vec; }
    Scale inverse() const { return Scale(1./vec[X], 1./vec[Y]); }

    Scale opBinary(string op)(Scale rhs) const if (op == "*")
    { return Scale(vec[X] * rhs[X], vec[Y] * rhs[Y]); }
    mixin self_assign;

    @property Affine toAffine() const
    { return Affine(vec[X], 0, 0, vec[Y], 0, 0); }
    
    alias toAffine this;

    private Point vec;
}

bool are_near(in Scale a, in Scale b, Coord eps = EPSILON)
{ return geom.coord.are_near(a[X], b[X], eps) && geom.coord.are_near(a[Y], b[Y], eps); }

/** Rotation around the origin.
 * Combine with translations to the origin and back to get a rotation around a different point.
 * @ingroup Transforms */
struct Rotate
{
    @disable this();

    /* Get a zero-degree rotation. */
    static Rotate identity() { return Rotate(1, 0); }

    /** Construct a rotation from its angle in degrees.
     * Positive arguments correspond to clockwise rotations if Y grows downwards. */
    static Rotate from_degrees(Coord deg)
    {
        Coord rad = (deg / 180.0) * math.PI;
        return Rotate(rad);
    }

    /** Construct a rotation from its angle in radians.
     * Positive arguments correspond to counter-clockwise rotations (if Y grows upwards). */
    this(in Coord theta) { vec = Point.polar(theta); }
    /* Construct a rotation from its characteristic vector. */
    this(in Point p) { vec = unit_vector(p); }
    /* Construct a rotation from the coordinates of its characteristic vector. */
    this(in Coord x, in Coord y) { this(Point(x, y)); }

    Rotate inverse() const
    { return Rotate(vec[X], -vec[Y]); }

    ref inout(Coord) opIndex(size_t i) inout { return vec[i]; }

    Rotate opBinary(string op)(Rotate rhs) const if (op == "*")
    { return Rotate(vec * cast(Affine)rhs); }
    mixin self_assign;

    @property Affine toAffine() const
    { return Affine(vec[X], vec[Y], -vec[Y], vec[X], 0, 0); }
    
    alias toAffine this;

    private Point vec;
}

bool are_near(in Rotate a, in Rotate b, Coord eps = EPSILON)
{ return geom.coord.are_near(a[X], b[X], eps) && geom.coord.are_near(a[Y], b[Y], eps); }

/** Horizontal shearing.
 * Points on the X axis will not move. Combine with translations to get a shear
 * with a different invariant line.
 * @ingroup Transforms */
struct HShear
{
    @disable this();
    this(Coord _f) { f = _f; }

    Coord factor() const { return f; }
    void setFactor(Coord nf) { f = nf; }
    HShear inverse() const { return HShear(-f); }
    static HShear identity() { return HShear(0); }

    HShear opBinary(string op)(HShear rhs) const if (op == "*")
    { return HShear(f + rhs.f); }
    mixin self_assign;

    @property Affine toAffine() const
    { return Affine(1, 0, f, 1, 0, 0); }

    alias toAffine this;

    private Coord f;
}

bool are_near(in HShear a, in HShear b, Coord eps = EPSILON)
{ return geom.coord.are_near(a.factor(), b.factor(), eps); }

/** Vertical shearing.
 * Points on the Y axis will not move. Combine with translations to get a shear
 * with a different invariant line.
 * @ingroup Transforms */
struct VShear
{
    @disable this();
    this(Coord _f) { f = _f; }

    Coord factor() const { return f; }
    void setFactor(Coord nf) { f = nf; }
    VShear inverse() const { return VShear(-f); }
    static VShear identity() { return VShear(0); }

    HShear opBinary(string op)(VShear rhs) const if (op == "*")
    { return VShear(f + rhs.f); }
    mixin self_assign;

    @property Affine toAffine() const
    { return Affine(1, f, 0, 1, 0, 0); }
    
    alias toAffine this;

    private Coord f;
}

bool are_near(in VShear a, in VShear b, Coord eps = EPSILON)
{ return geom.coord.are_near(a.factor(), b.factor(), eps); }


/** Combination of a translation and uniform scale.
 * The translation part is applied first, then the result is scaled from the new origin.
 * This way when the class is used to accumulate a zoom transform, trans always points
 * to the new origin in original coordinates.
 * @ingroup Transforms */
struct Zoom
{
    @disable this();
    /* Get a zoom that doesn't do anything. */
    static Zoom identity() { return Zoom(1); }

    /* Construct a zoom from a scaling factor. */
    this(in Coord s) { _scale = s; _trans = [0,0]; }
    /* Construct a zoom from a translation. */
    this(in Translate t) { _scale = 1; _trans = t.vector(); }
    /* Construct a zoom from a scaling factor and a translation. */
    this(in Coord s, in Translate t) { _scale = s; _trans = t.vector(); }

    Coord scale() const { return _scale; }
    void setScale(in Coord s) { _scale = s; }
    Point translation() const { return _trans; }
    void setTranslation(in Point p) { _trans = p; }

    Zoom inverse() const { return Zoom(1/_scale, Translate(-_trans*_scale)); }

    /** Zoom between rectangles.
     * Given two rectangles, compute a zoom that maps one to the other.
     * Rectangles are assumed to have the same aspect ratio. */
    import geom.rect;

    Zoom map_rect(in Rect old_r, in Rect new_r)
    {
        Zoom ret = Zoom(1);
        ret._scale = new_r.width() / old_r.width();
        ret._trans = new_r.min() - old_r.min();
        return ret;
    }

    Zoom opBinary(string op)(Zoom rhs) const if (op == "*")
    { return Zoom(_scale * rhs._scale, Translate(_trans + rhs._trans / _scale)); }
    mixin self_assign;

    @property Affine toAffine() const
    { return Affine(_scale, 0, 0, _scale, _trans[X] * _scale, _trans[Y] * _scale); }

    alias toAffine this;

private:
    Coord _scale;
    Point _trans;
}

/** Reflects objects about line.
 * The line, defined by a vector along the line and a point on it, acts as a mirror. */
Affine reflection(in Point vector, in Point origin)
{
    Point vec_norm = unit_vector(vector);
    Affine mirror = [vec_norm[X]*vec_norm[X] - vec_norm[Y]*vec_norm[Y], 2 * vec_norm[X] * vec_norm[Y] ,
                    2 * vec_norm[X] * vec_norm[Y], vec_norm[Y]*vec_norm[Y] - vec_norm[X]*vec_norm[X] ,
                    0 ,0];
    return Translate(-origin) * mirror * Translate(origin);
}

private mixin template self_assign()
{
    void opOpAssign(string op, T)(T rhs)
    { mixin("this = this "~op~" rhs;"); }
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
