/** 
 * 3x3 affine transformation matrix.
 *
 * Authors:
 *   Lauris Kaplinski <lauris@kaplinski.com> (Original NRAffine definition and related macros)
 *   Nathan Hurst <njh@mail.csse.monash.edu.au> (Geom::Affine class version of the above)
 *   Michael G. Sloan <mgsloan@gmail.com> (reorganization and additions)
 *   Krzysztof Kosiński <tweenk.pl@gmail.com> (removal of boilerplate, docs)
 *
 * This code is in public domain.
 */

module geom.affine;

public import geom.coord;
import math = std.math;
import geom.point;
import geom.transforms;

/**
 * 3x3 matrix representing an affine transformation.
 *
 * Affine transformations on elements of a vector space are transformations which can be
 * expressed in terms of matrix multiplication followed by addition
 * (\f$x \mapsto A x + b\f$). They can be thought of as generalizations of linear functions
 * (\f$y = a x + b\f$) to vector spaces. Affine transformations of points on a 2D plane preserve
 * the following properties:
 *
 * - Colinearity: if three points lie on the same line, they will still be colinear after
 *   an affine transformation.
 * - Ratios of distances between points on the same line are preserved
 * - Parallel lines remain parallel.
 *
 * All affine transformations on 2D points can be written as combinations of scaling, rotation,
 * shearing and translation. They can be represented as a combination of a vector and a 2x2 matrix,
 * but this form is inconvenient to work with. A better solution is to represent all affine
 * transformations on the 2D plane as 3x3 matrices, where the last column has fixed values.
 * \f[ A = \left[ \begin{array}{ccc}
    c_0 & c_1 & 0 \\
    c_2 & c_3 & 0 \\
    c_4 & c_5 & 1 \end{array} \right]\f]
 *
 * We then interpret points as row vectors of the form \f$[p_X, p_Y, 1]\f$. Applying a
 * transformation to a point can be written as right-multiplication by a 3x3 matrix
 * (\f$p' = pA\f$). This subset of matrices is closed under multiplication - combination
 * of any two transforms can be expressed as the multiplication of their matrices.
 * In this representation, the \f$c_4\f$ and \f$c_5\f$ coefficients represent
 * the translation component of the transformation.
 *
 * Matrices can be multiplied by other more specific transformations. When multiplying,
 * the transformations are applied from left to right, so for example <code>m = a * b</code>
 * means: @a m first transforms by a, then by b.
 *
 * @ingroup Transforms
 */
struct Affine
{
    @disable this();

    /** Create a matrix from its coefficient values.
     * It's rather inconvenient to directly create matrices in this way. Use transform classes
     * if your transformation has a geometric interpretation.
     * @see Translate
     * @see Scale
     * @see Rotate
     * @see HShear
     * @see VShear
     * @see Zoom */
    this(Coord c0, Coord c1, Coord c2, Coord c3, Coord c4, Coord c5)
    { _c = [c0, c1, c2, c3, c4, c5]; }
    
    this(const(Coord[6]) arr)
    { _c = arr; }

    /** Create an identity matrix.
     * @return The matrix
     *         \f$\left[\begin{array}{ccc}
               1 & 0 & 0 \\
               0 & 1 & 0 \\
               0 & 0 & 1 \end{array}\right]\f$.
     * @relates Affine */
    static Affine identity()
    { return Affine([1., 0, 0, 1, 0, 0]); }

    /** Access a coefficient by its index. */
    ref inout(Coord) opIndex(size_t i) inout
    { return _c[i]; }

    /+ Transform an Affine +/

    /** Combine this transformation with another one.
     * After this operation, the matrix will correspond to the transformation
     * obtained by first applying the original version of this matrix, and then
     * applying @a m. */
    Affine opBinary(string op)(in Affine rhs) const if (op == "*")
    {
        Affine ret = this;
        Coord nc[6];
        for(int a = 0; a < 5; a += 2) {
            for(int b = 0; b < 2; b++) {
                nc[a + b] = _c[a] * rhs._c[b] + _c[a + 1] * rhs._c[b + 2];
            }
        }
        ret._c = nc;
        ret[4] += rhs._c[4];
        ret[5] += rhs._c[5];
        return ret;
    }

    /+ Get the parameters of the matrix's transform +/

    Point xAxis() const
    { return Point(_c[0], _c[1]); }

    Point yAxis() const
    { return Point(_c[2], _c[3]); }

    /** Gets the translation imparted by the Affine. */
    Point translation() const
    { return Point(_c[4], _c[5]); }

    /** Calculates the amount of x-scaling imparted by the Affine.  This is the scaling applied to
     *  the original x-axis region.  It is \emph{not} the overall x-scaling of the transformation.
     *  Equivalent to L2(m.xAxis())
     */
    Coord expansionX() const
    { return math.sqrt(_c[0] * _c[0] + _c[1] * _c[1]); }

    /** Calculates the amount of y-scaling imparted by the Affine.  This is the scaling applied before
     *  the other transformations.  It is \emph{not} the overall y-scaling of the transformation. 
     *  Equivalent to L2(m.yAxis())
     */
    Coord expansionY() const
    { return math.sqrt(_c[2]*_c[2] + _c[3]*_c[3]); }
    
    Point expansion() const { return Point(expansionX(), expansionY()); }


    /+ Modify the matrix +/

    void setXAxis(in Point vec)
    { _c[0..2] = [vec[0], vec[1]]; }

    void setYAxis(in Point vec)
    { _c[2..4] = [vec[0], vec[1]]; }

    /** Sets the translation imparted by the Affine */
    void setTranslation(in Point loc)
    { _c[4..6] = [loc[0], loc[1]]; }

    void setExpansionX(Coord val)
    {
        auto exp_x = expansionX();
        if(!geom.coord.are_near(exp_x, 0.0)) {  // TODO: best way to deal with it is to skip op?
            auto coef = val / expansionX();
            for(uint i=0;i<2;i++) _c[i] *= coef;
        }
    }

    void setExpansionY(Coord val)
    {
        auto exp_y = expansionY();
        if(!geom.coord.are_near(exp_y, 0.0)) {  // TODO: best way to deal with it is to skip op?
            auto coef = val / expansionY();
            for(uint i=2; i<4; i++) _c[i] *= coef;
        }
    }

    /** Sets this matrix to be the Identity Affine. */
    void setIdentity()
    { _c = [1, 0, 0, 1, 0, 0]; }


    /+ Inspect the matrix's transform +/

    /** Check whether this matrix is an identity matrix.
     * @param eps Numerical tolerance
     * @return True iff the matrix is of the form
     *         \f$\left[\begin{array}{ccc}
               1 & 0 & 0 \\
               0 & 1 & 0 \\
               0 & 0 & 1 \end{array}\right]\f$ */
    bool isIdentity(Coord eps = EPSILON) const
    { return are_near(this, identity(), eps); }

    /** Check whether this matrix represents a pure translation.
     * Will return true for the identity matrix, which represents a zero translation.
     * @param eps Numerical tolerance
     * @return True iff the matrix is of the form
     *         \f$\left[\begin{array}{ccc}
               1 & 0 & 0 \\
               0 & 1 & 0 \\
               a & b & 1 \end{array}\right]\f$ */
    bool isTranslation(Coord eps = EPSILON) const
    {
        return geom.coord.are_near(_c[0], 1.0, eps) && geom.coord.are_near(_c[1], 0.0, eps) &&
               geom.coord.are_near(_c[2], 0.0, eps) && geom.coord.are_near(_c[3], 1.0, eps);
    }

    /** Check whether this matrix represents pure scaling.
     * @param eps Numerical tolerance
     * @return True iff the matrix is of the form
     *         \f$\left[\begin{array}{ccc}
               a & 0 & 0 \\
               0 & b & 0 \\
               0 & 0 & 1 \end{array}\right]\f$. */
    bool isScale(Coord eps = EPSILON) const
    {
        if (isSingular(eps)) return false;
        return geom.coord.are_near(_c[1], 0.0, eps) && geom.coord.are_near(_c[2], 0.0, eps) && 
               geom.coord.are_near(_c[4], 0.0, eps) && geom.coord.are_near(_c[5], 0.0, eps);
    }

    /** Check whether this matrix represents pure uniform scaling.
     * @param eps Numerical tolerance
     * @return True iff the matrix is of the form
     *         \f$\left[\begin{array}{ccc}
               a_1 & 0 & 0 \\
               0 & a_2 & 0 \\
               0 & 0 & 1 \end{array}\right]\f$ where \f$|a_1| = |a_2|\f$. */
    bool isUniformScale(Coord eps = EPSILON) const
    {
        if (isSingular(eps)) return false;
        return geom.coord.are_near(math.fabs(_c[0]), math.fabs(_c[3]), eps) &&
               geom.coord.are_near(_c[1], 0.0, eps) && geom.coord.are_near(_c[2], 0.0, eps) &&  
               geom.coord.are_near(_c[4], 0.0, eps) && geom.coord.are_near(_c[5], 0.0, eps);
    }

    /** Check whether this matrix represents a pure rotation.
     * @param eps Numerical tolerance
     * @return True iff the matrix is of the form
     *         \f$\left[\begin{array}{ccc}
               a & b & 0 \\
               -b & a & 0 \\
               0 & 0 & 1 \end{array}\right]\f$ and \f$a^2 + b^2 = 1\f$. */
    bool isRotation(Coord eps = EPSILON) const
    {
        return geom.coord.are_near(_c[0], _c[3], eps) && geom.coord.are_near(_c[1], -_c[2], eps) &&
               geom.coord.are_near(_c[4], 0.0, eps) && geom.coord.are_near(_c[5], 0.0, eps) &&
               geom.coord.are_near(_c[0]*_c[0] + _c[1]*_c[1], 1.0, eps);
    }

    /** Check whether this matrix represents pure horizontal shearing.
     * @param eps Numerical tolerance
     * @return True iff the matrix is of the form
     *         \f$\left[\begin{array}{ccc}
               1 & 0 & 0 \\
               k & 1 & 0 \\
               0 & 0 & 1 \end{array}\right]\f$. */
    bool isHShear(Coord eps = EPSILON) const
    {
        return geom.coord.are_near(_c[0], 1.0, eps) && geom.coord.are_near(_c[1], 0.0, eps) &&
               geom.coord.are_near(_c[3], 1.0, eps) && geom.coord.are_near(_c[4], 0.0, eps) &&
               geom.coord.are_near(_c[5], 0.0, eps);
    }

    /** Check whether this matrix represents pure vertical shearing.
     * @param eps Numerical tolerance
     * @return True iff the matrix is of the form
     *         \f$\left[\begin{array}{ccc}
               1 & k & 0 \\
               0 & 1 & 0 \\
               0 & 0 & 1 \end{array}\right]\f$. */
    bool isVShear(Coord eps = EPSILON) const
    {
        return geom.coord.are_near(_c[0], 1.0, eps) && geom.coord.are_near(_c[2], 0.0, eps) &&
               geom.coord.are_near(_c[3], 1.0, eps) && geom.coord.are_near(_c[4], 0.0, eps) &&
               geom.coord.are_near(_c[5], 0.0, eps);
    }

    /** Check whether this matrix represents a pure nonzero translation.
     * @param eps Numerical tolerance
     * @return True iff the matrix is of the form
     *         \f$\left[\begin{array}{ccc}
               1 & 0 & 0 \\
               0 & 1 & 0 \\
               a & b & 1 \end{array}\right]\f$ and \f$a, b \neq 0\f$ */
    bool isNonzeroTranslation(Coord eps = EPSILON) const
    {
        return geom.coord.are_near(_c[0], 1.0, eps) && geom.coord.are_near(_c[1], 0.0, eps) &&
               geom.coord.are_near(_c[2], 0.0, eps) && geom.coord.are_near(_c[3], 1.0, eps) &&
               (!geom.coord.are_near(_c[4], 0.0, eps) || !geom.coord.are_near(_c[5], 0.0, eps));
    }

    /** Check whether this matrix represents pure, nonzero scaling.
     * @param eps Numerical tolerance
     * @return True iff the matrix is of the form
     *         \f$\left[\begin{array}{ccc}
               a & 0 & 0 \\
               0 & b & 0 \\
               0 & 0 & 1 \end{array}\right]\f$ and \f$a, b \neq 1\f$. */
    bool isNonzeroScale(Coord eps = EPSILON) const
    {
        if (isSingular(eps)) return false;
        return (!geom.coord.are_near(_c[0], 1.0, eps) || !geom.coord.are_near(_c[3], 1.0, eps)) &&  //NOTE: these are the diags, and the next line opposite diags
               geom.coord.are_near(_c[1], 0.0, eps) && geom.coord.are_near(_c[2], 0.0, eps) && 
               geom.coord.are_near(_c[4], 0.0, eps) && geom.coord.are_near(_c[5], 0.0, eps);
    }

    /** Check whether this matrix represents pure, nonzero uniform scaling.
     * @param eps Numerical tolerance
     * @return True iff the matrix is of the form
     *         \f$\left[\begin{array}{ccc}
               a_1 & 0 & 0 \\
               0 & a_2 & 0 \\
               0 & 0 & 1 \end{array}\right]\f$ where \f$|a_1| = |a_2|\f$
     * and \f$a_1, a_2 \neq 1\f$. */
    bool isNonzeroUniformScale(Coord eps = EPSILON) const
    {
        if (isSingular(eps)) return false;
        // we need to test both c0 and c3 to handle the case of flips,
        // which should be treated as nonzero uniform scales
        return !(geom.coord.are_near(_c[0], 1.0, eps) && geom.coord.are_near(_c[3], 1.0, eps)) &&
               geom.coord.are_near(math.fabs(_c[0]), math.fabs(_c[3]), eps) &&
               geom.coord.are_near(_c[1], 0.0, eps) && geom.coord.are_near(_c[2], 0.0, eps) &&  
               geom.coord.are_near(_c[4], 0.0, eps) && geom.coord.are_near(_c[5], 0.0, eps);
    }

    /** Check whether this matrix represents a pure, nonzero rotation.
     * @param eps Numerical tolerance
     * @return True iff the matrix is of the form
     *         \f$\left[\begin{array}{ccc}
               a & b & 0 \\
               -b & a & 0 \\
               0 & 0 & 1 \end{array}\right]\f$, \f$a^2 + b^2 = 1\f$ and \f$a \neq 1\f$. */
    bool isNonzeroRotation(Coord eps = EPSILON) const
    {
        return !geom.coord.are_near(_c[0], 1.0, eps) &&
               geom.coord.are_near(_c[0], _c[3], eps) && geom.coord.are_near(_c[1], -_c[2], eps) &&
               geom.coord.are_near(_c[4], 0.0, eps) && geom.coord.are_near(_c[5], 0.0, eps) &&
               geom.coord.are_near(_c[0]*_c[0] + _c[1]*_c[1], 1.0, eps);
    }

    /** Check whether this matrix represents a non-zero rotation about any point.
     * @param eps Numerical tolerance
     * @return True iff the matrix is of the form
     *         \f$\left[\begin{array}{ccc}
               a & b & 0 \\
               -b & a & 0 \\
               c & d & 1 \end{array}\right]\f$, \f$a^2 + b^2 = 1\f$ and \f$a \neq 1\f$. */
    bool isNonzeroNonpureRotation(Coord eps = EPSILON) const
    {
        return !geom.coord.are_near(_c[0], 1.0, eps) &&
               geom.coord.are_near(_c[0], _c[3], eps) && geom.coord.are_near(_c[1], -_c[2], eps) &&
               geom.coord.are_near(_c[0]*_c[0] + _c[1]*_c[1], 1.0, eps);
    }

    /** For a (possibly non-pure) non-zero-rotation matrix, calculate the rotation center.
     * @pre The matrix must be a non-zero-rotation matrix to prevent division by zero, see isNonzeroNonpureRotation().
     * @return The rotation center x, the solution to the equation
     *         \f$A x = x\f$. */
    Point rotationCenter() const
    {
        Coord x = (_c[2]*_c[5]+_c[4]-_c[4]*_c[3]) / (1-_c[3]-_c[0]+_c[0]*_c[3]-_c[2]*_c[1]);
        Coord y = (_c[1]*x + _c[5]) / (1 - _c[3]);
        return Point(x,y);
    }

    /** Check whether this matrix represents pure, nonzero horizontal shearing.
     * @param eps Numerical tolerance
     * @return True iff the matrix is of the form
     *         \f$\left[\begin{array}{ccc}
               1 & 0 & 0 \\
               k & 1 & 0 \\
               0 & 0 & 1 \end{array}\right]\f$ and \f$k \neq 0\f$. */
    bool isNonzeroHShear(Coord eps = EPSILON) const
    {
        return geom.coord.are_near(_c[0], 1.0, eps) && geom.coord.are_near(_c[1], 0.0, eps) &&
              !geom.coord.are_near(_c[2], 0.0, eps) && geom.coord.are_near(_c[3], 1.0, eps) &&
               geom.coord.are_near(_c[4], 0.0, eps) && geom.coord.are_near(_c[5], 0.0, eps);
    }

    /** Check whether this matrix represents pure, nonzero vertical shearing.
     * @param eps Numerical tolerance
     * @return True iff the matrix is of the form
     *         \f$\left[\begin{array}{ccc}
               1 & k & 0 \\
               0 & 1 & 0 \\
               0 & 0 & 1 \end{array}\right]\f$ and \f$k \neq 0\f$. */
    bool isNonzeroVShear(Coord eps = EPSILON) const
    {
        return geom.coord.are_near(_c[0], 1.0, eps) && !geom.coord.are_near(_c[1], 0.0, eps) &&
               geom.coord.are_near(_c[2], 0.0, eps) && geom.coord.are_near(_c[3], 1.0, eps) &&
               geom.coord.are_near(_c[4], 0.0, eps) && geom.coord.are_near(_c[5], 0.0, eps);
    }

    /** Check whether this matrix represents zooming.
     * Zooming is any combination of translation and uniform non-flipping scaling.
     * It preserves angles, ratios of distances between arbitrary points
     * and unit vectors of line segments.
     * @param eps Numerical tolerance
     * @return True iff the matrix is invertible and of the form
     *         \f$\left[\begin{array}{ccc}
               a & 0 & 0 \\
               0 & a & 0 \\
               b & c & 1 \end{array}\right]\f$. */
    bool isZoom(Coord eps = EPSILON) const
    {
        if (isSingular(eps)) return false;
        return geom.coord.are_near(_c[0], _c[3], eps) &&
        geom.coord.are_near(_c[1], 0, eps) && geom.coord.are_near(_c[2], 0, eps);
    }

    /** Check whether the transformation preserves areas of polygons.
     * This means that the transformation can be any combination of translation, rotation,
     * shearing and squeezing (non-uniform scaling such that the absolute value of the product
     * of Y-scale and X-scale is 1).
     * @param eps Numerical tolerance
     * @return True iff \f$|\det A| = 1\f$. */
    bool preservesArea(Coord eps = EPSILON) const
    { return geom.coord.are_near(descrim2(), 1.0, eps); }

    /** Check whether the transformation preserves angles between lines.
     * This means that the transformation can be any combination of translation, uniform scaling,
     * rotation and flipping.
     * @param eps Numerical tolerance
     * @return True iff the matrix is of the form
     *         \f$\left[\begin{array}{ccc}
                 a & b & 0 \\
                 -b & a & 0 \\
                 c & d & 1 \end{array}\right]\f$ or
               \f$\left[\begin{array}{ccc}
                 -a & b & 0 \\
                 b & a & 0 \\
                 c & d & 1 \end{array}\right]\f$. */
    bool preservesAngles(Coord eps = EPSILON) const
    {
        if (isSingular(eps)) return false;
        return (geom.coord.are_near(_c[0], _c[3], eps) && geom.coord.are_near(_c[1], -_c[2], eps)) ||
               (geom.coord.are_near(_c[0], -_c[3], eps) && geom.coord.are_near(_c[1], _c[2], eps));
    }

    /** Check whether the transformation preserves distances between points.
     * This means that the transformation can be any combination of translation,
     * rotation and flipping.
     * @param eps Numerical tolerance
     * @return True iff the matrix is of the form
     *         \f$\left[\begin{array}{ccc}
               a & b & 0 \\
               -b & a & 0 \\
               c & d & 1 \end{array}\right]\f$ or
               \f$\left[\begin{array}{ccc}
               -a & b & 0 \\
               b & a & 0 \\
               c & d & 1 \end{array}\right]\f$ and \f$a^2 + b^2 = 1\f$. */
    bool preservesDistances(Coord eps = EPSILON) const
    {
        return ((geom.coord.are_near(_c[0], _c[3], eps) && geom.coord.are_near(_c[1], -_c[2], eps)) ||
                (geom.coord.are_near(_c[0], -_c[3], eps) && geom.coord.are_near(_c[1], _c[2], eps))) &&
               geom.coord.are_near(_c[0] * _c[0] + _c[1] * _c[1], 1.0, eps);
    }

    /** Check whether this transformation flips objects.
     * A transformation flips objects if it has a negative scaling component. */
    bool flips() const {
        // TODO shouldn't this be det() < 0?
        return cross(xAxis(), yAxis()) > 0;
    }


    /** Check whether this matrix is singular.
     * Singular matrices have no inverse, which means that applying them to a set of points
     * results in a loss of information.
     * @param eps Numerical tolerance
     * @return True iff the determinant is near zero. */
    bool isSingular(Coord eps = EPSILON) const
    { return geom.coord.are_near(det(), 0.0, eps); }


    /+ Compute other matrices +/

    Affine withoutTranslation() const
    {
        Affine ret = this;
        ret.setTranslation(Point(0,0));
        return ret;
    }

    /** Compute the inverse matrix.
     * Inverse is a matrix (denoted \f$A^{-1}) such that \f$AA^{-1} = A^{-1}A = I\f$.
     * Singular matrices have no inverse (for example a matrix that has two of its columns equal).
     * For such matrices, the identity matrix will be returned instead.
     * @param eps Numerical tolerance
     * @return Inverse of the matrix, or the identity matrix if the inverse is undefined.
     * @post (m * m.inverse()).isIdentity() == true */
    Affine inverse() const
    {
        Affine d = identity();
        
        Coord mx = math.fmax(math.fabs(_c[0]) + math.fabs(_c[1]), 
                             math.fabs(_c[2]) + math.fabs(_c[3])); // a random matrix norm (either l1 or linfty)
        if (mx > 0) {
            const Coord determ = det();
            if (!rel_error_bound(determ, mx*mx)) {
                const Coord ideterm = 1.0 / (determ);
                
                d._c[0] =  _c[3] * ideterm;
                d._c[1] = -_c[1] * ideterm;
                d._c[2] = -_c[2] * ideterm;
                d._c[3] =  _c[0] * ideterm;
                d._c[4] = (-_c[4] * d._c[0] - _c[5] * d._c[2]);
                d._c[5] = (-_c[4] * d._c[1] - _c[5] * d._c[3]);
            } else {
                d.setIdentity();
            }
        } else {
            d.setIdentity();
        }

        return d;
    }

    /+ Compute scalar values +/

    /** Calculate the determinant.
     * @return \f$\det A\f$. */
    Coord det() const // TODO this can overflow
    { return _c[0] * _c[3] - _c[1] * _c[2]; }

    /** Calculate the square of the descriminant.
     * This is simply the absolute value of the determinant.
     * @return \f$|\det A|\f$. */
    Coord descrim2() const
    { return math.fabs(det()); }

    /** Calculate the descriminant.
     * If the matrix doesn't contain a shearing or non-uniform scaling component, this value says
     * how will the length of any line segment change after applying this transformation
     * to arbitrary objects on a plane. The new length will be
     * @code line_seg.length() * m.descrim()) @endcode
     * @return \f$\sqrt{|\det A|}\f$. */
    Coord descrim() const
    { return math.sqrt(descrim2()); }

    private Coord[6] _c;
}

unittest
{
    /* Equality */
    Affine e = identity(); // identity
    Affine a = [1, 2, 3, 4, 5, 6];
    assert(e == e);
    assert(e == e.identity());
    assert(e != a);
    a = e;
    assert(e == a);
    
    /* Classification */
    assert(a.isIdentity());
    assert(a.isTranslation());
    assert(a.isScale());
    assert(a.isUniformScale());
    assert(a.isRotation());
    assert(a.isHShear());
    assert(a.isVShear());
    assert(a.isZoom());
    assert(!a.isNonzeroTranslation());
    assert(!a.isNonzeroScale());
    assert(!a.isNonzeroUniformScale());
    assert(!a.isNonzeroRotation());
    assert(!a.isNonzeroNonpureRotation());
    assert(!a.isNonzeroHShear());
    assert(!a.isNonzeroVShear());
    assert(a.preservesArea());
    assert(a.preservesAngles());
    assert(a.preservesDistances());
    assert(!a.flips());
    assert(!a.isSingular());
    
    a.setTranslation(Point(10, 15)); // pure translation
    assert(!a.isIdentity());
    assert(a.isTranslation());
    assert(!a.isScale());
    assert(!a.isUniformScale());
    assert(!a.isRotation());
    assert(!a.isHShear());
    assert(!a.isVShear());
    assert(a.isZoom());
    assert(a.isNonzeroTranslation());
    assert(!a.isNonzeroScale());
    assert(!a.isNonzeroUniformScale());
    assert(!a.isNonzeroRotation());
    assert(!a.isNonzeroNonpureRotation());
    assert(!a.isNonzeroHShear());
    assert(!a.isNonzeroVShear());
    assert(a.preservesArea());
    assert(a.preservesAngles());
    assert(a.preservesDistances());
    assert(!a.flips());
    assert(!a.isSingular());
}

void main() { }

/** Nearness predicate for affine transforms
 * Returns true if all entries of matrices are within eps of each other */
bool are_near(in Affine a, in Affine b, Coord eps = EPSILON)
{
    return geom.coord.are_near(a[0], b[0], eps) && geom.coord.are_near(a[1], b[1], eps) &&
           geom.coord.are_near(a[2], b[2], eps) && geom.coord.are_near(a[3], b[3], eps) &&
           geom.coord.are_near(a[4], b[4], eps) && geom.coord.are_near(a[5], b[5], eps);
}

/** Create an identity matrix.
 * This is a convenience function identical to Affine.identity(). */
deprecated("Use Affine.identity() instead") Affine identity()
{ return Affine.identity(); }

/+ Affine factories +/

/** Creates a Affine given an axis and origin point.
 *  The axis is represented as two vectors, which represent skew, rotation, and scaling in two dimensions.
 *  from_basis(Point(1, 0), Point(0, 1), Point(0, 0)) would return the identity matrix.

 @param x_basis the vector for the x-axis.
 @param y_basis the vector for the y-axis.
 @param offset the translation applied by the matrix.
 @return The new Affine.
 */
Affine from_basis(in Point x_basis, in Point y_basis, in Point offset = Point(0,0))
{
    return Affine(x_basis[X], x_basis[Y],
                  y_basis[X], y_basis[Y],
                  offset [X], offset [Y]);
}

// TODO: What's this!?!
/** Given a matrix m such that unit_circle = m*x, this returns the
 * quadratic form x*A*x = 1.
 * @relates Affine */
Affine elliptic_quadratic_form(in Affine m)
{
    Coord od = m[0] * m[1]  +  m[2] * m[3];
    Affine ret = [m[0]*m[0] + m[1]*m[1], od, od, m[2]*m[2] + m[3]*m[3], 0, 0];
    return ret;
}

/** Given a matrix (ignoring the translation) this returns the eigen
 * values and vectors. */
struct Eigen
{
    Point vectors[2];
    double values[2];

    @disable this();
    this(in Affine m)
    {
        const Coord B = -m[0] - m[3];
        const Coord C = m[0]*m[3] - m[1]*m[2];
        const Coord center = -B/2.0;
        const Coord delta = math.sqrt(B*B-4*C)/2.0;
        values[0] = center + delta;
        values[1] = center - delta;
        vectors[0] = unit_vector(rot90(Point(m[0]-values[0], m[1])));
        vectors[1] = unit_vector(rot90(Point(m[0]-values[1], m[1])));
    }

    this(in double[2][2] m)
    {
        vectors = [Point(0,0), Point(0,0)];

        const Coord B = -m[0][0] - m[1][1];
        const Coord C = m[0][0]*m[1][1] - m[1][0]*m[0][1];
        //double const desc = B*B-4*C;
        //double t = -0.5*(B+sgn(B)*desc);
        int n;
        values[0] = values[1] = 0;
        quadratic_roots(C, B, 1, n, values[0], values[1]);
        for (int i = 0; i < n; i++)
            vectors[i] = unit_vector(rot90(Point(m[0][0]-values[i], m[0][1])));
        for (int i = n; i < 2; i++) 
            vectors[i] = Point(0,0);
    }

    private void quadratic_roots(in double q0, in double q1, in double q2, ref int n, ref double r0, ref double r1)
    {
        if(q2 == 0) {
            if(q1 == 0) { // zero or infinite roots
                n = 0;
            } else {
                n = 1;
                r0 = -q0/q1;
            }
        } else {
            double desc = q1*q1 - 4*q2*q0;
            if (desc < 0)
                n = 0;
            else if (desc == 0) {
                n = 1;
                r0 = -q1/(2*q2);
            } else {
                n = 2;
                desc = math.sqrt(desc);
                double t = -0.5*(q1+math.sgn(q1)*desc);
                r0 = t/q2;
                r1 = q0/t;
            }
        }
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
