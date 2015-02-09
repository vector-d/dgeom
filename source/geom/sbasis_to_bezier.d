/*
 * Symmetric Power Basis - Bernstein Basis conversion routines
 *
 * Authors:
 *      Marco Cecchetti <mrcekets at gmail.com>
 *      Nathan Hurst <njh@mail.csse.monash.edu.au>
 *      Liam P. White
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

module geom.sbasis_to_bezier;

import geom.bezier;
import geom.choose;
import geom.coord;
import geom.d2;
import geom.point;
import geom.sbasis;

/*
 *  Symmetric Power Basis - Bernstein Basis conversion routines
 *
 *  some remark about precision:
 *  interval [0,1], subdivisions: 10^3
 *  - bezier_to_sbasis : up to degree ~72 precision is at least 10^-5
 *                       up to degree ~87 precision is at least 10^-3
 *  - sbasis_to_bezier : up to order ~63 precision is at least 10^-15
 *                       precision is at least 10^-14 even beyond order 200
 *
 *  interval [-1,1], subdivisions: 10^3
 *  - bezier_to_sbasis : up to degree ~21 precision is at least 10^-5
 *                       up to degree ~24 precision is at least 10^-3
 *  - sbasis_to_bezier : up to order ~11 precision is at least 10^-5
 *                       up to order ~13 precision is at least 10^-3
 *
 *  interval [-10,10], subdivisions: 10^3
 *  - bezier_to_sbasis : up to degree ~7 precision is at least 10^-5
 *                       up to degree ~8 precision is at least 10^-3
 *  - sbasis_to_bezier : up to order ~3 precision is at least 10^-5
 *                       up to order ~4 precision is at least 10^-3
 *
 *  references:
 *  this implementation is based on the following article:
 *  J.Sanchez-Reyes - The Symmetric Analogue of the Polynomial Power Basis
 */

alias binomial = choose!Coord;

/** Changes the basis of p to be bernstein.
 * @param p the Symmetric basis polynomial
 * @return the Bernstein basis polynomial
 *
 * if the degree is even q is the order in the symmetrical power basis,
 * if the degree is odd q is the order + 1
 * n is always the polynomial degree, i. e. the Bezier order
 * sz is the number of bezier handles.
 */
void sbasis_to_bezier(ref Bezier bz, in SBasis sb, size_t sz = 0)
{
    if (sb.size() == 0) {
        throw new Exception("size of sb is too small");
    }

    size_t q, n;
    bool even;
    if (sz == 0) {
        q = sb.size();
        if (sb[q-1][0] == sb[q-1][1]) {
            even = true;
            --q;
            n = 2*q;
        } else {
            even = false;
            n = 2*q-1;
        }
    } else {
        q = (sz > 2*sb.size()-1) ?  sb.size() : (sz+1)/2;
        n = sz-1;
        even = false;
    }
    bz.clear();
    bz.resize(n+1);
    double Tjk;
    for (size_t k = 0; k < q; ++k) {
        for (size_t j = k; j < n-k; ++j) { // j <= n-k-1
            Tjk = binomial(n-2*k-1, j-k);
            bz[j] += (Tjk * sb[k][0]);
            bz[n-j] += (Tjk * sb[k][1]); // n-k <-> [k][1]
        }
    } if (even) {
        bz[q] += sb[q][0];
    }
    // the resulting coefficients are with respect to the scaled Bernstein
    // basis so we need to divide them by (n, j) binomial coefficient
    for (size_t j = 1; j < n; ++j) {
        bz[j] /= binomial(n, j);
    }
    bz[0] = sb[0][0];
    bz[n] = sb[0][1];
}

/** Changes the basis of p to be Bernstein.
 * @param p the D2 Symmetric basis polynomial
 * @returns the D2 Bernstein basis polynomial
 *
 * sz is always the polynomial degree, i. e. the Bezier order
 */
void sbasis_to_bezier(ref Point[] bz, in D2!SBasis sb, size_t sz = 0)
{
    D2!Bezier bez;
    sbasis_to_bezier(bez, sb, sz);
    bz = bezier_points(bez);
}

void sbasis_to_bezier(ref D2!Bezier bz, in D2!SBasis sb, size_t sz)
{
    if (sz == 0) {
        sz = max(sb[X].size(), sb[Y].size())*2;
    }
    sbasis_to_bezier(bz[X], sb[X], sz);
    sbasis_to_bezier(bz[Y], sb[Y], sz);
}

/** Changes the basis of p to be sbasis.
 * @param p the Bernstein basis polynomial
 * @return the Symmetric basis polynomial
 *
 * if the degree is even q is the order in the symmetrical power basis,
 * if the degree is odd q is the order + 1
 * n is always the polynomial degree, i. e. the Bezier order
 */
void bezier_to_sbasis(ref SBasis sb, in Bezier bz)
{
    size_t n = bz.order();
    size_t q = (n+1) / 2;
    size_t even = (n & 1u) ? 0 : 1;
    sb.clear();
    sb.resize(q + even);
    double Tjk;
    for (size_t k = 0; k < q; ++k)
    {
        for (size_t j = k; j < q; ++j)
        {
            Tjk = sgn(j, k) * binomial(n-j-k, j-k) * binomial(n, k);
            sb[j][0] += (Tjk * bz[k]);
            sb[j][1] += (Tjk * bz[n-k]); // n-j <-> [j][1]
        }
        for (size_t j = k+1; j < q; ++j)
        {
            Tjk = sgn(j, k) * binomial(n-j-k-1, j-k-1) * binomial(n, k);
            sb[j][0] += (Tjk * bz[n-k]);
            sb[j][1] += (Tjk * bz[k]);   // n-j <-> [j][1]
        }
    }
    if (even)
    {
        for (size_t k = 0; k < q; ++k)
        {
            Tjk = sgn(q,k) * binomial(n, k);
            sb[q][0] += (Tjk * (bz[k] + bz[n-k]));
        }
        sb[q][0] += (binomial(n, q) * bz[q]);
        sb[q][1] = sb[q][0];
    }
    sb[0][0] = bz[0];
    sb[0][1] = bz[n];
}


/** Changes the basis of d2 p to be sbasis.
 * @param p the d2 Bernstein basis polynomial
 * @return the d2 Symmetric basis polynomial
 *
 * if the degree is even q is the order in the symmetrical power basis,
 * if the degree is odd q is the order + 1
 * n is always the polynomial degree, i. e. the Bezier order
 */
void bezier_to_sbasis(ref D2!SBasis sb, in Point[] bz)
{
    size_t n = bz.length - 1;
    size_t q = (n+1) / 2;
    size_t even = (n & 1u) ? 0 : 1;
    sb[X].clear();
    sb[Y].clear();
    sb[X].resize(q + even);
    sb[Y].resize(q + even);
    double Tjk;
    for (size_t k = 0; k < q; ++k)
    {
        for (size_t j = k; j < q; ++j)
        {
            Tjk = sgn(j, k) * binomial(n-j-k, j-k) * binomial(n, k);
            sb[X][j][0] += (Tjk * bz[k][X]);
            sb[X][j][1] += (Tjk * bz[n-k][X]);
            sb[Y][j][0] += (Tjk * bz[k][Y]);
            sb[Y][j][1] += (Tjk * bz[n-k][Y]);
        }
        for (size_t j = k+1; j < q; ++j)
        {
            Tjk = sgn(j, k) * binomial(n-j-k-1, j-k-1) * binomial(n, k);
            sb[X][j][0] += (Tjk * bz[n-k][X]);
            sb[X][j][1] += (Tjk * bz[k][X]);
            sb[Y][j][0] += (Tjk * bz[n-k][Y]);
            sb[Y][j][1] += (Tjk * bz[k][Y]);
        }
    }
    if (even)
    {
        for (size_t k = 0; k < q; ++k)
        {
            Tjk = sgn(q,k) * binomial(n, k);
            sb[X][q][0] += (Tjk * (bz[k][X] + bz[n-k][X]));
            sb[Y][q][0] += (Tjk * (bz[k][Y] + bz[n-k][Y]));
        }
        sb[X][q][0] += (binomial(n, q) * bz[q][X]);
        sb[X][q][1] = sb[X][q][0];
        sb[Y][q][0] += (binomial(n, q) * bz[q][Y]);
        sb[Y][q][1] = sb[Y][q][0];
    }
    sb[X][0][0] = bz[0][X];
    sb[X][0][1] = bz[n][X];
    sb[Y][0][0] = bz[0][Y];
    sb[Y][0][1] = bz[n][Y];
}

private int sgn(size_t j, size_t k)
{
    assert (j >= k);
    // we are sure that j >= k
    return ((j-k) &  1u) ? -1 : 1;
}

private T max(T)(T a, T b)
{ return a > b ? a : b; }

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
