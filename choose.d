/**
 *  Calculation of binomial cefficients
 *
 * Copyright 2006-2015 Nathan Hurst <njh@mail.csse.monash.edu.au>
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

module geom.choose;

// XXX: Can we keep only the left terms easily?
// this would more than halve the array
// row index becomes n2 = n/2, row2 = n2*(n2+1)/2, row = row2*2+(n&1)?n2:0
// we could also leave off the ones

T choose(T)(size_t n, size_t k)
{
    static T[] pascals_triangle;
    static size_t rows_done = 0;

    // indexing is (0,0,), (1,0), (1,1), (2, 0)...
    // to get (i, j) i*(i+1)/2 + j
    if (/*k < 0 ||*/ k > n) return 0;

    if (rows_done <= n) { // we haven't built that much of the triangle yet
        if (rows_done == 0) {
            pascals_triangle ~= 1;
            rows_done = 1;
        }

        while (rows_done <= n) {
            size_t p = pascals_triangle.length - rows_done;
            pascals_triangle ~= 1;
            for (size_t i = 0; i < rows_done - 1; ++i) {
                pascals_triangle ~= pascals_triangle[p] + pascals_triangle[p+1];
		p++;
            }

            pascals_triangle ~= 1;
            rows_done++;
        }
    }
    size_t row = (n*(n+1))/2;
    return pascals_triangle[row+k];
}

// Is it faster to store them or compute them on demand?
/+T choose(T)(size_t n, size_t k)
{
    T r = 1;
    for (size_t i = 1; i <= k; ++i)
        r = (r*(n-k+i))/i;
    return r;
}+/

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
