/*
 * Defines the Coord "real" type with sufficient precision for coordinates.
 *
 * Copyright 2006 Nathan Hurst <njh@mail.csse.monash.edu.au>
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

module geom.coord;

/**
 * Floating point type used to store coordinates.
 *
 * You may safely assume that double (or even float) provides enough precision for storing
 * on-canvas points, and hence that double provides enough precision for dot products of
 * differences of on-canvas points.
 */
alias Coord = double;
alias IntCoord = int;

const Coord EPSILON = 1e-5; //1e-18;

/* Dim2 */
const size_t X = 0;
const size_t Y = 1;

deprecated("Redundant") Coord infinity() { return Coord.infinity; }

/* IMPL: NearConcept */
bool are_near(in Coord a, in Coord b, in double eps = EPSILON) { return a-b <= eps && a-b >= -eps; }
bool rel_error_bound(in Coord a, in Coord b, in double eps = EPSILON) { return a <= eps*b && a >= -eps*b; }

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
