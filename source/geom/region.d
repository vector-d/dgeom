/*
 * Uncrossed path for boolean algorithms
 *
 * Authors:
 *   ? <?@?.?>
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

module geom.region;

import geom.affine;
import geom.path;
import geom.path_intersection;
import geom.path_sequence;
import geom.point;
import geom.rect;

class Region
{
    this() { boundary = new Path(); }
    this(in Path p) { boundary = new Path(p); fill = path_direction(p); }
    this(in Path p, bool dir) { boundary = new Path(p); fill = dir; }
    this(in Path p, in Rect b) { boundary = new Path(p); box = b; fill = path_direction(p); }
    this(in Path p, in Rect b, bool dir) { boundary = new Path(p); box = b; fill = dir; }

    size_t size() const { return boundary.size(); }
    bool isFill() const { return fill; }
    Path toPath() const { return new Path(boundary); }
    alias toPath this;

    Rect boundsFast() const
    {
        //if (box.isEmpty) box = boundary.boundsFast();  // TODO this doesn't look right at all...
        return box;
    }

    bool contains(in Point p) const
    {
        if (!box.isEmpty() && !box.contains(p)) return false;
        return geom.path_intersection.contains(boundary, p);
    }

    bool contains(in Region other) const { return contains(other.boundary.initialPoint()); }

    Region inverse() const { return new Region(boundary.reversed(), box, !fill); }

    void invariants()
    {
        assert(self_crossings(boundary).length == 0);
    }

    Region opBinary(string op, T : Affine)(in T m) const if (op == "*")
    {
        Region r = new Region((m.flips() ? boundary.reversed() : boundary) * m, fill);
        if (!box.isEmpty() && m.isZoom()) r.box = box * m;
        return r;
    }

private:
    Path boundary = null;
    Rect box = Rect.empty();
    bool fill = true;
}

alias Regions = Region[];

size_t outer_index(in Regions ps)
{
    if (ps.length <= 1 || ps[0].contains(ps[1])) {
        return 0;
    } else {
        /* Since we've already shown that chunks[0] is not outside
           it can be used as an exemplar inner. */
        Point exemplar = ps[0].boundary.initialPoint();
        foreach (i; 1 .. ps.length) {
            if (ps[i].contains(exemplar)) {
                return i;
            }
        }
    }
    return ps.length;
}

//assumes they're already sanitized somewhat
Regions regions_from_paths(in PathSequence ps)
{
    Regions res;
    foreach (i; 0 .. ps.size())
        res ~= new Region(ps[i]);
    return res;
}

PathSequence paths_from_regions(in Regions rs)
{
    PathSequence res = new PathSequence;
    foreach (i; rs) res ~= i;
    return res;
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
