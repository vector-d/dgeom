/*
 * Path intersection
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

module geom.crossing;

import geom.coord;
import geom.rect;

class Crossing
{
    bool dir = false;       // true: along a, a becomes outside
    Coord ta = 0, tb = 1;   // time on a and b of crossing
    size_t a = 0, b = 1;    // storage of indices

    this() {}
    this(Coord t_a, Coord t_b, bool direction) { ta = t_a; tb = t_b; dir = direction; }
    this(Coord t_a, Coord t_b, size_t ai, size_t bi, bool direction) { ta = t_a; tb = t_b; a = ai; b = bi; dir = direction; }

    size_t getOther(size_t cur) const { return a == cur ? b : a; }
    Coord getTime(size_t cur) const { return a == cur ? ta : tb; }
    Coord getOtherTime(size_t cur) const { return a == cur ? tb : ta; }
    bool onIx(size_t ix) const { return a == ix || b == ix; }
}

auto crossing_sorter(alias rev, alias ix)(in Crossing a, in Crossing b)
{
    if (rev) 
        return (ix == a.a ? a.ta : a.tb) <
               (ix == b.a ? b.ta : b.tb);
    else
        return (ix == a.a ? a.ta : a.tb) >
               (ix == a.a ? a.ta : a.tb);
};

alias Crossings = Crossing[];
alias CrossingSet = Crossings[];

void sort_crossings(ref Crossings cr, size_t ix)
{
    import std.algorithm : sort;

    sort!(crossing_sorter!(false, ix))(cr);
}

class Crosser(T)
{
    Crossings crossings(in T a, in T b) const { return crossings([a], [b])[0]; }

    CrossingSet crossings(in T[] a, in T[] b) const
    {
        CrossingSet results = new CrossingSet(a.length + b.length);

        size_t[][] cull = sweep_bounds(bounds(a), bounds(b));
        foreach (ix, i; cull) {
            foreach (jx, j; i) {
                auto jc = j + a.length;
                Crossings cr = crossings(a[ix], b[j]);
                foreach (ref k; cr) { k.a = ix; k.b = jc; }

                // Sort and add A-sorted crossings
                sort_crossings(cr, ix);
                auto n = merge!(crossing_sorter!(false, ix))(results[ix], cr);
                results[ix] = n;
                
                // Sort and add B-sorted crossings
                sort_crossings(cr, jc);
                n = merge!(crossing_sorter!(false, jc))(results[jc], cr);
                results[jc] = n;
            }
        }
        return results;
    }
}

/**
 * Make a list of pairs of self intersections in a list of Rects.
 * 
 * @param rs Rect[].
 * @param d dimension to sweep along
 *
 * [(A = rs[i], B = rs[j]) for i,J in enumerate(pairs) for j in J]
 * then A.left <= B.left
 */
size_t[][] sweep_bounds(in Rect[] rs, size_t d = X)
{
    import std.algorithm : sort, find, remove;
    import object : reserve;

    Event[] events;
    events.reserve(rs.length*2);

    size_t[][] pairs;
    pairs.reserve(rs.length);

    foreach (i, e; rs) {
        events ~= Event(e[d][X], i, false);
        events ~= Event(e[d][Y], i, true);
    }

    events.sort();

    size_t[] open;
    foreach (i, ei; events) {
        size_t ix = ei.ix;
        if (events[i].closing) {
            auto iter = open.find(ix);
            //if(iter != open.end())
            open.remove(iter);
        } else {
            foreach (jx; open) {
                if (rs[jx][1-d].intersects(rs[ix][1-d])) {
                    pairs[jx] ~= ix;
                }
            }
            open ~= ix;
        }
    }

    return pairs;
}

/**
 * @brief Make a list of pairs of red-blue intersections between two lists of Rects.
 * 
 * @param a: Rect[].
 * @param b: Rect[].
 * @param d: dimension to scan along
 *
 * [(A = rs[i], B = rs[j]) for i,J in enumerate(pairs) for j in J]
 * then A.left <= B.left, A in a, B in b
 */
size_t[][] sweep_bounds(in Rect[] a, in Rect[] b, size_t d = X)
{
    import std.algorithm : find, remove, sort;
    import object : reserve;

    size_t[][] pairs = new size_t[][a.length];
    if (a.length == 0 || b.length == 0) return pairs;

    Event[][2] events;
    events[X].reserve(a.length*2);
    events[Y].reserve(b.length*2);

    foreach (n; 0 .. 2) {
        auto sz = n ? b.length : a.length;
        events[n].reserve(sz*2);
        foreach (i; 0 .. sz) {
            Rect r = n ? b[i] : a[i];
            events[n] ~= Event(r[d][0], i, false);
            events[n] ~= Event(r[d][1], i, true);
        }
        events[n].sort();
    }

    size_t[][2] open;
    bool n = events[Y][0] < events[X][0];

    for (size_t[2] i = [0,0]; i[n] < events[n].length; ) {
        size_t ix = events[n][i[n]].ix;
        bool closing = events[n][i[n]].closing;
        if (closing) {
            open[n].remove(open[n].find(ix));
        } else {
            if(n) {
                //n = 1
                //opening a B, add to all open a
                foreach (j; 0 .. open[X].length) {
                    size_t jx = open[X][j];
                    if (a[jx][1-d].intersects(b[ix][1-d])) {
                        pairs[jx] ~= ix;
                    }
                }
            } else {
                //n = 0
                //opening an A, add all open b
                foreach (j; 0 .. open[Y].length) {
                    size_t jx = open[Y][j];
                    if (b[jx][1-d].intersects(a[ix][1-d])) {
                        pairs[ix] ~= jx;
                    }
                }
            }
            open[n] ~= ix;
        }
        i[n]++;

	if (i[n] >= events[n].length) break;

        n = (events[!n][i[!n]] < events[n][i[n]]) ? !n : n;
    }

    return pairs;
}

void merge_crossings(ref Crossings a, ref Crossings b, size_t i)
{
    sort_crossings(b, i);
    a = merge!(crossing_sorter!(false, i))(a, b);
}

void offset_crossings(ref Crossings cr, Coord a, Coord b)
{
    foreach (i; cr) {
        i.ta += a;
        i.tb += b;
    }
}

Crossings reverse_ta(in Crossings cr, in Coord[] max)
{
    Crossings ret;
    foreach (i; cr) {
        Coord mx = max[i.a];
        ret ~= new Crossing(i.ta > mx+0.01 ? (1 - (i.ta - mx) + mx) : mx - i.ta, i.tb, !i.dir);
    }
    return ret;
}

Crossings reverse_tb(in Crossings cr, size_t split, in Coord[] max)
{
    Crossings ret;
    foreach (i; cr) {
        double mx = max[i.b - split];
        ret ~= new Crossing(i.ta, i.tb > mx+0.01 ? (1 - (i.tb - mx) + mx) : mx - i.tb, !i.dir);
    }
    return ret;
}

CrossingSet reverse_ta(in CrossingSet cr, size_t split, in Coord[] max)
{
    import std.algorithm : reverse;

    CrossingSet ret;
    foreach (i; 0 .. cr.length) {
        Crossings res = reverse_ta(cr[i], max);
        if (i < split) res.reverse();
        ret ~= res;
    }
    return ret;
}

CrossingSet reverse_tb(in CrossingSet cr, size_t split, in Coord[] max)
{
    import std.algorithm : reverse;

    CrossingSet ret;
    foreach (i; 0 .. cr.length) {
        Crossings res = reverse_tb(cr[i], split, max);
        if (i >= split) res.reverse();
        ret ~= res;
    }
    return ret;
}

/** Delete any duplicates in an array of crossings.
 * A crossing is considered to be a duplicate when it has both t_a and t_b near to another crossing's t_a and t_b
 * For example, duplicates will be found when calculating the intersections of a linesegment with a polygon, if the
 * endpoint of that line coincides with a cusp node of the polygon. In that case, an intersection will be found of
 * the linesegment with each of the polygon's linesegments extending from the cusp node (i.e. two intersections).
 */
void delete_duplicates(ref Crossings crs)
{
    import std.algorithm : remove;

    auto i = crs.length - 1;
    for (; i != 0; --i) {
        auto j = i;
        while (--j != 0) {
            if (are_near(crs[i].ta, crs[j].ta) && are_near(crs[i].tb, crs[j].tb)) {
                crs.remove(i);
                break; // out of while loop, and continue with next iteration of for loop
            }
        }
    }
}

Rect[] bounds(C)(in C a)
{
    Rect[] rs;
    foreach (i; a) {
        Rect bb = i.boundsFast();
        if (!bb.isEmpty()) {
            rs ~= bb;
        }
    }
    return rs;
}

import geom.path;
// provide specific method for Paths because paths can be closed or open. Path.size is named somewhat wrong...
Rect[] bounds(C : Path)(in C a)
{
    Rect[] rs;
    foreach (i; 0 .. a.size_default) {
        auto bb = a[i].boundsFast();
        if (!bb.isEmpty()) {
            rs ~= bb;
        }
    }
    return rs;
}

T[] merge(alias Comparison, T)(T[] a, T[] b)
{
    import std.array;

    int idx = 0;
    T[] result = uninitializedArray!(T[])(a.length + b.length);

    while (a.length || b.length)
    {
        if (b.length == 0) {
            result[idx++] = a[0];
            a.popFront;
        } else if (a.length == 0) {
            result[idx++] = b[0];
            b.popFront;
        } else if (Comparison(a[0], b[0])) {
            result[idx++] = a[0];
            a.popFront;
        } else {
            result[idx++] = b[0];
            b.popFront;
        }
    }

    return result;
}

private struct Event
{
    Coord x;
    size_t ix;
    bool closing;

    @disable this();
    this(Coord pos, size_t i, bool c) { x = pos; ix = i; closing = c; }

    // Lexicographic ordering by x then closing
    int opCmp(in Event other) const
    {
        if (x < other.x) return 1;
        if (x > other.x) return -1;
        return closing < other.closing;
    }

    bool opEquals(in Event other) const
    {
        return other.x == x && other.ix == ix && other.closing == closing;
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
