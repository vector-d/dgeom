/*
 * Piecewise function class
 *
 * Copyright 2007-2015 Michael Sloan <mgsloan@gmail.com>
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

module geom.piecewise;

public import geom.coord;
import geom.d2;
import geom.interval;
import geom.sbasis;
import std.math;

/**
 * Piecewise function class.
 * The Piecewise class manages a sequence of elements of a type as segments and
 * the 'cuts' between them. These cuts are time values which separate the pieces.
 * This function representation allows for more interesting functions, as it provides
 * a viable output for operations such as inversion, which may require multiple
 * SBasis to properly invert the original.
 * As for technical details, while the actual SBasis segments begin on the ﬁrst
 * cut and end on the last, the function is deﬁned throughout all inputs by ex-
 * tending the ﬁrst and last segments. The exact switching between segments is
 * arbitrarily such that beginnings (t=0) have preference over endings (t=1). This
 * only matters if it is discontinuous at the location.
 * \f[
 *      f(t) \rightarrow \left\{ 
 *      \begin{array}{cc}
 *      s_1,& t <= c_2 \\
 *      s_2,& c_2 <= t <= c_3\\
 *      \ldots \\
 *      s_n,& c_n <= t
 *      \end{array}\right.
 * \f]
 */
struct Piecewise(T)
{
    Coord[] cuts;
    T[] segs;
    ref inout(T) opIndex(size_t i) inout { return segs[i]; }
    //segs[i] stretches from cuts[i] to cuts[i+1].

    this(in T s)
    {
        push_cut(0.);
        push_seg(s);
        push_cut(1.);
    }
    
    this(in Piecewise!T o)
    { opAssign(o); }
    
    void opAssign(in Piecewise!T o)
    {
        cuts = (cast(Coord[])(o.cuts)).dup;
        segs = (cast(T[])(o.segs)).dup;
    }
    
    alias output_type = T.output_type;
    size_t input_dim() { return 1; }

    this(in output_type v)
    {
        push_cut(0.);
        push_seg(T(v));
        push_cut(1.);
    }
    
    void reserve(size_t i)
    {
        import object : reserve;
        segs.reserve(i);
        cuts.reserve(i + 1);
    }

    // Convenience/implementation hiding function to add cuts.
    void push_cut(Coord c)
    {
        assert(cuts.length == 0 || c > cuts[$]);
        cuts ~= c;
    }

    // Convenience/implementation hiding function to add segments.
    void push_seg(in T s) { segs ~= T(s); }

    /** Convenience/implementation hiding function to add segment/cut pairs.
     * Asserts that basic size and order invariants are correct
     */
    void push(in T s, Coord to)
    {
        assert(cuts.length - segs.length == 1);
        push_seg(s);
        push_cut(to);
    }

    size_t size() const { return segs.length; }
    bool empty() const { return segs.length == 0; }
    void clear() { segs = []; cuts = []; }

    /**
    *  The size of the returned vector equals n_derivs+1.
    */
    output_type[] valueAndDerivatives(Coord t, uint n_derivs) const
    {
        uint n = segN(t);
        output_type[] ret, val = segs[n].valueAndDerivatives(segT(t, n), n_derivs);
        Coord mult = 1;
        for (size_t i = 0; i < val.length; i++) {
            ret ~= val[i] * mult;
            mult /= cuts[n+1] - cuts[n];
        }
        return ret;
    }

    /** Returns the segment index which corresponds to a 'global' piecewise time.
     * Also takes optional low/high parameters to expedite the search for the segment.
     */
    uint segN(Coord t, int low = 0, int high = -1) const
    {
        high = (high == -1) ? cast(int)size() : high;
        if (t < cuts[0]) return 0;
        if (t >= cuts[size()]) return cast(uint)size() - 1;
        while (low < high) {
            int mid = (high + low) / 2; // Let's not plan on having huge (> INT_MAX / 2) cut sequences
            Coord mv = cuts[mid];
            if (mv < t) {
                if (t < cuts[mid + 1]) return mid; else low = mid + 1;
            } else if(t < mv) {
                if (cuts[mid - 1] < t) return mid - 1; else high = mid - 1;
            } else {
                return mid;
            }
        }
        return low;
    }

    /** Returns the time within a segment, given the 'global' piecewise time.
     * Also takes an optional index parameter which may be used for efficiency or to find the time on a
     * segment outside its range.  If it is left to its default, -1, it will call segN to find the index.
     */
    Coord segT(Coord t, int i = -1) const
    {
        if (i == -1) i = segN(t);
        assert(i >= 0);
        return (t - cuts[i]) / (cuts[i+1] - cuts[i]);
    }

    output_type opCall(Coord t) const { return valueAt(t); }
    output_type valueAt(Coord t) const
    {
        uint n = segN(t);
        return segs[n](segT(t, n));
    }

    Coord mapToDomain(Coord t, uint i) const
    { return (1-t)*cuts[i] + t*cuts[i+1]; } // same as: t * (cuts[i+1] - cuts[i]) + cuts[i]

    // Offsets the piecewise domain
    void offsetDomain(Coord o)
    {
        assert(isFinite(o));
        if (o != 0)
            for (size_t i = 0; i <= size(); i++)
                cuts[i] += o;
    }

    // Scales the domain of the function by a value. 0 will result in an empty Piecewise.
    void scaleDomain(Coord s)
    {
        assert(s > 0);
        if (s == 0) {
            cuts = []; segs = [];
            return;
        }
        for (size_t i = 0; i <= size(); i++)
            cuts[i] *= s;
    }

    // Retrieves the domain in interval form
    Interval domain() const { return Interval(cuts[0], cuts[$]); }

    // Transforms the domain into another interval
    void setDomain(Interval dom)
    {
        if (empty()) return;
        /* dom can not be empty
        if(dom.isEmpty()) {
            cuts.clear(); segs.clear();
            return;
        }*/
        Coord cf = cuts[0];
        Coord o = dom.min() - cf, s = dom.extent() / (cuts[$] - cf);
        for (size_t i = 0; i <= size(); i++)
            cuts[i] = (cuts[i] - cf) * s + o;
        //fix floating point precision errors.
        cuts[0] = dom.min();
        cuts[size()] = dom.max();
    }

    // Concatenates this Piecewise function with another, offseting time of the other to match the end.
    void concat(Piecewise!T other)
    {
        import object : reserve;
        if (other.empty()) return;

        if (empty()) {
            opAssign(other);
            return;
        }

        segs ~= (cast(T[])(other.segs)).dup;
        //segs.insert(segs.end(), other.segs.begin(), other.segs.end());
        Coord t = cuts[$] - other.cuts[0];
        cuts.reserve(cuts.length + other.size());
        for (size_t i = 0; i < other.size(); i++)
            push_cut(other.cuts[i + 1] + t);
    }

    // Like concat, but ensures continuity.
    void continuousConcat(Piecewise!T other)
    {
        if (other.empty()) return;

        if (empty()) { opAssign(other); return; }

        T.output_type y = segs[$].at1() - other.segs[0].at0();
        Coord t = cuts[$] - other.cuts[0];
        reserve(size() + other.size());
        for (size_t i = 0; i < other.size(); i++)
            push(other[i] + y, other.cuts[i + 1] + t);
    }

    // returns true if the Piecewise!T meets some basic invariants.
    bool invariants() const
    {
        // segs between cuts
        if (!(segs.length + 1 == cuts.length || (segs.length == 0 && cuts.length == 0)))
            return false;

        // cuts in order
        for (size_t i = 0; i < segs.length; i++)
            if(cuts[i] >= cuts[i+1])
                return false;
        return true;
    }


    /** Further subdivides the Piecewise<T> such that there is a cut at every value in c.
     * Precondition: c sorted lower to higher.
     *
     * Given Piecewise!T a and b:
     * Piecewise!T ac = a.partition(b.cuts);
     * Piecewise!T bc = b.partition(a.cuts);
     * ac.cuts should be equivalent to bc.cuts
     *
     */
    Piecewise!T partition(in Coord[] c) const
    {
        Piecewise!T pw = Piecewise!T(this);

        assert(pw.invariants());
        if (c.length == 0) return pw;

        Piecewise!T ret;
        ret.reserve(c.length + pw.cuts.length - 1);

        if (pw.empty()) {
            ret.cuts = c.dup;
            for (size_t i = 0; i < c.length - 1; i++)
                ret.push_seg(T.init);
            return ret;
        }

        size_t si = 0, ci = 0; // Segment index, Cut index

        // if the cuts have something earlier than the Piecewise!T, add portions of the first segment
        while (c[ci] < pw.cuts[0] && ci < c.length) {
            bool isLast = (ci == c.length-1 || c[ci + 1] >= pw.cuts[0]);
            ret.push_cut(c[ci]);
            ret.push_seg( pw.elem_portion(0, c[ci], isLast ? pw.cuts[0] : c[ci + 1]) );
            ci++;
        }

        ret.push_cut(pw.cuts[0]);
        Coord prev = pw.cuts[0];  // previous cut
        // Loop which handles cuts within the Piecewise!T domain
        // Should have the cuts = segs + 1 invariant
        while(si < pw.size() && ci <= c.length) {
            if (ci == c.length && prev <= pw.cuts[si]) { // cuts exhausted, straight copy the rest
                ret.segs ~= pw.segs[si .. $];
                ret.cuts ~= pw.cuts[si + 1 .. $];
                //ret.segs.insert(ret.segs.end(), pw.segs.begin() + si, pw.segs.end());
                //ret.cuts.insert(ret.cuts.end(), pw.cuts.begin() + si + 1, pw.cuts.end());
                return ret;
            } else if(ci == c.length || c[ci] >= pw.cuts[si + 1]) { // no more cuts within this segment, finalize
                if (prev > pw.cuts[si]) {      // segment already has cuts, so portion is required
                    ret.push_seg(pw[si].portion(pw.segT(prev, cast(int)si), 1.0));
                } else {                       // plain copy is fine
                    ret.push_seg(pw[si]);
                }
                ret.push_cut(pw.cuts[si + 1]);
                prev = pw.cuts[si + 1];
                si++;
            } else if(c[ci] == pw.cuts[si]){                  // coincident
                // Already finalized the seg with the code immediately above
                ci++;
            } else {                                         // plain old subdivision
                ret.push(pw.elem_portion(si, prev, c[ci]), c[ci]);
                prev = c[ci];
                ci++;
            }
        }

        // input cuts extend further than this Piecewise!T, extend the last segment.
        while (ci < c.length) {
            if(c[ci] > prev) {
                ret.push(pw.elem_portion(pw.size() - 1, prev, c[ci]), c[ci]);
                prev = c[ci];
            }
            ci++;
        }
        return ret;
    }

    /**
     *  Returns a portion of a piece of a Piecewise!T, given the piece's index and a to/from time.
     */
    T elem_portion(size_t i, Coord from, Coord to) const
    {
        assert(i < size());
        Coord rwidth = 1 / (cuts[i+1] - cuts[i]);
        return opIndex(i).portion((from - cuts[i]) * rwidth, (to - cuts[i]) * rwidth);
    }

    /**
     *  Returns a Piecewise!T with a defined domain of [min(from, to), max(from, to)].
     */
    Piecewise!T portion(Coord from, Coord to) const
    {
        import std.math;
        Piecewise!T pw = Piecewise!T(this);

        if (pw.empty() || from == to) return Piecewise!(T).init;

        Piecewise!T ret;

        Coord temp = from;
        from = fmin(from, to);
        to = fmax(temp, to);

        uint i = pw.segN(from);
        ret.push_cut(from);
        if (i == pw.size() - 1 || to <= pw.cuts[i + 1]) { // to/from inhabit the same segment
            ret.push(pw.elem_portion(i, from, to), to);
            return ret;
        }
        ret.push_seg(pw[i].portion(pw.segT(from, i), 1.0 ));
        i++;
        uint fi = pw.segN(to, i);
        ret.reserve(fi - i + 1);
        if (to == pw.cuts[fi]) fi-=1;

        ret.segs ~= pw.segs[i .. fi];
        ret.cuts ~= pw.cuts[i .. fi + 1];
        //ret.segs.insert(ret.segs.end(), pw.segs.begin() + i, pw.segs.begin() + fi);  //copy segs
        //ret.cuts.insert(ret.cuts.end(), pw.cuts.begin() + i, pw.cuts.begin() + fi + 1);  //and their cuts

        ret.push_seg( pw[fi].portion(0.0, pw.segT(to, fi)));
        if (to != ret.cuts[$]) ret.push_cut(to);
        assert(ret.invariants());
        return ret;
    }
}


Piecewise!T remove_short_cuts(T)(Piecewise!T f, Coord tol)
{
    if (f.empty()) return f;
    Piecewise!T ret;
    ret.reserve(f.size());
    ret.push_cut(f.cuts[0]);
    for (size_t i=0; i<f.size(); i++) {
        if (f.cuts[i+1]-f.cuts[i] >= tol || i==f.size()-1) {
            ret.push(f[i], f.cuts[i+1]);
        }
    }
    return ret;
}

alias PWD2 = Piecewise!(D2!SBasis);
alias D2PW = D2!(Piecewise!SBasis);

unittest
{
    PWD2 test;
    test = test.partition([0, 1]);
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
