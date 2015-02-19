/** root finding for sbasis functions.
 * Copyright 2006 N Hurst
 * Copyright 2007 JF Barraud
 *
 * It is more efficient to find roots of f(t) = c_0, c_1, ... all at once, rather than iterating.
 *
 * Todo/think about:
 *  multi-roots using bernstein method, one approach would be:
    sort c
    take median and find roots of that
    whenever a segment lies entirely on one side of the median,
    find the median of the half and recurse.

    in essence we are implementing quicksort on a continuous function

 *  the gsl poly roots finder is faster than bernstein too, but we don't use it for 3 reasons:

 a) it requires convertion to poly, which is numerically unstable

 b) it requires gsl (which is currently not a dependency, and would bring in a whole slew of unrelated stuff)

 c) it finds all roots, even complex ones.  We don't want to accidently treat a nearly real root as a real root

From memory gsl poly roots was about 10 times faster than bernstein in the case where all the roots
are in [0,1] for polys of order 5.  I spent some time working out whether eigenvalue root finding
could be done directly in sbasis space, but the maths was too hard for me. -- njh

jfbarraud: eigenvalue root finding could be done directly in sbasis space ?

njh: I don't know, I think it should.  You would make a matrix whose characteristic polynomial was
correct, but do it by putting the sbasis terms in the right spots in the matrix.  normal eigenvalue
root finding makes a matrix that is a diagonal + a row along the top.  This matrix has the property
that its characteristic poly is just the poly whose coefficients are along the top row.

Now an sbasis function is a linear combination of the poly coeffs.  So it seems to me that you
should be able to put the sbasis coeffs directly into a matrix in the right spots so that the
characteristic poly is the sbasis.  You'll still have problems b) and c).

We might be able to lift an eigenvalue solver and include that directly into 2geom.  Eigenvalues
also allow you to find intersections of multiple curves but require solving n*m x n*m matrices.

 **/

module geom.sbasis_roots;

import geom.interval;
import geom.linear;
import geom.sbasis;
import std.math;


//-- multi_roots ------------------------------------
// goal: solve f(t)=c for several c at once.
/* algo: -compute f at both ends of the given segment [a,b].
         -compute bounds m<df(t)<M for df on the segment.
         let c and C be the levels below and above f(a):
         going from f(a) down to c with slope m takes at least time (f(a)-c)/m
         going from f(a)  up  to C with slope M takes at least time (C-f(a))/M
         From this we conclude there are no roots before a'=a+min((f(a)-c)/m,(C-f(a))/M).
         Do the same for b: compute some b' such that there are no roots in (b',b].
         -if [a',b'] is not empty, repeat the process with [a',(a'+b')/2] and [(a'+b')/2,b'].
  unfortunately, extra care is needed about rounding errors, and also to avoid the repetition of roots,
  making things tricky and unpleasant...
*/
//TODO: Make sure the code is "rounding-errors proof" and take care about repetition of roots!

private {

int upper_level(in Coord[] levels, Coord x, Coord tol = 0.)
{
    import std.range;
    auto r = assumeSorted(levels);
    auto p = r.upperBound(x - tol);
    return cast(int)(levels.length - p.length);
}

void multi_roots_internal(in SBasis f, in SBasis df, in Coord[] levels, ref Coord[][] roots, Coord htol, Coord vtol, Coord a, Coord fa, Coord b, Coord fb)
{
    if (f.size() == 0) {
        int idx = upper_level(levels, 0, vtol);

        if (idx < cast(int)levels.length && fabs(levels[idx]) <= vtol) {
            roots[idx] ~= [a, b];
        }
        return;
    }

    if ((b - a) < htol) {
        //TODO: use different tol for t and f ?
        //TODO: unsigned idx ? (remove int casts when fixed)

        int idx = cast(int)fmin(upper_level(levels, fa, vtol), upper_level(levels, fb, vtol));
        if (idx == cast(int)levels.length)
            idx--;

        Coord c = levels[idx];

        if ((fa-c)*(fb-c) <= 0 || fabs(fa-c) < vtol || fabs(fb-c) < vtol) {
            roots[idx] ~= (a+b)/2;
        }
        return;
    }

    int idxa = upper_level(levels, fa, vtol);
    int idxb = upper_level(levels, fb, vtol);

    Interval bs = bounds_local(df, Interval(a, b));

    // first times when a level (higher or lower) can be reached from a or b.
    Coord ta_hi, tb_hi, ta_lo, tb_lo;
    ta_hi = ta_lo = b+1; // default values => no root there.
    tb_hi = tb_lo = a-1; // default values => no root there.

    if (idxa < cast(int)levels.length && fabs(fa - levels[idxa]) < vtol) { // a can be considered a root.
        //ta_hi=ta_lo=a;
        roots[idxa] ~= a;
        ta_hi = ta_lo = a+htol;
    } else {
        if (bs.max()>0 && idxa < cast(int)levels.length)
            ta_hi=a+(levels[idxa  ]-fa)/bs.max();
        if (bs.min()<0 && idxa>0)
            ta_lo=a+(levels[idxa-1]-fa)/bs.min();
    }

    if (idxb < cast(int)levels.length && fabs(fb-levels[idxb]) < vtol) { // b can be considered a root.
        //tb_hi=tb_lo=b;
        roots[idxb] ~= b;
        tb_hi = tb_lo = b - htol;
    } else {
        if (bs.min()<0 && idxb< cast(int)levels.length)
            tb_hi=b+(levels[idxb  ]-fb)/bs.min();
        if (bs.max()>0 && idxb>0)
            tb_lo=b+(levels[idxb-1]-fb)/bs.max();
    }

    Coord t0 = fmin(ta_hi, ta_lo);
    Coord t1 = fmax(tb_hi, tb_lo);

    // hum, rounding errors frighten me! so I add this +tol...
    if (t0 > t1 + htol) return; // no root here.

    if (fabs(t1-t0) < htol) {
        multi_roots_internal(f,df,levels,roots,htol,vtol,t0,f(t0),t1,f(t1));
    } else {
        Coord t, t_left, t_right;
        t_left = t_right = t = (t0 + t1) / 2;

        Coord ft, ft_left, ft_right;
        ft_left = ft_right = ft = f(t);

        int idx = upper_level(levels, ft, vtol);

        if (idx < cast(int)levels.length && fabs(ft-levels[idx]) < vtol) { // t can be considered a root.
            roots[idx] ~= t;
            // we do not want to count it twice (from the left and from the right)
            t_left = t - htol/2;
            t_right = t + htol/2;
            ft_left = f(t_left);
            ft_right = f(t_right);
        }

        multi_roots_internal(f, df, levels, roots, htol, vtol, t0     , f(t0)   , t_left, ft_left);
        multi_roots_internal(f, df, levels, roots, htol, vtol, t_right, ft_right, t1    , f(t1)  );
    }
}

/** Solve f(t)=c for several c at once.
    \param f sbasis function
    \param levels vector of 'y' values
    \param htol, vtol 
    \param a, b left and right bounds
    \returns an array of arrays, one for each y giving roots

Effectively computes:
results = roots(f(y_i)) for all y_i

* algo: -compute f at both ends of the given segment [a,b].
         -compute bounds m<df(t)<M for df on the segment.
         let c and C be the levels below and above f(a):
         going from f(a) down to c with slope m takes at least time (f(a)-c)/m
         going from f(a)  up  to C with slope M takes at least time (C-f(a))/M
         From this we conclude there are no roots before a'=a+min((f(a)-c)/m,(C-f(a))/M).
         Do the same for b: compute some b' such that there are no roots in (b',b].
         -if [a',b'] is not empty, repeat the process with [a',(a'+b')/2] and [(a'+b')/2,b'].
  unfortunately, extra care is needed about rounding errors, and also to avoid the repetition of roots,
  making things tricky and unpleasant...

TODO: Make sure the code is "rounding-errors proof" and take care about repetition of roots!
*/
Coord[][] multi_roots(in SBasis f, in Coord[] levels, Coord htol, Coord vtol, Coord a, Coord b)
{
    Coord[][] roots;
    roots.length = levels.length;

    SBasis df = derivative(SBasis(f));
    multi_roots_internal(f, df, levels, roots, htol, vtol, a, f(a), b, f(b));

    return(roots);
}


bool compareIntervalMin(Interval I, Interval J)
{ return I.min()<J.min(); }

bool compareIntervalMax(Interval I, Interval J)
{ return I.max()<J.max(); }

//find the first interval whose max is >= x
uint upper_level(in Interval[] levels, double x)
{
    import std.range;
    auto r = assumeSorted!("a.max() < b.max()")(levels);
    auto p = r.lowerBound(Interval(x, x));
    return cast(uint)(levels.length - p.length);
}

Interval[] fuseContiguous(in Interval[] sets, double tol = 0.)
{
    Interval[] result;
    if (sets.length == 0) return result;

    result ~= sets[0];
    foreach (i; sets[1 .. $]) {
        if (result[$].max() + tol >= i.min())
            result[$].unionWith(i);
        else
            result ~= i;
    }
    return result;
}

/** level_sets internal method.
* algorithm: (~adaptation of Newton method versus 'accroissements finis')
   -compute f at both ends of the given segment [a,b].
   -compute bounds m<df(t)<M for df on the segment.
    Suppose f(a) is between two 'levels' c and C. Then
      f wont enter c before a + (f(a)-c.max())/m
      f wont enter C before a + (C.min()-f(a))/M
    From this we conclude nothing happens before a'=a+min((f(a)-c.max())/m,(C.min()-f(a))/M).
    We do the same for b: compute some b' such that nothing happens in (b',b].
    -if [a',b'] is not empty, repeat the process with [a',(a'+b')/2] and [(a'+b')/2,b'].

    If f(a) or f(b) belongs to some 'level' C, then use the same argument to find a' or b' such
    that f remains in C on [a,a'] or [b',b]. In case f is monotonic, we also know f won't enter another
    level before or after some time, allowing us to restrict the search a little more.

  unfortunately, extra care is needed about rounding errors, and also to avoid the repetition of roots,
  making things tricky and unpleasant...
*/

void level_sets_internal(in SBasis f, in SBasis df, in Interval[] levels, ref Interval[][] solsets, Coord a, Coord fa, Coord b, Coord fb, Coord tol = 1e-5)
{
    if (f.size() == 0) {
        uint idx = upper_level(levels, 0.);
        if (idx < levels.length && levels[idx].contains(0.)){
            solsets[idx] ~= Interval(a,b);
        }
        return;
    }

    uint idxa = upper_level(levels, fa);
    uint idxb = upper_level(levels, fb);

    Interval bs = bounds_local(df, Interval(a,b));

    // first times when a level (higher or lower) can be reached from a or b.
    Coord ta_hi; // f remains below next level for t<ta_hi
    Coord ta_lo; // f remains above prev level for t<ta_lo
    Coord tb_hi; // f remains below next level for t>tb_hi
    Coord tb_lo; // f remains above next level for t>tb_lo

    ta_hi = ta_lo = b+1; // default values => no root there.
    tb_hi = tb_lo = a-1; // default values => no root there.

    // --- if f(a) belongs to a level.-------
    if (idxa < levels.length && levels[idxa].contains(fa)) {
    	// find the first time when we may exit this level.
    	ta_lo = a + (levels[idxa].min() - fa)/bs.min();
    	ta_hi = a + (levels[idxa].max() - fa)/bs.max();
    	if (ta_lo < a || ta_lo > b) ta_lo = b;
    	if (ta_hi < a || ta_hi > b) ta_hi = b;
    	// move to that time for the next iteration.
    	solsets[idxa] ~= Interval(a, fmin(ta_lo, ta_hi));
    } else {
        // --- if f(b) does not belong to a level.-------
        if (idxa == 0) {
            ta_lo = b;
        } else {
            ta_lo = a + (levels[idxa - 1].max() - fa)/bs.min();
            if (ta_lo < a) ta_lo = b;
        }

        if (idxa == levels.length) {
            ta_hi = b;
        } else {
            ta_hi = a + (levels[idxa].min() - fa)/bs.max();
            if (ta_hi < a) ta_hi = b;
        }
    }

    // --- if f(b) belongs to a level.-------
    if (idxb < levels.length && levels[idxb].contains(fb)) {
    	// find the first time from b when we may exit this level.
    	tb_lo = b + (levels[idxb].min() - fb) / bs.max();
    	tb_hi = b + (levels[idxb].max() - fb) / bs.min();
    	if (tb_lo > b || tb_lo < a) tb_lo = a;
    	if (tb_hi > b || tb_hi < a) tb_hi = a;
    	// move to that time for the next iteration.
    	solsets[idxb] ~= Interval(fmax(tb_lo, tb_hi), b);
    } else {
        // --- if f(b) does not belong to a level.-------
        if (idxb == 0) {
            tb_lo = a;
        } else {
            tb_lo = b + (levels[idxb-1].max() - fb)/bs.max();
            if (tb_lo > b) tb_lo = a;
        }
        if (idxb == levels.length) {
            tb_hi = a;
        } else {
            tb_hi = b + (levels[idxb].min() - fb)/bs.min();
            if (tb_hi > b) tb_hi = a;
        }


    	if ( bs.min() < 0 && idxb < levels.length )
            tb_hi = b + ( levels[idxb  ].min() - fb ) / bs.min();
        if ( bs.max() > 0 && idxb > 0 )
            tb_lo = b + ( levels[idxb-1].max() - fb ) / bs.max();
    }

    // let [t0,t1] be the next interval where to search.
    Coord t0 = fmin(ta_hi, ta_lo);
    Coord t1 = fmax(tb_hi, tb_lo);

    if (t0 >= t1) return; // no root here.

    // if the interval is smaller than our resolution:
    // pretend f simultaneously meets all the levels between f(t0) and f(t1)...
    if (t1 - t0 <= tol) {
    	Interval f_t0t1 = Interval(f(t0), f(t1));
    	uint idxmin = cast(uint)fmin(idxa, idxb);
    	uint idxmax = cast(uint)fmax(idxa, idxb);
    	// push [t0,t1] into all crossed level. Cheat to avoid overlapping intervals on different levels?
    	if (idxmax > idxmin) {
            for (uint idx = idxmin; idx < idxmax; idx++) {
                solsets[idx] ~= Interval(t0, t1);
            }
    	}

    	if (idxmax < levels.length && f_t0t1.intersects(levels[idxmax])) {
            solsets[idxmax] ~= Interval(t0, t1);
    	}
    	return;
    }

    // To make sure we finally exit the level jump at least by tol:
    t0 = fmin(fmax(t0, a + tol), b);
    t1 = fmax(fmin(t1, b - tol), a);

    Coord t =(t0+t1)/2;
    Coord ft = f(t);
    level_sets_internal(f, df, levels, solsets, t0, f(t0), t, ft);
    level_sets_internal(f, df, levels, solsets, t, ft, t1, f(t1));
}

Interval[][] level_sets(in SBasis f, in Interval[] levels, Coord a, Coord b, Coord tol)
{
    import std.algorithm : sort;
    Interval[][] solsets;
    solsets.length = levels.length;

    SBasis df = derivative(SBasis(f));
    level_sets_internal(f, df, levels, solsets, a, f(a), b, f(b), tol);

    // Fuse overlapping intervals...
    foreach(i; solsets) {
        if (i.length == 0) continue;
        
    }
    for (size_t i = 0; i<solsets.length; i++){
    	if ( solsets[i].length == 0 ) continue;
    	sort!"a.min()<a.min()"(solsets[i]);
    	solsets[i] = fuseContiguous( solsets[i], tol );
    }
    return solsets;
}

Interval[] level_set(in SBasis f, Coord level, Coord vtol, Coord a, Coord b, Coord tol)
{
    Interval fat_level = Interval(level - vtol, level + vtol);
    return level_set(f, fat_level, a, b, tol);
}

Interval[] level_set(in SBasis f, in Interval level, Coord a, Coord b, Coord tol)
{
    Interval[] levels;
    levels ~= level;
    return level_sets(f,levels, a, b, tol)[0];
}

Interval[][] level_sets(in SBasis f, in Coord[] levels, Coord vtol, Coord a, Coord b, Coord tol)
{
    Interval[] fat_levels;
    fat_levels.length = levels.length;
    foreach(i, m; levels) {
        fat_levels[i] = Interval(m - vtol, m + vtol);
    }
    return level_sets(f, fat_levels, a, b, tol);
}


//-------------------------------------
//-------------------------------------


void subdiv_sbasis(in SBasis s, ref Coord[] roots, Coord left, Coord right)
{
    Interval bs = bounds_fast(s);

    if (bs.isEmpty() || bs.min() > 0 || bs.max() < 0)
        return; // no roots here

    if (s.tailError(1) < 1e-7) {
        Coord t = s[0][0] / (s[0][0] - s[0][1]);
        roots ~= left*(1-t) + t*right;
        return;
    }

    Coord middle = (left + right)/2;
    subdiv_sbasis(compose(s, SBasis(Linear(0, 0.5))), roots, left, middle);
    subdiv_sbasis(compose(s, SBasis(Linear(0.5, 1.))), roots, middle, right);
}

} // private

// It is faster to use the bernstein root finder for small degree polynomials (<100?.

Coord[] roots1(in SBasis s)
{
    Coord[] res;
    Coord d = s[0][0] - s[0][1];
    if (d != 0) {
        Coord r = s[0][0] / d;
        if (0 <= r && r <= 1)
            res ~= r;
    }
    return res;
}

Coord[] roots1(in SBasis s, in Interval ivl)
{
    Coord[] res;
    Coord d = s[0][0] - s[0][1];
    if (d != 0) {
        Coord r = s[0][0] / d;
        if (ivl.contains(r))
            res ~= r;
    }
    return res;
}

/** Find all t s.t s(t) = 0
 * @param a sbasis function
 * @return array of zeros (roots)
 */

Coord[] roots(in SBasis s)
{
    switch(s.size()) {
        case 0:
            return [];
        case 1:
            return roots1(s);
        default:
        {
            import geom.bezier, geom.sbasis_to_bezier;
            Bezier bz;
            sbasis_to_bezier(bz, s);
            return bz.roots();
        }
    }
}

Coord[] roots(in SBasis s, in Interval ivl)
{
    switch(s.size()) {
        case 0:
            return [];
        case 1:
            return roots1(s, ivl);
        default:
        {
            import geom.bezier, geom.sbasis_to_bezier;
            Bezier bz;
            sbasis_to_bezier(bz, s);
            return bz.roots(ivl);
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
