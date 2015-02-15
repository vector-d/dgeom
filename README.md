# dgeom :new:
Based on [lib2geom](http://lib2geom.sourceforge.net/), dgeom is completely rewritten in D, actively developed, and designed to build faster, run quicker, and take advantage of the many superior qualities D has to C++.

Written and maintained by [liamwhite](http://github.com/liamwhite) and [Mihail-K](http://github.com/Mihail-K).
<hr>
##Building
**Dependencies**
* [dub package manager](https://github.com/D-Programming-Language/dub)
* D compiler
  * Digital Mars D compiler (dmd) 2.xx or newer is **preferred**
  * <strike>LLVM D2 compiler (ldc2) 3.4x or newer</strike> (Current [ldc2 bug](https://github.com/ldc-developers/ldc/issues/837) renders uncompilable)
  * GNU D compiler (gdc) 4.8x or newer

> GNU make can be used to build dgeom if necessary, but requires different makefile settings for different compilers

**Steps**

1. Install dependencies (if needed)
2. Clone directory from GitHub (if needed)
3. `cd` into dgeom directory
4. `$ dub build`

<hr>
###License
dgeom is licensed under GNU Public License Version 3 (GPLv3)
