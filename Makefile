# Compiler config
DC ?= gdc
DCFLAGS ?= -g3 -O0 -I.. -Wall -funittest

LDFLAGS ?= 
OBJEXT ?= 

# Target
OBJ_NAME = dgeom

.PHONY: $(OBJ_NAME) clean

$(OBJ_NAME) : $(SOURCES)
	$(DC) $(DCFLAGS) $(LDFLAGS) -o $(OBJ_NAME) $(SOURCES)

clean:
	rm $(OBJ_NAME)

SOURCES = \
	source/app.d			\
	source/geom/affine.d		\
	source/geom/angle.d		\
	source/geom/bezier.d		\
	source/geom/bezier_curve.d	\
	source/geom/bezier_fitting.d	\
	source/geom/choose.d		\
	source/geom/coord.d		\
	source/geom/crossing.d		\
	source/geom/curve.d		\
	source/geom/d2.d		\
	source/geom/ellipse.d		\
	source/geom/elliptical_arc.d	\
	source/geom/interval.d		\
	source/geom/intpoint.d		\
	source/geom/linear.d		\
	source/geom/nearest_time.d	\
	source/geom/path.d		\
	source/geom/path_builder.d	\
	source/geom/path_intersection.d	\
	source/geom/path_sequence.d	\
	source/geom/path_outline.d	\
	source/geom/piecewise.d		\
	source/geom/point.d		\
	source/geom/rect.d		\
	source/geom/region.d		\
	source/geom/sbasis.d		\
	source/geom/sbasis_curve.d	\
	source/geom/sbasis_roots.d	\
	source/geom/solve_bezier.d	\
	source/geom/sbasis_to_bezier.d	\
	source/geom/svg_elliptical_arc.d\
	source/geom/transforms.d	
