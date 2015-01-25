# Compiler config
DC ?= gdc
DCFLAGS ?= -g3 -O0 -I.. -Wall -funittest

LDFLAGS ?= 
OBJEXT ?= 

# Target
OBJ_NAME = lib2geom

.PHONY: $(OBJ_NAME) clean

$(OBJ_NAME) : $(SOURCES)
	$(DC) $(DCFLAGS) $(LDFLAGS) -o $(OBJ_NAME) $(SOURCES)

clean:
	rm $(OBJ_NAME)

SOURCES = \
	affine.d	\
	d2.d		\
	bezier.d	\
	bezier_curve.d	\
	choose.d	\
	coord.d		\
	curve.d		\
	interval.d	\
	intpoint.d	\
	linear.d	\
	main.d		\
	piecewise.d	\
	point.d		\
	rect.d		\
	sbasis.d	\
	sbasis_roots.d	\
	solve_bezier.d	\
	sbasis_to_bezier.d	\
	transforms.d	
