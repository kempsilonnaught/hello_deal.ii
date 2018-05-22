#include <deal.II/grid/tria.h>

#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/grid/tria_boundary_lib.h>

#include <deal.II/grid/grid_out.h>

#include <fstream>
#include <cmath>

using namespace dealii;

void first_grid(){
	Trangulation<2> triangulation;
	GridGenerator::hyper_cube(triangulation); //Figure out what this line does.
	trangulation.refine_global(4); //What does this line do?

	std::ofstream out("grid-1.eps");
	GridOut grid_out; //ask Jef why GridOut isn't blue
	grid_out.write_eps(triangulation, out);
}

void second_grid(){
	Triangulation<2> triangulation;
	const Point<2> center(1,0);
	const double inner_radius = .5, outer_radius = 1.0;
	GridGenerator::hyper_shell(triangulation, center, inner_radius, outer_radius, 10); //what is hyper_shell

	const HyperShellBoundary<2> boundary_description(center);
	triangulation.set_boundary(0, boundary_description);
}