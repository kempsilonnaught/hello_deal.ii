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
	Triangulation<2> triangulation;
	GridGenerator::hyper_cube(triangulation); //Figure out what this line does.
	triangulation.refine_global(4); //What does this line do?

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

	for(unsigned int step = 0; step < 5; ++step){
		Triangulation<2>::active_cell_iterator //why no blue
		cell = triangulation.begin_active(),
		endc = triangulation.end();
		for(; cell!=endc; ++cell)
			for(unsigned int v = 0; v < GeometryInfo<2>::vertices_per_cell; ++v){
				const double distance_from_center = center.distance(cell->vertex(v));

				if(std::fabs(distance_from_center - inner_radius) < 1e-10){
					cell->set_refine_flag();
					break;
				}
			}
		triangulation.execute_coarsening_and_refinement();	
	}

	std::ofstream out("grid-2.eps");
	GridOut grid_out;
	grid_out.write_eps(triangulation, out);
	triangulation.set_boundary(0);
}

int main(){
	first_grid();
	second_grid();
}