#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_boundary.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/lac/sparse_matrix.h> //I don't get the comment on the tutorial for this one

#include <deal.II/lac/compressed_sparsity_pattern.h> //This too

#include <deal.II/dofs/dof_renumbering.h> //Why are the degrees of freedom being renumbered?

#include <fstream>

using namespace dealii;

void make_grid(Triangulation<2> &triangulation){ //Triangulation memory address? That's weird.
	
	const Point<2> center(1, 0);
	const double inner_radius = 0.5, outer_radius = 1.0;
	GridGenerator::hyper_shell(triangulation, center, inner_radius, outer_radius, 10);

	static const HyperShellBoundary<2> boundary_description(center);
	triangulation.set_boundary(0, boundary_description);

	for(unsigned int step = 0; step < 5; ++step){
		Triangulation<2>::active_cell_iterator 
		cell = triangulation.begin_active(), endc = triangulation.end();

		for(; cell!=endc; ++cell)
			for(unsigned int v = 0; v < GeometryInfo<2>::vertices_per_cell; ++v){
				const double distance_from_center = center.distance(cell->vertex(v));

				if(std::fabs(distance_from_center - inner_radius < 1e-10){
					cell->set_refine_flag();
					break;
				}
			}
		triangulation.execute_coarsening_and_refinement();	
	}
}


void distribute_dofs(DofHandler<2> &dof_handler){
	static const FE_Q<2> finite_element(1);
	dof_handler.distribut_dofs(finites_element);
	CompressedSparsityPatter compressed_sparsity_pattern(dof_handler.n_dofs(), dof_handler.ndofs());
	
	DoFTools::make_sparsity_pattern(dof_hanhdler, compressed_sparsity_pattern);

	SparsityPattern sparsity_pattern;
	sparsity_pattern.copy_from(compressed_sparsity_pattern);

	std::ofstream out ("sparsity_pattern.1");
	sparsity_pattern.print_gnuplot(out);		
}