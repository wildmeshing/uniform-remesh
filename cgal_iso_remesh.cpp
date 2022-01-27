#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_ratio_stop_predicate.h>
#include <CGAL/Surface_mesh/IO.h>
#include <chrono>
#include <fstream>
#include <iostream>
typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point_3;
typedef CGAL::Surface_mesh<Point_3> Surface_mesh;
typedef boost::graph_traits<Surface_mesh>::halfedge_descriptor halfedge_descriptor;
typedef boost::graph_traits<Surface_mesh>::edge_descriptor edge_descriptor;
namespace SMS = CGAL::Surface_mesh_simplification;
namespace PMP = CGAL::Polygon_mesh_processing;

struct halfedge2edge
{
    halfedge2edge(const Surface_mesh &m, std::vector<edge_descriptor> &edges)
        : m_mesh(m), m_edges(edges)
    {
    }
    void operator()(const halfedge_descriptor &h) const
    {
        m_edges.push_back(edge(h, m_mesh));
    }
    const Surface_mesh &m_mesh;
    std::vector<edge_descriptor> &m_edges;
};

int main(int argc, char *argv[])
{
    Surface_mesh surface_mesh;
    const char *filename = argv[1];
    std::cout << filename << std::endl;
    std::ifstream is(filename);
    if (!is || !(is >> surface_mesh))
    {
        std::cerr << "Failed to read input mesh: " << filename << std::endl;
        return EXIT_FAILURE;
    }
    if (!CGAL::is_triangle_mesh(surface_mesh))
    {
        std::cerr << "Input geometry is not triangulated." << std::endl;
        return EXIT_FAILURE;
    }
    double target_edge_length = (argc > 1) ? std::stod(std::string(argv[2])) : 0.04;
    unsigned int nb_iter = 2;

    std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now();
    std::vector<edge_descriptor> border;
    // PMP::border_halfedges(faces(surface_mesh), surface_mesh, boost::make_function_output_iterator(halfedge2edge(surface_mesh, border)));
    // PMP::split_long_edges(border, target_edge_length, surface_mesh);
    std::cout << "done." << std::endl;
    PMP::isotropic_remeshing(faces(surface_mesh), target_edge_length, surface_mesh,
                             PMP::parameters::number_of_iterations(nb_iter)
                                 .protect_constraints(true));

    // In this example, the simplification stops when the number of undirected edges
    // drops below 10% of the initial count
    //double stop_ratio = (argc > 2) ? std::stod(argv[2]) : 0.1;
    //SMS::Count_ratio_stop_predicate<Surface_mesh> stop(stop_ratio);
    //int r = SMS::edge_collapse(surface_mesh, stop);
    std::chrono::steady_clock::time_point end_time = std::chrono::steady_clock::now();
    std::cout << "\nFinished!\n"
              //   << r << " edges removed.\n"
              << surface_mesh.number_of_vertices() << " final vertices.\n"
              << surface_mesh.number_of_faces() << " final faces.\n";
    std::cout << "Time elapsed: " << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count() << "ms" << std::endl;
    //CGAL::IO::write_polygon_mesh((argc > 3) ? argv[3] : "out.off", surface_mesh, CGAL::parameters::stream_precision(17));
    std::ofstream out("remesh_out.off");
    if (!out || !(out << surface_mesh))
    {
        std::cerr << "Failed to write mesh: " << filename << std::endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
