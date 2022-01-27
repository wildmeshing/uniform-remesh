#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/border.h>
//#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <boost/iterator/function_output_iterator.hpp>
#include <fstream>
#include <vector>
#include <string>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Triangulation_utils_3.h>
#include <CGAL/Real_timer.h>

#include <CGAL/IO/read_off_points.h>
#include <CGAL/Surface_mesh/IO/OFF.h>
#include <CGAL/IO/read_ply_points.h>

#include "CLI11.hpp"

#include <iostream>
#include <fstream>
#include <vector>
//#include <set>
//#include <unordered_set>
//#include <bitset>
#include <queue>
//#include <map>
//#include <chrono>
//#include <random>
#include <algorithm>

//#include <boost/math/tools/minima.hpp>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::FT FT;

typedef CGAL::Triangulation_vertex_base_with_info_3<int, K> Vb;
typedef CGAL::Triangulation_cell_base_with_info_3<int, K> Cb;
typedef CGAL::Triangulation_data_structure_3<Vb, Cb> Tds;

typedef CGAL::Delaunay_triangulation_3<K, Tds> DT;
typedef CGAL::Triangulation_3<K, Tds> Triangulation;

typedef Triangulation::Finite_vertices_iterator Finite_vertices_iterator;
typedef Triangulation::Finite_edges_iterator Finite_edges_iterator;
typedef Triangulation::Finite_facets_iterator Finite_facets_iterator;
typedef Triangulation::Finite_cells_iterator Finite_cells_iterator;

typedef Triangulation::Point Point;
typedef Triangulation::Facet Facet;
typedef Triangulation::Edge Edge;
typedef Triangulation::Vertex_handle VH;
typedef Triangulation::Cell_handle CH;
typedef Triangulation::Cell_circulator CC;
typedef Triangulation::Facet_circulator FC;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3> Mesh;
typedef boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;
typedef boost::graph_traits<Mesh>::edge_descriptor edge_descriptor;

typedef Mesh::Vertex_index vertex_descriptor;
typedef Mesh::Face_index face_descriptor;

namespace PMP = CGAL::Polygon_mesh_processing;

struct halfedge2edge
{
    halfedge2edge(const Mesh &m, std::vector<edge_descriptor> &edges)
        : m_mesh(m), m_edges(edges)
    {
    }
    void operator()(const halfedge_descriptor &h) const
    {
        m_edges.push_back(edge(h, m_mesh));
    }
    const Mesh &m_mesh;
    std::vector<edge_descriptor> &m_edges;
};
int main(int argc, char *argv[])
{
    std::string infilename = "default.off";
    std::ifstream in(infilename);
    CGAL::IO::read_OFF(std::istream & in,
                       Mesh & sm);

    // std::vector<Point_3> P;
    // if (!CGAL::read_off_points(in, std::back_inserter(P)))
    // {
    //     std::cerr << "Error: cannot read off file " << infilename << std::endl;
    //     exit(0);
    // }

    // int n = P.size();
    // std::cout << "Read " << n << " points from file " << infilename << std::endl;

    double target_edge_length = (argc > 2) ? std::stod(std::string(argv[2])) : 0.04;
    unsigned int nb_iter = 3;
    // std::cout << "Split border...";
    // std::vector<edge_descriptor> border;
    // PMP::border_halfedges(faces(mesh), mesh, boost::make_function_output_iterator(halfedge2edge(mesh, border)));
    // PMP::split_long_edges(border, target_edge_length, mesh);
    // std::cout << "done." << std::endl;
    // std::cout << "Start remeshing of " << filename
    //           << " (" << num_faces(mesh) << " faces)..." << std::endl;
    PMP::isotropic_remeshing(faces(P), target_edge_length, P,
                             PMP::parameters::number_of_iterations(nb_iter)
                                 .protect_constraints(true)); //i.e. protect border, here
    std::cout << "Remeshing done." << std::endl;
    return 0;
}
