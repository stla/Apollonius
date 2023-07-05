#include <Rcpp.h>
#define CGAL_EIGEN3_ENABLED 1
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Apollonius_graph_2.h>
#include <CGAL/Apollonius_site_2.h>
#include <CGAL/Apollonius_graph_traits_2.h>
//#include <CGAL/Triangulation_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2                                          Point2;
typedef CGAL::Apollonius_graph_traits_2<K>                  Traits;
typedef CGAL::Apollonius_graph_2<Traits>                    ApolloniusGraph;
typedef ApolloniusGraph::Site_2                             Site2;
//typedef ApolloniusGraph::Triangulation_data_structure       Tds;
//typedef CGAL::Triangulation_2<K,Tds>                        Triangulation;
//typedef Triangulation::Vertex_handle                        Vertex_handle;

// [[Rcpp::export]]
void test() {

  Point2 p1(0, 0);
  Point2 p2(4, 1);
  Point2 p3(2, 4);
  Point2 p4(7, 4);

  Site2 s1(p1, 1);
  Site2 s2(p2, 1.5);
  Site2 s3(p3, 1);
  Site2 s4(p4, 1);

  ApolloniusGraph ag;
  ag.insert(s1);
  ag.insert(s2);
  ag.insert(s3);
  ag.insert(s4);
  // validate the Apollonius graph
  assert(ag.is_valid(true, 1));

  for(auto f = ag.finite_faces_begin(); f < ag.finite_faces_end(); f++) {
    Rcpp::Rcout << "\nface vertices:\n";
    Rcpp::Rcout << f->vertex(0)->site() << "\n";
    Rcpp::Rcout << f->vertex(1)->site() << "\n";
    Rcpp::Rcout << f->vertex(2)->site() << "\n";
    Rcpp::Rcout << "face dual:\n";
    Rcpp::Rcout << ag.dual(f) << "\n";
  }

}
