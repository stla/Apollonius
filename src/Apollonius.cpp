#include <Rcpp.h>
#define CGAL_EIGEN3_ENABLED 1
#include <CGAL/Apollonius_graph_2.h>
#include <CGAL/Apollonius_graph_traits_2.h>
#include <CGAL/Apollonius_graph_vertex_base_2.h>
#include <CGAL/Apollonius_site_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
//#include <CGAL/Triangulation_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel       K;
typedef K::Point_2                                                Point2;
typedef CGAL::Apollonius_graph_traits_2<K>                        Traits;
typedef CGAL::Triangulation_vertex_base_with_info_2<int, K>       Vbi2;
typedef CGAL::Apollonius_graph_vertex_base_2<Traits, false, Vbi2> Vb;
typedef CGAL::Triangulation_face_base_with_info_2<int, K>         Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb>              Tds;
//typedef CGAL::Triangulation_2<K,Tds>                        Triangulation;
typedef CGAL::Apollonius_graph_2<Traits, Tds>                     ApolloniusGraph;
typedef ApolloniusGraph::Site_2                                   Site2;
//typedef ApolloniusGraph::Triangulation_data_structure       Tds;
//typedef Triangulation::Vertex_handle                        Vertex_handle;

// [[Rcpp::export]]
Rcpp::List test(Rcpp::NumericMatrix sites, Rcpp::NumericVector radii){

  const int nsites = sites.nrow();

  ApolloniusGraph ag;
  for(int i = 0; i < nsites; i++) {
    Point2 p(sites(i, 0), sites(i, 1));
    Site2 s(p, radii(i));
    ApolloniusGraph::Vertex_handle v = ag.insert(s);
    v->info() = i + 1;
  }

  // validate the Apollonius graph
  if(!ag.is_valid(true, 1)) {
    Rcpp::stop("The Apollonius graph is not valid.");
  }

  int n = 1;
  for(auto f = ag.finite_faces_begin(); f < ag.finite_faces_end(); f++) {
    f->info() = n++;
  }

  const int nfaces = n - 1;

  Rcpp::IntegerMatrix Neighbors(nfaces, 3);
  Rcpp::IntegerMatrix CommonVertex1(nfaces, 3);
  Rcpp::IntegerMatrix CommonVertex2(nfaces, 3);
  Rcpp::NumericMatrix DualPoints(nfaces, 2);
  Rcpp::List          Vertices(nfaces);

  int fid = 0;
  for(auto f = ag.finite_faces_begin(); f < ag.finite_faces_end(); f++) {

    Rcpp::Rcout << "\nface vertices:\n";
    Rcpp::Rcout << "v1: " << f->vertex(0)->site() << "\n";
    Rcpp::Rcout << "v2: " << f->vertex(1)->site() << "\n";
    Rcpp::Rcout << "v3: " << f->vertex(2)->site() << "\n";

    Rcpp::NumericMatrix VS(3, 3);
    for(int i = 0; i < 3; i++) {
      Site2 site = f->vertex(i)->site();
      Point2 pt  = site.point();
      VS(i, 0) = pt.x();
      VS(i, 1) = pt.y();
      VS(i, 2) = site.weight();
    }
    // VS(0, 0) = f->vertex(0)->site().point().x();
    // VS(0, 1) = f->vertex(0)->site().point().y();
    // VS(0, 2) = f->vertex(0)->site().weight();
    // VS(1, 0) = f->vertex(1)->site().point().x();
    // VS(1, 1) = f->vertex(1)->site().point().y();
    // VS(1, 2) = f->vertex(1)->site().weight();
    // VS(2, 0) = f->vertex(2)->site().point().x();
    // VS(2, 1) = f->vertex(2)->site().point().y();
    // VS(2, 2) = f->vertex(2)->site().weight();
    Vertices(fid) = VS;

    Rcpp::Rcout << "face id: ";
    Rcpp::Rcout << f->info() << "\n";

    Rcpp::IntegerVector neighbors(3);
    for(int i = 0; i < 3; i++) {
      if(ag.is_infinite(f->neighbor(i))) {
        neighbors(i) = Rcpp::IntegerVector::get_na();
      } else {
        neighbors(i) = f->neighbor(i)->info();
      }
    }
    Neighbors(fid, Rcpp::_) = neighbors;

    for(int j = 0; j < 3; j++) {
      if(!ag.is_infinite(f->neighbor(j))) {
        Rcpp::Rcout << "face neighbor " << j + 1 << ": ";
        Rcpp::Rcout << f->neighbor(j)->info() << "\n";
        std::vector<int> commonVertices;
        commonVertices.reserve(2);
        for(int k = 0; k < 3; k++) {
          if(f->neighbor(j)->has_vertex(f->vertex(k))) {
            commonVertices.emplace_back(k + 1);
          }
        }
        Rcpp::Rcout << "common vertices: "
                    << commonVertices[0] << ", " << commonVertices[1] << "\n";
        CommonVertex1(fid, j) = commonVertices[0];
        CommonVertex2(fid, j) = commonVertices[1];
      }
    }

    Rcpp::Rcout << "face dual:\n";
    Rcpp::Rcout << ag.dual(f) << "\n";

    Site2 site = ag.dual(f);
    Point2 pt  = site.point();
    DualPoints(fid, 0) = pt.x();
    DualPoints(fid, 1) = pt.y();
    fid++;
  }

  return Rcpp::List::create(Rcpp::Named("vertices")  = Vertices,
                            Rcpp::Named("neighbors") = Neighbors,
                            Rcpp::Named("cvertex1")  = CommonVertex1,
                            Rcpp::Named("cvertex2")  = CommonVertex2,
                            Rcpp::Named("dpoints")   = DualPoints);
}
