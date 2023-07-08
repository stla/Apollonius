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
typedef CGAL::Apollonius_graph_vertex_base_2<Traits, false, Vbi2>  Vb;
typedef CGAL::Triangulation_face_base_with_info_2<int, K>         Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb>              Tds;
//typedef CGAL::Triangulation_2<K,Tds>                        Triangulation;
typedef CGAL::Apollonius_graph_2<Traits, Tds>                     ApolloniusGraph;
typedef ApolloniusGraph::Site_2                                   Site2;
typedef Traits::Line_2                                            Line2;
typedef Traits::Object_2                                          Object2;
typedef ApolloniusGraph::Triangulation_data_structure             Tds;
typedef Tds::Face_handle                                          Face_handle;

// [[Rcpp::export]]
int crash() {

  Point2 p1(-1, -1);
  Point2 p2(-1, 1);
  Point2 p3(1, 1);
  Point2 p4(1, -1);
  Point2 p5(0, 0);

  Site2 s1(p1, 2);
  Site2 s2(p2, 2);
  Site2 s3(p3, 2);
  Site2 s4(p4, 2);
  Site2 s5(p5, 2);

  // make Apollonius graph
  ApolloniusGraph ag;
  ApolloniusGraph::Vertex_handle v1 = ag.insert(s1);
  ApolloniusGraph::Vertex_handle v2 = ag.insert(s2);
  ApolloniusGraph::Vertex_handle v3 = ag.insert(s3);
  ApolloniusGraph::Vertex_handle v4 = ag.insert(s4);
  ApolloniusGraph::Vertex_handle v5 = ag.insert(s5);
  for(int i = 1; i <= 12; i++) {
    double a = 2.0 * M_PI * i / 12.0;
    Point2 pt(2 * cos(a), 2 * sin(a));
    Site2  st(pt, 1);
    ApolloniusGraph::Vertex_handle v = ag.insert(st);
    v->info() = i + 1;
  }
  // validate the Apollonius graph
  assert(ag.is_valid(true, 1));


  int n = 1;
  for(auto f = ag.all_faces_begin(); f != ag.all_faces_end(); f++) {
    f->info() = n++;
  }

  const int nfaces = n - 1;

  int fid = 0;
  for(auto f = ag.all_faces_begin(); f != ag.all_faces_end(); f++) {

    Rcpp::Rcout << "\nface vertices:\n";
    Rcpp::Rcout << "v1: " << f->vertex(0)->site() << "\n";
    Rcpp::Rcout << "v2: " << f->vertex(1)->site() << "\n";
    Rcpp::Rcout << "v3: " << f->vertex(2)->site() << "\n";

    Rcpp::Rcout << "face id: ";
    Rcpp::Rcout << f->info() << "\n";

    for(int i = 0; i < 3; i++) {
      Site2 site = f->vertex(i)->site();
      Point2 pt  = site.point();
      Rcpp::Rcout << site << "\n";
      Rcpp::Rcout << pt << "\n";
    }

    for(int i = 0; i < 3; i++) {
      if(ag.is_infinite(f->neighbor(i))) {
        Rcpp::Rcout << "face neighbor is infinite\n";
      } else {
        Rcpp::Rcout << "face neighbor is finite and its id is "
                    << f->neighbor(i)->info() << "\n";
      }
    }

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
      }
    }

    Rcpp::Rcout << "face dual:\n";
    Rcpp::Rcout << ag.dual(f) << "\n";

    // face dual
    Face_handle fh(f);
    Object2 obj = ag.dual(fh);
    Traits traits = ag.geom_traits();
    Traits::Assign_2 assgn = traits.assign_2_object();
    Site2 site;
    Line2 line;
    if(assgn(site, obj)) {
      // the dual is a site; we store its point and we assign NA to
      // the third element
      Rcpp::Rcout << site << "\n";
    } else if(assgn(line, obj)) {
      // the dual is a line; we store the three parameters of its
      // equation  ax + by + c = 0
      Rcpp::Rcout << line << "\n";
    }

    fid++;
  }


  return 1;
}


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
  Rcpp::NumericMatrix DualPoints(nfaces, 3);
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

    Traits traits = ag.geom_traits();
    Traits::Construct_object_2 cstrct = traits.construct_object_2_object();
    Object2 obj = cstrct(ag.dual(f));
    Traits::Assign_2 assgn = traits.assign_2_object();
    Site2 site;
    Line2 line;

    if(assgn(site, obj)) {
      Point2 pt  = site.point();
      DualPoints(fid, 0) = pt.x();
      DualPoints(fid, 1) = pt.y();
      DualPoints(fid, 2) = Rcpp::NumericVector::get_na();
    } else if(assgn(line, obj)) {
      DualPoints(fid, 0) = line.a();
      DualPoints(fid, 1) = line.b();
      DualPoints(fid, 2) = line.c();
    }

    fid++;
  }

  return Rcpp::List::create(Rcpp::Named("vertices")  = Vertices,
                            Rcpp::Named("neighbors") = Neighbors,
                            Rcpp::Named("cvertex1")  = CommonVertex1,
                            Rcpp::Named("cvertex2")  = CommonVertex2,
                            Rcpp::Named("dpoints")   = DualPoints);
}


// ------------------------------------------------------------------------- //
// ------------------------------------------------------------------------- //
// [[Rcpp::export]]
Rcpp::List test2(Rcpp::NumericMatrix sites, Rcpp::NumericVector radii){

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
  for(auto f = ag.all_faces_begin(); f != ag.all_faces_end(); f++) {
    f->info() = n++;
  }

  const int nfaces = n - 1;

  Rcpp::IntegerMatrix Neighbors(nfaces, 3);
  Rcpp::IntegerMatrix CommonVertex1(nfaces, 3);
  Rcpp::IntegerMatrix CommonVertex2(nfaces, 3);
  Rcpp::NumericMatrix DualPoints(nfaces, 3);
  Rcpp::List          Vertices(nfaces);

  int fid = 0;
  for(auto f = ag.all_faces_begin(); f != ag.all_faces_end(); f++) {

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

    // Rcpp::Rcout << "face dual:\n";
    // Rcpp::Rcout << ag.dual(f) << "\n";

    Face_handle fh(f);

    Traits traits = ag.geom_traits();
    // Traits::Construct_object_2 cstrct = traits.construct_object_2_object();
    Object2 obj = ag.dual(fh);
    Traits::Assign_2 assgn = traits.assign_2_object();
    Site2 site;
    Line2 line;

    if(assgn(site, obj)) {
      Rcpp::Rcout << "face dual:\n";
      Rcpp::Rcout << site << "\n";
      Point2 pt  = site.point();
      DualPoints(fid, 0) = pt.x();
      DualPoints(fid, 1) = pt.y();
      DualPoints(fid, 2) = Rcpp::NumericVector::get_na();
    } else if(assgn(line, obj)) {
      DualPoints(fid, 0) = line.a();
      DualPoints(fid, 1) = line.b();
      DualPoints(fid, 2) = line.c();
    }

    fid++;
  }

  return Rcpp::List::create(Rcpp::Named("vertices")  = Vertices,
                            Rcpp::Named("neighbors") = Neighbors,
                            Rcpp::Named("cvertex1")  = CommonVertex1,
                            Rcpp::Named("cvertex2")  = CommonVertex2,
                            Rcpp::Named("dpoints")   = DualPoints);
}


// -------------------------------------------------------------------------- //
// -------------------------------------------------------------------------- //
void Message(const std::string& msg) {
  SEXP rmsg = Rcpp::wrap(msg);
  Rcpp::message(rmsg);
}


// ------------------------------------------------------------------------- //
// ------------------------------------------------------------------------- //
// [[Rcpp::export]]
Rcpp::List ApolloniusCpp(Rcpp::NumericMatrix sites, Rcpp::NumericVector radii) {

  const int nsites = sites.nrow();

  // make the Apollonius graph
  ApolloniusGraph ag;
  for(int i = 0; i < nsites; i++) {
    Point2 p(sites(i, 0), sites(i, 1));
    Site2 s(p, radii(i));
    ag.insert(s);
  }

  int nfvertices = ag.number_of_vertices();

  if(nfvertices < nsites) {
    Rcpp::warning("Some vertices have been destroyed.");
    std::string msg =
      "Number of remaining vertices: " + std::to_string(nfvertices);
    Message(msg);
  }

  int n_all_vertices = 0;
  int vcounter = 1;
  for(auto v = ag.all_vertices_begin(); v != ag.all_vertices_end(); v++) {
    ApolloniusGraph::Vertex_handle vh(v);
    if(!ag.is_infinite(vh)) {
      vh->info() = vcounter++;
    }
    n_all_vertices++;
  }
  std::string msg =
    "Total number of vertices, including the infinite ones: " +
      std::to_string(n_all_vertices);
  Message(msg);

  // validate the Apollonius graph
  if(!ag.is_valid(true, 1)) {
    Rcpp::stop("The Apollonius graph is not valid.");
  }

  // assign an id to each face
  int n = 1;
  for(auto f = ag.all_faces_begin(); f != ag.all_faces_end(); f++) {
    f->info() = n++;
  }

  const int nfaces = n - 1;

  Rcpp::IntegerMatrix Faces(3, nfaces);
  Rcpp::IntegerMatrix Neighbors(3, nfaces);
  Rcpp::IntegerMatrix CommonVertex1(nfaces, 3);
  Rcpp::IntegerMatrix CommonVertex2(nfaces, 3);
  Rcpp::NumericMatrix Duals(3, nfaces);
  Rcpp::List          Vertices(nfaces);

  int fid = 0;
  for(auto f = ag.all_faces_begin(); f != ag.all_faces_end(); f++) {

    Rcpp::IntegerVector Face = Rcpp::IntegerVector::create(
      f->vertex(0)->info(), f->vertex(1)->info(), f->vertex(2)->info()
    );
    Faces(Rcpp::_, fid) = Face;

    Rcpp::NumericMatrix VS(3, 3);
    for(int i = 0; i < 3; i++) {
      Site2 site = f->vertex(i)->site();
      Point2 pt  = site.point();
      VS(0, i) = pt.x();
      VS(1, i) = pt.y();
      VS(2, i) = site.weight();
    }
    Vertices(fid) = Rcpp::transpose(VS);

    Rcpp::IntegerVector neighbors(3);
    for(int i = 0; i < 3; i++) {
      if(ag.is_infinite(f->neighbor(i))) {
        neighbors(i) = Rcpp::IntegerVector::get_na();
      } else {
        neighbors(i) = f->neighbor(i)->info();
      }
    }
    Neighbors(Rcpp::_, fid) = neighbors;

    for(int j = 0; j < 3; j++) {
      if(!ag.is_infinite(f->neighbor(j))) {
        std::vector<int> commonVertices;
        commonVertices.reserve(2);
        for(int k = 0; k < 3; k++) {
          if(f->neighbor(j)->has_vertex(f->vertex(k))) {
            commonVertices.emplace_back(k + 1);
          }
        }
        CommonVertex1(fid, j) = commonVertices[0];
        CommonVertex2(fid, j) = commonVertices[1];
      }
    }

    // face dual
    Face_handle fh(f);
    Object2 obj = ag.dual(fh);
    Traits traits = ag.geom_traits();
    Traits::Assign_2 assgn = traits.assign_2_object();
    Site2 site;
    Line2 line;
    Rcpp::NumericVector Dual(3);
    if(assgn(site, obj)) {
      // the dual is a site; we store its point and we assign NA to
      // the third element
      Point2 pt  = site.point();
      Dual(0) = pt.x();
      Dual(1) = pt.y();
      Dual(2) = Rcpp::NumericVector::get_na();
    } else if(assgn(line, obj)) {
      // the dual is a line; we store the three parameters of its
      // equation  ax + by + c = 0
      Dual(0) = line.a();
      Dual(1) = line.b();
      Dual(2) = line.c();
    }

    Duals(Rcpp::_, fid++) = Dual;
  }

  return Rcpp::List::create(
    Rcpp::Named("faces")     = Faces,
    Rcpp::Named("vertices")  = Vertices,
    Rcpp::Named("neighbors") = Rcpp::transpose(Neighbors),
    Rcpp::Named("cvertex1")  = CommonVertex1,
    Rcpp::Named("cvertex2")  = CommonVertex2,
    Rcpp::Named("duals")     = Rcpp::transpose(Duals)
  );
}
