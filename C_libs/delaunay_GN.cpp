
//============================================================================
// delaunay_GN.cpp -*- C++ -*-; generic implementation for reg. triangulations
//
// Copyright (C) 2007, 2008, 2009, 2010 Gerd Neudert
//
// This file is part of fconv.
// fconv is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
// version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
// more details.
//
// You should have received a copy of the GNU General Public License along
// with this program; if not, see <http://www.gnu.org/licenses/>.
//----------------------------------------------------------------------------
//
// author:      Gerd Neudert
//              gneudert(place_at_here)web.de
// supervisor:  Prof. Dr. Gerhard Klebe
//              klebe(place_at_here)staff.uni-marburg.de
// Department of Pharmaceutical Chemistry, Philipps-University of Marburg
//
// This library is part of the program fconv. Please see the
// documentation in the file fconv.cpp first!
// This library provides a generic interface for regular triangulations.
//============================================================================


#include"delaunay_GN.h"


#if defined (_SPEEDTEST)
#include"timer_GN.hpp"
TIMER_GN<double> t_triang;
TIMER_GN<double> t_locate;
TIMER_GN<double> t_flip;
TIMER_GN<double> t_regular;
TIMER_GN<double> t_regtest;
TIMER_GN<double> t_deltet;
TIMER_GN<double> t_proof;
#endif

//!***************************************************************************************************
//!***************************************************************************************************

ostream &operator<<(ostream &os,stl_ptr<TETRAHEDRON> const& t) {
    os << "Tetrahedron from generation " << t->generation << ": next.size() = "
       << t->next.size() << "   prev.size() = " << t->prev.size() << endl;
    os << "Coords = " << t->p0->coord << " -- " << t->p1->coord << " -- " << t->p2->coord << " -- " << t->p3->coord << endl;
    return os;
}

ostream &operator<<(ostream &os,stl_ptr<L_FACET> const& f) {
    os << "Facet " << f->id << ":" << endl;
    os << "t1: " << f->t1->id << "     t2: " << f->t2->id << endl;
    return os;
}


//!===============================================================================
//!Definition der Klasse DT_POINT:
DT_POINT::DT_POINT(int const& i,double const& x,double const& y,double const& z,double const& weight,
                   double const& radius,void *at) :
                   id(i),coord(vec3d<double>(x,y,z)),w(weight),r(radius),myobj(at) {calc_vol();}

DT_POINT::DT_POINT(DT_POINT const& ref) {
    id = ref.id;
    coord = ref.coord;
    w = ref.w;
    r = ref.r;
    myobj = ref.myobj;
    vol = ref.vol;
}

DT_POINT::~DT_POINT() {}

void DT_POINT::calc_vol() {
    if (w < 0.00001) vol = 0.;
    else {
        if (r < 0.00001) r = sqrt(w);
        vol = 4. * acos(-1.) * w * r / 3.;
    }
}
//!===============================================================================


//!===============================================================================
//!Definition der Klasse TETRAHEDRON:
TETRAHEDRON::TETRAHEDRON(stl_ptr<DT_POINT> const& v0,stl_ptr<DT_POINT> const& v1,
                         stl_ptr<DT_POINT> const& v2,stl_ptr<DT_POINT> const& v3,
                         int64_t const& currgen,int64_t const& i) :
                         p0(v0), p1(v1), p2(v2), p3(v3), generation(currgen), id(i) {
    if (!calc_sphere()) {sphere_square_radius = 999999999999.; sphere_middle = p0->coord;}
}

TETRAHEDRON::TETRAHEDRON(DT_POINT& v0,DT_POINT& v1,
                         DT_POINT& v2,DT_POINT& v3,int64_t const& currgen,int64_t const& i) :
                         p0(v0), p1(v1), p2(v2), p3(v3), generation(currgen), id(i) {
    if (!calc_sphere()) {sphere_square_radius = 999999999999.; sphere_middle = p0->coord;}
}

TETRAHEDRON::TETRAHEDRON(TETRAHEDRON const& ref) {
    p0 = ref.p0; p1 = ref.p1; p2=ref.p2; p3 = ref.p3;
    sphere_middle = ref.sphere_middle;
    sphere_square_radius = ref.sphere_square_radius;
    generation = ref.generation;
    id = ref.id;
}

TETRAHEDRON::~TETRAHEDRON() {}

bool TETRAHEDRON::operator<(TETRAHEDRON const& rechts) {
    if (sphere_square_radius < rechts.sphere_square_radius) return true;
    return false;
}

void TETRAHEDRON::visualize(int const& state,bool const& vdw_rad) {
    ofstream f_out;
    ostringstream os;
    os << "tetrahedron" << state << ".py";
    string fname = os.str();
    f_out.open(fname.c_str());
    
    float rad = 0.1;
    
    f_out << "p0=[7.0,"
        << p0->coord[0] << "," << p0->coord[1] << "," << p0->coord[2];
    if (vdw_rad) f_out << "," << p0->r << "]" << endl;
    else f_out << ",0.3]" << endl;
    f_out << "cmd.load_cgo(p0,'p0'," << state << ")" << endl;
    
    f_out << "p1=[7.0,"
        << p1->coord[0] << "," << p1->coord[1] << "," << p1->coord[2];
    if (vdw_rad) f_out << "," << p1->r << "]" << endl;
    else f_out << ",0.3]" << endl;
    f_out << "cmd.load_cgo(p1,'p1'," << state << ")" << endl;
    
    f_out << "p2=[7.0,"
        << p2->coord[0] << "," << p2->coord[1] << "," << p2->coord[2];
    if (vdw_rad) f_out << "," << p2->r << "]" << endl;
    else f_out << ",0.3]" << endl;
    f_out << "cmd.load_cgo(p2,'p2'," << state << ")" << endl;
    
    f_out << "p3=[7.0,"
        << p3->coord[0] << "," << p3->coord[1] << "," << p3->coord[2];
    if (vdw_rad) f_out << "," << p3->r << "]" << endl;
    else f_out << ",0.3]" << endl;
    f_out << "cmd.load_cgo(p3,'p3'," << state << ")" << endl;
    
    f_out << "self=[9.0,"
            << p0->coord[0] << "," << p0->coord[1] << "," << p0->coord[2] << ","
        << p1->coord[0] << "," << p1->coord[1] << "," << p1->coord[2] << ","
        << rad << ",1.0,1.0,0.0,1.0,1.0,0.0,9.0,"
        
        << p0->coord[0] << "," << p0->coord[1] << "," << p0->coord[2] << ","
        << p2->coord[0] << "," << p2->coord[1] << "," << p2->coord[2] << ","
        << rad << ",1.0,1.0,0.0,1.0,1.0,0.0,9.0,"
        
        << p0->coord[0] << "," << p0->coord[1] << "," << p0->coord[2] << ","
        << p3->coord[0] << "," << p3->coord[1] << "," << p3->coord[2] << ","
        << rad << ",1.0,1.0,0.0,1.0,1.0,0.0,9.0,"
        
        << p1->coord[0] << "," << p1->coord[1] << "," << p1->coord[2] << ","
        << p2->coord[0] << "," << p2->coord[1] << "," << p2->coord[2] << ","
        << rad << ",1.0,1.0,0.0,1.0,1.0,0.0,9.0,"
        
        << p2->coord[0] << "," << p2->coord[1] << "," << p2->coord[2] << ","
        << p3->coord[0] << "," << p3->coord[1] << "," << p3->coord[2] << ","
        << rad << ",1.0,1.0,0.0,1.0,1.0,0.0,9.0,"
        
        << p3->coord[0] << "," << p3->coord[1] << "," << p3->coord[2] << ","
        << p1->coord[0] << "," << p1->coord[1] << "," << p1->coord[2] << ","
        << rad << ",1.0,1.0,0.0,1.0,1.0,0.0]" << endl;
    f_out << "cmd.load_cgo(self,'self'," << state << ")" << endl;
    
    if (!p0n.zero()) {
    f_out << "p0n=[9.0,"
              << p0n->p0->coord[0] << "," << p0n->p0->coord[1] << "," << p0n->p0->coord[2] << ","
              << p0n->p1->coord[0] << "," << p0n->p1->coord[1] << "," << p0n->p1->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0,9.0,"
              
              << p0n->p0->coord[0] << "," << p0n->p0->coord[1] << "," << p0n->p0->coord[2] << ","
              << p0n->p2->coord[0] << "," << p0n->p2->coord[1] << "," << p0n->p2->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0,9.0,"
              
              << p0n->p0->coord[0] << "," << p0n->p0->coord[1] << "," << p0n->p0->coord[2] << ","
              << p0n->p3->coord[0] << "," << p0n->p3->coord[1] << "," << p0n->p3->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0,9.0,"
              
              << p0n->p1->coord[0] << "," << p0n->p1->coord[1] << "," << p0n->p1->coord[2] << ","
              << p0n->p2->coord[0] << "," << p0n->p2->coord[1] << "," << p0n->p2->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0,9.0,"
              
              << p0n->p2->coord[0] << "," << p0n->p2->coord[1] << "," << p0n->p2->coord[2] << ","
              << p0n->p3->coord[0] << "," << p0n->p3->coord[1] << "," << p0n->p3->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0,9.0,"
              
              << p0n->p3->coord[0] << "," << p0n->p3->coord[1] << "," << p0n->p3->coord[2] << ","
              << p0n->p1->coord[0] << "," << p0n->p1->coord[1] << "," << p0n->p1->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0]" << endl;
    f_out << "cmd.load_cgo(p0n,'p0n'," << state << ")" << endl;
    }
    
    if (!p1n.zero()) {
    f_out << "p1n=[9.0,"
              << p1n->p0->coord[0] << "," << p1n->p0->coord[1] << "," << p1n->p0->coord[2] << ","
              << p1n->p1->coord[0] << "," << p1n->p1->coord[1] << "," << p1n->p1->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0,9.0,"
              
              << p1n->p0->coord[0] << "," << p1n->p0->coord[1] << "," << p1n->p0->coord[2] << ","
              << p1n->p2->coord[0] << "," << p1n->p2->coord[1] << "," << p1n->p2->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0,9.0,"
              
              << p1n->p0->coord[0] << "," << p1n->p0->coord[1] << "," << p1n->p0->coord[2] << ","
              << p1n->p3->coord[0] << "," << p1n->p3->coord[1] << "," << p1n->p3->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0,9.0,"
              
              << p1n->p1->coord[0] << "," << p1n->p1->coord[1] << "," << p1n->p1->coord[2] << ","
              << p1n->p2->coord[0] << "," << p1n->p2->coord[1] << "," << p1n->p2->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0,9.0,"
              
              << p1n->p2->coord[0] << "," << p1n->p2->coord[1] << "," << p1n->p2->coord[2] << ","
              << p1n->p3->coord[0] << "," << p1n->p3->coord[1] << "," << p1n->p3->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0,9.0,"
              
              << p1n->p3->coord[0] << "," << p1n->p3->coord[1] << "," << p1n->p3->coord[2] << ","
              << p1n->p1->coord[0] << "," << p1n->p1->coord[1] << "," << p1n->p1->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0]" << endl;
    f_out << "cmd.load_cgo(p1n,'p1n'," << state << ")" << endl;
    }
    
    if (!p2n.zero()) {
    f_out << "p2n=[9.0,"
              << p2n->p0->coord[0] << "," << p2n->p0->coord[1] << "," << p2n->p0->coord[2] << ","
              << p2n->p1->coord[0] << "," << p2n->p1->coord[1] << "," << p2n->p1->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0,9.0,"
              
              << p2n->p0->coord[0] << "," << p2n->p0->coord[1] << "," << p2n->p0->coord[2] << ","
              << p2n->p2->coord[0] << "," << p2n->p2->coord[1] << "," << p2n->p2->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0,9.0,"
              
              << p2n->p0->coord[0] << "," << p2n->p0->coord[1] << "," << p2n->p0->coord[2] << ","
              << p2n->p3->coord[0] << "," << p2n->p3->coord[1] << "," << p2n->p3->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0,9.0,"
              
              << p2n->p1->coord[0] << "," << p2n->p1->coord[1] << "," << p2n->p1->coord[2] << ","
              << p2n->p2->coord[0] << "," << p2n->p2->coord[1] << "," << p2n->p2->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0,9.0,"
              
              << p2n->p2->coord[0] << "," << p2n->p2->coord[1] << "," << p2n->p2->coord[2] << ","
              << p2n->p3->coord[0] << "," << p2n->p3->coord[1] << "," << p2n->p3->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0,9.0,"
              
              << p2n->p3->coord[0] << "," << p2n->p3->coord[1] << "," << p2n->p3->coord[2] << ","
              << p2n->p1->coord[0] << "," << p2n->p1->coord[1] << "," << p2n->p1->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0]" << endl;
    f_out << "cmd.load_cgo(p2n,'p2n'," << state << ")" << endl;
    }
    
    if (!p3n.zero()) {
    f_out << "p3n=[9.0,"
              << p3n->p0->coord[0] << "," << p3n->p0->coord[1] << "," << p3n->p0->coord[2] << ","
              << p3n->p1->coord[0] << "," << p3n->p1->coord[1] << "," << p3n->p1->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0,9.0,"
              
              << p3n->p0->coord[0] << "," << p3n->p0->coord[1] << "," << p3n->p0->coord[2] << ","
              << p3n->p2->coord[0] << "," << p3n->p2->coord[1] << "," << p3n->p2->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0,9.0,"
              
              << p3n->p0->coord[0] << "," << p3n->p0->coord[1] << "," << p3n->p0->coord[2] << ","
              << p3n->p3->coord[0] << "," << p3n->p3->coord[1] << "," << p3n->p3->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0,9.0,"
              
              << p3n->p1->coord[0] << "," << p3n->p1->coord[1] << "," << p3n->p1->coord[2] << ","
              << p3n->p2->coord[0] << "," << p3n->p2->coord[1] << "," << p3n->p2->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0,9.0,"
              
              << p3n->p2->coord[0] << "," << p3n->p2->coord[1] << "," << p3n->p2->coord[2] << ","
              << p3n->p3->coord[0] << "," << p3n->p3->coord[1] << "," << p3n->p3->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0,9.0,"
              
              << p3n->p3->coord[0] << "," << p3n->p3->coord[1] << "," << p3n->p3->coord[2] << ","
              << p3n->p1->coord[0] << "," << p3n->p1->coord[1] << "," << p3n->p1->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0]" << endl;
    f_out << "cmd.load_cgo(p3n,'p3n'," << state << ")" << endl;
    }
    
    f_out.close();
}

void TETRAHEDRON::visualize_sphere(string const& sname) {
    ofstream f_out;
    string fname = sname + ".py";
    f_out.open(fname.c_str());
    
    f_out << sname << "=[";
        
    f_out << "7.0," << sphere_middle[0] << "," << sphere_middle[1] << "," << sphere_middle[2]
          << "," << sqrt(sphere_square_radius) << "]" << endl;
    
    f_out << "cmd.load_cgo(" << sname << ",'" << sname << "',1)" << endl;
    
    f_out.close();
}

bool TETRAHEDRON::is_regular(vector<stl_ptr<DT_POINT> >& dt_points) {
    #if defined (_SPEEDTEST)
    t_regtest.continue_t();
    #endif
    double testval = sphere_square_radius - (sphere_square_radius / 1000000.);
    for (dtp_vec it=dt_points.begin(); it!=dt_points.end(); ++it) {
        //ausschliessen, dass der Punkt zum tet gehoert -> gucken ob der Punkt in der umkugel liegt
        if ((*it).equal_addr(p0)) continue;
        if ((*it).equal_addr(p1)) continue;
        if ((*it).equal_addr(p2)) continue;
        if ((*it).equal_addr(p3)) continue;
        float d2 = get_square_distance((*it)->coord,sphere_middle)-(*it)->w;
        if (d2 < testval) {
            return false;
        }
    }
    #if defined (_SPEEDTEST)
    t_regtest.stop_t();
    #endif
    return true;
}

bool TETRAHEDRON::coplanar() {
    return is_coplanar(p0->coord,p1->coord,p2->coord,p3->coord);
}

bool TETRAHEDRON::shares_point(stl_ptr<TETRAHEDRON> const& tet) {
    if (p0.equal_addr(tet->p0)) return true;
    else if (p0.equal_addr(tet->p1)) return true;
    else if (p0.equal_addr(tet->p2)) return true;
    else if (p0.equal_addr(tet->p3)) return true;
    
    if (p1.equal_addr(tet->p0)) return true;
    else if (p1.equal_addr(tet->p1)) return true;
    else if (p1.equal_addr(tet->p2)) return true;
    else if (p1.equal_addr(tet->p3)) return true;
    
    if (p2.equal_addr(tet->p0)) return true;
    else if (p2.equal_addr(tet->p1)) return true;
    else if (p2.equal_addr(tet->p2)) return true;
    else if (p2.equal_addr(tet->p3)) return true;
    
    if (p3.equal_addr(tet->p0)) return true;
    else if (p3.equal_addr(tet->p1)) return true;
    else if (p3.equal_addr(tet->p2)) return true;
    else if (p3.equal_addr(tet->p3)) return true;
    
    return false;
}

bool TETRAHEDRON::has_infinit_point() {
    //! Die Punkte des allumfassenden Tetraeder haben einen Index <0
    if (p0->id < 0 || p1->id < 0 || p2->id < 0  || p3->id < 0) return true;
    else return false;
}

bool TETRAHEDRON::shares_facets(stl_ptr<TETRAHEDRON> const& tet) {
    if (p0n->id == tet->id) return true;
    if (p1n->id == tet->id) return true;
    if (p2n->id == tet->id) return true;
    if (p3n->id == tet->id) return true;
    return false;
}

bool TETRAHEDRON::calc_sphere() {
    //!Kugel um 4 Punkte (des Tetraeders) berechnen
    //!Nach "P.Bourke, Equation of a sphere from 4 points on the surface"
    //!Modifiziert fuer gewichtete Punkte (siehe linalg_GN.hpp)
    return weighted_points_to_sphere(p0->coord,p0->w,p1->coord,p1->w,p2->coord,p2->w,p3->coord,p3->w,
                                     sphere_middle,sphere_square_radius);
}

void TETRAHEDRON::calc_volume() {
    volume = (fabs(triple_product(p0->coord,p1->coord,p2->coord,p3->coord)) / 6.);
}

void TETRAHEDRON::calc_weighted_volume() {
    volume = (fabs(triple_product(p0->coord,p1->coord,p2->coord,p3->coord)) / 6.);
    double check = volume;
    
    //! Die Schitte der vdW-Kugeln mit dem Tetraeder abziehen:
    volume -= get_sphere_tetrahedron_intersection(p0->coord,p1->coord,p2->coord,p3->coord,p0->vol);
    volume -= get_sphere_tetrahedron_intersection(p1->coord,p0->coord,p2->coord,p3->coord,p1->vol);
    volume -= get_sphere_tetrahedron_intersection(p2->coord,p0->coord,p1->coord,p3->coord,p2->vol);
    volume -= get_sphere_tetrahedron_intersection(p3->coord,p0->coord,p1->coord,p2->coord,p3->vol);
    
    
    //!==============================================================
    //! Jetzt muss noch der eventuelle Overlap je zweier untereinander
    //! wieder dazuaddiert werden. (wurde ja doppelt abgezogen)
    //! Nur der Teil des Overlaps, der auch innerhalb des Tetraeders
    //! liegt darf beruecksichtigt werden!
    //!   => Ebenenwinkel der jeweiligen Kante ins Verhaeltnis zu 2Pi setzen:
    double pi = acos(-1.);
    double gamma = p2p_angle(p2->coord,p1->coord,p3->coord,p0->coord,p1->coord,p3->coord) / pi;
    volume += gamma * get_sphere_overlap(p1->r,p3->r,p1->coord,p3->coord);
    gamma = p2p_angle(p0->coord,p1->coord,p2->coord,p0->coord,p1->coord,p3->coord) / pi;
    volume += gamma * get_sphere_overlap(p1->r,p0->r,p1->coord,p0->coord);
    gamma = p2p_angle(p2->coord,p1->coord,p3->coord,p2->coord,p1->coord,p0->coord) / pi;
    volume += gamma * get_sphere_overlap(p1->r,p2->r,p1->coord,p2->coord);
    gamma = p2p_angle(p2->coord,p0->coord,p3->coord,p2->coord,p0->coord,p1->coord) / pi;
    volume += gamma * get_sphere_overlap(p0->r,p2->r,p0->coord,p2->coord);
    gamma = p2p_angle(p3->coord,p0->coord,p1->coord,p3->coord,p0->coord,p2->coord) / pi;
    volume += gamma * get_sphere_overlap(p0->r,p3->r,p0->coord,p3->coord);
    gamma = p2p_angle(p0->coord,p2->coord,p3->coord,p1->coord,p2->coord,p3->coord) / pi;
    volume += gamma * get_sphere_overlap(p2->r,p3->r,p2->coord,p3->coord);
    //!==============================================================
    
    
    //! Das Volumen kann jetzt negativ sein, weil bei gewichteten Punkten eine Kugel auch
    //! die gegenueberliegende Seite schneiden kann!
    //! Fuer eine exakte Volumenberechnung muss deshalb eine Menge an Tetraedern
    //! betrachett werden (siehe DL_SOLVER.calc_exact_volume):
    //! Das Volumen kann an dieser Stelle aber auch groesser sein, als das des leeren Tetraeders.
    //! Dies ist moeglich, wenn sich 3 Kugeln schneiden (dieses Schnittvolumen muesste eigentlich
    //! wieder abgezogen werden). Haeufig ist es in diesem Fall eine gute Naeherung das Volumen
    //! ebenfalls auf Null zu setzen. In der exakten Berechnung (DT_SOLVER) werden die 3er-Schnitte
    //! mit beruecksichtigt.
    if (volume < 0. || volume > check) volume = 0.;
}

void TETRAHEDRON::calc_exact_volume(double const&salpha) {
    //! Das folgende ist ein bissel Schwachsinn, weil so fuer einzelne Tetraeder keine
    //! exakten Volumen berechnet werden koennen. Das geht nur fuer eine Menge von Tetraedern
    if (sphere_square_radius <= salpha) {
        //! auch bei alpha = 0  ist sichergestellt, dass Tetraeder, bei denen die Kugeln voll
        //! ueberlappen ein Volumen von 0 zugewiesen bekommen, weil ihr sphere_square_radius
        //! einen negativen Wert hat!!!
        volume = 0.;
        return;
    }
    
    double a = get_square_distance(p0->coord,p1->coord);
    if (a < p0->w || a < p1->w) {volume = 0.; return;}
    a = get_square_distance(p0->coord,p2->coord);
    if (a < p0->w || a < p2->w) {volume = 0.; return;}
    a = get_square_distance(p0->coord,p3->coord);
    if (a < p0->w || a < p3->w) {volume = 0.; return;}
    a = get_square_distance(p1->coord,p2->coord);
    if (a < p1->w || a < p2->w) {volume = 0.; return;}
    a = get_square_distance(p1->coord,p3->coord);
    if (a < p1->w || a < p3->w) {volume = 0.; return;}
    a = get_square_distance(p2->coord,p3->coord);
    if (a < p2->w || a < p3->w) {volume = 0.; return;}
    
    volume = fabs(triple_product(p0->coord,p1->coord,p2->coord,p3->coord)) / 6.;
    
    double check1 = volume; /// DEBUG!
    
    //! 1.) Fuer alle 4 Punkte das Sektorvolumen abziehen:
    //!     Laut Edelsbrunner nur, wenn der Punkt zum alpha-Komplex gehoert. Ein Punkt, der frei im Void
    //!     liegt sollte aber auch abgezogen werden!
    //! ACHTUNG: Dringend mal ueberpruefen, ob get_sphere_tetrahedron_intersection() nicht schneller ist!
    volume -= sector_vol(p0->coord,p1->coord,p2->coord,p3->coord,p0->r);
    volume -= sector_vol(p1->coord,p0->coord,p2->coord,p3->coord,p1->r);
    volume -= sector_vol(p2->coord,p1->coord,p0->coord,p3->coord,p2->r);
    volume -= sector_vol(p3->coord,p1->coord,p2->coord,p0->coord,p3->r);
    
    double check2 = volume; /// DEBUG!
    
    //! 2.) Fuer alle Kanten des alpha-Komplex das wedge-Volumen aufaddieren (Schnitt von 2 Kugeln
    //!     und 2 Halbraeumen). ACHTUNG: Das darf wirklich nur fuer Kanten des Shapes gemacht werden,
    //!     deren Tetraeder NICHT zum alpha-Komplex gehoert. Nur so ist sichergestellt, dass der
    //!     zweiseitig begrenzte Kugelschnitt sich komplett im Tetraeder befindet!
    //!     (Wenn die Kante nicht zum alpha-Komplex gehoert, dann duerften auch die Kugeln nicht ueberlappen.)
    
    //! Ich lasse hier den Test weg, ob die entsprechende Kante zum alpha-Shape gehoert, denn wenn nicht, dann
    //! liefert wedge_vol Null, weil es zuerst prueft, ob die Kugeln ueberhaupt ueberlappen:
    volume += wedge_vol(p0->coord,p1->coord,p2->coord,p3->coord,p0->r,p1->r);
    volume += wedge_vol(p0->coord,p2->coord,p1->coord,p3->coord,p0->r,p2->r);
    volume += wedge_vol(p0->coord,p3->coord,p2->coord,p1->coord,p0->r,p3->r);
    volume += wedge_vol(p1->coord,p2->coord,p0->coord,p3->coord,p1->r,p2->r);
    volume += wedge_vol(p1->coord,p3->coord,p2->coord,p0->coord,p1->r,p3->r);
    volume += wedge_vol(p2->coord,p3->coord,p0->coord,p1->coord,p2->r,p3->r);
    
    double check3 = volume; /// DEBUG!
    
    if (check3 < 0.) {
        cerr << "c: " << id << endl;
        cerr << "check1 = " << check1 << "   check2 = " << check2 << "   check3 = " << check3 << endl;
        cerr << "sr = " << sphere_square_radius << "  p0n->sr = " << p0n->sphere_square_radius
             << "  p1n->sr = " << p1n->sphere_square_radius << "  p2n->sr = " << p2n->sphere_square_radius
             << "  p3n->sr = " << p3n->sphere_square_radius << endl;
        visualize(1);
        exit(1);
    }
    
    //! 3.) Fuer alle Facets des alpha-Komplex die pawn-Volumes abziehen (Haelfte des Schnittes
    //!     von 3 Kugeln). Auch hier duerfen wirklich nur Facets, die zum alpha-Shape gehoeren
    //!     betrachtet werden!
    
    if (isnan(volume)) {
        cerr << "a: " << id << endl;
        visualize(1);
        exit(1);
    }
    
    if (p3n->sphere_square_radius <= salpha) volume -= pawn_vol(p0->coord,p1->coord,p2->coord,p0->r,p1->r,p2->r);
    
    if (isnan(volume)) {
        cerr << "b: " << id << endl;
        visualize(1);
        exit(1);
    }
    
    if (p2n->sphere_square_radius <= salpha) volume -= pawn_vol(p0->coord,p1->coord,p3->coord,p0->r,p1->r,p3->r);
    
    if (isnan(volume)) {
        
        double os0 = get_sphere_tetrahedron_intersection(p0->coord,p0->r,p1->coord,p2->coord,p3->coord);
        double os1 = get_sphere_tetrahedron_intersection(p1->coord,p1->r,p0->coord,p2->coord,p3->coord);
        double os2 = get_sphere_tetrahedron_intersection(p2->coord,p2->r,p0->coord,p1->coord,p3->coord);
        double os3 = get_sphere_tetrahedron_intersection(p3->coord,p3->r,p0->coord,p1->coord,p2->coord);
    
        cerr << "c: " << id << endl;
        cerr << "check1 = " << check1 << "   check2 = " << check2 << "   check3 = " << check3 << endl;
        cerr << "s0 = " << sector_vol(p0->coord,p1->coord,p2->coord,p3->coord,p0->r) << "   os0 = " << os0 << endl;
        cerr << "s1 = " << sector_vol(p1->coord,p0->coord,p2->coord,p3->coord,p1->r) << "   os1 = " << os1 << endl;
        cerr << "s2 = " << sector_vol(p2->coord,p1->coord,p0->coord,p3->coord,p2->r) << "   os2 = " << os2 << endl;
        cerr << "s3 = " << sector_vol(p3->coord,p1->coord,p2->coord,p0->coord,p3->r) << "   os3 = " << os3 << endl;
        
        cerr << angle_solid(p0->coord,p1->coord,p2->coord,p3->coord) << endl;
        cerr << angle_solid(p1->coord,p0->coord,p2->coord,p3->coord) << endl;
        cerr << angle_solid(p2->coord,p1->coord,p0->coord,p3->coord) << endl;
        cerr << angle_solid(p3->coord,p1->coord,p2->coord,p0->coord) << endl;
        
        cerr << "rad = " << sqrt(sphere_square_radius) << endl;
        ofstream f_out;
        
        f_out.open("alpha_sphere.py");
        f_out << "aspheres = [";
        f_out << "7.0," << sphere_middle[0] << "," << sphere_middle[1] << "," << sphere_middle[2] << "," << sqrt(sphere_square_radius);
        f_out << "]" << endl;
        f_out << "cmd.load_cgo(aspheres,'aspheres',1)" << endl;
        f_out.close();
        
        visualize(1);
        exit(1);
    }
    
    if (p1n->sphere_square_radius <= salpha) volume -= pawn_vol(p0->coord,p2->coord,p3->coord,p0->r,p2->r,p3->r);
    
    if (isnan(volume)) {
        cerr << "d: " << id << endl;
        visualize(1);
        exit(1);
    }
    
    if (p0n->sphere_square_radius <= salpha) volume -= pawn_vol(p1->coord,p2->coord,p3->coord,p1->r,p2->r,p3->r);
    
    if (volume < 0. || volume > check1) {
        //! Es kommt leider vor, dass die Punkte fast koplanar sind und dann wird ein
        //! sehr hohes Volumen ausgerechnet
        
    //    cerr << "check1 = " << check1 << "   check2 = " << check2 << "   check3 = " << check3 << "   vol = " << volume << endl;
    //    if (volume > check1) {
        //    cerr << volume << " > " << check1 << endl;
        //    cerr << "check1 = " << check1 << "   check2 = " << check2 << "   check3 = " << check3 << endl;
        //    visualize(1);
        //    exit(1);
    //    }
        
        volume = 0.;
    }
    
    if (isnan(volume)) {
        cerr << "e: " << id << endl;
        visualize(1);
        exit(1);
    }
    
}

void TETRAHEDRON::reduce_sphere() {
    if (sphere_square_radius < 0.) { //!Kann bei gewichteten Punkten negativ sein!!! (siehe linalg_GN.hpp)
        sphere_square_radius = 0.;
        return;
    }
    double hlp = sqrt(sphere_square_radius);
    double od0 = sqrt(p0->w) + hlp - get_distance(sphere_middle,p0->coord);
    double od1 = sqrt(p1->w) + hlp - get_distance(sphere_middle,p1->coord);
    double od2 = sqrt(p2->w) + hlp - get_distance(sphere_middle,p2->coord);
    double od3 = sqrt(p3->w) + hlp - get_distance(sphere_middle,p3->coord);
    double hlp2 = od0;
    if (od1 > hlp2) hlp2 = od1;
    if (od2 > hlp2) hlp2 = od2;
    if (od3 > hlp2) hlp2 = od3;
    sphere_square_radius = (hlp - hlp2) * (hlp - hlp2);
}

bool TETRAHEDRON::is_positive_oriented() {
    if (triple_product(p0->coord,p1->coord,p2->coord,p3->coord) > 0.) return true;
    return false;
}

bool TETRAHEDRON::contains_p(stl_ptr<DT_POINT> const& point) {
    //!Wenn die Determinante Null ist liegt der neue Punkt exakt auf einer Tetraederflaeche.
    //!In diesem Fall wird er willkuerlich dem aktuellen Tetraeder zugeordnet, was zu einem
    //!unechten (4 koplanare Punkte) Tetraeder fuehrt. Die folgenden Flips sollten das Problem
    //!wieder beheben:
    if (triple_product(p0->coord,p1->coord,p2->coord,point->coord) < 0.) return false;
    if (triple_product(p0->coord,p3->coord,p1->coord,point->coord) < 0.) return false;
    if (triple_product(p0->coord,p2->coord,p3->coord,point->coord) < 0.) return false;
    if (triple_product(p1->coord,p3->coord,p2->coord,point->coord) < 0.) return false;
    
    return true;
}

void TETRAHEDRON::check_link(stl_ptr<TETRAHEDRON> const& t) {
    //! So ist die Funktion zwar sehr haesslich und lang, dafuer aber etwa 10mal schneller als
    //! die auskommentierte Variante mit den Sets!!!
    if (!p0.equal_addr(t->p0) && !p0.equal_addr(t->p1) && !p0.equal_addr(t->p2) && !p0.equal_addr(t->p3)) {
        if (p1.equal_addr(t->p0)) {
            if (p2.equal_addr(t->p1)) {
                if (p3.equal_addr(t->p2)) {
                    //! 123 -> 012
                    p0p = t->p3; t->p3n = this; t->p3p = p0;
                } else if (p3.equal_addr(t->p3)) {
                    //! 123 -> 013
                    p0p = t->p2; t->p2n = this; t->p2p = p0;
                }
            } else if (p2.equal_addr(t->p2)) {
                if (p3.equal_addr(t->p1)) {
                    //! 123 -> 021
                    p0p = t->p3; t->p3n = this; t->p3p = p0;
                } else if (p3.equal_addr(t->p3)) {
                    //! 123 -> 023
                    p0p = t->p1; t->p1n = this; t->p1p = p0;
                }
            } else if (p2.equal_addr(t->p3)) {
                if (p3.equal_addr(t->p1)) {
                    //! 123 -> 031
                    p0p = t->p2; t->p2n = this; t->p2p = p0;
                } else if (p3.equal_addr(t->p2)) {
                    //! 123 -> 032
                    p0p = t->p1; t->p1n = this; t->p1p = p0;
                }
            }
        } else if (p1.equal_addr(t->p1)) {
            if (p2.equal_addr(t->p0)) {
                if (p3.equal_addr(t->p2)) {
                    //! 123 -> 102
                    p0p = t->p3; t->p3n = this; t->p3p = p0;
                } else if (p3.equal_addr(t->p3)) {
                    //! 123 -> 103
                    p0p = t->p2; t->p2n = this; t->p2p = p0;
                }
            } else if (p2.equal_addr(t->p2)) {
                if (p3.equal_addr(t->p0)) {
                    //! 123 -> 120
                    p0p = t->p3; t->p3n = this; t->p3p = p0;
                } else if (p3.equal_addr(t->p3)) {
                    //! 123 -> 123
                    p0p = t->p0; t->p0n = this; t->p0p = p0;
                }
            } else if (p2.equal_addr(t->p3)) {
                if (p3.equal_addr(t->p0)) {
                    //! 123 -> 130
                    p0p = t->p2; t->p2n = this; t->p2p = p0;
                } else if (p3.equal_addr(t->p2)) {
                    //! 123 -> 132
                    p0p = t->p0; t->p0n = this; t->p0p = p0;
                }
            }
        } else if (p1.equal_addr(t->p2)) {
            if (p2.equal_addr(t->p1)) {
                if (p3.equal_addr(t->p0)) {
                    //! 123 -> 210
                    p0p = t->p3; t->p3n = this; t->p3p = p0;
                } else if (p3.equal_addr(t->p3)) {
                    //! 123 -> 213
                    p0p = t->p0; t->p0n = this; t->p0p = p0;
                }
            } else if (p2.equal_addr(t->p0)) {
                if (p3.equal_addr(t->p1)) {
                    //! 123 -> 201
                    p0p = t->p3; t->p3n = this; t->p3p = p0;
                } else if (p3.equal_addr(t->p3)) {
                    //! 123 -> 203
                    p0p = t->p1; t->p1n = this; t->p1p = p0;
                }
            } else if (p2.equal_addr(t->p3)) {
                if (p3.equal_addr(t->p1)) {
                    //! 123 -> 231
                    p0p = t->p0; t->p0n = this; t->p0p = p0;
                } else if (p3.equal_addr(t->p0)) {
                    //! 123 -> 230
                    p0p = t->p1; t->p1n = this; t->p1p = p0;
                }
            }
        } else if (p1.equal_addr(t->p3)) {
            if (p2.equal_addr(t->p1)) {
                if (p3.equal_addr(t->p2)) {
                    //! 123 -> 312
                    p0p = t->p0; t->p0n = this; t->p0p = p0;
                } else if (p3.equal_addr(t->p0)) {
                    //! 123 -> 310
                    p0p = t->p2; t->p2n = this; t->p2p = p0;
                }
            } else if (p2.equal_addr(t->p2)) {
                if (p3.equal_addr(t->p1)) {
                    //! 123 -> 321
                    p0p = t->p0; t->p0n = this; t->p0p = p0;
                } else if (p3.equal_addr(t->p0)) {
                    //! 123 -> 320
                    p0p = t->p1; t->p1n = this; t->p1p = p0;
                }
            } else if (p2.equal_addr(t->p0)) {
                if (p3.equal_addr(t->p1)) {
                    //! 123 -> 301
                    p0p = t->p2; t->p2n = this; t->p2p = p0;
                } else if (p3.equal_addr(t->p2)) {
                    //! 123 -> 302
                    p0p = t->p1; t->p1n = this; t->p1p = p0;
                }
            }
        }
        //!###########################################################################################################
    } else if (!p1.equal_addr(t->p0) && !p1.equal_addr(t->p1) && !p1.equal_addr(t->p2) && !p1.equal_addr(t->p3)) {
        if (p0.equal_addr(t->p0)) {
            if (p2.equal_addr(t->p1)) {
                if (p3.equal_addr(t->p2)) {
                    //! 023 -> 012
                    p1p = t->p3; t->p3n = this; t->p3p = p1;
                } else if (p3.equal_addr(t->p3)) {
                    //! 023 -> 013
                    p1p = t->p2; t->p2n = this; t->p2p = p1;
                }
            } else if (p2.equal_addr(t->p2)) {
                if (p3.equal_addr(t->p1)) {
                    //! 023 -> 021
                    p1p = t->p3; t->p3n = this; t->p3p = p1;
                } else if (p3.equal_addr(t->p3)) {
                    //! 023 -> 023
                    p1p = t->p1; t->p1n = this; t->p1p = p1;
                }
            } else if (p2.equal_addr(t->p3)) {
                if (p3.equal_addr(t->p1)) {
                    //! 023 -> 031
                    p1p = t->p2; t->p2n = this; t->p2p = p1;
                } else if (p3.equal_addr(t->p2)) {
                    //! 023 -> 032
                    p1p = t->p1; t->p1n = this; t->p1p = p1;
                }
            }
        } else if (p0.equal_addr(t->p1)) {
            if (p2.equal_addr(t->p0)) {
                if (p3.equal_addr(t->p2)) {
                    //! 023 -> 102
                    p1p = t->p3; t->p3n = this; t->p3p = p1;
                } else if (p3.equal_addr(t->p3)) {
                    //! 023 -> 103
                    p1p = t->p2; t->p2n = this; t->p2p = p1;
                }
            } else if (p2.equal_addr(t->p2)) {
                if (p3.equal_addr(t->p0)) {
                    //! 023 -> 120
                    p1p = t->p3; t->p3n = this; t->p3p = p1;
                } else if (p3.equal_addr(t->p3)) {
                    //! 023 -> 123
                    p1p = t->p0; t->p0n = this; t->p0p = p1;
                }
            } else if (p2.equal_addr(t->p3)) {
                if (p3.equal_addr(t->p0)) {
                    //! 023 -> 130
                    p1p = t->p2; t->p2n = this; t->p2p = p1;
                } else if (p3.equal_addr(t->p2)) {
                    //! 023 -> 132
                    p1p = t->p0; t->p0n = this; t->p0p = p1;
                }
            }
        } else if (p0.equal_addr(t->p2)) {
            if (p2.equal_addr(t->p1)) {
                if (p3.equal_addr(t->p0)) {
                    //! 023 -> 210
                    p1p = t->p3; t->p3n = this; t->p3p = p1;
                } else if (p3.equal_addr(t->p3)) {
                    //! 023 -> 213
                    p1p = t->p0; t->p0n = this; t->p0p = p1;
                }
            } else if (p2.equal_addr(t->p0)) {
                if (p3.equal_addr(t->p1)) {
                    //! 023 -> 201
                    p1p = t->p3; t->p3n = this; t->p3p = p1;
                } else if (p3.equal_addr(t->p3)) {
                    //! 023 -> 203
                    p1p = t->p1; t->p1n = this; t->p1p = p1;
                }
            } else if (p2.equal_addr(t->p3)) {
                if (p3.equal_addr(t->p1)) {
                    //! 023 -> 231
                    p1p = t->p0; t->p0n = this; t->p0p = p1;
                } else if (p3.equal_addr(t->p0)) {
                    //! 023 -> 230
                    p1p = t->p1; t->p1n = this; t->p1p = p1;
                }
            }
        } else if (p0.equal_addr(t->p3)) {
            if (p2.equal_addr(t->p1)) {
                if (p3.equal_addr(t->p2)) {
                    //! 023 -> 312
                    p1p = t->p0; t->p0n = this; t->p0p = p1;
                } else if (p3.equal_addr(t->p0)) {
                    //! 023 -> 310
                    p1p = t->p2; t->p2n = this; t->p2p = p1;
                }
            } else if (p2.equal_addr(t->p2)) {
                if (p3.equal_addr(t->p1)) {
                    //! 023 -> 321
                    p1p = t->p0; t->p0n = this; t->p0p = p1;
                } else if (p3.equal_addr(t->p0)) {
                    //! 023 -> 320
                    p1p = t->p1; t->p1n = this; t->p1p = p1;
                }
            } else if (p2.equal_addr(t->p0)) {
                if (p3.equal_addr(t->p1)) {
                    //! 023 -> 301
                    p1p = t->p2; t->p2n = this; t->p2p = p1;
                } else if (p3.equal_addr(t->p2)) {
                    //! 023 -> 302
                    p1p = t->p1; t->p1n = this; t->p1p = p1;
                }
            }
        }
        //!###########################################################################################################
    } else if (!p2.equal_addr(t->p0) && !p2.equal_addr(t->p1) && !p2.equal_addr(t->p2) && !p2.equal_addr(t->p3)) {
        if (p0.equal_addr(t->p0)) {
            if (p1.equal_addr(t->p1)) {
                if (p3.equal_addr(t->p2)) {
                    //! 013 -> 012
                    p2p = t->p3; t->p3n = this; t->p3p = p2;
                } else if (p3.equal_addr(t->p3)) {
                    //! 013 -> 013
                    p2p = t->p2; t->p2n = this; t->p2p = p2;
                }
            } else if (p1.equal_addr(t->p2)) {
                if (p3.equal_addr(t->p1)) {
                    //! 013 -> 021
                    p2p = t->p3; t->p3n = this; t->p3p = p2;
                } else if (p3.equal_addr(t->p3)) {
                    //! 013 -> 023
                    p2p = t->p1; t->p1n = this; t->p1p = p2;
                }
            } else if (p1.equal_addr(t->p3)) {
                if (p3.equal_addr(t->p1)) {
                    //! 013 -> 031
                    p2p = t->p2; t->p2n = this; t->p2p = p2;
                } else if (p3.equal_addr(t->p2)) {
                    //! 013 -> 032
                    p2p = t->p1; t->p1n = this; t->p1p = p2;
                }
            }
        } else if (p0.equal_addr(t->p1)) {
            if (p1.equal_addr(t->p0)) {
                if (p3.equal_addr(t->p2)) {
                    //! 013 -> 102
                    p2p = t->p3; t->p3n = this; t->p3p = p2;
                } else if (p3.equal_addr(t->p3)) {
                    //! 013 -> 103
                    p2p = t->p2; t->p2n = this; t->p2p = p2;
                }
            } else if (p1.equal_addr(t->p2)) {
                if (p3.equal_addr(t->p0)) {
                    //! 013 -> 120
                    p2p = t->p3; t->p3n = this; t->p3p = p2;
                } else if (p3.equal_addr(t->p3)) {
                    //! 013 -> 123
                    p2p = t->p0; t->p0n = this; t->p0p = p2;
                }
            } else if (p1.equal_addr(t->p3)) {
                if (p3.equal_addr(t->p0)) {
                    //! 013 -> 130
                    p2p = t->p2; t->p2n = this; t->p2p = p2;
                } else if (p3.equal_addr(t->p2)) {
                    //! 013 -> 132
                    p2p = t->p0; t->p0n = this; t->p0p = p2;
                }
            }
        } else if (p0.equal_addr(t->p2)) {
            if (p1.equal_addr(t->p1)) {
                if (p3.equal_addr(t->p0)) {
                    //! 013 -> 210
                    p2p = t->p3; t->p3n = this; t->p3p = p2;
                } else if (p3.equal_addr(t->p3)) {
                    //! 013 -> 213
                    p2p = t->p0; t->p0n = this; t->p0p = p2;
                }
            } else if (p1.equal_addr(t->p0)) {
                if (p3.equal_addr(t->p1)) {
                    //! 013 -> 201
                    p2p = t->p3; t->p3n = this; t->p3p = p2;
                } else if (p3.equal_addr(t->p3)) {
                    //! 013 -> 203
                    p2p = t->p1; t->p1n = this; t->p1p = p2;
                }
            } else if (p1.equal_addr(t->p3)) {
                if (p3.equal_addr(t->p1)) {
                    //! 013 -> 231
                    p2p = t->p0; t->p0n = this; t->p0p = p2;
                } else if (p3.equal_addr(t->p0)) {
                    //! 013 -> 230
                    p2p = t->p1; t->p1n = this; t->p1p = p2;
                }
            }
        } else if (p0.equal_addr(t->p3)) {
            if (p1.equal_addr(t->p1)) {
                if (p3.equal_addr(t->p2)) {
                    //! 013 -> 312
                    p2p = t->p0; t->p0n = this; t->p0p = p2;
                } else if (p3.equal_addr(t->p0)) {
                    //! 013 -> 310
                    p2p = t->p2; t->p2n = this; t->p2p = p2;
                }
            } else if (p1.equal_addr(t->p2)) {
                if (p3.equal_addr(t->p1)) {
                    //! 013 -> 321
                    p2p = t->p0; t->p0n = this; t->p0p = p2;
                } else if (p3.equal_addr(t->p0)) {
                    //! 013 -> 320
                    p2p = t->p1; t->p1n = this; t->p1p = p2;
                }
            } else if (p1.equal_addr(t->p0)) {
                if (p3.equal_addr(t->p1)) {
                    //! 013 -> 301
                    p2p = t->p2; t->p2n = this; t->p2p = p2;
                } else if (p3.equal_addr(t->p2)) {
                    //! 013 -> 302
                    p2p = t->p1; t->p1n = this; t->p1p = p2;
                }
            }
        }
        //!###########################################################################################################
    } else if (!p3.equal_addr(t->p0) && !p3.equal_addr(t->p1) && !p3.equal_addr(t->p2) && !p3.equal_addr(t->p3)) {
        if (p0.equal_addr(t->p0)) {
            if (p1.equal_addr(t->p1)) {
                if (p2.equal_addr(t->p2)) {
                    //! 012 -> 012
                    p3p = t->p3; t->p3n = this; t->p3p = p3;
                } else if (p2.equal_addr(t->p3)) {
                    //! 012 -> 013
                    p3p = t->p2; t->p2n = this; t->p2p = p3;
                }
            } else if (p1.equal_addr(t->p2)) {
                if (p2.equal_addr(t->p1)) {
                    //! 012 -> 021
                    p3p = t->p3; t->p3n = this; t->p3p = p3;
                } else if (p2.equal_addr(t->p3)) {
                    //! 012 -> 023
                    p3p = t->p1; t->p1n = this; t->p1p = p3;
                }
            } else if (p1.equal_addr(t->p3)) {
                if (p2.equal_addr(t->p1)) {
                    //! 012 -> 031
                    p3p = t->p2; t->p2n = this; t->p2p = p3;
                } else if (p2.equal_addr(t->p2)) {
                    //! 012 -> 032
                    p3p = t->p1; t->p1n = this; t->p1p = p3;
                }
            }
        } else if (p0.equal_addr(t->p1)) {
            if (p1.equal_addr(t->p0)) {
                if (p2.equal_addr(t->p2)) {
                    //! 012 -> 102
                    p3p = t->p3; t->p3n = this; t->p3p = p3;
                } else if (p2.equal_addr(t->p3)) {
                    //! 012 -> 103
                    p3p = t->p2; t->p2n = this; t->p2p = p3;
                }
            } else if (p1.equal_addr(t->p2)) {
                if (p2.equal_addr(t->p0)) {
                    //! 012 -> 120
                    p3p = t->p3; t->p3n = this; t->p3p = p3;
                } else if (p2.equal_addr(t->p3)) {
                    //! 012 -> 123
                    p3p = t->p0; t->p0n = this; t->p0p = p3;
                }
            } else if (p1.equal_addr(t->p3)) {
                if (p2.equal_addr(t->p0)) {
                    //! 012 -> 130
                    p3p = t->p2; t->p2n = this; t->p2p = p3;
                } else if (p2.equal_addr(t->p2)) {
                    //! 012 -> 132
                    p3p = t->p0; t->p0n = this; t->p0p = p3;
                }
            }
        } else if (p0.equal_addr(t->p2)) {
            if (p1.equal_addr(t->p1)) {
                if (p2.equal_addr(t->p0)) {
                    //! 012 -> 210
                    p3p = t->p3; t->p3n = this; t->p3p = p3;
                } else if (p2.equal_addr(t->p3)) {
                    //! 012 -> 213
                    p3p = t->p0; t->p0n = this; t->p0p = p3;
                }
            } else if (p1.equal_addr(t->p0)) {
                if (p2.equal_addr(t->p1)) {
                    //! 012 -> 201
                    p3p = t->p3; t->p3n = this; t->p3p = p3;
                } else if (p2.equal_addr(t->p3)) {
                    //! 012 -> 203
                    p3p = t->p1; t->p1n = this; t->p1p = p3;
                }
            } else if (p1.equal_addr(t->p3)) {
                if (p2.equal_addr(t->p1)) {
                    //! 012 -> 231
                    p3p = t->p0; t->p0n = this; t->p0p = p3;
                } else if (p2.equal_addr(t->p0)) {
                    //! 012 -> 230
                    p3p = t->p1; t->p1n = this; t->p1p = p3;
                }
            }
        } else if (p0.equal_addr(t->p3)) {
            if (p1.equal_addr(t->p1)) {
                if (p2.equal_addr(t->p2)) {
                    //! 012 -> 312
                    p3p = t->p0; t->p0n = this; t->p0p = p3;
                } else if (p2.equal_addr(t->p0)) {
                    //! 012 -> 310
                    p3p = t->p2; t->p2n = this; t->p2p = p3;
                }
            } else if (p1.equal_addr(t->p2)) {
                if (p2.equal_addr(t->p1)) {
                    //! 012 -> 321
                    p3p = t->p0; t->p0n = this; t->p0p = p3;
                } else if (p2.equal_addr(t->p0)) {
                    //! 012 -> 320
                    p3p = t->p1; t->p1n = this; t->p1p = p3;
                }
            } else if (p1.equal_addr(t->p0)) {
                if (p2.equal_addr(t->p1)) {
                    //! 012 -> 301
                    p3p = t->p2; t->p2n = this; t->p2p = p3;
                } else if (p2.equal_addr(t->p2)) {
                    //! 012 -> 302
                    p3p = t->p1; t->p1n = this; t->p1p = p3;
                }
            }
        }
        //!###########################################################################################################
    }
}
/*
void TETRAHEDRON::check_link(stl_ptr<TETRAHEDRON> const& t) {
    tr1::unordered_set<int> t1,t2;
    
    if (p0.equal_addr(t->p0)) {t1.insert(0); t2.insert(0);}
    else if (p0.equal_addr(t->p1)) {t1.insert(0); t2.insert(1);}
    else if (p0.equal_addr(t->p2)) {t1.insert(0); t2.insert(2);}
    else if (p0.equal_addr(t->p3)) {t1.insert(0); t2.insert(3);}
    
    if (p1.equal_addr(t->p0)) {t1.insert(1); t2.insert(0);}
    else if (p1.equal_addr(t->p1)) {t1.insert(1); t2.insert(1);}
    else if (p1.equal_addr(t->p2)) {t1.insert(1); t2.insert(2);}
    else if (p1.equal_addr(t->p3)) {t1.insert(1); t2.insert(3);}
    
    if (t1.size() == 0) return;
    
    if (p2.equal_addr(t->p0)) {t1.insert(2); t2.insert(0);}
    else if (p2.equal_addr(t->p1)) {t1.insert(2); t2.insert(1);}
    else if (p2.equal_addr(t->p2)) {t1.insert(2); t2.insert(2);}
    else if (p2.equal_addr(t->p3)) {t1.insert(2); t2.insert(3);}
    
    if (t1.size() < 2) return;
    
    if (p3.equal_addr(t->p0)) {t1.insert(3); t2.insert(0);}
    else if (p3.equal_addr(t->p1)) {t1.insert(3); t2.insert(1);}
    else if (p3.equal_addr(t->p2)) {t1.insert(3); t2.insert(2);}
    else if (p3.equal_addr(t->p3)) {t1.insert(3); t2.insert(3);}
    
    if (t1.size() == 3) {
        if (t1.find(0) == t1.end()) {
            p0n = t;
            if (t2.find(0) == t2.end()) {p0p = t->p0; t->p0n = this; t->p0p = p0;}
            else if (t2.find(1) == t2.end()) {p0p = t->p1; t->p1n = this; t->p1p = p0;}
            else if (t2.find(2) == t2.end()) {p0p = t->p2; t->p2n = this; t->p2p = p0;}
            else if (t2.find(3) == t2.end()) {p0p = t->p3; t->p3n = this; t->p3p = p0;}
        } else if (t1.find(1) == t1.end()) {
            p1n = t;
            if (t2.find(0) == t2.end()) {p1p = t->p0; t->p0n = this; t->p0p = p1;}
            else if (t2.find(1) == t2.end()) {p1p = t->p1; t->p1n = this; t->p1p = p1;}
            else if (t2.find(2) == t2.end()) {p1p = t->p2; t->p2n = this; t->p2p = p1;}
            else if (t2.find(3) == t2.end()) {p1p = t->p3; t->p3n = this; t->p3p = p1;}
        } else if (t1.find(2) == t1.end()) {
            p2n = t;
            if (t2.find(0) == t2.end()) {p2p = t->p0; t->p0n = this; t->p0p = p2;}
            else if (t2.find(1) == t2.end()) {p2p = t->p1; t->p1n = this; t->p1p = p2;}
            else if (t2.find(2) == t2.end()) {p2p = t->p2; t->p2n = this; t->p2p = p2;}
            else if (t2.find(3) == t2.end()) {p2p = t->p3; t->p3n = this; t->p3p = p2;}
        } else if (t1.find(3) == t1.end()) {
            p3n = t;
            if (t2.find(0) == t2.end()) {p3p = t->p0; t->p0n = this; t->p0p = p3;}
            else if (t2.find(1) == t2.end()) {p3p = t->p1; t->p1n = this; t->p1p = p3;}
            else if (t2.find(2) == t2.end()) {p3p = t->p2; t->p2n = this; t->p2p = p3;}
            else if (t2.find(3) == t2.end()) {p3p = t->p3; t->p3n = this; t->p3p = p3;}
        }
    }
}
*/

void TETRAHEDRON::relink_neighbours() {
    //!Diese Methode ist recht umstaendlich! Besser ist es die Tetraeder direkt dort umzulinken, wo sie erzeugt werden.
    //!Dies hier ist nur zu Debug-Zwecken, um festzustellen, ob das "effiziente" Umlinken an Ort und Stelle
    //!einen Fehler hat.
    if (!p0n.zero()) check_link(p0n);
    if (!p1n.zero()) check_link(p1n);
    if (!p2n.zero()) check_link(p2n);
    if (!p3n.zero()) check_link(p3n);
}

void TETRAHEDRON::generate_facets(stack<stl_ptr<L_FACET> > &facets,int64_t curr_id) {
    stl_ptr<TETRAHEDRON> tz(this);
    if (!p0n.zero()) {
        if (generation != p0n->generation) {
            facets.push(new L_FACET(p1,p2,p3,
                                        tz,p0n,
                                            p0,p0p,curr_id));
            ++curr_id;
        }
    }
    
    if (!p1n.zero()) {
        if (generation != p1n->generation) {
            facets.push(new L_FACET(p2,p0,p3,
                                        tz,p1n,
                                            p1,p1p,curr_id));
            ++curr_id;
        }
    }
    
    if (!p2n.zero()) {
        if (generation != p2n->generation) {
            facets.push(new L_FACET(p0,p1,p3,
                                        tz,p2n,
                                            p2,p2p,curr_id));
            ++curr_id;
        }
    }
    
    if (!p3n.zero()) {
        if (generation != p3n->generation) {
            facets.push(new L_FACET(p0,p2,p1,
                                        tz,p3n,
                                            p3,p3p,curr_id));
            ++curr_id;
        }
    }
}
//!===============================================================================


//!===============================================================================
//!Definition der Klasse L_FACET:
L_FACET::L_FACET() {
    id = -1;
}

L_FACET::L_FACET(stl_ptr<DT_POINT> const& vp0,stl_ptr<DT_POINT> const& vp1,
                 stl_ptr<DT_POINT> const& vp2) : p0(vp0),p1(vp1),p2(vp2) {}

L_FACET::L_FACET(stl_ptr<DT_POINT> const& vp0,stl_ptr<DT_POINT> const& vp1,
                 stl_ptr<DT_POINT> const& vp2,stl_ptr<TETRAHEDRON> const& vt1,
                 stl_ptr<TETRAHEDRON> const& vt2,stl_ptr<DT_POINT> const& vps1,
                 stl_ptr<DT_POINT> const& vps2,int64_t const& i) : p0(vp0), p1(vp1),
                 p2(vp2), t1(vt1), t2(vt2), ps1(vps1), ps2(vps2), id(i) {}

L_FACET::~L_FACET() {}

L_FACET& L_FACET::operator=(L_FACET const& facet) {
    p0 = facet.p0;
    p1 = facet.p1;
    p2 = facet.p2;
    t1 = facet.t1;
    t2 = facet.t2;
    ps1 = facet.ps1;
    ps2 = facet.ps2;
    id = facet.id;
    return *this;
}

bool L_FACET::is_regular() {
    if (t1->next.size() > 0 || t2->next.size() > 0) return true;
    
    if ((get_square_distance(ps1->coord,t2->sphere_middle)-ps1->w) < t2->sphere_square_radius ||
        (get_square_distance(ps2->coord,t1->sphere_middle)-ps2->w) < t1->sphere_square_radius) return false;
    return true;
}

void L_FACET::calc_n() {
    if (t1->p0.equal_addr(p0)) {
        t1n0 = t1->p0n;
        t1p0 = t1->p0p;
    } else if (t1->p1.equal_addr(p0)) {
        t1n0 = t1->p1n;
        t1p0 = t1->p1p;
    } else if (t1->p2.equal_addr(p0)) {
        t1n0 = t1->p2n;
        t1p0 = t1->p2p;
    } else {
        t1n0 = t1->p3n;
        t1p0 = t1->p3p;
    }
    if (t1->p0.equal_addr(p1)) {
        t1n1 = t1->p0n;
        t1p1 = t1->p0p;
    } else if (t1->p1.equal_addr(p1)) {
        t1n1 = t1->p1n;
        t1p1 = t1->p1p;
    } else if (t1->p2.equal_addr(p1)) {
        t1n1 = t1->p2n;
        t1p1 = t1->p2p;
    } else {
        t1n1 = t1->p3n;
        t1p1 = t1->p3p;
    }
    if (t1->p0.equal_addr(p2)) {
        t1n2 = t1->p0n;
        t1p2 = t1->p0p;
    } else if (t1->p1.equal_addr(p2)) {
        t1n2 = t1->p1n;
        t1p2 = t1->p1p;
    } else if (t1->p2.equal_addr(p2)) {
        t1n2 = t1->p2n;
        t1p2 = t1->p2p;
    } else {
        t1n2 = t1->p3n;
        t1p2 = t1->p3p;
    }
    
    if (t2->p0.equal_addr(p0)) {
        t2n0 = t2->p0n;
        t2p0 = t2->p0p;
    } else if (t2->p1.equal_addr(p0)) {
        t2n0 = t2->p1n;
        t2p0 = t2->p1p;
    } else if (t2->p2.equal_addr(p0)) {
        t2n0 = t2->p2n;
        t2p0 = t2->p2p;
    } else {
        t2n0 = t2->p3n;
        t2p0 = t2->p3p;
    }
    if (t2->p0.equal_addr(p1)) {
        t2n1 = t2->p0n;
        t2p1 = t2->p0p;
    } else if (t2->p1.equal_addr(p1)) {
        t2n1 = t2->p1n;
        t2p1 = t2->p1p;
    } else if (t2->p2.equal_addr(p1)) {
        t2n1 = t2->p2n;
        t2p1 = t2->p2p;
    } else {
        t2n1 = t2->p3n;
        t2p1 = t2->p3p;
    }
    if (t2->p0.equal_addr(p2)) {
        t2n2 = t2->p0n;
        t2p2 = t2->p0p;
    } else if (t2->p1.equal_addr(p2)) {
        t2n2 = t2->p1n;
        t2p2 = t2->p1p;
    } else if (t2->p2.equal_addr(p2)) {
        t2n2 = t2->p2n;
        t2p2 = t2->p2p;
    } else {
        t2n2 = t2->p3n;
        t2p2 = t2->p3p;
    }
}

void L_FACET::visualize() {
    ofstream f_out;
    f_out.open("facet.py");
    
    float rad = 0.1;
    
    f_out << "p0=[7.0,"
        << p0->coord[0] << "," << p0->coord[1] << "," << p0->coord[2] << ",0.3]" << endl;
    f_out << "cmd.load_cgo(p0,'p0',1)" << endl;
    
    f_out << "p1=[7.0,"
        << p1->coord[0] << "," << p1->coord[1] << "," << p1->coord[2] << ",0.3]" << endl;
    f_out << "cmd.load_cgo(p1,'p1',1)" << endl;
    
    f_out << "p2=[7.0,"
        << p2->coord[0] << "," << p2->coord[1] << "," << p2->coord[2] << ",0.3]" << endl;
    f_out << "cmd.load_cgo(p2,'p2',1)" << endl;
    /*
    f_out << "ps1=[7.0,"
        << ps1->coord[0] << "," << ps1->coord[1] << "," << ps1->coord[2] << ",0.3]" << endl;
    f_out << "cmd.load_cgo(ps1,'ps1',1)" << endl;
    
    f_out << "ps2=[7.0,"
        << ps2->coord[0] << "," << ps2->coord[1] << "," << ps2->coord[2] << ",0.3]" << endl;
    f_out << "cmd.load_cgo(ps2,'ps2',1)" << endl;
    */
    f_out << "facet=[9.0,"
          << p0->coord[0] << "," << p0->coord[1] << "," << p0->coord[2] << ","
          << p1->coord[0] << "," << p1->coord[1] << "," << p1->coord[2] << ","
          << rad << ",1.0,0.0,0.0,1.0,0.0,0.0,9.0,"
          
          << p0->coord[0] << "," << p0->coord[1] << "," << p0->coord[2] << ","
          << p2->coord[0] << "," << p2->coord[1] << "," << p2->coord[2] << ","
          << rad << ",1.0,0.0,0.0,1.0,0.0,0.0,9.0,"
    
          << p2->coord[0] << "," << p2->coord[1] << "," << p2->coord[2] << ","
          << p1->coord[0] << "," << p1->coord[1] << "," << p1->coord[2] << ","
          << rad << ",1.0,0.0,0.0,1.0,0.0,0.0]" << endl;
    
    f_out << "cmd.load_cgo(facet,'facet',1)" << endl;
    /*
    f_out << "t1=[9.0,"
              << t1->p0->coord[0] << "," << t1->p0->coord[1] << "," << t1->p0->coord[2] << ","
              << t1->p1->coord[0] << "," << t1->p1->coord[1] << "," << t1->p1->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0,9.0,"
              
              << t1->p0->coord[0] << "," << t1->p0->coord[1] << "," << t1->p0->coord[2] << ","
              << t1->p2->coord[0] << "," << t1->p2->coord[1] << "," << t1->p2->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0,9.0,"
              
              << t1->p0->coord[0] << "," << t1->p0->coord[1] << "," << t1->p0->coord[2] << ","
              << t1->p3->coord[0] << "," << t1->p3->coord[1] << "," << t1->p3->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0,9.0,"
              
              << t1->p1->coord[0] << "," << t1->p1->coord[1] << "," << t1->p1->coord[2] << ","
              << t1->p2->coord[0] << "," << t1->p2->coord[1] << "," << t1->p2->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0,9.0,"
              
              << t1->p2->coord[0] << "," << t1->p2->coord[1] << "," << t1->p2->coord[2] << ","
              << t1->p3->coord[0] << "," << t1->p3->coord[1] << "," << t1->p3->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0,9.0,"
              
              << t1->p3->coord[0] << "," << t1->p3->coord[1] << "," << t1->p3->coord[2] << ","
              << t1->p1->coord[0] << "," << t1->p1->coord[1] << "," << t1->p1->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0]" << endl;
    f_out << "cmd.load_cgo(t1,'t1',1)" << endl;
    
    f_out << "t2=[9.0,"
              << t2->p0->coord[0] << "," << t2->p0->coord[1] << "," << t2->p0->coord[2] << ","
              << t2->p1->coord[0] << "," << t2->p1->coord[1] << "," << t2->p1->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0,9.0,"
              
              << t2->p0->coord[0] << "," << t2->p0->coord[1] << "," << t2->p0->coord[2] << ","
              << t2->p2->coord[0] << "," << t2->p2->coord[1] << "," << t2->p2->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0,9.0,"
              
              << t2->p0->coord[0] << "," << t2->p0->coord[1] << "," << t2->p0->coord[2] << ","
              << t2->p3->coord[0] << "," << t2->p3->coord[1] << "," << t2->p3->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0,9.0,"
              
              << t2->p1->coord[0] << "," << t2->p1->coord[1] << "," << t2->p1->coord[2] << ","
              << t2->p2->coord[0] << "," << t2->p2->coord[1] << "," << t2->p2->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0,9.0,"
              
              << t2->p2->coord[0] << "," << t2->p2->coord[1] << "," << t2->p2->coord[2] << ","
              << t2->p3->coord[0] << "," << t2->p3->coord[1] << "," << t2->p3->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0,9.0,"
              
              << t2->p3->coord[0] << "," << t2->p3->coord[1] << "," << t2->p3->coord[2] << ","
              << t2->p1->coord[0] << "," << t2->p1->coord[1] << "," << t2->p1->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0]" << endl;
    f_out << "cmd.load_cgo(t2,'t2',1)" << endl;
    
    if (!t1n0.zero()) {
    f_out << "t1n0=[9.0,"
              << t1n0->p0->coord[0] << "," << t1n0->p0->coord[1] << "," << t1n0->p0->coord[2] << ","
              << t1n0->p1->coord[0] << "," << t1n0->p1->coord[1] << "," << t1n0->p1->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0,9.0,"
              
              << t1n0->p0->coord[0] << "," << t1n0->p0->coord[1] << "," << t1n0->p0->coord[2] << ","
              << t1n0->p2->coord[0] << "," << t1n0->p2->coord[1] << "," << t1n0->p2->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0,9.0,"
              
              << t1n0->p0->coord[0] << "," << t1n0->p0->coord[1] << "," << t1n0->p0->coord[2] << ","
              << t1n0->p3->coord[0] << "," << t1n0->p3->coord[1] << "," << t1n0->p3->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0,9.0,"
              
              << t1n0->p1->coord[0] << "," << t1n0->p1->coord[1] << "," << t1n0->p1->coord[2] << ","
              << t1n0->p2->coord[0] << "," << t1n0->p2->coord[1] << "," << t1n0->p2->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0,9.0,"
              
              << t1n0->p2->coord[0] << "," << t1n0->p2->coord[1] << "," << t1n0->p2->coord[2] << ","
              << t1n0->p3->coord[0] << "," << t1n0->p3->coord[1] << "," << t1n0->p3->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0,9.0,"
              
              << t1n0->p3->coord[0] << "," << t1n0->p3->coord[1] << "," << t1n0->p3->coord[2] << ","
              << t1n0->p1->coord[0] << "," << t1n0->p1->coord[1] << "," << t1n0->p1->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0]" << endl;
    f_out << "cmd.load_cgo(t1n0,'t1n0',1)" << endl;
    }
    if (!t1n1.zero()) {
    f_out << "t1n1=[9.0,"
              << t1n1->p0->coord[0] << "," << t1n1->p0->coord[1] << "," << t1n1->p0->coord[2] << ","
              << t1n1->p1->coord[0] << "," << t1n1->p1->coord[1] << "," << t1n1->p1->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0,9.0,"
              
              << t1n1->p0->coord[0] << "," << t1n1->p0->coord[1] << "," << t1n1->p0->coord[2] << ","
              << t1n1->p2->coord[0] << "," << t1n1->p2->coord[1] << "," << t1n1->p2->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0,9.0,"
              
              << t1n1->p0->coord[0] << "," << t1n1->p0->coord[1] << "," << t1n1->p0->coord[2] << ","
              << t1n1->p3->coord[0] << "," << t1n1->p3->coord[1] << "," << t1n1->p3->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0,9.0,"
              
              << t1n1->p1->coord[0] << "," << t1n1->p1->coord[1] << "," << t1n1->p1->coord[2] << ","
              << t1n1->p2->coord[0] << "," << t1n1->p2->coord[1] << "," << t1n1->p2->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0,9.0,"
              
              << t1n1->p2->coord[0] << "," << t1n1->p2->coord[1] << "," << t1n1->p2->coord[2] << ","
              << t1n1->p3->coord[0] << "," << t1n1->p3->coord[1] << "," << t1n1->p3->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0,9.0,"
              
              << t1n1->p3->coord[0] << "," << t1n1->p3->coord[1] << "," << t1n1->p3->coord[2] << ","
              << t1n1->p1->coord[0] << "," << t1n1->p1->coord[1] << "," << t1n1->p1->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0]" << endl;
    f_out << "cmd.load_cgo(t1n1,'t1n1',1)" << endl;
    }
    if (!t1n2.zero()) {
    f_out << "t1n2=[9.0,"
              << t1n2->p0->coord[0] << "," << t1n2->p0->coord[1] << "," << t1n2->p0->coord[2] << ","
              << t1n2->p1->coord[0] << "," << t1n2->p1->coord[1] << "," << t1n2->p1->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0,9.0,"
              
              << t1n2->p0->coord[0] << "," << t1n2->p0->coord[1] << "," << t1n2->p0->coord[2] << ","
              << t1n2->p2->coord[0] << "," << t1n2->p2->coord[1] << "," << t1n2->p2->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0,9.0,"
              
              << t1n2->p0->coord[0] << "," << t1n2->p0->coord[1] << "," << t1n2->p0->coord[2] << ","
              << t1n2->p3->coord[0] << "," << t1n2->p3->coord[1] << "," << t1n2->p3->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0,9.0,"
              
              << t1n2->p1->coord[0] << "," << t1n2->p1->coord[1] << "," << t1n2->p1->coord[2] << ","
              << t1n2->p2->coord[0] << "," << t1n2->p2->coord[1] << "," << t1n2->p2->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0,9.0,"
              
              << t1n2->p2->coord[0] << "," << t1n2->p2->coord[1] << "," << t1n2->p2->coord[2] << ","
              << t1n2->p3->coord[0] << "," << t1n2->p3->coord[1] << "," << t1n2->p3->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0,9.0,"
              
              << t1n2->p3->coord[0] << "," << t1n2->p3->coord[1] << "," << t1n2->p3->coord[2] << ","
              << t1n2->p1->coord[0] << "," << t1n2->p1->coord[1] << "," << t1n2->p1->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0]" << endl;
    f_out << "cmd.load_cgo(t1n2,'t1n2',1)" << endl;
    }
    if (!t2n0.zero()) {
    f_out << "t2n0=[9.0,"
              << t2n0->p0->coord[0] << "," << t2n0->p0->coord[1] << "," << t2n0->p0->coord[2] << ","
              << t2n0->p1->coord[0] << "," << t2n0->p1->coord[1] << "," << t2n0->p1->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0,9.0,"
              
              << t2n0->p0->coord[0] << "," << t2n0->p0->coord[1] << "," << t2n0->p0->coord[2] << ","
              << t2n0->p2->coord[0] << "," << t2n0->p2->coord[1] << "," << t2n0->p2->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0,9.0,"
              
              << t2n0->p0->coord[0] << "," << t2n0->p0->coord[1] << "," << t2n0->p0->coord[2] << ","
              << t2n0->p3->coord[0] << "," << t2n0->p3->coord[1] << "," << t2n0->p3->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0,9.0,"
              
              << t2n0->p1->coord[0] << "," << t2n0->p1->coord[1] << "," << t2n0->p1->coord[2] << ","
              << t2n0->p2->coord[0] << "," << t2n0->p2->coord[1] << "," << t2n0->p2->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0,9.0,"
              
              << t2n0->p2->coord[0] << "," << t2n0->p2->coord[1] << "," << t2n0->p2->coord[2] << ","
              << t2n0->p3->coord[0] << "," << t2n0->p3->coord[1] << "," << t2n0->p3->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0,9.0,"
              
              << t2n0->p3->coord[0] << "," << t2n0->p3->coord[1] << "," << t2n0->p3->coord[2] << ","
              << t2n0->p1->coord[0] << "," << t2n0->p1->coord[1] << "," << t2n0->p1->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0]" << endl;
    f_out << "cmd.load_cgo(t2n0,'t2n0',1)" << endl;
    }
    if (!t2n1.zero()) {
    f_out << "t2n1=[9.0,"
              << t2n1->p0->coord[0] << "," << t2n1->p0->coord[1] << "," << t2n1->p0->coord[2] << ","
              << t2n1->p1->coord[0] << "," << t2n1->p1->coord[1] << "," << t2n1->p1->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0,9.0,"
              
              << t2n1->p0->coord[0] << "," << t2n1->p0->coord[1] << "," << t2n1->p0->coord[2] << ","
              << t2n1->p2->coord[0] << "," << t2n1->p2->coord[1] << "," << t2n1->p2->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0,9.0,"
              
              << t2n1->p0->coord[0] << "," << t2n1->p0->coord[1] << "," << t2n1->p0->coord[2] << ","
              << t2n1->p3->coord[0] << "," << t2n1->p3->coord[1] << "," << t2n1->p3->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0,9.0,"
              
              << t2n1->p1->coord[0] << "," << t2n1->p1->coord[1] << "," << t2n1->p1->coord[2] << ","
              << t2n1->p2->coord[0] << "," << t2n1->p2->coord[1] << "," << t2n1->p2->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0,9.0,"
              
              << t2n1->p2->coord[0] << "," << t2n1->p2->coord[1] << "," << t2n1->p2->coord[2] << ","
              << t2n1->p3->coord[0] << "," << t2n1->p3->coord[1] << "," << t2n1->p3->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0,9.0,"
              
              << t2n1->p3->coord[0] << "," << t2n1->p3->coord[1] << "," << t2n1->p3->coord[2] << ","
              << t2n1->p1->coord[0] << "," << t2n1->p1->coord[1] << "," << t2n1->p1->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0]" << endl;
    f_out << "cmd.load_cgo(t2n1,'t2n1',1)" << endl;
    }
    if (!t2n2.zero()) {
    f_out << "t2n2=[9.0,"
              << t2n2->p0->coord[0] << "," << t2n2->p0->coord[1] << "," << t2n2->p0->coord[2] << ","
              << t2n2->p1->coord[0] << "," << t2n2->p1->coord[1] << "," << t2n2->p1->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0,9.0,"
              
              << t2n2->p0->coord[0] << "," << t2n2->p0->coord[1] << "," << t2n2->p0->coord[2] << ","
              << t2n2->p2->coord[0] << "," << t2n2->p2->coord[1] << "," << t2n2->p2->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0,9.0,"
              
              << t2n2->p0->coord[0] << "," << t2n2->p0->coord[1] << "," << t2n2->p0->coord[2] << ","
              << t2n2->p3->coord[0] << "," << t2n2->p3->coord[1] << "," << t2n2->p3->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0,9.0,"
              
              << t2n2->p1->coord[0] << "," << t2n2->p1->coord[1] << "," << t2n2->p1->coord[2] << ","
              << t2n2->p2->coord[0] << "," << t2n2->p2->coord[1] << "," << t2n2->p2->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0,9.0,"
              
              << t2n2->p2->coord[0] << "," << t2n2->p2->coord[1] << "," << t2n2->p2->coord[2] << ","
              << t2n2->p3->coord[0] << "," << t2n2->p3->coord[1] << "," << t2n2->p3->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0,9.0,"
              
              << t2n2->p3->coord[0] << "," << t2n2->p3->coord[1] << "," << t2n2->p3->coord[2] << ","
              << t2n2->p1->coord[0] << "," << t2n2->p1->coord[1] << "," << t2n2->p1->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0]" << endl;
    f_out << "cmd.load_cgo(t2n2,'t2n2',1)" << endl;
    }
    */
    f_out.close();
}
//!===============================================================================


//!===============================================================================
//!Definition der Klasse DT_SOLVER:
DT_SOLVER::DT_SOLVER(vector<stl_ptr<DT_POINT> > &dt_pnt) {
    for (dtp_vec it=dt_pnt.begin(); it!=dt_pnt.end(); ++it) {
        stl_ptr<DT_POINT> np(new DT_POINT(**it));
        dt_points.push_back(np);
    }
    get_initial_tetrahedron();
    curr_id = 0;
    curr_tet_id = 1;
}

DT_SOLVER::DT_SOLVER(vector<DT_POINT> &dt_pnt) {
    for (vector<DT_POINT>::iterator it=dt_pnt.begin(); it!=dt_pnt.end(); ++it) {
        stl_ptr<DT_POINT> np(new DT_POINT(*it));
        dt_points.push_back(np);
    }
    get_initial_tetrahedron();
    curr_id = 0;
    curr_tet_id = 1;
}

DT_SOLVER::~DT_SOLVER() {
    //!zunaechst die Punkte des initialen Tetraeders killen:
    tetrahedrons[0]->p0.kill();
    tetrahedrons[0]->p1.kill();
    tetrahedrons[0]->p2.kill();
    tetrahedrons[0]->p3.kill();
    for (th_vec it=tetrahedrons.begin(); it!=tetrahedrons.end(); ++it) {
        it->kill();
    }
    for (dtp_vec it=dt_points.begin(); it!=dt_points.end(); ++it) {
        it->kill();
    }
}

void DT_SOLVER::visualize_tetrahedrons() {
    ofstream f_out;
    f_out.open("tetrahedrons.py");
    
    float rad = 0.04;
    
    f_out << "tet=[";
    bool first = true;
    for (th_vec it=tetrahedrons.begin(); it!=tetrahedrons.end(); ++it) {
        if ((*it)->next.size() > 0) continue;
        if ((*it)->p0.equal_addr(tetrahedrons[0]->p0)) continue;
        else if ((*it)->p0.equal_addr(tetrahedrons[0]->p1)) continue;
        else if ((*it)->p0.equal_addr(tetrahedrons[0]->p2)) continue;
        else if ((*it)->p0.equal_addr(tetrahedrons[0]->p3)) continue;
        else if ((*it)->p1.equal_addr(tetrahedrons[0]->p0)) continue;
        else if ((*it)->p1.equal_addr(tetrahedrons[0]->p1)) continue;
        else if ((*it)->p1.equal_addr(tetrahedrons[0]->p2)) continue;
        else if ((*it)->p1.equal_addr(tetrahedrons[0]->p3)) continue;
        else if ((*it)->p2.equal_addr(tetrahedrons[0]->p0)) continue;
        else if ((*it)->p2.equal_addr(tetrahedrons[0]->p1)) continue;
        else if ((*it)->p2.equal_addr(tetrahedrons[0]->p2)) continue;
        else if ((*it)->p2.equal_addr(tetrahedrons[0]->p3)) continue;
        else if ((*it)->p3.equal_addr(tetrahedrons[0]->p0)) continue;
        else if ((*it)->p3.equal_addr(tetrahedrons[0]->p1)) continue;
        else if ((*it)->p3.equal_addr(tetrahedrons[0]->p2)) continue;
        else if ((*it)->p3.equal_addr(tetrahedrons[0]->p3)) continue;
        
        if (first) {
            f_out << "9.0,";
            first = false;
        } else f_out << ",9.0,";
        
        f_out << (*it)->p0->coord[0] << "," << (*it)->p0->coord[1] << "," << (*it)->p0->coord[2] << ","
              << (*it)->p1->coord[0] << "," << (*it)->p1->coord[1] << "," << (*it)->p1->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0,9.0,"
              
              << (*it)->p0->coord[0] << "," << (*it)->p0->coord[1] << "," << (*it)->p0->coord[2] << ","
              << (*it)->p2->coord[0] << "," << (*it)->p2->coord[1] << "," << (*it)->p2->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0,9.0,"
              
              << (*it)->p0->coord[0] << "," << (*it)->p0->coord[1] << "," << (*it)->p0->coord[2] << ","
              << (*it)->p3->coord[0] << "," << (*it)->p3->coord[1] << "," << (*it)->p3->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0,9.0,"
              
              << (*it)->p1->coord[0] << "," << (*it)->p1->coord[1] << "," << (*it)->p1->coord[2] << ","
              << (*it)->p2->coord[0] << "," << (*it)->p2->coord[1] << "," << (*it)->p2->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0,9.0,"
              
              << (*it)->p2->coord[0] << "," << (*it)->p2->coord[1] << "," << (*it)->p2->coord[2] << ","
              << (*it)->p3->coord[0] << "," << (*it)->p3->coord[1] << "," << (*it)->p3->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0,9.0,"
              
              << (*it)->p3->coord[0] << "," << (*it)->p3->coord[1] << "," << (*it)->p3->coord[2] << ","
              << (*it)->p1->coord[0] << "," << (*it)->p1->coord[1] << "," << (*it)->p1->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0";
    }
    
    f_out << "]" << endl;
    
    f_out << "cmd.load_cgo(tet,'tet',1)" << endl;
    
    f_out.close();
}

void DT_SOLVER::visualize_tetrahedrons(vector<stl_ptr<TETRAHEDRON> > &tets,string fname,string pname) {
    ofstream f_out;
    if (fname == "X") f_out.open("alpha_tets.py");
    else f_out.open(fname.c_str());
    
    float rad = 0.04;
    
    f_out << "tet=[";
    bool first = true;
    for (th_vec it=tets.begin(); it!=tets.end(); ++it) {
        if ((*it)->next.size() > 0) continue;
        if ((*it)->p0.equal_addr(tetrahedrons[0]->p0)) continue;
        else if ((*it)->p0.equal_addr(tetrahedrons[0]->p1)) continue;
        else if ((*it)->p0.equal_addr(tetrahedrons[0]->p2)) continue;
        else if ((*it)->p0.equal_addr(tetrahedrons[0]->p3)) continue;
        else if ((*it)->p1.equal_addr(tetrahedrons[0]->p0)) continue;
        else if ((*it)->p1.equal_addr(tetrahedrons[0]->p1)) continue;
        else if ((*it)->p1.equal_addr(tetrahedrons[0]->p2)) continue;
        else if ((*it)->p1.equal_addr(tetrahedrons[0]->p3)) continue;
        else if ((*it)->p2.equal_addr(tetrahedrons[0]->p0)) continue;
        else if ((*it)->p2.equal_addr(tetrahedrons[0]->p1)) continue;
        else if ((*it)->p2.equal_addr(tetrahedrons[0]->p2)) continue;
        else if ((*it)->p2.equal_addr(tetrahedrons[0]->p3)) continue;
        else if ((*it)->p3.equal_addr(tetrahedrons[0]->p0)) continue;
        else if ((*it)->p3.equal_addr(tetrahedrons[0]->p1)) continue;
        else if ((*it)->p3.equal_addr(tetrahedrons[0]->p2)) continue;
        else if ((*it)->p3.equal_addr(tetrahedrons[0]->p3)) continue;
        
        if (first) {
            f_out << "9.0,";
            first = false;
        } else f_out << ",9.0,";
        
        f_out << (*it)->p0->coord[0] << "," << (*it)->p0->coord[1] << "," << (*it)->p0->coord[2] << ","
              << (*it)->p1->coord[0] << "," << (*it)->p1->coord[1] << "," << (*it)->p1->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0,9.0,"
              
              << (*it)->p0->coord[0] << "," << (*it)->p0->coord[1] << "," << (*it)->p0->coord[2] << ","
              << (*it)->p2->coord[0] << "," << (*it)->p2->coord[1] << "," << (*it)->p2->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0,9.0,"
              
              << (*it)->p0->coord[0] << "," << (*it)->p0->coord[1] << "," << (*it)->p0->coord[2] << ","
              << (*it)->p3->coord[0] << "," << (*it)->p3->coord[1] << "," << (*it)->p3->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0,9.0,"
              
              << (*it)->p1->coord[0] << "," << (*it)->p1->coord[1] << "," << (*it)->p1->coord[2] << ","
              << (*it)->p2->coord[0] << "," << (*it)->p2->coord[1] << "," << (*it)->p2->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0,9.0,"
              
              << (*it)->p2->coord[0] << "," << (*it)->p2->coord[1] << "," << (*it)->p2->coord[2] << ","
              << (*it)->p3->coord[0] << "," << (*it)->p3->coord[1] << "," << (*it)->p3->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0,9.0,"
              
              << (*it)->p3->coord[0] << "," << (*it)->p3->coord[1] << "," << (*it)->p3->coord[2] << ","
              << (*it)->p1->coord[0] << "," << (*it)->p1->coord[1] << "," << (*it)->p1->coord[2] << ","
              << rad << ",1.0,1.0,0.0,1.0,1.0,0.0";
    }
    
    f_out << "]" << endl;
    
    if (pname != "X") f_out << "cmd.load_cgo(tet,'" << pname << "',1)" << endl;
    else f_out << "cmd.load_cgo(tet,'tet',1)" << endl;
    
    f_out.close();
}

void DT_SOLVER::visualize_convex_hull(string fname) {
    ofstream f_out;
    
    if (fname == "X") f_out.open("convex_hull.py");
    else f_out.open(fname.c_str());
    
    float rad = 0.04;
    
    vector<stl_ptr<TETRAHEDRON> > hull_tet;
    
    for (th_vec it=tetrahedrons.begin(); it!=tetrahedrons.end(); ++it) {
        if ((*it)->next.size() > 0) continue;
        if ((*it)->p0.equal_addr(tetrahedrons[0]->p0) ||
            (*it)->p0.equal_addr(tetrahedrons[0]->p1) ||
            (*it)->p0.equal_addr(tetrahedrons[0]->p2) ||
            (*it)->p0.equal_addr(tetrahedrons[0]->p3) ||
            (*it)->p1.equal_addr(tetrahedrons[0]->p0) ||
            (*it)->p1.equal_addr(tetrahedrons[0]->p1) ||
            (*it)->p1.equal_addr(tetrahedrons[0]->p2) ||
            (*it)->p1.equal_addr(tetrahedrons[0]->p3) ||
            (*it)->p2.equal_addr(tetrahedrons[0]->p0) ||
            (*it)->p2.equal_addr(tetrahedrons[0]->p1) ||
            (*it)->p2.equal_addr(tetrahedrons[0]->p2) ||
            (*it)->p2.equal_addr(tetrahedrons[0]->p3) ||
            (*it)->p3.equal_addr(tetrahedrons[0]->p0) ||
            (*it)->p3.equal_addr(tetrahedrons[0]->p1) ||
            (*it)->p3.equal_addr(tetrahedrons[0]->p2) ||
            (*it)->p3.equal_addr(tetrahedrons[0]->p3)) hull_tet.push_back(*it);
    }
    
    f_out << "hull=[";
    bool first = true;
    
    for (th_vec it=hull_tet.begin(); it!=hull_tet.end(); ++it) {
        int sp = 0;
        vec3d<float> p[3];
        if (!(*it)->p0.equal_addr(tetrahedrons[0]->p0) &&
            !(*it)->p0.equal_addr(tetrahedrons[0]->p1) &&
            !(*it)->p0.equal_addr(tetrahedrons[0]->p2) &&
            !(*it)->p0.equal_addr(tetrahedrons[0]->p3)) {
            p[sp] = (*it)->p0->coord;
            sp++;
        }
        
        if (!(*it)->p1.equal_addr(tetrahedrons[0]->p0) &&
            !(*it)->p1.equal_addr(tetrahedrons[0]->p1) &&
            !(*it)->p1.equal_addr(tetrahedrons[0]->p2) &&
            !(*it)->p1.equal_addr(tetrahedrons[0]->p3)) {
            p[sp] = (*it)->p1->coord;
            sp++;
        }
        
        if (!(*it)->p2.equal_addr(tetrahedrons[0]->p0) &&
            !(*it)->p2.equal_addr(tetrahedrons[0]->p1) &&
            !(*it)->p2.equal_addr(tetrahedrons[0]->p2) &&
            !(*it)->p2.equal_addr(tetrahedrons[0]->p3)) {
            p[sp] = (*it)->p2->coord;
            sp++;
        }
        
        if (!(*it)->p3.equal_addr(tetrahedrons[0]->p0) &&
            !(*it)->p3.equal_addr(tetrahedrons[0]->p1) &&
            !(*it)->p3.equal_addr(tetrahedrons[0]->p2) &&
            !(*it)->p3.equal_addr(tetrahedrons[0]->p3)) {
            p[sp] = (*it)->p3->coord;
            sp++;
        }
        
        if (sp == 3) {
            if (first) {
                f_out << "9.0,";
                first = false;
            } else f_out << ",9.0,";
            f_out << p[0][0] << "," << p[0][1] << "," << p[0][2] << ","
                      << p[1][0] << "," << p[1][1] << "," << p[1][2] << ","
                      << rad << ",1.0,0.0,0.0,1.0,0.0,0.0,9.0,"
                      
                      << p[0][0] << "," << p[0][1] << "," << p[0][2] << ","
                      << p[2][0] << "," << p[2][1] << "," << p[2][2] << ","
                      << rad << ",1.0,0.0,0.0,1.0,0.0,0.0,9.0,"
                      
                      << p[2][0] << "," << p[2][1] << "," << p[2][2] << ","
                      << p[1][0] << "," << p[1][1] << "," << p[1][2] << ","
                      << rad << ",1.0,0.0,0.0,1.0,0.0,0.0";
        } else if (sp == 2) {
            if (first) {
                f_out << "9.0,";
                first = false;
            } else f_out << ",9.0,";
            f_out << p[0][0] << "," << p[0][1] << "," << p[0][2] << ","
                      << p[1][0] << "," << p[1][1] << "," << p[1][2] << ","
                      << rad << ",1.0,0.0,0.0,1.0,0.0,0.0";
        }
    }
    
    f_out << "]" << endl;
    
    f_out << "cmd.load_cgo(hull,'hull',1)" << endl;
    
    f_out.close();
}

void DT_SOLVER::visualize_alpha_spheres_new(string fname,double alpha) {
    ofstream f_out;
    
    if (fname == "X") f_out.open("alpha_spheres.py");
    else f_out.open(fname.c_str());
    
    f_out << "aspheres = [";
    
    bool first = true;
    
    for (vector<vec3d<double> >::iterator it=alpha_spheres.begin(); it!=alpha_spheres.end(); ++it) {
        
        if (first) {
            first = false;
        } else f_out << ",";
        
        f_out << "7.0," << (*it)[0] << "," << (*it)[1] << "," << (*it)[2] << "," << alpha;
    }
    
    f_out << "]" << endl;
    
    f_out << "cmd.load_cgo(aspheres,'aspheres',1)" << endl;
    
    f_out.close();
}

void DT_SOLVER::visualize_alpha_facets(string fname) {
    ofstream f_out;
    float rad = 0.04;
    if (fname == "X") f_out.open("alpha_facets.py");
    else f_out.open(fname.c_str());
    
    f_out << "afacets = [";
    bool first = true;
    for (vector<stl_ptr<L_FACET> >::iterator it=alpha_facets.begin(); it!=alpha_facets.end(); ++it) {
        if (first) {
            first = false;
        } else f_out << ",";
        f_out << "9.0,";
        f_out << (*it)->p0->coord[0] << "," << (*it)->p0->coord[1] << "," << (*it)->p0->coord[2] << ","
              << (*it)->p1->coord[0] << "," << (*it)->p1->coord[1] << "," << (*it)->p1->coord[2] << ","
              << rad << ",1.0,0.0,0.0,1.0,0.0,0.0,9.0,";
        f_out << (*it)->p0->coord[0] << "," << (*it)->p0->coord[1] << "," << (*it)->p0->coord[2] << ","
              << (*it)->p2->coord[0] << "," << (*it)->p2->coord[1] << "," << (*it)->p2->coord[2] << ","
              << rad << ",1.0,0.0,0.0,1.0,0.0,0.0,9.0,";
        f_out << (*it)->p2->coord[0] << "," << (*it)->p2->coord[1] << "," << (*it)->p2->coord[2] << ","
              << (*it)->p1->coord[0] << "," << (*it)->p1->coord[1] << "," << (*it)->p1->coord[2] << ","
              << rad << ",1.0,0.0,0.0,1.0,0.0,0.0";
    }
    f_out << "]" << endl;
    
    f_out << "cmd.load_cgo(afacets,'afacets',1)" << endl;
    
    f_out.close();
}

void DT_SOLVER::visualize_facets(vector<stl_ptr<L_FACET> > &fcs,string fname) {
    ofstream f_out;
    float rad = 0.04;
    if (fname == "X") f_out.open("facets.py");
    else f_out.open(fname.c_str());
    
    f_out << "facets = [";
    bool first = true;
    for (vector<stl_ptr<L_FACET> >::iterator it=fcs.begin(); it!=fcs.end(); ++it) {
        if (first) {
            first = false;
        } else f_out << ",";
        f_out << "9.0,";
        f_out << (*it)->p0->coord[0] << "," << (*it)->p0->coord[1] << "," << (*it)->p0->coord[2] << ","
              << (*it)->p1->coord[0] << "," << (*it)->p1->coord[1] << "," << (*it)->p1->coord[2] << ","
              << rad << ",1.0,0.0,0.0,1.0,0.0,0.0,9.0,";
        f_out << (*it)->p0->coord[0] << "," << (*it)->p0->coord[1] << "," << (*it)->p0->coord[2] << ","
              << (*it)->p2->coord[0] << "," << (*it)->p2->coord[1] << "," << (*it)->p2->coord[2] << ","
              << rad << ",1.0,0.0,0.0,1.0,0.0,0.0,9.0,";
        f_out << (*it)->p2->coord[0] << "," << (*it)->p2->coord[1] << "," << (*it)->p2->coord[2] << ","
              << (*it)->p1->coord[0] << "," << (*it)->p1->coord[1] << "," << (*it)->p1->coord[2] << ","
              << rad << ",1.0,0.0,0.0,1.0,0.0,0.0";
    }
    f_out << "]" << endl;
    
    f_out << "cmd.load_cgo(facets,'facets',1)" << endl;
    
    f_out.close();
}

void DT_SOLVER::visualize_alpha_spheres() {
    ofstream f_out;
    f_out.open("alpha_spheres.py");
    
    float max_srad = 100.;
    float min_srad = 48.;
    
    f_out << "aspheres = [";
    
    bool first = true;
    
    for (th_vec it=delaunay_tetrahedrons.begin(); it!=delaunay_tetrahedrons.end(); ++it) {
        
        
        if ((*it)->sphere_square_radius > max_srad) continue;
        if ((*it)->sphere_square_radius < min_srad) continue;
        
        
        if (!(is_in_hull((*it)->sphere_middle))) continue;
        
        if (first) {
            first = false;
        } else f_out << ",";
        
        f_out << "7.0," << (*it)->sphere_middle[0] << "," << (*it)->sphere_middle[1] << "," << (*it)->sphere_middle[2] 
              << "," << sqrt((*it)->sphere_square_radius);
    }
    
    f_out << "]" << endl;
    
    f_out << "cmd.load_cgo(aspheres,'aspheres',1)" << endl;
    
    f_out.close();
}

void DT_SOLVER::visualize_spheres(vector<stl_ptr<TETRAHEDRON> > &tets) {
    ofstream f_out;
    f_out.open("spheres.py");
    
    f_out << "spheres = [";
    
    bool first = true;
    for (th_vec it=tets.begin(); it!=tets.end(); ++it) {
        if (first) {
            first = false;
        } else f_out << ",";
        f_out << "7.0," << (*it)->sphere_middle[0] << "," << (*it)->sphere_middle[1] << "," << (*it)->sphere_middle[2] 
              << "," << sqrt((*it)->sphere_square_radius);
    }
    f_out << "]" << endl;
    
    f_out << "cmd.load_cgo(spheres,'spheres',1)" << endl;
    f_out.close();
}

void DT_SOLVER::visualize_sphere_maps(map<int,vector<stl_ptr<TETRAHEDRON> > > &cluster,string fname) {
    ofstream f_out;
    
    if (fname == "X") f_out.open("cluster.py");
    else f_out.open(fname.c_str());
    
    int state = 1;
    
    for (map<int,vector<stl_ptr<TETRAHEDRON> > >::iterator jt=cluster.begin(); jt!=cluster.end(); ++jt) {
        f_out << "spheres" << state << " = [";
        
        bool first = true;
        for (th_vec it=jt->second.begin(); it!=jt->second.end(); ++it) {
            if (first) {
                first = false;
            } else f_out << ",";
            f_out << "7.0," << (*it)->sphere_middle[0] << "," << (*it)->sphere_middle[1] << "," << (*it)->sphere_middle[2] 
            << "," << sqrt((*it)->sphere_square_radius);
        }
        f_out << "]" << endl;
        
        f_out << "cmd.load_cgo(spheres" << state << ",'spheres" << state << "',1)" << endl;
        
        ++state;
    }
    
    f_out.close();
}

void DT_SOLVER::get_initial_tetrahedron() {
    double const unlimited = 30000.; //!Wenn diese Zahl zu gross gewaehlt wird (relativ zu den vorkommenden Koordinaten,
                                     //!dann sind die dem ersten 1->4er Flip folgenden Punkte fast immer koplanar mit
                                     //!einer der neuen Tetraderseiten (oder mehreren)
    vec3d<double> p1(-unlimited,-unlimited,-unlimited);
    vec3d<double> p2(unlimited,-unlimited,-unlimited);
    vec3d<double> p3(0.,unlimited,-unlimited);
    vec3d<double> p4(0.,0.,unlimited);
    stl_ptr<DT_POINT> v1(new DT_POINT(-4,p1));
    stl_ptr<DT_POINT> v2(new DT_POINT(-3,p2));
    stl_ptr<DT_POINT> v3(new DT_POINT(-2,p3));
    stl_ptr<DT_POINT> v4(new DT_POINT(-1,p4));
    stl_ptr<TETRAHEDRON> ith(new TETRAHEDRON(v1,v2,v3,v4,0,0));
    tetrahedrons.push_back(ith);
}

stl_ptr<TETRAHEDRON>& DT_SOLVER::locate_tetrahedron(stl_ptr<DT_POINT> const& point,stl_ptr<TETRAHEDRON> &root) {
    //!Rekursive Funktion um zu ermitteln in welchem Tetraeder point liegt
    //!=> mit tetrahedrons[0] aufrufen
    if (root->next.empty()) return root; //!Das Tetraeder wurde gefunden
    for (th_vec it=root->next.begin(); it!=root->next.end(); ++it) {
        if ((*it)->contains_p(point)) return locate_tetrahedron(point,*it);
    }
    
    cerr << "error: point localization failed for " << point->coord << " !" << endl;
    cerr << "       trying hardcore method..." << endl;
    
    //!Wurde kein Tetraeder gefunden, so ist entweder das allumfassende Tetraeder zu klein gewaehlt oder
    //!es liegt ein Problem mit der numerischen Genauigkeit vor.
    //!Jetzt die brutale Methode:  Alle aktuellen Tetraeder durchgehen und gucken in welchem es am ehesten liegt:
    //root->visualize(1); int vi = 2;
    cerr << root << endl;
    for (th_vec it=root->next.begin(); it!=root->next.end(); ++it) {
        cerr << triple_product((*it)->p0->coord,(*it)->p1->coord,(*it)->p2->coord,point->coord) << endl;
        cerr << triple_product((*it)->p0->coord,(*it)->p3->coord,(*it)->p1->coord,point->coord) << endl;
        cerr << triple_product((*it)->p0->coord,(*it)->p2->coord,(*it)->p3->coord,point->coord) << endl;
        cerr << triple_product((*it)->p1->coord,(*it)->p3->coord,(*it)->p2->coord,point->coord) << endl;
        cerr << endl;
    //    (*it)->visualize(vi); vi++;
        cerr << *it << endl;
    }
    
    for (th_vec it=tetrahedrons.begin(); it!=tetrahedrons.end(); ++it) {
        if ((*it)->next.size() > 0) continue;
        bool hit = true;
        if (triple_product((*it)->p0->coord,(*it)->p1->coord,(*it)->p2->coord,point->coord) <= 0.) hit = false;
        if (triple_product((*it)->p0->coord,(*it)->p3->coord,(*it)->p1->coord,point->coord) <= 0.) hit = false;
        if (triple_product((*it)->p0->coord,(*it)->p2->coord,(*it)->p3->coord,point->coord) <= 0.) hit = false;
        if (triple_product((*it)->p1->coord,(*it)->p3->coord,(*it)->p2->coord,point->coord) <= 0.) hit = false;
        
        if (hit) return *it;
    }
    
    cerr << "error: point localization failed again!" << endl;
    exit(1);
    return zero;
}

void DT_SOLVER::flip14(stl_ptr<DT_POINT> const& point) {
    //! Ein neuer Punkt wurde der aktuellen Triangulierung hinzugefuegt und das Tetraeder in
    //! welchem dieser Punkt liegt wurde lokalisiert => 'curr_tet' zeigt auf dieses Tetraeder.
    //! Als initialer Schritt zur Wiederherstellung der Triangulierung wird nun ein 1->4 Flip
    //! ausgefuehrt. Hierbei entstehen aus dem urspruenglichen Tetraeder mit den Punkten 0 bis 3
    //! und dem neuen Punkt 'point' 4 neue Tetraeder.
    //! Um die Triangulierung in weiteren Schritten wieder regulaer zu machen werden die
    //! entstandenen Link-Facets (Dreiecksflaechen, die neu entstandenen Tetraeder an andere
    //! Tetraeder angrenzen) zunaechst auf einem Stapel 'facets' abgelegt.
    
    ++curr_generation;
    
    //! 1.) Mit point 4 neue Tetraeder innerhalb von curr_tet erzeugen:
    stl_ptr<TETRAHEDRON> nt0(new TETRAHEDRON(curr_tet->p0,curr_tet->p3,curr_tet->p1,point,curr_generation,curr_tet_id)); curr_tet_id++;
    stl_ptr<TETRAHEDRON> nt1(new TETRAHEDRON(curr_tet->p0,curr_tet->p2,curr_tet->p3,point,curr_generation,curr_tet_id)); curr_tet_id++;
    stl_ptr<TETRAHEDRON> nt2(new TETRAHEDRON(curr_tet->p1,curr_tet->p3,curr_tet->p2,point,curr_generation,curr_tet_id)); curr_tet_id++;
    stl_ptr<TETRAHEDRON> nt3(new TETRAHEDRON(curr_tet->p0,curr_tet->p1,curr_tet->p2,point,curr_generation,curr_tet_id)); curr_tet_id++;
    
    
    /*
    if (nt0->coplanar()) cerr << "nt0 ist coplanar" << endl;
    if (nt1->coplanar()) cerr << "nt1 ist coplanar" << endl;
    if (nt2->coplanar()) cerr << "nt2 ist coplanar" << endl;
    if (nt3->coplanar()) cerr << "nt3 ist coplanar" << endl;
    cerr << "nt_id = " << nt0->id << "   det = " << get_4point_det(nt0->p0->coord,nt0->p1->coord,nt0->p2->coord,nt0->p3->coord) << endl;
    cerr << "nt_id = " << nt1->id << "   det = " << get_4point_det(nt1->p0->coord,nt1->p1->coord,nt1->p2->coord,nt1->p3->coord) << endl;
    cerr << "nt_id = " << nt2->id << "   det = " << get_4point_det(nt2->p0->coord,nt2->p1->coord,nt2->p2->coord,nt2->p3->coord) << endl;
    cerr << "nt_id = " << nt3->id << "   det = " << get_4point_det(nt3->p0->coord,nt3->p1->coord,nt3->p2->coord,nt3->p3->coord) << endl;
    */
    
    
    //! 2.) Die neuen Tetraeder zufuegen:
    tetrahedrons.push_back(nt0);
    tetrahedrons.push_back(nt1);
    tetrahedrons.push_back(nt2);
    tetrahedrons.push_back(nt3);
    
    //! 3.) next und prev entsprechend fuer alle Tetraeder anpassen:
    curr_tet->next.push_back(nt0);
    curr_tet->next.push_back(nt1);
    curr_tet->next.push_back(nt2);
    curr_tet->next.push_back(nt3);
    nt0->prev.push_back(curr_tet);
    nt1->prev.push_back(curr_tet);
    nt2->prev.push_back(curr_tet);
    nt3->prev.push_back(curr_tet);
    
    //! 4.) p0n bis p3n fuer alle 4 Tetraeder aendern:
    nt0->p0n = nt2; nt0->p1n = nt3; nt0->p2n = nt1; nt0->p3n = curr_tet->p2n; //gegenueberliegende Tetraeder von nt0
    nt0->p0p = curr_tet->p2; nt0->p1p = curr_tet->p2; nt0->p2p = curr_tet->p2; nt0->p3p = curr_tet->p2p; //gegenueberliegende Punkte
    
    nt1->p0n = nt2; nt1->p1n = nt0; nt1->p2n = nt3; nt1->p3n = curr_tet->p1n;
    nt1->p0p = curr_tet->p1; nt1->p1p = curr_tet->p1; nt1->p2p = curr_tet->p1; nt1->p3p = curr_tet->p1p;
    
    nt2->p0n = nt1; nt2->p1n = nt3; nt2->p2n = nt0; nt2->p3n = curr_tet->p0n;
    nt2->p0p = curr_tet->p0; nt2->p1p = curr_tet->p0; nt2->p2p = curr_tet->p0; nt2->p3p = curr_tet->p0p;
    
    nt3->p0n = nt2; nt3->p1n = nt1; nt3->p2n = nt0; nt3->p3n = curr_tet->p3n;
    nt3->p0p = curr_tet->p3; nt3->p1p = curr_tet->p3; nt3->p2p = curr_tet->p3; nt3->p3p = curr_tet->p3p;
    
    
    nt0->relink_neighbours();    
    nt1->relink_neighbours();    
    nt2->relink_neighbours();    
    nt3->relink_neighbours();    
    
    //! 5.) Die 0 bis 4 Link-Facets auf den Stapel schieben:
    nt0->generate_facets(facets,curr_id);
    nt1->generate_facets(facets,curr_id);
    nt2->generate_facets(facets,curr_id);
    nt3->generate_facets(facets,curr_id);
}

void DT_SOLVER::flip32(stl_ptr<DT_POINT> const& kp1,stl_ptr<DT_POINT> const& kp2,stl_ptr<DT_POINT> const& kp3,
                       stl_ptr<TETRAHEDRON> const& t1, stl_ptr<TETRAHEDRON> const& t2,
                       stl_ptr<TETRAHEDRON> const& lt1n0, stl_ptr<TETRAHEDRON> const& lt2n0,
                       stl_ptr<TETRAHEDRON> const& lt1n1, stl_ptr<TETRAHEDRON> const& lt2n1,
                       stl_ptr<DT_POINT> const& lt1p0, stl_ptr<DT_POINT> const& lt2p0,
                       stl_ptr<DT_POINT> const& lt1p1, stl_ptr<DT_POINT> const& lt2p1) {
    //! Ein Facet der neu entstandenen Tetraeder hat eine reflexe Kante 'kp1','kp2'
    //! Wenn 3 Tetraeder an diese Kante angrenzen, so werden 2 Tetraeder daraus gemacht.
    //! Gegeben seien 3 Tetraeder mit den Punkten  1,2,3,4   5,2,4,3   1,2,4,5
    //! Die Punkte  2,4  definieren eine reflexe Kante => Die 3 Tetraeder werden zu
    //! 2 neuen vereint mit den Punkten  2,1,5,3   4,1,3,5
    
    ++curr_generation;
    
    //! 1.) Pruefen, ob genau 3 Tetraeder an die reflexe Kante angrenzen:
    if (t1.zero()) return;
    if (!t1.equal_addr(t2)) return;
    
    //! 2.) Aus den 3 Tetraedern 2 neue machen:
    stl_ptr<TETRAHEDRON> nt0(new TETRAHEDRON(kp1,cf->ps1,cf->ps2,kp3,curr_generation,curr_tet_id)); curr_tet_id++;
    stl_ptr<TETRAHEDRON> nt1(new TETRAHEDRON(kp2,cf->ps2,cf->ps1,kp3,curr_generation,curr_tet_id)); curr_tet_id++;
    
    tetrahedrons.push_back(nt0);
    tetrahedrons.push_back(nt1);
    
    //! 3.) next und prev anpassen:
    cf->t1->next.push_back(nt0);
    cf->t1->next.push_back(nt1);
    cf->t2->next.push_back(nt0);
    cf->t2->next.push_back(nt1);
    t1->next.push_back(nt0);
    t1->next.push_back(nt1);
    nt0->prev.push_back(cf->t1);
    nt1->prev.push_back(cf->t1);
    nt0->prev.push_back(cf->t2);
    nt1->prev.push_back(cf->t2);
    nt0->prev.push_back(t1);
    nt1->prev.push_back(t1);
    
    
    nt0->p0n = nt1; nt0->p1n = lt2n1; nt0->p2n = lt1n1;
    nt0->p0p = kp2; nt0->p1p = lt2p1; nt0->p2p = lt1p1;
    if (t1->p0.equal_addr(kp2)) {nt0->p3n = t1->p0n; nt0->p3p = t1->p0p;}
    else if (t1->p1.equal_addr(kp2)) {nt0->p3n = t1->p1n; nt0->p3p = t1->p1p;}
    else if (t1->p2.equal_addr(kp2)) {nt0->p3n = t1->p2n; nt0->p3p = t1->p2p;}
    else {nt0->p3n = t1->p3n; nt0->p3p = t1->p3p;}
    
    nt1->p0n = nt0; nt1->p1n = lt1n0; nt1->p2n = lt2n0;
    nt1->p0p = kp1; nt1->p1p = lt1p0; nt1->p2p = lt2p0;

    if (t1->p0.equal_addr(kp1)) {nt1->p3n = t1->p0n; nt1->p3p = t1->p0p;}
    else if (t1->p1.equal_addr(kp1)) {nt1->p3n = t1->p1n; nt1->p3p = t1->p1p;}
    else if (t1->p2.equal_addr(kp1)) {nt1->p3n = t1->p2n; nt1->p3p = t1->p2p;}
    else {nt1->p3n = t1->p3n; nt1->p3p = t1->p3p;}
    
    
    nt0->relink_neighbours();
    nt1->relink_neighbours();

    nt0->generate_facets(facets,curr_id);
    nt1->generate_facets(facets,curr_id);
}

void DT_SOLVER::flip23() {
    //! Hat das Facet  2,3,4  der Tetraeder  1,2,3,4   5,2,4,3  keine reflexe Kante, so
    //! wird ein  2->3 Flip  ausgefuehrt. Aus den beiden Tetraedern entstehen die drei
    //! neuen  5,1,2,4   5,1,4,3   5,1,3,2
    
    ++curr_generation;
    
    //! 1.) Aus den 2 Tetraedern des aktuellen Facets 3 neue erzeugen:
    stl_ptr<TETRAHEDRON> nt0(new TETRAHEDRON(cf->ps1,cf->p0,cf->ps2,cf->p2,curr_generation,curr_tet_id)); curr_tet_id++;
    stl_ptr<TETRAHEDRON> nt1(new TETRAHEDRON(cf->ps1,cf->p1,cf->ps2,cf->p0,curr_generation,curr_tet_id)); curr_tet_id++;
    stl_ptr<TETRAHEDRON> nt2(new TETRAHEDRON(cf->ps1,cf->p2,cf->ps2,cf->p1,curr_generation,curr_tet_id)); curr_tet_id++;
    
    //! 2.) Die neuen Tetraeder zufuegen:
    tetrahedrons.push_back(nt0);
    tetrahedrons.push_back(nt1);
    tetrahedrons.push_back(nt2);
    
    //! 3.) next und prev entsprechend fuer alle Tetraeder anpassen:
    cf->t1->next.push_back(nt0);
    cf->t1->next.push_back(nt1);
    cf->t1->next.push_back(nt2);
    cf->t2->next.push_back(nt0);
    cf->t2->next.push_back(nt1);
    cf->t2->next.push_back(nt2);
    nt0->prev.push_back(cf->t1);
    nt0->prev.push_back(cf->t2);
    nt1->prev.push_back(cf->t1);
    nt1->prev.push_back(cf->t2);
    nt2->prev.push_back(cf->t1);
    nt2->prev.push_back(cf->t2);
    
    
    nt0->p0n = cf->t2n1; nt0->p1n = nt2; nt0->p2n = cf->t1n1; nt0->p3n = nt1;
    nt0->p0p = cf->t2p1; nt0->p1p = cf->p1; nt0->p2p = cf->t1p1; nt0->p3p = cf->p1;
    
    nt1->p0n = cf->t2n2; nt1->p1n = nt0; nt1->p2n = cf->t1n2; nt1->p3n = nt2;
    nt1->p0p = cf->t2p2; nt1->p1p = cf->p2; nt1->p2p = cf->t1p2; nt1->p3p = cf->p2;
    
    nt2->p0n = cf->t2n0; nt2->p1n = nt1; nt2->p2n = cf->t1n0; nt2->p3n = nt0;
    nt2->p0p = cf->t2p0; nt2->p1p = cf->p0; nt2->p2p = cf->t1p0; nt2->p3p = cf->p0;
    
    
    nt0->relink_neighbours();
    nt1->relink_neighbours();
    nt2->relink_neighbours();
    
    
    nt0->generate_facets(facets,curr_id);
    nt1->generate_facets(facets,curr_id);
    nt2->generate_facets(facets,curr_id);    
}

void DT_SOLVER::make_regular() {
    //!Die DT durch eine Reihe von 2->3  und 3->2 Flips wieder herstellen:
    while (!facets.empty()) {
        
        if (facets.size() > 1000000) break;
        
        //! 1.) Facet vom Stapel holen (top und pop):
        cf = facets.top(); facets.pop();
        
        //! 2.) Regularitaet pruefen:
        //!     Ein Tetraeder ist regulaer, wenn kein anderer Punkt innerhalb des Tetraeders liegt. Es ist
        //!     ausreichend hierzu die Facets zu betrachten, die an alte Tetraeder grenzen:
        //!     (genaugenommen wird hier nur die Regularitaet in bezug auf dieses Link-Facet geprueft)
        
        if (!cf->is_regular()) {
            //! 3.) Auf reflexe Kanten pruefen:
            //!     Sei a,b,c ein Facet des Tetraeders a,b,c,d welches entsprechend an ein Tetraeder a*,b*,c*,d* angrenzt.
            //!     (a = a*     b = b*     c = c*)
            //!     Die Kante a,b ist reflex, wenn c und d* auf unterschiedlichen Seiten der durch a,b,d definierten Flaeche
            //!     liegen.
            //!     Wenn reflexe Kante, dann auf Flipbarkeit pruefen und gegebenenfalls 3->2 Flip (sonst das Facet verwerfen)
            //!     Wenn keine reflexe Kante, dann einen 2->3 Flip
            
            if (triple_product(cf->ps1->coord,cf->p1->coord,cf->p0->coord,cf->ps2->coord) > 0.) {
                //im Linksumlauf bezogen auf cf->t1:
                //p0 der reflexen Kante, p1 der reflexen Kante, dritter Punkt des Facets,
                //Nachbartetraeder an der reflexen Kante bezogen auf cf->t1, Nachbar bezogen auf cf->t2,
                //Nachbar zu p0 der reflexen Kante in Bezug auf cf->t1, ... in Bezug auf cf->t2,   #A
                //Nachbar zu p1 der reflexen Kante in Bezug auf cf->t1, ... in Bezug auf cf->t2,   #B
                //entsprechende Punkte gegenueber A
                //entsprechende Punkte gegenueber B
                cf->calc_n();
                flip32(cf->p1,cf->p0,cf->p2,
                       cf->t1n2,cf->t2n2,
                       cf->t1n1,cf->t2n1,
                       cf->t1n0,cf->t2n0,
                       cf->t1p1,cf->t2p1,
                       cf->t1p0,cf->t2p0); //prueft selbst, ob Flip moeglich
            } else if (triple_product(cf->ps1->coord,cf->p2->coord,cf->p1->coord,cf->ps2->coord) > 0.) {
                cf->calc_n();
                flip32(cf->p2,cf->p1,cf->p0,
                       cf->t1n0,cf->t2n0,
                       cf->t1n2,cf->t2n2,
                       cf->t1n1,cf->t2n1,
                       cf->t1p2,cf->t2p2,
                       cf->t1p1,cf->t2p1);
            } else if (triple_product(cf->ps1->coord,cf->p0->coord,cf->p2->coord,cf->ps2->coord) > 0.) {    
                cf->calc_n();
                flip32(cf->p0,cf->p2,cf->p1,
                       cf->t1n1,cf->t2n1,
                       cf->t1n0,cf->t2n0,
                       cf->t1n2,cf->t2n2,
                       cf->t1p0,cf->t2p0,
                       cf->t1p2,cf->t2p2);
            } else { //keine reflexe Kante
                cf->calc_n();
                flip23();
            }
        }
        
        cf.kill();
    }
    
    if (facets.size() > 1000000) {
        cerr << "error: DT_SOLVER::make_regular() -> too many facets on stack => triangulation will not be regular!" << endl;
        cerr << "                                    -> please make sure there are no 2 points with identical coordinates!" << endl;
        while (!facets.empty()) {
            cf = facets.top(); facets.pop();
            cf.kill();
        }
    }
}

void DT_SOLVER::add_dtpoint(stl_ptr<DT_POINT> const& point) {
    //!Fuegt der aktuellen Triangulierung einen Punkt hinzu und stellt die DT wieder her:
    
    //! 1.) Tetraeder lokalisieren in dem der neue Punkt liegt:
    #if defined (_SPEEDTEST)
    t_locate.continue_t();
    #endif
    curr_tet = locate_tetrahedron(point,tetrahedrons[0]);
    #if defined (_SPEEDTEST)
    t_locate.stop_t();
    #endif
    
    //! 2.) Einen 1->4 Flip ausfuehren:
    #if defined (_SPEEDTEST)
    t_flip.continue_t();
    #endif
    flip14(point); //wird immer fuer curr_tet ausgefuehrt
    #if defined (_SPEEDTEST)
    t_flip.stop_t();
    #endif
    
    //! 3.) wieder eine regulaere Triangulation herstellen:
    #if defined (_SPEEDTEST)
    t_regular.continue_t();
    #endif
    make_regular();
    #if defined (_SPEEDTEST)
    t_regular.stop_t();
    #endif
}

void DT_SOLVER::solve(bool check_regularity) {
    //!Stellt die DT (bzw. regulaere T.) ueber alle DT_POINTs her:
    //!Nach jedem eingefuegten Punkt wird zunaechst wieder eine regulaere Triangulierung hergestellt
    curr_generation = 0;
    
    #if defined (_SPEEDTEST)
    t_triang.start_t();
    #endif
    for (dtp_vec it=dt_points.begin(); it!=dt_points.end(); ++it) {
        add_dtpoint(*it);
    }
    #if defined (_SPEEDTEST)
    t_triang.stop_t();
    #endif
    
    #if defined (_SPEEDTEST)
    t_deltet.start_t();
    #endif
    for (th_vec it=tetrahedrons.begin(); it!=tetrahedrons.end(); ++it) {
        if ((*it)->next.size() > 0) continue;
        if ((*it)->p0.equal_addr(tetrahedrons[0]->p0) ||
            (*it)->p0.equal_addr(tetrahedrons[0]->p1) ||
            (*it)->p0.equal_addr(tetrahedrons[0]->p2) ||
            (*it)->p0.equal_addr(tetrahedrons[0]->p3) ||
            (*it)->p1.equal_addr(tetrahedrons[0]->p0) ||
            (*it)->p1.equal_addr(tetrahedrons[0]->p1) ||
            (*it)->p1.equal_addr(tetrahedrons[0]->p2) ||
            (*it)->p1.equal_addr(tetrahedrons[0]->p3) ||
            (*it)->p2.equal_addr(tetrahedrons[0]->p0) ||
            (*it)->p2.equal_addr(tetrahedrons[0]->p1) ||
            (*it)->p2.equal_addr(tetrahedrons[0]->p2) ||
            (*it)->p2.equal_addr(tetrahedrons[0]->p3) ||
            (*it)->p3.equal_addr(tetrahedrons[0]->p0) ||
            (*it)->p3.equal_addr(tetrahedrons[0]->p1) ||
            (*it)->p3.equal_addr(tetrahedrons[0]->p2) ||
            (*it)->p3.equal_addr(tetrahedrons[0]->p3)) {
            hull_tetrahedrons.push_back(*it);
        } else delaunay_tetrahedrons.push_back(*it);
    }
    #if defined (_SPEEDTEST)
    t_deltet.stop_t();
    #endif
    
    //!Ueberpruefen, ob das Ergenis wirklich eine regulaere Triangulierung ist:
//    cout << "checking all tetrahedrons of final triangulation..." << endl;
    #if defined (_SPEEDTEST)
    t_proof.start_t();
    #endif
    if (check_regularity) {
        //! kostet unheimlich viel Zeit, weil sinnloser Weise in is_regular alle Punkte geprueft werden
        //! -> mal nacharbeiten!!!!!!!!
        bool not_regular = false;
        for (th_vec it=delaunay_tetrahedrons.begin(); it!=delaunay_tetrahedrons.end(); ++it) {
            if (!(*it)->is_regular(dt_points)) {
                cerr << "error: tetrahedron " << (*it)->id << " is not regular!!!" << *it << endl;
                not_regular = true;
            }
        }
        if (not_regular) {
            cerr << "error: final triangulation is not regular!" << endl;
            cerr << "make sure you have no two points with the same coordinates!" << endl;
            //exit(1);
        }
    }
    #if defined (_SPEEDTEST)
    t_proof.stop_t();
    #endif
    
    #if defined (_SPEEDTEST)
    cout << "t_triang   = " << t_triang.get_time() << endl;
    cout << "t_locate   = " << t_locate.get_time() << endl;
    cout << "t_flip     = " << t_flip.get_time() << endl;
    cout << "t_regular  = " << t_regular.get_time() << endl;
    cout << "t_regtest  = " << t_regtest.get_time() << endl;
    cout << "t_deltet   = " << t_deltet.get_time() << endl;
    cout << "t_proof    = " << t_proof.get_time() << endl;
    #endif
}

bool DT_SOLVER::has_hull_point(stl_ptr<TETRAHEDRON> const& t) {
    if (hull_set.find(t->p0->id) != hull_set.end() ||
        hull_set.find(t->p1->id) != hull_set.end() ||
        hull_set.find(t->p2->id) != hull_set.end() ||
        hull_set.find(t->p3->id) != hull_set.end()) return true;
    else return false;
}

void DT_SOLVER::get_hull_set() {
    hull_set.clear();
    for (th_vec it=hull_tetrahedrons.begin(); it!=hull_tetrahedrons.end(); ++it) {
        hull_set.insert((*it)->p0->id);
        hull_set.insert((*it)->p1->id);
        hull_set.insert((*it)->p2->id);
        hull_set.insert((*it)->p3->id);
    }
}

void DT_SOLVER::get_hull_facets() {
    //!derzeit kann es doppelte hull_facets geben!!! => noch ein extra set<id> machen
    /*
    get_hull_points();
    hull_facets.clear();
    for (th_vec it=tetrahedrons.begin(); it!=tetrahedrons.end(); ++it) {
        if ((*it)->next.size() > 0) continue;
        set<int> sp;
        if (hull_set.find((*it)->p0->id) != hull_set.end()) sp.insert(0);
        if (hull_set.find((*it)->p1->id) != hull_set.end()) sp.insert(1);
        if (hull_set.find((*it)->p2->id) != hull_set.end()) sp.insert(2);
        if (hull_set.find((*it)->p3->id) != hull_set.end()) sp.insert(3);
        if (sp.size() != 3) continue;
        if (sp.find(0) == sp.end()) {
            hull_facets.push_back(new L_FACET((*it)->p1,(*it)->p2,(*it)->p3));
        } else if (sp.find(1) == sp.end()) {
            hull_facets.push_back(new L_FACET((*it)->p0,(*it)->p3,(*it)->p2));
        } else if (sp.find(2) == sp.end()) {
            hull_facets.push_back(new L_FACET((*it)->p0,(*it)->p1,(*it)->p3));
        } else if (sp.find(3) == sp.end()) {
            hull_facets.push_back(new L_FACET((*it)->p0,(*it)->p2,(*it)->p1));
        }
    }
    */
    //! obige Variante am 02.10.2009 ersetzt, weil es zu einem falschen hull_facet bei der 12A Tasche von 1phg kam:
    for (lf_vec lt=hull_facets.begin(); lt!=hull_facets.end(); ++lt) lt->kill();
    hull_facets.clear();
//    for (th_vec it=delaunay_tetrahedrons.begin(); it!=delaunay_tetrahedrons.end(); ++it) {
    for (th_vec it=hull_tetrahedrons.begin(); it!=hull_tetrahedrons.end(); ++it) {
        if ((*it)->p0.equal_addr(tetrahedrons[0]->p0) ||
            (*it)->p0.equal_addr(tetrahedrons[0]->p1) ||
            (*it)->p0.equal_addr(tetrahedrons[0]->p2) ||
            (*it)->p0.equal_addr(tetrahedrons[0]->p3)) {
            if ((*it)->p1.equal_addr(tetrahedrons[0]->p0) ||
                    (*it)->p1.equal_addr(tetrahedrons[0]->p1) ||
                    (*it)->p1.equal_addr(tetrahedrons[0]->p2) ||
                    (*it)->p1.equal_addr(tetrahedrons[0]->p3) ||
                    (*it)->p2.equal_addr(tetrahedrons[0]->p0) ||
                    (*it)->p2.equal_addr(tetrahedrons[0]->p1) ||
                    (*it)->p2.equal_addr(tetrahedrons[0]->p2) ||
                    (*it)->p2.equal_addr(tetrahedrons[0]->p3) ||
                    (*it)->p3.equal_addr(tetrahedrons[0]->p0) ||
                    (*it)->p3.equal_addr(tetrahedrons[0]->p1) ||
                    (*it)->p3.equal_addr(tetrahedrons[0]->p2) ||
                    (*it)->p3.equal_addr(tetrahedrons[0]->p3)) continue;
            hull_facets.push_back(new L_FACET((*it)->p1,(*it)->p2,(*it)->p3));
        } else if ((*it)->p1.equal_addr(tetrahedrons[0]->p0) ||
                   (*it)->p1.equal_addr(tetrahedrons[0]->p1) ||
                   (*it)->p1.equal_addr(tetrahedrons[0]->p2) ||
                   (*it)->p1.equal_addr(tetrahedrons[0]->p3)) {
            if ((*it)->p0.equal_addr(tetrahedrons[0]->p0) ||
                    (*it)->p0.equal_addr(tetrahedrons[0]->p1) ||
                    (*it)->p0.equal_addr(tetrahedrons[0]->p2) ||
                    (*it)->p0.equal_addr(tetrahedrons[0]->p3) ||
                    (*it)->p2.equal_addr(tetrahedrons[0]->p0) ||
                    (*it)->p2.equal_addr(tetrahedrons[0]->p1) ||
                    (*it)->p2.equal_addr(tetrahedrons[0]->p2) ||
                    (*it)->p2.equal_addr(tetrahedrons[0]->p3) ||
                    (*it)->p3.equal_addr(tetrahedrons[0]->p0) ||
                    (*it)->p3.equal_addr(tetrahedrons[0]->p1) ||
                    (*it)->p3.equal_addr(tetrahedrons[0]->p2) ||
                    (*it)->p3.equal_addr(tetrahedrons[0]->p3)) continue;
            hull_facets.push_back(new L_FACET((*it)->p0,(*it)->p3,(*it)->p2));
        } else if ((*it)->p2.equal_addr(tetrahedrons[0]->p0) ||
                   (*it)->p2.equal_addr(tetrahedrons[0]->p1) ||
                   (*it)->p2.equal_addr(tetrahedrons[0]->p2) ||
                   (*it)->p2.equal_addr(tetrahedrons[0]->p3)) {
            if ((*it)->p1.equal_addr(tetrahedrons[0]->p0) ||
                    (*it)->p1.equal_addr(tetrahedrons[0]->p1) ||
                    (*it)->p1.equal_addr(tetrahedrons[0]->p2) ||
                    (*it)->p1.equal_addr(tetrahedrons[0]->p3) ||
                    (*it)->p0.equal_addr(tetrahedrons[0]->p0) ||
                    (*it)->p0.equal_addr(tetrahedrons[0]->p1) ||
                    (*it)->p0.equal_addr(tetrahedrons[0]->p2) ||
                    (*it)->p0.equal_addr(tetrahedrons[0]->p3) ||
                    (*it)->p3.equal_addr(tetrahedrons[0]->p0) ||
                    (*it)->p3.equal_addr(tetrahedrons[0]->p1) ||
                    (*it)->p3.equal_addr(tetrahedrons[0]->p2) ||
                    (*it)->p3.equal_addr(tetrahedrons[0]->p3)) continue;
            hull_facets.push_back(new L_FACET((*it)->p0,(*it)->p1,(*it)->p3));
        } else if ((*it)->p3.equal_addr(tetrahedrons[0]->p0) ||
                   (*it)->p3.equal_addr(tetrahedrons[0]->p1) ||
                   (*it)->p3.equal_addr(tetrahedrons[0]->p2) ||
                   (*it)->p3.equal_addr(tetrahedrons[0]->p3)) {
            if ((*it)->p1.equal_addr(tetrahedrons[0]->p0) ||
                    (*it)->p1.equal_addr(tetrahedrons[0]->p1) ||
                    (*it)->p1.equal_addr(tetrahedrons[0]->p2) ||
                    (*it)->p1.equal_addr(tetrahedrons[0]->p3) ||
                    (*it)->p2.equal_addr(tetrahedrons[0]->p0) ||
                    (*it)->p2.equal_addr(tetrahedrons[0]->p1) ||
                    (*it)->p2.equal_addr(tetrahedrons[0]->p2) ||
                    (*it)->p2.equal_addr(tetrahedrons[0]->p3) ||
                    (*it)->p0.equal_addr(tetrahedrons[0]->p0) ||
                    (*it)->p0.equal_addr(tetrahedrons[0]->p1) ||
                    (*it)->p0.equal_addr(tetrahedrons[0]->p2) ||
                    (*it)->p0.equal_addr(tetrahedrons[0]->p3)) continue;
            hull_facets.push_back(new L_FACET((*it)->p0,(*it)->p2,(*it)->p1));
        }
    }
}

void DT_SOLVER::get_alpha_shape(vector<stl_ptr<L_FACET> >&afc,double alpha) {
    //! Ueber alle Tetraeder T der finalen Triangulierung DT gehen:
    //!    Wenn T->sphere_square_radius <= alpha*alpha: // T gehoert zum alpha-Komplex
    //!       Ueber alle Facets F von T:
    //!          Wenn fuer das an F angrenzende Tetraeder N gilt N->sphere_square_radius > alpha: // N gehoert nicht zum alpha-Komplex
    //!             F gehoert zum alpha-Shape
    double srad = alpha * alpha;
    tr1::unordered_set<int> visited;
    for (th_vec it=delaunay_tetrahedrons.begin(); it!=delaunay_tetrahedrons.end(); ++it) {
        //! ACHTUNG: Hier fehlt hier fehlen die Tetraeder, die mit dem allumfassenden Tetraeder gebildet werden.
        //!          Ist aber kein Problem, weil die ja generell nicht zum alpha-Komplex gehoeren koennen
        if ((*it)->sphere_square_radius > srad) continue;
        if ((*it)->p0n->sphere_square_radius > srad) {
                afc.push_back(new L_FACET((*it)->p1,(*it)->p2,(*it)->p3,
                                          (*it)->p0n,*it,
                                          (*it)->p0p,(*it)->p0,-1));
        }
        if ((*it)->p1n->sphere_square_radius > srad) {
                afc.push_back(new L_FACET((*it)->p0,(*it)->p3,(*it)->p2,
                                          (*it)->p1n,*it,
                                          (*it)->p1p,(*it)->p1,-1));
        }
        if ((*it)->p2n->sphere_square_radius > srad) {
                afc.push_back(new L_FACET((*it)->p0,(*it)->p1,(*it)->p3,
                                          (*it)->p2n,*it,
                                          (*it)->p2p,(*it)->p2,-1));
        }
        if ((*it)->p3n->sphere_square_radius > srad) {
                afc.push_back(new L_FACET((*it)->p0,(*it)->p2,(*it)->p1,
                                          (*it)->p3n,*it,
                                          (*it)->p3p,(*it)->p3,-1));
        }
    }
}

double DT_SOLVER::get_surface_area(vector<stl_ptr<L_FACET> >& fv) {
    double result = 0.;
    for (lf_vec lt=fv.begin(); lt!=fv.end(); ++lt) {
        result += get_triangle_area((*lt)->p0->coord,(*lt)->p1->coord,(*lt)->p2->coord);
    }
    return result;
}

double DT_SOLVER::get_exact_void_volume(vector<stl_ptr<TETRAHEDRON> >&tets) {
    //! Nach Inclusion-Exclusion Prinzip:
    double volume = 0.;
    //! 1.) Die Volumina der Tetraeder bestimmen:
    for (th_vec it=tets.begin(); it!=tets.end(); ++it) {
        volume += fabs(triple_product((*it)->p0->coord,(*it)->p1->coord,(*it)->p2->coord,(*it)->p3->coord)) / 6.;
    }
    
    
    //! 2.) Ueber alle Tetraeder gehen und die Sektoren vom Volumen abziehen (Schnitte der
    //!     Kugeln mit dem Tetraeder:
    for (th_vec it=tets.begin(); it!=tets.end(); ++it) {
        volume -= get_sphere_tetrahedron_intersection((*it)->p0->coord,(*it)->p1->coord,(*it)->p2->coord,(*it)->p3->coord,(*it)->p0->vol);
        volume -= get_sphere_tetrahedron_intersection((*it)->p1->coord,(*it)->p0->coord,(*it)->p2->coord,(*it)->p3->coord,(*it)->p1->vol);
        volume -= get_sphere_tetrahedron_intersection((*it)->p2->coord,(*it)->p0->coord,(*it)->p1->coord,(*it)->p3->coord,(*it)->p2->vol);
        volume -= get_sphere_tetrahedron_intersection((*it)->p3->coord,(*it)->p0->coord,(*it)->p1->coord,(*it)->p2->coord,(*it)->p3->vol);
    }
    
    
    //! 3.) Jetzt ueber alle gehen und die Schnittvolumen (je zweier Kugeln) wieder aufaddieren:
    double pi = acos(-1.);
    for (th_vec it=tets.begin(); it!=tets.end(); ++it) {
        double gamma = p2p_angle((*it)->p2->coord,(*it)->p1->coord,(*it)->p3->coord,(*it)->p0->coord,(*it)->p1->coord,(*it)->p3->coord) / pi;
        volume += gamma * get_sphere_overlap((*it)->p1->r,(*it)->p3->r,(*it)->p1->coord,(*it)->p3->coord);
        gamma = p2p_angle((*it)->p0->coord,(*it)->p1->coord,(*it)->p2->coord,(*it)->p0->coord,(*it)->p1->coord,(*it)->p3->coord) / pi;
        volume += gamma * get_sphere_overlap((*it)->p1->r,(*it)->p0->r,(*it)->p1->coord,(*it)->p0->coord);
        gamma = p2p_angle((*it)->p2->coord,(*it)->p1->coord,(*it)->p3->coord,(*it)->p2->coord,(*it)->p1->coord,(*it)->p0->coord) / pi;
        volume += gamma * get_sphere_overlap((*it)->p1->r,(*it)->p2->r,(*it)->p1->coord,(*it)->p2->coord);
        gamma = p2p_angle((*it)->p2->coord,(*it)->p0->coord,(*it)->p3->coord,(*it)->p2->coord,(*it)->p0->coord,(*it)->p1->coord) / pi;
        volume += gamma * get_sphere_overlap((*it)->p0->r,(*it)->p2->r,(*it)->p0->coord,(*it)->p2->coord);
        gamma = p2p_angle((*it)->p3->coord,(*it)->p0->coord,(*it)->p1->coord,(*it)->p3->coord,(*it)->p0->coord,(*it)->p2->coord) / pi;
        volume += gamma * get_sphere_overlap((*it)->p0->r,(*it)->p3->r,(*it)->p0->coord,(*it)->p3->coord);
        gamma = p2p_angle((*it)->p0->coord,(*it)->p2->coord,(*it)->p3->coord,(*it)->p1->coord,(*it)->p2->coord,(*it)->p3->coord) / pi;
        volume += gamma * get_sphere_overlap((*it)->p2->r,(*it)->p3->r,(*it)->p2->coord,(*it)->p3->coord);
    }
    
    
    //! 4.) Jetzt das Schnittvolumen je dreier Kugeln abziehen:
    for (th_vec it=tets.begin(); it!=tets.end(); ++it) {
        volume -= 0.5 * get_sphere3_overlap((*it)->p0->r,(*it)->p1->r,(*it)->p2->r,(*it)->p0->coord,(*it)->p1->coord,(*it)->p2->coord);
        volume -= 0.5 * get_sphere3_overlap((*it)->p0->r,(*it)->p1->r,(*it)->p3->r,(*it)->p0->coord,(*it)->p1->coord,(*it)->p3->coord);
        volume -= 0.5 * get_sphere3_overlap((*it)->p0->r,(*it)->p2->r,(*it)->p3->r,(*it)->p0->coord,(*it)->p2->coord,(*it)->p3->coord);
        volume -= 0.5 * get_sphere3_overlap((*it)->p1->r,(*it)->p2->r,(*it)->p3->r,(*it)->p1->coord,(*it)->p2->coord,(*it)->p3->coord);
    }
    
    return volume;
}

double DT_SOLVER::get_exact_occ_volume(vector<stl_ptr<TETRAHEDRON> >&tets) {
    //! Da hier auch Tetraeder mitgegeben werden koennen, die einen gemeinsamen Punkt
    //! mit dem allumfassenden Tet haben, muessen die Punkte darauf gechecked werden
    double volume = 0.;
    
    //! 1.) Volumina der einzelnen Kugeln addieren:
    tr1::unordered_set<int> visited;
    for (th_vec it=tets.begin(); it!=tets.end(); ++it) {
        if ((*it)->p0->id >= 0) {
            if (visited.find((*it)->p0->id) == visited.end()) {
                visited.insert((*it)->p0->id);
                volume += (*it)->p0->vol;
            }
        }
        if ((*it)->p1->id >= 0) {
            if (visited.find((*it)->p1->id) == visited.end()) {
                visited.insert((*it)->p1->id);
                volume += (*it)->p1->vol;
            }
        }
        if ((*it)->p2->id >= 0) {
            if (visited.find((*it)->p2->id) == visited.end()) {
                visited.insert((*it)->p2->id);
                volume += (*it)->p2->vol;
            }
        }
        if ((*it)->p3->id >= 0) {
            if (visited.find((*it)->p3->id) == visited.end()) {
                visited.insert((*it)->p3->id);
                volume += (*it)->p3->vol;
            }
        }
    }
    
//    cerr << "single spheres = " << volume << endl;
//    double deb1 = volume;
    
    //! 2.) 2er-Schnitte abziehen, wenn gemeinsame Kante.
    //!     Dabei darf jede Kante nur einmal beruecksichtigt werden!
    tr1::unordered_set<uint64_t> edge_hash;
    uint64_t curr_hash;
    for (th_vec it=tets.begin(); it!=tets.end(); ++it) {
        if ((*it)->p0->id >= 0) {
            if ((*it)->p1->id >= 0) {
                if ((*it)->p0->id < (*it)->p1->id) curr_hash = (*it)->p0->id + 10000000*(*it)->p1->id;
                else curr_hash = (*it)->p1->id + 10000000*(*it)->p0->id;
                if (edge_hash.find(curr_hash) == edge_hash.end()) {
                    volume -= get_sphere_overlap((*it)->p1->r,(*it)->p0->r,(*it)->p1->coord,(*it)->p0->coord);
                    edge_hash.insert(curr_hash);
                }
            }
            if ((*it)->p2->id >= 0) {
                if ((*it)->p0->id < (*it)->p2->id) curr_hash = (*it)->p0->id + 10000000*(*it)->p2->id;
                else curr_hash = (*it)->p2->id + 10000000*(*it)->p0->id;
                if (edge_hash.find(curr_hash) == edge_hash.end()) {
                    volume -= get_sphere_overlap((*it)->p0->r,(*it)->p2->r,(*it)->p0->coord,(*it)->p2->coord);
                    edge_hash.insert(curr_hash);
                }
            }
            if ((*it)->p3->id >= 0) {
                if ((*it)->p0->id < (*it)->p3->id) curr_hash = (*it)->p0->id + 10000000*(*it)->p3->id;
                else curr_hash = (*it)->p3->id + 10000000*(*it)->p0->id;
                if (edge_hash.find(curr_hash) == edge_hash.end()) {
                    volume -= get_sphere_overlap((*it)->p0->r,(*it)->p3->r,(*it)->p0->coord,(*it)->p3->coord);
                    edge_hash.insert(curr_hash);
                }
            }
        }
        if ((*it)->p1->id >= 0) {
            if ((*it)->p2->id >= 0) {
                if ((*it)->p1->id < (*it)->p2->id) curr_hash = (*it)->p1->id + 10000000*(*it)->p2->id;
                else curr_hash = (*it)->p2->id + 10000000*(*it)->p1->id;
                if (edge_hash.find(curr_hash) == edge_hash.end()) {
                    volume -= get_sphere_overlap((*it)->p1->r,(*it)->p2->r,(*it)->p1->coord,(*it)->p2->coord);
                    edge_hash.insert(curr_hash);
                }
            }
            if ((*it)->p3->id >= 0) {
                if ((*it)->p1->id < (*it)->p3->id) curr_hash = (*it)->p1->id + 10000000*(*it)->p3->id;
                else curr_hash = (*it)->p3->id + 10000000*(*it)->p1->id;
                if (edge_hash.find(curr_hash) == edge_hash.end()) {
                    volume -= get_sphere_overlap((*it)->p1->r,(*it)->p3->r,(*it)->p1->coord,(*it)->p3->coord);
                    edge_hash.insert(curr_hash);
                }
            }
        }
        if ((*it)->p2->id >= 0) {
            if ((*it)->p3->id >= 0) {
                if ((*it)->p2->id < (*it)->p3->id) curr_hash = (*it)->p2->id + 10000000*(*it)->p3->id;
                else curr_hash = (*it)->p3->id + 10000000*(*it)->p2->id;
                if (edge_hash.find(curr_hash) == edge_hash.end()) {
                    volume -= get_sphere_overlap((*it)->p2->r,(*it)->p3->r,(*it)->p2->coord,(*it)->p3->coord);
                    edge_hash.insert(curr_hash);
                }
            }
        }
    }
    
//    cerr << "double intersections = " << (deb1-volume) << endl;
//    deb1 = volume;
    
    //! 3.) 3er-Schnitte addieren, wenn gemeinsames Facet:
    //!     Dabei darf jedes Facet nur einmal beruecksichtigt werden!
    edge_hash.clear();
    uint64_t o1 = 2000000;
    uint64_t o2 = o1*o1;
    for (th_vec it=tets.begin(); it!=tets.end(); ++it) {
        if ((*it)->p0->id >= 0) {
            if ((*it)->p1->id >= 0 && (*it)->p2->id >= 0) {
                if ((*it)->p0->id < (*it)->p1->id) {
                    if ((*it)->p1->id < (*it)->p2->id) curr_hash = (*it)->p0->id + o1*(*it)->p1->id + o2*(*it)->p2->id;
                    else if ((*it)->p0->id < (*it)->p2->id) curr_hash = (*it)->p0->id + o1*(*it)->p2->id + o2*(*it)->p1->id;
                    else curr_hash = (*it)->p2->id + o1*(*it)->p0->id + o2*(*it)->p1->id;
                } else {
                    if ((*it)->p0->id < (*it)->p2->id) curr_hash = (*it)->p1->id + o1*(*it)->p0->id + o2*(*it)->p2->id;
                    else if ((*it)->p1->id < (*it)->p2->id) curr_hash = (*it)->p1->id + o1*(*it)->p2->id + o2*(*it)->p0->id;
                    else curr_hash = (*it)->p2->id + o1*(*it)->p1->id + o2*(*it)->p0->id;
                }
                if (edge_hash.find(curr_hash) == edge_hash.end()) {
                    volume += get_sphere3_overlap((*it)->p0->r,(*it)->p1->r,(*it)->p2->r,
                                                  (*it)->p0->coord,(*it)->p1->coord,(*it)->p2->coord);
                    edge_hash.insert(curr_hash);
                }
            }
            if ((*it)->p1->id >= 0 && (*it)->p3->id >= 0) {
                if ((*it)->p0->id < (*it)->p1->id) {
                    if ((*it)->p1->id < (*it)->p3->id) curr_hash = (*it)->p0->id + o1*(*it)->p1->id + o2*(*it)->p3->id;
                    else if ((*it)->p0->id < (*it)->p3->id) curr_hash = (*it)->p0->id + o1*(*it)->p3->id + o2*(*it)->p1->id;
                    else curr_hash = (*it)->p3->id + o1*(*it)->p0->id + o2*(*it)->p1->id;
                } else {
                    if ((*it)->p0->id < (*it)->p3->id) curr_hash = (*it)->p1->id + o1*(*it)->p0->id + o2*(*it)->p3->id;
                    else if ((*it)->p1->id < (*it)->p3->id) curr_hash = (*it)->p1->id + o1*(*it)->p3->id + o2*(*it)->p0->id;
                    else curr_hash = (*it)->p3->id + o1*(*it)->p1->id + o2*(*it)->p0->id;
                }
                if (edge_hash.find(curr_hash) == edge_hash.end()) {
                    volume += get_sphere3_overlap((*it)->p0->r,(*it)->p1->r,(*it)->p3->r,
                                                  (*it)->p0->coord,(*it)->p1->coord,(*it)->p3->coord);
                    edge_hash.insert(curr_hash);
                }
            }
            if ((*it)->p2->id >= 0 && (*it)->p3->id >= 0) {
                if ((*it)->p0->id < (*it)->p2->id) {
                    if ((*it)->p2->id < (*it)->p3->id) curr_hash = (*it)->p0->id + o1*(*it)->p2->id + o2*(*it)->p3->id;
                    else if ((*it)->p0->id < (*it)->p3->id) curr_hash = (*it)->p0->id + o1*(*it)->p3->id + o2*(*it)->p2->id;
                    else curr_hash = (*it)->p3->id + o1*(*it)->p0->id + o2*(*it)->p2->id;
                } else {
                    if ((*it)->p0->id < (*it)->p3->id) curr_hash = (*it)->p2->id + o1*(*it)->p0->id + o2*(*it)->p3->id;
                    else if ((*it)->p2->id < (*it)->p3->id) curr_hash = (*it)->p2->id + o1*(*it)->p3->id + o2*(*it)->p0->id;
                    else curr_hash = (*it)->p3->id + o1*(*it)->p2->id + o2*(*it)->p0->id;
                }
                if (edge_hash.find(curr_hash) == edge_hash.end()) {
                    volume += get_sphere3_overlap((*it)->p0->r,(*it)->p2->r,(*it)->p3->r,
                                                  (*it)->p0->coord,(*it)->p2->coord,(*it)->p3->coord);
                    edge_hash.insert(curr_hash);
                }
            }
        }
        if ((*it)->p1->id >= 0) {
            if ((*it)->p2->id >= 0 && (*it)->p3->id >= 0) {
                if ((*it)->p1->id < (*it)->p2->id) {
                    if ((*it)->p2->id < (*it)->p3->id) curr_hash = (*it)->p1->id + o1*(*it)->p2->id + o2*(*it)->p3->id;
                    else if ((*it)->p1->id < (*it)->p3->id) curr_hash = (*it)->p1->id + o1*(*it)->p3->id + o2*(*it)->p2->id;
                    else curr_hash = (*it)->p3->id + o1*(*it)->p1->id + o2*(*it)->p2->id;
                } else {
                    if ((*it)->p1->id < (*it)->p3->id) curr_hash = (*it)->p2->id + o1*(*it)->p1->id + o2*(*it)->p3->id;
                    else if ((*it)->p2->id < (*it)->p3->id) curr_hash = (*it)->p2->id + o1*(*it)->p3->id + o2*(*it)->p1->id;
                    else curr_hash = (*it)->p3->id + o1*(*it)->p2->id + o2*(*it)->p1->id;
                }
                if (edge_hash.find(curr_hash) == edge_hash.end()) {
                    volume += get_sphere3_overlap((*it)->p1->r,(*it)->p2->r,(*it)->p3->r,
                                                  (*it)->p1->coord,(*it)->p2->coord,(*it)->p3->coord);
                    edge_hash.insert(curr_hash);
                }
            }
        }
    }
    
//    cerr << "triple intersections = " << (volume-deb1) << endl;
//    deb1 = volume;
    
    //! 4.) 4er-Schnitte abziehen:
    for (th_vec it=tets.begin(); it!=tets.end(); ++it) {
        if ((*it)->p0->id < 0 || (*it)->p1->id < 0 || (*it)->p2->id < 0 || (*it)->p3->id < 0) continue;
        volume -= get_sphere4_overlap((*it)->p0->r,(*it)->p1->r,(*it)->p2->r,(*it)->p3->r,
                              (*it)->p0->coord,(*it)->p1->coord,(*it)->p2->coord,(*it)->p3->coord);
    }
    
//    cerr << "quad intersections = " << (deb1-volume) << endl;
    
    return volume;
}

void DT_SOLVER::get_tet_shape(vector<stl_ptr<L_FACET> >&afc,vector<stl_ptr<TETRAHEDRON> >&tets) {
    //! Die Facets, welche die Aussenhuelle von tets bilden in afc ablegen:
    tr1::unordered_set<int64_t> partof;
    for (th_vec it=tets.begin(); it!=tets.end(); ++it) partof.insert((*it)->id);
    for (th_vec it=tets.begin(); it!=tets.end(); ++it) {
        if (partof.find((*it)->p0n->id) == partof.end()) {
            afc.push_back(new L_FACET((*it)->p1,(*it)->p2,(*it)->p3));
        }
        if (partof.find((*it)->p1n->id) == partof.end()) {
            afc.push_back(new L_FACET((*it)->p0,(*it)->p3,(*it)->p2));
        }
        if (partof.find((*it)->p2n->id) == partof.end()) {
            afc.push_back(new L_FACET((*it)->p0,(*it)->p1,(*it)->p3));
        }
        if (partof.find((*it)->p3n->id) == partof.end()) {
            afc.push_back(new L_FACET((*it)->p0,(*it)->p2,(*it)->p1));
        }
    }
}

void DT_SOLVER::facets2stl(char const* fname,vector<stl_ptr<L_FACET> >& fv,vec3d<double> const& shift) {
    //! Ein STL File mit den Facets rausschreiben (veraltetes Format -> besser facets2wrl benutzen):
    ofstream fout;
    fout.open(fname);
    
    fout << "solid " << fname << "\n";
    for (lf_vec lt=fv.begin(); lt!=fv.end(); ++lt) {
        vec3d<double> normal((*lt)->p1->coord); normal -= (*lt)->p0->coord;
        vec3d<double> v2((*lt)->p2->coord); v2 -= (*lt)->p0->coord;
        normal *= v2; normal.norm();
        fout << "facet normal " << normal[0] << " " << normal[1] << " " << normal[2] << "\n";
        fout << "outer loop\n";
        fout << "vertex " << (*lt)->p0->coord[0]-shift[0] << " " << (*lt)->p0->coord[1]-shift[1] << " " << (*lt)->p0->coord[2]-shift[2] << "\n";
        fout << "vertex " << (*lt)->p1->coord[0]-shift[0] << " " << (*lt)->p1->coord[1]-shift[1] << " " << (*lt)->p1->coord[2]-shift[2] << "\n";
        fout << "vertex " << (*lt)->p2->coord[0]-shift[0] << " " << (*lt)->p2->coord[1]-shift[1] << " " << (*lt)->p2->coord[2]-shift[2] << "\n";
        fout << "endloop\n";
        fout << "endfacet\n";
    }
    fout << "endsolid " << fname << endl;
    
    fout.close();
}

void DT_SOLVER::facets2wrl(char const* fname,vector<stl_ptr<L_FACET> >& fv,vec3d<double> const& s) {
    //! Ein VRML97 File mit den Facets generieren:
    ofstream fout;
    fout.open(fname);
    
    fout << "Shape {\n";
    fout << " appearance Appearance {\n";
    fout << "  material Material { diffuseColor 1.0 1.0 0.0 }\n";
    fout << "}\n";
    fout << "geometry IndexedFaceSet {\n";
    fout << " coord Coordinate {\n";
    fout << "  point [\n";
    for (lf_vec lt=fv.begin(); lt!=fv.end(); ++lt) {
        fout << (*lt)->p0->coord[0]-s[0] << " " << (*lt)->p0->coord[1]-s[1] << " " << (*lt)->p0->coord[2]-s[2] << ",\n";
        fout << (*lt)->p1->coord[0]-s[0] << " " << (*lt)->p1->coord[1]-s[1] << " " << (*lt)->p1->coord[2]-s[2] << ",\n";
        fout << (*lt)->p2->coord[0]-s[0] << " " << (*lt)->p2->coord[1]-s[1] << " " << (*lt)->p2->coord[2]-s[2] << ",\n";
    }
    fout << "]\n}\n";
    fout << "coordIndex [\n";
    unsigned int idx = 0;
    for (lf_vec lt=fv.begin(); lt!=fv.end(); ++lt) {
        fout << idx << " ";
        ++idx;
        fout << idx << " ";
        ++idx;
        fout << idx << " -1,\n";
        ++idx;
    }
    fout << "]\n";
    fout << "normalPerVertex TRUE\n";
    fout << " normal Normal {\n";
    fout << "  vector [\n";
    for (lf_vec lt=fv.begin(); lt!=fv.end(); ++lt) {
        vec3d<double> normal((*lt)->p1->coord); normal -= (*lt)->p0->coord;
        vec3d<double> v2((*lt)->p2->coord); v2 -= (*lt)->p0->coord;
        normal *= v2; normal.norm();
        fout << normal[0] << " " << normal[1] << " " << normal[2] << ",\n";
    }
    fout << "]\n}\n";
    fout << "normalIndex [\n";
    idx = 0;
    for (lf_vec lt=fv.begin(); lt!=fv.end(); ++lt) {
        fout << idx << " ";
        ++idx;
        fout << idx << " ";
        ++idx;
        fout << idx << " -1,\n";
        ++idx;
    }
    fout << "]\n}\n}" << endl;
    
    fout.close();
}

void DT_SOLVER::facets2off(char const* fname,vector<stl_ptr<L_FACET> >& fv) {
    //! Ein OFF File fuer Mesh Viewer schreiben:
    ofstream fout;
    fout.open(fname);
    
    int idx = 0;
    tr1::unordered_map<int,int> id2id;
    for (lf_vec lt=fv.begin(); lt!=fv.end(); ++lt) {
        if (id2id.find((*lt)->p0->id) == id2id.end()) {
            id2id[(*lt)->p0->id] = idx; ++idx;
        }
        if (id2id.find((*lt)->p1->id) == id2id.end()) {
            id2id[(*lt)->p1->id] = idx; ++idx;
        }
        if (id2id.find((*lt)->p2->id) == id2id.end()) {
            id2id[(*lt)->p2->id] = idx; ++idx;
        }
    }
    fout << "OFF\n";
    fout << id2id.size() << " " << fv.size() << " 0\n";
    idx = 0;
    id2id.clear();
    for (lf_vec lt=fv.begin(); lt!=fv.end(); ++lt) {
        if (id2id.find((*lt)->p0->id) == id2id.end()) {
            fout << (*lt)->p0->coord[0] << " " << (*lt)->p0->coord[1] << " " << (*lt)->p0->coord[2] << "\n";
            id2id[(*lt)->p0->id] = idx; ++idx;
        }
        if (id2id.find((*lt)->p1->id) == id2id.end()) {
            fout << (*lt)->p1->coord[0] << " " << (*lt)->p1->coord[1] << " " << (*lt)->p1->coord[2] << "\n";
            id2id[(*lt)->p1->id] = idx; ++idx;
        }
        if (id2id.find((*lt)->p2->id) == id2id.end()) {
            fout << (*lt)->p2->coord[0] << " " << (*lt)->p2->coord[1] << " " << (*lt)->p2->coord[2] << "\n";
            id2id[(*lt)->p2->id] = idx; ++idx;
        }
    }
    for (lf_vec lt=fv.begin(); lt!=fv.end(); ++lt) {
        fout << 3 << " " << id2id[(*lt)->p0->id] << " ";
        fout << id2id[(*lt)->p1->id] << " ";
        fout << id2id[(*lt)->p2->id] << "\n";
    }
    
    fout.close();
}

void DT_SOLVER::get_volume_cluster(float const& alpha1,float const& alpha2,float const& min_vol,
                                   map<int,vector<stl_ptr<TETRAHEDRON> > > &final_tets,
                                   map<int,double> &final_vol) {
    vector<stl_ptr<L_FACET> > alpha_facets2;
    get_alpha_shape(alpha_facets2,alpha2);

    tr1::unordered_set<int> fset;
    for (lf_vec lt=alpha_facets2.begin(); lt!=alpha_facets2.end(); ++lt) {
        fset.insert((*lt)->p0->id);
        fset.insert((*lt)->p1->id);
        fset.insert((*lt)->p2->id);
    }

    vector<stl_ptr<TETRAHEDRON> > pt;
    tr1::unordered_map<int64_t,int> id2idx;
    int cindex = 0;
    double srad2 = alpha2 * alpha2;
    double srad1 = alpha1 * alpha1;
    tr1::unordered_set<int64_t> aforbidden;

    for (th_vec it=delaunay_tetrahedrons.begin(); it!=delaunay_tetrahedrons.end(); ++it) {
        id2idx[(*it)->id] = cindex; ++cindex;
        if ((*it)->sphere_square_radius > srad2) continue; // gehoert nicht zum alpha2-Komplex
        if ((*it)->sphere_square_radius <= srad1) continue; // gehoert zum alpha1-Komplex

        //! ACHTUNG: Das gewichtete Volumen der Tetraeder ist nur ein approximierter Wert
        //!          Die exakte Berechnung erfolgt spaeter fuer eine definierte Menge von
        //!          Tetraedern ueber calc_exact_volume
        (*it)->calc_weighted_volume();

        //! An dieser Stelle muss das Tetraeder zum alpha2-Komplex gehoeren, ist
        //! aber nicht Teil des alpha1-Komplex
        pt.push_back(*it);
        //! In forbidden packen, wenn gemeinsamer Punkt mit alpha2-Shape:
        if (fset.find((*it)->p0->id) != fset.end() ||
            fset.find((*it)->p1->id) != fset.end() ||
            fset.find((*it)->p2->id) != fset.end() ||
            fset.find((*it)->p3->id) != fset.end()) aforbidden.insert((*it)->id);
    }

    tr1::unordered_set<int64_t> avisited;
    list<ALPHACLUSTER*> clusters;
    for (th_vec it=pt.begin(); it!=pt.end(); ++it) {
        if (avisited.find((*it)->id) == avisited.end() && aforbidden.find((*it)->id) == aforbidden.end()) {
            avisited.insert((*it)->id);
            ALPHACLUSTER* clust = new ALPHACLUSTER();
            clust->elements.insert((*it)->id);
            clusters.push_back(clust);
            stack<stl_ptr<TETRAHEDRON> > tstack;
            tstack.push((*it)->p0n);
            tstack.push((*it)->p1n);
            tstack.push((*it)->p2n);
            tstack.push((*it)->p3n);
            while (tstack.size() > 0) {
                stl_ptr<TETRAHEDRON> t = tstack.top(); tstack.pop();
                if (avisited.find(t->id) != avisited.end()) {
                    if (clust->elements.find(t->id) == clust->elements.end()) {
                        for (list<ALPHACLUSTER*>::iterator jt=clusters.begin(); jt!=clusters.end(); ++jt) {
                            if (*jt == clust) continue;
                            if ((*jt)->elements.find(t->id) != (*jt)->elements.end()) {
                                for (tr1::unordered_set<int64_t>::iterator iit=(*jt)->elements.begin();
                                                                           iit!=(*jt)->elements.end(); ++iit) {
                                    clust->elements.insert(*iit);
                                }
                                delete *jt;
                                clusters.erase(jt);
                                break;
                            }
                        }
                    }
                } else {
                    if (t->sphere_square_radius <= srad2 && t->sphere_square_radius > srad1 &&
                        aforbidden.find(t->id) == aforbidden.end()) {
                        clust->elements.insert(t->id);
                        avisited.insert(t->id);
                        tstack.push(t->p0n);
                        tstack.push(t->p1n);
                        tstack.push(t->p2n);
                        tstack.push(t->p3n);
                    }
                }
            }
        }
    }

    int n_clust = 0;
    for (list<ALPHACLUSTER*>::iterator it=clusters.begin(); it!=clusters.end(); ++it) {
        double vol = 0.;
        for (tr1::unordered_set<int64_t>::iterator jt=(*it)->elements.begin(); jt!=(*it)->elements.end(); ++jt) {
            vol += delaunay_tetrahedrons[id2idx[*jt]]->volume;
        }
        if (vol < min_vol) continue;
        final_tets[n_clust] = vector<stl_ptr<TETRAHEDRON> >();
        for (tr1::unordered_set<int64_t>::iterator jt = (*it)->elements.begin(); jt!=(*it)->elements.end(); ++jt) {
            final_tets[n_clust].push_back(delaunay_tetrahedrons[id2idx[*jt]]);
        }
        final_vol[n_clust] = get_exact_void_volume(final_tets[n_clust]);
        ++n_clust;
    }
    for (list<ALPHACLUSTER*>::iterator it=clusters.begin(); it!=clusters.end(); ++it) delete *it;
}


void DT_SOLVER::get_voro_points(vector<stl_ptr<TETRAHEDRON> > &vt,vector<vec3d<float> > &vpoints,
                                float simplify_radius) {
    if (simplify_radius <= 0.) {
        for (th_vec it=vt.begin(); it!=vt.end(); ++it) {
            vpoints.push_back(vec3d<float>((*it)->sphere_middle));
        }
        return;
    }

    // Wenn der Sphere-Mittelpunkt eines Nachbartetraeders naeher als
    // simplify radius liegt werden die beiden zusammengefasst
    float srad = simplify_radius * simplify_radius;
    tr1::unordered_map<int64_t,int> id2idx;
    int cindex = 0;
    for (th_vec it=vt.begin(); it!=vt.end(); ++it) {
        id2idx[(*it)->id] = cindex; ++cindex;
    }

    tr1::unordered_set<int64_t> avisited;
    list<ALPHACLUSTER*> clusters;
    for (th_vec it=vt.begin(); it!=vt.end(); ++it) {
        if (avisited.find((*it)->id) == avisited.end()) {
            avisited.insert((*it)->id);
            ALPHACLUSTER* clust = new ALPHACLUSTER();
            clust->elements.insert((*it)->id);
            clusters.push_back(clust);
            stack<stl_ptr<TETRAHEDRON> > tstack;
            if (id2idx.find((*it)->p0n->id) != id2idx.end()) {
                if (get_square_distance((*it)->p0n->sphere_middle,(*it)->sphere_middle) < srad) tstack.push((*it)->p0n);
            }
            if (id2idx.find((*it)->p1n->id) != id2idx.end()) {
                if (get_square_distance((*it)->p1n->sphere_middle,(*it)->sphere_middle) < srad) tstack.push((*it)->p1n);
            }
            if (id2idx.find((*it)->p2n->id) != id2idx.end()) {
                if (get_square_distance((*it)->p2n->sphere_middle,(*it)->sphere_middle) < srad) tstack.push((*it)->p2n);
            }
            if (id2idx.find((*it)->p3n->id) != id2idx.end()) {
                if (get_square_distance((*it)->p3n->sphere_middle,(*it)->sphere_middle) < srad) tstack.push((*it)->p3n);
            }
            while (tstack.size() > 0) {
                stl_ptr<TETRAHEDRON> t = tstack.top(); tstack.pop();
                if (avisited.find(t->id) != avisited.end()) {
                    if (clust->elements.find(t->id) == clust->elements.end()) {
                        for (list<ALPHACLUSTER*>::iterator jt=clusters.begin(); jt!=clusters.end(); ++jt) {
                            if (*jt == clust) continue;
                            if ((*jt)->elements.find(t->id) != (*jt)->elements.end()) {
                                for (tr1::unordered_set<int64_t>::iterator iit=(*jt)->elements.begin();
                                                                           iit!=(*jt)->elements.end(); ++iit) {
                                    clust->elements.insert(*iit);
                                }
                                delete *jt;
                                clusters.erase(jt);
                                break;
                            }
                        }
                    }
                } else {
                    clust->elements.insert(t->id);
                    avisited.insert(t->id);
                    if (id2idx.find(t->p0n->id) != id2idx.end()) {
                        if (get_square_distance(t->p0n->sphere_middle,t->sphere_middle) < srad) tstack.push(t->p0n);
                    }
                    if (id2idx.find(t->p1n->id) != id2idx.end()) {
                        if (get_square_distance(t->p1n->sphere_middle,t->sphere_middle) < srad) tstack.push(t->p1n);
                    }
                    if (id2idx.find(t->p2n->id) != id2idx.end()) {
                        if (get_square_distance(t->p2n->sphere_middle,t->sphere_middle) < srad) tstack.push(t->p2n);
                    }
                    if (id2idx.find(t->p3n->id) != id2idx.end()) {
                        if (get_square_distance(t->p3n->sphere_middle,t->sphere_middle) < srad) tstack.push(t->p3n);
                    }
                }
            }
        }
    }

    for (list<ALPHACLUSTER*>::iterator it=clusters.begin(); it!=clusters.end(); ++it) {
        vec3d<float> cp(0.,0.,0.);
        for (tr1::unordered_set<int64_t>::iterator jt=(*it)->elements.begin(); jt!=(*it)->elements.end(); ++jt) {
            cp += vt[id2idx[*jt]]->sphere_middle;
        }
        cp /= (*it)->elements.size();
        vpoints.push_back(cp);
    }
    
    for (list<ALPHACLUSTER*>::iterator it=clusters.begin(); it!=clusters.end(); ++it) delete *it;
}
//!===============================================================================
