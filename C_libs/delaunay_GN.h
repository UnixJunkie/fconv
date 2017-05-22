
//============================================================================
// delaunay_GN.h -*- C++ -*-; generic implementation for regular triangulation
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


#ifndef __DELAUNAYGN
#define __DELAUNAYGN

#include<stdint.h>
#include<stdlib.h>
#include<iostream>
#include<fstream>
#include<sstream>
#include<vector>
#include<set>
#include<map>
#include<list>
#include<tr1/unordered_map>
#include<tr1/unordered_set>
#include<stack>
#include"linalg_GN.hpp"
#include"stl_ptr_GN.hpp"

using namespace std;

//!===============================================================================
//!Forward Deklarationen:
class DT_POINT;     //!Diese Klasse repraesentiert einen Punkt der DT
class TETRAHEDRON;  //!Ein Simplex der DT im R3
class L_FACET;      //!Repraesentiert ein Link-Facet (R3 => Dreieck)
class DT_SOLVER;    //!Diese Klasse errechnet die DT und ist das Interface fuer den user
//!===============================================================================


//!===============================================================================
//!Typen und Konstanten:
typedef stl_ptr<DT_POINT> dtp_pointer;                  //!Zeiger auf einen DT_POINT
typedef stl_ptr<TETRAHEDRON> th_pointer;                //!Zeiger auf ein TETRAHEDRON
typedef stl_ptr<L_FACET> lf_pointer;                    //!Zeiger auf ein L_FACET
typedef vector<stl_ptr<DT_POINT> > dtp_container;       //!vector fuer dtp_pointer
typedef vector<stl_ptr<TETRAHEDRON> > th_container;     //!vector fuer th_pointer
typedef vector<stl_ptr<L_FACET> > lf_container;         //!vector fuer lf_pointer
typedef vector<stl_ptr<DT_POINT> >::iterator dtp_vec;   //!Iterator ueber DT_POINTs
typedef vector<stl_ptr<TETRAHEDRON> >::iterator th_vec; //!Iterator ueber TETRAHEDRONs
typedef vector<stl_ptr<L_FACET> >::iterator lf_vec;     //!Iterator ueber L_FACETs
//!===============================================================================

//!===============================================================================
//!Deklaration der Klasse DT_POINT:
class DT_POINT {
    public:
        int id;
        vec3d<double> coord;
        double w; //! Gewichtung des Punktes (handelt es sich um eine Kugel ist r^2 das Gewicht)
        double r; //! Radius ( == sqrt(w) )
        double vol; //! Volumen der Kugel
        void* myobj;
        
        DT_POINT(int const& i,double const& x,double const& y,double const& z,double const& weight = 0.,
                 double const& radius = 0.,void *at = NULL);
        template<class T> DT_POINT(int const& i,vec3d<T> const& pos,double const& weight = 0.,
                                   double const& radius = 0.,void *at = NULL):id(i),
                                   coord(pos),w(weight),r(radius),myobj(at) {calc_vol();}
        DT_POINT(DT_POINT const& ref);
        ~DT_POINT();
        
        void calc_vol();
};
//!===============================================================================

//!===============================================================================
//!Deklaration der Klasse TETRAHEDRON:
class  TETRAHEDRON{
    public:
        stl_ptr<DT_POINT> p0;
        stl_ptr<DT_POINT> p1;
        stl_ptr<DT_POINT> p2;
        stl_ptr<DT_POINT> p3;
        
        stl_ptr<TETRAHEDRON> p0n; //!Tetraeder welches an die Flaeche p1,p2,p3 grenzt
        stl_ptr<TETRAHEDRON> p1n;
        stl_ptr<TETRAHEDRON> p2n;
        stl_ptr<TETRAHEDRON> p3n;
        
        stl_ptr<DT_POINT> p0p; //!nicht gemeinsamer Punkt von p0n
        stl_ptr<DT_POINT> p1p;
        stl_ptr<DT_POINT> p2p;
        stl_ptr<DT_POINT> p3p;
        
        vector<stl_ptr<TETRAHEDRON> > prev; //!prev und next verweisen auf die entsprechenden Vorgaenger
        vector<stl_ptr<TETRAHEDRON> > next; //!und Nachfolger ==> sie ersetzen den His_DAG
                                            //!==> alte Tetraeder duerfen nicht geloescht werden!!!
                                            //!==> Tetraeder mit next.empty() gehoeren zur aktuellen DT
        
        vec3d<double> sphere_middle; //!Mittelpunkt der Umkugel (bzw. orthogonalen Kugel bei Gewichten != 0)
        double sphere_square_radius; //!Quadrat des Radius der Umkugel
        double volume;               //!Volumen des Tetraeders (steht erst nach calc_volume() zur Verfuegung!)
        
        int64_t generation;
        
        int64_t id;
        
        TETRAHEDRON(stl_ptr<DT_POINT> const& v0,stl_ptr<DT_POINT> const& v1,
                    stl_ptr<DT_POINT> const& v2,stl_ptr<DT_POINT> const& v3,
                    int64_t const& currgen,int64_t const& i);
        TETRAHEDRON(DT_POINT& v0,DT_POINT& v1,
                    DT_POINT& v2,DT_POINT& v3,
                    int64_t const& currgen,int64_t const& i);
        TETRAHEDRON(TETRAHEDRON const& ref);
        ~TETRAHEDRON();
        
        inline bool operator<(TETRAHEDRON const& rechts);
        
        void visualize(int const& state,bool const& vdw_rad = false);
        void visualize_sphere(string const& sname);
        
        inline bool is_regular(vector<stl_ptr<DT_POINT> > &dt_points);
        inline bool coplanar();
        bool shares_point(stl_ptr<TETRAHEDRON> const& tet);
        bool has_infinit_point();
        bool shares_facets(stl_ptr<TETRAHEDRON> const& tet);
        inline bool calc_sphere(); //!Die Umkugel (bzw. orthogonale Kugel) berechnen
        void calc_volume(); //!Das Volumen des Tetraeders berechnen
        void calc_weighted_volume(); //!Das Volumen des Tetraeders berechnen und die vdW-Volumina abziehen
        void calc_exact_volume(double const& alpha); //!
        void reduce_sphere(); //!Die orthogonale Kugel verkleinern, so dass sie komplett innerhalb der vdW-Kugeln liegt
        inline bool is_positive_oriented(); //!bilden p0p1  p0p2  p0p3  ein Rechtssystem?
        inline bool contains_p(stl_ptr<DT_POINT> const& point); //!Ist point im Tetraeder?
        template<class T> inline bool contains_p(vec3d<T> const& point) {
            //!Wenn die Determinante Null ist liegt der neue Punkt exakt auf einer Tetraederflaeche.
            //!In diesem Fall wird er willkuerlich dem aktuellen Tetraeder zugeordnet, was zu einem
            //!unechten (4 koplanare Punkte) Tetraeder fuehrt. Die folgenden Flips sollten das Problem
            //!wieder beheben:
            if (triple_product(p0->coord,p1->coord,p2->coord,point) < 0.) return false;
            if (triple_product(p0->coord,p3->coord,p1->coord,point) < 0.) return false;
            if (triple_product(p0->coord,p2->coord,p3->coord,point) < 0.) return false;
            if (triple_product(p1->coord,p3->coord,p2->coord,point) < 0.) return false;
            return true;
        }
        inline void check_link(stl_ptr<TETRAHEDRON> const& t);
        inline void relink_neighbours();
        inline void generate_facets(stack<stl_ptr<L_FACET> > &facets,int64_t curr_id);
};
//!===============================================================================


//!===============================================================================
//!Deklaration der Klasse L_FACET:
class  L_FACET{
    public:
        stl_ptr<DT_POINT> p0; //!Die 3 Punkte des Facets
        stl_ptr<DT_POINT> p1;
        stl_ptr<DT_POINT> p2;
        
        stl_ptr<TETRAHEDRON> t1; //!Die beiden beteiligten Tetraeder
        stl_ptr<TETRAHEDRON> t2;
        
        stl_ptr<DT_POINT> ps1; //!Die jeweiligen Punkte von t1 und t2 die
        stl_ptr<DT_POINT> ps2; //!nicht am Link-Facet beteiligt sind
        
        stl_ptr<TETRAHEDRON> t1n0; //!Tetraeder gegenueber p0 auf t1 bezogen
        stl_ptr<TETRAHEDRON> t1n1;
        stl_ptr<TETRAHEDRON> t1n2;
        stl_ptr<TETRAHEDRON> t2n0; //!Tetraeder gegenueber p0 auf t2 bezogen
        stl_ptr<TETRAHEDRON> t2n1;
        stl_ptr<TETRAHEDRON> t2n2;
        
        stl_ptr<DT_POINT> t1p0; //!nicht gemeinsamer Punkt von t1n0
        stl_ptr<DT_POINT> t1p1;
        stl_ptr<DT_POINT> t1p2;
        stl_ptr<DT_POINT> t2p0; //!nicht gemeinsamer Punkt von t2n0
        stl_ptr<DT_POINT> t2p1;
        stl_ptr<DT_POINT> t2p2;
        
        int id;
        
        L_FACET();
        L_FACET(stl_ptr<DT_POINT> const& vp0,stl_ptr<DT_POINT> const& vp1,stl_ptr<DT_POINT> const& vp2);
        L_FACET(stl_ptr<DT_POINT> const& vp0,stl_ptr<DT_POINT> const& vp1,stl_ptr<DT_POINT> const& vp2,
                stl_ptr<TETRAHEDRON> const& vt1,stl_ptr<TETRAHEDRON> const& vt2,
                stl_ptr<DT_POINT> const& vps1,stl_ptr<DT_POINT> const& vps2,int64_t const& i);
        ~L_FACET();
        
        inline L_FACET& operator=(L_FACET const& facet);
        
        inline bool is_regular();
        inline void calc_n();
        
        template<class T> bool is_inside(vec3d<T> const& point) {
            if (triple_product(p0->coord,p1->coord,p2->coord,point) > 0.) return true;
            else return false;
        }
        
        void visualize();
};
//!===============================================================================


struct ALPHACLUSTER {
    tr1::unordered_set<int64_t> elements;
};


//!===============================================================================
//!Deklaration der Klasse DT_SOLVER:
class DT_SOLVER {
    private:
        
        vector<stl_ptr<DT_POINT> > dt_points; //!alle Punkte, die in die Triangulierung einbezogen werden sollen
        
        stl_ptr<TETRAHEDRON> zero; //Hilfsvariable
        
        stl_ptr<TETRAHEDRON> curr_tet; //!aktuelles Tetraeder
        
        stack<stl_ptr<L_FACET> > facets; //!abzuarbeitender Stack mit Link-Facets
        
        stl_ptr<L_FACET> cf; //!aktuelles Link-Facet
        
        int64_t curr_generation; //!aktuelle Generation von Punkten -> zum ermitteln, ob ein Facet an ein altes
                                  //!Tetraeder grenzt und somit ein Link-Facet ist
        
        int64_t curr_id; //!laufende ID fuer Facets
        int64_t curr_tet_id; //!laufende ID fuer Tetraeder
        
        template<class T> inline void points2dt_points(vector<vec3d<T> > &points) { //!vec3d's in DT_POINTs umwandeln
            //!ACHTUNG: Hier bekommen alle Punkte das Gewicht Null !!!
            int pid = 0;
            for (typename vector<vec3d<T> >::iterator it=points.begin(); it!=points.end(); ++it) {
                stl_ptr<DT_POINT> np(new DT_POINT(pid,*it));
                ++pid;
                dt_points.push_back(np);
            }
        }
        inline void get_initial_tetrahedron(); //!das allumfassende Tetraeder erzeugen
        
        stl_ptr<TETRAHEDRON>& locate_tetrahedron(stl_ptr<DT_POINT> const& point,
                                                 stl_ptr<TETRAHEDRON> &root);
        template<class T> stl_ptr<TETRAHEDRON>& locate_tetrahedron(vec3d<T> const& point,stl_ptr<TETRAHEDRON>& root) {
            //!Rekursive Funktion um zu ermitteln in welchem Tetraeder point liegt
            //!=> mit tetrahedrons[0] aufrufen
            if (root->next.empty()) return root; //!Das Tetraeder wurde gefunden
            for (th_vec it=root->next.begin(); it!=root->next.end(); ++it) {
                if ((*it)->contains_p(point)) return locate_tetrahedron(point,*it);
            }
            
            cerr << "error: point localization failed for " << point << " !" << endl;
            cerr << "       trying hardcore method..." << endl;
            
            //!Wurde kein Tetraeder gefunden, so ist entweder das allumfassende Tetraeder zu klein gewaehlt oder
            //!es liegt ein Problem mit der numerischen Genauigkeit vor.
            //!Jetzt die brutale Methode:  Alle aktuellen Tetraeder durchgehen und gucken in welchem es am ehesten liegt:
            for (th_vec it=tetrahedrons.begin(); it!=tetrahedrons.end(); ++it) {
                if ((*it)->next.size() > 0) continue;
                bool hit = true;
                if (triple_product((*it)->p0->coord,(*it)->p1->coord,(*it)->p2->coord,point) <= 0.) hit = false;
                if (triple_product((*it)->p0->coord,(*it)->p3->coord,(*it)->p1->coord,point) <= 0.) hit = false;
                if (triple_product((*it)->p0->coord,(*it)->p2->coord,(*it)->p3->coord,point) <= 0.) hit = false;
                if (triple_product((*it)->p1->coord,(*it)->p3->coord,(*it)->p2->coord,point) <= 0.) hit = false;
                
                if (hit) return *it;
            }
            cerr << "error: point localization failed again!" << endl;
            exit(1);
            return zero;
        }
        inline void flip14(stl_ptr<DT_POINT> const& point); //nutzt curr_tet
        inline void flip32(stl_ptr<DT_POINT> const& kp1,stl_ptr<DT_POINT> const& kp2,stl_ptr<DT_POINT> const& kp3,
                           stl_ptr<TETRAHEDRON> const& t1,stl_ptr<TETRAHEDRON> const& t2,
                           stl_ptr<TETRAHEDRON> const& lt1n0,stl_ptr<TETRAHEDRON> const& lt2n0,
                           stl_ptr<TETRAHEDRON> const& lt1n1,stl_ptr<TETRAHEDRON> const& lt2n1,
                           stl_ptr<DT_POINT> const& lt1p0,stl_ptr<DT_POINT> const& lt2p0,
                           stl_ptr<DT_POINT> const& lt1p1,stl_ptr<DT_POINT> const& lt2p1); //nutzt cf
        inline void flip23();
        
        inline void make_regular();
    public: //!Es folgt das Interface dieser Bibliothek:
        //!Die Initialisierung kann mit vec3d-Objekten oder mit DT_POINTs erfolgen
        //!Sollen die Punkte gewichtet werden ist es ZWINGEND notwendig direkt mit
        //!gewichteten DT_POINTs zu Initialisieren!!! ( DT_POINT(double x,double y, double z, double weight) )
        vector<stl_ptr<TETRAHEDRON> > tetrahedrons; //!alle Tetraeder des HisDAG  //war frueher private
        vector<stl_ptr<TETRAHEDRON> > delaunay_tetrahedrons; //!Tetraeder der finalen Triangulierung OHNE allumfassendes Tetraeder
        vector<stl_ptr<TETRAHEDRON> > hull_tetrahedrons; //! NUR Tetraeder, die sich MINDESTENS einen Punkt mit dem allumfassenden teilen
        vector<stl_ptr<DT_POINT> > hull_points; //!get_hull_points() schreibt hier alle Punkte rein, die zur convexen Huelle gehoeren
        vector<stl_ptr<L_FACET> > hull_facets; //!get_hull_facets() legt hier alle Facets ab aus denen die convexe Huelle besteht
        
        tr1::unordered_set<int> hull_set; //! IDs der Punkte, welche zur convexen Huelle gehoeren
        
        vector<stl_ptr<L_FACET> > alpha_facets; //!get_alpha_hull(double alpha) legt hier die Facets der alpha-shape ab
        vector<vec3d<double> > alpha_spheres;
        
        template<class T> DT_SOLVER(vector<vec3d<T> > &points) {
            //!ACHTUNG: Hier bekommen alle Punkte das Gewicht Null !!!
            points2dt_points(points);
            get_initial_tetrahedron();
            curr_id = 0;
            curr_tet_id = 1;
        }
        DT_SOLVER(vector<stl_ptr<DT_POINT> > &dt_pnt);
        DT_SOLVER(vector<DT_POINT> &dt_pnt);
        ~DT_SOLVER();
        
        void visualize_tetrahedrons();
        void visualize_tetrahedrons(vector<stl_ptr<TETRAHEDRON> > &tets,string fname = "X",string pname = "X");
        void visualize_convex_hull(string fname = "X");
        void visualize_alpha_spheres_new(string fname = "X",double alpha = 3.);
        void visualize_alpha_facets(string fname = "X");
        void visualize_alpha_spheres(); //! Das sind die falschen alpha_spheres!!!
        static void visualize_spheres(vector<stl_ptr<TETRAHEDRON> > &tets);
        static void visualize_facets(vector<stl_ptr<L_FACET> > &fcs,string fname = "X");
        static void visualize_sphere_maps(map<int,vector<stl_ptr<TETRAHEDRON> > > &cluster,string fname = "X");
        void add_dtpoint(stl_ptr<DT_POINT> const& point); //!Der aktuellen DT einen Punkt zufuegen und die DT wieder herstellen
        void solve(bool check_regularity = false); //!Fuer alle Punkte die DT herstellen
        
        bool has_hull_point(stl_ptr<TETRAHEDRON> const& t);
        void get_hull_set();
        void get_hull_facets();
        
        void get_alpha_shape(vector<stl_ptr<L_FACET> >& afc,double alpha);
        void get_tet_shape(vector<stl_ptr<L_FACET> >& afc,vector<stl_ptr<TETRAHEDRON> >& tets);
        void facets2stl(char const* fname,vector<stl_ptr<L_FACET> >& fv,vec3d<double> const& shift);
        void facets2wrl(char const* fname,vector<stl_ptr<L_FACET> >& fv,vec3d<double> const& s);
        void facets2off(char const* fname,vector<stl_ptr<L_FACET> >& fv);
        double get_surface_area(vector<stl_ptr<L_FACET> >& fv);
        double get_exact_void_volume(vector<stl_ptr<TETRAHEDRON> >&tets);
        double get_exact_occ_volume(vector<stl_ptr<TETRAHEDRON> >&tets);

        void get_volume_cluster(float const& alpha1,float const& alpha2,float const& min_vol,
                                map<int,vector<stl_ptr<TETRAHEDRON> > > &final_tets,
                                map<int,double> &final_vol);
        void get_voro_points(vector<stl_ptr<TETRAHEDRON> > &vt,vector<vec3d<float> > &vpoints,
                             float simplify_radius = -1.);
        
        template<class T> inline bool is_in_hull(vec3d<T> const& p) {
            /*
            curr_tet = locate_tetrahedron(p,tetrahedrons[0]);
            if (curr_tet->shares_point(tetrahedrons[0])) return false;
            return true;
            */
            
            //! 02.10.09:
            if (hull_facets.size() == 0) {
                cerr << "warning: using slow function for is_in_hull() -- call get_hull_facets() before to enable fast version" 
                     << endl;
                curr_tet = locate_tetrahedron(p,tetrahedrons[0]);
                if (curr_tet->shares_point(tetrahedrons[0])) return false;
                return true;
            } else {
                for (lf_vec lt=hull_facets.begin(); lt!=hull_facets.end(); ++lt) {
                    if (triple_product((*lt)->p0->coord,(*lt)->p1->coord,(*lt)->p2->coord,p) < 0.) return false;
                }
                return true;
            }
        }
};

ostream &operator<<(ostream &os,stl_ptr<TETRAHEDRON> const& t);
ostream &operator<<(ostream &os,stl_ptr<L_FACET> const& f);
//!===============================================================================

#endif
