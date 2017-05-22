
//============================================================================
// structure_GN.h -*- C++ -*-; high level representation of molecular data
//
// Copyright (C) 2006, 2007, 2008, 2009, 2010, 2011 Gerd Neudert
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
// This library implements structure objects as high level abstraction of
// molecular data.
//============================================================================


#ifndef __STRUCTUREGN
#define __STRUCTUREGN


#include<stdlib.h>
#include<iostream>
#include<string>
#include<fstream>
#include<sstream>
#include<algorithm>
#include<vector>
#include<map>
#include<tr1/unordered_map>
#include<set>
#include<tr1/unordered_set>
#include<queue>
#include<stack>
#include<list>

#include"linalg_GN.hpp"
#include"stl_ptr_GN.hpp"
#include"grid_GN.hpp"
#include"message_GN.hpp"
#include"optalign_GN.hpp"
#include"seqalign_GN.hpp"
#include"bron_kerbosch_GN.hpp"
#include"string_fu_GN.hpp"
#include"delaunay_GN.h"

#include"atom_properties_GN.h"
#include"atom_GN.h"
#include"molecule_GN.h"
#include"protein_GN.h"
#include"structure_additional_GN.h"


using namespace std;
using namespace TXT;


//==================================================================================================
//Forward-Deklarationen:
//==================================================================================================

class ATOM;        //atom_GN.h
class ATOM_EXT;        //atom_GN.h
class COMPARENODE;    //atom_GN.h
class BOND;        //atom_GN.h
class RING;        //atom_GN.h
class COMMENT;        //atom_GN.h
class PDBCONECT;    //atom_GN.h
class MOLECULE;        //molecule_GN.h
class NODE;        //atom_GN.h
class LIGAND;        //molecule_GN.h
class HET;        //atom_GN.h
class WATER;        //atom_GN.h
class METAL;        //atom_GN.h
class AACID;        //molecule_GN.h
class CHAIN;        //molecule_GN.h
class PROTEIN;        //
class CAVITY;        //molecule_GN.h
class CRYSIN_POSITION;    //structure_additional_GN.h
class CRYSIN;        //structure_additional_GN.h
class HOLE;        //
class CAVDUM;        //
class RELIBASE_CENTER;    //molecule_GN.h
class STRUCTURE;    //hier

class CELL;        //
class STRUCTURE_GRID;    //

class COORD_GRID;    //

class CRYSTALIZER;    //

//==================================================================================================
//Typedefs:
//==================================================================================================

typedef vector< stl_ptr<HOLE> >::iterator holes_vec;
typedef map<int,vector< stl_ptr<HOLE> > >::iterator holes_map;
typedef vector< stl_ptr<CAVDUM> >::iterator cavdums_vec;
typedef vector< stl_ptr<STRUCTURE> >::iterator alt_struct_vec;
typedef vector< stl_ptr<CELL> >::iterator cells_vec;


//==================================================================================================
//Deklaration der Klasse HOLE:
//==================================================================================================

class HOLE {
    public:
        vector<stl_ptr<CELL> > cells;
        HOLE();
        ~HOLE();
};


//==================================================================================================
//Deklaration der Klasse CAVDUM:
//==================================================================================================

class CAVDUM {
    public:
        float total_volume;
        float deepness;
        map<int,vector<stl_ptr<HOLE> > > holes;  //! schale:holes
        CAVDUM();
        ~CAVDUM();
        bool operator<(CAVDUM &rechts);
        void get_volume();
        void pymol_vis(int j,float spacing);
};


//==================================================================================================
//Deklaration der Klasse STRUCTURE:
//==================================================================================================

class STRUCTURE {
    private:
        tr1::unordered_set<int64_t> aforbidden;
        tr1::unordered_set<int64_t> avisited;
    public:
        string name;
        bool has_peptides;
        vector< stl_ptr<PDBCONECT> > conects;
        vector< stl_ptr<LIGAND> > ligands;
        vector< stl_ptr<LIGAND> > splitted_ligands;
        vector< stl_ptr<LIGAND> > mol2_protein; //siehe full_pdb2mol2()
        vector< stl_ptr<HET> > het_entrys;
        vector< stl_ptr<WATER> > waters;
        vector< stl_ptr<METAL> > metals;
        stl_ptr<LIGAND> water_ligand;
        stl_ptr<LIGAND> metal_ligand;
        PROTEIN* protein;
        vector< stl_ptr<CAVITY> > cavities;
        int verb;
        vector< stl_ptr<STRUCTURE> > alt_struct; //zeigt auf alternative Strukturen (z.B. fr NMR-pdb's)
        int model_number; //fuer alt_structs
        STRUCTURE(int verbosity = 1);
        ~STRUCTURE(); //muss conects, ligands, het_entrys, waters, metals, protein, alt_struct loeschen
        void clear(); //loescht alle Daten
        void clear_ligands(); //loescht nur ligands, waters und metals
        bool pdb2mol2(int mode,const char *def_file = "X",bool fill_X = false,bool bond_building = true); //pdb-Protein in mol2-Format
                                                                                //umwandeln (kann dann mit 
                                                                                //write_mol2_protein vom Parser geschrieben werden)
        void full_pdb2mol2(int mode = 1,bool get_ele = true,const char *def_file = "X",bool get_bonds = true,int verb = 1,
                           bool kill_ext = true,bool fill_X = false,const char *alt_def_file = "X",int max_ring_members = 10,
                           bool no_free_rot_planar = true,LIGAND* at_ref_mol = 0,LIGAND* at_ref_mol2 = 0,bool prot_acids = false,
                           bool prot_guanidin = true,bool prot_amidin = true,bool prot_amin = false,bool prot_phosphate = false,
                           bool prot_sulfate = false); //benutzt das sybyl atom typing fuer das ganze Protein samt Liganden
                                                       //und schreibt das Protein und die Liganden in mol2_protein
        void get_lig_cavity(float dist,LIGAND* lig,bool full_res=true, //erzeugt ein Cavity-Objekt unter Bercksichtigung des
                            bool no_hyd=false);                        //angegebenen Abstandes zum Liganden  -- Wird kein Ligand
                                                                       //mitgegeben, so wird der groesste Ligand genommen
        void get_buriedness();
        void get_cavity_auto(float radius = 6.,int min_buried = 18,bool debug_mode = false,unsigned int vis_number = 1); 
                                                                                                                      //Versucht automatisch
                                                                                                                      //Cavities zu finden
        void get_cavdum_cavity(float dist, stl_ptr<CAVDUM> &cavdum);
        void get_alpha_cavity(vector<stl_ptr<ATOM> >& pv,map<int,vector<stl_ptr<ATOM> > >& p_clust,map<int,double> &output_vol,
                              vector<stl_ptr<ATOM> >& lig_atoms,double alpha1=3.2,double alpha2=10.0,bool refine=true,double min_vol=300.0,
                              bool visualize=false,double def_rad=4.0,double surr_radius = 0.0,bool get_cavity_obj = false,
                              bool profile=false);
        void get_cavity_profile(DT_SOLVER& mysolver,vector<stl_ptr<TETRAHEDRON> >& pt,double alpha1,double alpha2,string& fname);
        void get_cavity_profile2(DT_SOLVER& mysolver,tr1::unordered_map<int64_t,int>& id2idx,vector<stl_ptr<TETRAHEDRON> >& pt,
                                 double alpha1,double alpha2,string& fname);
        void get_alt_loc_ligs(bool prot_mode = false); //neue Ligand-Objekte fuer jede alternative Location
        void merge_ligands(); //!wenn im pdb ein Peptidteil zwischen den HETATMs war, dann pruefen auf Verbindungen
        void merge_peptides(); //!sollte merge_ligands ersetzen
        void shorten_lig_names(); //!um ewig lange Namen von Peptiden zu kuerzen
        void water2lig();
        void metal2lig();
        void align_by_protein(PROTEIN *ref,bool spatial = true,bool debug=false,bool only_C_alpha=true);
        bool align_by_protein2(PROTEIN *ref,bool spatial = true,float const& max_displacement=1.4);
        float get_shape_sim(stl_ptr<STRUCTURE> const& ref);
        void align_by_cavity(STRUCTURE *ref); //!Fuer beide Strukturen muss VORHER get_cavity_auto() aufgerufen werden!!!
        vec3d<float> get_ligand_diversity();
};


//==================================================================================================
//Deklaration der Klasse CELL:
//==================================================================================================
class CELL {
    public:
        bool empty;
        bool surf;
        bool real_surf;
        bool next_to_surf;
        int burial;
        int friends;
        int index;
        vec3d<float> middle;
        map<int,stl_ptr<CELL> > n_cell_map;
        CELL();
        ~CELL();
};


//==================================================================================================
//Deklaration der Klasse BURIAL_CELL:
//==================================================================================================
class BURIAL_CELL {
    public:
        bool empty;
        bool surf;
        bool next_to_surf;
        bool real_surf;
        int burial;
        int friends;
        vec3d<float> middle;
        vector<stl_ptr<ATOM> > atoms;
        BURIAL_CELL();
        ~BURIAL_CELL();
};

//==================================================================================================
//Deklaration der Klasse STRUCTURE_GRID:
//==================================================================================================
class STRUCTURE_GRID {
    protected:
        vector<stl_ptr<ATOM> > atoms;
        float cell_size;
        int overadd;
        float half_size;
        vec3d<int> sizes; //x,y,z-Ausdehnungen (Anzahl der Zellen)
        vec3d<float> min; //Eckpunkt
        vec3d<float> max; //Eckpunkt
        vec3d<float> center;
        vec3d<float> buf; //fuer Zwischenrechnungen
        int grid_pointer[3]; //zeigt auf eine bestimmte Gitter-Stelle
        
        inline void calc_properties(); //!Ausgelagert aus dem Konstruktor
        inline void generate_grid(); //!eigentliches Gitter erzeugen
        inline void get_cell_middle();
    public:
        typedef GRID<3,CELL>::iterator g_iterator;
        stl_ptr<GRID<3,CELL> > grid; //!Das eigentliche Gitter
        
        STRUCTURE_GRID(vector< stl_ptr<ATOM> > &atms, float c_size, int add = 0);
        STRUCTURE_GRID(PROTEIN *protein, float c_size, int add = 0);
        ~STRUCTURE_GRID();
        
        inline vec3d<int> get_sizes();
        inline void get_index(int *index,vec3d<float> &coord);
        inline void check_empty(float add_to_vdw = 0.);
        inline void check_empty_buried(float add_to_vdw = 0.);
        inline void get_next_to_surf(); //!vorher check_empty aufrufen!!!
        inline void get_next_to_surf_buried();
        inline void get_real_surf(int min_buried = 18,int min_clust = 8); //!vorher get_next_to_surf aufrufen!!!
        inline void get_buriedness();
};


//==================================================================================================
//Deklaration der Klasse BURIAL_STRUCTURE_GRID:
//==================================================================================================
class BURIAL_STRUCTURE_GRID {
    protected:
        vector<stl_ptr<ATOM> > atoms;
        float cell_size;
        int overadd;
        float half_size;
        vec3d<int> sizes; //x,y,z-Ausdehnungen (Anzahl der Zellen)
        vec3d<float> min; //Eckpunkt
        vec3d<float> max; //Eckpunkt
        vec3d<float> center;
        vec3d<float> buf; //fuer Zwischenrechnungen
        int grid_pointer[3]; //zeigt auf eine bestimmte Gitter-Stelle
        
        inline void calc_properties(); //!Ausgelagert aus dem Konstruktor
        inline void generate_grid(); //!eigentliches Gitter erzeugen
        inline void get_cell_middle();
    public:
        typedef GRID<3,BURIAL_CELL>::iterator g_iterator;
        stl_ptr<GRID<3,BURIAL_CELL> > grid; //!Das eigentliche Gitter
        
        BURIAL_STRUCTURE_GRID(vector< stl_ptr<ATOM> > &atms, float c_size, int add = 0);
        BURIAL_STRUCTURE_GRID(PROTEIN *protein, float c_size, int add = 0);
        ~BURIAL_STRUCTURE_GRID();
        
        inline vec3d<int> get_sizes();
        inline void get_index(int *index,vec3d<float> &coord);
        inline void check_empty(float add_to_vdw = 0.);
        inline void get_next_to_surf(); //!vorher check_empty aufrufen!!!
        inline void get_real_surf();
        inline void get_buriedness();
};


//==================================================================================================
//Deklaration der Klasse COORD_GRID:
//==================================================================================================

class COORD_GRID {
    public:
        GRID<3,vec3d<float> > grid;
        float spacing;
        vec3d<int> sizes;
        vec3d<float> min; //Eckpunkt
        vec3d<float> max; //Eckpunkt
        COORD_GRID(vector<stl_ptr<ATOM> > &atoms, float spa, float add_size); //die Koordinaten der uebergebenen Atome zum
                                                                              //Bestimmen der Bindetasche nehmen und das
                                              //Gitter noch um add_size in jede Richtung vergroessern
        COORD_GRID(vec3d<float> &mi,vec3d<float> &ma,float spa,float add_size);
        ~COORD_GRID();
};

#endif
