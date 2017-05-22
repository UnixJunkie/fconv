
//============================================================================
// protein_GN.h -*- C++ -*-; representation of molecular protein data
//
// Copyright (C) 2006, 2007, 2008, 2009, 2010 Gerd Neudert
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
// This library implements protein objects, which are used
// as high-level abstraction for molecular data
//============================================================================


#ifndef PROTEINGN_H
#define PROTEINGN_H

#include<stdlib.h>
#include<iostream>
#include<string>
#include<sstream>
#include<vector>
#include<map>

#include"linalg_GN.hpp"
#include"stl_ptr_GN.hpp"
#include"message_GN.hpp"
#include"atom_GN.h"
#include"molecule_GN.h"
#include"bron_kerbosch_GN.hpp"

using namespace std;
using namespace TXT;


class STRUCTURE;

class PROTEIN;


//==================================================================================================
//Deklaration der Klasse PROTEIN:
//==================================================================================================

class PROTEIN {
    public:
        stl_ptr<STRUCTURE> main_structure;
        vector< stl_ptr<CHAIN> > chains;
        stl_ptr<ATOM> last_C;
        string last_name;
        int bnd_id;
        float aligned_rmsd;
        matrix<float> optalign_rotm;
        vec3d<float> optalign_trans[2];
        PROTEIN();
        ~PROTEIN();
        void get_aacids();
        void get_atom_types(int mode,const char *def_file = "X",bool fill_X = false); //weist allen Proteinatomen Sybylatomtypen zu
        void build_bonds(); //achtung: BOND-Objekte aus den CONECT-Eintraegen gibt es schon
                            //         fehlen also nur noch die Verknuepfungen in den Ketten
        void ter_correct(); //fuer Umwandlung nach mol2 muessen die intern_id's neu berechnet werden
                            //weil die TER-Eintraege wegfallen
        void get_reli_centers();
        void visualize_reli_centers();
        void align(PROTEIN *ref,bool debug=false,bool only_C_alpha=true); //! Strukturalignment basierend auf Sequenzalignment
//        bool align2(PROTEIN *ref,unsigned int const& early_term=100); //! Strukturalignment basierend auf C-alpha Graphmatching
        float align2(PROTEIN *ref,float const& max_displace=1.4,bool const& sim_only=false);
        void align2b(PROTEIN *ref,vector<stl_ptr<ATOM> >&atms_a,vector<stl_ptr<ATOM> >&atms_b);
        void get_CA_poc_atoms(vector< stl_ptr<ATOM> > &cv);
        pair<int,float> get_similarity(vector< stl_ptr<ATOM> > &ref);
        void search_spattern(PROTEIN *ref,multimap<float,vector<stl_ptr<ATOM> > > &mmap,
                             vector<vec3d<float> > &t1,vector<vec3d<float> > &t2,vector<matrix<float> > &rm);
        
        void search_spattern2(PROTEIN *ref,multimap<float,vector<stl_ptr<ATOM> > > &mmap,
                              vector<vec3d<float> > &t1,vector<vec3d<float> > &t2,vector<matrix<float> > &rm);
        
        void search_spattern3(PROTEIN *ref,tr1::unordered_set<int> &set1,tr1::unordered_set<int> &set2);
};


#endif
