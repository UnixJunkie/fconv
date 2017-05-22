
//============================================================================
// crystallizer_GN.h -*- C++ -*-; building crystal packings
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
// This library provides routines for building crystal packings from
// information contained in CIF files or in CRYSIN entries of MOL2 files.
//============================================================================

#ifndef CRYSTALIZERGN_H
#define CRYSTALIZERGN_H


#include<stdlib.h>
#include<iostream>
#include<string>
#include<fstream>
#include<sstream>
#include<vector>
#include<map>

#include"linalg_GN.hpp"
#include"message_GN.hpp"
#include"atom_GN.h"
#include"molecule_GN.h"
#include"structure_GN.h"


using namespace std;
using namespace TXT;

class STRUCTURE;


//==================================================================================================
//Deklaration der Klasse CRYSTALIZER:
//==================================================================================================

class CRYSTALIZER {
    private:
        ofstream f_out; //fuer ein eventuelles pymol-visualisierungsfile
        int sym_ele; //Anzahl der Symmetrieelemente
        vector<string> sym_name; //namen der Symmetrieelemente
        int n_splitted;
        int n_mol; //Anzahl der Molekle pro unit_cell
        int n_ligs; //Molekle insgesamt
        bool visualize; //pymol file fuer Symmetrieelemente
        bool verbosity;
        float max_d; //maximal erlaubter Abstand eines Atoms zum mittleren Liganden
        float smax_d; //max_d^2
        stl_ptr<LIGAND> l;
        vec3d<float> cart_x;
        vec3d<float> cart_y;
        vec3d<float> cart_z;
    //    vector<int> non_asym_ids; //ID's der Atome die nicht zur asymmetrischen Einheit gehoeren
        bool is_copy(stl_ptr<LIGAND> &lig); //prueft, ob die Koordinaten von lig bereits in structure->ligands sind
        void correct_setting(vec3d<float> &v);
    //    void change_origin_choice(vec3d<float> &v,vec3d<float> &sym1,vec3d<float> &sym2);
        bool origin_one();
        void kill_oc1();
        void correct_pos(stl_ptr<LIGAND> &lig);
        void gl(vec3d<float> v1,vec3d<float> v2,vec3d<float> posi,vec3d<float> direction); //glide
        void sc(vec3d<float> axis,float angle,vec3d<float> posi,vec3d<float> direction); //screw
        void sci(vec3d<float> axis,float angle,vec3d<float> posi,vec3d<float> direction,vec3d<float> point); //Rotoinversion
        void in(vec3d<float> point); //inversion
        void tr(vec3d<float> trans); //Translation
        void t(stl_ptr<CRYSIN_POSITION> &p); //Position
        bool build_unit_cell_from_positions();
        bool build_unit_cell(); //Einheitszelle aufbauen
    //    void correct_unit_cell();
    //    void duplicate_unit_cell(); //Einheitszelle symmetrisch vervielfachen
        void duplicate_unit_cell_splitted();
    //    void duplicate_unit_cell2(int number);
        void duplicate_unit_cell2_splitted(int number);
    //    void get_asym_ids();
    public:
        STRUCTURE *structure; //!muss selbst an parser zum schreiben uebergeben werden
        CRYSTALIZER(LIGAND &ligand,bool vis = false,bool verb = true);
        CRYSTALIZER(stl_ptr<LIGAND> &ligand,bool vis = false,bool verb = true);
        ~CRYSTALIZER();
    //    bool build_crystal(float max_dst = 8.);
        int build_crystal_splitted(float max_dst = 8.);
    //    bool build_crystal2(int number);
        int build_crystal2_splitted(int number);
};

#endif
