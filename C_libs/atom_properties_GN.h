
//============================================================================
// atom_properties_GN.h -*- C++ -*-; contains type specific atom properties
//
// Copyright (C) 2007, 2008, 2009, 2010, 2011 Gerd Neudert
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
// This library contains some mappings of atom types and atom type specific
// properties like vdW-radii, molecular weights, etc.
//============================================================================

#ifndef ATOMPROPERTIESGN_H
#define ATOMPROPERTIESGN_H

#include<iostream>
#include<string>
#include<vector>
#include<tr1/unordered_map>
#include<tr1/unordered_set>

using namespace std;
//==================================================================================================
//Atomtypen:
//==================================================================================================

//!Anzahl der internen Atomtypen:
extern const int n_intern_types;

//!Versionsstring fuer die atom types:
extern const string a_t_version;

extern const int n_metal_list;
extern const string metal_list[19];
extern const string metal_list_ra[19];

//!es folgen die internen Atomtypen (mode 0):
extern const string i_t[160];

//!es folgen die entsprechenden Sybyltypen (mode 1):
extern const string mode_1[160];

//!es folgen die modifizierten Sybyltypen (mode 2):
extern const string mode_2[160];

//!Jetzt die Anzahl der Valenzen fuer jeden Atomtyp:
extern const int atom_valence[160];

//!Jetzt die Hybridisierungen aller Typen:
extern const int atom_hybridizations[160];

//!Bindungsgeometrien:
extern const int atom_geometrie[160];

//! X--H  Bindungslaengen:
extern const float atom_hb_length[160];

//! Defaults fuer einige Flags:
extern const bool prot_acids_default;
extern const bool prot_guanidin_default;
extern const bool prot_amidin_default;
extern const bool prot_amin_default;
extern const bool prot_phosphate_default;
extern const bool prot_sulfate_default;
extern const bool get_bonds_default;
extern const bool kekulize_aromatics_default;
extern const bool kekulize_charged_default;
extern const bool allow_charged_aromatics_default;
extern const int max_ring_members_default;



namespace def_map {
    //!und hier die Zuweisungsmap fuer die Atomtypen:
    extern tr1::unordered_map<string,string> a_t_map;
    
    //!eine zweite map fuer anwendungen in denen immer zwischen 2 maps gewechselt wird:
    extern tr1::unordered_map<string,string> a_t_map2;
    
    //!ein flag, mit welchem def_file die map bereits befuellt ist:
    extern string curr_a_t_map; //X / def_file / m0 / m1 / m2
    
    extern string curr_a_t_map2;

    extern bool curr_prot_acids_default;
    extern bool curr_prot_guanidin_default;
    extern bool curr_prot_amidin_default;
    extern bool curr_prot_amin_default;
    extern bool curr_prot_phosphate_default;
    extern bool curr_prot_sulfate_default;
    extern bool curr_get_bonds_default;
    extern bool curr_kekulize_aromatics_default;
    extern bool curr_kekulize_charged_default;
    extern bool curr_allow_charged_aromatics_default;
    extern int curr_max_ring_members_default;
    extern bool curr_prot_acids_default2;
    extern bool curr_prot_guanidin_default2;
    extern bool curr_prot_amidin_default2;
    extern bool curr_prot_amin_default2;
    extern bool curr_prot_phosphate_default2;
    extern bool curr_prot_sulfate_default2;
    extern bool curr_get_bonds_default2;
    extern bool curr_kekulize_aromatics_default2;
    extern bool curr_kekulize_charged_default2;
    extern bool curr_allow_charged_aromatics_default2;
    extern int curr_max_ring_members_default2;
    
    extern int last_a_t_map;
}

namespace atom_properties {
    //! Muss hier erstmal alles extern deklariert werden!!!
    //! Wenn hier schon definiert wird, so bekommen alle Module die atom_properties einbinden eigene Versionen
    //! der Variablen => Overhead UND auch nach initialize() UNGEFUELLTE maps!!!
    extern const string normal_aacids[31];
    extern const string modified_aacids[11];
    extern const string nucleic_acids[11];
    extern const string known_elements[29];
    extern const int number_of_elements;
    extern const string vdW_elements[25];//[number_of_elements];
    extern const float vdW_radii[25];//[number_of_elements];
    extern const float clash_radii[25];
    extern const float mweights[25];//[number_of_elements];
    extern tr1::unordered_map<string,float> vdW_map; //!vdW-Radien
    extern tr1::unordered_map<string,float> square_vdW_map;
    extern tr1::unordered_map<string,float> clash_map; //!Atom-Radien
    extern tr1::unordered_map<string,float> square_clash_map;
    extern tr1::unordered_map<string,float> mol_weight;
    extern tr1::unordered_set<string> known_aacids;
    extern tr1::unordered_set<string> known_mod_aacids;
    extern tr1::unordered_set<string> known_nucleic_acids;
    void initialize();
}

#endif

