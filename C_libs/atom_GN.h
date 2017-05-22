
//============================================================================
// atom_GN.h -*- C++ -*-; objects representing atoms and related stuff
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
// This library provides objects representing atoms, bonds, rings and some
// PDB specific entries.
//============================================================================

#ifndef ATOMGN_H
#define ATOMGN_H

#include<stdlib.h>

#include"linalg_GN.hpp"
#include"stl_ptr_GN.hpp"
#include"message_GN.hpp"

#include<string>
#include<sstream>
#include<vector>
#include<iostream>
#include<tr1/unordered_set>
#include<tr1/unordered_map>

class ATOM;
class ATOM_EXT;
class COMPARENODE;
class RING;
class NODE;
class BOND;
class PDBCONECT;
class HET;
class COMMENT;
class WATER;
class METAL;

using namespace std;
using namespace TXT;

typedef vector< stl_ptr<ATOM> >::iterator atoms_vec;
typedef vector< stl_ptr<ATOM> >::const_iterator const_atoms_vec;
typedef vector< stl_ptr<BOND> >::iterator bonds_vec;
typedef vector< stl_ptr<BOND> >::const_iterator const_bonds_vec;
typedef vector< stl_ptr<COMPARENODE> >::iterator compares_vec;
typedef vector< stl_ptr<COMPARENODE> >::const_iterator const_compares_vec;
typedef vector< stl_ptr<RING> >::iterator rings_vec;
typedef vector< stl_ptr<RING> >::const_iterator const_rings_vec;
typedef vector< stl_ptr<NODE> >::iterator nodes_vec;
typedef vector< stl_ptr<NODE> >::const_iterator const_nodes_vec;
typedef vector< stl_ptr<PDBCONECT> >::iterator conects_vec;
typedef vector< stl_ptr<PDBCONECT> >::const_iterator const_conects_vec;
typedef vector< stl_ptr<WATER> >::iterator waters_vec;
typedef vector< stl_ptr<WATER> >::const_iterator const_waters_vec;
typedef vector< stl_ptr<METAL> >::iterator metals_vec;
typedef vector< stl_ptr<METAL> >::const_iterator const_metals_vec;
typedef vector< stl_ptr<HET> >::iterator het_entrys_vec;
typedef vector< stl_ptr<HET> >::const_iterator const_het_entrys_vec;
typedef vector< stl_ptr<COMMENT> >::iterator comments_vec;
typedef vector< stl_ptr<COMMENT> >::const_iterator const_comments_vec;

//==================================================================================================
//Deklaration der Klasse ATOM:
//==================================================================================================

class ATOM {
    public:
        int id; //6-10 im pdb
        int intern_id; //interne fortlaufende Nummer (id soll spaeter ohne h_correction auskommen!)
        string name; //12-15 im pdb
        int type; //1: pdb-ATOM / 2: pdb-HETATM / 3: pdb-metal / 4: pdb-water / 5:mol2-Atom
        string sybyl_type; //Der mol2-Atomtyp
        string intern_type; // X:nicht zugewiesen
        char alt_loc_id; //16 im pdb
        string res_name; //17-19 im pdb / sub_name vom mol2
        char chain_id; //21 im pdb
        int res_number; // 22-25 im pdb / sub_id vom mol2
        char alt_res; //26 fuer alternativen Rest
        vec3d<float> coord; //30-53 im pdb
        float occupancy; //54-59 im pdb
        float b_factor; //60-65 im pdb
        string element; //76-77 im pdb
        string charge; //78-79 im pdb
        bool is_ter; //folgt ein TER-Entry ?
        bool is_model; //folgt ein MODEL ?
        bool is_endmdl; //folgt ein ENDMDL ?
        int model_number; //0 wenn keine MODEL-Entrys sonst 1..n
        string dict_type; //fuer dict-Eintrag in mol2 file
        int bond_ind; //fuer Verknpfungsregeln bei Aminosaeuren
        string full_res_name; //res_name+res_number+alt_res
        vector< stl_ptr<ATOM> > bonded_atoms; //Zeiger auf gebundene Atome (werden nicht vom Destruktor geloescht)
        stl_ptr<ATOM_EXT> ext; //erweiterte Eigenschaften (wird nicht mit Copykonstruktor kopiert)
        ATOM();
        ATOM(ATOM const& atom); //uebernimmt keine bonded_atoms und keinen ext Zeiger
        ~ATOM();
        ATOM& operator=(ATOM const& atom);
        /*inline*/ void remove_ext(); //um Speicherplatz freizugeben
        bool operator<(ATOM const& rechts);
        inline bool operator>(ATOM const& rechts);
        bool operator==(ATOM const& rechts);
        /*inline*/ bool operator!=(ATOM const& rechts);
        /*inline*/ void get_compare_number();
        /*inline*/ void get_compare_number2();
        void get_element(bool reset_names = true);
        
        bool is_equal(ATOM &ref,vector<stl_ptr<ATOM> > &a_roots,vector<stl_ptr<ATOM> > &b_roots,
                      tr1::unordered_map<int,tr1::unordered_set<int> > &known,
                      tr1::unordered_map<int,tr1::unordered_set<int> > &negknown);
        
        bool is_in_same_ring(stl_ptr<ATOM>& ref);
        
        float const& get_smallest_bond_angle();
        float const& get_trigo_angle();

        friend ostream &operator<<(ostream &os,ATOM const& ref);
};


//==================================================================================================
//Deklaration der Klasse ATOM_EXT:
//==================================================================================================
// Diese Klasse dient der Erweiterung normaler ATOM Objekte um weitere Merkmale:
class ATOM_EXT {
    public:
        //!mal checken, ob hier noch alles gebraucht wird!!!
        unsigned int compare_number; //war frueher in COMPARENODE
        unsigned int compare_number2;
        
        stl_ptr<COMPARENODE> root_node;
        
        int hybridization; //! 0:unbekannt; 1:sp; 2:sp2; 3:sp3; 
        int n_heavy_bonded; //! 0:nicht ermittelt
        float sp; //Summe der Wahrscheinlichkeiten fuer sp-Hybridisierung
        float sp2;
        float sp3;
        float buriedness;
        float smallest_ba;
        float trigo_angle;
        bool check_hyb;
        bool sec_check;
        bool is_aromatic;
        bool is_ring;
        bool is_planar_ring;
        bool not_full_sp2_ring;
        bool pos_metal; //!bei dem Atom koennte es sich aufgrund des namens auch um ein Metall handeln
        vector<stl_ptr<ATOM> > existing_bonds; //!um feststellen zu koennen, ob es ein BOND Object schon gibt
        tr1::unordered_set<RING*> ring_ptrs; //! neu 05.08.2010: um festzustellen, ob 2 Atome im gleichen Ring sind
        stl_ptr<ATOM> prev;
        int level;
        ATOM_EXT();
        ATOM_EXT(ATOM_EXT const& ext);
        ~ATOM_EXT();
        void clear();
};


//==================================================================================================
//Deklaration der Klasse COMPARENODE:
//==================================================================================================

class COMPARENODE {
    public:
        stl_ptr<ATOM> root_atom;
        vector<stl_ptr<COMPARENODE> > sub_nodes;
        vector<stl_ptr<ATOM> > prev_roots;
        inline void make_subnodes();
        inline bool is_equal(COMPARENODE& cnode);
        COMPARENODE(stl_ptr<ATOM> &atom,vector<stl_ptr<ATOM> > &prev);
        ~COMPARENODE();
};


//==================================================================================================
//Deklaration der Klasse RING:
//==================================================================================================

class RING {
    public:
        int n_hetero;
        int n_members;
        int id;
        int sq_id;
        int n_pi_ele;
        bool is_full_sp2;
        bool is_planar;
        bool is_aromatic;
        bool is_pos;
        bool is_prot;
        vector<stl_ptr<ATOM> > ring_atoms;
        tr1::unordered_set<int> atom_set; // nur bei Bedarf befuellen
        RING();
        ~RING();
        bool fused(stl_ptr<RING> const& rg);
};


//Hilfsklasse fuer get_rings:
class NODE {
    public:
        bool is_end;
        int level;
        int node_id;
        stl_ptr<ATOM> atm;
        stl_ptr<NODE> prev;
        NODE(stl_ptr<ATOM> &at);
        ~NODE();
};


//==================================================================================================
//Deklaration der Klasse BOND:
//==================================================================================================

class BOND {
    public:
        int id; //Bond-id
        bool free_rot; //freie Drehbarkeit -- per default false
        stl_ptr<ATOM> from; //werden nicht vom Destruktor geloescht (=> BOND loeschen, wenn die ATOMs geloescht werden)
        stl_ptr<ATOM> to;
        string dict_type; //fuer Protein-mol2
        string type; //1 : single  2 : double  3 : triple  am: amid  ar: aromatisch  un: unknown
        bool operator<(BOND const& rechts);
        BOND();
        ~BOND();
};


//==================================================================================================
//Deklaration der Klasse PDBCONECT:
//==================================================================================================

class PDBCONECT { //fuer je eine CONECT-Line im pdb-File
    public:
        int from_id;
        stl_ptr<ATOM> from; //wird nicht vom Destruktor geloescht
        int cov_bonds[4];
        int h_bonds[4];
        int salt_bonds[2];
        PDBCONECT();
        ~PDBCONECT();
};


//==================================================================================================
//Deklaration der Klasse HET:
//==================================================================================================

class HET {
    public:
        string res_name; //7-10 im pdb  --  the het_id
        char chain_id; //12 im pdb
        HET();
        ~HET();
};


//==================================================================================================
//Deklaration der Klasse COMMENT:
//==================================================================================================

class COMMENT {
    public:
        string text;
        COMMENT();
        ~COMMENT();
};


//==================================================================================================
//Deklaration der Klasse WATER:
//==================================================================================================

class WATER {
    public:
        stl_ptr<ATOM> atom;
        vec3d<float> v1;
        vec3d<float> v2;
        bool calculated;
        WATER(bool calc = false);
        ~WATER(); //muss das ATOM-Objekt zerstoeren
};


//==================================================================================================
//Deklaration der Klasse METAL:
//==================================================================================================

class METAL {
    public:
        stl_ptr<ATOM> atom;
        bool is_covalent;
        METAL(bool cov = false);
        ~METAL(); //muss das ATOM-Objekt zerstoeren
};


ostream &operator<<(ostream &os,ATOM const& ref);

#endif //atom_GN.h

