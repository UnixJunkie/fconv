
//============================================================================
// molecule_GN.h -*- C++ -*-; representation of molecules and related stuff
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
// This library implements the most important features of fconv. It supplies
// objects for molecules and related molecular data and it implements the
// complete atom type perception routines.
//============================================================================


#ifndef MOLECULEGN_H
#define MOLECULEGN_H

#include<stdlib.h>
#include<iostream>
#include<string>
#include<sstream>
#include<fstream>
#include<vector>
#include<map>
#include<tr1/unordered_map>
#include<set>
#include<tr1/unordered_set>
#include<queue>

#include"ELEMENT_DATA_new_GN.h"
#include"linalg_GN.hpp"
#include"stl_ptr_GN.hpp"
#include"message_GN.hpp"
#include"optalign_GN.hpp"
#include"seqalign_GN.hpp"
#include"bron_kerbosch_GN.hpp"
#include"string_fu_GN.hpp"
#include"atom_GN.h" //ATOM BOND RING usw.
#include"structure_additional_GN.h" //CRYSIN usw.
#include"atom_properties_GN.h"
#include"octree_GN.hpp"
#include"mw_match_GN.h"


using namespace std;
using namespace TXT;

class STRUCTURE;
class PROTEIN;
class CRYSIN;
class CRYSIN_POSITION;


class AACID;
class CAVITY;
class RELIBASE_CENTER;
class MOLECULE;
class LIGAND;
class CHAIN;


typedef vector< stl_ptr<LIGAND> >::iterator ligands_vec;
typedef vector< stl_ptr<LIGAND> >::const_iterator const_ligands_vec;
typedef vector< stl_ptr<CHAIN> >::iterator chains_vec;
typedef vector< stl_ptr<CHAIN> >::const_iterator const_chains_vec;
typedef vector< stl_ptr<CAVITY> >::iterator cavities_vec;
typedef vector< stl_ptr<CAVITY> >::const_iterator const_cavities_vec;
typedef vector< stl_ptr<AACID> >::iterator aacids_vec;
typedef vector< stl_ptr<AACID> >::const_iterator const_aacids_vec;
typedef map<int,stl_ptr<AACID> >::iterator aacids_map;
typedef map<int,stl_ptr<AACID> >::const_iterator const_aacids_map;
typedef vector< stl_ptr<RELIBASE_CENTER> >::iterator reli_centers_vec;
typedef vector< stl_ptr<RELIBASE_CENTER> >::const_iterator const_reli_centers_vec;
typedef BK_NODE<stl_ptr<ATOM> > bkn_node;
typedef BK_NODE<stl_ptr<ATOM> >* bkn_pointer;
typedef vector<BK_NODE<stl_ptr<ATOM> >* > bkn_container;
typedef vector<BK_NODE<stl_ptr<ATOM> >* >::iterator bkn_vec;
typedef vector<BK_NODE<stl_ptr<ATOM> >* >::const_iterator const_bkn_vec;
typedef map<unsigned int,vector<BK_NODE<stl_ptr<ATOM> >* > > clique_container;
typedef map<unsigned int,vector<BK_NODE<stl_ptr<ATOM> >* > >::iterator clique_map;
typedef map<unsigned int,vector<BK_NODE<stl_ptr<ATOM> >* > >::const_iterator const_clique_map;


void write_def_file(int const& mode,const char *name,vector<string> *alt_mode = 0,
                    bool prot_acids = prot_acids_default,
                    bool prot_guanidin = prot_guanidin_default,
                    bool prot_amidin = prot_amidin_default,
                    bool prot_amin = prot_amin_default,
                    bool prot_phosphate = prot_phosphate_default,
                    bool prot_sulfate = prot_sulfate_default,
                    bool get_bonds = get_bonds_default,
                    int max_ring_members = max_ring_members_default,
                    bool kekulize_aromatics = kekulize_aromatics_default,
                    bool kekulize_charged = kekulize_charged_default,
                    bool allow_charged_aromatics = allow_charged_aromatics_default);
void write_def_file(int const& mode,const char *name,vector<string> *alt_mode,
                    tr1::unordered_map<string,string> new_flags);
void load_def_file_global(char const* name,
                          bool &prot_acids,
                          bool &prot_guanidin,
                          bool &prot_amidin,
                          bool &prot_amin,
                          bool &prot_phosphate,
                          bool &prot_sulfate,
                          bool &get_bonds,
                          int &max_ring_members,
                          bool &kekulize_aromatics,
                          bool &kekulize_charged,
                          bool &allow_charged_aromatics);


//==================================================================================================
//Deklaration der Klasse MOLECULE:
//==================================================================================================

class MOLECULE {
    protected:
        tr1::unordered_map<string,int> key_map;
        int freerot_bonds;
        int n_heavy_atoms;
        int compare_sum;
        int verbosity;
        int curr_int_id; //!eingefuehrt fuer set_standard_protonation

        bool check_trigo(stl_ptr<ATOM> const& base,float lim = 0.78);
        inline bool is_linear(stl_ptr<ATOM> const& base);
        inline void delete_worst_bonded(stl_ptr<ATOM> const& base);
        inline bool check_db_tors(stl_ptr<ATOM> const& base,float const& lim = 0.36); // war mal auf 0.36
        inline bool check_planar_tors(stl_ptr<ATOM> const& from,stl_ptr<ATOM> const& to,float const& lim = 0.36); // war mal auf 0.36
        inline void get_connections();
        void get_hybridization(bool const& get_bonds,int const& max_ring_members = 10,
                               bool const& allow_protonated = true,bool const& kekulize_aromatics = false);
        inline void check_ring_plan();
        void get_rings(int const& max_members = 10);
        bool check_SO(stl_ptr<ATOM> const& base,int const& O_count,int const& OH_count,int const& O_2heavy,int const& N_count);
        bool check_PO(stl_ptr<ATOM> const& base,int const& O_count,int const& OH_count,int const& O_2heavy);
        void set_intern_types(bool const& prot_acids,bool const& prot_guanidin,bool const& prot_amidin,
                              bool const& prot_amin,bool const& prot_phosphate,bool const& prot_sulfate,
                              bool const& kekulize_charged);
        void generate_new_def_file(char const* name);
        void make_bonds(bool const& no_free_rot_planar,
                        bool const& prot_acids,
                        bool const& prot_guanidin,
                        bool const& prot_amidin,
                        bool const& prot_amin,
                        bool const& prot_phosphate,
                        bool const& prot_sulfate,
                        bool const& kekulize_aromatics,
                        bool const& kekulize_charged);
        void add_bond_object(stl_ptr<ATOM> &at1, stl_ptr<ATOM> &at2,string const& tpy,bool test = true);
        void load_def_file(char const* name,
                           bool &prot_acids,
                           bool &prot_guanidin,
                           bool &prot_amidin,
                           bool &prot_amin,
                           bool &prot_phosphate,
                           bool &prot_sulfate,
                           bool &get_bonds,
                           int &max_ring_members,
                           bool &kekulize_aromatics,
                           bool &kekulize_charged,
                           bool &allow_charged_aromatics);
        inline void add_H(stl_ptr<ATOM> &atom,vec3d<float> &hvec);
        float get_hpos_score(vec3d<float> &hvec,stl_ptr<ATOM> &atom);
        void get_optimal_hpos(vec3d<float> &hvec,stl_ptr<ATOM> &atom,vec3d<float> &axis);
        void get_optimal_hpos2(vec3d<float> &hvec,vec3d<float> &hvec2,vec3d<float> &hvec3,stl_ptr<ATOM> &atom,vec3d<float> &axis);
        int get_closest_ring_id(stl_ptr<ATOM> const& from);
        int get_compare_sum();
    public:
        bool ele_already_set;
        bool atom_types_already_set;
        bool has_stereo;
        string name;
        int** SP_map;
        void calc_SP_map();
        void release_SP_map();
        
        matrix<float> optalign_rotm;    //!Rotationsmatrix fuer optimales Alignment
        vec3d<float> optalign_trans[2]; //!Verschiebungsvektoren fuer optalign
                                        //!Beides wird in get_rmsd berechnet (muss also vorher aufgerufen werden!)
        
        stl_ptr<STRUCTURE> main_structure;
        map< int,stl_ptr<ATOM> > sub_map; //haelt die root-atoms fr die substructures (mol2)
        vector< stl_ptr<ATOM> > atoms;
        vector< stl_ptr<ATOM> > atom_buf;
        vector< stl_ptr<BOND> > bonds;
        vector<stl_ptr<RING> > rings;
        vector<stl_ptr<COMMENT> > comments;
        MOLECULE();
        MOLECULE(MOLECULE const& lig);
        MOLECULE(LIGAND const& lig);
        bool copy_without_coords(LIGAND *lig);
        virtual ~MOLECULE(); //muss die ATOM-Obj. und die BOND-Obj. und die RING-Obj. zerstoeren
        void change_ele_set(bool const& new_stat) {ele_already_set = new_stat;}
        void dlg2mol2(int mode = 1,const char *def_file = "X",bool get_bonds = true,int verb = 1,
                      bool kill_ext = true,bool fill_X = false,const char *alt_def_file = "X",int max_ring_members = 10,
                      bool no_free_rot_planar = true,LIGAND* at_ref_mol = 0,LIGAND* at_ref_mol2 = 0,
                      bool prot_acids = prot_acids_default,
                      bool prot_guanidin = prot_guanidin_default,
                      bool prot_amidin = prot_amidin_default,
                      bool prot_amin = prot_amin_default,
                      bool prot_phosphate = prot_phosphate_default,
                      bool prot_sulfate = prot_sulfate_default,
                      bool kekulize_aromatics = kekulize_aromatics_default,
                      bool kekulize_charged = kekulize_charged_default,
                      bool allow_charged_aromatics = allow_charged_aromatics_default);
        void dlg2mol2(stl_ptr<LIGAND> &ref_mol);
        int get_freerot_bonds(); //Anzahl frei drehbarer Bindungen liefern
        bool check_equality(stl_ptr<LIGAND> &ref_mol,bool recalc_compare_numbers,bool stereo);
        
        template<class T>
        float get_rmsd(stl_ptr<T> &ref_mol,bool with_hyd = false,bool debug_mode = false,bool optalign = false);
        
        template<class T>
        float get_bk_rmsd(stl_ptr<T> &ref_mol,bool with_hyd = false,bool debug_mode = false,bool optalign = false,bool recalc_hybrid = true);
        
        int has_substructure(stl_ptr<LIGAND> &ref_mol,bool consider_hybrid = false,bool type_based = false,bool own_types = false);
        
        void opt_align(); //!get_rmsd muss vorher aufgerufen werden!!!

        void get_elements(bool reset_names = true);
        void get_bonded_only();
        void get_hybrid_only(bool prot_acids = prot_acids_default,
                             bool prot_guanidin = prot_guanidin_default,
                             bool prot_amidin = prot_amidin_default,
                             bool prot_amin = prot_amin_default,
                             bool prot_phosphate = prot_phosphate_default,
                             bool prot_sulfate = prot_sulfate_default, //! prot_sulfate betrifft auch Sulfonamide !!!
                             bool kekulize_aromatics = kekulize_aromatics_default,
                             bool kekulize_charged = kekulize_charged_default,
                             bool allow_charged_aromatics = allow_charged_aromatics_default);
        void get_atom_typing(int mode = 1,
                             bool get_ele = true,
                             const char *def_file = "X",
                             bool get_bonds = get_bonds_default, // Defaults werden in atom_properties_GN.cpp gesetzt
                             int verb = 1,
                             bool kill_ext = true,
                             bool fill_X = false,
                             const char *alt_def_file = "X",
                             int max_ring_members = max_ring_members_default,
                             bool no_free_rot_planar = true,
                             LIGAND* at_ref_mol = 0,
                             LIGAND* at_ref_mol2 = 0,
                             bool prot_acids = prot_acids_default,
                             bool prot_guanidin = prot_guanidin_default,
                             bool prot_amidin = prot_amidin_default,
                             bool prot_amin = prot_amin_default,
                             bool prot_phosphate = prot_phosphate_default,
                             bool prot_sulfate = prot_sulfate_default, //! prot_sulfate betrifft auch Sulfonamide !!!
                             bool kekulize_aromatics = kekulize_aromatics_default,
                             bool kekulize_charged = kekulize_charged_default,
                             bool allow_charged_aromatics = allow_charged_aromatics_default);
        void set_standard_protonation(int verb = 1); //!neu: Standardprotonierungen fuer CNOS setzen
        void rename_atoms(bool get_ele = true);
        void reorder_atoms(bool get_ele = true);
        int get_n_heavy();
        void kick_water(bool intern_is_set = false);
        void get_compare_numbers();
        void intern2mode(bool fill_X = false,const char *alt_def_file = "X");
        void try_for_def_file(int const& mode,const char *def_file,const char *alt_def_file,
                              bool &prot_acids,
                              bool &prot_guanidin,
                              bool &prot_amidin,
                              bool &prot_amin,
                              bool &prot_phosphate,
                              bool &prot_sulfate,
                              bool &get_bonds,
                              int &max_ring_members,
                              bool &kekulize_aromatics,
                              bool &kekulize_charged,
                              bool &allow_charged_aromatics);
        float get_mol_weight();
        float get_vdw_volume();
        float get_vdw_volume_grid();
        int get_number_of_rings(int max_ring_members = 10,bool debug_mode = false);
        
        void rmsd_by_mol2_residues(vector<string> &vrn,vector<float> &vtr,vector<float> &vbr,vector<float> &vrr,
                                   vector<float> &otr,vector<float> &obr,vector<float> &orr,stl_ptr<LIGAND> &ref,bool debug_mode);
};

//==================================================================================================
//Deklaration der Klasse LIGAND:
//==================================================================================================

class LIGAND : public MOLECULE {
    private:
        void get_connected_atoms(tr1::unordered_set<int> &con,stl_ptr<ATOM> &base);
    public:
        stl_ptr<CHAIN> next_chain;
        bool is_peptide;
        string type; //SMALL, BIOPOLYMER, PROTEIN, NUCLEIC_ACID, SACCHARIDE (noch nicht fuer pdb's)
        string charge_type;
        stl_ptr<CRYSIN> crysin;
        LIGAND();
        LIGAND(LIGAND const& lig);
        LIGAND(LIGAND const& lig,tr1::unordered_set<int> &at2take);
        ~LIGAND();
        int split(bool keep_biggest_fragment = false); //!falls mehrere Molekuele im Liganden sind werden diese in extra Ligand-Objekte
                     //!geschrieben. Die Anzahl der neuen Liganden wird zurueckgegeben. Diese werden
                     //!in main_structure->splitted_ligands abgelegt.
};


//==================================================================================================
//Deklaration der Klasse CHAIN:
//==================================================================================================

class CHAIN : public MOLECULE {
    public:
        string last_name;
        int n_aacids;
        stl_ptr<LIGAND> prev_lig;
        stl_ptr<PROTEIN> protein;
        map<int,stl_ptr<AACID> > aacids;  //res_number als key
    //    vector< stl_ptr<CENTER> > centers;
        vector< stl_ptr<RELIBASE_CENTER> > reli_centers;
        CHAIN();
        ~CHAIN(); //muss die AACIDs und die reli_centers loeschen
        void get_atom_types(int mode,const char *def_file = "X",bool fill_X = false);
        void get_aacids(int const& preset_n_aacids = 0);
        void build_bonds();
        void get_reli_centers();
};



//==================================================================================================
//Deklaration der Klasse AACID:
//==================================================================================================

class AACID {
    private:
        inline void create_bond(stl_ptr<ATOM> &from, stl_ptr<ATOM> &to, string type, string &d_type);
    public:
        bool last_res;
        string res_name; //ARG PHE ...
        int res_number; //original res_number aus dem pdb
        stl_ptr<CHAIN> chain; //Zeiger auf die entsprechende Kette zu der die AS gehoert
        vector< stl_ptr<ATOM> > atoms; //nur Zeiger auf die atoms von CHAIN => werden von CHAIN
        AACID();                       //geloescht !
        ~AACID();
        inline void get_intern_types();
        inline void build_bonds();
};


//==================================================================================================
//Deklaration der Klasse CAVITY:
//==================================================================================================

class CAVITY {
    public:
        vector< stl_ptr<AACID> > aacids;
        float total_volume;
        matrix<float> optalign_rotm;
        vec3d<float> optalign_trans[2];
        CAVITY();
        ~CAVITY();
        void align(stl_ptr<CAVITY> ref);
        bool operator<(CAVITY &rechts);
};

//==================================================================================================
//Deklaration der Klasse RELIBASE_CENTER:
//==================================================================================================

class RELIBASE_CENTER {
    public:
        vec3d<float> coord;
        string type;
        stl_ptr<AACID> aacid;
        RELIBASE_CENTER();
        ~RELIBASE_CENTER();
        
};


template<class T>
float MOLECULE::get_rmsd(stl_ptr<T> &ref_mol,bool with_hyd,bool debug_mode,bool optalign) {
    get_compare_numbers();
    int ref_n_heavy = 0;

    for (atoms_vec at=ref_mol->atoms.begin(); at!=ref_mol->atoms.end(); ++at) {
        (*at)->get_compare_number();
        if ((*at)->element == "H") continue;
        ++ref_n_heavy;
    }
    for (atoms_vec at=ref_mol->atoms.begin(); at!=ref_mol->atoms.end(); ++at) {
        (*at)->get_compare_number2();
    }

    map<stl_ptr<ATOM>,set<stl_ptr<ATOM> > > clist; //!Vergleichsmap die zu jedem Atom eine Liste mit gleichwertigen Atomen
                                                   //!der Referenzstruktur enthaelt
    tr1::unordered_map<int,tr1::unordered_set<int> > known;
    tr1::unordered_map<int,tr1::unordered_set<int> > negknown;

    int soll_clique_size = 0;
    for (atoms_vec at=atoms.begin(); at!=atoms.end(); ++at) {
        if ((*at)->element == "H") continue;
        ++soll_clique_size;
        clist[*at] = set<stl_ptr<ATOM> >();
        for (atoms_vec bt=ref_mol->atoms.begin(); bt!=ref_mol->atoms.end(); ++bt) {
            if ((*at)->element == "H") continue;
            if ((*at)->element == (*bt)->element) {
                vector<stl_ptr<ATOM> > prevs1;
                vector<stl_ptr<ATOM> > prevs2;
                if ((*at)->is_equal(**bt,prevs1,prevs2,known,negknown)) clist[*at].insert(*bt);
            }
        }
    }

    if (soll_clique_size != ref_n_heavy) {
        if (debug_mode) cerr << "-> different number of heavy atoms in both molecules" << endl;
        return -1.;
    }

    //! 2.) alle moeglichen Atom-Atom-Paare erzeugen (nur gleiche Elemente zusammen):
    bkn_container test_list;
    for (atoms_vec at=atoms.begin(); at!=atoms.end(); ++at) {
        if ((*at)->element == "H") continue;
        for (set<stl_ptr<ATOM> >::iterator bt=clist[*at].begin(); bt!=clist[*at].end(); ++bt) {
            bkn_pointer new_node = new bkn_node(*at,*bt);
            test_list.push_back(new_node);
        }
    }

    if (test_list.size() == 0) {
        if (debug_mode) cerr << "-> no nodes created" << endl;
        return -1.;
    }

    if (ref_mol->SP_map == 0) ref_mol->calc_SP_map();
    if (SP_map == 0) calc_SP_map();

    //! 3.) Jetzt die Kanten erzeugen:
    int edges = 0;
    for (bkn_vec it=test_list.begin(); it!=test_list.end(); ++it) {
        for (bkn_vec jt=it+1; jt!=test_list.end(); ++jt) {
            if ((*it)->obj1.equal_addr((*jt)->obj1) || (*it)->obj2.equal_addr((*jt)->obj2)) continue;
            bool bonded = false;
            if ((*it)->obj1->intern_id > (*jt)->obj1->intern_id) {
                if ((*it)->obj2->intern_id > (*jt)->obj2->intern_id) {
                    if (SP_map[(*it)->obj1->intern_id][(*jt)->obj1->intern_id] == ref_mol->SP_map[(*it)->obj2->intern_id][(*jt)->obj2->intern_id]) {
                        if (SP_map[(*it)->obj1->intern_id][(*jt)->obj1->intern_id] != 0) bonded = true;
                        else if (ref_mol->SP_map[(*it)->obj2->intern_id][(*jt)->obj2->intern_id] == 0) bonded = true;
                    }
                } else {
                    if (SP_map[(*it)->obj1->intern_id][(*jt)->obj1->intern_id] == ref_mol->SP_map[(*jt)->obj2->intern_id][(*it)->obj2->intern_id]) {
                        if (SP_map[(*it)->obj1->intern_id][(*jt)->obj1->intern_id] != 0) bonded = true;
                        else if (ref_mol->SP_map[(*jt)->obj2->intern_id][(*it)->obj2->intern_id] == 0) bonded = true;
                    }
                }
            } else {
                if ((*it)->obj2->intern_id > (*jt)->obj2->intern_id) {
                    if (SP_map[(*jt)->obj1->intern_id][(*it)->obj1->intern_id] == ref_mol->SP_map[(*it)->obj2->intern_id][(*jt)->obj2->intern_id]) {
                       if (SP_map[(*jt)->obj1->intern_id][(*it)->obj1->intern_id] != 0) bonded = true;
                       else if (ref_mol->SP_map[(*it)->obj2->intern_id][(*jt)->obj2->intern_id] == 0) bonded = true;
                    }
                } else {
                    if (SP_map[(*jt)->obj1->intern_id][(*it)->obj1->intern_id] == ref_mol->SP_map[(*jt)->obj2->intern_id][(*it)->obj2->intern_id]) {
                        if (SP_map[(*jt)->obj1->intern_id][(*it)->obj1->intern_id] != 0) bonded = true;
                        else if (ref_mol->SP_map[(*jt)->obj2->intern_id][(*it)->obj2->intern_id] == 0) bonded = true;
                    }
                }
            }

            if (bonded) {
                (*it)->edges.insert(*jt);
                (*jt)->edges.insert(*it);
                ++edges;
            }
        }
    }
    if (debug_mode) cerr << "-> " << edges << " edges and " << test_list.size() << " nodes" << endl;

    if (edges == 0) return -1.;

    //! 4.) Jetzt die Cliquen suchen:
    BK_SOLVER<stl_ptr<ATOM> > mysolver(test_list);
    mysolver.solve(soll_clique_size); //! hier kann optional auch eine Mindestgroesse fuer die Cliquen mitgegeben werden
    //! => mysolver.cliques  ist eine multimap, die nun alle gefundenen Cliquen enthaelt (Cliquengroesse als Keys)
    //!    mysolver.max_clique  ist ein unsigned int der die groesste Clique angibt

    if (mysolver.cliques.size() == 0) {
        if (debug_mode) cerr << "-> no cliques found" << endl;
        return -1.;
    }

    //! 5.) Es koennen mehrere groesste Cliquen gefunden worden sein => nun muss die Clique ermittelt werden, die
    //!     den besten RMSD liefert:
    has_stereo = false;
    float rmsd = 999999999999.;
    bkn_container* best_clique = 0;
    vector<vector<stl_ptr<ATOM> > > stereo_bonds;
    for (bonds_vec bnt=bonds.begin(); bnt!=bonds.end(); ++bnt) {
        if ((*bnt)->type != "2") continue;
        if ((*bnt)->from->ext->n_heavy_bonded < 2 ||
            (*bnt)->to->ext->n_heavy_bonded < 2) continue;
        // Wenn auf BEIDEN Seiten der DB jeweils 2 gleiche Substituenten sind
        // liegt keine Isomerie vor:
        if ((*bnt)->from->ext->n_heavy_bonded == 3 &&
            (*bnt)->to->ext->n_heavy_bonded == 3) {
            int spairs = 0;
            for (atoms_vec at=(*bnt)->from->bonded_atoms.begin(); at!=(*bnt)->from->bonded_atoms.end(); ++at) {
                if ((*at)->intern_id == (*bnt)->to->intern_id) continue;
                for (atoms_vec bt=at+1; bt!=(*bnt)->from->bonded_atoms.end(); ++bt) {
                    if ((*bt)->intern_id == (*bnt)->to->intern_id) continue;
                    for (set<stl_ptr<ATOM> >::iterator cta=clist[*at].begin(); cta!=clist[*at].end(); ++cta) {
                        if (clist[*bt].find(*cta) != clist[*bt].end()) {
                            ++spairs;
                            break;
                        }
                    }
                    if (spairs > 0) break;
                }
                if (spairs > 0) break;
            }
            for (atoms_vec at=(*bnt)->to->bonded_atoms.begin(); at!=(*bnt)->to->bonded_atoms.end(); ++at) {
                if ((*at)->intern_id == (*bnt)->from->intern_id) continue;
                for (atoms_vec bt=at+1; bt!=(*bnt)->to->bonded_atoms.end(); ++bt) {
                    if ((*bt)->intern_id == (*bnt)->from->intern_id) continue;
                    for (set<stl_ptr<ATOM> >::iterator cta=clist[*at].begin(); cta!=clist[*at].end(); ++cta) {
                        if (clist[*bt].find(*cta) != clist[*bt].end()) {
                            ++spairs;
                            break;
                        }
                    }
                    if (spairs > 1) break;
                }
                if (spairs > 1) break;
            }
            if (spairs > 1) continue;
        }
        vector<stl_ptr<ATOM> > a_vecs;
        for (atoms_vec at=(*bnt)->from->bonded_atoms.begin(); at!=(*bnt)->from->bonded_atoms.end(); ++at) {
            if ((*at)->element == "H") continue;
            if ((*at)->intern_id == (*bnt)->to->intern_id) continue;
            a_vecs.push_back(*at);
            break;
        }
        a_vecs.push_back((*bnt)->from);
        a_vecs.push_back((*bnt)->to);
        for (atoms_vec at=(*bnt)->to->bonded_atoms.begin(); at!=(*bnt)->to->bonded_atoms.end(); ++at) {
            if ((*at)->element == "H") continue;
            if ((*at)->intern_id == (*bnt)->from->intern_id) continue;
            a_vecs.push_back(*at);
            break;
        }
        if (a_vecs.size() != 4) continue;
        has_stereo = true;
        stereo_bonds.push_back(a_vecs);
    }

    for (clique_map it=mysolver.cliques.begin(); it!=mysolver.cliques.end(); ++it) {
        if (it->first < mysolver.max_clique) continue;

        //! Jetzt muss noch getestet werden, ob die Stereochemie der Clique stimmt
        //! Zunaechst asymmetrische Cs pruefen:
        bool correct_stereo = true;
        for (bkn_vec jt=it->second.begin(); jt!=it->second.end(); ++jt) {
            if ((*jt)->obj1->ext->n_heavy_bonded < 3) continue;
            if (check_trigo((*jt)->obj1)) continue;
            // Bei 4 heavy atoms muss noch gechecked werden, ob mindestens 3
            // verschieden sind:
            if ((*jt)->obj1->ext->n_heavy_bonded == 4) {
                int spairs = 0;
                for (atoms_vec at=(*jt)->obj1->bonded_atoms.begin(); at!=(*jt)->obj1->bonded_atoms.end(); ++at) {
                    if ((*at)->element == "H") continue;
                    for (atoms_vec bt=at+1; bt!=(*jt)->obj1->bonded_atoms.end(); ++bt) {
                        for (set<stl_ptr<ATOM> >::iterator cta=clist[*at].begin(); cta!=clist[*at].end(); ++cta) {
                            if (clist[*bt].find(*cta) != clist[*bt].end()) {
                                ++spairs;
                                break;
                            }
                        }
                    }
                    if (spairs > 1) break;
                }
                if (spairs > 1) continue;
            }

            vector<stl_ptr<ATOM> > a_vecs;
            vector<stl_ptr<ATOM> > b_vecs;
            for (atoms_vec at=(*jt)->obj1->bonded_atoms.begin(); at!=(*jt)->obj1->bonded_atoms.end(); ++at) {
                if ((*at)->element == "H") continue;
                a_vecs.push_back(*at);
                for (bkn_vec jjt=it->second.begin(); jjt!=it->second.end(); ++jjt) {
                    if ((*jjt)->obj1->intern_id == (*at)->intern_id) {
                        b_vecs.push_back((*jjt)->obj2);
                        break;
                    }
                }
                if (a_vecs.size() == 3) break;
            }
            if (a_vecs.size() != 3 || b_vecs.size() != 3) {
                if (debug_mode) cerr << "-> error in R/S symmetrie check for " << (*jt)->obj1->name << endl;
                break;
            }
            has_stereo = true;
            if (triple_product((*jt)->obj1->coord,a_vecs[0]->coord,a_vecs[1]->coord,a_vecs[2]->coord) > 0) {
                if (triple_product((*jt)->obj2->coord,b_vecs[0]->coord,b_vecs[1]->coord,b_vecs[2]->coord) < 0) {
                    correct_stereo = false;
                    /*
                    if (debug_mode) {
                        cerr << "no correct stereo for:" << endl;
                        cerr << (*jt)->obj1->name << " " << (*jt)->obj2->name << endl;
                        cerr << a_vecs[0]->name << " " << b_vecs[0]->name << endl;
                        cerr << a_vecs[1]->name << " " << b_vecs[1]->name << endl;
                        cerr << a_vecs[2]->name << " " << b_vecs[2]->name << endl << endl;
                    }
                    break;
                    */
                }
            } else {
                if (triple_product((*jt)->obj2->coord,b_vecs[0]->coord,b_vecs[1]->coord,b_vecs[2]->coord) > 0) {
                    correct_stereo = false;
                    /*
                    if (debug_mode) {
                        cerr << "no correct stereo for:" << endl;
                        cerr << (*jt)->obj1->name << " " << (*jt)->obj2->name << endl;
                        cerr << a_vecs[0]->name << " " << b_vecs[0]->name << endl;
                        cerr << a_vecs[1]->name << " " << b_vecs[1]->name << endl;
                        cerr << a_vecs[2]->name << " " << b_vecs[2]->name << endl << endl;
                    }
                    */
                    break;
                }
            }
        }

        if (!correct_stereo) continue;

        //! Jetzt E/Z Isomerie pruefen (cis/trans bei Ringen wird bereits bei
        //!                             den asym Cs mit erfasst):
        for (vector<vector<stl_ptr<ATOM> > >::iterator bnt=stereo_bonds.begin(); bnt!=stereo_bonds.end(); ++bnt) {
            vector<stl_ptr<ATOM> > b_vecs;
            for (atoms_vec at=bnt->begin(); at!=bnt->end(); ++at) {
                for (bkn_vec jt=it->second.begin(); jt!=it->second.end(); ++jt) {
                    if ((*jt)->obj1->intern_id == (*at)->intern_id) {
                        b_vecs.push_back((*jt)->obj2);
                        break;
                    }
                }
            }
            if (b_vecs.size() != 4) {
                if (debug_mode) cerr << "-> error in E/Z symmetrie check for "
                                     << (*bnt)[1]->name << " === " << (*bnt)[2]->name << endl;
                break;
            }
            if (dihedral((*bnt)[0]->coord,(*bnt)[1]->coord,(*bnt)[2]->coord,(*bnt)[3]->coord) < 1.571) {
                if (dihedral(b_vecs[0]->coord,b_vecs[1]->coord,b_vecs[2]->coord,b_vecs[3]->coord) > 1.571) {
                    correct_stereo = false;
                    /*
                    if (debug_mode) {
                        cerr << "no correct stereo for "
                             << a_vecs[0]->name << "--"
                             << a_vecs[1]->name << "=="
                             << a_vecs[2]->name << "--"
                             << a_vecs[3]->name << endl << endl;
                    }
                    break;
                    */
                }
            } else {
                if (dihedral(b_vecs[0]->coord,b_vecs[1]->coord,b_vecs[2]->coord,b_vecs[3]->coord) < 1.571) {
                    correct_stereo = false;
                    /*
                    if (debug_mode) {
                        cerr << "no correct stereo for "
                             << a_vecs[0]->name << "--"
                             << a_vecs[1]->name << "=="
                             << a_vecs[2]->name << "--"
                             << a_vecs[3]->name << endl << endl;
                    }
                    break;
                    */
                }
            }
        }

        if (!correct_stereo) continue;

        float test_rmsd = 0.;
        if (optalign) {
            vector<vec3d<float> > A_comp_list;
            vector<vec3d<float> > B_comp_list;
            for (bkn_vec jt=it->second.begin(); jt!=it->second.end(); ++jt) {
                A_comp_list.push_back((*jt)->obj1->coord);
                B_comp_list.push_back((*jt)->obj2->coord);
            }
            get_align_matrix(B_comp_list,A_comp_list,optalign_rotm,optalign_trans[0],optalign_trans[1]);
            for (vector<vec3d<float> >::iterator at=A_comp_list.begin(); at!=A_comp_list.end(); ++at) {
                *at -= optalign_trans[1];
                *at *= optalign_rotm;
                *at += optalign_trans[0];
            }
            for (unsigned int i=0; i<A_comp_list.size(); ++i) test_rmsd += get_square_distance(A_comp_list[i],B_comp_list[i]);
        } else for (bkn_vec jt=it->second.begin(); jt!=it->second.end(); ++jt) {
            test_rmsd += get_square_distance((*jt)->obj1->coord,(*jt)->obj2->coord);
        }
        if (test_rmsd < rmsd) {
            rmsd = test_rmsd;
            best_clique = &(it->second);
        }
    }

    ref_mol->has_stereo = has_stereo;

    if (best_clique == 0) {
        for (bkn_vec it=test_list.begin(); it!=test_list.end(); ++it) delete *it;
        return -2; // -2, wenn falsche Stereochemie
    }

    if (rmsd > 999999999998.) {
        for (bkn_vec it=test_list.begin(); it!=test_list.end(); ++it) delete *it;
        return -1.; //! -1 zeigt an, dass hier was nicht geklappt hat
    }

    rmsd /= mysolver.max_clique;

    if (with_hyd) {
        if (optalign) cerr << c_message<cWARNING>("MOLECULE::get_rmsd --> currently rmsd with hydrogens ")
                           << "is not available in combination with an optimal alignment" << endl;
        //! Jetzt die Wasserstoffe noch den gematchten Atomen zuordnen und den rmsd entsprechend modifizieren:
        int n_match = best_clique->size();
        rmsd *= n_match;
        for (bkn_vec it=best_clique->begin(); it!=best_clique->end(); ++it) {
            tr1::unordered_set<int> matched;
            for (atoms_vec at=(*it)->obj1->bonded_atoms.begin(); at != (*it)->obj1->bonded_atoms.end(); ++at) {
                if ((*at)->element == "H") {
                    float sd_tmp;
                    float sd_best = 999999999.;
                    int match_id = -1;
                    for (atoms_vec bt=(*it)->obj2->bonded_atoms.begin(); bt != (*it)->obj2->bonded_atoms.end(); ++bt) {
                        if ((*bt)->element != "H") continue;
                        if (matched.find((*bt)->intern_id) != matched.end()) continue;
                        sd_tmp = get_square_distance((*at)->coord,(*bt)->coord);
                        if (sd_tmp < sd_best) {
                            sd_best = sd_tmp;
                            match_id = (*bt)->intern_id;
                        }
                    }
                    if (match_id != -1) {

                        //cerr << "matched " << (*at)->intern_id << " with " << match_id << " d = " << sqrt(sd_best) << endl;

                        rmsd += sd_best;
                        ++n_match;
                        matched.insert(match_id);
                    }
                }
            }
        }
        rmsd /= n_match;
    }

    if (debug_mode) {
        for (bkn_vec it=best_clique->begin(); it!=best_clique->end(); ++it) {
            cerr << "compared  " << (*it)->obj1->name << " (" << (*it)->obj1->id << ")  with  "
                 << (*it)->obj2->name << " (" << (*it)->obj2->id << ")" << endl;
        }
    }

    //! 6.) Aufraeumen!:
    for (bkn_vec it=test_list.begin(); it!=test_list.end(); ++it) delete *it;

    release_SP_map(); // Die von der Referenz wird noch gebraucht

    return sqrt(rmsd);
}


template<class T>
float MOLECULE::get_bk_rmsd(stl_ptr<T> &ref_mol,bool with_hyd,bool debug_mode,bool optalign,bool recalc_hybrid) {
    //! 1.) Elemente und Konnektivitaeten beider Molekuele bestimmen
    if (recalc_hybrid) get_hybrid_only();
    //! Fuer die Referenz unbedingt vorher einmal aufrufen!!! (damits nur einmal aufgerufen wird)
//    ref_mol->get_hybrid_only();
    
    //!Neu: 04.11.2009:
    get_compare_numbers();
    for (atoms_vec at=ref_mol->atoms.begin(); at!=ref_mol->atoms.end(); ++at) {
        (*at)->get_compare_number();
    }
    
    //! 2.) alle moeglichen Atom-Atom-Paare erzeugen (nur gleiche Elemente zusammen):
    bkn_container test_list;
    for (atoms_vec at=atoms.begin(); at!=atoms.end(); ++at) {
        if ((*at)->element == "H") continue;
        for (atoms_vec bt=ref_mol->atoms.begin(); bt!=ref_mol->atoms.end(); ++bt) {
            if ((*bt)->element == "H") continue;
            if ((*at)->element != (*bt)->element) continue;
            
            if ((*at)->ext->compare_number != (*bt)->ext->compare_number) continue; //!Neu: 04.11.2009
            
            if ((*at)->ext->is_ring) {
                if (!((*bt)->ext->is_ring)) continue;
                if ((*at)->ext->is_planar_ring != (*bt)->ext->is_planar_ring) continue;
                if ((*at)->ext->is_aromatic != (*bt)->ext->is_aromatic) continue;
            } else if ((*bt)->ext->is_ring) continue;
            bkn_pointer new_node = new bkn_node(*at,*bt);
            test_list.push_back(new_node);
        }
    }
    
    //! 3.) Konnektivitaeten:
    //! Eine simple Variante die Kanten im assoziativen Graphen zu erzeugen waere die,
    //! jeweils einen Knoten (at1,bt1) mit einem Knoten (at2,bt2) zu verbinden, wenn eines
    //! der folgenden Kriterien erfuellt ist:
    //! - at1 ist direkt mit at2 verbunden UND bt1 ist direkt mit bt2 verbunden
    //! - at1 ist nicht direkt mit at2 verbunden UND bt1 ist nicht direkt mit bt2 verbunden
    //!
    //! Die Laufzeit des BK steigt exponentiell mit der Anzahl der Kanten im Graphen, weswegen
    //! schon bei Molekuelen mit mehr als 20 Atomen die Laufzeit erheblich wird mit den obigen
    //! einfachen Kriterien zum Kanten setzen. Eine andere Variante waere, dass dist(at1,at2)
    //! und dist(bt1,bt2) einen aehnlichen Wert haben. Dies ist aber fuer eine 3D-Substruktur-
    //! suche unvorteilhaft, weil sehr unterschiedliche Konformere so nicht als gleich erkannt werden.
    //! Eine exakte Loesung waere es nicht nur auf direkte Bindungen zu gucken, sondern die komplette
    //! Pfadlaenge zwischen at1 und at2. Das Problem hierbei ist, dass bei Beteiligung von Ringsystemen
    //! mehrere Pfade zwischen at1 und at2 moeglich sind. Fuer eine exakte Loesung muesste nur jeweils
    //! ein moeglicher Pfad(at1,at2) mit einem moeglichen Pfad(bt1,bt2) uebereinstimmen. Saemtliche Pfade
    //! in einem Graphen zu bestimmen ist aber leider ein NP-schweres Problem. Stattdessen nehme ich als
    //! nicht ganz exakte Loesung lieber nur den jeweils kuerzesten Pfad. Zusaetzlich bestimme ich noch
    //! eine Pfadsumme, deren Einzelwerte sich aus den Elementen der einzelnen Verbindungen ergeben:
    
    
//    cerr << "-> Ermittle Pfade:" << endl;

    if (ref_mol->SP_map == 0) ref_mol->calc_SP_map();
    if (SP_map == 0) calc_SP_map();

    //! 3.3.) Jetzt die Kanten erzeugen:
    int edges = 0;
    for (bkn_vec it=test_list.begin(); it!=test_list.end(); ++it) {
        for (bkn_vec jt=it+1; jt!=test_list.end(); ++jt) {
            if ((*it)->obj1.equal_addr((*jt)->obj1) || (*it)->obj2.equal_addr((*jt)->obj2)) continue;
            bool bonded = false;
            if ((*it)->obj1->intern_id > (*jt)->obj1->intern_id) {
                if ((*it)->obj2->intern_id > (*jt)->obj2->intern_id) {
                    if (SP_map[(*it)->obj1->intern_id][(*jt)->obj1->intern_id] == ref_mol->SP_map[(*it)->obj2->intern_id][(*jt)->obj2->intern_id]) {
                        if (SP_map[(*it)->obj1->intern_id][(*jt)->obj1->intern_id] != 0) bonded = true;
                        else if (ref_mol->SP_map[(*it)->obj2->intern_id][(*jt)->obj2->intern_id] == 0) bonded = true;
                    }
                } else {
                    if (SP_map[(*it)->obj1->intern_id][(*jt)->obj1->intern_id] == ref_mol->SP_map[(*jt)->obj2->intern_id][(*it)->obj2->intern_id]) {
                        if (SP_map[(*it)->obj1->intern_id][(*jt)->obj1->intern_id] != 0) bonded = true;
                        else if (ref_mol->SP_map[(*jt)->obj2->intern_id][(*it)->obj2->intern_id] == 0) bonded = true;
                    }
                }
            } else {
                if ((*it)->obj2->intern_id > (*jt)->obj2->intern_id) {
                    if (SP_map[(*jt)->obj1->intern_id][(*it)->obj1->intern_id] == ref_mol->SP_map[(*it)->obj2->intern_id][(*jt)->obj2->intern_id]) {
                       if (SP_map[(*jt)->obj1->intern_id][(*it)->obj1->intern_id] != 0) bonded = true;
                       else if (ref_mol->SP_map[(*it)->obj2->intern_id][(*jt)->obj2->intern_id] == 0) bonded = true;
                    }
                } else {
                    if (SP_map[(*jt)->obj1->intern_id][(*it)->obj1->intern_id] == ref_mol->SP_map[(*jt)->obj2->intern_id][(*it)->obj2->intern_id]) {
                        if (SP_map[(*jt)->obj1->intern_id][(*it)->obj1->intern_id] != 0) bonded = true;
                        else if (ref_mol->SP_map[(*jt)->obj2->intern_id][(*it)->obj2->intern_id] == 0) bonded = true;
                    }
                }
            }
            if (bonded) {
                (*it)->edges.insert(*jt);
                (*jt)->edges.insert(*it);
                ++edges;
            }
        }
    }
    if (debug_mode) cerr << "-> " << edges << " edges and " << test_list.size() << " nodes" << endl;

    //! 4.) Jetzt die Cliquen suchen:
    BK_SOLVER<stl_ptr<ATOM> > mysolver(test_list);
    mysolver.solve(); //! hier kann optional auch eine Mindestgroesse fuer die Cliquen mitgegeben werden
    //! => mysolver.cliques  ist eine multimap, die nun alle gefundenen Cliquen enthaelt (Cliquengroesse als Keys)
    //!    mysolver.max_clique  ist ein unsigned int der die groesste Clique angibt

    //! 5.) Es koennen mehrere groesste Cliquen gefunden worden sein => nun muss die Clique ermittelt werden, die
    //!     den besten RMSD liefert:
    float rmsd = 999999999999.;
    bkn_container * best_clique = 0;
    for (clique_map it=mysolver.cliques.begin(); it!=mysolver.cliques.end(); ++it) {
        if (it->first < mysolver.max_clique) continue;
        float test_rmsd = 0.;
        if (optalign) {
            vector<vec3d<float> > A_comp_list;
            vector<vec3d<float> > B_comp_list;
            for (bkn_vec jt=it->second.begin(); jt!=it->second.end(); ++jt) {
                A_comp_list.push_back((*jt)->obj1->coord);
                B_comp_list.push_back((*jt)->obj2->coord);
            }
            get_align_matrix(B_comp_list,A_comp_list,optalign_rotm,optalign_trans[0],optalign_trans[1]);
            for (vector<vec3d<float> >::iterator at=A_comp_list.begin(); at!=A_comp_list.end(); ++at) {
                *at -= optalign_trans[1];
                *at *= optalign_rotm;
                *at += optalign_trans[0];
            }
            for (unsigned int i=0; i<A_comp_list.size(); ++i) test_rmsd += get_square_distance(A_comp_list[i],B_comp_list[i]);
        } else for (bkn_vec jt=it->second.begin(); jt!=it->second.end(); ++jt) {
            test_rmsd += get_square_distance((*jt)->obj1->coord,(*jt)->obj2->coord);
        }
        if (test_rmsd < rmsd) {
            rmsd = test_rmsd;
            best_clique = &(it->second);
        }
    }

    if (rmsd > 999999999998.) {
        for (bkn_vec it=test_list.begin(); it!=test_list.end(); ++it) delete *it;
        return -1.; //! -1 zeigt an, dass hier was nicht geklappt hat
    }

    rmsd /= mysolver.max_clique;
    
    if (with_hyd) {
        if (optalign) cerr << c_message<cWARNING>("MOLECULE::get_bk_rmsd --> currently rmsd with hydrogens ")
                           << "is not available in combination with an optimal alignment" << endl;
        //! Jetzt die Wasserstoffe noch den gematchten Atomen zuordnen und den rmsd entsprechend modifizieren:
        int n_match = best_clique->size();
        rmsd *= n_match;
        for (bkn_vec it=best_clique->begin(); it!=best_clique->end(); ++it) {
            tr1::unordered_set<int> matched;
            for (atoms_vec at=(*it)->obj1->bonded_atoms.begin(); at != (*it)->obj1->bonded_atoms.end(); ++at) {
                if ((*at)->element == "H") {
                    float sd_tmp;
                    float sd_best = 999999999.;
                    int match_id = -1;
                    for (atoms_vec bt=(*it)->obj2->bonded_atoms.begin(); bt != (*it)->obj2->bonded_atoms.end(); ++bt) {
                        if ((*bt)->element != "H") continue;
                        if (matched.find((*bt)->intern_id) != matched.end()) continue;
                        sd_tmp = get_square_distance((*at)->coord,(*bt)->coord);
                        if (sd_tmp < sd_best) {
                            sd_best = sd_tmp;
                            match_id = (*bt)->intern_id;
                        }
                    }
                    if (match_id != -1) {
                        rmsd += sd_best;
                        ++n_match;
                        matched.insert(match_id);
                    }
                }
            }
        }
        rmsd /= n_match;
    }

    if (debug_mode) {
        for (bkn_vec it=best_clique->begin(); it!=best_clique->end(); ++it) {
            cerr << "compared  " << (*it)->obj1->name << " (" << (*it)->obj1->id << ")  with  "
                 << (*it)->obj2->name << " (" << (*it)->obj2->id << ")" << endl;
        }
    }
    
    //! 6.) Aufraeumen!:
    for (bkn_vec it=test_list.begin(); it!=test_list.end(); ++it) delete *it;

    release_SP_map(); // Die von der Referenz wird noch gebraucht

    return sqrt(rmsd);
}


#endif

