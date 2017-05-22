
//============================================================================
// molecule_GN.cpp -*- C++ -*-; representation of molecules and related stuff
//
// Copyright (C) 2006, 2007, 2008, 2009, 2010, 2011, 2012 Gerd Neudert
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


#include"molecule_GN.h"
#include"structure_GN.h"


class BOND_SORT {
public:
    bool operator()(stl_ptr<BOND> const& a,stl_ptr<BOND> const& b) {
        if (a->from->intern_id == b->from->intern_id) return (a->to->intern_id < b->to->intern_id);
        else return (a->from->intern_id < b->from->intern_id);
    }
} bond_sort;


RELIBASE_CENTER::RELIBASE_CENTER() {}

RELIBASE_CENTER::~RELIBASE_CENTER() {}



MOLECULE::MOLECULE() {
    freerot_bonds = -1;
    ele_already_set = false;
    has_stereo = false;
    atom_types_already_set = false;
    n_heavy_atoms = -1;
    compare_sum = -1;
    SP_map = 0;
}


MOLECULE::MOLECULE(MOLECULE const& lig): freerot_bonds(lig.freerot_bonds),n_heavy_atoms(lig.n_heavy_atoms),
                                         compare_sum(lig.compare_sum),verbosity(lig.verbosity),ele_already_set(lig.ele_already_set),
                                         atom_types_already_set(lig.atom_types_already_set),has_stereo(lig.has_stereo),
                                         name(lig.name),SP_map(lig.SP_map) {
    for (const_atoms_vec at=lig.atoms.begin(); at!=lig.atoms.end(); ++at) {
        stl_ptr<ATOM> atm(new ATOM(**at));
        atoms.push_back(atm);
    }
    
    for (const_bonds_vec bt=lig.bonds.begin(); bt!=lig.bonds.end(); ++bt) {
        stl_ptr<BOND> bnd(new BOND());
        bnd->id = (*bt)->id;
        bnd->free_rot = (*bt)->free_rot;
        bnd->type = (*bt)->type;
        for (atoms_vec at=atoms.begin(); at!=atoms.end(); ++at) {
            if ((*at)->intern_id == (*bt)->from->intern_id) bnd->from = (*at);
            else if ((*at)->intern_id == (*bt)->to->intern_id) bnd->to = (*at);
        }
        bonds.push_back(bnd);
    }
    
    for (map<int,stl_ptr<ATOM> >::const_iterator it=lig.sub_map.begin(); it!=lig.sub_map.end(); ++it) {
        for (atoms_vec at=atoms.begin(); at!=atoms.end(); ++at) {
            if (it->first == (*at)->intern_id) sub_map[it->first] = *at;
        }
    }
    
    for (const_comments_vec ct=lig.comments.begin(); ct!=lig.comments.end(); ++ct) {
        stl_ptr<COMMENT> cnt(new COMMENT());
        cnt->text = (*ct)->text;
        comments.push_back(cnt);
    }
    
    main_structure = lig.main_structure; //!neu: 11.10.07
    
    bool b_flag;
    
    for (bonds_vec it=bonds.begin(); it!=bonds.end();++it) {
        b_flag = false;
        for (atoms_vec wt=(*it)->from->bonded_atoms.begin(); wt!=(*it)->from->bonded_atoms.end();++wt) {
            if ((*it)->to->intern_id == (*wt)->intern_id) {
                b_flag = true;
                break;
            }
        }    
        if (!b_flag) {
            (*it)->from->bonded_atoms.push_back((*it)->to);
            (*it)->to->bonded_atoms.push_back((*it)->from);
        }
    }
}


MOLECULE::MOLECULE(LIGAND const& lig): freerot_bonds(lig.freerot_bonds),n_heavy_atoms(lig.n_heavy_atoms),
                                 compare_sum(lig.compare_sum),verbosity(lig.verbosity),ele_already_set(lig.ele_already_set),
                                 atom_types_already_set(lig.atom_types_already_set),has_stereo(lig.has_stereo),
                                 name(lig.name),SP_map(lig.SP_map) {
    for (const_atoms_vec at=lig.atoms.begin(); at!=lig.atoms.end(); ++at) {
        stl_ptr<ATOM> atm(new ATOM(**at));
        atoms.push_back(atm);
    }
    
    for (const_bonds_vec bt=lig.bonds.begin(); bt!=lig.bonds.end(); ++bt) {
        stl_ptr<BOND> bnd(new BOND());
        bnd->id = (*bt)->id;
        bnd->free_rot = (*bt)->free_rot;
        bnd->type = (*bt)->type;
        for (atoms_vec at=atoms.begin(); at!=atoms.end(); ++at) {
            if ((*at)->intern_id == (*bt)->from->intern_id) bnd->from = (*at);
            else if ((*at)->intern_id == (*bt)->to->intern_id) bnd->to = (*at);
        }
        bonds.push_back(bnd);
    }
    
    for (map<int,stl_ptr<ATOM> >::const_iterator it=lig.sub_map.begin(); it!=lig.sub_map.end(); ++it) {
        for (atoms_vec at=atoms.begin(); at!=atoms.end(); ++at) {
            if (it->first == (*at)->intern_id) sub_map[it->first] = *at;
        }
    }
    
    for (const_comments_vec ct=lig.comments.begin(); ct!=lig.comments.end(); ++ct) {
        stl_ptr<COMMENT> cnt(new COMMENT());
        cnt->text = (*ct)->text;
        comments.push_back(cnt);
    }
    
    main_structure = lig.main_structure; //!neu: 11.10.07
    
    bool b_flag;
    
    for (bonds_vec it=bonds.begin(); it!=bonds.end();++it) {
        b_flag = false;
        for (atoms_vec wt=(*it)->from->bonded_atoms.begin(); wt!=(*it)->from->bonded_atoms.end();++wt) {
            if ((*it)->to->intern_id == (*wt)->intern_id) {
                b_flag = true;
                break;
            }
        }    
        if (!b_flag) {
            (*it)->from->bonded_atoms.push_back((*it)->to);
            (*it)->to->bonded_atoms.push_back((*it)->from);
        }
    }
}


bool MOLECULE::copy_without_coords(LIGAND* lig) {
    if (atoms.size() != lig->atoms.size()) return false;
    ele_already_set = lig->ele_already_set;
    atom_types_already_set = lig->atom_types_already_set;
    has_stereo = lig->has_stereo;
    freerot_bonds = lig->freerot_bonds;
    verbosity = lig->verbosity;
    n_heavy_atoms = lig->n_heavy_atoms;
    compare_sum = lig->compare_sum;
    name = lig->name;
    SP_map = lig->SP_map;
    for (unsigned int i=0; i<lig->atoms.size(); ++i) {
        atoms[i]->bonded_atoms.clear();
        atoms[i]->id = lig->atoms[i]->id;
        atoms[i]->intern_id = lig->atoms[i]->intern_id;
        atoms[i]->name = lig->atoms[i]->name;
        atoms[i]->type = lig->atoms[i]->type;
        atoms[i]->sybyl_type = lig->atoms[i]->sybyl_type;
        atoms[i]->intern_type = lig->atoms[i]->intern_type;
        atoms[i]->res_name = lig->atoms[i]->res_name;
        atoms[i]->res_number = lig->atoms[i]->res_number;
        atoms[i]->element = lig->atoms[i]->element;
        atoms[i]->charge = lig->atoms[i]->charge;
        atoms[i]->bond_ind = lig->atoms[i]->bond_ind;
        atoms[i]->full_res_name = lig->atoms[i]->full_res_name;
        if (!atoms[i]->ext.zero()) atoms[i]->ext.kill();
        if (!lig->atoms[i]->ext.zero()) atoms[i]->ext = new ATOM_EXT(*(lig->atoms[i]->ext));
    }
    
    for (bonds_vec it=bonds.begin(); it!=bonds.end();++it) {
        it->kill();
    }
    
    bonds.clear();
    for (bonds_vec bt=lig->bonds.begin(); bt!=lig->bonds.end(); ++bt) {
        stl_ptr<BOND> bnd(new BOND());
        bnd->id = (*bt)->id;
        bnd->free_rot = (*bt)->free_rot;
        bnd->type = (*bt)->type;
        for (atoms_vec at=atoms.begin(); at!=atoms.end(); ++at) {
            if ((*at)->intern_id == (*bt)->from->intern_id) bnd->from = (*at);
            else if ((*at)->intern_id == (*bt)->to->intern_id) bnd->to = (*at);
        }
        bonds.push_back(bnd);
    }
    
    bool b_flag;
    for (bonds_vec it=bonds.begin(); it!=bonds.end();++it) {
        b_flag = false;
        for (atoms_vec wt=(*it)->from->bonded_atoms.begin(); wt!=(*it)->from->bonded_atoms.end();++wt) {
            if ((*it)->to->intern_id == (*wt)->intern_id) {
                b_flag = true;
                break;
            }
        }    
        if (!b_flag) {
            (*it)->from->bonded_atoms.push_back((*it)->to);
            (*it)->to->bonded_atoms.push_back((*it)->from);
        }
    }
    
    for (map<int,stl_ptr<ATOM> >::iterator it=lig->sub_map.begin(); it!=lig->sub_map.end(); ++it) {
        for (atoms_vec at=atoms.begin(); at!=atoms.end(); ++at) {
            if (it->first == (*at)->intern_id) sub_map[it->first] = *at;
        }
    }
    
    for (rings_vec it=rings.begin(); it!=rings.end(); ++it) {
        it->kill();
    }
    rings.clear();
    
    for (rings_vec it=lig->rings.begin(); it!=lig->rings.end(); ++it) {
        stl_ptr<RING> rng(new RING());
        rng->n_hetero = (*it)->n_hetero;
        rng->n_members = (*it)->n_members;
        rng->is_planar = (*it)->is_planar;
        rng->is_aromatic = (*it)->is_aromatic;
        rng->is_pos = (*it)->is_pos;
        rng->is_prot = (*it)->is_prot;
        for (atoms_vec rt=(*it)->ring_atoms.begin(); rt!=(*it)->ring_atoms.end(); ++rt) {
            for (atoms_vec at=atoms.begin(); at!=atoms.end(); ++at) {
                if ((*at)->intern_id == (*rt)->intern_id) {
                    rng->ring_atoms.push_back(*at);
                    break;
                }
            }
        }
        rings.push_back(rng);
    }

    return true;
}


MOLECULE::~MOLECULE() {
    release_SP_map();
    for (bonds_vec it=bonds.begin(); it!=bonds.end();++it) {
        it->kill();
    }
    
    for (atoms_vec it=atoms.begin(); it!=atoms.end();++it) {
        it->kill();
    }
    
    for (rings_vec it=rings.begin(); it!=rings.end();++it) {
        it->kill();
    }
    
    for (comments_vec it=comments.begin(); it!=comments.end();++it) {
        it->kill();
    }
}


void MOLECULE::get_compare_numbers() {
    for (atoms_vec at=atoms.begin(); at!=atoms.end(); ++at) {
        (*at)->get_compare_number();
    }
    for (atoms_vec at=atoms.begin(); at!=atoms.end(); ++at) {
        (*at)->get_compare_number2();
    }
}


int MOLECULE::get_compare_sum() {
    if (compare_sum > -1) return compare_sum;
    compare_sum = 0;
    for (atoms_vec at=atoms.begin(); at!=atoms.end(); ++at) {
        if ((*at)->ext.zero()) return -1;
        if ((*at)->sybyl_type[0] == 'H') continue;
        compare_sum += (*at)->ext->compare_number;
    }
    return compare_sum;
}


int MOLECULE::get_freerot_bonds() {
    if (freerot_bonds > -1) return freerot_bonds;
    int res = 0;
    for (bonds_vec it=bonds.begin(); it!=bonds.end();++it) {
        if ((*it)->free_rot) ++res;
        
    }
    freerot_bonds = res;
    return res;
}


void MOLECULE::dlg2mol2(int mode,const char *def_file,bool get_bonds,int verb,bool kill_ext,bool fill_X,const char *alt_def_file,
                        int max_ring_members,bool no_free_rot_planar,LIGAND* at_ref_mol,LIGAND* at_ref_mol2,bool prot_acids,
                        bool prot_guanidin,bool prot_amidin,bool prot_amin,bool prot_phosphate,bool prot_sulfate,
                        bool kekulize_aromatics,bool kekulize_charged,bool allow_charged_aromatics) {
    get_atom_typing(mode,false,def_file,get_bonds,verb,kill_ext,fill_X,alt_def_file,max_ring_members,no_free_rot_planar,at_ref_mol,
                    at_ref_mol2,prot_acids,prot_guanidin,prot_amidin,prot_amin,prot_phosphate,prot_sulfate,
                    kekulize_aromatics,kekulize_charged,allow_charged_aromatics);
}


void MOLECULE::dlg2mol2(stl_ptr<LIGAND> &ref_mol) {
    //!Fuer ref_mol vorher die bonded_atoms erzeugen!!!
    //zunaechst eigene Typen und Bonds setzen (Elemente sind schon gesetzt!):
    
    get_atom_typing(1,false,"X",true,1,false,false,
                    "X",10,true,0,0,prot_acids_default,
                    prot_guanidin_default,prot_amidin_default,prot_amin_default,prot_phosphate_default,
                    prot_sulfate_default,kekulize_aromatics_default,kekulize_charged_default,allow_charged_aromatics_default);

    get_compare_numbers();
    ref_mol->get_compare_numbers();

    tr1::unordered_map<int,tr1::unordered_set<int> > known;
    tr1::unordered_map<int,tr1::unordered_set<int> > negknown;
    
    //Jetzt gucken was zugeordnet werden kann:
    for (atoms_vec at=atoms.begin(); at!=atoms.end(); ++at) {
        bool found_match = false;
        
        for (atoms_vec bt=ref_mol->atoms.begin(); bt!=ref_mol->atoms.end(); ++bt) {
            
            if ((*at)->element == (*bt)->element) {
                
                vector<stl_ptr<ATOM> > prevs1;
                vector<stl_ptr<ATOM> > prevs2;
                if ((*at)->is_equal(**bt,prevs1,prevs2,known,negknown)) {
                    //alles relevante ausser den Koordinaten von der Referenz uebernehmen:
                    (*at)->name = (*bt)->name;
                    (*at)->type = (*bt)->type;
                    (*at)->sybyl_type = (*bt)->sybyl_type;
                    (*at)->res_name = (*bt)->res_name;
                    (*at)->res_number = (*bt)->res_number;
                    (*at)->charge = (*bt)->charge;
                    (*at)->full_res_name = (*bt)->full_res_name;
                    found_match = true;
                    (*at)->bond_ind = (*bt)->intern_id; //!um unten die bonds vergleichen zu koennen
                    break;
                }
            }
        }
        if (!found_match) if (verbosity) {
            cerr << c_message<cWARNING>("MOLECULE::dlg2mol2 --> found no match for ") << (*at)->name 
                 << " in the reference structure! Using own atom typing!" << endl;
        }
    }
    
    for (bonds_vec it=bonds.begin(); it!=bonds.end();++it) {
        for (bonds_vec jt=ref_mol->bonds.begin(); jt!=ref_mol->bonds.end(); ++jt) {
            if (((*it)->from->element == (*jt)->from->element) && ((*it)->to->element == (*jt)->to->element)) {
                if (((*it)->from->bond_ind != (*jt)->from->intern_id) || ((*it)->to->bond_ind != (*jt)->to->intern_id)) continue;
                (*it)->type = (*jt)->type;
            } else if (((*it)->from->element == (*jt)->to->element) && ((*it)->to->element == (*jt)->from->element)) {
                if (((*it)->from->bond_ind != (*jt)->to->intern_id) || ((*it)->to->bond_ind != (*jt)->from->intern_id)) continue;
                (*it)->type = (*jt)->type;
            }
        }
    }
    
    for (atoms_vec at=atoms.begin(); at!=atoms.end(); ++at) {
        (*at)->remove_ext();
    }
}


bool MOLECULE::check_equality(stl_ptr<LIGAND> &ref_mol,bool recalc_compare_numbers,bool stereo) {
    if (recalc_compare_numbers) {
        get_compare_numbers();
        for (atoms_vec at=ref_mol->atoms.begin(); at!=ref_mol->atoms.end(); ++at) {
            (*at)->get_compare_number();
        }
        for (atoms_vec at=ref_mol->atoms.begin(); at!=ref_mol->atoms.end(); ++at) {
            (*at)->get_compare_number2();
        }
    }

    if (get_compare_sum() != ref_mol->get_compare_sum()) return false;
    
    tr1::unordered_map<int,tr1::unordered_set<int> > known;
    tr1::unordered_map<int,tr1::unordered_set<int> > negknown;

    bool equality = false;
    for (atoms_vec at=atoms.begin(); at!=atoms.end(); ++at) {
        if ((*at)->element == "H") continue;
        for (atoms_vec bt=ref_mol->atoms.begin(); bt!=ref_mol->atoms.end(); ++bt) {
            if ((*at)->element == "H") continue;
            if ((*at)->element == (*bt)->element) {
                vector<stl_ptr<ATOM> > prevs1;
                vector<stl_ptr<ATOM> > prevs2;
                if ((*at)->is_equal(**bt,prevs1,prevs2,known,negknown)) {
                    equality = true;
                    break;
                }
            }
        }
        break;
    }

    if (equality && stereo && ref_mol->has_stereo) {
        if (get_rmsd(ref_mol,false,false) > -0.5) return true;
        else return false;
    } else return equality;
}


int MOLECULE::get_closest_ring_id(stl_ptr<ATOM> const& from) {
    //!ACHTUNG: diese Funktion ist veraltet! => Bei Gelegenheit mal ueberarbeiten mit Hilfe von
    //!         ATOM::is_in_same_ring
    
    //! Nur zwei Konnektivitaetsebenen weit pruefen (reicht fuer den Zweck)
    for (atoms_vec at=from->bonded_atoms.begin(); at!=from->bonded_atoms.end(); ++at) {
        if ((*at)->ext->is_ring) {
            for (rings_vec rt=rings.begin(); rt!=rings.end(); ++rt) {
                for (atoms_vec rat=(*rt)->ring_atoms.begin(); rat!=(*rt)->ring_atoms.end(); ++rat) {
                    if (*at == *rat) return (*rt)->id;
                }
            }
        }
    }
    for (atoms_vec at=from->bonded_atoms.begin(); at!=from->bonded_atoms.end(); ++at) {
        for (atoms_vec bt=(*at)->bonded_atoms.begin(); bt!=(*at)->bonded_atoms.end(); ++bt) {
            if (*bt == from) continue;
            if ((*bt)->ext->is_ring) {
                for (rings_vec rt=rings.begin(); rt!=rings.end(); ++rt) {
                    for (atoms_vec rat=(*rt)->ring_atoms.begin(); rat!=(*rt)->ring_atoms.end(); ++rat) {
                        if (*bt == *rat) return (*rt)->id;
                    }
                }
            }
        }
    }
    return -1;
}


void MOLECULE::release_SP_map() {
    if (SP_map != 0) {
        for (unsigned int i=2; i<(atoms.size()+1); ++i) delete[] SP_map[i];
        delete[] SP_map;
        SP_map = 0;
    }
}

void MOLECULE::calc_SP_map() {
    release_SP_map();

    SP_map = new int*[atoms.size()+1];
    for (unsigned int i=2; i<=atoms.size(); ++i) {
        SP_map[i] = new int[i];
        for (unsigned int j=0; j<i; ++j) SP_map[i][j] = 0;
    }

    for (atoms_vec at=atoms.begin(); at!=atoms.end(); ++at) {
        if ((*at)->element == "H") continue;
        queue<vector<stl_ptr<ATOM> > > bat;
        tr1::unordered_set<int> visited;
        visited.insert((*at)->intern_id);
        vector<stl_ptr<ATOM> > vbat;
        for (atoms_vec ct=(*at)->bonded_atoms.begin(); ct!=(*at)->bonded_atoms.end(); ++ct) {
            if ((*ct)->element == "H") continue;
            vbat.push_back(*ct);
        }
        bat.push(vbat);
        int pl = 1;
        while (!bat.empty()) {
            vector<stl_ptr<ATOM> > vb;
            for (atoms_vec avbat=bat.front().begin(); avbat!=bat.front().end(); ++avbat) {
                if ((*avbat)->element == "H") continue;
                if (visited.find((*avbat)->intern_id) != visited.end()) continue;
                visited.insert((*avbat)->intern_id);
                if ((*at)->intern_id > (*avbat)->intern_id) SP_map[(*at)->intern_id][(*avbat)->intern_id] = pl;
                else SP_map[(*avbat)->intern_id][(*at)->intern_id] = pl;
                for (atoms_vec gt=(*avbat)->bonded_atoms.begin(); gt!=(*avbat)->bonded_atoms.end(); ++gt) {
                    if ((*gt)->element == "H") continue;
                    if (visited.find((*gt)->intern_id) != visited.end()) continue;
                    vb.push_back(*gt);
                }
            }
            if (vb.size() > 0) bat.push(vb);
            ++pl;
            bat.pop();
        }
    }
}


int MOLECULE::has_substructure(stl_ptr<LIGAND> &ref_mol,bool consider_hybrid,bool type_based,bool own_types) {
    unsigned int min_size = 0;
    for (atoms_vec bt=ref_mol->atoms.begin(); bt!=ref_mol->atoms.end(); ++bt) {
            if ((*bt)->element == "H") continue;
            ++min_size;
    }
    
    bkn_container test_list;
    
    if (type_based) for (atoms_vec at=atoms.begin(); at!=atoms.end(); ++at) {
        if ((*at)->element == "H") continue;
        for (atoms_vec bt=ref_mol->atoms.begin(); bt!=ref_mol->atoms.end(); ++bt) {
            if ((*bt)->element == "H") continue;
            if ((*at)->sybyl_type != (*bt)->sybyl_type) {
                if (own_types) {
                    if ((*at)->sybyl_type != "X" && (*bt)->sybyl_type != "X") continue;
                } else continue;
            }
            if ((*at)->ext->is_ring) {
                if (!((*bt)->ext->is_ring)) continue;
                if ((*at)->ext->is_planar_ring != (*bt)->ext->is_planar_ring) continue;
                if ((*at)->ext->is_aromatic != (*bt)->ext->is_aromatic) continue;
            }
            bkn_pointer new_node = new bkn_node(*at,*bt);
            test_list.push_back(new_node);
        }
    } else for (atoms_vec at=atoms.begin(); at!=atoms.end(); ++at) {
        if ((*at)->element == "H") continue;
        for (atoms_vec bt=ref_mol->atoms.begin(); bt!=ref_mol->atoms.end(); ++bt) {
            if ((*bt)->element == "H") continue;
            if ((*at)->element != (*bt)->element) continue;
        //    if (consider_hybrid && ((*at)->ext->hybridization != (*bt)->ext->hybridization)) continue;
            if (consider_hybrid) {
                if (((*bt)->charge != "free_hyb" && (*bt)->charge != "exact_bonded") && ((*at)->ext->hybridization != (*bt)->ext->hybridization)) continue;
            //    if (((*bt)->charge != "free_hyb" && (*bt)->charge != "exact_match") && ((*at)->ext->hybridization != (*bt)->ext->hybridization)) continue;
            }
            if ((*bt)->charge == "exact_match" || (*bt)->charge == "exact_bonded") {
        //    if ((*bt)->charge == "exact_match") {
                if ((*at)->ext->n_heavy_bonded != (*bt)->ext->n_heavy_bonded) continue;
            }

            //---------------------------------------
            // 19.01.2011: GN:  at und bt vertauscht:
            if ((*bt)->ext->is_ring) {
                if (!((*at)->ext->is_ring)) {
                    if ((*bt)->ext->is_planar_ring) continue;
                }
                if ((*bt)->ext->is_planar_ring != (*at)->ext->is_planar_ring) continue;
                if ((*bt)->ext->is_aromatic != (*at)->ext->is_aromatic) continue;
            }
            //---------------------------------------

            bkn_pointer new_node = new bkn_node(*at,*bt);
            test_list.push_back(new_node);
        }
    }

    if (ref_mol->SP_map == 0) ref_mol->calc_SP_map();
    if (SP_map == 0) calc_SP_map();

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
            }
        }
    }

    BK_SOLVER<stl_ptr<ATOM> > mysolver(test_list);
    mysolver.solve(min_size);
    
    //!-----------------------------------------------------------------------------
    //! DEBUG:
    /*
    for (clique_map jt=mysolver.cliques.begin(); jt!=mysolver.cliques.end(); ++jt) {
        cerr << endl << "=> new clique" << endl;
        for (bkn_vec it=jt->second.begin(); it!=jt->second.end(); ++it) {
            cerr << "compared  " << (*it)->obj1->name << " (" << (*it)->obj1->id << ")  with  "
                << (*it)->obj2->name << " (" << (*it)->obj2->id << ")" << endl;
        //    if ((*it)->obj2->id == 4) cerr << (*it)->obj2 << endl;
        }
    }
    */
    //!-----------------------------------------------------------------------------
    
    /*
    if (mysolver.cliques.size() == 0) {
        for (bkn_vec it=test_list.begin(); it!=test_list.end(); ++it) delete *it;
        return false;
    } else {
        for (bkn_vec it=test_list.begin(); it!=test_list.end(); ++it) delete *it;
        return true;
    }
    */

    for (bkn_vec it=test_list.begin(); it!=test_list.end(); ++it) delete *it;

    release_SP_map(); // Die von der Referenz wird noch gebraucht

    return mysolver.cliques.size();
}


void MOLECULE::opt_align() {
    for (atoms_vec at=atoms.begin(); at!=atoms.end(); ++at) {
        (*at)->coord -= optalign_trans[1];
        (*at)->coord *= optalign_rotm;
        (*at)->coord += optalign_trans[0];
    }
}

bool MOLECULE::check_trigo(stl_ptr<ATOM> const& base,float lim) {
    if (base->bonded_atoms.size() < 3) return true;

    if (lim < 0.9) {
        if (base->bonded_atoms[0]->element == "S"
            || base->bonded_atoms[1]->element == "S"
            || base->bonded_atoms[2]->element == "S") lim = 0.9; //!neu: 06.09.06
    }
    if (base->get_trigo_angle() < lim) return true;
    return false;

//    //Spatprodukt muss unter einem bestimmten threshold liegen:
//    float val;
//    vec3d<float> v1 = base->bonded_atoms[0]->coord; v1 -= base->coord;
//    vec3d<float> v2 = base->bonded_atoms[1]->coord; v2 -= base->coord;
//    vec3d<float> v3 = base->bonded_atoms[2]->coord; v3 -= base->coord;
//    v1 *= v2;
//    val = v1.skalar_product(v3);
//
//    if (lim < 0.9) {
//        if (base->bonded_atoms[0]->element == "S"
//            || base->bonded_atoms[1]->element == "S"
//            || base->bonded_atoms[2]->element == "S") lim = 0.9; //!neu: 06.09.06
//    }
//
//    if (fabs(val) < lim) return true;
//    return false;
}

void MOLECULE::delete_worst_bonded(stl_ptr<ATOM> const& base) {
    stl_ptr<ATOM> wa;
    float wp = 1.1;
    for (atoms_vec bt=base->bonded_atoms.begin(); bt!=base->bonded_atoms.end(); ++bt) {
        float mod = 1.;
        if ((*bt)->bonded_atoms.size() > 4) {
            mod = 0.5;
        } else if ((*bt)->element == "O" && (*bt)->bonded_atoms.size() > 2) mod = 0.5;

        float dst = get_distance(base->coord,(*bt)->coord);
        float bp = EDnew::EDAccess::get_max_probability(base->element,(*bt)->element,dst);

        bp *= mod;

        if (bp < wp) {
            wp = bp;
            wa = *bt;
        }
    }

    //!wa sollte jetzt das Atom mit der geringsten Bindungswahrscheinlichkeit sein
    for (atoms_vec bt=base->bonded_atoms.begin(); bt!=base->bonded_atoms.end(); ++bt) {
        if (*bt == wa) {
            base->bonded_atoms.erase(bt);
            if (wa->element != "H") base->ext->n_heavy_bonded--;
            break;
        }
    }
    for (atoms_vec bt=wa->bonded_atoms.begin(); bt!=wa->bonded_atoms.end(); ++bt) {
        if (*bt == base) {
            wa->bonded_atoms.erase(bt);
            if (base->element != "H") wa->ext->n_heavy_bonded--;
            break;
        }
    }

    if (verbosity > 1) {
        cerr << c_message<cWARNING>("MOLECULE::get_hybridizations --> close contact warning for ") << name << ": " << base->name << "(" << base->res_name
             << base->res_number << ") and " << wa->name << " (" << wa->res_name << wa->res_number << ")" << endl;
        cerr << "                                          (disconnecting because valence is exceeded)" << endl;
    }
}

bool MOLECULE::is_linear(stl_ptr<ATOM> const& base) {
    if (base->get_smallest_bond_angle() < 2.7925) return false; // ca. 160° Grad
    return true;
}

bool MOLECULE::check_db_tors(stl_ptr<ATOM> const& base,float const& lim) {
    //!liefert true, wenn die Torsion   'atom1---base===atom3---atom4' oder
    //! die Torsion              atom0---atom1===base---atom3          planar ist
    //! Je nachdem, welche DB eine hoehere Wahrscheinlichkeit hat
    stl_ptr<ATOM> p1;
    stl_ptr<ATOM> p3;
    stl_ptr<ATOM> p4;

    //!zunaechst den DB-Partner (p3) bestimmen:
    base->ext->sp2 = 0.;
    float sd = 0.;
    for (atoms_vec bt=base->bonded_atoms.begin(); bt!=base->bonded_atoms.end(); ++bt) {
        if ((*bt)->element == "H") continue;

    //    if ((*bt)->ext->hybridization == 3 && !(*bt)->ext->sec_check) continue;

        float dst = get_distance(base->coord,(*bt)->coord); //Abstand der Atome

        float prob2 = EDnew::EDAccess::get_probability(base->element,(*bt)->element,2,dst);

        if (prob2 > base->ext->sp2 || base->ext->sp2 < 0.00001 || (dst < sd && fabs(prob2-base->ext->sp2) < 0.01)) {
            base->ext->sp2 = prob2;
            p3 = *bt;
            sd = dst; //!03.11.09, weil fuer 1.3A und 1.4A gleiche Wahrscheinlichkeit
        } else if (base->ext->sp2 > 0.001) {
            if (p3->ext->hybridization == 3 && (*bt)->ext->hybridization == 2) {
                base->ext->sp2 = prob2;
                p3 = *bt;
                sd = dst; //!03.11.09, weil fuer 1.3A und 1.4A gleiche Wahrscheinlichkeit
            }
        }
    }

    //!jetzt ein moegliches p4:
    if (p3->ext->n_heavy_bonded < 2) return true; //! => eventuell falsch positiv
    for (atoms_vec bt=p3->bonded_atoms.begin(); bt!=p3->bonded_atoms.end(); ++bt) {
        if (((*bt)->intern_id != base->intern_id) && ((*bt)->element != "H")) {
            p4 = *bt;
        }
    }

    //!jetzt p1:
    for (atoms_vec bt=base->bonded_atoms.begin(); bt!=base->bonded_atoms.end(); ++bt) {
        if (((*bt)->intern_id != p3->intern_id) && ((*bt)->element != "H")) {
            p1 = *bt;
        }
    }

    float val = dihedral(p1->coord,base->coord,p3->coord,p4->coord);

    if ((val < lim) || (val > (3.14159-lim))) return true; // bis lim zulassen //bis ca. 20,6 Grad Torsion zulassen

    return false;
}

bool MOLECULE::check_planar_tors(stl_ptr<ATOM> const& from,stl_ptr<ATOM> const& to,float const& lim) {
    //!liefert true, wenn die Torsion   'atom1---from---to---atom4' planar ist
    if (from->ext->n_heavy_bonded < 2 || to->ext->n_heavy_bonded < 2) return true;
    stl_ptr<ATOM> p1;
    stl_ptr<ATOM> p4;
    
    //!jetzt ein moegliches p4:
    for (atoms_vec bt=to->bonded_atoms.begin(); bt!=to->bonded_atoms.end(); ++bt) {
        if (((*bt)->intern_id != from->intern_id) && ((*bt)->element != "H")) {
            p4 = *bt;
            break;
        }
    }
    
    //!jetzt p1:
    for (atoms_vec bt=from->bonded_atoms.begin(); bt!=from->bonded_atoms.end(); ++bt) {
        if (((*bt)->intern_id != to->intern_id) && ((*bt)->element != "H")) {
            p1 = *bt;
            break;
        }
    }
    
    float val = dihedral(p1->coord,from->coord,to->coord,p4->coord);
    
    if ((val < lim) || (val > (3.14159-lim))) return true; // bis lim zulassen //bis ca. 20,6 Grad Torsion zulassen
    
    return false;
}

void MOLECULE::get_elements(bool reset_names) {
    if (ele_already_set) return;
    for (atoms_vec at=atoms.begin(); at!=atoms.end(); ++at) {
        (*at)->get_element(reset_names);
    }
    ele_already_set = true;
}


//!====================================================================================================
//! ab hier die neue Atomtypzuweisung
//!====================================================================================================

void MOLECULE::get_connections() {
    EDnew::initialize();
    float dst;
    string key;

    OCTREE<ATOM>* pwp = new OCTREE<ATOM>(atoms,2.5);
    for (atoms_vec at=atoms.begin(); at!=atoms.end(); ++at) {
        for (vector<ATOM*>::iterator bt=pwp->begin(2.,(*at)->coord); bt!=pwp->end(); ++bt) {
            if ((*bt)->intern_id <= (*at)->intern_id) continue;

            dst = (*at)->coord[0] - (*bt)->coord[0]; dst *= dst;
            if (dst > 6.75) continue;
            dst += ((*at)->coord[1] - (*bt)->coord[1]) * ((*at)->coord[1] - (*bt)->coord[1]);
            if (dst > 6.75) continue;
            dst += ((*at)->coord[2] - (*bt)->coord[2]) * ((*at)->coord[2] - (*bt)->coord[2]);
            if (dst > 6.75) continue;
            if ((*at)->element == "H" || (*bt)->element == "H") {
                if (dst > 1.8225) continue;
            }

            if (!EDnew::EDAccess::connection_is_possible((*at)->element,(*bt)->element)) continue;

            dst = sqrt(dst);

            if (dst < EDnew::ED_min_dist || dst > EDnew::ED_max_dist) continue;

            float prob[5] = {EDnew::EDAccess::get_probability(0,dst),  // ar
                             EDnew::EDAccess::get_probability(1,dst),  // 1
                             EDnew::EDAccess::get_probability(2,dst),  // 2
                             EDnew::EDAccess::get_probability(3,dst)}; // 3

            if (prob[1] > (*at)->ext->sp3) (*at)->ext->sp3 = prob[1];
            if (prob[1] > (*bt)->ext->sp3) (*bt)->ext->sp3 = prob[1];
            if (prob[2] > (*at)->ext->sp2) (*at)->ext->sp2 = prob[2];
            if (prob[2] > (*bt)->ext->sp2) (*bt)->ext->sp2 = prob[2];
            if (prob[0] > (*at)->ext->sp2) (*at)->ext->sp2 = prob[0];
            if (prob[0] > (*bt)->ext->sp2) (*bt)->ext->sp2 = prob[0];
            if (prob[3] > (*at)->ext->sp) (*at)->ext->sp = prob[3];
            if (prob[3] > (*bt)->ext->sp) (*bt)->ext->sp = prob[3];

            if (prob[0] > 0. || prob[1] > 0. || prob[2] > 0. || prob[3] > 0.) {
                (*at)->bonded_atoms.push_back(*bt);
                (*bt)->bonded_atoms.push_back(*at);
            }
        }
    }
    delete pwp;


    for (atoms_vec at=atoms.begin(); at!=atoms.end(); ++at) {
        if ((*at)->element == "H") {
            if ((*at)->bonded_atoms.size() == 0 && (*at)->ext->pos_metal) {
                (*at)->type = 3;
                (*at)->element = "Hg";
                continue;
            } else while ((*at)->bonded_atoms.size() > 1) { //!H mit mehr als einem gebundenen Atom
                delete_worst_bonded(*at);
                continue;
            }
        } else if ((*at)->element == "C" || (*at)->element == "N" || (*at)->element == "Si" || 
                   (*at)->element == "Se" || (*at)->element == "As") {
            if ((*at)->bonded_atoms.size() > 4) {
                do {
                    delete_worst_bonded(*at);
                } while ((*at)->bonded_atoms.size() > 4);
            } else if ((*at)->bonded_atoms.size() == 0) {
                if ((*at)->ext->pos_metal) {
                    if ((*at)->element == "C") {
                        (*at)->type = 3;
                        (*at)->element = "Ca";
                    } else if ((*at)->element == "N") {
                        (*at)->type = 3;
                        if ((*at)->name.size() > 1) {
                            if ((*at)->name[1] ==  'I' || (*at)->name[1] ==  'i') (*at)->element = "Ni";
                            else if ((*at)->name[1] ==  'A' || (*at)->name[1] ==  'a') (*at)->element = "Na";
                        }
                    }
                }
            }
        } else if ((*at)->element == "S") {
            if ((*at)->bonded_atoms.size() > 6) {
                do {
                    delete_worst_bonded(*at);
                } while ((*at)->bonded_atoms.size() > 6);
            }
        } else if ((*at)->element == "P") {
            if ((*at)->bonded_atoms.size() > 5) {
                do {
                    delete_worst_bonded(*at);
                } while ((*at)->bonded_atoms.size() > 5);
            }
        } else if ((*at)->element == "O") {
            while ((*at)->bonded_atoms.size() > 2) delete_worst_bonded(*at);
        } else if ((*at)->element == "B") {
            while ((*at)->bonded_atoms.size() > 3) delete_worst_bonded(*at);
        } else while ((*at)->bonded_atoms.size() > 1) delete_worst_bonded(*at);
    }
    for (atoms_vec at=atoms.begin(); at!=atoms.end(); ++at) {
        (*at)->ext->n_heavy_bonded = 0;
        for (atoms_vec bt=(*at)->bonded_atoms.begin(); bt!=(*at)->bonded_atoms.end(); ++bt) {
            if ((*bt)->element != "H") (*at)->ext->n_heavy_bonded++;
        }
    }


    //!Jetzt die initialen Hybridisierungen festlegen:
    for (atoms_vec at=atoms.begin(); at!=atoms.end(); ++at) {
        (*at)->ext->hybridization = 0;
        (*at)->ext->check_hyb = false;  //Flag fuer Plausibilitaetstest
        (*at)->ext->sec_check = false;

        if ((*at)->type == 4) { //!Wasser
            (*at)->ext->hybridization = 3;
            continue;
        }

        if ((*at)->ext->sp > 0. && ((*at)->ext->sp2 / (*at)->ext->sp) < 50.) (*at)->ext->hybridization = 1;
        else if ((*at)->ext->n_heavy_bonded == 2) {
            if (is_linear(*at)) (*at)->ext->hybridization = 1;
            else if ((*at)->ext->sp2 > 0. && ((*at)->ext->sp3 / (*at)->ext->sp2) < 50.) (*at)->ext->hybridization = 2;
            else (*at)->ext->hybridization = 3;
        } else if ((*at)->ext->sp2 > 0. && ((*at)->ext->sp3 / (*at)->ext->sp2) < 50.) (*at)->ext->hybridization = 2;
        else (*at)->ext->hybridization = 3;
    }
}


void MOLECULE::get_hybridization(bool const& get_bonds,int const& max_ring_members,bool const& allow_protonated,bool const& kekulize_aromatics) {
    //! Relevant Rings bestimmen:
    get_rings(max_ring_members);

    //! Ringe auf Planaritaet pruefen:
    check_ring_plan();

//    for (atoms_vec at=atoms.begin(); at!=atoms.end(); ++at) {
//        if ((*at)->name == "CA_" || (*at)->name == "CB_") cerr << *at << endl;
//    } cerr << endl;


    vector<stl_ptr<ATOM> > atoms2check;
    vector<stl_ptr<ATOM> > sp2check;
    tr1::unordered_set<int> do_not_check;
    for (atoms_vec at=atoms.begin(); at!=atoms.end(); ++at) {
        if ((*at)->element == "H") continue;
        if ((*at)->ext->is_planar_ring && (!(*at)->ext->not_full_sp2_ring)) {
            (*at)->ext->sec_check = false;
            if ((*at)->element == "C") atoms2check.push_back(*at);
            else if ((*at)->element == "N") {
                if ((*at)->bonded_atoms.size() == 2) atoms2check.push_back(*at);
                else if (allow_protonated) {
                    atoms2check.push_back(*at);
                    (*at)->ext->sec_check = true;
                }
            }
        } else if ((*at)->element == "C") {
            //!***************** Kohlenstoff ******************
            if ((*at)->ext->hybridization == 1) { //linearer Kohlenstoff
                if ((*at)->bonded_atoms.size() > 2) { //kann nicht sein => auf sp2 pruefen
                    (*at)->ext->hybridization = 2; // durchlaeuft noch die 2er Pruefung
                } else if ((*at)->bonded_atoms.size() == 2) {
                    if (is_linear(*at)) {
                        if ((*at)->bonded_atoms[0]->ext->hybridization == 1) {
                            if ((*at)->bonded_atoms[0]->bonded_atoms.size() == 2) {
                                if (is_linear((*at)->bonded_atoms[0])) {
                                    sp2check.push_back(*at);
                                    continue;
                                }
                            } else if ((*at)->bonded_atoms[0]->bonded_atoms.size() == 1) {
                                if ((*at)->bonded_atoms[0]->element == "N") {
                                    if (get_bonds) add_bond_object((*at)->bonded_atoms[0],*at,string("3"),true);
                                    continue;
                                }
                            }
                        }
                        if ((*at)->bonded_atoms[1]->ext->hybridization == 1) {
                            if ((*at)->bonded_atoms[1]->bonded_atoms.size() == 2) {
                                if (is_linear((*at)->bonded_atoms[1])) {
                                    sp2check.push_back(*at);
                                    continue;
                                }
                            } else if ((*at)->bonded_atoms[1]->bonded_atoms.size() == 1) {
                                if ((*at)->bonded_atoms[1]->element == "N") {
                                    if (get_bonds) add_bond_object((*at)->bonded_atoms[1],*at,string("3"),true);
                                    continue;
                                }
                            }
                        }
                        if ((*at)->bonded_atoms[0]->ext->hybridization < 3 && (*at)->bonded_atoms[0]->bonded_atoms.size() < 4 &&
                            (*at)->bonded_atoms[1]->ext->hybridization < 3 && (*at)->bonded_atoms[1]->bonded_atoms.size() < 4){ // ==C==
                            do_not_check.insert((*at)->bonded_atoms[0]->intern_id);
                            do_not_check.insert((*at)->bonded_atoms[1]->intern_id);
                            if (get_bonds) {
                                add_bond_object((*at)->bonded_atoms[0],*at,string("2"),true);
                                add_bond_object((*at)->bonded_atoms[1],*at,string("2"),true);
                            }
                            continue;
                        }
                        sp2check.push_back(*at);
                        continue;
                    } else { //kann nicht sein => auf sp2 pruefen
                        (*at)->ext->hybridization = 2; // durchlaeuft noch die 2er Pruefung
                    }
                } else {
                    if ((*at)->bonded_atoms.size() == 1) {
                        if ((*at)->bonded_atoms[0]->ext->hybridization == 1) {
                            if ((*at)->bonded_atoms[0]->bonded_atoms.size() == 2) {
                                if (is_linear((*at)->bonded_atoms[0])) { // endstaendiger sp, der an anderes lineares sp-Atom gebunden ist
                                    sp2check.push_back(*at);             // Der einzige Fall, wo eine endst. Dreifachbindung erlaubt ist
                                    continue;
                                }
                            }
                        }
                        (*at)->ext->hybridization = 2; // durchlaeuft noch die 2er Pruefung
                    }
                    // else Methan
                }
            }
            if ((*at)->ext->hybridization == 2) { //trigonal planarer Kohlenstoff
                if ((*at)->bonded_atoms.size() > 3) {
                    (*at)->ext->hybridization = 3;
                    continue;
                } else if ((*at)->bonded_atoms.size() == 3) {
                    if (check_trigo(*at)) {
                        int N_count = 0;
                        for (atoms_vec bt=(*at)->bonded_atoms.begin(); bt!=(*at)->bonded_atoms.end(); ++bt) {
                            if ((*bt)->element == "N") ++N_count;
                        }
                        if (N_count > 1) {
                            for (atoms_vec bt=(*at)->bonded_atoms.begin(); bt!=(*at)->bonded_atoms.end(); ++bt) {
                                if ((*bt)->element == "N") { // OK, weil hier nur NICHT-Ring-Ns betroffen sind
                                    do_not_check.insert((*bt)->intern_id);
                                }
                            }
                        } else atoms2check.push_back(*at);
                        continue; //ok - ist wirklich sp2
                    } else {
                        (*at)->ext->hybridization = 3; //nicht planar => auf sp3 setzen
                        continue;
                    }
                } else if ((*at)->ext->n_heavy_bonded == 2) {
                    //Torsion um vermeintliche DB pruefen
                    atoms2check.push_back(*at);
                    continue;
                } else {
                    atoms2check.push_back(*at);
                    continue;
                }
            }


        } else if ((*at)->element == "N") {
            //!***************** Stickstoff ******************
            if ((*at)->bonded_atoms.size() == 4) {
                (*at)->ext->hybridization = 3; //positiv geladener Stickstoff
                continue;
            }
            if ((*at)->ext->n_heavy_bonded == 0) {(*at)->ext->hybridization = 3; continue;} //Ammoniak
            //!auf  sp  pruefen:
            //wird derzeit nur in Cyano Gruppen zugelassen
            if ((*at)->ext->hybridization == 1) {
                if ((*at)->bonded_atoms.size() > 1) { //geht in meiner Definition nicht
                    (*at)->ext->hybridization = 2; //auf sp2 pruefen
                } else if ((*at)->bonded_atoms.size() == 1) {
                    if ((*at)->bonded_atoms[0]->element == "C" && (*at)->bonded_atoms[0]->bonded_atoms.size() == 2) {
                        if (is_linear((*at)->bonded_atoms[0])) {
                            do_not_check.insert((*at)->bonded_atoms[0]->intern_id);
                            if (get_bonds) add_bond_object((*at)->bonded_atoms[0],*at,string("3"),true);
                            continue; //ist wirklich sp , weil an sp Kohlenstoff gebunden
                        }
                    }
                    (*at)->ext->hybridization = 2; //ist keine Cyano-Gruppe => auf sp2 pruefen
                }
            }
            //!auf sp2 pruefen:
            if ((*at)->ext->hybridization == 2) { //trigonal planarer Stickstoff
                if ((*at)->bonded_atoms.size() == 3) {
                    if (check_trigo(*at)) {
                        atoms2check.push_back(*at);
                        continue; //ok - ist wirklich sp2
                    } else {
                        (*at)->ext->hybridization = 3; //nicht planar => sp3
                        continue;
                    }
                } else {
                    if ((*at)->ext->n_heavy_bonded == 1) { //endstaendiger Stickstoff
                        if ((*at)->bonded_atoms.size() > 2) {
                            (*at)->ext->hybridization = 3;
                            continue;
                        } else {
                            atoms2check.push_back(*at);
                            continue;
                        }
                    } else if ((*at)->ext->n_heavy_bonded == 2) {
                        atoms2check.push_back(*at);
                        continue;
                    } else { // sollte nicht vorkommen:
                        atoms2check.push_back(*at);
                        continue;
                    }
                }
            }


        } else if ((*at)->element == "O") {
            //!***************** Sauerstoff ******************
            //!auf sp  pruefen:
            if ((*at)->ext->hybridization == 1) (*at)->ext->hybridization = 2; //gibt kein sp Sauerstoff
            //!auf sp2  pruefen:
            if ((*at)->ext->hybridization == 2) {
                if ((*at)->ext->n_heavy_bonded == 0) {(*at)->ext->hybridization = 3; (*at)->intern_type = "O.h2o"; continue;} // Wasser
                else if ((*at)->ext->n_heavy_bonded == 1) {
                    if ((*at)->bonded_atoms.size() == 1) {
                        if ((*at)->bonded_atoms[0]->ext->n_heavy_bonded == 3) {
                            if (check_trigo((*at)->bonded_atoms[0])) { //wenn an planarem C => ist Carbonyl-O
                                                                       //!Enole sollen auch als Keton gesetzt werden
                                if ((*at)->bonded_atoms[0]->ext->is_planar_ring) {
                                    if ((*at)->bonded_atoms[0]->ext->not_full_sp2_ring) {
                                        do_not_check.insert((*at)->bonded_atoms[0]->intern_id);
                                        continue;
                                    } else {
                                        (*at)->bonded_atoms[0]->ext->check_hyb = true;
                                        atoms2check.push_back(*at);
                                        continue;
                                    }
                                } else {
                                    do_not_check.insert((*at)->bonded_atoms[0]->intern_id);
                                    continue;
                                }
                            }
                            (*at)->ext->hybridization = 3; //auf sp3 setzen, weil nicht an planarem Atom
                            continue;
                        }
                        atoms2check.push_back(*at);
                        continue;
                    } else {
                        (*at)->ext->hybridization = 3; //ist wohl ne OH-Gruppe => sp3 setzen
                        continue;
                    }
                } else { // 2 heavy bonded
                    (*at)->ext->hybridization = 3;
                    continue;
                }
            }
            //!auf sp3  pruefen:
            if ((*at)->ext->hybridization == 3) {
                if ((*at)->ext->n_heavy_bonded == 0) {(*at)->intern_type = "O.h2o"; continue;} // Wasser
                if ((*at)->bonded_atoms.size() == 1) {
                    if ((*at)->bonded_atoms[0]->ext->n_heavy_bonded == 3) {
                        if (check_trigo((*at)->bonded_atoms[0])) { //!Enole sollen auch als Keton gesetzt werden
                            if ((*at)->bonded_atoms[0]->ext->is_planar_ring) {
                                if ((*at)->bonded_atoms[0]->ext->not_full_sp2_ring) {
                                    do_not_check.insert((*at)->bonded_atoms[0]->intern_id);
                                    continue;
                                }
                            } else {
                                do_not_check.insert((*at)->bonded_atoms[0]->intern_id);
                                continue;
                            }
                        } else continue;
                    } else if ((*at)->bonded_atoms[0]->ext->n_heavy_bonded == 1) { // O=O
                        (*at)->ext->hybridization = 2;
                        continue;
                    } else continue;
                }
            }


        } else if ((*at)->element == "S") {
            //!***************** Schwefel ******************
            if ((*at)->ext->hybridization == 2) {
                bool NO_flag = false;
                for (atoms_vec bt=(*at)->bonded_atoms.begin(); bt!=(*at)->bonded_atoms.end(); ++bt) {
                    if ((*bt)->element == "O" || (*bt)->element == "N") {
                        NO_flag = true;
                        break;
                    }
                }
                if (NO_flag) continue;
                else atoms2check.push_back(*at);
            }
        }
    }


    // Dreifachbindungen:
    if (sp2check.size() > 0 && get_bonds) {
        if (sp2check.size() == 1) {
            sp2check[0]->ext->hybridization = 2;
            atoms2check.push_back(sp2check[0]);
        } else {
            vector<MW_EDGE*> edges;
            tr1::unordered_map<ATOM*,unsigned int> at2id;
            for (unsigned int i=0; i<sp2check.size(); ++i) at2id[sp2check[i].get_pointer()] = i;
            for (unsigned int i=0; i<sp2check.size(); ++i) {
                if (sp2check[i]->bonded_atoms.size() > 2) continue;
                for (atoms_vec bt=sp2check[i]->bonded_atoms.begin(); bt!=sp2check[i]->bonded_atoms.end(); ++bt) {
                    if (at2id.find(bt->get_pointer()) == at2id.end()) continue;
                    if (at2id[bt->get_pointer()] < i) continue;
                    if ((*bt)->bonded_atoms.size() > 2) continue;
                    float dst = get_distance(sp2check[i]->coord,(*bt)->coord);
                    float weight = EDnew::EDAccess::get_probability(sp2check[i]->element,(*bt)->element,3,dst);

                    if (weight < 0.03) continue;

//                    cerr << "push " << atoms2check[i]->name << "," << (*bt)->name << ",w=" << weight
//                         << "   (" << i << "," << at2id[bt->get_pointer()] << ")" << endl;

                    edges.push_back(new MW_EDGE(i,at2id[bt->get_pointer()],weight));
                }
            }

            MW_MATCH matcher;
            matcher.solve(sp2check.size(),edges);

            for (vector<MW_EDGE*>::const_iterator it=matcher.get_max_match().begin(); it!=matcher.get_max_match().end(); ++it) {
                add_bond_object(sp2check[(*it)->node_i],sp2check[(*it)->node_j],string("3"),true);
            }
            for (vector<MW_EDGE*>::const_iterator it=matcher.get_non_max().begin(); it!=matcher.get_non_max().end(); ++it) {
                add_bond_object(sp2check[(*it)->node_i],sp2check[(*it)->node_j],string("1"),true);
            }
            for (vector<MW_EDGE*>::iterator it=edges.begin(); it!=edges.end(); ++it) delete *it;
        }
    }


    // Jetzt max. weighted match fuer die atoms2check:
    if (atoms2check.size() > 0) {
        for (unsigned int i=0; i<atoms2check.size(); ++i) {
            if (atoms2check[i]->element == "N") {
                for (tr1::unordered_set<RING*>::iterator rit=atoms2check[i]->ext->ring_ptrs.begin();
                                                         rit!=atoms2check[i]->ext->ring_ptrs.end(); ++rit) {
                    //if ((*rit)->n_members % 2 == 0) { // In Ringen mit gerader Anzahl DBs zum N bevorzugen
                    if ((*rit)->n_members == 6) {
                        atoms2check[i]->ext->check_hyb = true;
                        break;
                    }
                }
                for (atoms_vec bt=atoms2check[i]->bonded_atoms.begin(); bt!=atoms2check[i]->bonded_atoms.end(); ++bt) {
                    if ((*bt)->element == "C" && (*bt)->ext->check_hyb) {
                        (*bt)->ext->sec_check = true; //Amidcarbonyl kennzeichnen
                    }
                }
            }
        }

        vector<MW_EDGE*> edges;
        tr1::unordered_map<ATOM*,unsigned int> at2id;
        for (unsigned int i=0; i<atoms2check.size(); ++i) at2id[atoms2check[i].get_pointer()] = i;
        for (unsigned int i=0; i<atoms2check.size(); ++i) {
            if (do_not_check.find(atoms2check[i]->intern_id) != do_not_check.end()) continue;
            for (atoms_vec bt=atoms2check[i]->bonded_atoms.begin(); bt!=atoms2check[i]->bonded_atoms.end(); ++bt) {
                if (at2id.find(bt->get_pointer()) == at2id.end()) continue;
                if (at2id[bt->get_pointer()] < i) continue;
                if (do_not_check.find((*bt)->intern_id) != do_not_check.end()) continue;

                float weight = 0.;
                if (atoms2check[i]->ext->is_planar_ring && (!atoms2check[i]->ext->not_full_sp2_ring) &&
                    (*bt)->ext->is_planar_ring && (!(*bt)->ext->not_full_sp2_ring)) {
                    float dst = get_distance(atoms2check[i]->coord,(*bt)->coord);
                    float w1 = EDnew::EDAccess::get_probability(atoms2check[i]->element,(*bt)->element,2,dst);
                    float w2 = EDnew::EDAccess::get_probability(atoms2check[i]->element,(*bt)->element,0,dst);
                    if (w2 > w1) weight = w2;
                    else weight = w1;
                    if (!atoms2check[i]->is_in_same_ring(*bt)) {
                        if (!check_planar_tors(atoms2check[i],*bt)) continue;
                        else weight += 1.5;
                    } else {
//                        if (atoms2check[i]->element == "C") {
//                            if ((*bt)->element == "C") {
//                                if (atoms2check[i]->ext->check_hyb || (*bt)->ext->check_hyb) weight += 4.5;//4.0
//                                else weight += 6.5;
//                                if (atoms2check[i]->ext->sec_check || (*bt)->ext->sec_check) weight -= 2.5;
//                            } else if ((*bt)->element == "N") {
//                                if (atoms2check[i]->ext->check_hyb) continue; // Amid
//                                if ((*bt)->ext->sec_check) weight += 1.5;
//                                else weight += 2.5;
//                                if ((*bt)->ext->check_hyb) weight += 1.5; // gerade Anzahl Ringatome
//                            }
//                        } else if (atoms2check[i]->element == "N") {
//                            if ((*bt)->element == "C") {
//                                if ((*bt)->ext->check_hyb) continue; // Amid
//                                if (atoms2check[i]->ext->sec_check) weight += 1.5;
//                                else weight += 2.5;
//                                if (atoms2check[i]->ext->check_hyb) weight += 1.5; // gerade Anzahl Ringatome
//                            } else if ((*bt)->element == "N") {
//                                if ((*bt)->ext->sec_check || atoms2check[i]->ext->sec_check) {
//                                    if ((*bt)->ext->sec_check && atoms2check[i]->ext->sec_check) continue;
//                                    else weight += 1.3;
//                                } else weight += 2.3;
//                                if (atoms2check[i]->ext->check_hyb && (*bt)->ext->check_hyb) weight += 1.5; // gerade Anzahl Ringatome
//                            }
//                        }



                        if (atoms2check[i]->element == "C") {
                            if (atoms2check[i]->ext->check_hyb) weight += 3.;
                            else weight += 4.;
                            if (atoms2check[i]->ext->sec_check) weight -= 2.;
                            else if ((*bt)->element == "C") weight += 2.;
                        }
                        else if (atoms2check[i]->element == "N") {
                            if (atoms2check[i]->ext->sec_check) weight += 1.;
                            else weight += 2.;
                            if ((*bt)->ext->check_hyb) weight -= 2.1;
                        }
                        if ((*bt)->element == "C") {
                            if ((*bt)->ext->check_hyb) weight += 3.;
                            else weight += 4.;
                            if ((*bt)->ext->sec_check) weight -= 2.;
                            else if (atoms2check[i]->element == "C") weight += 2.;
                        }
                        else if ((*bt)->element == "N") {
                            if ((*bt)->ext->sec_check) weight += 1.;
                            else weight += 2.;
                            if (atoms2check[i]->ext->check_hyb) weight -= 2.1;
                        }

                    }

                    if (weight < 0.03) continue;

//                    cerr << "push " << atoms2check[i]->name << "," << (*bt)->name << ",w=" << weight
//                         << "   (" << i << "," << at2id[bt->get_pointer()] << ")" << endl;

                    edges.push_back(new MW_EDGE(i,at2id[bt->get_pointer()],weight));
                } else {
                    if (!check_planar_tors(atoms2check[i],*bt)) continue;
                    float dst = get_distance(atoms2check[i]->coord,(*bt)->coord);
                    float w1 = EDnew::EDAccess::get_probability(atoms2check[i]->element,(*bt)->element,2,dst);
                    float w2 = EDnew::EDAccess::get_probability(atoms2check[i]->element,(*bt)->element,0,dst);
                    if (w2 > w1) weight = w2;
                    else weight = w1;

                    float aw1 = 0.;
                    if (atoms2check[i]->bonded_atoms.size() > 1) aw1 = EDnew::EDAccess::get_angle_probability(atoms2check[i]->element,1,atoms2check[i]->get_smallest_bond_angle());
                    if (aw1 < 0.) aw1 = 0.;
                    float aw2 = 0.;
                    if ((*bt)->bonded_atoms.size() > 1) aw2 = EDnew::EDAccess::get_angle_probability((*bt)->element,1,(*bt)->get_smallest_bond_angle());
                    if (aw2 < 0.) aw2 = 0.;
                    if (aw1 < aw2) weight += aw1;
                    else weight += aw2;

                    if (atoms2check[i]->ext->is_planar_ring && (!atoms2check[i]->ext->not_full_sp2_ring)) {
                        if ((*bt)->element == "O") weight += 1.;
                        else if ((*bt)->element == "C" || (*bt)->element == "N") {
                            if ((*bt)->ext->n_heavy_bonded > 1) {
                                bool increase = true;
                                for (unsigned int c=0; c<(*bt)->bonded_atoms.size(); ++c) {
                                    if (!((*bt)->bonded_atoms[c]->ext->is_planar_ring) ||
                                        (*bt)->bonded_atoms[c]->ext->not_full_sp2_ring) {
                                        increase = false;
                                        break;
                                    }
                                }
                                if (increase) {
                                    if ((*bt)->element == "C" && (*bt)->bonded_atoms.size() < 4) weight += 12.;//6.0;
                                    else if ((*bt)->element == "N" && (*bt)->bonded_atoms.size() < 3) weight += 6.;//4.0;
                                }
                            }
                        }
                    } else if ((*bt)->ext->is_planar_ring && (!(*bt)->ext->not_full_sp2_ring)) {
                        if (atoms2check[i]->element == "O") weight += 1.;
                        else if (atoms2check[i]->element == "C" || atoms2check[i]->element == "N") {
                            if (atoms2check[i]->ext->n_heavy_bonded > 1) {
                                bool increase = true;
                                for (unsigned int c=0; c<atoms2check[i]->bonded_atoms.size(); ++c) {
                                    if (!(atoms2check[i]->bonded_atoms[c]->ext->is_planar_ring) ||
                                        atoms2check[i]->bonded_atoms[c]->ext->not_full_sp2_ring) {
                                        increase = false;
                                        break;
                                    }
                                }
                                if (increase) {
                                    if (atoms2check[i]->element == "C" && atoms2check[i]->bonded_atoms.size() < 4) weight += 12.;//6.0;
                                    else if (atoms2check[i]->element == "N" && atoms2check[i]->bonded_atoms.size() < 3) weight += 6.;//4.0;
                                }
                            }
                        }
                    } else if (atoms2check[i]->ext->n_heavy_bonded < 2 || (*bt)->ext->n_heavy_bonded < 2) weight *= 0.7; // endstaendige runter gewichten

                    if ((atoms2check[i]->element == "N" && (*bt)->element == "O") ||
                        (atoms2check[i]->element == "O" && (*bt)->element == "N")) weight *= 0.35;
                    else if (atoms2check[i]->element == "O" && atoms2check[i]->bonded_atoms[0]->ext->n_heavy_bonded == 2) weight *= 0.7;
                    else if ((*bt)->element == "O" && (*bt)->bonded_atoms[0]->ext->n_heavy_bonded == 2) weight *= 0.7;

                    if (weight < 0.03) continue;

//                    cerr << "push " << atoms2check[i]->name << "," << (*bt)->name << ",w=" << weight << "   aw1=" << aw1
//                         << "   aw2=" << aw2 << "   (" << i << "," << at2id[bt->get_pointer()] << ")" << endl;

                    edges.push_back(new MW_EDGE(i,at2id[bt->get_pointer()],weight));
                }
            }
        }

        MW_MATCH matcher;
        matcher.solve(atoms2check.size(),edges);


//        cerr << "Match:" << endl;
//        for (vector<MW_EDGE*>::const_iterator it=matcher.get_max_match().begin(); it!=matcher.get_max_match().end(); ++it) {
//            cerr << "(" << atoms2check[(*it)->node_i]->name << "," << atoms2check[(*it)->node_j]->name << ") ";
//        }
//        cerr << endl;


        for (unsigned int i=0; i<atoms2check.size(); ++i) {
            if (do_not_check.find(atoms2check[i]->intern_id) != do_not_check.end()) continue;
            if (matcher.is_in_max(i) == false) {
                atoms2check[i]->ext->hybridization = 3;
                if (atoms2check[i]->ext->is_planar_ring && !atoms2check[i]->ext->not_full_sp2_ring) {
                    if (atoms2check[i]->element == "N") {
                        atoms2check[i]->ext->hybridization = 2;
                    } else for (tr1::unordered_set<RING*>::iterator it=atoms2check[i]->ext->ring_ptrs.begin();
                                                             it!=atoms2check[i]->ext->ring_ptrs.end(); ++it) {
                        for (atoms_vec rat=(*it)->ring_atoms.begin(); rat!=(*it)->ring_atoms.end(); ++rat) {
                            (*rat)->ext->not_full_sp2_ring = true;
                        }
                        (*it)->is_full_sp2 = false;
                    }
                } else if (atoms2check[i]->element == "N" && atoms2check[i]->bonded_atoms.size() == 3) {
                    atoms2check[i]->ext->hybridization = 2;
                }
            }
        }


        //#######################################################################
        // Neue Aro-Zuweisung:
        int* em = new int[atoms2check.size()+1];
        for (int i=0; i<(atoms2check.size()+1); ++i) em[i] = -1;
        for (vector<MW_EDGE*>::const_iterator it=matcher.get_max_match().begin(); it!=matcher.get_max_match().end(); ++it) {
            em[(*it)->node_i] = (*it)->node_j;
            em[(*it)->node_j] = (*it)->node_i;
        }
        tr1::unordered_map<RING*,int> pi_count;
        tr1::unordered_map<RING*,tr1::unordered_set<RING*> > exo_DBs; // zeigt auf Nachbarringe in denen exo-DBs liegen
        tr1::unordered_map<RING*,int> aro_status; // -1: unprozessiert / 0: unmoeglich ein Aromat / 1: candidate / 2: fused_candidate
        for (rings_vec rit=rings.begin(); rit!=rings.end(); ++rit) {
            (*rit)->is_aromatic = false;
            aro_status[rit->get_pointer()] = -1;
            exo_DBs[rit->get_pointer()] = tr1::unordered_set<RING*>();
            pi_count[rit->get_pointer()] = 0; // Pi-Elektronen fuer den Ring zaehlen
            if ((*rit)->is_full_sp2) { // nicht planar und komplett sp2-Hyb
                if ((*rit)->atom_set.size() == 0) { // Atom-Set fuer Pruefungen anlegen
                    for (atoms_vec rat=(*rit)->ring_atoms.begin(); rat!=(*rit)->ring_atoms.end(); ++rat) {
                        (*rit)->atom_set.insert((*rat)->intern_id);
                    }
                }
            } else aro_status[rit->get_pointer()] = 0;
        }
        for (rings_vec rit=rings.begin(); rit!=rings.end(); ++rit) {
            if (aro_status[rit->get_pointer()] == 0) continue;
            for (atoms_vec rat=(*rit)->ring_atoms.begin(); rat!=(*rit)->ring_atoms.end(); ++rat) {
                if (at2id.find(rat->get_pointer()) != at2id.end()) {
                    if (matcher.is_in_max(at2id[rat->get_pointer()])) {
                        if ((*rit)->atom_set.find(atoms2check[em[at2id[rat->get_pointer()]]]->intern_id) != (*rit)->atom_set.end()) {
                            pi_count[rit->get_pointer()] += 1; // DB ist in (*rit)
                            continue;
                        } else { // DB ist NICHT in (*rit)
                            aro_status[rit->get_pointer()] = 0;
                            for (rings_vec rit2=rings.begin(); rit2!=rings.end(); ++rit2) {
                                if (aro_status[rit2->get_pointer()] == 0) continue;
                                if (rit->get_pointer() == rit2->get_pointer()) continue;
                                if ((*rit2)->atom_set.find(atoms2check[em[at2id[rat->get_pointer()]]]->intern_id) != (*rit2)->atom_set.end()) {
                                    if ((*rit)->fused(*rit2)) {
                                        exo_DBs[rit->get_pointer()].insert(rit2->get_pointer());
                                        aro_status[rit->get_pointer()] = 2;
                                        break;
                                    }
                                }
                            }
                            if (aro_status[rit->get_pointer()] == 0) break;
                        }
                    } else if ((*rat)->element == "N" || (*rat)->element == "O" || (*rat)->element == "S") pi_count[rit->get_pointer()] += 2;
                    else {aro_status[rit->get_pointer()] = 0; break;}
                } else {
                    if ((*rat)->element == "N" || (*rat)->element == "O" || (*rat)->element == "S") pi_count[rit->get_pointer()] += 2;
                    else {aro_status[rit->get_pointer()] = 0; break;}
                }
            }
            if (aro_status[rit->get_pointer()] == -1) aro_status[rit->get_pointer()] = 1;
        }
        delete[] em;

        // Kandidaten muessen 4n+2 Pi-Elektronen haben:
        for (rings_vec rit=rings.begin(); rit!=rings.end(); ++rit) {
            if (aro_status[rit->get_pointer()] != 1) continue;
            if (!(pi_count[rit->get_pointer()] == 6 || pi_count[rit->get_pointer()] == 10 || pi_count[rit->get_pointer()] == 14 ||
                  pi_count[rit->get_pointer()] == 18 || pi_count[rit->get_pointer()] == 22 || pi_count[rit->get_pointer()] == 26)) {
                aro_status[rit->get_pointer()] = 0;
            } else (*rit)->is_aromatic = true;
        }

        // Fused_Candidates (Kandidaten mit exocyklischer DB) als ganzes System pruefen:
        for (rings_vec rit=rings.begin(); rit!=rings.end(); ++rit) {
            if ((*rit)->is_aromatic || aro_status[rit->get_pointer()] != 2) continue;
            tr1::unordered_set<RING*> fused_system;
            stack<RING*> fuse_stack;
            fuse_stack.push(rit->get_pointer());
            fused_system.insert(rit->get_pointer());
            while (!fuse_stack.empty()) {
                RING* cr = fuse_stack.top(); fuse_stack.pop();
                for (tr1::unordered_set<RING*>::iterator rit2=exo_DBs[cr].begin(); rit2!=exo_DBs[cr].end(); ++rit2) {
                    if (fused_system.find(*rit2) != fused_system.end()) continue;
                    fused_system.insert(*rit2);
                    fuse_stack.push(*rit2);
                }
            }
            int total_pi = 0;
            for (tr1::unordered_set<RING*>::iterator rjt=fused_system.begin(); rjt!=fused_system.end(); ++rjt) {
                total_pi += pi_count[*rjt];
            }
            if (!(total_pi == 6 || total_pi == 10 || total_pi == 14 ||
                  total_pi == 18 || total_pi == 22 || total_pi == 26 ||
                  total_pi == 30 || total_pi == 34 || total_pi == 38)) {
                aro_status[rit->get_pointer()] = 0;
            } else { // Alles aromatisch
                for (tr1::unordered_set<RING*>::iterator rjt=fused_system.begin(); rjt!=fused_system.end(); ++rjt) {
                    (*rjt)->is_aromatic = true;
                }
            }
        }

        for (rings_vec rit=rings.begin(); rit!=rings.end(); ++rit) {
            if ((*rit)->is_aromatic) {
                for (atoms_vec rat=(*rit)->ring_atoms.begin(); rat!=(*rit)->ring_atoms.end(); ++rat) {
                    (*rat)->ext->is_aromatic = true;
                    if ((*rat)->element == "N") {
                        if (matcher.is_in_max(at2id[rat->get_pointer()])) {
                            if ((*rat)->bonded_atoms.size() > 2) {
                                if ((*rit)->n_members == 6) {
                                    (*rat)->intern_type = "N.ar6p";
                                    (*rit)->is_pos = true;
                                    if ((*rat)->ext->n_heavy_bonded == 2) (*rit)->is_prot = true;
                                } else {
                                    (*rat)->intern_type = "N.arp";
                                    (*rit)->is_pos = true;
                                    if ((*rat)->ext->n_heavy_bonded == 2) (*rit)->is_prot = true;
                                }
                            } else {
                                if ((*rit)->n_members == 6) (*rat)->intern_type = "N.ar6";
                                else {
                                    (*rat)->intern_type = "N.ar2";
                                    (*rat)->sybyl_type = "X"; //!spaeter nicht auf N.ar3h setzen
                                }
                            }
                        } else {
                            if ((*rat)->ext->n_heavy_bonded == 2) (*rat)->intern_type = "N.ar3h";
                            else (*rat)->intern_type = "N.ar3";
                        }
                    }
                }
            }
        }
        // Ende neue Aro-Zuweisung
        //#######################################################################


        for (vector<MW_EDGE*>::const_iterator it=matcher.get_max_match().begin(); it!=matcher.get_max_match().end(); ++it) {
            bool conti = true;
            if ((*it)->weight < 0.06) conti = false;
            else {
                float w1 = EDnew::EDAccess::get_probability(atoms2check[(*it)->node_i]->element,
                                                            atoms2check[(*it)->node_j]->element,2,
                                                            get_distance(atoms2check[(*it)->node_i]->coord,atoms2check[(*it)->node_j]->coord));
                if (w1 < 0.06) {
                    if (atoms2check[(*it)->node_i]->element == "C" && atoms2check[(*it)->node_j]->element == "C") {
                        float w2 = EDnew::EDAccess::get_probability(atoms2check[(*it)->node_i]->element,
                                                                    atoms2check[(*it)->node_j]->element,0,
                                                                    get_distance(atoms2check[(*it)->node_i]->coord,atoms2check[(*it)->node_j]->coord));
                        if (w2 < 0.05) conti = false;
                    } else conti = false;
                }
            }
            if (conti) continue;

            bool is_isolated = true;
            for (atoms_vec bt=atoms2check[(*it)->node_i]->bonded_atoms.begin(); bt!=atoms2check[(*it)->node_i]->bonded_atoms.end(); ++bt) {
                if ((*bt)->intern_id == atoms2check[(*it)->node_j]->intern_id) continue;
                if (at2id.find(bt->get_pointer()) == at2id.end() || do_not_check.find((*bt)->intern_id) != do_not_check.end()) {
                    if ((*bt)->ext->hybridization == 2) {
                        is_isolated = false;
                        break;
                    } else continue;
                } else if (matcher.is_in_max(at2id[bt->get_pointer()])) {
                    is_isolated = false;
                    break;
                }
            }
            if (!is_isolated) continue;
            for (atoms_vec bt=atoms2check[(*it)->node_j]->bonded_atoms.begin(); bt!=atoms2check[(*it)->node_j]->bonded_atoms.end(); ++bt) {
                if ((*bt)->intern_id == atoms2check[(*it)->node_i]->intern_id) continue;
                if (at2id.find(bt->get_pointer()) == at2id.end() || do_not_check.find((*bt)->intern_id) != do_not_check.end()) {
                    if ((*bt)->ext->hybridization == 2) {
                        is_isolated = false;
                        break;
                    } else continue;
                } else if (matcher.is_in_max(at2id[bt->get_pointer()])) {
                    is_isolated = false;
                    break;
                }
            }
            if (is_isolated) {
                //cerr << "[" << atoms2check[(*it)->node_i]->name << "," << atoms2check[(*it)->node_j]->name << "] is isolated" << endl;
                atoms2check[(*it)->node_i]->ext->hybridization = 3;
                atoms2check[(*it)->node_j]->ext->hybridization = 3;
            }
        }

        if (get_bonds) {
            // Die bestimmten Single und DB soweit schonmal setzen (noch keine 'ar' und keine zu 'ar'):
            for (vector<MW_EDGE*>::const_iterator it=matcher.get_max_match().begin(); it!=matcher.get_max_match().end(); ++it) {
                if (atoms2check[(*it)->node_i]->ext->is_aromatic && !kekulize_aromatics) continue;
                if (atoms2check[(*it)->node_i]->ext->hybridization == 3 || atoms2check[(*it)->node_j]->ext->hybridization == 3) { //  Isoliert und auf sp3 gesetzt
                    add_bond_object(atoms2check[(*it)->node_i],atoms2check[(*it)->node_j],string("1"),false);
                } else add_bond_object(atoms2check[(*it)->node_i],atoms2check[(*it)->node_j],string("2"),false);
            }
            for (vector<MW_EDGE*>::const_iterator it=matcher.get_non_max().begin(); it!=matcher.get_non_max().end(); ++it) {
                if (!kekulize_aromatics && (atoms2check[(*it)->node_i]->ext->is_aromatic ||
                                            atoms2check[(*it)->node_j]->ext->is_aromatic)) continue;
                add_bond_object(atoms2check[(*it)->node_i],atoms2check[(*it)->node_j],string("1"),false);
            }
        }

        for (vector<MW_EDGE*>::iterator it=edges.begin(); it!=edges.end(); ++it) delete *it;
    }
}


void MOLECULE::get_rings(int const& max_members) {
    for (atoms_vec ct=atoms.begin(); ct!=atoms.end(); ++ct) (*ct)->ext->prev = 0;
    for (atoms_vec at=atoms.begin(); at!=atoms.end(); ++at) {

        if ((*at)->element == "H" || (*at)->ext->n_heavy_bonded < 2) continue;

        if ((*at)->ext->is_ring) {
            if ((*at)->ext->n_heavy_bonded < 3) continue;
            else if ((*at)->ext->n_heavy_bonded < 4 && (*at)->ext->ring_ptrs.size() > 2) continue;
            else if ((*at)->ext->ring_ptrs.size() > 5) continue;
        }

        (*at)->ext->prev = 0;

        (*at)->ext->level = 1;
        list<stl_ptr<ATOM> > pq;
        pq.push_back(*at);

        bool next_root = false;
        int shortest = max_members;

        vector<stl_ptr<ATOM> > reset_prevs;
        reset_prevs.push_back(*at);

        while (!pq.empty()) {
            if (next_root) break;
            stl_ptr<ATOM> v = pq.front(); pq.pop_front();

            for (atoms_vec bt=v->bonded_atoms.begin(); bt!=v->bonded_atoms.end(); ++bt) {

                if (!v->ext->prev.zero()) {
                    if ((*bt)->intern_id == v->ext->prev->intern_id) continue;
                }

                if ((*bt)->element == "H" || (*bt)->ext->n_heavy_bonded < 2) continue;

                if (!(*bt)->ext->prev.zero()) { // schon aus anderer Richtung bekannt

                    if (v->ext->level + (*bt)->ext->level -1 > shortest) break;

                    // Pruefen, ob ein gemeinsamer nicht-root Vorgaenger vorhanden ist:
                    bool conti = false;
                    stl_ptr<ATOM> tat = v->ext->prev;
                    while (!tat->ext->prev.zero()) {
                        stl_ptr<ATOM> tbt = *bt;
                        do {
                            if (tat->intern_id == tbt->intern_id) {
                                conti = true;
                                break;
                            }
                            tbt = tbt->ext->prev;
                        } while (!tbt->ext->prev.zero());
                        if (conti) break;
                        tat = tat->ext->prev;
                    }
                    if (conti) continue;

                    // potentiell neuer kleinster Ring:
                    stl_ptr<RING> rng = new RING();
                    rng->n_members = v->ext->level + (*bt)->ext->level - 1;
                    rng->is_planar = false; rng->n_hetero = 0; rng->is_aromatic = false;
                    rng->id = 0; rng->sq_id = 0;

                    tat = v;
                    do {
                        rng->ring_atoms.push_back(tat);
                        rng->id += tat->intern_id;
                        rng->sq_id += (tat->intern_id * tat->intern_id);
                        tat->ext->is_ring = true;
                        if (tat->element != "C") rng->n_hetero++;
                        tat = tat->ext->prev;
                    } while (tat.zero() == false);
                    tat = *bt;
                    vector<stl_ptr<ATOM> > helpi;
                    do {
                        if (tat->intern_id == (*at)->intern_id) break;
                        helpi.push_back(tat);
                        rng->id += tat->intern_id;
                        rng->sq_id += (tat->intern_id * tat->intern_id);
                        tat->ext->is_ring = true;
                        if (tat->element != "C") rng->n_hetero++;
                        tat = tat->ext->prev;
                    } while (tat.zero() == false);

                    for (vector<stl_ptr<ATOM> >::reverse_iterator ht=helpi.rbegin(); ht!=helpi.rend(); ++ht) rng->ring_atoms.push_back(*ht);

                    bool new_ring = true;
                    for (rings_vec rit=rings.begin(); rit!=rings.end(); ++rit) {
                        if ((*rit)->n_members != rng->n_members) continue;
                        if ((*rit)->id != rng->id) continue;
                        if ((*rit)->sq_id != rng->sq_id) continue;
                        new_ring = false;
                        break;
                    }
                    if (new_ring) {
                        rng->is_full_sp2 = true;
                        for (atoms_vec rat=rng->ring_atoms.begin();
                            rat!=rng->ring_atoms.end(); ++rat) {
                            (*rat)->ext->is_ring = true;
                            (*rat)->ext->ring_ptrs.insert(rng.get_pointer());
                        }
                        rings.push_back(rng);
                        next_root = true;
                        break;
                    } else {
                        shortest = rng->n_members;
                        rng.kill();
                    }

                } else {
                    if ((2 * v->ext->level -1)  >= shortest) continue;
                    pq.push_back(*bt);
                    reset_prevs.push_back(*bt);
                    (*bt)->ext->level = v->ext->level + 1;
                    (*bt)->ext->prev = v;
                }
            }
        }

        for (atoms_vec rp=reset_prevs.begin(); rp!=reset_prevs.end(); ++rp) {
            (*rp)->ext->prev = 0;
        }

    }
}


void MOLECULE::check_ring_plan() {
    for (rings_vec rit=rings.begin(); rit!=rings.end(); ++rit) {
        if (((*rit)->n_members > 4)) {
            float testa = dihedral((*rit)->ring_atoms[0]->coord,(*rit)->ring_atoms[1]->coord,
                                   (*rit)->ring_atoms[2]->coord,(*rit)->ring_atoms[3]->coord);
            testa += dihedral((*rit)->ring_atoms[1]->coord,(*rit)->ring_atoms[2]->coord,
                             (*rit)->ring_atoms[3]->coord,(*rit)->ring_atoms[4]->coord);
            if ((*rit)->n_members == 5) {
                testa += dihedral((*rit)->ring_atoms[2]->coord,(*rit)->ring_atoms[3]->coord,
                                  (*rit)->ring_atoms[4]->coord,(*rit)->ring_atoms[0]->coord);
                testa += dihedral((*rit)->ring_atoms[3]->coord,(*rit)->ring_atoms[4]->coord,
                                  (*rit)->ring_atoms[0]->coord,(*rit)->ring_atoms[1]->coord);
                testa += dihedral((*rit)->ring_atoms[4]->coord,(*rit)->ring_atoms[0]->coord,
                                  (*rit)->ring_atoms[1]->coord,(*rit)->ring_atoms[2]->coord);
                if ((testa/5.) > 0.195) continue;
            } else if ((*rit)->n_members == 6) {
                testa += dihedral((*rit)->ring_atoms[2]->coord,(*rit)->ring_atoms[3]->coord,
                                  (*rit)->ring_atoms[4]->coord,(*rit)->ring_atoms[5]->coord);
                testa += dihedral((*rit)->ring_atoms[3]->coord,(*rit)->ring_atoms[4]->coord,
                                  (*rit)->ring_atoms[5]->coord,(*rit)->ring_atoms[0]->coord);
                testa += dihedral((*rit)->ring_atoms[4]->coord,(*rit)->ring_atoms[5]->coord,
                                  (*rit)->ring_atoms[0]->coord,(*rit)->ring_atoms[1]->coord);
                testa += dihedral((*rit)->ring_atoms[5]->coord,(*rit)->ring_atoms[0]->coord,
                                  (*rit)->ring_atoms[1]->coord,(*rit)->ring_atoms[2]->coord);
                if ((testa/6.) > 0.262) continue; //bis ca. 15 Grad Torsion zulassen
            } else if ((*rit)->n_members > 6) {
                testa = dihedral((*rit)->ring_atoms[2]->coord,(*rit)->ring_atoms[3]->coord,
                                 (*rit)->ring_atoms[4]->coord,(*rit)->ring_atoms[5]->coord);
                testa = dihedral((*rit)->ring_atoms[3]->coord,(*rit)->ring_atoms[4]->coord,
                                 (*rit)->ring_atoms[5]->coord,(*rit)->ring_atoms[6]->coord);
                if ((testa/4.) > 0.26 && (testa/4.) < 2.88) continue;
            }

            (*rit)->is_planar = true;

        } else {
            (*rit)->is_full_sp2 = false; //3 und 4 ringe zwar auch planar, aber nicht sp2
        }
        if ((*rit)->is_planar) {
            bool no_full = false;
            for (atoms_vec rat=(*rit)->ring_atoms.begin(); rat!=(*rit)->ring_atoms.end(); ++rat) {
                (*rat)->ext->is_planar_ring = true;
                if ((*rat)->bonded_atoms.size() > 3 || ((*rat)->element != "C" && (*rat)->element != "N" && // full_sp2 nur fuer CNOSP
                    (*rat)->element != "O" && (*rat)->element != "S" && (*rat)->element != "P")) {          // siehe LIYHIK
                    no_full = true;
                    break;
                } else {
                    if (!check_trigo(*rat,1.2)) { // hoher Threshold ist noetig (siehe z.B. 1TNK)
                        no_full = true;
                        break;
                    }
                }
            }
            if (no_full) {
                for (atoms_vec rat=(*rit)->ring_atoms.begin(); rat!=(*rit)->ring_atoms.end(); ++rat) {
                    (*rit)->is_full_sp2 = false;
                    (*rat)->ext->not_full_sp2_ring = true; //wird in hyb2 gebraucht
                }
            } else {
                for (atoms_vec rat=(*rit)->ring_atoms.begin(); rat!=(*rit)->ring_atoms.end(); ++rat) {
                    if ((*rat)->ext->hybridization != 2) (*rat)->ext->hybridization = 2;
                }
            }
        }
    }
}


bool MOLECULE::check_SO(stl_ptr<ATOM> const& base,int const& O_count,int const& OH_count,int const& O_2heavy,int const& N_count) {
    int diff = O_count + N_count - O_2heavy;
    bool depro = false;
    if (diff > 0) {
        if ((diff-OH_count) > 2) { // deprotoniert
            for (atoms_vec ct=base->bonded_atoms.begin(); ct!=base->bonded_atoms.end(); ++ct) {
                if ((*ct)->element == "O") {
                    if ((*ct)->bonded_atoms.size() == 2) {
                        if ((*ct)->ext->n_heavy_bonded == 1) {
                            (*ct)->intern_type = "O.3ac";
                            (*ct)->ext->hybridization = 3;
                        } else {
                            (*ct)->intern_type = "O.3so";
                            (*ct)->ext->hybridization = 3;
                        }
                    } else {
                        (*ct)->intern_type = "O.2s";
                        (*ct)->ext->hybridization = 2;
                    }
                } else if ((*ct)->element == "S") {
                    if ((*ct)->bonded_atoms.size() == 2) {
                        if ((*ct)->ext->n_heavy_bonded == 1) {
                            (*ct)->intern_type = "S.sh";
                            (*ct)->ext->hybridization = 3;
                        } else {
                            (*ct)->intern_type = "S.3";
                            (*ct)->ext->hybridization = 3;
                        }
                    } else {
                        (*ct)->intern_type = "S.2";
                        (*ct)->ext->hybridization = 2;
                    }
                } // N Typen sollten bereits korrekt gesetzt sein!
            }
            depro = true;
        } else {
            for (atoms_vec ct=base->bonded_atoms.begin(); ct!=base->bonded_atoms.end(); ++ct) {
                if ((*ct)->element == "O") {
                    if ((*ct)->bonded_atoms.size() == 2) {
                        if ((*ct)->ext->n_heavy_bonded == 1) {
                            (*ct)->intern_type = "O.3oh";
                            (*ct)->ext->hybridization = 3;
                        } else {
                            (*ct)->intern_type = "O.3so";
                            (*ct)->ext->hybridization = 3;
                        }
                    } else {
                        (*ct)->intern_type = "O.2so";
                        (*ct)->ext->hybridization = 2;
                    }
                } else if ((*ct)->element == "S") {
                    if ((*ct)->bonded_atoms.size() == 2) {
                        if ((*ct)->ext->n_heavy_bonded == 1) {
                            (*ct)->intern_type = "S.sh";
                            (*ct)->ext->hybridization = 3;
                        } else {
                            (*ct)->intern_type = "S.3";
                            (*ct)->ext->hybridization = 3;
                        }
                    } else {
                        (*ct)->intern_type = "S.2";
                        (*ct)->ext->hybridization = 2;
                    }
                }
            }
        }
    } else {
        for (atoms_vec ct=base->bonded_atoms.begin(); ct!=base->bonded_atoms.end(); ++ct) {
            if ((*ct)->element == "O") {
                if ((*ct)->bonded_atoms.size() == 2) {
                    (*ct)->intern_type = "O.3so";
                    (*ct)->ext->hybridization = 3;
                } else {
                    (*ct)->intern_type = "O.2so";
                    (*ct)->ext->hybridization = 2;
                }
            } else if ((*ct)->element == "S") {
                if ((*ct)->bonded_atoms.size() == 2) {
                    if ((*ct)->ext->n_heavy_bonded == 1) {
                        (*ct)->intern_type = "S.sh";
                        (*ct)->ext->hybridization = 3;
                    } else {
                        (*ct)->intern_type = "S.3";
                        (*ct)->ext->hybridization = 3;
                    }
                } else {
                    (*ct)->intern_type = "S.2";
                    (*ct)->ext->hybridization = 2;
                }
            }
        }
    }

    if (O_count == 4) {
        if (depro) base->intern_type = "S.o4";
        else base->intern_type = "S.o4h";
        base->ext->hybridization = 2;
        return true;
    } else if (O_count == 3) {
        if (depro) base->intern_type = "S.o3";
        else base->intern_type = "S.o3h";
        base->ext->hybridization = 2;
        return true;
    } else if (O_count == 2) {
        if (depro) base->intern_type = "S.o2";
        else base->intern_type = "S.o2h";
        base->ext->hybridization = 2;
        return true;
    } else if (O_count == 1) {
        base->intern_type = "S.o";
        base->ext->hybridization = 2;
        return true;
    } else return false;
}


bool MOLECULE::check_PO(stl_ptr<ATOM> const& base,int const& O_count,int const& OH_count,int const& O_2heavy) {
    bool depro = false;
    int diff = O_count - O_2heavy;
    if (diff > 0) {
        if ((diff-OH_count) > 1) { // deprotoniert
            for (atoms_vec ct=base->bonded_atoms.begin(); ct!=base->bonded_atoms.end(); ++ct) {
                if ((*ct)->element == "O") {
                    if ((*ct)->bonded_atoms.size() == 2) {
                        if ((*ct)->ext->n_heavy_bonded == 1) {
                            (*ct)->intern_type = "O.3ac";
                            (*ct)->ext->hybridization = 3;
                        } else {
                            (*ct)->intern_type = "O.3po";
                            (*ct)->ext->hybridization = 3;
                        }
                    } else {
                        (*ct)->intern_type = "O.2p";
                        (*ct)->ext->hybridization = 2;
                    }
                } else if ((*ct)->element == "S") {
                    if ((*ct)->bonded_atoms.size() == 2) {
                        if ((*ct)->ext->n_heavy_bonded == 1) {
                            (*ct)->intern_type = "S.sh";
                            (*ct)->ext->hybridization = 3;
                        } else {
                            (*ct)->intern_type = "S.3";
                            (*ct)->ext->hybridization = 3;
                        }
                    } else {
                        (*ct)->intern_type = "S.2";
                        (*ct)->ext->hybridization = 2;
                    }
                }
            }
            depro = true;
        } else {
            for (atoms_vec ct=base->bonded_atoms.begin(); ct!=base->bonded_atoms.end(); ++ct) {
                if ((*ct)->element == "O") {
                    if ((*ct)->bonded_atoms.size() == 2) {
                        if ((*ct)->ext->n_heavy_bonded == 1) {
                            (*ct)->intern_type = "O.3oh";
                            (*ct)->ext->hybridization = 3;
                        } else {
                            (*ct)->intern_type = "O.3po";
                            (*ct)->ext->hybridization = 3;
                        }
                    } else {
                        (*ct)->intern_type = "O.2po";
                        (*ct)->ext->hybridization = 2;
                    }
                } else if ((*ct)->element == "S") {
                    if ((*ct)->bonded_atoms.size() == 2) {
                        if ((*ct)->ext->n_heavy_bonded == 1) {
                            (*ct)->intern_type = "S.sh";
                            (*ct)->ext->hybridization = 3;
                        } else {
                            (*ct)->intern_type = "S.3";
                            (*ct)->ext->hybridization = 3;
                        }
                    } else {
                        (*ct)->intern_type = "S.2";
                        (*ct)->ext->hybridization = 2;
                    }
                }
            }
        }
    } else {
        for (atoms_vec ct=base->bonded_atoms.begin(); ct!=base->bonded_atoms.end(); ++ct) {
            if ((*ct)->element == "O") {
                if ((*ct)->bonded_atoms.size() == 2) {
                    (*ct)->intern_type = "O.3po";
                    (*ct)->ext->hybridization = 3;
                } else {
                    (*ct)->intern_type = "O.2po";
                    (*ct)->ext->hybridization = 2;
                }
            } else if ((*ct)->element == "S") {
                if ((*ct)->bonded_atoms.size() == 2) {
                    if ((*ct)->ext->n_heavy_bonded == 1) {
                        (*ct)->intern_type = "S.sh";
                        (*ct)->ext->hybridization = 3;
                    } else {
                        (*ct)->intern_type = "S.3";
                        (*ct)->ext->hybridization = 3;
                    }
                } else {
                    (*ct)->intern_type = "S.2";
                    (*ct)->ext->hybridization = 2;
                }
            }
        }
    }

    if (O_count == 4) {
        if (depro) base->intern_type = "P.o4";
        else base->intern_type = "P.o4h";
        return true;
    } else if (O_count == 3) {
        if (depro) base->intern_type = "P.o3";
        else base->intern_type = "P.o3h";
        return true;
    } else if (O_count == 2) {
        if (depro) base->intern_type = "P.o2";
        else base->intern_type = "P.o2h";
        return true;
    } else if (O_count == 1) {
        base->intern_type = "P.o";
        return true;
    } else return false;
}


void MOLECULE::set_intern_types(bool const& prot_acids,bool const& prot_guanidin,bool const& prot_amidin,
                                bool const& prot_amin,bool const& prot_phosphate,bool const& prot_sulfate,
                                bool const& kekulize_charged) {
        //!Die internen Typen koennen spaeter auf einfachere Typen runtergebrochen werden
    
    bool breaker;
    
    //!************************************************************************************************
    //!zunaechst alle Typen fuer die Ringatome setzen:
    //!-----------------------------------------------
    //!folgende Typen sind hier moeglich (gleichzeitig Prioritaetenreihenfolge):
    //! C.ar6p    6-Ring aro. C mit geladener Resonanzstruktur (ortho oder para zu N.ar6p)
    //! C.ar6x    6-Ring heteroaro. C
    //! C.ar6    Benzolring
    //! C.arp    nicht-6-Ring aro. C mit geladener Resonanzstruktur
    //! C.arx    nicht-6-Ring heteroaro. C
    //! C.ar    andere nicht-6-Ring Aromaten
    //! C.2r3o    Carbonylkohlenstoff im 3-Ring
    //! C.2r3x    sp2 im Hetero-3-Ring
    //! C.3r3x    sp3 im Hetero-3-Ring
    //! C.2r3    sp2 im 3-Ring
    //! C.3r3    sp3 im 3-Ring
    //! N.ar6p    geladener N im 6-Ring Aromaten (z.B. Pyridinium oder NAD+)
    //! N.ar6    N im 6-Ring Aromaten
    //! N.arp    geladener N im nicht-6-Ring Aromaten (z.B. protoniertes Imidazol)
    //! N.ar2    N im nicht-6-Ring Aromaten mit 2 bonded_atoms (dem Sybyltyp N.2 entsprechend)
    //! N.ar3    N im nicht-6-Ring Aromaten mit 3 heavy atoms (dem Sybyltyp N.pl3 entsprechend)
    //! N.ar3h    N im nicht-6-Ring Aromaten mit 2 heavy atoms (dem Sybyltyp N.pl3 entsprechend)
    //! N.r3    N im Aziridin oder Aziren
    //! O.ar    aromatischer Sauerstoff (z.B. Furan)
    //! O.r3    O im Oxiran
    //! S.ar    aromatischer Schwefel (z.B. Thiophen)
    //! S.r3    S im Thiiran
    //! P.r3    P im Phosphiran
    //!
    //!folgende Typen muessen spaeter noch gesetzt werden:
//    //! N.agu1    aro. N (N.2) welches Teil einer Guanidinstruktur ist
//    //! N.agu2    aro  NH (N.pl3) welches Teil einer Guanidinstruktur ist
//    //! N.aguh    aro. N welches Teil einer protonierten Guanidinstruktur ist
    
    for (rings_vec rit=rings.begin(); rit!=rings.end(); ++rit) {
        if ((*rit)->is_aromatic) {
            if ((*rit)->n_members == 6) {
                if ((*rit)->n_hetero > 0) { //!6-Ring Heteroaromaten
                    if ((*rit)->is_pos) { //positiv geladen
                        for (atoms_vec at=(*rit)->ring_atoms.begin(); at!=(*rit)->ring_atoms.end(); ++at) {
                            if ((*at)->intern_type != "X") continue;
                            //!Der geladene Stickstoff wurde bereits in check_aro() gesetzt!
                            if ((*at)->element == "C") {
                                //wenn mit pos. Resonanzstruktur, dann C.ar6p sonst C.ar6x:
                                if ((*at)->bonded_atoms[0]->intern_type == "N.ar6p" ||
                                    (*at)->bonded_atoms[1]->intern_type == "N.ar6p") {
                                    (*at)->intern_type = "C.ar6p";
                                } else {
                                    for (atoms_vec bt=(*at)->bonded_atoms.begin();
                                                   bt!=(*at)->bonded_atoms.end(); ++bt) {
                                        for (atoms_vec ct=(*bt)->bonded_atoms.begin();
                                                       ct!=(*bt)->bonded_atoms.end(); ++ct) {
                                            if ((*ct)->intern_type == "N.ar6p")
                                                (*at)->intern_type = "C.ar6x";
                                        }
                                    }
                                }
                                if ((*at)->intern_type != "C.ar6x") (*at)->intern_type = "C.ar6p";
                            }
                            else if ((*at)->element == "N") {
                                if ((*at)->intern_type != "N.ar6p") (*at)->intern_type = "N.ar6";
                            }
                            else if ((*at)->element == "O") (*at)->intern_type = "O.ar";
                            else if ((*at)->element == "S") (*at)->intern_type = "S.ar";
                        }
                    } else {
                        for (atoms_vec at=(*rit)->ring_atoms.begin(); at!=(*rit)->ring_atoms.end(); ++at) {
                            if ((*at)->intern_type == "C.ar6p" || (*at)->intern_type == "N.ar6p") continue;
                            if ((*at)->element == "C") (*at)->intern_type = "C.ar6x";
                            else if ((*at)->element == "N") (*at)->intern_type = "N.ar6";
                        }
                    }
                } else { //!Benzolring
                    for (atoms_vec at=(*rit)->ring_atoms.begin(); at!=(*rit)->ring_atoms.end(); ++at) {
                        if ((*at)->intern_type == "C.ar6x" || (*at)->intern_type == "C.ar6p") continue;
                        (*at)->intern_type = "C.ar6";
                    }
                }
            } else {
                if ((*rit)->n_hetero > 0) { //!andere Heteroaromaten
                    if ((*rit)->is_pos) { //positiv geladen
                        for (atoms_vec at=(*rit)->ring_atoms.begin(); at!=(*rit)->ring_atoms.end(); ++at) {
                            if ((*at)->intern_type != "X") continue;
                            if ((*at)->intern_type == "C.ar6x" || (*at)->intern_type == "C.ar6p" ||
                                (*at)->intern_type == "C.ar6") continue;
                            //!Die N-Typen fuer protonierte 5-Ringe wurden bereits in check_aro() gesetzt !!!
                            // (N.arp und N.ar2)
                            if ((*at)->element == "C") {
                                //C.arp wenn zwischen den geladenen N's, sonst C.arx:
                                //wenn die geladenen N's nebeneinander liegen, dann die benachbarten C's auf C.arp:
                                bool N_N = false;
                                for (atoms_vec bt=(*rit)->ring_atoms.begin(); bt!=(*rit)->ring_atoms.end(); ++bt) {
                                    if ((*bt)->intern_type == "N.arp") {
                                        if ((*bt)->bonded_atoms[0]->intern_type == "N.arp" ||
                                            (*bt)->bonded_atoms[1]->intern_type == "N.arp") {
                                                N_N = true;
                                                break;
                                            }
                                    }
                                }
                                if (N_N) {
                                    if ((*at)->bonded_atoms[0]->intern_type == "N.arp" ||
                                        (*at)->bonded_atoms[1]->intern_type == "N.arp") {
                                        (*at)->intern_type = "C.arp";
                                    } else (*at)->intern_type = "C.arx";
                                } else {
                                    if ((*at)->bonded_atoms[0]->intern_type == "N.arp" &&
                                        (*at)->bonded_atoms[1]->intern_type == "N.arp") {
                                        (*at)->intern_type = "C.arp";
                                    } else (*at)->intern_type = "C.arx";
                                }
                            }
                            else if ((*at)->element == "O") (*at)->intern_type = "O.ar";
                            else if ((*at)->element == "S") (*at)->intern_type = "S.ar";
                        }
                    } else {
                        for (atoms_vec at=(*rit)->ring_atoms.begin(); at!=(*rit)->ring_atoms.end(); ++at) {
                            if ((*at)->intern_type != "X") continue;
                            if ((*at)->intern_type == "C.ar6x" || (*at)->intern_type == "C.ar6p" ||
                                (*at)->intern_type == "C.ar6") continue;
                            if ((*at)->element == "C") (*at)->intern_type = "C.arx";
                            else if ((*at)->element == "N") {
                                if ((*at)->sybyl_type == "N.pl3" ||
                                    (*at)->bonded_atoms.size() == 3) {
                                    if ((*at)->ext->n_heavy_bonded == 3) (*at)->intern_type = "N.ar3";
                                    else (*at)->intern_type = "N.ar3h";
                                } else (*at)->intern_type = "N.ar2";
                            }
                            else if ((*at)->element == "O") (*at)->intern_type = "O.ar";
                            else if ((*at)->element == "S") (*at)->intern_type = "S.ar";
                        }
                    }
                } else { //!nicht-6-Ring-nicht-Heteroaromat
                    for (atoms_vec at=(*rit)->ring_atoms.begin();
                         at!=(*rit)->ring_atoms.end(); ++at) {
                             if ((*at)->intern_type == "C.ar6x" || (*at)->intern_type == "C.ar6p" ||
                            (*at)->intern_type == "C.ar6" || (*at)->intern_type == "C.arx" ||
                            (*at)->intern_type == "C.arp") continue;
                             (*at)->intern_type = "C.ar";
                    }
                }
            }
        } else if ((*rit)->n_members == 3) { //! 3-Ringe
            for (atoms_vec at=(*rit)->ring_atoms.begin(); at!=(*rit)->ring_atoms.end(); ++at) {
                if ((*at)->element == "N") (*at)->intern_type = "N.r3";
                else if ((*at)->element == "O") (*at)->intern_type = "O.r3";
                else if ((*at)->element == "S") {
                    (*at)->intern_type = "S.r3";
                    (*at)->ext->hybridization = 3;
                } else if ((*at)->element == "P") {
                    (*at)->intern_type = "P.r3";
                    (*at)->ext->hybridization = 3;
                } else if ((*at)->element == "C") { //auf Carbonyl pruefen
                    for (atoms_vec bt=(*at)->bonded_atoms.begin(); bt!=(*at)->bonded_atoms.end(); ++bt) {
                        if ((*bt)->element == "O" && (*bt)->bonded_atoms.size() == 1) {
                            (*at)->intern_type = "C.2r3o"; //sp2 im Cyclopropanon
                            (*bt)->intern_type = "O.carb";
                            break;
                        }
                    }
                    if ((*at)->intern_type == "C.2r3o") continue;
                    if ((*rit)->n_hetero > 0) {
                        if ((*at)->ext->hybridization == 2) (*at)->intern_type = "C.2r3x";
                        else (*at)->intern_type = "C.3r3x";
                    } else {
                        if ((*at)->ext->hybridization == 2) (*at)->intern_type = "C.2r3";
                        else (*at)->intern_type = "C.3r3";
                    }
                } else (*at)->intern_type = (*at)->element;
            }
        } else if ((*rit)->is_planar) {
            for (atoms_vec at=(*rit)->ring_atoms.begin(); at!=(*rit)->ring_atoms.end(); ++at) {
                if ((*at)->element == "O") (*at)->ext->hybridization = 3; //wurde urspr. auf sp2 gesetzt fuer aro. Erk.
                //auch N wurde lieber auf sp2 gelassen, kann aber nicht generell zurueckgesetzt werden!!!
            }
        }
    }
    //!************************************************************************************************
    
    for (atoms_vec at=atoms.begin(); at!=atoms.end(); ++at) {
        if ((*at)->intern_type != "X") continue;
        
        if ((*at)->ext->hybridization == 0) { //sonstige Elemente
            (*at)->intern_type = (*at)->element; //!Achtung: Hier wird auch der intern_type von Wasserstoff gesetzt
            continue;
        }
        
        //!****************************************************************************************
        //!Kohlenstoff  1. Durchlauf:
        //!--------------------------
        //!folgende Typen sind hier moeglich:
        //! C.1n    sp in Cyano-Gruppe
        //! C.1p    sp mit einem heavy atom
        //! C.1s    sp mit 2 heavy atoms
        //! C.co2h    in protonierter COOH Gruppe (nur mit explizitem Proton)
        //! C.co2    in deprotonierter COO-  Gruppe
        //! C.es    Carbonyl-C im Ester oder Anhydrid
        //! C.hal    Carbonyl-C im Saeurehalogenid
        //! C.am    Carbonyl-C im Amid
        //! C.o        sonstiger Carbonylkohlenstoff
        //! C.s        sp2 in Thionylgruppe
        //! C.gu    sp2 in unprotonierter Guanidinogruppe
        //! C.guh    sp2 in protonierter Guanidinogruppe (bzw. bei unbekanntem Protonierungsstatus)
        //! C.mi    sp2 in unprotonierter Amidinogruppe
        //! C.mih    sp2 in protonierter Amidinogruppe (bzw. bei unbekanntem Protonierungsstatus)
        //! C.2p    sonstige sp2 mit einem heavy atom
        //! C.2s    sonstige sp2 mit 2 heavy atoms
        //! C.2t    sonstige sp2 mit 3 heavy atoms
        //! C.et    sp3 im Ether
        //! C.ohp    sp3 im primaeren Alkohol
        //! C.ohs    sp3 im sekundaeren Alkohol
        //! C.oht    sp3 im tertiaeren Alkohol
        //! C.3n    sp3 an Stickstoff gebunden
        //! C.3p    sonstiger sp3 mit einem heavy atom
        //! C.3s    sonstiger sp3 mit 2 heavy atoms
        //! C.3t    sonstiger sp3 mit 3 heavy atoms
        //! C.3q    sonstiger sp3 mit 4 heavy atoms
        //!
        //!folgende Typen muessen spaeter noch gesetzt werden:
        //! C.n        sp2 im Imin
        if ((*at)->element == "C") {
            if ((*at)->ext->hybridization == 1) { //! sp  Kohlenstoff
                if ((*at)->ext->n_heavy_bonded == 1) {
                    if ((*at)->bonded_atoms[0]->element == "N") {
                        (*at)->intern_type = "C.1n";
                        continue;
                    }
                    (*at)->intern_type = "C.1p";
                    continue;
                }
                (*at)->intern_type = "C.1s";
                continue;
            } else if ((*at)->ext->hybridization == 2) { //! sp2  Kohlenstoff
                int O_sp2 = 0;
                int O_sp3 = 0;
                int O_oh = 0;
                int O_et = 0;
                int n_N = 0;
                int N_sp2 = 0;
                int N_single = 0;
                int N_pl3 = 0;
                int S_single = 0;
                int n_Hal = 0;
                for (atoms_vec bt=(*at)->bonded_atoms.begin(); bt!=(*at)->bonded_atoms.end(); ++bt) {
                    if ((*bt)->element == "N") {
                        ++n_N;
                        if ((*bt)->bonded_atoms.size() == 3 && (*bt)->ext->hybridization == 2) ++N_pl3;
                        else if ((*bt)->ext->hybridization == 2) {
                            ++N_sp2;
                            if ((*bt)->ext->n_heavy_bonded == 1) ++N_single;
                        } else if ((*bt)->ext->n_heavy_bonded == 1) ++N_single;
                    }
                    else if ((*bt)->element == "O") {
                        if ((*bt)->ext->hybridization == 2) ++O_sp2;
                        else if ((*bt)->ext->hybridization == 3) {
                            ++O_sp3;
                            if ((*bt)->ext->n_heavy_bonded == 2) ++O_et;
                            else if ((*bt)->bonded_atoms.size() == 2 && (*bt)->ext->n_heavy_bonded == 1) ++O_oh;
                        }
                    }
                    else if ((*bt)->element == "S") {
                        if ((*bt)->bonded_atoms.size() == 1) ++S_single;
                    }
                    else if ((*bt)->element == "F" || (*bt)->element == "Cl" 
                             || (*bt)->element == "Br" || (*bt)->element == "I") ++n_Hal;
                }

                if (O_sp2 == 1 && O_oh == 1) {
                    (*at)->intern_type = "C.co2h"; //protoniert
                    continue;
                }
                if (O_sp2 == 1 && O_et > 0) {
                    (*at)->intern_type = "C.es"; 
                    continue;
                }
                if (O_sp2 == 1 && n_Hal == 1) {
                    (*at)->intern_type = "C.hal";
                    continue;
                }
                if (O_sp2 == 2 || (O_sp2 == 1 && O_sp3 == 1)) {
                    (*at)->intern_type = "C.co2"; //deprotoniert
                    continue;
                }
                if (O_sp2 == 1 && n_N > 0) {
                    (*at)->intern_type = "C.am";
                    continue;
                }
                if (O_sp2 == 1) {
                    (*at)->intern_type = "C.o";
                    continue;
                }
                if (S_single == 1) {
                    (*at)->intern_type = "C.s";
                    continue;
                }
                if (n_N == 3) { //kann eine Guanidinogruppe sein
                    if (N_pl3 == 2) (*at)->intern_type = "C.gu";
                    else (*at)->intern_type = "C.guh"; //!wenn Protonierungsgrad nicht explizit vorgegeben immer C.guh
                    continue;
                }
                if (n_N == 2) { //koennte Amidinogruppe sein --- muss aber nicht!
                    if (N_pl3 == 1 && N_sp2 == 1) {
                        (*at)->intern_type = "C.mi";
                        continue;
                    }
                    else if (N_pl3 == 2 || N_sp2 == 2 || N_single == 2) {
                        (*at)->intern_type = "C.mih"; //!wenn Protonierungsgrad nicht explizit vorgegeben immer C.mih
                        continue;
                    }
                }
                if ((*at)->ext->n_heavy_bonded == 1) {
                    bool to_aro = false;
                    for (atoms_vec bt=(*at)->bonded_atoms.begin(); bt!=(*at)->bonded_atoms.end(); ++bt) {
                        if ((*bt)->ext->is_aromatic) to_aro = true;
                    }
                    if (to_aro) {
                        (*at)->intern_type = "C.3p";
                        (*at)->ext->hybridization = 3;
                    }
                    else if ((*at)->ext->sp2 > (*at)->ext->sp3) {(*at)->intern_type = "C.2p";}
                    else { //wahrschl. doch eher ein sp3
                        bool sp3_flag = false;
                        for (atoms_vec bt=(*at)->bonded_atoms[0]->bonded_atoms.begin(); 
                             bt!=(*at)->bonded_atoms[0]->bonded_atoms.end(); ++bt) {
                            if ((*bt != *at) && (*bt)->ext->hybridization == 2) sp3_flag = true;
                        }
                        if (sp3_flag) {
                            (*at)->intern_type = "C.3p";
                            (*at)->ext->hybridization = 3;
                        } else (*at)->intern_type = "C.2p";
                    }
                    continue;
                } else if ((*at)->ext->n_heavy_bonded == 2) {
                    (*at)->intern_type = "C.2s";
                    continue;
                } else if ((*at)->ext->n_heavy_bonded == 3) {
                    (*at)->intern_type = "C.2t";
                    continue;
                }
                //!Die konjugierten Typen werden erst spaeter gesetzt!!!
            } else if ((*at)->ext->hybridization == 3) { //! sp3  Kohlenstoff
                int n_N = 0;
                int n_O = 0;
                int O_et = 0;
                int n_OH = 0;
                for (atoms_vec bt=(*at)->bonded_atoms.begin(); bt!=(*at)->bonded_atoms.end(); ++bt) {
                    if ((*bt)->element == "N") ++n_N;
                    else if ((*bt)->element == "O") {
                        ++n_O;
                        if ((*bt)->ext->n_heavy_bonded == 2) ++O_et;
                        else ++n_OH;
                    }
                }
                if (O_et > 0) {
                    if (n_O > O_et && check_trigo(*at,1.0)) {
                        bool no_ester = true;
                        for (atoms_vec bt=(*at)->bonded_atoms.begin(); bt!=(*at)->bonded_atoms.end(); ++bt) {
                            if ((*bt)->element == "O" && (*bt)->bonded_atoms.size() == 1) {
                                (*at)->intern_type = "C.es";
                                (*at)->ext->hybridization = 2;
                                (*bt)->intern_type = "O.2es";
                                (*bt)->ext->hybridization = 2;
                                no_ester = false;
                                break;
                            }
                        }
                        if (no_ester) {
                            (*at)->intern_type = "C.et";
                        }
                    } else {
                        (*at)->intern_type = "C.et";
                    }
                    continue;
                }
                if (n_OH > 0) {
                    if ((*at)->ext->n_heavy_bonded == 2) (*at)->intern_type = "C.ohp";
                    else if ((*at)->ext->n_heavy_bonded == 3) (*at)->intern_type = "C.ohs";
                    else if ((*at)->ext->n_heavy_bonded == 4) (*at)->intern_type = "C.oht";
                    continue;
                }
                if (n_N > 0) {
                    (*at)->intern_type = "C.3n";
                    continue;
                }
                if ((*at)->ext->n_heavy_bonded == 1) {
                    (*at)->intern_type = "C.3p";
                    continue;
                } else if ((*at)->ext->n_heavy_bonded == 2) {
                    (*at)->intern_type = "C.3s";
                    continue;
                } else if ((*at)->ext->n_heavy_bonded == 3) {
                    (*at)->intern_type = "C.3t";
                    continue;
                } else if ((*at)->ext->n_heavy_bonded == 4) {
                    (*at)->intern_type = "C.3q";
                    continue;
                }
            }
        }
        //!****************************************************************************************
        
        //!****************************************************************************************
        //!Stickstoff  1. Durchlauf:
        //!--------------------------
        //!folgende Typen sind hier moeglich:
        //! N.az    mittlers N im Azid
        //! N.1        sonstiger sp Stickstoff
        //! N.aap    primaeres aromatisches Amin (Hybridisierung kann nicht eindeutig bestimmt werden)
        //! N.aas2    sekundaeres aro. Amin planar
        //! N.aas3    sekundaeres aro. Amin nicht planar
        //! N.aat2    sp2 tertiaeres aro. Amin
        //! N.aat3    sp3 tertiaeres aro. Amin
        //! N.o2    Nitrogruppe
        //! N.ohac    Hydroxamsaeure
        //! N.oh    Hydroxylamin
        //! N.ims    Imid mit 2 heavy atoms
        //! N.imt    Imid mit 3 heavy atoms
        //! N.samp    Sulfonamid mit einem heavy atom
        //! N.sams    Sulfonamid mit 2 heavy atoms
        //! N.samt    Sulfonamid mit 3 heavy atoms
        //! N.amp    Carbonsaeure- oder Thioamid mit einem heavy atom
        //! N.ams    Carbonsaeure- oder Thioamid mit 2 heavy atoms
        //! N.amt    Carbonsaeure- oder Thioamid mit 3 heavy atoms
        //! N.gu1    NH in unprotonierter Guanidinogruppe (nur bei expliziter Protonierung)
        //! N.gu2    NH2 in unprotonierter Guanidinogruppe (nur bei expliziter Protonierung)
        //! N.guh    N in protonierter Guanidinogruppe (bzw. bei unklarem Protonierungsstatus)
        //! N.mi1    NH in unprotonierter Amidinogruppe (nur bei expliziter Protonierung)
        //! N.mi2    NH2 in unprotonierter Amidinogruppe (nur bei expliziter Protonierung)
        //! N.mih    N in protonierter Amidinogruppe (bzw. bei unklarem Protonierungsstatus)
        //! N.2n    sp2 N dass an anderes N gebunden ist
        //! N.2p    sonstige sp2 mit einem heavy atom
        //! N.2s    sonstige sp2 mit 2 heavy atoms
        //! N.4q    sp3 mit 4 bonded heavy atoms
        //! N.4h    sp3 mit 4 bonded atoms, dabei min. 1 Wasserstoff
        //! N.3n    sp3 N an anderes N gebunden
        //! N.3p    sp3 mit einem heavy atom
        //! N.3s    sp3 mit 2 heavy atoms
        //! N.3t    sp3 mit 3 heavy atoms
        if ((*at)->element == "N") {
            if ((*at)->ext->hybridization == 1) { //! sp  Stickstoff
                if ((*at)->ext->n_heavy_bonded == 2) {
                    if ((*at)->bonded_atoms[0]->element == "N" &&
                        (*at)->bonded_atoms[1]->element == "N") {
                        (*at)->intern_type = "N.az";
                        continue;
                    }
                }
                (*at)->intern_type = "N.1";
                continue;
            }
            for (atoms_vec bt=(*at)->bonded_atoms.begin(); bt!=(*at)->bonded_atoms.end(); ++bt) {
                if ((*bt)->ext->is_aromatic) {
                    if ((*at)->ext->n_heavy_bonded == 1) (*at)->intern_type = "N.aap";
                    else if ((*at)->ext->n_heavy_bonded == 2) {
                        if (check_db_tors(*at)) {
                            (*at)->intern_type = "N.aas2";
                            (*at)->ext->hybridization = 2;
                        } else {
                            (*at)->intern_type = "N.aas3";
                            (*at)->ext->hybridization = 3;
                        }
                    } else if ((*at)->ext->n_heavy_bonded == 3) {
                        if (check_trigo(*at)) {
                            (*at)->intern_type = "N.aat2";
                            (*at)->ext->hybridization = 2;
                        } else {
                            (*at)->intern_type = "N.aat3";
                            (*at)->ext->hybridization = 3;
                        }
                    }
                    break;
                }
            }
            //!Trotzdem noch nicht abbrechen, da andere Typen wie z.B. Amid vorrang haben
            //!Erst beim setzen von N.3 und N.2 Typen abbrechen, wenn schon ein N.aa Typ gesetzt ist
            if ((*at)->ext->hybridization == 2) { //! sp2  Stickstoff
                //!Hier hab ich erstmal die alte Variante gelassen, obwohl mit den neuen C-Typen in einem
                //!Extralauf eine leichtere Typisierung moeglich waere.
                //!Andererseits bleibt so die Option fuer einen weiteren Plausibilitaetscheck offen.
                int O_sp2_count = 0;
                int O_count = 0;
                int NO_count = 0;
                int S_o2_count = 0;
                bool N_flag = false;
                bool reset_mi_to_n3 = false; //!03.04.07
                for (atoms_vec bt=(*at)->bonded_atoms.begin(); bt!=(*at)->bonded_atoms.end(); ++bt) {
                    if ((*bt)->element == "C") {
                        int N_count = 0;
                        int Npl3_count = 0;
                        bool no_mi = false; //!03.04.07
                        if ((*bt)->ext->hybridization == 3) no_mi = true; //!03.04.07
                        for (atoms_vec ct=(*bt)->bonded_atoms.begin(); ct!=(*bt)->bonded_atoms.end(); ++ct) {
                            if ((*ct)->element == "N") {
                                ++N_count;
                                if ((*ct)->bonded_atoms.size() == 3) ++Npl3_count;
                            }
                            if ((*ct)->element == "O" && (*ct)->bonded_atoms.size() == 1) {
                                ++O_sp2_count; //Carbonylgruppe
                            }
                            if ((*ct)->element == "S" && (*ct)->bonded_atoms.size() == 1) {
                                ++O_sp2_count; //Thionylgruppe wird hier wie Carbonylgruppe behandelt
                            }
                        }
                        if (N_count == 3) {
                            if (Npl3_count == 2) {
                                if ((*at)->bonded_atoms.size() == 3) (*at)->intern_type = "N.gu2"; //! NH2 im Guan.
                                else (*at)->intern_type = "N.gu1"; //! NH im Guanidin
                            }
                            else {
                                (*at)->intern_type = "N.guh"; //!auch bei unbekanntem Prot.-Status
                            }
                            //!kein break, weil es auch noch ein Amid sein koennte
                        }
                        if (N_count == 2 && (*at)->intern_type == "X") { //!Amidin nicht setzen, wenn schon Guanidin
                            if (Npl3_count == 1) {
                                if ((*at)->bonded_atoms.size() == 3) (*at)->intern_type = "N.mi2"; //! NH2 im Ami.
                                else (*at)->intern_type = "N.mi1"; //! NH im Amidin
                            }
                            else (*at)->intern_type = "N.mih"; //!auch bei unbekanntem Prot.-Status
                            //!03.04.2007: Pruefen, ob wirklich ein planares C zwischen den N's
                            //!-> sonst ist es eher ein N.3-Typ:
                            if (no_mi) {
                                reset_mi_to_n3 = true;
                                (*at)->intern_type = "X";
                            } else reset_mi_to_n3 = false;
                        }
                        
                    }
                    else if ((*bt)->element == "S") {
                        for (atoms_vec ct=(*bt)->bonded_atoms.begin(); ct!=(*bt)->bonded_atoms.end(); ++ct) {
                            if ((*ct)->element == "O" && (*ct)->bonded_atoms.size() == 1) {
                                ++S_o2_count;
                                break;
                            }
                        }
                    }
                    else if ((*bt)->element == "O") {
                        ++O_count;
                        if ((*bt)->bonded_atoms.size() == 1) ++NO_count;
                    }
                    else if ((*bt)->element == "N") {
                        N_flag = true;
                    }
                }
                if (NO_count == 2) {
                    (*at)->intern_type = "N.o2"; //!Nitrogruppe (hat hoechste Prioritaet)
                    continue;
                }
                if (O_count == 1) {
                    if (O_sp2_count == 1) {
                        (*at)->intern_type = "N.ohac"; //!Hydroxamsaeure
                        for (atoms_vec bt=(*at)->bonded_atoms.begin(); bt!=(*at)->bonded_atoms.end(); ++bt) {
                            if ((*bt)->element == "C") {
                                for (atoms_vec ct=(*bt)->bonded_atoms.begin(); ct!=(*bt)->bonded_atoms.end(); ++ct) {
                                    if ((*ct)->element == "O" && (*ct)->bonded_atoms.size() == 1) {
                                        (*ct)->intern_type = "O.am";
                                        (*bt)->intern_type = "C.am";
                                        break;
                                    }
                                }
                            }
                        }
                    }
                    else (*at)->intern_type = "N.oh"; //!Hydroxylamin
                    continue;
                }
                if ((O_sp2_count + S_o2_count) > 1) {
                    if ((*at)->ext->n_heavy_bonded == 2) (*at)->intern_type = "N.ims"; //!Imid
                    else if ((*at)->ext->n_heavy_bonded == 3) (*at)->intern_type = "N.imt";
                    continue;
                }
                if (S_o2_count > 0) {
                    //!Sulfonamidtyp hat Vorrang gegenueber Carbonsaeureamid
                    if ((*at)->ext->n_heavy_bonded == 1) (*at)->intern_type = "N.samp"; //!primaeres Sulfonamid
                    else if ((*at)->ext->n_heavy_bonded == 2) (*at)->intern_type = "N.sams";
                    else if ((*at)->ext->n_heavy_bonded == 3) (*at)->intern_type = "N.samt";
                    continue;
                }
                if (O_sp2_count == 1) {
                    //!Amidtyp hat Vorrang gegenueber Guanidin und Amidin => entspr. Typ wird ueberschrieben
                    if ((*at)->ext->n_heavy_bonded == 1) (*at)->intern_type = "N.amp"; //!primaeres Amid
                    else if ((*at)->ext->n_heavy_bonded == 2) (*at)->intern_type = "N.ams";
                    else if ((*at)->ext->n_heavy_bonded == 3) (*at)->intern_type = "N.amt";
                    continue;
                }
                
                if ((*at)->intern_type != "X") continue;
                if (reset_mi_to_n3) { //!03.04.07
                    (*at)->ext->hybridization = 3;
                } else {
                    if (N_flag) {
                        (*at)->intern_type = "N.2n"; //!sp2 N, dass an anderes N gebunden ist
                        continue;
                    }
                    if ((*at)->ext->n_heavy_bonded == 1) {
                        (*at)->intern_type = "N.2p"; //!primaeres sp2 N
                        continue;
                    }
                    if ((*at)->ext->n_heavy_bonded == 2) {
                        (*at)->intern_type = "N.2s"; //!sekundaeres sp2 N
                        continue;
                    }
                    if ((*at)->ext->n_heavy_bonded == 3) {
                        (*at)->intern_type = "N.2t"; //!gibt es auch an nicht-Aromaten (z.B. IDIHIN), aber so kommt es der Sache Nahe
                        continue;
                    }
                }
            }
            if ((*at)->ext->hybridization == 3) { //! sp3  Stickstoff
                if ((*at)->bonded_atoms.size() == 4) {
                    if ((*at)->ext->n_heavy_bonded == 4) (*at)->intern_type = "N.4q";
                    else (*at)->intern_type = "N.4h";
                    continue;
                }
                //!und jetzt nochmal pruefen auf N.am, N.sam, N.im, N.ohac, N.oh:
                int O_sp2_count = 0;
                int O_count = 0;
                int S_o2_count = 0;
                bool N_flag = false;
                for (atoms_vec bt=(*at)->bonded_atoms.begin(); bt!=(*at)->bonded_atoms.end(); ++bt) {
                    if ((*bt)->element == "C") {
                        int N_count = 0;
                        for (atoms_vec ct=(*bt)->bonded_atoms.begin(); ct!=(*bt)->bonded_atoms.end(); ++ct) {
                            if ((*ct)->element == "N") {
                                ++N_count;
                            }
                            if ((*ct)->element == "O" && (*ct)->bonded_atoms.size() == 1) {
                                ++O_sp2_count; //Carbonylgruppe
                            }
                            if ((*ct)->element == "S" && (*ct)->bonded_atoms.size() == 1) {
                                ++O_sp2_count; //Thionylgruppe wird hier wie Carbonylgruppe behandelt
                            }
                        }
                    }
                    else if ((*bt)->element == "S") {
                        for (atoms_vec ct=(*bt)->bonded_atoms.begin(); ct!=(*bt)->bonded_atoms.end(); ++ct) {
                            if ((*ct)->element == "O" && (*ct)->bonded_atoms.size() == 1) {
                                ++S_o2_count;
                                break;
                            }
                        }
                    }
                    else if ((*bt)->element == "N") {
                        N_flag = true;
                    }
                    else if ((*bt)->element == "O") { //!neu: 14.11.2008
                        ++O_count;
                    }
                }
                if (O_count == 1) {
                    if (O_sp2_count == 1) {
                        (*at)->intern_type = "N.ohac"; //!Hydroxamsaeure
                        for (atoms_vec bt=(*at)->bonded_atoms.begin(); bt!=(*at)->bonded_atoms.end(); ++bt) {
                            if ((*bt)->element == "C") {
                                for (atoms_vec ct=(*bt)->bonded_atoms.begin(); ct!=(*bt)->bonded_atoms.end(); ++ct) {
                                    if ((*ct)->element == "O" && (*ct)->bonded_atoms.size() == 1) {
                                        (*ct)->intern_type = "O.am";
                                        (*bt)->intern_type = "C.am";
                                        break;
                                    }
                                }
                            }
                        }
                    }
                    else (*at)->intern_type = "N.oh"; //!Hydroxylamin
                    continue;
                }
                if ((O_sp2_count + S_o2_count) > 1) {
                    if ((*at)->ext->n_heavy_bonded == 2) (*at)->intern_type = "N.ims"; //!Imid
                    else if ((*at)->ext->n_heavy_bonded == 3) (*at)->intern_type = "N.imt";
                    continue;
                }
                if (S_o2_count > 0) {
                    //!Sulfonamidtyp hat Vorrang gegenueber Carbonsaeureamid
                    if ((*at)->ext->n_heavy_bonded == 1) (*at)->intern_type = "N.samp"; //!primaeres Sulfonamid
                    else if ((*at)->ext->n_heavy_bonded == 2) (*at)->intern_type = "N.sams";
                    else if ((*at)->ext->n_heavy_bonded == 3) (*at)->intern_type = "N.samt";
                    continue;
                }
                if (O_sp2_count == 1) {
                    //!Amidtyp hat Vorrang gegenueber Guanidin und Amidin => entspr. Typ wird ueberschrieben
                    if ((*at)->ext->n_heavy_bonded == 1) (*at)->intern_type = "N.amp"; //!primaeres Amid
                    else if ((*at)->ext->n_heavy_bonded == 2) (*at)->intern_type = "N.ams";
                    else if ((*at)->ext->n_heavy_bonded == 3) (*at)->intern_type = "N.amt";
                    continue;
                }
                if ((*at)->intern_type != "X") continue;
                if (N_flag) {
                    (*at)->intern_type = "N.3n"; //!sp3 N, dass an anderes N gebunden ist
                    continue;
                }
                if ((*at)->ext->n_heavy_bonded == 1) {
                    (*at)->intern_type = "N.3p"; //!primaeres sp3 N
                    continue;
                }
                if ((*at)->ext->n_heavy_bonded == 2) {
                    (*at)->intern_type = "N.3s"; //!sekundaeres sp3 N
                    continue;
                }
                if ((*at)->ext->n_heavy_bonded == 3) {
                    (*at)->intern_type = "N.3t"; //!tertiaeres sp3 N
                    continue;
                }
            }
        }
        //!****************************************************************************************
    }
    
    //!Jetzt wo alle N und C gesetzt sind faellt das setzen der restlichen Atome deutlich leichter:
    for (atoms_vec at=atoms.begin(); at!=atoms.end(); ++at) {
        //!**************Iminkohlenstoff***********************************************************
        //!folgender Typ wird gesetzt:
        //! C.n        Iminkohlenstoff
        bool im_flag = false;
        bool no_other = true;
        if ((*at)->element == "C") {
            if ((*at)->intern_type == "C.2p" || (*at)->intern_type == "C.2s" || (*at)->intern_type == "C.2t") {
                for (atoms_vec bt=(*at)->bonded_atoms.begin(); bt!=(*at)->bonded_atoms.end(); ++bt) {
                    if ((*bt)->intern_type == "N.2p" || (*bt)->intern_type == "N.2s") {
                        //!pruefen ob es sich um eine DB handelt -> wenn ja, dann auf C.n setzen
                        vec3d<float> vv = (*at)->coord;
                        vv -= (*bt)->coord;
                        float dst = vv.value(); 
                        if (dst < 1.40) im_flag = true;
                    } else if ((*bt)->ext->hybridization == 2) {
                        if ((*bt)->element == "O" && (*bt)->bonded_atoms.size() > 1) continue;
                        vec3d<float> vv = (*at)->coord;
                        vv -= (*bt)->coord;
                        float dst = vv.value(); 
                        if (dst < 1.37) no_other = false;
                        else if (dst < 1.43 && !((*bt)->ext->is_aromatic)) for (atoms_vec ct=(*bt)->bonded_atoms.begin(); 
                                                  ct!=(*bt)->bonded_atoms.end(); ++ct) {
                            if ((*ct) != (*at)) {
                                if ((*ct)->ext->hybridization == 2) no_other = false;
                            }
                        }
                    }
                }
                if (im_flag && no_other) (*at)->intern_type = "C.n";
            }
        }
        
        int O_count, OH_count, N_count;
        //!****************************************************************************************
        
        //!****************************************************************************************
        //!Sauerstoff  2. Durchlauf:
        //!--------------------------
        //!folgende Typen sind hier moeglich:
        //! O.h2o    Wasser
        //! O.n        O in Nitrogruppe
        //! O.noh    O im Hydroxylamin oder in Hydroxamsaeure
        //! O.2es    sp2 O im Ester oder Anhydrid
        //! O.2hal    sp2 O im Saeurehalogenid
        //! O.am    Carbonsaeureamid
        //! O.2co2    sp2 in COOH (explizit protoniert)
        //! O.co2    sp2 in COO-
        //! O.carb    sonstige Carbonylgruppen
        //! O.2po    P=O in nicht deprotonierter Gruppe
        //! O.2p    P~O in deprotonierter Gruppe
        //! O.3po    sp3  in P--O--heavy  (Ester)
        //! O.2so    S=O in nicht deprotonierter Gruppe
        //! O.2s    S~O in deprotonierter Gruppe
        //! O.3so    sp3  in S--O--heavy  (Ester)
        //! O.o        im Peroxid
        //! O.3ac    OH-Gruppe in COOH, CSOH, POOHOH, POOH oder SOOOH
        //! O.ph    Phenolische OH-Gruppe
        //! O.3oh    normaler Alkohol
        //! O.3es    sp3 im Ester
        //! O.3eta    aromatischer Ether
        //! O.3et    Ether (aliphatisch)
        if ((*at)->element == "O") {
            
            if ((*at)->ext->n_heavy_bonded == 0) {
                (*at)->intern_type = "O.h2o"; //! Holy water !!!
                continue;
            }
            
            if ((*at)->intern_type != "X") continue;
            
            if ((*at)->ext->hybridization == 2 && (*at)->bonded_atoms.size() > 1) (*at)->ext->hybridization = 3;
            
            if ((*at)->ext->hybridization == 2) { //! sp2  Sauerstoff
                breaker = false;
                for (atoms_vec bt=(*at)->bonded_atoms.begin(); bt!=(*at)->bonded_atoms.end(); ++bt) {
                    if ((*bt)->element == "N") {
                        if ((*bt)->intern_type == "N.o2") {
                            (*at)->intern_type = "O.n"; //!Nitrogruppe
                            breaker = true;
                            break;
                        }
                        if ((*bt)->intern_type == "N.oh" || (*bt)->intern_type == "N.ohac") {
                            (*at)->intern_type = "O.noh"; //!Hydroxylamin oder Hydroxamsaeure
                            (*at)->ext->hybridization = 3;
                            breaker = true;
                            break;
                        }
                    }
                    if ((*bt)->element == "C") { //anders bei sp3
                        if ((*bt)->intern_type == "C.es") {
                            (*at)->intern_type = "O.2es";
                            breaker = true;
                            break;
                        }
                        if ((*bt)->intern_type == "C.hal") {
                            (*at)->intern_type = "O.2hal";
                            breaker = true;
                            break;
                        }
                        if ((*bt)->intern_type == "C.am") {
                            (*at)->intern_type = "O.am";
                            breaker = true;
                            break;
                        }
                        if ((*bt)->intern_type == "C.co2h") {
                            (*at)->intern_type = "O.2co2";
                            breaker = true;
                            break;
                        }
                        if ((*bt)->intern_type == "C.co2" || (*bt)->intern_type == "C.s") {
                            (*at)->intern_type = "O.co2";
                            breaker = true;
                            break;
                        }
                        else {
                            if ((*bt)->ext->is_aromatic) {
                                (*at)->intern_type = "O.ph";
                                (*at)->ext->hybridization = 3;
                                breaker = true;
                                break;
                            } else {
                                (*at)->intern_type = "O.carb";
                                breaker = true;
                                break;
                            }
                        }
                        
                    }
                    O_count = 0;
                    OH_count = 0;
                    N_count = 0;
                    int O_3s_count = 0;
                    if ((*bt)->element == "S") {
                        for (atoms_vec ct=(*bt)->bonded_atoms.begin(); ct!=(*bt)->bonded_atoms.end(); ++ct) {
                            if ((*ct)->element == "O" || (*ct)->element == "S") {
                                ++O_count;
                                if ((*ct)->bonded_atoms.size() == 2) {
                                    if ((*ct)->ext->n_heavy_bonded == 1) ++OH_count;
                                    else ++O_3s_count;
                                }
                            } else if ((*ct)->element == "N") {
                                ++N_count;
                                if ((*ct)->bonded_atoms.size() == 3) {
                                    if ((*ct)->ext->n_heavy_bonded == 3) ++O_3s_count; // NR3
                                    else ++OH_count; // NHR oder NH2
                                }
                            }
                        }
                        if (check_SO(*bt,O_count,OH_count,O_3s_count,N_count)) {
                            breaker = true;
                            break;
                        }
                    }

                    O_count = 0;
                    OH_count = 0;
                    int O_2heavy = 0;
                    if ((*bt)->element == "P") {
                        for (atoms_vec ct=(*bt)->bonded_atoms.begin(); ct!=(*bt)->bonded_atoms.end(); ++ct) {
                            if ((*ct)->element == "O" || (*ct)->element == "S") { // S = 200911 (BUHZEJ)
                                ++O_count;
                                if ((*ct)->bonded_atoms.size() == 2) {
                                    if ((*ct)->ext->n_heavy_bonded == 1) ++OH_count; // at -> ct 200911
                                    else ++O_2heavy;
                                }
                            }
                        }
                        if (check_PO(*bt,O_count,OH_count,O_2heavy)) {
                            breaker = true;
                            break;
                        }
                    } else if ((*bt)->element == "O") {
                        (*at)->intern_type = "O.o"; //!Peroxid
                        breaker = true;
                        break;
                    }
                }
                if (breaker) continue;
            }
            if ((*at)->ext->hybridization == 3) { //! sp3  Sauerstoff    
                if ((*at)->ext->n_heavy_bonded == 1) { //muss denn wohl ein OH, S=O oder P=O sein
                    breaker = false;
                    for (atoms_vec bt=(*at)->bonded_atoms.begin(); bt!=(*at)->bonded_atoms.end(); ++bt) {
                        if ((*bt)->element == "N") {
                            if ((*bt)->intern_type == "N.oh" || (*bt)->intern_type == "N.ohac") {
                                (*at)->intern_type = "O.noh"; //!Hydroxylamin oder Hydroxamsaeure
                                breaker = true;
                                break;
                            }
                        }
                        if ((*bt)->element == "C") {
                            if ((*bt)->intern_type == "C.co2") {
                                if ((*at)->bonded_atoms.size() > 1) {
                                    O_count = 0;
                                    OH_count = 0;
                                    int O_3s_count = 0;
                                    for (atoms_vec ct=(*bt)->bonded_atoms.begin(); ct!=(*bt)->bonded_atoms.end(); ++ct) {
                                        if ((*ct)->element == "O" || (*ct)->element == "S") {
                                            ++O_count;
                                            if ((*ct)->bonded_atoms.size() == 2) {
                                                if ((*ct)->ext->n_heavy_bonded == 1) ++OH_count;
                                                else ++O_3s_count;
                                            }
                                        }
                                    }
                                    if ((O_count - OH_count -O_3s_count) > 1) {
                                        if ((*at)->ext->n_heavy_bonded < 2) (*at)->intern_type = "O.3oh";
                                        else (*at)->intern_type = "O.3es";
                                    } else {
                                        (*bt)->intern_type == "C.co2h";
                                        if ((*at)->ext->n_heavy_bonded < 2) (*at)->intern_type = "O.3oh";
                                        else (*at)->intern_type = "O.3es";
                                    }
                                } else {
                                    (*at)->intern_type = "O.co2"; //!in COO-
                                    (*at)->ext->hybridization = 2;
                                }
                            } else if ((*bt)->intern_type == "C.o" || (*bt)->intern_type == "C.s" ||
                                (*bt)->intern_type == "C.co2h") {
                                bool carb_flag = false;
                                O_count = 0;
                                for (atoms_vec ct=(*bt)->bonded_atoms.begin(); ct!=(*bt)->bonded_atoms.end(); ++ct) {
                                    if (((*ct)->element == "O" || (*ct)->element == "S")) {//!auch Thiolsaeuren
                                        ++O_count;
                                        if ((*ct)->bonded_atoms.size() == 2 &&
                                            (*ct)->ext->n_heavy_bonded == 1) carb_flag = true;
                                    }
                                }
                                if ((O_count < 2) || carb_flag) (*at)->intern_type = "O.3ac"; //!in COOH
                                else {
                                    (*at)->intern_type = "O.co2"; //!in COO-
                                    (*at)->ext->hybridization = 2;
                                }
                                breaker = true;
                                break;
                            } else {
                                if ((*bt)->ext->is_aromatic) (*at)->intern_type = "O.ph";
                                else (*at)->intern_type = "O.3oh";
                                breaker = true;
                                break;
                            }
                        }
                        O_count = 0;
                        OH_count = 0;
                        N_count = 0;
                        int O_3s_count = 0;
                        if ((*bt)->element == "S") {
                                for (atoms_vec ct=(*bt)->bonded_atoms.begin(); ct!=(*bt)->bonded_atoms.end(); ++ct) {
                                if ((*ct)->element == "O" || (*ct)->element == "S") {
                                    ++O_count;
                                    if ((*ct)->bonded_atoms.size() == 2) {
                                        if ((*ct)->ext->n_heavy_bonded == 1) ++OH_count;
                                        else ++O_3s_count;
                                    }
                                } else if ((*ct)->element == "N") {
                                    ++N_count;
                                    if ((*ct)->bonded_atoms.size() == 3) {
                                        if ((*ct)->ext->n_heavy_bonded == 3) ++O_3s_count; // NR3
                                        else ++OH_count; // NHR oder NH2
                                    }
                                }
                            }
                            if (check_SO(*bt,O_count,OH_count,O_3s_count,N_count)) {
                                breaker = true;
                                break;
                            }
                        }
                        
                        O_count = 0;
                        OH_count = 0;
                        int O_2heavy = 0;
                        if ((*bt)->element == "P") {
                            for (atoms_vec ct=(*bt)->bonded_atoms.begin(); ct!=(*bt)->bonded_atoms.end(); ++ct) {
                                if ((*ct)->element == "O" || (*ct)->element == "S") {
                                    ++O_count;
                                    if ((*ct)->bonded_atoms.size() == 2) {
                                        if ((*ct)->ext->n_heavy_bonded == 1) ++OH_count;
                                        else ++O_2heavy;
                                    }
                                }
                            }
                            if (check_PO(*bt,O_count,OH_count,O_2heavy)) {
                                breaker = true;
                                break;
                            }
                        }
                    }
                    if (breaker) continue;
                } else { //! > 1 heavy atom
                    breaker = false;
                    for (atoms_vec bt=(*at)->bonded_atoms.begin(); bt!=(*at)->bonded_atoms.end(); ++bt) {
                        if ((*bt)->element == "C") {
                            if ((*bt)->intern_type == "C.es") {
                                (*at)->intern_type = "O.3es"; //!Ester oder Anhydrid
                                breaker = true;
                                break;
                            }
                            if ((*bt)->intern_type == "C.et") {
                                bool ar_flag = false;
                                for (atoms_vec ct=(*at)->bonded_atoms.begin(); 
                                     ct!=(*at)->bonded_atoms.end(); ++ct) {
                                    if ((*ct)->ext->is_aromatic) {
                                        ar_flag = true;
                                        break;
                                    }
                                }
                                if (ar_flag) (*at)->intern_type = "O.3eta";
                                else (*at)->intern_type = "O.3et";
                                //!Wenn sich spaeter herausstellt, dass an S.o3 oder S.o4 oder
                                //!P.o3 oder P.o4 gebunden, dann
                                //!auf O.3po bzw. O.3so setzen!!!
                                breaker = true;
                                continue;
                            }
                        }
                    }
                    if ((*at)->intern_type == "X") (*at)->intern_type = "O.3et";
                    if ((*at)->intern_type == "O.3et" || (*at)->intern_type == "O.3eta" || (*at)->intern_type == "O.3es") {
                        for (atoms_vec bt=(*at)->bonded_atoms.begin(); bt!=(*at)->bonded_atoms.end(); ++bt) {
                            if ((*bt)->element == "P") {
                                (*at)->intern_type = "O.3po";
                                breaker = true;
                                break;
                            } else if ((*bt)->element == "S") {
                                (*at)->intern_type = "O.3so";
                                breaker = true;
                                break;
                            }
                        }
                    }
                    if (breaker) continue;
                }
            }
        }
        //!****************************************************************************************
    }
    //!Jetzt noch die Schwefel- und Phosphor-Typen:
    for (atoms_vec at=atoms.begin(); at!=atoms.end(); ++at) {
        if ((*at)->intern_type != "X") continue;
        
        //!****************************************************************************************
        //!Schwefel  3. Durchlauf:
        //!--------------------------
        //!folgende Typen sind hier moeglich:
        //! S.thi    sp2 in Thionylgruppe
        //! S.o        in SO            bereits gesetzt
        //! S.o2h    in SO2NH        bereits gesetzt
        //! S.o2    in SO2            bereits gesetzt
        //! S.o3h    in SO3h            bereits gesetzt
        //! S.o3    in SO3            bereits gesetzt
        //! S.o4h    in OSO3h        bereits gesetzt
        //! S.o4    in OSO3            bereits gesetzt
        //! S.2        sonstige sp2 S
        //! S.sh    in SH-Gruppen
        //! S.s        in Disulfidbruecken
        //! S.3        sonstige sp3 S
        if ((*at)->element == "S") {
            if ((*at)->ext->n_heavy_bonded == 1) {
                if ((*at)->bonded_atoms[0]->element == "C") {
                    (*at)->intern_type = "S.sh"; //!SH-Gruppe
                    (*at)->ext->hybridization = 3;
                } else {
                    (*at)->intern_type = "S.3"; //! sp3  S- 
                    (*at)->ext->hybridization = 3;
                }
                continue;
            }
            int O_count = 0;
            int N_count = 0;
            int O_neg = 0;
            bool S_flag = false;
            bool thi_flag = false;
            for (atoms_vec bt=(*at)->bonded_atoms.begin(); bt!=(*at)->bonded_atoms.end(); ++bt) {
                if ((*bt)->element == "O") {
                    ++O_count;
                    if ((*bt)->ext->n_heavy_bonded == 1 && (*bt)->bonded_atoms.size() == 1) ++O_neg;
                } else if ((*bt)->element == "N") {
                    if ((*bt)->ext->n_heavy_bonded == 1) ++N_count;
                } else if ((*bt)->element == "S") S_flag = true;
                else if ((*bt)->intern_type == "C.s") thi_flag = true;
            }
            if (thi_flag) {
                if (O_neg > 0) {
                    (*at)->intern_type = "S.2"; //! sp2 S-
                    (*at)->ext->hybridization = 2;
                } else {
                    (*at)->intern_type = "S.thi"; //!Thionylgruppe
                    (*at)->ext->hybridization = 2;
                }
                continue;
            }
            
            //!eigentlich ueberfluessiger Teil:
            if (O_count == 3) {
                (*at)->intern_type = "S.o3";
                for (atoms_vec bt=(*at)->bonded_atoms.begin(); bt!=(*at)->bonded_atoms.end(); ++bt) {
                    if ((*bt)->intern_type == "O.3et" || (*bt)->intern_type == "O.3eta") {
                        (*bt)->intern_type = "O.3es";
                        break;
                    }
                }
                continue;
            } else if (O_count == 2) {
                (*at)->intern_type = "S.o2";
                continue;
            } else if (O_count == 1) {
                (*at)->intern_type = "S.o";
                continue;
            } else if (N_count == 2) {
                (*at)->intern_type = "S.2";
                (*at)->ext->hybridization = 2;
                for (atoms_vec bt=(*at)->bonded_atoms.begin(); bt!=(*at)->bonded_atoms.end(); ++bt) {
                    if ((*bt)->element == "N") {
                        (*bt)->intern_type = "N.2p";
                        (*bt)->ext->hybridization = 2;
                    }
                }
                continue;
            }
            //!=================================
            
            if (S_flag) {
                (*at)->intern_type = "S.s"; //!Disulfid
                (*at)->ext->hybridization = 3;
                continue;
            }

            if ((*at)->ext->hybridization == 2) (*at)->intern_type = "S.2";
            else (*at)->intern_type = "S.3"; //!Der Rest
            continue;
        }
        //!****************************************************************************************
        
        //!****************************************************************************************
        //!Phosphor  3. Durchlauf:
        //!--------------------------
        //!folgende Typen sind hier moeglich: (alle bis auf P.3 bereits gesetzt!!!)
        //! P.o        in PO
        //! P.o2h    in protonierter PO2
        //! P.o3h    in protonierter PO3
        //! P.o4h    in protonierter OPO3
        //! P.o2    in deprotonierter PO2 bzw. unklarer Protonierungsstatus
        //! P.o3    in deprotonierter PO3  bzw. unklarer Protonierungsstatus
        //! P.o4    in deprotonierter OPO3  bzw. unklarer Protonierungsstatus
        //! P.3        sonstige sp3 P
        if ((*at)->element == "P") {
            (*at)->intern_type = "P.3"; //! Alle anderen wurden schon gesetzt
            (*at)->ext->hybridization = 3;
            continue;
        }
    }
    //!************************************************************************************************

    //! Jetzt die Halogene:****************************************************************************
    for (atoms_vec at=atoms.begin(); at!=atoms.end(); ++at) {
        if ((*at)->element == "F") {
            if ((*at)->bonded_atoms.size() > 0) (*at)->intern_type = "F.0";
            else (*at)->intern_type = "F.i";
        } else if ((*at)->element == "Cl") {
            if ((*at)->bonded_atoms.size() > 0) (*at)->intern_type = "Cl.0";
            else (*at)->intern_type = "Cl.i";
        } else if ((*at)->element == "Br") {
            if ((*at)->bonded_atoms.size() > 0) (*at)->intern_type = "Br.0";
            else (*at)->intern_type = "Br.i";
        } else if ((*at)->element == "I") {
            if ((*at)->bonded_atoms.size() > 0) (*at)->intern_type = "I.0";
            else (*at)->intern_type = "I.i";
        }
    }
    //!**********************************************************************************************
    
    
    for (atoms_vec at=atoms.begin(); at!=atoms.end(); ++at) {
        if ((*at)->element != "H") continue;
        
        //!****************************************************************************************
        //!Wasserstoff  4. Durchlauf:
        //!--------------------------
        //!folgende Typen sind hier moeglich:
        //! H.ac    acider Wasserstoff (an O.3ac gebunden oder N.im oder N.sam oder N.ohac)
        //! H.onh    Amid  H
        //! H.n        an sonstigem Stickstoff
        //! H.o        an sonstigem Sauerstoff
        //! H.0        sonstige
        if ((*at)->bonded_atoms.size() == 0) {
            (*at)->intern_type = "H.0";
            continue;
        }
        if ((*at)->bonded_atoms[0]->intern_type == "O.3ac") {
            (*at)->intern_type = "H.ac";
            continue;
        }
        if ((*at)->bonded_atoms[0]->element == "N") {
            if ((*at)->bonded_atoms[0]->intern_type == "N.ims" ||
                (*at)->bonded_atoms[0]->intern_type == "N.imt" ||
                (*at)->bonded_atoms[0]->intern_type == "N.samp" ||
                (*at)->bonded_atoms[0]->intern_type == "N.sams" ||
                (*at)->bonded_atoms[0]->intern_type == "N.samt" ||
                (*at)->bonded_atoms[0]->intern_type == "N.ohac") {
                (*at)->intern_type = "H.ac";
                continue;
            }
            if ((*at)->bonded_atoms[0]->intern_type == "N.amp" ||
                (*at)->bonded_atoms[0]->intern_type == "N.ams" ||
                (*at)->bonded_atoms[0]->intern_type == "N.amt") {
                (*at)->intern_type = "H.onh";
                continue;
            }
            (*at)->intern_type = "H.n";
            continue;
        }
        if ((*at)->bonded_atoms[0]->element == "O") {
            (*at)->intern_type = "H.o";
            continue;
        }
        (*at)->intern_type = "H.0";
        continue;
    }
    //!************************************************************************************************
    
    //!NEU: 03.04.2007:  Fuer C.2p und C.2s nochmal pruefen, ob diese wirklich an ein C.2 oder N.2 gebunden sind
    for (atoms_vec at=atoms.begin(); at!=atoms.end(); ++at) {
        if ((*at)->element == "C") {
        if ((*at)->intern_type == "C.guh") {
            for (atoms_vec bt=(*at)->bonded_atoms.begin(); bt!=(*at)->bonded_atoms.end(); ++bt) {
                if ((*bt)->element == "N") {
                    if ((*bt)->intern_type != "N.guh" && (!(*bt)->ext->is_aromatic)) {
                        if ((*bt)->intern_type != "N.ams" && (*bt)->intern_type != "N.amt" &&
                            (*bt)->intern_type != "N.sams" && (*bt)->intern_type != "N.samt") {
                            (*bt)->intern_type = "N.guh";
                            (*bt)->ext->hybridization = 2;
                        }
                    }
                }
            }
        } else if ((*at)->intern_type == "C.mih") {
            for (atoms_vec bt=(*at)->bonded_atoms.begin(); bt!=(*at)->bonded_atoms.end(); ++bt) {
                if ((*bt)->element == "N") {
                    if ((*bt)->intern_type != "N.mih" && (!(*bt)->ext->is_aromatic)) {
                        if ((*bt)->intern_type != "N.ams" && (*bt)->intern_type != "N.amt" &&
                            (*bt)->intern_type != "N.sams" && (*bt)->intern_type != "N.samt") {
                            (*bt)->intern_type = "N.mih";
                            (*bt)->ext->hybridization = 2;
                        }
                    }
                }
            }
        }
        } else if ((*at)->element == "N") {
        if ((*at)->intern_type == "N.amp" || (*at)->intern_type == "N.ams" || (*at)->intern_type == "N.amt") {
            for (atoms_vec bt=(*at)->bonded_atoms.begin(); bt!=(*at)->bonded_atoms.end(); ++bt) {
                if ((*bt)->element == "C" && (*bt)->ext->hybridization == 3) {
                    if (!check_trigo(*bt,1.2)) continue;
                    for (atoms_vec ct=(*bt)->bonded_atoms.begin(); ct!=(*bt)->bonded_atoms.end(); ++ct) {
                        if ((*ct)->element == "O" && (*ct)->bonded_atoms.size() < 2) {
                            (*ct)->intern_type = "O.am";
                            (*ct)->ext->hybridization = 2;
                            (*bt)->intern_type = "C.am";
                            (*bt)->ext->hybridization = 2;
                            break;
                        }
                    }
                }
            }
        }
        }
    }


    for (atoms_vec at=atoms.begin(); at!=atoms.end(); ++at) {
        if ((*at)->intern_type == "X") {
            if ((*at)->element == "C") (*at)->intern_type = "C.3p";
            else if ((*at)->element == "N") (*at)->intern_type = "N.3p";
            else if ((*at)->element == "O") (*at)->intern_type = "O.3oh";
            else if ((*at)->element == "S") (*at)->intern_type = "S.3";
            else if ((*at)->element == "P") (*at)->intern_type = "P.3";
            else if ((*at)->element == "B") (*at)->intern_type = "B";
            else (*at)->intern_type = (*at)->element;
        }
    }
    
    
    //! Jetzt nochmal alle Hybridisierungen aktualisieren. Dies ist wichtig fuer das
    //! Setzen der Bindungen und anderer Routinen, die auf den Hybridisierungen aufbauen.
    tr1::unordered_map<string,int> a2h;
    for (int i=0; i<n_intern_types; ++i) a2h[i_t[i]] = atom_hybridizations[i];
    for (atoms_vec at=atoms.begin(); at!=atoms.end(); ++at) {
        if ((*at)->element == "C" || (*at)->element == "N") {
            if (a2h.find((*at)->intern_type) != a2h.end()) {
                (*at)->ext->hybridization = a2h[(*at)->intern_type];
            }
        }
    }
    
    //! Checken, ob 2 O.3oh an einem C:  wenn ja, dann auf O.co2 oder O.co2h:
    for (atoms_vec at=atoms.begin(); at!=atoms.end(); ++at) {
        if ((*at)->intern_type == "O.3oh") {
            if ((*at)->bonded_atoms[0]->element == "C") {
                for (atoms_vec bt=(*at)->bonded_atoms[0]->bonded_atoms.begin();
                               bt!=(*at)->bonded_atoms[0]->bonded_atoms.end(); ++bt) {
                    if ((*at)->intern_id == (*bt)->intern_id) continue;
                    if ((*bt)->intern_type == "O.3oh") {
                        if ((*at)->ext->n_heavy_bonded < int((*at)->bonded_atoms.size())) {
                            if ((*bt)->ext->n_heavy_bonded == int((*bt)->bonded_atoms.size())) {
                                (*at)->intern_type = "O.3ac";
                                (*bt)->intern_type = "O.2co2"; (*bt)->ext->hybridization = 2;
                                (*at)->bonded_atoms[0]->intern_type = "C.co2h";
                                (*at)->bonded_atoms[0]->ext->hybridization = 2;
                            }
                        } else if ((*bt)->ext->n_heavy_bonded < int((*bt)->bonded_atoms.size())) {
                            if ((*at)->ext->n_heavy_bonded == int((*at)->bonded_atoms.size())) {
                                (*bt)->intern_type = "O.3ac";
                                (*at)->intern_type = "O.2co2"; (*at)->ext->hybridization = 2;
                                (*at)->bonded_atoms[0]->intern_type = "C.co2h";
                                (*at)->bonded_atoms[0]->ext->hybridization = 2;
                            }
                        } else {
                            (*at)->intern_type = "O.co2"; (*at)->ext->hybridization = 2;
                            (*bt)->intern_type = "O.co2"; (*bt)->ext->hybridization = 2;
                            (*at)->bonded_atoms[0]->intern_type = "C.co2";
                            (*at)->bonded_atoms[0]->ext->hybridization = 2;
                        }
                        break;
                    }
                }
            } else if ((*at)->bonded_atoms.size() > 1) {
                for (atoms_vec bt=(*at)->bonded_atoms[1]->bonded_atoms.begin();
                               bt!=(*at)->bonded_atoms[1]->bonded_atoms.end(); ++bt) {
                    if ((*at)->intern_id == (*bt)->intern_id) continue;
                    if ((*bt)->intern_type == "O.3oh") {
                        if ((*at)->ext->n_heavy_bonded < int((*at)->bonded_atoms.size())) {
                            if ((*bt)->ext->n_heavy_bonded == int((*bt)->bonded_atoms.size())) {
                                (*at)->intern_type = "O.3ac";
                                (*bt)->intern_type = "O.2co2"; (*bt)->ext->hybridization = 2;
                                (*at)->bonded_atoms[1]->intern_type = "C.co2h";
                                (*at)->bonded_atoms[1]->ext->hybridization = 2;
                            }
                        } else if ((*bt)->ext->n_heavy_bonded < int((*bt)->bonded_atoms.size())) {
                            if ((*at)->ext->n_heavy_bonded == int((*at)->bonded_atoms.size())) {
                                (*bt)->intern_type = "O.3ac";
                                (*at)->intern_type = "O.2co2"; (*at)->ext->hybridization = 2;
                                (*at)->bonded_atoms[1]->intern_type = "C.co2h";
                                (*at)->bonded_atoms[1]->ext->hybridization = 2;
                            }
                        } else {
                            (*at)->intern_type = "O.co2"; (*at)->ext->hybridization = 2;
                            (*bt)->intern_type = "O.co2"; (*bt)->ext->hybridization = 2;
                            (*at)->bonded_atoms[1]->intern_type = "C.co2";
                            (*at)->bonded_atoms[1]->ext->hybridization = 2;
                        }
                        break;
                    }
                }
            }
        }
    }
    
    
    //! gegebenenfalls Protonierungszustaende aendern: ==========================================
    for (atoms_vec at=atoms.begin(); at!=atoms.end(); ++at) {
        if (prot_acids && (*at)->intern_type == "C.co2") {
            (*at)->intern_type = "C.co2h";
            //! O mit kuerzerer Bindung auf O.2co2 und das andere auf O.3ac:
            stl_ptr<ATOM> o1 = 0;
            stl_ptr<ATOM> o2 = 0;
            for (atoms_vec bt=(*at)->bonded_atoms.begin(); bt!=(*at)->bonded_atoms.end(); ++bt) {
                if ((*bt)->element == "O") {
                    if (o1.zero()) o1 = *bt;
                    else {
                        o2 = *bt;
                        break;
                    }
                }
            }
            if (get_square_distance((*at)->coord,o1->coord) < get_square_distance((*at)->coord,o2->coord)) {
                o1->intern_type = "O.2co2"; o1->ext->hybridization = 2;
                o2->intern_type = "O.3ac"; o2->ext->hybridization = 3;
            } else {
                o2->intern_type = "O.2co2"; o2->ext->hybridization = 2;
                o1->intern_type = "O.3ac"; o1->ext->hybridization = 3;
            }
        }
        
        if ((!prot_guanidin) && (*at)->intern_type == "C.guh") {
            (*at)->intern_type = "C.gu";
            //! N mit kuerzester Bindung auf N.gu1, WENN NICHT 2 H gebunden:
            stl_ptr<ATOM> n1 = 0; float d1;
            stl_ptr<ATOM> n2 = 0; float d2;
            stl_ptr<ATOM> n3 = 0; float d3;
            for (atoms_vec bt=(*at)->bonded_atoms.begin(); bt!=(*at)->bonded_atoms.end(); ++bt) {
                if ((*bt)->element == "N") {
                    if (n1.zero()) n1 = *bt;
                    else if (n2.zero()) n2 = *bt;
                    else {
                        n3 = *bt;
                        break;
                    }
                }
            }
            d1 = get_square_distance((*at)->coord,n1->coord);
            d2 = get_square_distance((*at)->coord,n2->coord);
            d3 = get_square_distance((*at)->coord,n3->coord);
            if (n1->bonded_atoms.size() < 3) {
                if (n2->bonded_atoms.size() < 3) {
                    if (n3->bonded_atoms.size() < 3) {
                        if (d1 < d2 && d1 < d3) {
                            n1->intern_type = "N.gu1"; n1->ext->hybridization = 2;
                            n2->intern_type = "N.gu2"; n1->ext->hybridization = 3;
                            n3->intern_type = "N.gu2"; n1->ext->hybridization = 3;
                        } else if (d2 < d1 && d2 < d3) {
                            n1->intern_type = "N.gu2"; n1->ext->hybridization = 3;
                            n2->intern_type = "N.gu1"; n1->ext->hybridization = 2;
                            n3->intern_type = "N.gu2"; n1->ext->hybridization = 3;
                        } else {
                            n1->intern_type = "N.gu2"; n1->ext->hybridization = 3;
                            n2->intern_type = "N.gu2"; n1->ext->hybridization = 3;
                            n3->intern_type = "N.gu1"; n1->ext->hybridization = 2;
                        }
                    } else {
                        if (d1 < d2) {
                            n1->intern_type = "N.gu1"; n1->ext->hybridization = 2;
                            n2->intern_type = "N.gu2"; n1->ext->hybridization = 3;
                            n3->intern_type = "N.gu2"; n1->ext->hybridization = 3;
                        } else {
                            n1->intern_type = "N.gu2"; n1->ext->hybridization = 3;
                            n2->intern_type = "N.gu1"; n1->ext->hybridization = 2;
                            n3->intern_type = "N.gu2"; n1->ext->hybridization = 3;
                        }
                    }
                } else if (n3->bonded_atoms.size() < 3) {
                    if (d1 < d3) {
                        n1->intern_type = "N.gu1"; n1->ext->hybridization = 2;
                        n2->intern_type = "N.gu2"; n1->ext->hybridization = 3;
                        n3->intern_type = "N.gu2"; n1->ext->hybridization = 3;
                    } else {
                        n1->intern_type = "N.gu2"; n1->ext->hybridization = 3;
                        n2->intern_type = "N.gu2"; n1->ext->hybridization = 3;
                        n3->intern_type = "N.gu1"; n1->ext->hybridization = 2;
                    }
                } else {
                    n1->intern_type = "N.gu1"; n1->ext->hybridization = 2;
                    n2->intern_type = "N.gu2"; n1->ext->hybridization = 3;
                    n3->intern_type = "N.gu2"; n1->ext->hybridization = 3;
                }
            } else if (n2->bonded_atoms.size() < 3) {
                if (n3->bonded_atoms.size() < 3) {
                    if (d2 < d3) {
                        n1->intern_type = "N.gu2"; n1->ext->hybridization = 3;
                        n2->intern_type = "N.gu1"; n1->ext->hybridization = 2;
                        n3->intern_type = "N.gu2"; n1->ext->hybridization = 3;
                    } else {
                        n1->intern_type = "N.gu2"; n1->ext->hybridization = 3;
                        n2->intern_type = "N.gu2"; n1->ext->hybridization = 3;
                        n3->intern_type = "N.gu1"; n1->ext->hybridization = 2;
                    }
                } else {
                    n1->intern_type = "N.gu2"; n1->ext->hybridization = 3;
                    n2->intern_type = "N.gu1"; n1->ext->hybridization = 2;
                    n3->intern_type = "N.gu2"; n1->ext->hybridization = 3;
                }
            } else if (n3->bonded_atoms.size() < 3) {
                n1->intern_type = "N.gu2"; n1->ext->hybridization = 3;
                n2->intern_type = "N.gu2"; n1->ext->hybridization = 3;
                n3->intern_type = "N.gu1"; n1->ext->hybridization = 2;
            }
        }
        
        if ((!prot_amidin) && (*at)->intern_type == "C.mih") {
            (*at)->intern_type = "C.mi";
            //! N mit kuerzester Bindung auf N.mi1, WENN NICHT 2 H gebunden:
            stl_ptr<ATOM> n1 = 0; float d1;
            stl_ptr<ATOM> n2 = 0; float d2;
            for (atoms_vec bt=(*at)->bonded_atoms.begin(); bt!=(*at)->bonded_atoms.end(); ++bt) {
                if ((*bt)->element == "N") {
                    if (n1.zero()) n1 = *bt;
                    else {
                        n2 = *bt;
                        break;
                    }
                }
            }
            d1 = get_square_distance((*at)->coord,n1->coord);
            d2 = get_square_distance((*at)->coord,n2->coord);
            if (n1->bonded_atoms.size() < 3) {
                if (n2->bonded_atoms.size() < 3) {
                    if (d1 < d2) {
                        n1->intern_type = "N.mi1"; n1->ext->hybridization = 2;
                        n2->intern_type = "N.mi2"; n1->ext->hybridization = 3;
                    } else {
                        n1->intern_type = "N.mi2"; n1->ext->hybridization = 3;
                        n2->intern_type = "N.mi1"; n1->ext->hybridization = 2;
                    }
                } else {
                    n1->intern_type = "N.mi1"; n1->ext->hybridization = 2;
                    n2->intern_type = "N.mi2"; n1->ext->hybridization = 3;
                }
            } else if (n2->bonded_atoms.size() < 3) {
                n1->intern_type = "N.mi2"; n1->ext->hybridization = 3;
                n2->intern_type = "N.mi1"; n1->ext->hybridization = 2;
            }
        }
        
        if (prot_amin && ((*at)->intern_type == "N.3p" || (*at)->intern_type == "N.3s" || (*at)->intern_type == "N.3t")) {
            (*at)->intern_type = "N.4h";
        }
        
        if (prot_phosphate) {
            if ((*at)->intern_type == "P.o2") {
                stl_ptr<ATOM> o1 = 0; float d1;
                stl_ptr<ATOM> o2 = 0; float d2;
                for (atoms_vec bt=(*at)->bonded_atoms.begin(); bt!=(*at)->bonded_atoms.end(); ++bt) {
                    if ((*bt)->intern_type == "O.2p") {
                        if (o1.zero()) o1 = *bt;
                        else {
                            o2 = *bt;
                            break;
                        }
                    } else if ((*bt)->intern_type == "O.2po") break;
                }
                if (o2.zero()) {
                    if (!o1.zero()) o1->intern_type = "O.2po";
                } else {
                    d1 = get_square_distance((*at)->coord,o1->coord);
                    d2 = get_square_distance((*at)->coord,o2->coord);
                    if (d1 < d2) {
                        o2->intern_type = "O.3ac"; o2->ext->hybridization = 3;
                        o1->intern_type = "O.2po";
                        (*at)->intern_type = "P.o2h";
                    } else {
                        o1->intern_type = "O.3ac"; o1->ext->hybridization = 3;
                        o2->intern_type = "O.2po";
                        (*at)->intern_type = "P.o2h";
                    }
                }
            } else if ((*at)->intern_type == "P.o3") {
                stl_ptr<ATOM> o1 = 0; float d1;
                stl_ptr<ATOM> o2 = 0; float d2;
                stl_ptr<ATOM> o3 = 0; float d3;
                for (atoms_vec bt=(*at)->bonded_atoms.begin(); bt!=(*at)->bonded_atoms.end(); ++bt) {
                    if ((*bt)->intern_type == "O.2p") {
                        if (o1.zero()) o1 = *bt;
                        else if (o2.zero()) o2 = *bt;
                        else {
                            o3 = *bt;
                            break;
                        }
                    } else if ((*bt)->intern_type == "O.2po") break;
                }
                if (o3.zero()) {
                    if (o2.zero()) {
                        if (!o1.zero()) o1->intern_type = "O.2po";
                    } else {
                        d1 = get_square_distance((*at)->coord,o1->coord);
                        d2 = get_square_distance((*at)->coord,o2->coord);
                        if (d1 < d2) {
                            o2->intern_type = "O.3ac"; o2->ext->hybridization = 3;
                            o1->intern_type = "O.2po";
                        } else {
                            o1->intern_type = "O.3ac"; o1->ext->hybridization = 3;
                            o2->intern_type = "O.2po";
                        }
                        (*at)->intern_type = "P.o3h";
                    }
                } else {
                    d1 = get_square_distance((*at)->coord,o1->coord);
                    d2 = get_square_distance((*at)->coord,o2->coord);
                    d3 = get_square_distance((*at)->coord,o3->coord);
                    if (d1 < d2 && d1 < d3) {
                        o3->intern_type = "O.3ac"; o3->ext->hybridization = 3;
                        o2->intern_type = "O.3ac"; o2->ext->hybridization = 3;
                        o1->intern_type = "O.2po";
                    } else if (d2 < d1 && d2 < d3) {
                        o3->intern_type = "O.3ac"; o3->ext->hybridization = 3;
                        o1->intern_type = "O.3ac"; o1->ext->hybridization = 3;
                        o2->intern_type = "O.2po";
                    } else {
                        o2->intern_type = "O.3ac"; o2->ext->hybridization = 3;
                        o1->intern_type = "O.3ac"; o1->ext->hybridization = 3;
                        o3->intern_type = "O.2po";
                    }
                    (*at)->intern_type = "P.o3h";
                }
            } else if ((*at)->intern_type == "P.o4") {
                stl_ptr<ATOM> o1 = 0; float d1;
                stl_ptr<ATOM> o2 = 0; float d2;
                stl_ptr<ATOM> o3 = 0; float d3;
                stl_ptr<ATOM> o4 = 0; float d4;
                for (atoms_vec bt=(*at)->bonded_atoms.begin(); bt!=(*at)->bonded_atoms.end(); ++bt) {
                    if ((*bt)->intern_type == "O.2p") {
                        if (o1.zero()) o1 = *bt;
                        else if (o2.zero()) o2 = *bt;
                        else if (o3.zero()) o3 = *bt;
                        else {
                            o4 = *bt;
                            break;
                        }
                    } else if ((*bt)->intern_type == "O.2po") break;
                }
                if (o4.zero()) {
                    if (o3.zero()) {
                        if (o2.zero()) {
                            if (!o1.zero()) o1->intern_type = "O.2po";
                        } else {
                            d1 = get_square_distance((*at)->coord,o1->coord);
                            d2 = get_square_distance((*at)->coord,o2->coord);
                            if (d1 < d2) {
                                o2->intern_type = "O.3ac"; o2->ext->hybridization = 3;
                                o1->intern_type = "O.2po";
                            } else {
                                o1->intern_type = "O.3ac"; o1->ext->hybridization = 3;
                                o2->intern_type = "O.2po";
                            }
                            (*at)->intern_type = "P.o4h";
                        }
                    } else {
                        d1 = get_square_distance((*at)->coord,o1->coord);
                        d2 = get_square_distance((*at)->coord,o2->coord);
                        d3 = get_square_distance((*at)->coord,o3->coord);
                        if (d1 < d2 && d1 < d3) {
                            o3->intern_type = "O.3ac"; o3->ext->hybridization = 3;
                            o2->intern_type = "O.3ac"; o2->ext->hybridization = 3;
                            o1->intern_type = "O.2po";
                        } else if (d2 < d1 && d2 < d3) {
                            o3->intern_type = "O.3ac"; o3->ext->hybridization = 3;
                            o1->intern_type = "O.3ac"; o1->ext->hybridization = 3;
                            o2->intern_type = "O.2po";
                        } else {
                            o2->intern_type = "O.3ac"; o2->ext->hybridization = 3;
                            o1->intern_type = "O.3ac"; o1->ext->hybridization = 3;
                            o3->intern_type = "O.2po";
                        }
                        (*at)->intern_type = "P.o4h";
                    }
                } else {
                    d1 = get_square_distance((*at)->coord,o1->coord);
                    d2 = get_square_distance((*at)->coord,o2->coord);
                    d3 = get_square_distance((*at)->coord,o3->coord);
                    d4 = get_square_distance((*at)->coord,o4->coord);
                    if (d1 < d2 && d1 < d3 && d1 < d4) {
                        o4->intern_type = "O.3ac"; o4->ext->hybridization = 3;
                        o3->intern_type = "O.3ac"; o3->ext->hybridization = 3;
                        o2->intern_type = "O.3ac"; o2->ext->hybridization = 3;
                        o1->intern_type = "O.2po";
                    } else if (d2 < d1 && d2 < d3 && d2 < d4) {
                        o4->intern_type = "O.3ac"; o4->ext->hybridization = 3;
                        o3->intern_type = "O.3ac"; o3->ext->hybridization = 3;
                        o1->intern_type = "O.3ac"; o1->ext->hybridization = 3;
                        o2->intern_type = "O.2po";
                    } else if (d3 < d1 && d3 < d2 && d3 < d4) {
                        o4->intern_type = "O.3ac"; o4->ext->hybridization = 3;
                        o1->intern_type = "O.3ac"; o1->ext->hybridization = 3;
                        o2->intern_type = "O.3ac"; o2->ext->hybridization = 3;
                        o3->intern_type = "O.2po";
                    } else {
                        o3->intern_type = "O.3ac"; o3->ext->hybridization = 3;
                        o1->intern_type = "O.3ac"; o1->ext->hybridization = 3;
                        o2->intern_type = "O.3ac"; o2->ext->hybridization = 3;
                        o4->intern_type = "O.2po";
                    }
                    (*at)->intern_type = "P.o4h";
                }
            }
        }
        
        if (prot_sulfate) {
            if ((*at)->intern_type == "S.o2") {
                stl_ptr<ATOM> o1 = 0;
                stl_ptr<ATOM> o2 = 0;
                stl_ptr<ATOM> n1 = 0;
                stl_ptr<ATOM> n2 = 0;
                for (atoms_vec bt=(*at)->bonded_atoms.begin(); bt!=(*at)->bonded_atoms.end(); ++bt) {
                    if ((*bt)->intern_type == "N.samp" || (*bt)->intern_type == "N.sams") {
                        if (n1.zero()) n1 = *bt;
                        else n2 = *bt;
                    } else if ((*bt)->intern_type == "O.2s") {
                        if (o1.zero()) o1 = *bt;
                        else o2 = *bt;
                    } else if ((*bt)->intern_type == "N.samt" || (*bt)->intern_type == "O.2so") break;
                }
                if (o2.zero()) {
                    //! alles gut
                } else {
                    if (!n1.zero()) {
                        o1->intern_type = "O.2so";
                        o2->intern_type = "O.2so";
                        (*at)->intern_type = "S.o2h";
                    } else {
                        o1->intern_type = "O.2so";
                        o2->intern_type = "O.2so";
                    }
                }
            } else if ((*at)->intern_type == "S.o3") {
                stl_ptr<ATOM> o1 = 0; float d1;
                stl_ptr<ATOM> o2 = 0; float d2;
                stl_ptr<ATOM> o3 = 0; float d3;
                stl_ptr<ATOM> n1 = 0;
                for (atoms_vec bt=(*at)->bonded_atoms.begin(); bt!=(*at)->bonded_atoms.end(); ++bt) {
                    if ((*bt)->intern_type == "N.samp" || (*bt)->intern_type == "N.sams") {
                        if (n1.zero()) n1 = *bt;
                    } else if ((*bt)->intern_type == "O.2s") {
                        if (o1.zero()) o1 = *bt;
                        else if (o2.zero()) o2 = *bt;
                        else o3 = *bt;
                    } else if ((*bt)->intern_type == "N.samt" || (*bt)->intern_type == "O.2so") break;
                }
                if (o3.zero()) {
                    if (o2.zero()) {
                        //! alles gut
                    } else {
                        if (!n1.zero()) {
                            o1->intern_type = "O.2so";
                            o2->intern_type = "O.2so";
                            (*at)->intern_type = "S.o3h";
                        } else {
                            o1->intern_type = "O.2so";
                            o2->intern_type = "O.2so";
                        }
                    }
                } else {
                    d1 = get_square_distance((*at)->coord,o1->coord);
                    d2 = get_square_distance((*at)->coord,o2->coord);
                    d3 = get_square_distance((*at)->coord,o3->coord);
                    if (d1 > d2 && d1 > d3) {
                        o3->intern_type = "O.2so";
                        o2->intern_type = "O.2so";
                        o1->intern_type = "O.3ac"; o1->ext->hybridization = 3;
                    } else if (d2 > d1 && d2 > d3) {
                        o3->intern_type = "O.2so";
                        o1->intern_type = "O.2so";
                        o2->intern_type = "O.3ac"; o2->ext->hybridization = 3;
                    } else {
                        o2->intern_type = "O.2so";
                        o1->intern_type = "O.2so";
                        o3->intern_type = "O.3ac"; o3->ext->hybridization = 3;
                    }
                    (*at)->intern_type = "S.o3h";
                }
            } else if ((*at)->intern_type == "S.o4") {
                stl_ptr<ATOM> o1 = 0; float d1;
                stl_ptr<ATOM> o2 = 0; float d2;
                stl_ptr<ATOM> o3 = 0; float d3;
                stl_ptr<ATOM> o4 = 0; float d4;
                for (atoms_vec bt=(*at)->bonded_atoms.begin(); bt!=(*at)->bonded_atoms.end(); ++bt) {
                    if ((*bt)->intern_type == "O.2s") {
                        if (o1.zero()) o1 = *bt;
                        else if (o2.zero()) o2 = *bt;
                        else if (o3.zero()) o3 = *bt;
                        else o4 = *bt;
                    } else if ((*bt)->intern_type == "O.2so") break;
                }
                if (o4.zero()) {
                    if (o3.zero()) {
                        if (o2.zero()) {
                            //! alles gut
                        } else {
                            o1->intern_type = "O.2so";
                            o2->intern_type = "O.2so";
                        }
                    } else {
                        d1 = get_square_distance((*at)->coord,o1->coord);
                        d2 = get_square_distance((*at)->coord,o2->coord);
                        d3 = get_square_distance((*at)->coord,o3->coord);
                        if (d1 > d2 && d1 > d3) {
                            o3->intern_type = "O.2so";
                            o2->intern_type = "O.2so";
                            o1->intern_type = "O.3ac"; o1->ext->hybridization = 3;
                        } else if (d2 > d1 && d2 > d3) {
                            o3->intern_type = "O.2so";
                            o1->intern_type = "O.2so";
                            o2->intern_type = "O.3ac"; o2->ext->hybridization = 3;
                        } else {
                            o2->intern_type = "O.2so";
                            o1->intern_type = "O.2so";
                            o3->intern_type = "O.3ac"; o3->ext->hybridization = 3;
                        }
                        (*at)->intern_type = "S.o4h";
                    }
                } else {
                    d1 = get_square_distance((*at)->coord,o1->coord);
                    d2 = get_square_distance((*at)->coord,o2->coord);
                    d3 = get_square_distance((*at)->coord,o3->coord);
                    d4 = get_square_distance((*at)->coord,o4->coord);
                    float w1 = d1+d2;
                    float w2 = d1+d3;
                    float w3 = d1+d4;
                    float w4 = d2+d3;
                    float w5 = d2+d4;
                    float w6 = d3+d4;
                    if (w1<w2 && w1<w3 && w1<w4 && w1<w5 && w1<w6) {
                        o1->intern_type = "O.2so";
                        o2->intern_type = "O.2so";
                        o3->intern_type = "O.3ac"; o3->ext->hybridization = 3;
                        o4->intern_type = "O.3ac"; o4->ext->hybridization = 3;
                    } else if (w2<w1 && w2<w3 && w2<w4 && w2<w5 && w2<w6) {
                        o1->intern_type = "O.2so";
                        o2->intern_type = "O.3ac"; o2->ext->hybridization = 3;
                        o3->intern_type = "O.2so";
                        o4->intern_type = "O.3ac"; o4->ext->hybridization = 3;
                    } else if (w3<w2 && w3<w1 && w3<w4 && w3<w5 && w3<w6) {
                        o1->intern_type = "O.2so";
                        o2->intern_type = "O.3ac"; o2->ext->hybridization = 3;
                        o3->intern_type = "O.3ac"; o3->ext->hybridization = 3;
                        o4->intern_type = "O.2so";
                    } else if (w4<w2 && w4<w3 && w4<w1 && w4<w5 && w4<w6) {
                        o1->intern_type = "O.3ac"; o1->ext->hybridization = 3;
                        o2->intern_type = "O.2so";
                        o3->intern_type = "O.2so";
                        o4->intern_type = "O.3ac"; o4->ext->hybridization = 3;
                    } else if (w5<w2 && w5<w3 && w5<w4 && w5<w1 && w5<w6) {
                        o1->intern_type = "O.3ac"; o1->ext->hybridization = 3;
                        o2->intern_type = "O.2so";
                        o3->intern_type = "O.3ac"; o3->ext->hybridization = 3;
                        o4->intern_type = "O.2so";
                    } else {
                        o1->intern_type = "O.3ac"; o1->ext->hybridization = 3;
                        o2->intern_type = "O.3ac"; o2->ext->hybridization = 3;
                        o3->intern_type = "O.2so";
                        o4->intern_type = "O.2so";
                    }
                    (*at)->intern_type = "S.o4h";
                }
            }
        }
    }
}


void MOLECULE::add_bond_object(stl_ptr<ATOM> &at1, stl_ptr<ATOM> &at2,string const& tpy,bool test) {
    if (test) for (bonds_vec bv=bonds.begin(); bv!=bonds.end(); ++bv) {
        if ((*bv)->from->intern_id != at1->intern_id) continue;
        else if ((*bv)->to->intern_id != at2->intern_id) continue;
        if (tpy != "1") (*bv)->type = tpy;
        return; //diese Bindung gibt es schon
    }
    
    //!neues Bond Object erzeugen:
    stl_ptr<BOND> bnd(new BOND());
    if (at2->intern_id > at1->intern_id) {
        bnd->from = at1;
        bnd->to = at2;
    } else {
        bnd->from = at2;
        bnd->to = at1;
    }
    bnd->type = tpy;
    bonds.push_back(bnd);
    //!jeweils zu den existing_bonds zufuegen:
    at1->ext->existing_bonds.push_back(at2);
    at2->ext->existing_bonds.push_back(at1);
}


void MOLECULE::make_bonds(bool const& no_free_rot_planar,bool const& prot_acids,bool const& prot_guanidin,bool const& prot_amidin,
                          bool const& prot_amin,bool const& prot_phosphate,bool const& prot_sulfate,bool const& kekulize_aromatics,
                          bool const& kekulize_charged) {
    string curr_type;

    for (atoms_vec at=atoms.begin(); at!=atoms.end(); ++at) {
        for (atoms_vec bt=(*at)->bonded_atoms.begin(); bt!=(*at)->bonded_atoms.end(); ++bt) {
            if ((*at)->intern_id > (*bt)->intern_id) continue; //! Jede Bindung nur einmal   => test == false fuer add_bond_object
            if ((*at)->element == "H" || (*bt)->element == "H") {
                curr_type = "1";
                add_bond_object(*at,*bt,curr_type,false);
            } else if ((((*at)->element == "S" || (*at)->element == "P") && ((*bt)->element == "O" && (*bt)->ext->hybridization == 2)) ||
                       (((*bt)->element == "S" || (*bt)->element == "P") && ((*at)->element == "O" && (*at)->ext->hybridization == 2))) {
                if ((*at)->intern_type == "O.2p" || (*at)->intern_type == "O.2s" || (*bt)->intern_type == "O.2p" || (*bt)->intern_type == "O.2s") {
                    curr_type = "ar";
                    add_bond_object(*at,*bt,curr_type,false);
                } else {
                    curr_type = "2";
                    add_bond_object(*at,*bt,curr_type,true);
                }
            } else if ((*at)->ext->hybridization == 2) {
                if ((*bt)->ext->hybridization == 3) {
                    curr_type = "1";
                    add_bond_object(*at,*bt,curr_type,true);
                } else if ((*bt)->ext->hybridization == 2) {
                    //! "2" oder "am" oder "ar" oder "1"
                    if ((*at)->ext->is_aromatic && (*bt)->ext->is_aromatic) {
                        if ((*at)->is_in_same_ring(*bt)) {
                            curr_type = "ar";
                            if (kekulize_aromatics) {
                                if (((*at)->element != "C" && (*at)->element != "N") ||
                                    ((*bt)->element != "C" && (*bt)->element != "N")) {
                                    add_bond_object(*at,*bt,string("1"),true);
                                }
                            } else add_bond_object(*at,*bt,curr_type,false);
                        } else {
                            curr_type = "1";
                            add_bond_object(*at,*bt,curr_type,false);
                        }
                    } else if ((*at)->ext->is_aromatic || (*bt)->ext->is_aromatic) {
                        curr_type = "1";
                        add_bond_object(*at,*bt,curr_type,false);
                    } else if (((*at)->intern_type == "C.am" || (*at)->intern_type == "C.es" || (*at)->intern_type == "S.o2" || (*at)->intern_type == "S.o3") &&
                               (*bt)->element == "N") {
                        if (((*bt)->intern_type == "N.samp" || (*bt)->intern_type == "N.sams") && prot_sulfate == false) {
                            curr_type = "ar";
                            add_bond_object(*at,*bt,curr_type,true);
                        } else {
                            curr_type = "am";
                            add_bond_object(*at,*bt,curr_type,true);
                        }
                    } else if (((*bt)->intern_type == "C.am" || (*bt)->intern_type == "C.es" || (*bt)->intern_type == "S.o2" || (*bt)->intern_type == "S.o3") &&
                               (*at)->element == "N") {
                        if (((*at)->intern_type == "N.samp" || (*at)->intern_type == "N.sams") && prot_sulfate == false) {
                            curr_type = "ar";
                            add_bond_object(*at,*bt,curr_type,true);
                        } else {
                            curr_type = "am";
                            add_bond_object(*at,*bt,curr_type,true);
                        }
                    } else if ((*at)->intern_type == "S.o2h" && ((*bt)->intern_type == "N.samp" || (*bt)->intern_type == "N.sams" ||
                               (*bt)->intern_type == "N.samt")) {
                        curr_type = "am";
                        add_bond_object(*at,*bt,curr_type,true);
                    } else if ((*bt)->intern_type == "S.o2h" && ((*at)->intern_type == "N.samp" || (*at)->intern_type == "N.sams" ||
                               (*at)->intern_type == "N.samt")) {
                        curr_type = "am";
                        add_bond_object(*at,*bt,curr_type,true);
                    } else if ((*at)->intern_type == "O.co2" || (*at)->intern_type == "O.2p" || (*at)->intern_type == "O.2s" ||
                               (*at)->intern_type == "O.n" ||
                               (*bt)->intern_type == "O.co2" || (*bt)->intern_type == "O.2p" || (*bt)->intern_type == "O.2s" ||
                               (*bt)->intern_type == "O.n") {
                        curr_type = "ar";
                        add_bond_object(*at,*bt,curr_type,true);
                    } else if (((*at)->intern_type == "C.guh" || (*at)->intern_type == "C.mih") &&
                               (*bt)->element == "N") {
                        curr_type = "ar";
                        add_bond_object(*at,*bt,curr_type,true);
                    } else if (((*bt)->intern_type == "C.guh" || (*bt)->intern_type == "C.mih") &&
                               (*at)->element == "N") {
                        curr_type = "ar";
                        add_bond_object(*at,*bt,curr_type,true);
                    } else if ((*at)->intern_type == "O.carb" || (*at)->intern_type == "O.am" ||
                               (*at)->intern_type == "O.2po" || (*at)->intern_type == "O.2so" ||
                               (*at)->intern_type == "O.2co2" || (*at)->intern_type == "O.2es" ||
                               (*at)->intern_type == "O.hal" || (*at)->intern_type == "S.thi" ||
                               (*bt)->intern_type == "O.carb" || (*bt)->intern_type == "O.am" ||
                               (*bt)->intern_type == "O.2po" || (*bt)->intern_type == "O.2so" ||
                               (*bt)->intern_type == "O.2co2" || (*bt)->intern_type == "O.2es" ||
                               (*bt)->intern_type == "O.hal" || (*bt)->intern_type == "S.thi") {
                        curr_type = "2";
                        add_bond_object(*at,*bt,curr_type,true);
                    } else if (((*at)->intern_type == "N.ams" && ((*bt)->intern_type == "C.2p" ||
                               (*bt)->intern_type == "C.2s" || (*bt)->intern_type == "C.2t")) ||
                               ((*at)->intern_type == "N.amt" && ((*bt)->intern_type == "C.2p" ||
                               (*bt)->intern_type == "C.2s" || (*bt)->intern_type == "C.2t")) ||
                               ((*bt)->intern_type == "N.ams" && ((*at)->intern_type == "C.2p" ||
                               (*at)->intern_type == "C.2s" || (*at)->intern_type == "C.2t")) ||
                               ((*bt)->intern_type == "N.amt" && ((*at)->intern_type == "C.2p" ||
                               (*at)->intern_type == "C.2s" || (*at)->intern_type == "C.2t"))) {
                        curr_type = "1";
                        add_bond_object(*at,*bt,curr_type,true);
                    } else if (((*at)->intern_type == "N.2p" && (*bt)->intern_type == "S.2") ||
                               ((*bt)->intern_type == "N.2p" && (*at)->intern_type == "S.2")) {
                        curr_type = "2";
                        add_bond_object(*at,*bt,curr_type,true);
                    } else {
                        curr_type = "1";
                        add_bond_object(*at,*bt,curr_type,true);
                    }
                } else {
                    curr_type = "1";
                    add_bond_object(*at,*bt,curr_type,true);
                }
            } else {
                curr_type = "1";
                add_bond_object(*at,*bt,curr_type,true);
            }
        }
    }

    // Bindungen nochmal sortieren:
    sort(bonds.begin(),bonds.end(),bond_sort);

    // Drehbarkeit bestimmen:
    int bnd_id = 1;
    for (bonds_vec bit=bonds.begin(); bit!=bonds.end(); ++bit) {
        (*bit)->id = bnd_id;
        ++bnd_id;
        if ((*bit)->type == "am") { // Carbonsaeureamide nicht frei drehbar
            //!Sulfonamide sind frei drehbar (naja, zumindest eher als normale Amide)!!!
            if ((*bit)->from->intern_type != "S.o2" && (*bit)->to->intern_type != "S.o2") continue;
        } else if ((*bit)->type != "1") continue;
        if ((*bit)->from->ext->n_heavy_bonded < 2 || (*bit)->to->ext->n_heavy_bonded < 2) continue; // 08.04.09 (keine
        if (((*bit)->from->intern_type == "C.es" && (*bit)->to->intern_type == "O.3es") ||          //           Endstaendigen)
            ((*bit)->from->intern_type == "O.3es" && (*bit)->to->intern_type == "C.es")) continue; // 08.04.09 (Ester nicht drehbar)
        if ((*bit)->from->ext->is_ring && (*bit)->to->ext->is_ring) {
            if ((*bit)->from->ext->is_aromatic && (*bit)->to->ext->is_aromatic) {// 08.04.09 (konjugierte Aromaten)
                if (check_planar_tors((*bit)->from,(*bit)->to)) {
                    //(*bit)->type = "ar"; // nicht mehr machen
                    if (no_free_rot_planar) continue;
                }
            }
            bool contin = false;
            for (rings_vec rit=rings.begin(); rit!=rings.end(); ++rit) {
                bool ffrom = false;
                bool fto = false;
                for (atoms_vec at=(*rit)->ring_atoms.begin(); at!=(*rit)->ring_atoms.end(); ++at) {
                    if (*at == (*bit)->from) ffrom = true;
                    if (*at == (*bit)->to) fto = true;
                }
                if (ffrom && fto) { // beide im gleichen Ring => nicht drehbar (ring---ring ist evtl drehbar)
                    contin = true;
                    break;
                }
            }
            if (contin) continue;
        }
        if (no_free_rot_planar) {
            //! Neu: 08.04.09
            //! einfach-Bindungen in planaren, konjugierten Systemen nicht als frei drehbar markieren:
            if ((*bit)->from->ext->hybridization == 2 && (*bit)->to->ext->hybridization == 2) {
                if (check_planar_tors((*bit)->from,(*bit)->to)) continue;
            }
        }
        (*bit)->free_rot = true;
    }
}


void MOLECULE::generate_new_def_file(char const* name) {
    //! Ein veraltetes deffile einlesen und anpassen, um ein neues rauszuschreiben:
    tr1::unordered_map<string,string> old_map;
    tr1::unordered_map<string,string> new_map;
    tr1::unordered_map<string,string> old_flags;
    tr1::unordered_map<string,string> new_flags;
    //! Namen fuer das neue deffile erfragen:
    string fname = ""; string buf;
    cout << " -> Please type in a name for your new definition file: ";
    cin >> fname;
    getline(cin,buf);
    fname += buf;

    cout << endl;

    //! Das alte einlesen:
    istringstream is;
    string row,key,vers,tot;
    ifstream f_in;
    f_in.open(name);
    if (!f_in) {
        cerr << c_message<cERROR>("could not open ") << name << endl;
        exit(1);
    }
    while (!f_in.eof()) {
        getline(f_in,row);
        is.clear(); is.str(row);
        is >> key;
        if (key == "<HEADER>") break;
    }
    while (!f_in.eof()) {
        getline(f_in,row);
        is.clear(); is.str(row);
        is >> key;
        if (key == "internal_set_type") is >> vers;
        else if (key == "total") {is >> tot; break;}
    }
    int n_types = 0;
    bool in_flags = false;
    while (!f_in.eof()) {
        getline(f_in,row);
        is.clear(); is.str(row);
        is >> key;
        if (key == "*") {
            is >> vers;
            is >> tot;
            old_map[vers] = tot;
            ++n_types;
            vers = "UNK";
            tot = "UNK";
            key = "X";
        } else if (key == "<FLAGS>") in_flags = true;
        else if (in_flags) {
            if (key == "protonate_acids") {
                string ti;
                is >> ti;
                old_flags[key] = ti;
            } else if (key == "protonate_amidin") {
                string ti;
                is >> ti;
                old_flags[key] = ti;
            } else if (key == "protonate_guanidin") {
                string ti;
                is >> ti;
                old_flags[key] = ti;
            } else if (key == "protonate_amine") {
                string ti;
                is >> ti;
                old_flags[key] = ti;
            } else if (key == "protonate_phosphate") {
                string ti;
                is >> ti;
                old_flags[key] = ti;
            } else if (key == "protonate_sulfate") {
                string ti;
                is >> ti;
                old_flags[key] = ti;
            } else if (key == "set_bonds") {
                string ti;
                is >> ti;
                old_flags[key] = ti;
            } else if (key == "max_ring_members") {
                string ti;
                is >> ti;
                old_flags[key] = ti;
            }
        }
    }
    f_in.close();
    
    if (old_flags.find("protonate_acids") == old_flags.end()) {
        new_flags["protonate_acids"] = "-1";
        cout << " -> There is a new flag: " << "protonate_acids" << " -1" << endl;
    } else new_flags["protonate_acids"] = old_flags["protonate_acids"];

    if (old_flags.find("protonate_amidin") == old_flags.end()) {
        new_flags["protonate_amidin"] = "-1";
        cout << " -> There is a new flag: " << "protonate_amidin" << " -1" << endl;
    } else new_flags["protonate_amidin"] = old_flags["protonate_amidin"];

    if (old_flags.find("protonate_guanidin") == old_flags.end()) {
        new_flags["protonate_guanidin"] = "-1";
        cout << " -> There is a new flag: " << "protonate_guanidin" << " -1" << endl;
    } else new_flags["protonate_guanidin"] = old_flags["protonate_guanidin"];

    if (old_flags.find("protonate_amine") == old_flags.end()) {
        new_flags["protonate_amine"] = "-1";
        cout << " -> There is a new flag: " << "protonate_amine" << " -1" << endl;
    } else new_flags["protonate_amine"] = old_flags["protonate_amine"];

    if (old_flags.find("protonate_phosphate") == old_flags.end()) {
        new_flags["protonate_phosphate"] = "-1";
        cout << " -> There is a new flag: " << "protonate_phosphate" << " -1" << endl;
    } else new_flags["protonate_phosphate"] = old_flags["protonate_phosphate"];

    if (old_flags.find("protonate_sulfate") == old_flags.end()) {
        new_flags["protonate_sulfate"] = "-1";
        cout << " -> There is a new flag: " << "protonate_sulfate" << " -1" << endl;
    } else new_flags["protonate_sulfate"] = old_flags["protonate_sulfate"];

    if (old_flags.find("set_bonds") == old_flags.end()) {
        new_flags["set_bonds"] = "-1";
        cout << " -> There is a new flag: " << "set_bonds" << " -1" << endl;
    } else new_flags["set_bonds"] = old_flags["set_bonds"];

    if (old_flags.find("kekulize_aromatics") == old_flags.end()) {
        new_flags["kekulize_aromatics"] = "-1";
        cout << " -> There is a new flag: " << "kekulize_aromatics" << " -1" << endl;
    } else new_flags["kekulize_aromatics"] = old_flags["kekulize_aromatics"];

//    if (old_flags.find("kekulize_charged") == old_flags.end()) {
//        new_flags["kekulize_charged"] = "-1";
//        cout << " -> There is a new flag: " << "kekulize_charged" << " -1" << endl;
//    } else new_flags["kekulize_charged"] = old_flags["kekulize_charged"];

    if (old_flags.find("allow_charged_aromatics") == old_flags.end()) {
        new_flags["allow_charged_aromatics"] = "-1";
        cout << " -> There is a new flag: " << "allow_charged_aromatics" << " -1" << endl;
    } else new_flags["allow_charged_aromatics"] = old_flags["allow_charged_aromatics"];

    if (old_flags.find("max_ring_members") == old_flags.end()) {
        new_flags["max_ring_members"] = "-1";
        cout << " -> There is a new flag: " << "max_ring_members" << " -1" << endl;
    } else new_flags["max_ring_members"] = old_flags["max_ring_members"];
    
    //! Jetzt das SOLL generieren:
    for (int i=0; i<n_intern_types; ++i) new_map[i_t[i]] = "X";
    cout << endl << " -> Types that are no longer present:" << endl;
    for (tr1::unordered_map<string,string>::iterator it=old_map.begin(); it!=old_map.end(); ++it) {
        if (new_map.find(it->first) != new_map.end()) new_map[it->first] = it->second;
        else cout << "     internal type = " << it->first << "   mapped type = " << it->second << endl;
    }
    cout << endl << " -> Types that are new: Please type in the correct mapping" << endl;
    cout << "     (If you don't know, just type something. Afterwards you can watch the description "
         << "in the new file and edit it)" << endl;
    for (tr1::unordered_map<string,string>::iterator it=new_map.begin(); it!=new_map.end(); ++it) {
        if (old_map.find(it->first) == old_map.end()) {
            string mt = ""; string mbuf;
            cout << "     internal type = " << it->first << "   mapped type = ";
            cin >> mt;
            getline(cin,mbuf);
            mt += mbuf;
            it->second = mt;
        }
    }
    cout << endl;
    
    //! Jetzt nach mode_2 uebertragen und rausschreiben:
    vector<string> wmode;
    for (int i=0; i<n_intern_types; ++i) wmode.push_back(new_map[i_t[i]]);
    write_def_file(0,fname.c_str(),&wmode,new_flags);
    cout << " -> The program will exit now. Please recall it with the newly generated file." << endl;
}

void MOLECULE::load_def_file(char const* name,
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
                             bool &allow_charged_aromatics) {
    istringstream is;
    string row,key,vers,tot;
    ifstream f_in;
    f_in.open(name);
    if (!f_in) {
        cerr << c_message<cERROR>("could not open ") << name << endl;
        exit(1);
    }
    while (!f_in.eof()) {
        getline(f_in,row);
        is.clear(); is.str(row);
        is >> key;
        if (key == "<HEADER>") break;
    }
    while (!f_in.eof()) {
        getline(f_in,row);
        is.clear(); is.str(row);
        is >> key;
        if (key == "internal_set_type") is >> vers;
        else if (key == "total") {is >> tot; break;}
    }
    if (vers != a_t_version) {
        f_in.close();
        cerr << c_message<cERROR>("MOLECULE::load_def_file --> your definition file has a wrong version number") << endl;
        cerr << " -> your version:     " << vers << endl;
        cerr << " -> demanded version: " << a_t_version << endl;
        cerr << " -> A new file may be generated automatically with some user interventions." << endl;
        cerr << "    Please type 'y' if you want to generate a new file now: ";
        char dummy = 'n'; cin >> dummy;
        if (dummy == 'y' || dummy == 'Y') generate_new_def_file(name);
        exit(1);
    }
    int n_types = 0;
    bool in_flags = false;
    while (!f_in.eof()) {
        getline(f_in,row);
        is.clear(); is.str(row);
        is >> key;
        if (key == "*") {
            is >> vers;
            is >> tot;
            if (def_map::last_a_t_map == 1) def_map::a_t_map[vers] = tot;
            else def_map::a_t_map2[vers] = tot;
            ++n_types;
            vers = "UNK";
            tot = "UNK";
            key = "X";
        } else if (key == "<FLAGS>") in_flags = true;
        else if (in_flags) {
            if (key == "protonate_acids") {
                int ti = -1;
                is >> ti;
                if (ti == 0) prot_acids = false;
                else if (ti > 0) prot_acids = true;
                else continue;
            } else if (key == "protonate_amidin") {
                int ti = -1;
                is >> ti;
                if (ti == 0) prot_amidin = false;
                else if (ti > 0) prot_amidin = true;
                else continue;
            } else if (key == "protonate_guanidin") {
                int ti = -1;
                is >> ti;
                if (ti == 0) prot_guanidin = false;
                else if (ti > 0) prot_guanidin = true;
                else continue;
            } else if (key == "protonate_amine") {
                int ti = -1;
                is >> ti;
                if (ti == 0) prot_amin = false;
                else if (ti > 0) prot_amin = true;
                else continue;
            } else if (key == "protonate_phosphate") {
                int ti = -1;
                is >> ti;
                if (ti == 0) prot_phosphate = false;
                else if (ti > 0) prot_phosphate = true;
                else continue;
            } else if (key == "protonate_sulfate") {
                int ti = -1;
                is >> ti;
                if (ti == 0) prot_sulfate = false;
                else if (ti > 0) prot_sulfate = true;
                else continue;
            } else if (key == "set_bonds") {
                int ti = -1;
                is >> ti;
                if (ti == 0) get_bonds = false;
                else if (ti > 0) get_bonds = true;
                else continue;
            } else if (key == "kekulize_aromatics") {
                int ti = -1;
                is >> ti;
                if (ti == 0) kekulize_aromatics = false;
                else if (ti > 0) kekulize_aromatics = true;
                else continue;
//            } else if (key == "kekulize_charged") {
//                int ti = -1;
//                is >> ti;
//                if (ti == 0) kekulize_charged = false;
//                else if (ti > 0) kekulize_charged = true;
//                else continue;
            } else if (key == "allow_charged_aromatics") {
                int ti = -1;
                is >> ti;
                if (ti == 0) allow_charged_aromatics = false;
                else if (ti > 0) allow_charged_aromatics = true;
                else continue;
            } else if (key == "max_ring_members") {
                int ti = -1;
                is >> ti;
                if (ti < 1) continue;
                else max_ring_members = ti;
            }
        }
    }
    if (n_intern_types != n_types) {
        if (verbosity) cerr << c_message<cWARNING>("MOLECULE::load_def_file --> read ") << n_types << " atom type definitions from " << name << " , but " << n_intern_types << " were expected!" << endl;
    }
    f_in.close();
}

void MOLECULE::get_bonded_only() {
    for (atoms_vec at=atoms.begin(); at!=atoms.end(); ++at) {
        //!neues ATOM_EXT Objekt:
        if ((*at)->ext.zero()) (*at)->ext = new ATOM_EXT();
        else (*at)->ext->clear();
        //!eventuell vorhandene Bindungen loeschen:
        (*at)->bonded_atoms.clear();
    }
    get_elements();
    
    get_connections();
}

void MOLECULE::get_hybrid_only(bool prot_acids,bool prot_guanidin,bool prot_amidin,bool prot_amin,bool prot_phosphate,bool prot_sulfate,
                               bool kekulize_aromatics,bool kekulize_charged,bool allow_charged_aromatics) {
    verbosity = 0;
    for (atoms_vec at=atoms.begin(); at!=atoms.end(); ++at) {
        //!neues ATOM_EXT Objekt:
        if ((*at)->ext.zero()) (*at)->ext = new ATOM_EXT();
        else (*at)->ext->clear();
        //!eventuell vorhandene Bindungen loeschen:
        (*at)->bonded_atoms.clear();
    }
    get_elements();
    
    get_connections();

    for (rings_vec it=rings.begin(); it!=rings.end();++it) {
        it->kill();
    }
    rings.clear();

    get_hybridization(false,8,allow_charged_aromatics);

    set_intern_types(prot_acids,prot_guanidin,prot_amidin,prot_amin,prot_phosphate,prot_sulfate,  //! noetig, weil hier auch noch
                     kekulize_charged);                                                           //! Hybridisierungen umgesetzt werden
}

void MOLECULE::intern2mode(bool fill_X,const char *alt_def_file) {
    for (atoms_vec at=atoms.begin(); at!=atoms.end(); ++at) {
        if (def_map::last_a_t_map == 1) {
            if (def_map::a_t_map.find((*at)->intern_type) == def_map::a_t_map.end()) {
                if (!fill_X) (*at)->sybyl_type = (*at)->intern_type;
                else (*at)->sybyl_type = "X";
            }
            else (*at)->sybyl_type = def_map::a_t_map[(*at)->intern_type];
        } else {
            if (def_map::a_t_map2.find((*at)->intern_type) == def_map::a_t_map2.end()) {
                if (!fill_X) (*at)->sybyl_type = (*at)->intern_type;
                else (*at)->sybyl_type = "X";
            }
            else (*at)->sybyl_type = def_map::a_t_map2[(*at)->intern_type];
        }
    }
    
    //!Fuer den Fall, dass ein alternatives atom typing erfolgen soll:
    //!Die entsprechenden Typen werden dann in 'dict_type' abgelegt!!!:
    string compare2 = alt_def_file;
    if (compare2 != "X") {
        for (atoms_vec at=atoms.begin(); at!=atoms.end(); ++at) {
            if (def_map::last_a_t_map == 1) {
                if (def_map::a_t_map.find((*at)->intern_type) == def_map::a_t_map.end()) {
                    if (!fill_X) (*at)->dict_type = (*at)->intern_type;
                    else (*at)->dict_type = "X";
                }
                else (*at)->dict_type = def_map::a_t_map[(*at)->intern_type];
            } else {
                if (def_map::a_t_map2.find((*at)->intern_type) == def_map::a_t_map2.end()) {
                    if (!fill_X) (*at)->dict_type = (*at)->intern_type;
                    else (*at)->dict_type = "X";
                }
                else (*at)->dict_type = def_map::a_t_map2[(*at)->intern_type];
            }
        }
    }
}

void MOLECULE::try_for_def_file(int const& mode,const char *def_file,const char *alt_def_file,
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
                                bool &allow_charged_aromatics) {
    string compare(def_file);
    //!  a_t_map erzeugen:
    if (compare != "X") {
        if ((def_map::curr_a_t_map != compare) and (def_map::curr_a_t_map2 != compare)) {
            //!gefordertes Setting ist in keiner der beiden maps
            if (def_map::last_a_t_map == 1) { //in map 2 laden
                def_map::last_a_t_map = 2;
                load_def_file(def_file,prot_acids,prot_guanidin,prot_amidin,prot_amin,
                              prot_phosphate,prot_sulfate,get_bonds,max_ring_members,
                              kekulize_aromatics,kekulize_charged,allow_charged_aromatics);
                def_map::curr_a_t_map2 = compare;

                def_map::curr_prot_acids_default2 = prot_acids;
                def_map::curr_prot_guanidin_default2 = prot_guanidin;
                def_map::curr_prot_amidin_default2 = prot_amidin;
                def_map::curr_prot_amin_default2 = prot_amin;
                def_map::curr_prot_phosphate_default2 = prot_phosphate;
                def_map::curr_prot_sulfate_default2 = prot_sulfate;
                def_map::curr_get_bonds_default2 = get_bonds;
                def_map::curr_kekulize_aromatics_default2 = kekulize_aromatics;
                def_map::curr_kekulize_charged_default2 = kekulize_charged;
                def_map::curr_allow_charged_aromatics_default2 = allow_charged_aromatics;
                def_map::curr_max_ring_members_default2 = max_ring_members;
            } else { //in map 1 laden
                def_map::last_a_t_map = 1;
                load_def_file(def_file,prot_acids,prot_guanidin,prot_amidin,prot_amin,
                              prot_phosphate,prot_sulfate,get_bonds,max_ring_members,
                              kekulize_aromatics,kekulize_charged,allow_charged_aromatics);
                def_map::curr_a_t_map = compare;

                def_map::curr_prot_acids_default = prot_acids;
                def_map::curr_prot_guanidin_default = prot_guanidin;
                def_map::curr_prot_amidin_default = prot_amidin;
                def_map::curr_prot_amin_default = prot_amin;
                def_map::curr_prot_phosphate_default = prot_phosphate;
                def_map::curr_prot_sulfate_default = prot_sulfate;
                def_map::curr_get_bonds_default = get_bonds;
                def_map::curr_kekulize_aromatics_default = kekulize_aromatics;
                def_map::curr_kekulize_charged_default = kekulize_charged;
                def_map::curr_allow_charged_aromatics_default = allow_charged_aromatics;
                def_map::curr_max_ring_members_default = max_ring_members;
            }
        } else if (def_map::curr_a_t_map != compare) {
            //!das geforderte setting ist bereits in map 2
            def_map::last_a_t_map = 2;
            prot_acids = def_map::curr_prot_acids_default2;
            prot_guanidin = def_map::curr_prot_guanidin_default2;
            prot_amidin = def_map::curr_prot_amidin_default2;
            prot_amin = def_map::curr_prot_amin_default2;
            prot_phosphate = def_map::curr_prot_phosphate_default2;
            prot_sulfate = def_map::curr_prot_sulfate_default2;
            get_bonds = def_map::curr_get_bonds_default2;
            kekulize_aromatics = def_map::curr_kekulize_aromatics_default2;
            kekulize_aromatics = def_map::curr_kekulize_aromatics_default2;
            allow_charged_aromatics = def_map::curr_allow_charged_aromatics_default2;
            max_ring_members = def_map::curr_max_ring_members_default2;
        } else {
            //!das geforderte setting ist bereits in map 1
            def_map::last_a_t_map = 1;
            prot_acids = def_map::curr_prot_acids_default;
            prot_guanidin = def_map::curr_prot_guanidin_default;
            prot_amidin = def_map::curr_prot_amidin_default;
            prot_amin = def_map::curr_prot_amin_default;
            prot_phosphate = def_map::curr_prot_phosphate_default;
            prot_sulfate = def_map::curr_prot_sulfate_default;
            get_bonds = def_map::curr_get_bonds_default;
            kekulize_aromatics = def_map::curr_kekulize_aromatics_default;
            kekulize_aromatics = def_map::curr_kekulize_aromatics_default;
            allow_charged_aromatics = def_map::curr_allow_charged_aromatics_default;
            max_ring_members = def_map::curr_max_ring_members_default;
        }
    }
    else if (mode == 0) {
        if ((def_map::curr_a_t_map != "m0") and (def_map::curr_a_t_map2 != "m0")) {
            //!gefordertes Setting ist in keiner der beiden maps
            if (def_map::last_a_t_map == 1) { //in map 2 laden
                for (int i=0; i<n_intern_types; ++i) {
                    def_map::a_t_map2[i_t[i]] = i_t[i];
                }
                def_map::curr_a_t_map2 = "m0";
                def_map::last_a_t_map = 2;
            } else { //in map 1 laden
                for (int i=0; i<n_intern_types; ++i) {
                    def_map::a_t_map[i_t[i]] = i_t[i];
                }
                def_map::curr_a_t_map = "m0";
                def_map::last_a_t_map = 1;
            }
        } else if (def_map::curr_a_t_map != "m0") {
            //!das geforderte setting ist bereits in map 2
            def_map::last_a_t_map = 2;
        } else {
            //!das geforderte setting ist bereits in map 1
            def_map::last_a_t_map = 1;
        }
    } else if (mode == 1) {
        if ((def_map::curr_a_t_map != "m1") and (def_map::curr_a_t_map2 != "m1")) {
            //!gefordertes Setting ist in keiner der beiden maps
            if (def_map::last_a_t_map == 1) { //in map 2 laden
                for (int i=0; i<n_intern_types; ++i) {
                    def_map::a_t_map2[i_t[i]] = mode_1[i];
                }
                def_map::curr_a_t_map2 = "m1";
                def_map::last_a_t_map = 2;
            } else { //in map 1 laden
                for (int i=0; i<n_intern_types; ++i) {
                    def_map::a_t_map[i_t[i]] = mode_1[i];
                }
                def_map::curr_a_t_map = "m1";
                def_map::last_a_t_map = 1;
            }
        } else if (def_map::curr_a_t_map != "m1") {
            //!das geforderte setting ist bereits in map 2
            def_map::last_a_t_map = 2;
        } else {
            //!das geforderte setting ist bereits in map 1
            def_map::last_a_t_map = 1;
        }
    } else if (mode == 2) {
        if ((def_map::curr_a_t_map != "m2") and (def_map::curr_a_t_map2 != "m2")) {
            //!gefordertes Setting ist in keiner der beiden maps
            if (def_map::last_a_t_map == 1) { //in map 2 laden
                for (int i=0; i<n_intern_types; ++i) {
                    def_map::a_t_map2[i_t[i]] = mode_2[i];
                }
                def_map::curr_a_t_map2 = "m2";
                def_map::last_a_t_map = 2;
            } else { //in map 1 laden
                for (int i=0; i<n_intern_types; ++i) {
                    def_map::a_t_map[i_t[i]] = mode_2[i];
                }
                def_map::curr_a_t_map = "m2";
                def_map::last_a_t_map = 1;
            }
        } else if (def_map::curr_a_t_map != "m2") {
            //!das geforderte setting ist bereits in map 2
            def_map::last_a_t_map = 2;
        } else {
            //!das geforderte setting ist bereits in map 1
            def_map::last_a_t_map = 1;
        }
    }

    //!Fuer den Fall, dass ein alternatives atom typing erfolgen soll:
    //!Die entsprechenden Typen werden dann in 'dict_type' abgelegt!!!:
    string compare2 = alt_def_file;
    if (compare2 != "X") {
        if ((def_map::curr_a_t_map != compare2) and (def_map::curr_a_t_map2 != compare2)) {
            //!gefordertes Setting ist in keiner der beiden maps
            if (def_map::last_a_t_map == 1) { //in map 2 laden
                def_map::last_a_t_map = 2;
                bool dummy; int idummy;
                load_def_file(alt_def_file,dummy,dummy,dummy,dummy,
                              dummy,dummy,dummy,idummy,dummy,dummy,dummy);
                def_map::curr_a_t_map2 = compare2;
            } else { //in map 1 laden
                def_map::last_a_t_map = 1;
                bool dummy; int idummy;
                load_def_file(alt_def_file,dummy,dummy,dummy,dummy,
                              dummy,dummy,dummy,idummy,dummy,dummy,dummy);
                def_map::curr_a_t_map = compare2;
            }
        } else if (def_map::curr_a_t_map != compare2) {
            //!das geforderte setting ist bereits in map 2
            def_map::last_a_t_map = 2;
        } else {
            //!das geforderte setting ist bereits in map 1
            def_map::last_a_t_map = 1;
        }
    }
}

void MOLECULE::get_atom_typing(int mode,bool get_ele,const char *def_file,bool get_bonds,int verb,bool kill_ext,bool fill_X,
                               const char *alt_def_file,int max_ring_members,bool no_free_rot_planar,LIGAND* at_ref_mol,LIGAND* at_ref_mol2,
                               bool prot_acids,bool prot_guanidin,bool prot_amidin,bool prot_amin,bool prot_phosphate,bool prot_sulfate,
                               bool kekulize_aromatics,bool kekulize_charged,bool allow_charged_aromatics) {
    //! Komplette Uebernahme aller Eigenschaften von einer Referenz (nur Koordinaten anpassen):
    if (at_ref_mol != 0) {
        //! Referenz wird vorher selbst typisert:
        if (atoms.size() != at_ref_mol->atoms.size()) {
            cerr << c_message<cWARNING>("MOLECULE::get_atom_typing --> different number of atoms compared to the given reference")
                 << "\n" << "    => using normal atom typing without reference" << endl;
            get_atom_typing(mode,get_ele,def_file,get_bonds,verb,kill_ext,fill_X,
                            alt_def_file,max_ring_members,no_free_rot_planar,0,0,prot_acids_default,
                            prot_guanidin_default,prot_amidin_default,prot_amin_default,prot_phosphate_default,
                            prot_sulfate_default,kekulize_aromatics_default,kekulize_charged_default,allow_charged_aromatics_default);
            return;
        }
        
        at_ref_mol->get_atom_typing(mode,get_ele,def_file,get_bonds,verb,kill_ext,fill_X,
                                            alt_def_file,max_ring_members,no_free_rot_planar,0,0,prot_acids_default,
                            prot_guanidin_default,prot_amidin_default,prot_amin_default,prot_phosphate_default,
                            prot_sulfate_default,kekulize_aromatics_default,kekulize_charged_default,allow_charged_aromatics_default);
        
        if (!copy_without_coords(at_ref_mol)) {
            cerr << c_message<cWARNING>("MOLECULE::get_atom_typing --> copying reference was not possible")
                 << "\n" << "    => using normal atom typing without reference" << endl;
            get_atom_typing(mode,get_ele,def_file,get_bonds,verb,kill_ext,fill_X,
                                        alt_def_file,max_ring_members,no_free_rot_planar,0,0,prot_acids_default,
                            prot_guanidin_default,prot_amidin_default,prot_amin_default,prot_phosphate_default,
                            prot_sulfate_default,kekulize_aromatics_default,kekulize_charged_default,allow_charged_aromatics_default);
            return;
        }
    } else if (at_ref_mol2 != 0) {
        //! Referenz wird NICHT vorher typisiert:
        if (atoms.size() != at_ref_mol2->atoms.size()) {
            cerr << c_message<cWARNING>("MOLECULE::get_atom_typing --> different number of atoms compared to the given reference")
                 << "\n" << "    => using normal atom typing without reference" << endl;
            get_atom_typing(mode,get_ele,def_file,get_bonds,verb,kill_ext,fill_X,
                                        alt_def_file,max_ring_members,no_free_rot_planar,0,0,prot_acids_default,
                            prot_guanidin_default,prot_amidin_default,prot_amin_default,prot_phosphate_default,
                            prot_sulfate_default,kekulize_aromatics_default,kekulize_charged_default,allow_charged_aromatics_default);
            return;
        }
        if (!copy_without_coords(at_ref_mol2)) {
            cerr << c_message<cWARNING>("MOLECULE::get_atom_typing --> copying reference was not possible")
                 << "\n" << "    => using normal atom typing without reference" << endl;
            get_atom_typing(mode,get_ele,def_file,get_bonds,verb,kill_ext,fill_X,
                                        alt_def_file,max_ring_members,no_free_rot_planar,0,0,prot_acids_default,
                            prot_guanidin_default,prot_amidin_default,prot_amin_default,prot_phosphate_default,
                            prot_sulfate_default,kekulize_aromatics_default,kekulize_charged_default,allow_charged_aromatics_default);
            return;
        }
    } else {
        if(atom_types_already_set) return; //! neu: 23.05.2009
        
        verbosity = verb;

        try_for_def_file(mode,def_file,alt_def_file,prot_acids,prot_guanidin,prot_amidin,
                         prot_amin,prot_phosphate,prot_sulfate,get_bonds,max_ring_members,
                         kekulize_aromatics,kekulize_charged,allow_charged_aromatics);
        
        //cerr << "def_acids_default = " << prot_acids << endl;
        //cerr << "def_guanidin_default = " << prot_guanidin << endl;
        //cerr << "def_amidin_default = " << prot_amidin << endl;
        //cerr << "def_amin_default = " << prot_amin << endl;
        //cerr << "def_phosphate_default = " << prot_phosphate << endl;
        //cerr << "def_sulfate_default = " << prot_sulfate << endl;
        //cerr << "def_bonds_default = " << get_bonds << endl;
        //cerr << "def_ring_members_default = " << max_ring_members << endl << endl;

        for (atoms_vec at=atoms.begin(); at!=atoms.end(); ++at) {
            //!neues ATOM_EXT Objekt:
            if ((*at)->ext.zero()) (*at)->ext = new ATOM_EXT();
            else (*at)->ext->clear();
            //!eventuell vorhandene Bindungen loeschen:
            (*at)->bonded_atoms.clear();
        }
        
        //!eventuell vorhandene Bindungen loeschen:
        for (bonds_vec it=bonds.begin(); it!=bonds.end();++it) {
            it->kill();
        }
        bonds.clear();
        
        //!eventuell vorhandene Ringe loeschen:
        for (rings_vec it=rings.begin(); it!=rings.end();++it) {
            it->kill();
        }
        rings.clear();

        //!Elementtyp zuweisen (zur Sicherheit, falls noch nicht geschehen):
        if (get_ele) get_elements();

        get_connections();

        get_hybridization(get_bonds,max_ring_members,allow_charged_aromatics,kekulize_aromatics);

        set_intern_types(prot_acids,prot_guanidin,prot_amidin,prot_amin,prot_phosphate,prot_sulfate,kekulize_charged);

        if (get_bonds) {
            make_bonds(no_free_rot_planar,prot_acids,prot_guanidin,prot_amidin,
                       prot_amin,prot_phosphate,prot_sulfate,kekulize_aromatics,kekulize_charged);
        }

        intern2mode(fill_X,alt_def_file);
        
        if (kill_ext) {
            for (atoms_vec at=atoms.begin(); at!=atoms.end(); ++at) {
                (*at)->remove_ext();
            }
        }
        atom_types_already_set = true;
    }
}

void MOLECULE::add_H(stl_ptr<ATOM> &atom,vec3d<float> &hvec) {
    //!neues H-Atom:
    stl_ptr<ATOM> hatom(new ATOM(*atom));
    hatom->intern_id = curr_int_id;
    curr_int_id++;
    //!hier muss eigentlich noch der intern_type gesetzt werden!
    hatom->sybyl_type = "H";
    hatom->element = "H";
    ostringstream os;
    os << "H" << hatom->intern_id;
    hatom->name = os.str();
    //hatom->name = "H";
    hatom->coord = hvec;
    hatom->bonded_atoms.push_back(atom);
    atom->bonded_atoms.push_back(hatom);
    atom_buf.push_back(hatom);
    
    //!neues Bond Object erzeugen:
    stl_ptr<BOND> bnd(new BOND());
    bnd->from = hatom;
    bnd->to = atom;
    bnd->type = "1";
    bonds.push_back(bnd);
}

float MOLECULE::get_hpos_score(vec3d<float> &hvec,stl_ptr<ATOM> &atom) {
    //!Die Funktion soll einen Score liefern, der um so niedriger ist je besser die
    //!Position   hvec+atom->coord    fuer ein Wasserstoffatom geeignet ist.
    //!Hierzu werden zunaechst nur die restlichen Atome des Molekuels selbst beruecksichtigt.
    //!Spaeter soll noch die weitere Umgebung des Molekuels beruecksichtigt werden koennen.
    float score = 0.;
    float dbuf;
    vec3d<float> buf(hvec);
    buf += atom->coord;
    
    //!Erstmal alle Clashes bestrafen:
    //!score += d^4 - 81    //Nullpunkt bei 3 Angstrom
    
    //for (ligands_vec lig=main_structure->ligands.begin(); lig!=main_structure->ligands.end(); ++lig) {
    //for (atoms_vec at=(*lig)->atoms.begin(); at!=(*lig)->atoms.end(); ++at) {
    for (atoms_vec at=atoms.begin(); at!=atoms.end(); ++at) {
        if ((*at)->intern_id == atom->intern_id) continue;
        dbuf = get_square_distance(buf,(*at)->coord);
        if (dbuf > 12.96) continue;
        dbuf = (-1. * (dbuf * dbuf)) + 168.;
        score += dbuf;
    }
    //}
    
    //!Ein Abstand zwischen 1.6 und 2.3 angstrom kann auch ein Wasserstoffbruecke sein.
    //!Pruefen, ob moegliche H-Bond und wenn ja guenstig bewerten:
    //for (ligands_vec lig=main_structure->ligands.begin(); lig!=main_structure->ligands.end(); ++lig) {
    //for (atoms_vec at=(*lig)->atoms.begin(); at!=(*lig)->atoms.end(); ++at) {
    for (atoms_vec at=atoms.begin(); at!=atoms.end(); ++at) {
        if ((*at)->intern_id == atom->intern_id) continue;
        if ((*at)->element == "C" || (*at)->element == "H") continue;
        dbuf = get_square_distance(buf,(*at)->coord);
        if (dbuf > 5.29) continue;
        if (dbuf < 2.56) continue;
        dbuf = (2. - fabs(dbuf - 3.61)) * 220.;
        score -= dbuf;
        dbuf = angle(atom->coord,buf,(*at)->coord);
        if (dbuf < 1.7453) score += 240.; //!kleiner 100 Grad
        else if (dbuf < 1.8675) score += 100.; //!kleiner 107 Grad
    }
    //}
    
    return score;
}

void MOLECULE::get_optimal_hpos(vec3d<float> &hvec,stl_ptr<ATOM> &atom,vec3d<float> &axis) {
    //!Diese Funktion wird aufgerufen, um die optimale Position fuer EIN Wasserstoffatom zu ermitteln.
    //!Die aktuelle Position ist:   hvec + atom->coord
    //!Diese Position soll durch Rotation um    axis   geaendert werden.
    //!Hierzu wird in 5� Schritten um 355� gedreht und jeweils ein Score fuer die entsprechende Position ermittelt.
    //!Am Ende erhaelt hvec die Position mit dem besten Score.
    //!In dieser Version sollen zunaechst nur die restlichen Atome des Molekuels selbst beruecksichtigt werden.
    //!(Protonierungen unter Beruecksichtigung der Umgebung spaeter ueber extra Option)
    vec3d<float> buf(hvec);
    float bestval = get_hpos_score(hvec,atom);
    float bufval;
    
    /*
    matrix<float> rm = rotmatrix(axis,0.0872664);
    for (int i=0; i<72; ++i) {
        hvec *= rm;
        bufval = get_hpos_score(hvec,atom);
        if (bufval < bestval) {
            bestval = bufval;
            buf = hvec;
        }
    }
    hvec = buf;
    */
    
    matrix<float> rm = rotmatrix(axis,0.2617993);
    for (int i=0; i<24; ++i) {
        hvec *= rm;
        bufval = get_hpos_score(hvec,atom);
        if (bufval < bestval) {
            bestval = bufval;
            buf = hvec;
        }
    }
    hvec = buf;
    rm = rotmatrix(axis,0.1308996);
    hvec *= rm;
    bufval = get_hpos_score(hvec,atom);
    if (bufval < bestval) {
        bestval = bufval;
        buf = hvec;
    }
    rm = rotmatrix(axis,-0.2617993);
    hvec *= rm;
    bufval = get_hpos_score(hvec,atom);
    if (bufval < bestval) {
        bestval = bufval;
        buf = hvec;
    }
    
    hvec = buf;
    rm = rotmatrix(axis,0.0654498);
    hvec *= rm;
    bufval = get_hpos_score(hvec,atom);
    if (bufval < bestval) {
        bestval = bufval;
        buf = hvec;
    }
    rm = rotmatrix(axis,-0.1308996);
    hvec *= rm;
    bufval = get_hpos_score(hvec,atom);
    if (bufval < bestval) {
        bestval = bufval;
        buf = hvec;
    }
    
    hvec = buf;
    
}

void MOLECULE::get_optimal_hpos2(vec3d<float> &hvec,vec3d<float> &hvec2,vec3d<float> &hvec3,stl_ptr<ATOM> &atom,vec3d<float> &axis) {
    //!Diese Funktion wird aufgerufen, um die optimale Position fuer 3 Wasserstoffatome zu ermitteln.
    //!Die aktuelle Position ist:   hvec(2/3) + atom->coord
    //!Diese Position soll durch Rotation um    axis   geaendert werden.
    //!Hierzu wird in 5� Schritten um 355� gedreht und jeweils ein Score fuer die entsprechende Position ermittelt.
    //!Am Ende erhaelt hvec die Position mit dem besten Score.
    //!In dieser Version sollen zunaechst nur die restlichen Atome des Molekuels selbst beruecksichtigt werden.
    //!(Protonierungen unter Beruecksichtigung der Umgebung spaeter ueber extra Option)
    vec3d<float> buf(hvec);
    vec3d<float> buf2(hvec2);
    vec3d<float> buf3(hvec3);
    float bestval = get_hpos_score(hvec,atom);
    bestval += get_hpos_score(hvec2,atom);
    bestval += get_hpos_score(hvec3,atom);
    float bufval;
    /*
    matrix<float> rm = rotmatrix(axis,0.0872664);
    for (int i=0; i<72; ++i) {
        hvec *= rm;
        hvec2 *= rm;
        hvec3 *= rm;
        bufval = get_hpos_score(hvec,atom);
        bufval += get_hpos_score(hvec2,atom);
        bufval += get_hpos_score(hvec3,atom);
        if (bufval < bestval) {
            
            //if (atom->name == "C6") cout << "taken (" << bufval << " < " << bestval << endl << endl;
        
            bestval = bufval;
            buf = hvec;
            buf2 = hvec2;
            buf3 = hvec3;
        }
    }
    hvec = buf;
    hvec2 = buf2;
    hvec3 = buf3;
    */
    
    matrix<float> rm = rotmatrix(axis,0.2617993);
    for (int i=0; i<24; ++i) {
        hvec *= rm;
        hvec2 *= rm;
        hvec3 *= rm;
        bufval = get_hpos_score(hvec,atom);
        bufval += get_hpos_score(hvec2,atom);
        bufval += get_hpos_score(hvec3,atom);
        if (bufval < bestval) {
            bestval = bufval;
            buf = hvec;
            buf2 = hvec2;
            buf3 = hvec3;
        }
    }
    hvec = buf;
    hvec2 = buf2;
    hvec3 = buf3;
    rm = rotmatrix(axis,0.1308996);
    hvec *= rm;
    hvec2 *= rm;
    hvec3 *= rm;
    bufval = get_hpos_score(hvec,atom);
    bufval += get_hpos_score(hvec2,atom);
    bufval += get_hpos_score(hvec3,atom);
    if (bufval < bestval) {
        bestval = bufval;
        buf = hvec;
        buf2 = hvec2;
        buf3 = hvec3;
    }
    rm = rotmatrix(axis,-0.2617993);
    hvec *= rm;
    hvec2 *= rm;
    hvec3 *= rm;
    bufval = get_hpos_score(hvec,atom);
    bufval += get_hpos_score(hvec2,atom);
    bufval += get_hpos_score(hvec3,atom);
    if (bufval < bestval) {
        bestval = bufval;
        buf = hvec;
        buf2 = hvec2;
        buf3 = hvec3;
    }
    
    hvec = buf;
    hvec2 = buf2;
    hvec3 = buf3;
    rm = rotmatrix(axis,0.0654498);
    hvec *= rm;
    hvec2 *= rm;
    hvec3 *= rm;
    bufval = get_hpos_score(hvec,atom);
    bufval += get_hpos_score(hvec2,atom);
    bufval += get_hpos_score(hvec3,atom);
    if (bufval < bestval) {
        bestval = bufval;
        buf = hvec;
        buf2 = hvec2;
        buf3 = hvec3;
    }
    rm = rotmatrix(axis,-0.1308996);
    hvec *= rm;
    hvec2 *= rm;
    hvec3 *= rm;
    bufval = get_hpos_score(hvec,atom);
    bufval += get_hpos_score(hvec2,atom);
    bufval += get_hpos_score(hvec3,atom);
    if (bufval < bestval) {
        bestval = bufval;
        buf = hvec;
        buf2 = hvec2;
        buf3 = hvec3;
    }
    
    hvec = buf;
    hvec2 = buf2;
    hvec3 = buf3;
    
}

void MOLECULE::set_standard_protonation(int verb) {
    if (atoms.size() == 0) return;
    if (!atom_types_already_set) {
        if (verbosity) cerr << c_message<cERROR>("MOLECULE::set_standard_protonation --> internal atom types are not set!") << endl;
        return;
    }
    atom_buf.clear();
    //!Map erstellen mit den Valenzen:
    tr1::unordered_map<string,int> v_map;
    for (int i=0; i<n_intern_types; ++i) v_map[i_t[i]] = atom_valence[i];
    //!Map erstellen mit Geometrie:
    tr1::unordered_map<string,int> g_map;
    for (int i=0; i<n_intern_types; ++i) g_map[i_t[i]] = atom_geometrie[i];
    //!Map erstellen mit H-Bindungslaengen:
    tr1::unordered_map<string,float> b_map;
    for (int i=0; i<n_intern_types; ++i) b_map[i_t[i]] = atom_hb_length[i];
    
    //!hoechste intern_id ermitteln:
    curr_int_id = 0;
    for (atoms_vec at=atoms.begin(); at!=atoms.end(); ++at) {
        if ((*at)->intern_id > curr_int_id) curr_int_id = (*at)->intern_id;
    }
    curr_int_id++;
    
    atom_buf.clear();
    
    //!Erstmal nur die NICHT frei drehbaren Wasserstoffe setzen:
    //!lineare Geometrie
    //!trigonal planar: ein H zu setzen und 2 andere Atome bereits gebunden
    //!tetraedrisch: ein H zu setzen und 3 andere Atome bereits gebunden
    //!tetraedrisch: zwei H zu setzen und 2 andere Atome bereits gebunden
    for (atoms_vec at=atoms.begin(); at!=atoms.end(); ++at) {
        if (!((*at)->element == "C" || (*at)->element == "N" || (*at)->element == "O" || (*at)->element == "S")) continue;
        if (v_map.find((*at)->intern_type) == v_map.end()) {
            cerr << c_message<cWARNING>("MOLECULE::set_standard_protonation --> atom ") << (*at)->name 
                 << " has no valid internal atom type!" << endl;
            cerr << "No protonation state will be set for this atom!" << endl;
            continue;
        }
        
        //!DEBUG: Schuwi will kein S1 protonieren:
        /*
        string na = (*at)->name;
        string_fu::remove_char(na);
        if (na == "S1") {
        //    cerr << " * ";
            continue;
        }
        */
        
        int h_to_set = v_map[(*at)->intern_type] - (*at)->bonded_atoms.size(); //!Anzahl der zu setzenden Wasserstoffe
        
        if (h_to_set <= 0) continue;
        
        //!An dieser Stelle duerften nur noch CNOS Atome mit mindestens einem zu setzenden Proton sein:
        vec3d<float> hvec(0.,0.,0.);
        if (g_map[(*at)->intern_type] == 1) { //!lineare Geometrie
            if ((*at)->bonded_atoms.size() < 1) continue;
            hvec = (*at)->coord;
            hvec -= (*at)->bonded_atoms[0]->coord;
            hvec += (*at)->coord;
            hvec.norm();
            hvec *= b_map[(*at)->intern_type]; //!Jetzt entspricht hvec den Koordinaten fuer den Wasserstoff
            add_H(*at,hvec);
        } else if (g_map[(*at)->intern_type] == 2) { //!trigonal planare Geometrie
            if (h_to_set == 1) { //!nur ein H anhaengen
                if ((*at)->bonded_atoms.size() == 2) {
                    float rangle = 6.2831853 - angle((*at)->bonded_atoms[0]->coord,(*at)->coord,(*at)->bonded_atoms[1]->coord);
                    rangle /= 2.;
                    vec3d<float> vb1((*at)->bonded_atoms[0]->coord);
                    vb1 -= (*at)->coord;
                    vec3d<float> vb2((*at)->bonded_atoms[1]->coord);
                    vb2 -= (*at)->coord;
                    vec3d<float> rotaxe(vb1);
                    rotaxe *= vb2;
                    hvec = vb2;
                    rotaxe.norm();
                    matrix<float> rm = rotmatrix(rotaxe,rangle);
                    hvec *= rm;
                    hvec.norm();
                    hvec *= b_map[(*at)->intern_type];
                    hvec += (*at)->coord;
                    add_H(*at,hvec);
                } else continue; //erst im naechsten Durchlauf
            } else if (h_to_set == 2) { //! 2 Wasserstoffe anhaengen (z.B.: an sp2-Kohlenstoff)
                vec3d<float> hvec2;
                vec3d<float> xvec((*at)->bonded_atoms[0]->coord);
                vec3d<float> yvec;
                if ((*at)->bonded_atoms[0]->bonded_atoms.size() < 2) {
                    cerr << c_message<cWARNING>("MOLECULE::set_standard_protonation --> atom ") 
                        << (*at)->name << " has unknown bond geometry" << endl;
                    cerr << "No protonation state will be set for this atom!" << endl;
                    continue;
                }
                if ((*at)->bonded_atoms[0]->bonded_atoms[0]->intern_id == (*at)->intern_id) {
                    yvec = (*at)->bonded_atoms[0]->bonded_atoms[1]->coord;
                } else yvec = (*at)->bonded_atoms[0]->bonded_atoms[0]->coord;
                vec3d<float> zvec((*at)->coord);
                vec3d<float> avec(yvec);
                avec -= xvec;
                vec3d<float> bvec(xvec);
                bvec -= zvec;
                vec3d<float> axis(bvec);
                axis *= avec;
                axis.norm();
                matrix<float> rm = rotmatrix(axis,2.0943951);
                hvec = bvec; hvec *= rm; hvec.norm(); hvec *= b_map[(*at)->intern_type];
                hvec2 = hvec; hvec2 *= rm; hvec2.norm(); hvec2 *= b_map[(*at)->intern_type];
                hvec += (*at)->coord;
                add_H(*at,hvec);
                hvec2 += (*at)->coord;
                add_H(*at,hvec2);
            }
        } else if (g_map[(*at)->intern_type] == 3) { //!tetraedrische Geometrie
            if (h_to_set == 1) { //!nur ein H anhaengen
                if ((*at)->bonded_atoms.size() == 3) {
                    vec3d<float> vb1((*at)->bonded_atoms[0]->coord); vb1 -= (*at)->coord; vb1.norm();
                    vec3d<float> vb2((*at)->bonded_atoms[1]->coord); vb2 -= (*at)->coord; vb2.norm();
                    vec3d<float> vb3((*at)->bonded_atoms[2]->coord); vb3 -= (*at)->coord; vb3.norm();
                    hvec = vb1; hvec += vb2; hvec += vb3;
                    hvec.norm();
                    hvec *= -1. * b_map[(*at)->intern_type];
                    hvec += (*at)->coord;
                    add_H(*at,hvec);
                } else if ((*at)->bonded_atoms.size() == 2) { //!z.B. an einen sp3 Stickstoff
                    continue; //erst im naechsten Durchlauf
                } else if ((*at)->bonded_atoms.size() == 1) { //!z.B. an einen sp3 Sauerstoff
                    continue; //erst im naechsten Durchlauf
                }
            } else if (h_to_set == 2) { //! 2 Wasserstoffe anhaengen
                if ((*at)->bonded_atoms.size() == 2) {
                    vec3d<float> vb1((*at)->bonded_atoms[0]->coord); vb1 -= (*at)->coord; vb1.norm();
                    vec3d<float> vb2((*at)->bonded_atoms[1]->coord); vb2 -= (*at)->coord; vb2.norm();
                    vec3d<float> vb3(vb1); vb3 += vb2; vb3.norm();
                    vec3d<float> vb4(vb1); vb4 *= vb2; vb4.norm();
                    hvec = vb1;
                    matrix<float> rm1 = rotmatrix(vb4,3.1415927);
                    hvec *= rm1;
                    matrix<float> rm2 = rotmatrix(vb3,1.5707963);
                    hvec *= rm2;
                    hvec *= b_map[(*at)->intern_type];
                    hvec += (*at)->coord;
                    add_H(*at,hvec);
                    hvec = vb2;
                    hvec *= rm1;
                    hvec *= rm2;
                    hvec *= b_map[(*at)->intern_type];
                    hvec += (*at)->coord;
                    add_H(*at,hvec);
                } else if ((*at)->bonded_atoms.size() == 1) {
                    continue; //erst im naechsten Durchlauf
                } else if ((*at)->bonded_atoms.size() == 0) { //!z.B.: Wasser
                    continue; //erst im naechsten Durchlauf
                }
            } else if (h_to_set == 3) { //! 3 Wasserstoffe anhaengen
                if ((*at)->bonded_atoms.size() == 1) {
                    continue; //erst im naechsten Durchlauf
                } else if ((*at)->bonded_atoms.size() == 0) { //!z.B.: Ammoniak
                    continue; //erst im naechsten Durchlauf
                }
            }
        }
    }
    
    for (atoms_vec at=atom_buf.begin(); at!=atom_buf.end(); ++at) {
        atoms.push_back(*at);
    }
    atom_buf.clear();
    
    //!Jetzt die drehbaren Wasserstoffe setzen:
    //!Hier nach jedem gesetzten diesen erstmal dem Molekuel zufuegen und die Schleife neu starten
    bool goon = true;
    while (goon) {
        goon = false;
        for (atoms_vec at=atoms.begin(); at!=atoms.end(); ++at) {
            if (!((*at)->element == "C" || (*at)->element == "N" || (*at)->element == "O" || (*at)->element == "S")) continue;
            if (v_map.find((*at)->intern_type) == v_map.end()) {
                continue;
            }
            
            
            //!DEBUG: Schuwi will kein S1 protonieren:
            /*
            string na = (*at)->name;
            string_fu::remove_char(na);
            if (na == "S1") {
            //    cerr << " * ";
                continue;
            }
            */
            
            int h_to_set = v_map[(*at)->intern_type] - (*at)->bonded_atoms.size(); //!Anzahl der zu setzenden Wasserstoffe
            
            if (h_to_set <= 0) continue;
            
            //!An dieser Stelle duerften nur noch CNOS Atome mit mindestens einem zu setzenden Proton sein:
            vec3d<float> hvec(0.,0.,0.);
            vec3d<float> hvec2(0.,0.,0.);
            if (g_map[(*at)->intern_type] == 2) { //!trigonal planare Geometrie
                if (h_to_set == 1) { //!nur ein H anhaengen
                    if ((*at)->bonded_atoms.size() == 1) { //!z.B. Imin
                        //Es kommen 2 Positionen in Frage von denen die bessere zu ermitteln ist:
                        vec3d<float> xvec((*at)->bonded_atoms[0]->coord);
                        vec3d<float> yvec;
                        if ((*at)->bonded_atoms[0]->bonded_atoms.size() < 2) {
                            cerr << c_message<cWARNING>("MOLECULE::set_standard_protonation --> atom ") 
                                 << (*at)->name << " has unknown bond geometry" << endl;
                            cerr << "No protonation state will be set for this atom!" << endl;
                            continue;
                        }
                        if ((*at)->bonded_atoms[0]->bonded_atoms[0]->intern_id == (*at)->intern_id) {
                            yvec = (*at)->bonded_atoms[0]->bonded_atoms[1]->coord;
                        } else yvec = (*at)->bonded_atoms[0]->bonded_atoms[0]->coord;
                        vec3d<float> zvec((*at)->coord);
                        vec3d<float> avec(yvec);
                        avec -= xvec;
                        vec3d<float> bvec(xvec);
                        bvec -= zvec;
                        vec3d<float> axis(bvec);
                        axis *= avec;
                        axis.norm();
                        matrix<float> rm = rotmatrix(axis,2.0943951);
                        hvec = bvec; hvec *= rm; hvec.norm(); hvec *= b_map[(*at)->intern_type];
                        hvec2 = hvec; hvec2 *= rm; hvec2.norm(); hvec2 *= b_map[(*at)->intern_type];
                        if (get_hpos_score(hvec,*at) < get_hpos_score(hvec2,*at)) {
                            hvec += (*at)->coord;
                            add_H(*at,hvec);
                        } else {
                            hvec2 += (*at)->coord;
                            add_H(*at,hvec2);
                        }
                        goon = true;
                        break;
                    }
                } else continue;
            } else if (g_map[(*at)->intern_type] == 3) { //!tetraedrische Geometrie
                if (h_to_set == 1) { //!nur ein H anhaengen
                    if ((*at)->bonded_atoms.size() == 2) { //!z.B. an einen sp3 Stickstoff
                        //Es kommen 2 Positionen in Frage von denen die bessere zu ermitteln ist:
                        vec3d<float> vb1((*at)->bonded_atoms[0]->coord); vb1 -= (*at)->coord; vb1.norm();
                        vec3d<float> vb2((*at)->bonded_atoms[1]->coord); vb2 -= (*at)->coord; vb2.norm();
                        vec3d<float> vb3(vb1); vb3 += vb2; vb3.norm();
                        vec3d<float> vb4(vb1); vb4 *= vb2; vb4.norm();
                        hvec = vb1;
                        matrix<float> rm1 = rotmatrix(vb4,3.1415927);
                        hvec *= rm1;
                        matrix<float> rm2 = rotmatrix(vb3,1.5707963);
                        hvec *= rm2;
                        hvec *= b_map[(*at)->intern_type];
                        hvec2 = vb2;
                        hvec2 *= rm1;
                        hvec2 *= rm2;
                        hvec2 *= b_map[(*at)->intern_type];
                        if (get_hpos_score(hvec,*at) < get_hpos_score(hvec2,*at)) {
                            hvec += (*at)->coord;
                            add_H(*at,hvec);
                        } else {
                            hvec2 += (*at)->coord;
                            add_H(*at,hvec2);
                        }
                        goon = true;
                        break;
                    } else if ((*at)->bonded_atoms.size() == 1) { //!z.B. an einen sp3 Sauerstoff
                        //Erstmal irgendwie im 109,5� Winkel setzen und dann auf "optimale" Position drehen
                        vec3d<float> vb1((*at)->bonded_atoms[0]->coord);
                        vb1 -= (*at)->coord; vb1.norm();
                        hvec = vb1;
                        // 1.) Vektor willkuerlich setzen:
                        vec3d<float> vb2(1.,0.,0.);
                        if (fabs(vb1[1]) < 0.1 && fabs(vb1[2]) < 0.1) {vb2[0]  = 0.; vb2[1] = 1.;}
                        // 2.) Rotationsachse bestimmen:
                        vec3d<float> vb3(vb1); vb3 *= vb2;
                        matrix<float> rm = rotmatrix(vb3,1.911135);
                        // 3.) hvec so drehen, dass er 109,5� zum "at---bonded[0]" Vektor hat
                        hvec *= rm;
                        hvec *= b_map[(*at)->intern_type];
                        // 4.) Jetzt um vb1 drehen bis gute Position:
                        get_optimal_hpos(hvec,*at,vb1);
                        hvec += (*at)->coord;
                        add_H(*at,hvec);
                        goon = true;
                        break;
                    }
                } else if (h_to_set == 2) { //! 2 Wasserstoffe anhaengen
                    if ((*at)->bonded_atoms.size() == 1) { //!z.B.: NH2
                        //!erstmal nur eins anhaengen im 109� Winkel
                        //!das zweite sollte dann automatisch im naechsten Durchlauf gesetzt werden
                        //! => wie fuer ein zu setzendes H bei einem bereits gebundenen Schweratom
                        vec3d<float> vb1((*at)->bonded_atoms[0]->coord);
                        vb1 -= (*at)->coord; vb1.norm();
                        hvec = vb1;
                        vec3d<float> vb2(1.,0.,0.);
                        if (fabs(vb1[1]) < 0.1 && fabs(vb1[2]) < 0.1) {vb2[0]  = 0.; vb2[1] = 1.;}
                        vec3d<float> vb3(vb1); vb3 *= vb2; vb3.norm();
                        matrix<float> rm = rotmatrix(vb3,1.911135);
                        hvec *= rm;
                        hvec *= b_map[(*at)->intern_type];
                        get_optimal_hpos(hvec,*at,vb1);
                        hvec += (*at)->coord;
                        add_H(*at,hvec);
                        goon = true;
                        break;
                    } else if ((*at)->bonded_atoms.size() == 0) { //!z.B.: Wasser
                        //!auch hier erstmal ein H setzen und das zweite wird dann automatisch im
                        //!naechsten Durchlauf gesetzt
                        //!spaeter mal...
                        
                        hvec = vec3d<float>(1.,0.,0.);
                        hvec *= b_map[(*at)->intern_type];
                        hvec += (*at)->coord;
                        add_H(*at,hvec);
                        goon = true;
                        break;
                        
                    ///    cerr << c_message<cWARNING>("MOLECULE::set_standard_protonation --> atom ") 
                    ///         << (*at)->name << " has unknown bond geometry" << endl;
                    ///    cerr << "No protonation state will be set for this atom!" << endl;
                    ///    //goon = true;
                    ///    break;
                    }
                } else if (h_to_set == 3) { //! 3 Wasserstoffe anhaengen
                    if ((*at)->bonded_atoms.size() == 1) { //!z.B.: CH3
                        //Erstmal irgendwie im 109,5� Winkel setzen und dann auf "optimale" Position drehen
                        vec3d<float> vb1((*at)->bonded_atoms[0]->coord);
                        vb1 -= (*at)->coord; vb1.norm();
                        hvec = vb1;
                        // 1.) Vektor willkuerlich setzen:
                        vec3d<float> vb2(1.,0.,0.);
                        if (fabs(vb1[1]) < 0.1 && fabs(vb1[2]) < 0.1) {vb2[0]  = 0.; vb2[1] = 1.;}
                        // 2.) Rotationsachse bestimmen:
                        vec3d<float> vb3(vb1); vb3 *= vb2;
                        matrix<float> rm = rotmatrix(vb3,1.911135);
                        // 3.) hvec so drehen, dass er 109,5� zum "at---bonded[0]" Vektor hat
                        hvec *= rm;
                        hvec *= b_map[(*at)->intern_type];
                        // 4.) Nun das Tetraeder vervollstaendigen:
                        vec3d<float> vb5((*at)->bonded_atoms[0]->coord); vb5 -= (*at)->coord; vb5.norm();
                        vec3d<float> vb6(hvec); vb6.norm();
                        vec3d<float> vb7(vb1); vb7 += vb6; vb7.norm();
                        vec3d<float> vb4(vb1); vb4 *= vb6; vb4.norm();
                        hvec2 = vb5;
                        matrix<float> rm1 = rotmatrix(vb4,3.1415927);
                        hvec2 *= rm1;
                        matrix<float> rm2 = rotmatrix(vb7,1.5707963);
                        hvec2 *= rm2;
                        hvec2 *= b_map[(*at)->intern_type];
                        vec3d<float> hvec3(vb6);
                        hvec3 *= rm1;
                        hvec3 *= rm2;
                        hvec3 *= b_map[(*at)->intern_type];
                        // 5.) Jetzt um vb1 drehen bis gute Position:
                        get_optimal_hpos2(hvec,hvec2,hvec3,*at,vb1);
                        hvec += (*at)->coord;
                        hvec2 += (*at)->coord;
                        hvec3 += (*at)->coord;
                        add_H(*at,hvec);
                        add_H(*at,hvec2);
                        add_H(*at,hvec3);
                        goon = true;
                        break;
                    } else if ((*at)->bonded_atoms.size() == 0) { //!z.B.: Ammoniak
                        //!erstmal nur ein H setzen
                        //!spaeter mal...
                        cerr << c_message<cWARNING>("MOLECULE::set_standard_protonation --> atom ") 
                             << (*at)->name << " has unknown bond geometry" << endl;
                        cerr << "No protonation state will be set for this atom!" << endl;
                        //goon = true;
                        break;
                    }
                }
            }
        }
        if (goon) {
            for (atoms_vec at=atom_buf.begin(); at!=atom_buf.end(); ++at) {
                atoms.push_back(*at);
            }
            atom_buf.clear();
        }
    }
    
    map<int,vector<stl_ptr<ATOM> > > res_buf;
    for (atoms_vec at=atoms.begin(); at!=atoms.end(); ++at) {
        if (res_buf.find((*at)->res_number) == res_buf.end()) res_buf[(*at)->res_number] = vector<stl_ptr<ATOM> >();
        res_buf[(*at)->res_number].push_back(*at);
    }
    atoms.clear();
    for (map<int,vector<stl_ptr<ATOM> > >::iterator it=res_buf.begin(); it!=res_buf.end(); ++it) {
        for (atoms_vec at=it->second.begin(); at!=it->second.end(); ++at) {
            atoms.push_back(*at);
        }
    }
}

void MOLECULE::kick_water(bool intern_is_set) {
    if (atoms.size() == 0) return;
    if (!intern_is_set) { //noch keine internen Typen
        for (atoms_vec at=atoms.begin(); at!=atoms.end(); ++at) {
            //!neues ATOM_EXT Objekt:
            if ((*at)->ext.zero()) (*at)->ext = new ATOM_EXT();
            else (*at)->ext->clear();
            //!eventuell vorhandene Bindungen loeschen:
            (*at)->bonded_atoms.clear();
        }
        get_elements();
        get_connections();
        for (atoms_vec at=atoms.begin(); at!=atoms.end(); ++at) {
            if ((*at)->element == "O" && (*at)->ext->n_heavy_bonded == 0) (*at)->intern_type = "O.h2o";
            (*at)->remove_ext();
        }
    }
    //interne Typen sind gesetzt => einfach O.h2o rausschmeissen
    for (atoms_vec at=atoms.begin(); at!=atoms.end(); ++at) {
        if ((*at)->intern_type == "O.h2o") {
            //beteiligte bonds killen:
            for (bonds_vec bnd=bonds.begin(); bnd!=bonds.end(); ++bnd) {
                if ((*at == (*bnd)->from) || (*at == (*bnd)->to)) {
                    stl_ptr<BOND> bon = *bnd;
                    bonds.erase(bnd);
                    bnd--;
                    bon.kill();
                }
            }
            for (atoms_vec bt=(*at)->bonded_atoms.begin(); bt!=(*at)->bonded_atoms.end(); ++bt) {
                (*bt)->intern_type ="kill";
            }
            (*at)->intern_type ="kill";
        }
    }
    for (atoms_vec at=atoms.begin(); at!=atoms.end(); ++at) {
        if ((*at)->intern_type == "kill") {
            stl_ptr<ATOM> atm = *at;
            atoms.erase(at);
            at--;
            atm.kill();
            continue;
        }
        if (!intern_is_set) (*at)->intern_type == "X";
    }
}

void MOLECULE::rename_atoms(bool get_ele) {
    if (get_ele) get_elements();
    tr1::unordered_map<string,int> ele_count;
    for (atoms_vec at=atoms.begin(); at!=atoms.end(); ++at) {
        if (ele_count.find((*at)->element) == ele_count.end()) ele_count[(*at)->element] = 1;
        else ele_count[(*at)->element] +=1;
        ostringstream os;
        os << (*at)->element << ele_count[(*at)->element];
        (*at)->name = os.str();
    }
}

void MOLECULE::reorder_atoms(bool get_ele) {
    if (get_ele) get_elements();
    vector<stl_ptr<ATOM> > other;
    vector<stl_ptr<ATOM> > hydrogens;
    for (atoms_vec at=atoms.begin(); at!=atoms.end(); ++at) {
        if ((*at)->element == "H") hydrogens.push_back(*at);
        else other.push_back(*at);
    }
    atoms.clear();
    for (atoms_vec at=other.begin(); at!=other.end(); ++at) atoms.push_back(*at);
    for (atoms_vec at=hydrogens.begin(); at!=hydrogens.end(); ++at) atoms.push_back(*at);
}

int MOLECULE::get_n_heavy() {
    if (n_heavy_atoms > -1) return n_heavy_atoms;
    n_heavy_atoms = 0;
    for (atoms_vec at=atoms.begin(); at!=atoms.end(); ++at) {
        if ((*at)->sybyl_type[0] != 'H') ++n_heavy_atoms;
    }
    return n_heavy_atoms;
}

int MOLECULE::get_number_of_rings(int max_ring_members,bool debug_mode) {
    if (!atom_types_already_set) {
        verbosity = 0;
        
        for (atoms_vec at=atoms.begin(); at!=atoms.end(); ++at) {
            //!neues ATOM_EXT Objekt:
            if ((*at)->ext.zero()) (*at)->ext = new ATOM_EXT();
            else (*at)->ext->clear();
            //!eventuell vorhandene Bindungen loeschen:
            (*at)->bonded_atoms.clear();
        }
        
        for (rings_vec it=rings.begin(); it!=rings.end();++it) {
            it->kill();
        }
        rings.clear();
    
        get_elements();
        
        get_connections();
        get_hybridization(false,max_ring_members);
    }
    
    if (debug_mode) {
        cout << "N_atoms N_hetero_atoms planar aromatic at1_name(at1_id) ... atN_name(atN_id)" << endl;
        for (rings_vec rt=rings.begin(); rt!=rings.end(); ++rt) {
            cout << (*rt)->n_members << " " << (*rt)->n_hetero << " " << (*rt)->is_planar << " " << (*rt)->is_aromatic << " ";
            for (atoms_vec at=(*rt)->ring_atoms.begin(); at!=(*rt)->ring_atoms.end(); ++at) {
                cout << (*at)->name << "(" << (*at)->id << ") ";
            }
            cout << endl;
        }
    }
    
    return rings.size();
}

float MOLECULE::get_mol_weight() {
    get_atom_typing(0,true,"X",false,0,false,true);
    atom_properties::initialize();
    float weight = 0.;
    for (atoms_vec at=atoms.begin(); at!=atoms.end(); ++at) {
        if ((*at)->intern_type[0] == 'X' || (*at)->intern_type[0] == 'H') continue;
        if (atom_properties::mol_weight.find((*at)->element) == atom_properties::mol_weight.end()) continue;
        weight += atom_properties::mol_weight[(*at)->element];
        //!Jetzt noch die Wasserstoffe dazuzaehlen:
        if ((*at)->element == "C") {
            int modi = ((*at)->ext->hybridization + 1 - (*at)->ext->n_heavy_bonded);
            if (modi > 0) weight += atom_properties::mol_weight["H"] * modi;
        } else if ((*at)->element == "N") {
            if ((*at)->intern_type == "N.ar6p" || (*at)->intern_type == "N.arp" || (*at)->intern_type == "N.ar3" ||
                (*at)->intern_type == "N.ar3h" || (*at)->intern_type == "N.amp" || (*at)->intern_type == "N.ams" ||
                (*at)->intern_type == "N.amt" || (*at)->intern_type == "N.samp" || (*at)->intern_type == "N.sams" ||
                (*at)->intern_type == "N.samt" || (*at)->intern_type == "N.gu2" || (*at)->intern_type == "N.guh" ||
                (*at)->intern_type == "N.mi2" || (*at)->intern_type == "N.mih" || (*at)->intern_type == "N.aap" ||
                (*at)->intern_type == "N.aas2" || (*at)->intern_type == "N.aas3" || (*at)->intern_type == "N.aat2" ||
                (*at)->intern_type == "N.aat3" || (*at)->intern_type == "N.3n" || (*at)->intern_type == "N.3p" ||
                (*at)->intern_type == "N.3s" || (*at)->intern_type == "N.3t" || (*at)->intern_type == "N.ims" ||
                (*at)->intern_type == "N.imt" || (*at)->intern_type == "N.2t") {
                int modi = 3 - (*at)->ext->n_heavy_bonded;
                if (modi > 0) weight += atom_properties::mol_weight["H"] * modi;
            } else if ((*at)->intern_type == "N.4q" || (*at)->intern_type == "N.4h") {
                int modi = 4 - (*at)->ext->n_heavy_bonded;
                if (modi > 0) weight += atom_properties::mol_weight["H"] * modi;
            } else if ((*at)->ext->hybridization > 1) {
                int modi = 2 - (*at)->ext->n_heavy_bonded;
                if (modi > 0) weight +=atom_properties:: mol_weight["H"] * modi;
            }
        } else if ((*at)->element == "O" || (*at)->element == "S") {
            if ((*at)->ext->hybridization == 3) {
                int modi = 2 - (*at)->ext->n_heavy_bonded;
                if (modi > 0) weight += atom_properties::mol_weight["H"] * modi;
            }
        }
    }
    return weight;
}

float MOLECULE::get_vdw_volume() {
    //! Nach Inclusion-Exclusion Prinzip:
    
    get_elements();
    //! 1.) Triangulieren:
    atom_properties::initialize();
    vector<DT_POINT> dtpoints;
    int id = 1;
    for (atoms_vec at=atoms.begin(); at!=atoms.end(); ++at) {
        if ((*at)->element == "H") continue;
        double w = atom_properties::square_vdW_map[(*at)->element];
        double r = atom_properties::vdW_map[(*at)->element];
        DT_POINT tmp(id,(*at)->coord[0],(*at)->coord[1],(*at)->coord[2],w,r,(ATOM*)&(**at));
        ++id;
        dtpoints.push_back(tmp);
    }
    
    DT_SOLVER mysolver(dtpoints);
    mysolver.solve();
    
    //! 2.) Volumen berechnen:
    vector<stl_ptr<TETRAHEDRON> > pt;
    for (th_vec it=mysolver.delaunay_tetrahedrons.begin(); it!=mysolver.delaunay_tetrahedrons.end(); ++it) {
        pt.push_back(*it);
    }
    for (th_vec it=mysolver.hull_tetrahedrons.begin(); it!=mysolver.hull_tetrahedrons.end(); ++it) {
        pt.push_back(*it);
    }
    
    float vol = mysolver.get_exact_occ_volume(pt);
    
    return vol;
}

float MOLECULE::get_vdw_volume_grid() {
    //! Grid-basierte Variante zur Ueberpruefung:
    
    get_elements();
    
    vec3d<float> min(atoms[0]->coord);
    vec3d<float> max(atoms[0]->coord);
    for (atoms_vec at=atoms.begin(); at!=atoms.end(); ++at) {
        if ((*at)->element == "H") continue;
        if ((*at)->coord[0] < min[0]) min[0] = (*at)->coord[0];
        else if ((*at)->coord[0] > max[0]) max[0] = (*at)->coord[0];
        if ((*at)->coord[1] < min[1]) min[1] = (*at)->coord[1];
        else if ((*at)->coord[1] > max[1]) max[1] = (*at)->coord[1];
        if ((*at)->coord[2] < min[2]) min[2] = (*at)->coord[2];
        else if ((*at)->coord[2] > max[2]) max[2] = (*at)->coord[2];
    }
    
    min[0] -= 3.; min[1] -= 3.; min[2] -= 3.;
    max[0] += 3.; max[1] += 3.; max[2] += 3.;
    
    double vol = 0.;
    
    double step_size = 0.1;//0.08;
    double ele_vol = step_size*step_size*step_size;
    
    int sx = ((max[0]-min[0]) / step_size) + 1;
    int sy = ((max[1]-min[1]) / step_size) + 1;
    int sz = ((max[2]-min[2]) / step_size) + 1;
    
    int64_t n_p = 0;

    OCTREE<ATOM>* pwp = new OCTREE<ATOM>(atoms,2.2);

    for (int x=0; x<sx; ++x) {
        for (int y=0; y<sy; ++y) {
            for (int z=0; z<sz; ++z) {
                vec3d<float> tv(min[0]+(x*step_size),
                                min[1]+(y*step_size),
                                min[2]+(z*step_size));
                for (vector<ATOM*>::iterator at=pwp->begin(1.0,tv); at!=pwp->end(); ++at) {
                    if ((*at)->element == "H") continue;
                    if (get_square_distance(tv,(*at)->coord) < atom_properties::square_vdW_map[(*at)->element]) {
                        ++n_p;
                        break;
                    }
                }
            }
        }
    }

    delete pwp;
    
    vol = n_p * ele_vol;
    
    return vol;
}

void MOLECULE::rmsd_by_mol2_residues(vector<string> &vrn,vector<float> &vtr,vector<float> &vbr,vector<float> &vrr,
                                     vector<float> &otr,vector<float> &obr,vector<float> &orr,stl_ptr<LIGAND> &ref,bool debug_mode) {
    //! 22.04.09
    //! Hier ist das Sequenzalignment in einer weniger umstaendlichen Form mit einem
    //! SA_CONTAINER Objekt anstatt eines String => Kein muehseliges zurueck-mappen
    //! auf die eigentlichen ATOM-Objekte!!!
    //! Die anderen Sequenz-Alignment-Funktionen sollten dringend entsprechend angepasst werden!!!
    
    vector<SA_CONTAINER<stl_ptr<ATOM>,string> > Aseq;
    vector<SA_CONTAINER<stl_ptr<ATOM>,string> > Bseq;
    
    SA_CONTAINER<stl_ptr<ATOM>,string>  gap_val;
    ATOM gap_at;
    gap_at.res_name = "XXX";
    stl_ptr<ATOM> tatmp;
    tatmp = &gap_at;
    gap_val.at = tatmp;
    gap_val.comparer = "XXX";
    
    tr1::unordered_map<string,vector<stl_ptr<ATOM> > > full_resB;
    int last_res_num = 1;
    for (atoms_vec at=atoms.begin(); at!=atoms.end(); ++at) {
        if ((*at)->res_name.size() > 3) (*at)->res_name = (*at)->res_name.substr(0,3);
        ostringstream hos;
        hos << (*at)->res_name << (*at)->res_number;
        (*at)->full_res_name = hos.str();
        if (full_resB.find((*at)->full_res_name) == full_resB.end()) {
            int n_gaps = abs((*at)->res_number - last_res_num);
            if (n_gaps > 9) n_gaps = 10;
            for (int i=1; i<n_gaps; ++i) {
                Bseq.push_back(gap_val);
            //    break; //! => also nur ein gap_val
            }
            last_res_num = (*at)->res_number;
            full_resB[(*at)->full_res_name] = vector<stl_ptr<ATOM> >();
            SA_CONTAINER<stl_ptr<ATOM>,string> csa;
            csa.at = *at;
            csa.comparer = (*at)->res_name;
            Bseq.push_back(csa);
        }
        full_resB[(*at)->full_res_name].push_back(*at);
    }
    tr1::unordered_map<string,vector<stl_ptr<ATOM> > > full_resA;
    last_res_num = 1;
    for (atoms_vec at=ref->atoms.begin(); at!=ref->atoms.end(); ++at) {
        if ((*at)->res_name.size() > 3) (*at)->res_name = (*at)->res_name.substr(0,3);
        ostringstream hos;
        hos << (*at)->res_name << (*at)->res_number;
        (*at)->full_res_name = hos.str();
        if (full_resA.find((*at)->full_res_name) == full_resA.end()) {
            int n_gaps = abs((*at)->res_number - last_res_num);
            if (n_gaps > 9) n_gaps = 10;
            for (int i=1; i<n_gaps; ++i) {
                Aseq.push_back(gap_val);
            //    break; //! => also nur ein gap_val
            }
            last_res_num = (*at)->res_number;
            full_resA[(*at)->full_res_name] = vector<stl_ptr<ATOM> >();
            SA_CONTAINER<stl_ptr<ATOM>,string> csa;
            csa.at = *at;
            csa.comparer = (*at)->res_name;
            Aseq.push_back(csa);
        }
        full_resA[(*at)->full_res_name].push_back(*at);
    }
    
    SEQALIGN<SA_CONTAINER<stl_ptr<ATOM>,string> > my_sa(Aseq,Bseq,gap_val);
    my_sa.needleman_wunsch1(); //!Arbeitet hier eher schlecht => lieber Smith-Waterman nehmen
                               //!Dann allerdings nur je einen kuenstlichen Gap fuer unterbrochene
                               //!Ketten einsetzen!
    //my_sa.smith_waterman();
    
    if (debug_mode) {
        for (unsigned int i=0; i<my_sa.ires.size(); ++i) {
            if (i >= my_sa.jres.size()) break;
            if (my_sa.ires[i].at->res_name != "XXX") {
                cout << my_sa.ires[i].at->res_name << " " << my_sa.ires[i].at->res_number << " <--> ";
            } else cout << "XXX <--> ";
            if (my_sa.jres[i].at->res_name != "XXX") {
                cout << my_sa.jres[i].at->res_name << " " << my_sa.jres[i].at->res_number << endl;
            } else cout << "XXX" << endl;
        }
    }
    
    vector<vec3d<float> > A_list;
    vector<vec3d<float> > B_list;
    for (unsigned int i=0; i<my_sa.ires.size(); ++i) {
        if (i >= my_sa.jres.size()) break;
        if (my_sa.ires[i].at->res_name == my_sa.jres[i].at->res_name && 
            !(my_sa.ires[i].at->res_name == "XXX" && my_sa.jres[i].at->res_name == "XXX")) {
            float res_rmsd = 0.;
            float back_rmsd = 0.;
            float total_rmsd = 0.;
            int n_back = 0;
            int n_res = 0;
            int n_tot = 0;
            unsigned max_atm = full_resA[my_sa.ires[i].at->full_res_name].size();
            if (full_resB[my_sa.jres[i].at->full_res_name].size() < max_atm) {
                max_atm = full_resB[my_sa.jres[i].at->full_res_name].size();
            }
            for (unsigned int j=0; j<max_atm; ++j) {
                istringstream is;
                string name;
                is.str(full_resA[my_sa.ires[i].at->full_res_name][j]->name);
                is >> name;
                if (name == "N" || name == "CA" || name == "C" || name == "O") {
                    ++n_tot;
                    ++n_back;
                    back_rmsd += get_square_distance(full_resA[my_sa.ires[i].at->full_res_name][j]->coord,
                                                         full_resB[my_sa.jres[i].at->full_res_name][j]->coord);
                } else {
                    ++n_tot;
                    ++n_res;
                    res_rmsd += get_square_distance(full_resA[my_sa.ires[i].at->full_res_name][j]->coord,
                                                        full_resB[my_sa.jres[i].at->full_res_name][j]->coord);
                }
                A_list.push_back(full_resA[my_sa.ires[i].at->full_res_name][j]->coord);
                B_list.push_back(full_resB[my_sa.jres[i].at->full_res_name][j]->coord);
            }
            if (n_tot == 0) total_rmsd = -1.;
            else {
                total_rmsd = back_rmsd + res_rmsd;
                total_rmsd /= n_tot;
                total_rmsd = sqrt(total_rmsd);
            }
            if (n_back == 0) back_rmsd = -1.;
            else {
                back_rmsd /= n_back;
                back_rmsd = sqrt(back_rmsd);
            }
            if (n_res == 0) res_rmsd = -1.;
            else {
                res_rmsd /= n_res;
                res_rmsd = sqrt(res_rmsd);
            }
            vrn.push_back(my_sa.jres[i].at->full_res_name);
            vtr.push_back(total_rmsd);
            vbr.push_back(back_rmsd);
            vrr.push_back(res_rmsd);
        }
    }
    
    //! Und jetzt den ganzen Unsinn nochmal mit vorhergehendem geometrischen Alignment:
    matrix<float> optalign_rotm;
    vec3d<float> optalign_trans[2];
    get_align_matrix(A_list,B_list,optalign_rotm,optalign_trans[0],optalign_trans[1]);
    
    for (unsigned int i=0; i<my_sa.ires.size(); ++i) {
        if (i >= my_sa.jres.size()) break;
        if (my_sa.ires[i].at->res_name == my_sa.jres[i].at->res_name && 
            !(my_sa.ires[i].at->res_name == "XXX" && my_sa.jres[i].at->res_name == "XXX")) {
            float res_rmsd = 0.;
            float back_rmsd = 0.;
            float total_rmsd = 0.;
            int n_back = 0;
            int n_res = 0;
            int n_tot = 0;
            unsigned max_atm = full_resA[my_sa.ires[i].at->full_res_name].size();
            if (full_resB[my_sa.jres[i].at->full_res_name].size() < max_atm) {
                max_atm = full_resB[my_sa.jres[i].at->full_res_name].size();
            }
            for (unsigned int j=0; j<max_atm; ++j) {
                istringstream is;
                string name;
                is.str(full_resA[my_sa.ires[i].at->full_res_name][j]->name);
                is >> name;
                
                full_resB[my_sa.jres[i].at->full_res_name][j]->coord -= optalign_trans[1];
                full_resB[my_sa.jres[i].at->full_res_name][j]->coord *= optalign_rotm;
                full_resB[my_sa.jres[i].at->full_res_name][j]->coord += optalign_trans[0];
                
                if (name == "N" || name == "CA" || name == "C" || name == "O") {
                    ++n_tot;
                    ++n_back;
                    back_rmsd += get_square_distance(full_resA[my_sa.ires[i].at->full_res_name][j]->coord,
                                                         full_resB[my_sa.jres[i].at->full_res_name][j]->coord);
                } else {
                    ++n_tot;
                    ++n_res;
                    res_rmsd += get_square_distance(full_resA[my_sa.ires[i].at->full_res_name][j]->coord,
                                                        full_resB[my_sa.jres[i].at->full_res_name][j]->coord);
                }
            }
            if (n_tot == 0) total_rmsd = -1.;
            else {
                total_rmsd = back_rmsd + res_rmsd;
                total_rmsd /= n_tot;
                total_rmsd = sqrt(total_rmsd);
                if (total_rmsd < 0.0009) total_rmsd = 0.;
            }
            if (n_back == 0) back_rmsd = -1.;
            else {
                back_rmsd /= n_back;
                back_rmsd = sqrt(back_rmsd);
                if (back_rmsd < 0.0009) back_rmsd = 0.;
            }
            if (n_res == 0) res_rmsd = -1.;
            else {
                res_rmsd /= n_res;
                res_rmsd = sqrt(res_rmsd);
                if (res_rmsd < 0.0009) res_rmsd = 0.;
            }
            otr.push_back(total_rmsd);
            obr.push_back(back_rmsd);
            orr.push_back(res_rmsd);
        }
    }
}


//==============================================================================================
//Definitionen fuer LIGAND:
//==============================================================================================

LIGAND::LIGAND():MOLECULE(), is_peptide(false), charge_type("NO_CHARGES") {}

LIGAND::LIGAND(LIGAND const& lig):MOLECULE(lig), type(lig.type), charge_type(lig.charge_type) { //Copy-Konstruktor
    if (!lig.crysin.zero()) crysin = new CRYSIN(*(lig.crysin));
}

LIGAND::LIGAND(LIGAND const& lig,tr1::unordered_set<int> &at2take):MOLECULE(), type(lig.type), charge_type(lig.charge_type) { //Copy-Konstruktor
    //! Nur die Atome uebernehmen, deren intern_id in at2take steht!!!
    for (const_atoms_vec at=lig.atoms.begin(); at!=lig.atoms.end(); ++at) {
        if (at2take.find((*at)->intern_id) == at2take.end()) continue;
        stl_ptr<ATOM> atm(new ATOM(**at));
        atoms.push_back(atm);
    }
    
    //!Keine BOND-Objekte kopieren
    
    for (map<int,stl_ptr<ATOM> >::const_iterator it=lig.sub_map.begin(); it!=lig.sub_map.end(); ++it) {
        for (atoms_vec at=atoms.begin(); at!=atoms.end(); ++at) {
            if (it->first == (*at)->intern_id) sub_map[it->first] = *at;
        }
    }
    
    for (const_comments_vec ct=lig.comments.begin(); ct!=lig.comments.end(); ++ct) {
        stl_ptr<COMMENT> cnt(new COMMENT());
        cnt->text = (*ct)->text;
        comments.push_back(cnt);
    }
    
    if (!lig.crysin.zero()) crysin = new CRYSIN(*(lig.crysin));
    
    main_structure = lig.main_structure;
}

LIGAND::~LIGAND() {
    if (!crysin.zero()) crysin.kill();
}

void LIGAND::get_connected_atoms(tr1::unordered_set<int> &con,stl_ptr<ATOM> &base) {
    for (atoms_vec at=base->bonded_atoms.begin(); at!=base->bonded_atoms.end(); ++at) {
        if (con.find((*at)->intern_id) != con.end()) continue;
        con.insert((*at)->intern_id);
        get_connected_atoms(con,*at);
    }
}

int LIGAND::split(bool keep_biggest_fragment) {
    //! Vorher get_atom_typing aufrufen!!!
    if (!atom_types_already_set) {
        if (verbosity) {cerr << c_message<cERROR>("LIGAND::split --> no atom types set for ") 
                                 << name << endl; return 0;}
    }
    
    int max = 0;
    
    map<int,tr1::unordered_set<int> > mol_at;
    mol_at[0] = tr1::unordered_set<int>();
    mol_at[0].insert(atoms[0]->intern_id);
    get_connected_atoms(mol_at[0],atoms[0]);
    stl_ptr<ATOM> unvisited = 0;
    for (atoms_vec at=atoms.begin(); at!=atoms.end(); ++at) {
        if (mol_at[0].find((*at)->intern_id) == mol_at[0].end()) {
            unvisited = *at;
            break;
        }
    }
    
    int map_n = 0;
    while (!unvisited.zero()) {
        ++map_n;
        mol_at[map_n] = tr1::unordered_set<int>();
        mol_at[map_n].insert(unvisited->intern_id);
        get_connected_atoms(mol_at[map_n],unvisited);
        
        if (mol_at[map_n].size() > mol_at[max].size()) max = map_n;
        
        unvisited = 0;
        for (atoms_vec at=atoms.begin(); at!=atoms.end(); ++at) {
            bool take_it = true;
            for (map<int,tr1::unordered_set<int> >::iterator it=mol_at.begin(); it!=mol_at.end(); ++it) {
                if (it->second.find((*at)->intern_id) != it->second.end()) {
                    take_it = false;
                    break;
                }
            }
            if (take_it) {
                unvisited = *at;
                break;
            }
        }
    }
    
    if (keep_biggest_fragment) {
        stl_ptr<LIGAND> lig(new LIGAND(*this,mol_at[max]));
        main_structure->splitted_ligands.push_back(lig);
    } else {
        int split_num = 0;
        for (map<int,tr1::unordered_set<int> >::iterator it=mol_at.begin(); it!=mol_at.end(); ++it) {
            stl_ptr<LIGAND> lig(new LIGAND(*this,it->second));
            ostringstream hos; hos << split_num; lig->name += "_S"; lig->name += hos.str();
            ++split_num;
            main_structure->splitted_ligands.push_back(lig);
        }
    }
    return mol_at.size();
}


//==============================================================================================
//Definitionen fuer AACID:
//==============================================================================================

AACID::AACID() : last_res(false) {
//    is_cav = false;
}

AACID::~AACID() {}

//--------------------------------------------------------------------------
//Methoden fuer AACID:

void AACID::get_intern_types() { //Have fun!
    istringstream is;
    string name;
    for (atoms_vec it=atoms.begin(); it!=atoms.end();++it) {
        is.clear(); //streamobjekt zurcksetzen
        is.str((*it)->name);
        is >> name;
        //Backbone:
        if (name == "N") {
            (*it)->element = "N";
            (*it)->intern_type = "N.ams";
            (*it)->dict_type = "BACKBONE | DICT | DIRECT";
            (*it)->bond_ind = 1;
            if ((*it)->res_name == "PRO") {
                (*it)->intern_type = "N.amt";
                (*it)->bond_ind = 60;
            }
            continue;
        }
        if (name == "H") {
            (*it)->element = "H";
            (*it)->intern_type = "H.onh";
            (*it)->dict_type = "BACKBONE | DICT | ESSENTIAL | DIRECT";
            (*it)->bond_ind = 80;
            continue;
        }
        if (name == "CA") {
            (*it)->element = "C";
            (*it)->intern_type = "C.3n";
            (*it)->dict_type = "BACKBONE | DICT | DIRECT";
            (*it)->bond_ind = 2;
            continue;
        }
        
        //!GN: 09.05.06   neben der HH12-Nomenklatur gibt es auch noch die 2HH1-Nomenklatur
        if (name[0] == '1' || name[0] == '2' || name[0] == '3') {
            if (name[1] == 'H') {
                string hbuf(name,1,string::npos);
                name = hbuf;
            }
        }
        
        if (name[0] == 'H') {
            if (name[1] == 'A') {
                (*it)->element = "H";
                (*it)->intern_type = "H.0";
                (*it)->dict_type = "BACKBONE | DICT";
                (*it)->bond_ind = 81;
                continue;
            }
            if (name[1] == 'O') {
                (*it)->element = "H";
                (*it)->intern_type = "H.ac";
                (*it)->dict_type = " ";
                (*it)->bond_ind = 101;
                continue;
            }
            if (name[1] == 'N') { //!GN: 09.05.06:  Terminaler N
                (*it)->element = "H";
                (*it)->intern_type = "H.n";
                (*it)->dict_type = "BACKBONE | DICT";
                (*it)->bond_ind = 80;
                //!Hier muss das entsprechende N auf N.3p setzen!!!
                //!Dafuer werden bonded_atoms gebraucht!!!
                            //!==========ToDo==========!//
                continue;
            }
        }
        if (name == "C") {
            (*it)->element = "C";
            (*it)->intern_type = "C.am";
            (*it)->dict_type = "BACKBONE | DICT | DIRECT";
            (*it)->bond_ind = 3;
            continue;
        }
        if (name == "O") {
            (*it)->element = "O";
            (*it)->intern_type = "O.am";
            (*it)->dict_type = "BACKBONE | DICT | DIRECT";
            (*it)->bond_ind = 4;
            continue;
        }
        if (name.size() > 1) {
            if (name[0] == 'O' && (name[1] == 'X' || name[1] == 'T')) { //terminaler Sauerstoff
                last_res = true;
                (*it)->element = "O";
                (*it)->intern_type = "O.co2"; //!wird von den meisten Programmen auf O.3 gesetzt
                (*it)->dict_type = "CAP | DICT";
                (*it)->bond_ind = 61;
                continue;
            }
        }
        
        (*it)->dict_type = "DICT"; //wenn kein Backbone
        //Alanin:
        if ((*it)->res_name == "ALA") {
            if (name == "CB") {(*it)->element = "C"; (*it)->intern_type = "C.3p"; (*it)->bond_ind = 5; continue;}
            if (name[0] == 'H') {(*it)->element = "H"; (*it)->intern_type = "H.0"; (*it)->bond_ind = 82; continue;}
        }
        
        //Arginin:
        if ((*it)->res_name == "ARG") {
            if (name == "CB") {(*it)->element = "C"; (*it)->intern_type = "C.3s"; (*it)->bond_ind = 5; continue;}
            if (name == "CG") {(*it)->element = "C"; (*it)->intern_type = "C.3s"; (*it)->bond_ind = 6; continue;}
            if (name == "CD") {(*it)->element = "C"; (*it)->intern_type = "C.3s"; (*it)->bond_ind = 8; continue;}
            if (name == "NE") {(*it)->element = "N"; (*it)->intern_type = "N.guh"; (*it)->bond_ind = 11; continue;}
            if (name == "NH1") {(*it)->element = "N"; (*it)->intern_type = "N.guh"; (*it)->bond_ind = 26; continue;}
            if (name == "NH2") {(*it)->element = "N"; (*it)->intern_type = "N.guh"; (*it)->bond_ind = 27; continue;}
            if (name == "CZ") {(*it)->element = "C"; (*it)->intern_type = "C.guh"; (*it)->bond_ind = 25; continue;}
            if (name[0] == 'H') {
                (*it)->element = "H";
                (*it)->intern_type = "H.0";
                if (name[1] == 'B') {(*it)->bond_ind = 82; continue;}
                if (name[1] == 'G') {(*it)->bond_ind = 83; continue;}
                if (name[1] == 'D') {(*it)->bond_ind = 85; continue;}
                if (name[1] == 'E') {(*it)->bond_ind = 99; (*it)->intern_type = "H.n"; continue;}
                if (name[1] == 'H') {
                    if (name[2] == '1') {(*it)->bond_ind = 97; (*it)->intern_type = "H.n"; continue;}
                    if (name[2] == '2') {(*it)->bond_ind = 102; (*it)->intern_type = "H.n"; continue;}
                }
            }
        }
        
        //Asparagin:
        if ((*it)->res_name == "ASN") {
            if (name == "CB") {(*it)->element = "C"; (*it)->intern_type = "C.3s"; (*it)->bond_ind = 5; continue;}
            if (name == "CG") {(*it)->element = "C"; (*it)->intern_type = "C.am"; (*it)->bond_ind = 6; continue;}
            if (name == "OD1") {(*it)->element = "O"; (*it)->intern_type = "O.am"; (*it)->bond_ind = 38; continue;}
            if (name == "ND2") {(*it)->element = "N"; (*it)->intern_type = "N.amp"; (*it)->bond_ind = 49; continue;}
            if (name[0] == 'H') {
                (*it)->element = "H";
                (*it)->intern_type = "H.0";
                if (name[1] == 'B') {(*it)->bond_ind = 82; continue;}
                if (name[1] == 'D') {(*it)->bond_ind = 93; (*it)->intern_type = "H.onh"; continue;}
            }
        }
        
        //Asparaginsaeure:
        if ((*it)->res_name == "ASP") {
            if (name == "CB") {(*it)->element = "C"; (*it)->intern_type = "C.3s"; (*it)->bond_ind = 5; continue;}
            if (name == "CG") {(*it)->element = "C"; (*it)->intern_type = "C.co2"; (*it)->bond_ind = 6; continue;}
            if (name == "OD1") {(*it)->element = "O"; (*it)->intern_type = "O.co2"; (*it)->bond_ind = 18; continue;}
            if (name == "OD2") {(*it)->element = "O"; (*it)->intern_type = "O.co2"; (*it)->bond_ind = 19; continue;}
            if (name[0] == 'H') {
                (*it)->element = "H";
                (*it)->intern_type = "H.0";
                if (name[1] == 'B') {(*it)->bond_ind = 82; continue;}
            }
        }
        
        //Asparaginsaeure neutral:
        if ((*it)->res_name == "ASH") {
            if (name == "CB") {(*it)->element = "C"; (*it)->intern_type = "C.3s"; (*it)->bond_ind = 5; continue;}
            if (name == "CG") {(*it)->element = "C"; (*it)->intern_type = "C.co2h"; (*it)->bond_ind = 6; continue;}
            if (name == "OD1") {(*it)->element = "O"; (*it)->intern_type = "O.2co2"; (*it)->bond_ind = 18; continue;}
            if (name == "OD2") {(*it)->element = "O"; (*it)->intern_type = "O.3ac"; (*it)->bond_ind = 19; continue;}
            if (name[0] == 'H') {
                (*it)->element = "H";
                (*it)->intern_type = "H.0";
                if (name[1] == 'B') {(*it)->bond_ind = 82; continue;}
            }
        }
        
        //unklar ob ASN oder ASP
        if ((*it)->res_name == "ASX") {
            if (name == "CB") {(*it)->element = "C"; (*it)->intern_type = "C.3s"; (*it)->bond_ind = 5; continue;}
            if (name == "CG") {(*it)->element = "C"; (*it)->intern_type = "C.2t"; (*it)->bond_ind = 6; continue;}
            if (name == "OD1") {
                (*it)->element = "O";
                (*it)->intern_type = "O.carb"; //akzeptable Wahl fuer ASN und ASP
                (*it)->bond_ind = 18;
                continue;
            }
            if (name[0] == 'H') {
                (*it)->element = "H";
                (*it)->intern_type = "H.0";
                if (name[1] == 'B') {(*it)->bond_ind = 82; continue;}
            }
            if (chain->protein->main_structure->verb > 0) {
                cerr << c_message<cWARNING>("AACID::get_sybyl_types --> ") << name << "  from  " 
                     << (*it)->res_name << " " << (*it)->res_number << " remain unconsidered" << endl;
            }
        }
        
        //Cystein:
        if ((*it)->res_name == "CYS") {
            if (name == "CB") {(*it)->element = "C"; (*it)->intern_type = "C.3s"; (*it)->bond_ind = 5; continue;}
            if (name == "SG") {(*it)->element = "S"; (*it)->intern_type = "S.sh"; (*it)->bond_ind = 6; continue;}
            if (name[0] == 'H') {
                if (name[1] == 'B') {(*it)->element = "H"; (*it)->bond_ind = 82; (*it)->intern_type = "H.0"; continue;}
                if (name[1] == 'G') {(*it)->element = "H"; (*it)->bond_ind = 83; (*it)->intern_type = "H.ac"; continue;}
            }
        }
        
        //Cystein in Disulfidbruecke:
        if ((*it)->res_name == "CYX") {
            if (name == "CB") {(*it)->element = "C"; (*it)->intern_type = "C.3s"; (*it)->bond_ind = 5; continue;}
            if (name == "SG") {(*it)->element = "S"; (*it)->intern_type = "S.s"; (*it)->bond_ind = 6; continue;}
            if (name[0] == 'H') {
                if (name[1] == 'B') {(*it)->element = "H"; (*it)->bond_ind = 82; (*it)->intern_type = "H.0"; continue;}
                if (name[1] == 'G') {(*it)->element = "H"; (*it)->bond_ind = 83; (*it)->intern_type = "H.ac"; continue;} //!duerfte nicht vorkommen
            }
        }
        
        //Cystein deprotoniert:
        if ((*it)->res_name == "CYM") {
            if (name == "CB") {(*it)->element = "C"; (*it)->intern_type = "C.3s"; (*it)->bond_ind = 5; continue;}
            if (name == "SG") {(*it)->element = "S"; (*it)->intern_type = "S.3"; (*it)->bond_ind = 6; continue;}
            if (name[0] == 'H') {
                if (name[1] == 'B') {(*it)->element = "H"; (*it)->bond_ind = 82; (*it)->intern_type = "H.0"; continue;}
                if (name[1] == 'G') {(*it)->element = "H"; (*it)->bond_ind = 83; (*it)->intern_type = "H.ac"; continue;} //!duerfte nicht vorkommen
            }
        }
        
        //Glutamin:
        if ((*it)->res_name == "GLN") {
            if (name == "CB") {(*it)->element = "C"; (*it)->intern_type = "C.3s"; (*it)->bond_ind = 5; continue;}
            if (name == "CG") {(*it)->element = "C"; (*it)->intern_type = "C.3s"; (*it)->bond_ind = 6; continue;}
            if (name == "CD") {(*it)->element = "C"; (*it)->intern_type = "C.am"; (*it)->bond_ind = 8; continue;}
            if (name == "OE1") {(*it)->element = "O"; (*it)->intern_type = "O.am"; (*it)->bond_ind = 40; continue;}
            if (name == "NE2") {(*it)->element = "N"; (*it)->intern_type = "N.amp"; (*it)->bond_ind = 51; continue;}
            if (name[0] == 'H') {
                (*it)->element = "H";
                (*it)->intern_type = "H.0";
                if (name[1] == 'B') {(*it)->bond_ind = 82; continue;}
                if (name[1] == 'G') {(*it)->bond_ind = 83; continue;}
                if (name[1] == 'E') {(*it)->bond_ind = 94; (*it)->intern_type = "H.onh"; continue;}
            }
        }
        
        //Glutaminsaeure:
        if ((*it)->res_name == "GLU") {
            if (name == "CB") {(*it)->element = "C"; (*it)->intern_type = "C.3s"; (*it)->bond_ind = 5; continue;}
            if (name == "CG") {(*it)->element = "C"; (*it)->intern_type = "C.3s"; (*it)->bond_ind = 6; continue;}
            if (name == "CD") {(*it)->element = "C"; (*it)->intern_type = "C.co2"; (*it)->bond_ind = 8; continue;}
            if (name == "OE1") {(*it)->element = "O"; (*it)->intern_type = "O.co2"; (*it)->bond_ind = 20; continue;}
            if (name == "OE2") {(*it)->element = "O"; (*it)->intern_type = "O.co2"; (*it)->bond_ind = 21; continue;}
            if (name[0] == 'H') {
                (*it)->element = "H";
                (*it)->intern_type = "H.0";
                if (name[1] == 'B') {(*it)->bond_ind = 82; continue;}
                if (name[1] == 'G') {(*it)->bond_ind = 83; continue;}
            }
        }
        
        //Glutaminsaeure neutral:
        if ((*it)->res_name == "GLH") {
            if (name == "CB") {(*it)->element = "C"; (*it)->intern_type = "C.3s"; (*it)->bond_ind = 5; continue;}
            if (name == "CG") {(*it)->element = "C"; (*it)->intern_type = "C.3s"; (*it)->bond_ind = 6; continue;}
            if (name == "CD") {(*it)->element = "C"; (*it)->intern_type = "C.co2h"; (*it)->bond_ind = 8; continue;}
            if (name == "OE1") {(*it)->element = "O"; (*it)->intern_type = "O.2co2"; (*it)->bond_ind = 20; continue;}
            if (name == "OE2") {(*it)->element = "O"; (*it)->intern_type = "O.3ac"; (*it)->bond_ind = 21; continue;}
            if (name[0] == 'H') {
                (*it)->element = "H";
                (*it)->intern_type = "H.0";
                if (name[1] == 'B') {(*it)->bond_ind = 82; continue;}
                if (name[1] == 'G') {(*it)->bond_ind = 83; continue;}
            }
        }
        
        //unklar ob GLN oder GLU:
        if ((*it)->res_name == "GLX") {
            if (name == "CB") {(*it)->element = "C"; (*it)->intern_type = "C.3s"; (*it)->bond_ind = 5; continue;}
            if (name == "CG") {(*it)->element = "C"; (*it)->intern_type = "C.3s"; (*it)->bond_ind = 6; continue;}
            if (name == "CD") {(*it)->element = "C"; (*it)->intern_type = "C.2t"; (*it)->bond_ind = 8; continue;}
            if (name == "OE1" || name == "OE2") {
                (*it)->element = "O";
                (*it)->intern_type = "O.carb";
                (*it)->bond_ind = 20;
                continue;
            }
            if (name[0] == 'H') {
                (*it)->element = "H";
                (*it)->intern_type = "H.0";
                if (name[1] == 'B') {(*it)->bond_ind = 82; continue;}
                if (name[1] == 'G') {(*it)->bond_ind = 83; continue;}
            }
            if (chain->protein->main_structure->verb > 0) {
                cerr << c_message<cWARNING>("AACID::get_sybyl_types --> ") << name << "  from  " 
                     << (*it)->res_name << " " << (*it)->res_number << " remains unconsidered" << endl;
            }
        }
        
        //Isoleucin:
        if ((*it)->res_name == "ILE") {
            if (name == "CB") {(*it)->element = "C"; (*it)->intern_type = "C.3t"; (*it)->bond_ind = 5; continue;}
            if (name == "CG1") {(*it)->element = "C"; (*it)->intern_type = "C.3s"; (*it)->bond_ind = 6; continue;}
            if (name == "CG2") {(*it)->element = "C"; (*it)->intern_type = "C.3p"; (*it)->bond_ind = 7; continue;}
            //GN: 04.07.06 :
            if (name == "CD1" || name == "CD") {(*it)->element = "C"; (*it)->intern_type = "C.3p"; (*it)->bond_ind = 8; continue;}
            if (name[0] == 'H') {
                (*it)->element = "H";
                (*it)->intern_type = "H.0";
                if (name[1] == 'B') {(*it)->bond_ind = 82; continue;}
                if (name[1] == 'G') {
                    if (name[2] == '1') {(*it)->bond_ind = 83; continue;}
                    if (name[2] == '2') {(*it)->bond_ind = 84; continue;}
                }
                if (name[1] == 'D') {(*it)->bond_ind = 85; continue;}
            }
        }
        //!====================TEST==================================
        //Leucin:
        if ((*it)->res_name == "LEU") {
            if (name == "CB") {(*it)->element = "C"; (*it)->intern_type = "C.3s"; (*it)->bond_ind = 5; continue;}
            if (name == "CG") {(*it)->element = "C"; (*it)->intern_type = "C.3t"; (*it)->bond_ind = 6; continue;}
            if (name == "CD1") {(*it)->element = "C"; (*it)->intern_type = "C.3p"; (*it)->bond_ind = 8; continue;}
            if (name == "CD2") {(*it)->element = "C"; (*it)->intern_type = "C.3p"; (*it)->bond_ind = 9; continue;}
            if (name[0] == 'H') {
                (*it)->element = "H"; (*it)->intern_type = "H.0";
                if (name[1] == 'B') {(*it)->bond_ind = 82; continue;}
                if (name[1] == 'D') {
                    if (name[2] == '1') {(*it)->bond_ind = 85; continue;}
                    if (name[2] == '2') {(*it)->bond_ind = 86; continue;}
                }
                if (name[1] == 'G') {(*it)->bond_ind = 83; continue;}
            }
        }
        
        //Lysin:
        if ((*it)->res_name == "LYS") {
            if (name == "CB") {(*it)->element = "C"; (*it)->intern_type = "C.3s"; (*it)->bond_ind = 5; continue;}
            if (name == "CG") {(*it)->element = "C"; (*it)->intern_type = "C.3s"; (*it)->bond_ind = 6; continue;}
            if (name == "CD") {(*it)->element = "C"; (*it)->intern_type = "C.3s"; (*it)->bond_ind = 8; continue;}
            if (name == "CE") {(*it)->element = "C"; (*it)->intern_type = "C.3n"; (*it)->bond_ind = 11; continue;}
            if (name == "NZ") {(*it)->element = "N"; (*it)->intern_type = "N.4h"; (*it)->bond_ind = 13; continue;}
            if (name[0] == 'H') {
                (*it)->element = "H"; (*it)->intern_type = "H.0";
                if (name[1] == 'B') {(*it)->bond_ind = 82; continue;}
                if (name[1] == 'G') {(*it)->bond_ind = 83; continue;}
                if (name[1] == 'D') {(*it)->bond_ind = 85; continue;}
                if (name[1] == 'E') {(*it)->bond_ind = 99; continue;}
                if (name[1] == 'Z') {(*it)->bond_ind = 92; (*it)->intern_type = "H.n"; continue;}
            }
        }
        
        //Lysin neutral:
        if ((*it)->res_name == "LYN") {
            if (name == "CB") {(*it)->element = "C"; (*it)->intern_type = "C.3s"; (*it)->bond_ind = 5; continue;}
            if (name == "CG") {(*it)->element = "C"; (*it)->intern_type = "C.3s"; (*it)->bond_ind = 6; continue;}
            if (name == "CD") {(*it)->element = "C"; (*it)->intern_type = "C.3s"; (*it)->bond_ind = 8; continue;}
            if (name == "CE") {(*it)->element = "C"; (*it)->intern_type = "C.3n"; (*it)->bond_ind = 11; continue;}
            if (name == "NZ") {(*it)->element = "N"; (*it)->intern_type = "N.3p"; (*it)->bond_ind = 13; continue;}
            if (name[0] == 'H') {
                (*it)->element = "H"; (*it)->intern_type = "H.0";
                if (name[1] == 'B') {(*it)->bond_ind = 82; continue;}
                if (name[1] == 'G') {(*it)->bond_ind = 83; continue;}
                if (name[1] == 'D') {(*it)->bond_ind = 85; continue;}
                if (name[1] == 'E') {(*it)->bond_ind = 99; continue;}
                if (name[1] == 'Z') {(*it)->bond_ind = 92; (*it)->intern_type = "H.n"; continue;}
            }
        }
        
        //Methionin:
        if ((*it)->res_name == "MET") {
            if (name == "CB") {(*it)->element = "C"; (*it)->intern_type = "C.3s"; (*it)->bond_ind = 5; continue;}
            if (name == "CG") {(*it)->element = "C"; (*it)->intern_type = "C.3s"; (*it)->bond_ind = 6; continue;}
            if (name == "SD") {(*it)->element = "S"; (*it)->intern_type = "S.3"; (*it)->bond_ind = 8; continue;}
            if (name == "CE") {(*it)->element = "C"; (*it)->intern_type = "C.3p"; (*it)->bond_ind = 11; continue;}
            if (name[0] == 'H') {
                (*it)->element = "H"; (*it)->intern_type = "H.0";
                if (name[1] == 'B') {(*it)->bond_ind = 82; continue;}
                if (name[1] == 'G') {(*it)->bond_ind = 83; continue;}
                if (name[1] == 'E') {(*it)->bond_ind = 99; continue;}
            }
        }
        
        //Phenylalanin:
        if ((*it)->res_name == "PHE") {
            if (name == "CB") {(*it)->element = "C"; (*it)->intern_type = "C.3s"; (*it)->bond_ind = 5; continue;}
            if (name == "CG") {(*it)->element = "C"; (*it)->intern_type = "C.ar6"; (*it)->bond_ind = 6; continue;}
            if (name == "CD1") {(*it)->element = "C"; (*it)->intern_type = "C.ar6"; (*it)->bond_ind = 18; continue;}
            if (name == "CD2") {(*it)->element = "C"; (*it)->intern_type = "C.ar6"; (*it)->bond_ind = 19; continue;}
            if (name == "CE1") {(*it)->element = "C"; (*it)->intern_type = "C.ar6"; (*it)->bond_ind = 20; continue;}
            if (name == "CE2") {(*it)->element = "C"; (*it)->intern_type = "C.ar6"; (*it)->bond_ind = 21; continue;}
            if (name == "CZ") {(*it)->element = "C"; (*it)->intern_type = "C.ar6"; (*it)->bond_ind = 22; continue;}
            if (name[0] == 'H') {
                (*it)->element = "H"; (*it)->intern_type = "H.0";
                if (name[1] == 'B') {(*it)->bond_ind = 82; continue;}
                if (name[1] == 'D') {
                    if (name[2] == '1') {(*it)->bond_ind = 87; continue;}
                    if (name[2] == '2') {(*it)->bond_ind = 88; continue;}
                }
                if (name[1] == 'E') {
                    if (name[2] == '1') {(*it)->bond_ind = 89; continue;}
                    if (name[2] == '2') {(*it)->bond_ind = 90; continue;}
                }
                if (name[1] == 'Z') {(*it)->bond_ind = 91; continue;}
            }
        }
        
        //Prolin:
        if ((*it)->res_name == "PRO") {
            if (name == "CB") {(*it)->element = "C"; (*it)->intern_type = "C.3s"; (*it)->bond_ind = 5; continue;}
            if (name == "CG") {(*it)->element = "C"; (*it)->intern_type = "C.3s"; (*it)->bond_ind = 6; continue;}
            if (name == "CD") {(*it)->element = "C"; (*it)->intern_type = "C.3s"; (*it)->bond_ind = 8; continue;}
            if (name[0] == 'H') {
                (*it)->element = "H"; (*it)->intern_type = "H.0";
                if (name[1] == 'B') {(*it)->bond_ind = 82; continue;}
                if (name[1] == 'G') {(*it)->bond_ind = 83; continue;}
                if (name[1] == 'D') {(*it)->bond_ind = 85; continue;}
            }
        }
        
        //Serin:
        if ((*it)->res_name == "SER") {
            if (name == "CB") {(*it)->element = "C"; (*it)->intern_type = "C.ohp"; (*it)->bond_ind = 5; continue;}
            if (name == "OG") {(*it)->element = "O"; (*it)->intern_type = "O.3oh"; (*it)->bond_ind = 6; continue;}
            if (name[0] == 'H') {
                if (name[1] == 'B') {(*it)->bond_ind = 82; (*it)->element = "H"; (*it)->intern_type = "H.0"; continue;}
                if (name[1] == 'G') {(*it)->bond_ind = 83; (*it)->element = "H"; (*it)->intern_type = "H.o"; continue;}
            }
        }
        
        //Threonin:
        if ((*it)->res_name == "THR") {
            if (name == "CB") {(*it)->element = "C"; (*it)->intern_type = "C.ohs"; (*it)->bond_ind = 5; continue;}
            if (name == "OG1") {(*it)->element = "O"; (*it)->intern_type = "O.3oh"; (*it)->bond_ind = 6; continue;}
            if (name == "CG2") {(*it)->element = "C"; (*it)->intern_type = "C.3p"; (*it)->bond_ind = 7; continue;}
            if (name[0] == 'H') {
                (*it)->element = "H"; (*it)->intern_type = "H.0";
                if (name[1] == 'B') {(*it)->bond_ind = 82; continue;}
                if (name[1] == 'G') {
                    if (name[2] == '1') {(*it)->bond_ind = 83; (*it)->intern_type = "H.o"; continue;}
                    if (name[2] == '2') {(*it)->bond_ind = 84; continue;}
                }
            }
        }
        
        //Tryptophan:
        if ((*it)->res_name == "TRP") {
            if (name == "CB") {(*it)->element = "C"; (*it)->intern_type = "C.3s"; (*it)->bond_ind = 5; continue;}
            if (name == "CG") {(*it)->element = "C"; (*it)->intern_type = "C.arx"; (*it)->bond_ind = 6; continue;}
            if (name == "CD1") {(*it)->element = "C"; (*it)->intern_type = "C.arx"; (*it)->bond_ind = 38; continue;}
            if (name == "CD2") {(*it)->element = "C"; (*it)->intern_type = "C.ar6"; (*it)->bond_ind = 9; continue;}
            if (name == "NE1") {(*it)->element = "N"; (*it)->intern_type = "N.ar3h"; (*it)->bond_ind = 10; continue;}
            if (name == "CE2") {(*it)->element = "C"; (*it)->intern_type = "C.ar6"; (*it)->bond_ind = 21; continue;}
            if (name == "CE3") {(*it)->element = "C"; (*it)->intern_type = "C.ar6"; (*it)->bond_ind = 25; continue;}
            if (name == "CZ2") {(*it)->element = "C"; (*it)->intern_type = "C.ar6"; (*it)->bond_ind = 22; continue;}
            if (name == "CZ3") {(*it)->element = "C"; (*it)->intern_type = "C.ar6"; (*it)->bond_ind = 26; continue;}
            if (name == "CH2") {(*it)->element = "C"; (*it)->intern_type = "C.ar6"; (*it)->bond_ind = 20; continue;}
            if (name[0] == 'H') {
                (*it)->element = "H"; (*it)->intern_type = "H.0";
                if (name[1] == 'B') {(*it)->bond_ind = 82; continue;}
                if (name[1] == 'Z') {
                    if (name[2] == '2') {(*it)->bond_ind = 91; continue;}
                    if (name[2] == '3') {(*it)->bond_ind = 97; continue;}
                }
                if (name[1] == 'E') {
                    if (name[2] == '1') {(*it)->bond_ind = 96; (*it)->intern_type = "H.n"; continue;}
                    if (name[2] == '3') {(*it)->bond_ind = 98; continue;}
                }
                if (name[1] == 'D') {(*it)->bond_ind = 95; continue;}
                if (name[1] == 'H') {(*it)->bond_ind = 89; continue;}
            }
        }
        
        //Tyrosin:
        if ((*it)->res_name == "TYR") {
            if (name == "CB") {(*it)->element = "C"; (*it)->intern_type = "C.3s"; (*it)->bond_ind = 5; continue;}
            if (name == "CG") {(*it)->element = "C"; (*it)->intern_type = "C.ar6"; (*it)->bond_ind = 6; continue;}
            if (name == "CD1") {(*it)->element = "C"; (*it)->intern_type = "C.ar6"; (*it)->bond_ind = 18; continue;}
            if (name == "CD2") {(*it)->element = "C"; (*it)->intern_type = "C.ar6"; (*it)->bond_ind = 19; continue;}
            if (name == "CE1") {(*it)->element = "C"; (*it)->intern_type = "C.ar6"; (*it)->bond_ind = 20; continue;}
            if (name == "CE2") {(*it)->element = "C"; (*it)->intern_type = "C.ar6"; (*it)->bond_ind = 21; continue;}
            if (name == "CZ") {(*it)->element = "C"; (*it)->intern_type = "C.ar6"; (*it)->bond_ind = 22; continue;}
            if (name == "OH") {(*it)->element = "O"; (*it)->intern_type = "O.ph"; (*it)->bond_ind = 13; continue;}
            if (name[0] == 'H') {
                (*it)->element = "H"; (*it)->intern_type = "H.0";
                if (name[1] == 'B') {(*it)->bond_ind = 82; continue;}
                if (name[1] == 'D') {
                    if (name[2] == '1') {(*it)->bond_ind = 87; continue;}
                    if (name[2] == '2') {(*it)->bond_ind = 88; continue;}
                }
                if (name[1] == 'E') {
                    if (name[2] == '1') {(*it)->bond_ind = 89; continue;}
                    if (name[2] == '2') {(*it)->bond_ind = 90; continue;}
                }
                if (name[1] == 'H') {(*it)->bond_ind = 92; (*it)->intern_type = "H.o"; continue;}
            }
        }
        
        //Tyrosin deprotoniert:
        if ((*it)->res_name == "TYM") {
            if (name == "CB") {(*it)->element = "C"; (*it)->intern_type = "C.3s"; (*it)->bond_ind = 5; continue;}
            if (name == "CG") {(*it)->element = "C"; (*it)->intern_type = "C.ar6"; (*it)->bond_ind = 6; continue;}
            if (name == "CD1") {(*it)->element = "C"; (*it)->intern_type = "C.ar6"; (*it)->bond_ind = 18; continue;}
            if (name == "CD2") {(*it)->element = "C"; (*it)->intern_type = "C.ar6"; (*it)->bond_ind = 19; continue;}
            if (name == "CE1") {(*it)->element = "C"; (*it)->intern_type = "C.ar6"; (*it)->bond_ind = 20; continue;}
            if (name == "CE2") {(*it)->element = "C"; (*it)->intern_type = "C.ar6"; (*it)->bond_ind = 21; continue;}
            if (name == "CZ") {(*it)->element = "C"; (*it)->intern_type = "C.ar6"; (*it)->bond_ind = 22; continue;}
            if (name == "OH") {(*it)->element = "O"; (*it)->intern_type = "O.ph"; (*it)->bond_ind = 13; continue;}
            if (name[0] == 'H') {
                (*it)->element = "H"; (*it)->intern_type = "H.0";
                if (name[1] == 'B') {(*it)->bond_ind = 82; continue;}
                if (name[1] == 'D') {
                    if (name[2] == '1') {(*it)->bond_ind = 87; continue;}
                    if (name[2] == '2') {(*it)->bond_ind = 88; continue;}
                }
                if (name[1] == 'E') {
                    if (name[2] == '1') {(*it)->bond_ind = 89; continue;}
                    if (name[2] == '2') {(*it)->bond_ind = 90; continue;}
                }
                if (name[1] == 'H') {(*it)->bond_ind = 92; (*it)->intern_type = "H.o"; continue;} //!duerfte es hier nicht geben
            }
        }
        
        //Valin:
        if ((*it)->res_name == "VAL") {
            if (name == "CB") {(*it)->element = "C"; (*it)->intern_type = "C.3t"; (*it)->bond_ind = 5; continue;}
            if (name == "CG1") {(*it)->element = "C"; (*it)->intern_type = "C.3p"; (*it)->bond_ind = 6; continue;}
            if (name == "CG2") {(*it)->element = "C"; (*it)->intern_type = "C.3p"; (*it)->bond_ind = 7; continue;}
            if (name[0] == 'H') {
                (*it)->element = "H"; (*it)->intern_type = "H.0";
                if (name[1] == 'B') {(*it)->bond_ind = 82; continue;}
                if (name[1] == 'G') {
                    if (name[2] == '1') {(*it)->bond_ind = 83; continue;}
                    if (name[2] == '2') {(*it)->bond_ind = 84; continue;}
                }
            }
        }
        
        //Histidin:
        if ((*it)->res_name == "HIS" || (*it)->res_name == "HID" || (*it)->res_name == "HIE") { //Histidin Tautomere
            //! Problem: sowohl ND1, als auch NE2 kann das N.pl3 sein
            //! => wenn HE, oder HD vorkommt entsprechend setzen, sonst die default-Variante
            //! mit Proton an ND
            if (name == "CB") {(*it)->element = "C"; (*it)->intern_type = "C.3s"; (*it)->bond_ind = 5; continue;}
            if (name == "CG") {(*it)->element = "C"; (*it)->intern_type = "C.arx"; (*it)->bond_ind = 6; continue;}
            if (name == "CD2") {(*it)->element = "C"; (*it)->intern_type = "C.arx"; (*it)->bond_ind = 38; continue;}
            if (name == "ND1" && (*it)->intern_type != "N.ar2") {(*it)->element = "N"; (*it)->intern_type = "N.ar3h"; (*it)->bond_ind = 8; continue;} 
            if (name == "CE1") {(*it)->element = "C"; (*it)->intern_type = "C.arx"; (*it)->bond_ind = 11; continue;}
            if (name == "NE2" && (*it)->intern_type != "N.ar3h") {(*it)->element = "N"; (*it)->intern_type = "N.ar2"; (*it)->bond_ind = 12; continue;}
            if (name[0] == 'H') {
                (*it)->element = "H"; (*it)->intern_type = "H.0";
                if (name[1] == 'B') {(*it)->bond_ind = 82; continue;}
                if (name[1] == 'D') {
                    if (name[2] == '1') {(*it)->bond_ind = 85; (*it)->intern_type = "H.n"; continue;} //zu ND1
                    if (name[2] == '2') {(*it)->bond_ind = 95; continue;} //zu CD2
                }
                if (name[1] == 'E') {
                    if (name[2] == '1') {(*it)->bond_ind = 99; continue;} //zu CE1
                    if (name[2] == '2') {
                        (*it)->intern_type = "H.n";
                        (*it)->bond_ind = 100;
                        //Proton ist also an NE2 => Typ von ND1 und NE2 entsprechend setzen
                        istringstream is2;
                        string name2;
                        for (atoms_vec it2=atoms.begin(); it2!=atoms.end();++it2) {
                            is2.clear(); //streamobjekt zurcksetzen
                            is2.str((*it2)->name);
                            is2 >> name2;
                            if (name2 == "ND1") {(*it2)->element = "N"; (*it2)->intern_type = "N.ar2"; (*it2)->bond_ind = 8; continue;}
                            if (name2 == "NE2") {(*it2)->element = "N"; (*it2)->intern_type = "N.ar3h"; (*it2)->bond_ind = 12; continue;}
                        }
                        continue;
                    }
                }
            }
        }
        
        //Histidin protoniert
        if ((*it)->res_name == "HIP") {
            if (name == "CB") {(*it)->element = "C"; (*it)->intern_type = "C.3s"; (*it)->bond_ind = 5; continue;}
            if (name == "CG") {(*it)->element = "C"; (*it)->intern_type = "C.arx"; (*it)->bond_ind = 6; continue;}
            if (name == "CD2") {(*it)->element = "C"; (*it)->intern_type = "C.arx"; (*it)->bond_ind = 38; continue;}
            if (name == "ND1") {(*it)->element = "N"; (*it)->intern_type = "N.arp"; (*it)->bond_ind = 8; continue;} 
            if (name == "CE1") {(*it)->element = "C"; (*it)->intern_type = "C.arp"; (*it)->bond_ind = 11; continue;}
            if (name == "NE2") {(*it)->element = "N"; (*it)->intern_type = "N.arp"; (*it)->bond_ind = 12; continue;}
            if (name[0] == 'H') {
                (*it)->element = "H"; (*it)->intern_type = "H.0";
                if (name[1] == 'B') {(*it)->bond_ind = 82; continue;}
                if (name[1] == 'D') {
                    if (name[2] == '1') {(*it)->bond_ind = 85; (*it)->intern_type = "H.n"; continue;}
                    if (name[2] == '2') {(*it)->bond_ind = 95; continue;}
                }
                if (name[1] == 'E') {
                    if (name[2] == '1') {(*it)->bond_ind = 99; continue;}
                    if (name[2] == '2') {(*it)->bond_ind = 100; (*it)->intern_type = "H.n"; continue;}
                }
            }
        }
        
        //! 07.10.2009: modifizierte Aminosaeuren: =================
        //! ACHTUNG: Diese AS muessen als HETATM deklariert sein!!!
        //Selenomethionin
        if ((*it)->res_name == "MSE") {
            if (name == "CB") {(*it)->element = "C"; (*it)->intern_type = "C.3s"; (*it)->bond_ind = 5; continue;}
            if (name == "CG") {(*it)->element = "C"; (*it)->intern_type = "C.3s"; (*it)->bond_ind = 6; continue;}
            if (name == "SE") {(*it)->element = "Se"; (*it)->intern_type = "Se"; (*it)->bond_ind = 8; continue;}
            if (name == "CE") {(*it)->element = "C"; (*it)->intern_type = "C.3p"; (*it)->bond_ind = 11; continue;}
            if (name[0] == 'H') {
                (*it)->element = "H"; (*it)->intern_type = "H.0";
                if (name[1] == 'B') {(*it)->bond_ind = 82; continue;}
                if (name[1] == 'G') {(*it)->bond_ind = 83; continue;}
                if (name[1] == 'E') {(*it)->bond_ind = 99; continue;}
            }
        }
        
        //3-Sulfinoalanin
        if ((*it)->res_name == "CSD") {
            if (name == "CB") {(*it)->element = "C"; (*it)->intern_type = "C.3s"; (*it)->bond_ind = 5; continue;}
            if (name == "SG") {(*it)->element = "S"; (*it)->intern_type = "S.o2"; (*it)->bond_ind = 6; continue;}
            if (name == "OD1") {(*it)->element = "O"; (*it)->intern_type = "O.2s"; (*it)->bond_ind = 18; continue;}
            if (name == "OD2") {(*it)->element = "O"; (*it)->intern_type = "O.2s"; (*it)->bond_ind = 19; continue;}
            if (name[0] == 'H') {
                (*it)->element = "H"; (*it)->intern_type = "H.0";
                if (name[1] == 'B') {(*it)->bond_ind = 82; continue;}
            }
        }
        
        //N-Dimethyl-Lysin
        if ((*it)->res_name == "LYN") {
            if (name == "CB") {(*it)->element = "C"; (*it)->intern_type = "C.3s"; (*it)->bond_ind = 5; continue;}
            if (name == "CG") {(*it)->element = "C"; (*it)->intern_type = "C.3s"; (*it)->bond_ind = 6; continue;}
            if (name == "CD") {(*it)->element = "C"; (*it)->intern_type = "C.3s"; (*it)->bond_ind = 8; continue;}
            if (name == "CE") {(*it)->element = "C"; (*it)->intern_type = "C.3n"; (*it)->bond_ind = 11; continue;}
            if (name == "NZ") {(*it)->element = "N"; (*it)->intern_type = "N.3t"; (*it)->bond_ind = 13; continue;}
            if (name == "CH1") {(*it)->element = "C"; (*it)->intern_type = "C.3n"; (*it)->bond_ind = 92; continue;}
            if (name == "CH2") {(*it)->element = "C"; (*it)->intern_type = "C.3n"; (*it)->bond_ind = 92; continue;}
            if (name[0] == 'H') {
                (*it)->element = "H"; (*it)->intern_type = "H.0";
                if (name[1] == 'B') {(*it)->bond_ind = 82; continue;}
                if (name[1] == 'G') {(*it)->bond_ind = 83; continue;}
                if (name[1] == 'D') {(*it)->bond_ind = 85; continue;}
                if (name[1] == 'E') {(*it)->bond_ind = 99; continue;}
                //! Hier fehlen noch die Wasserstoffe an den beiden Methylgruppen !!!!!
            }
        }
        
        //N-Mehtyl-Asparagin
        if ((*it)->res_name == "MEN") {
            if (name == "CB") {(*it)->element = "C"; (*it)->intern_type = "C.3s"; (*it)->bond_ind = 5; continue;}
            if (name == "CG") {(*it)->element = "C"; (*it)->intern_type = "C.am"; (*it)->bond_ind = 6; continue;}
            if (name == "OD1") {(*it)->element = "O"; (*it)->intern_type = "O.am"; (*it)->bond_ind = 38; continue;}
            if (name == "ND2") {(*it)->element = "N"; (*it)->intern_type = "N.ams"; (*it)->bond_ind = 49; continue;}
            if (name == "CE2") {(*it)->element = "C"; (*it)->intern_type = "C.3p"; (*it)->bond_ind = 93; continue;}
            if (name[0] == 'H') {
                (*it)->element = "H";
                (*it)->intern_type = "H.0";
                if (name[1] == 'B') {(*it)->bond_ind = 82; continue;}
                if (name[1] == 'D') {(*it)->bond_ind = 93; (*it)->intern_type = "H.onh"; continue;}
                //! Hier fehlen noch die Wasserstoffe an der Methylgruppe !!!!!
            }
        }
        
        //! Hier fehlen noch viele wichtige modifizierte AS, insbesondere PCA, PTR, SEP,
        //! TPO, CME, CGU, ABA
        //! -> siehe 'http://deposit.pdb.org/het_dictionary.txt'
        //!=========================================================
        
        
        //!neu: 06.02.2007:
        //Terminaler Stickstoff hat noch extra H-Atome:
        if (name[0] == 'H') {
            (*it)->element = "H"; (*it)->intern_type = "H.n";
            (*it)->dict_type = "BACKBONE | DICT | ESSENTIAL | DIRECT";
            (*it)->bond_ind = 80;
            continue;
        }
        
        if (chain->protein->main_structure->verb > 0) {
            cerr << c_message<cWARNING>("AACID::get_sybyl_types --> could not process  ") << name << "  from  " 
                 << (*it)->res_name << " " << (*it)->res_number << endl;
        }
    }
}

void AACID::create_bond(stl_ptr<ATOM> &from, stl_ptr<ATOM> &to, string type, string &d_type) {
    stl_ptr<BOND> bnd(new BOND());
    chain->bonds.push_back(bnd);
    bnd->from = from;
    bnd->to = to;
    bnd->type = type;
    bnd->id = chain->protein->bnd_id;
    chain->protein->bnd_id += 1;
    bnd->dict_type = d_type;
}

void AACID::build_bonds() { //Have fun again!
    string dt;
    
    if (!(chain->protein->last_C.zero())) { //Bindung zum vorherigen Rest
    //    if (!(chain->protein->last_C->is_ter)) { //!auf TER checken, weil auch bei fehlenden AS Ter-Eintraege
            for (atoms_vec it=atoms.begin(); it!=atoms.end();++it) {
                
                if ((*it)->alt_loc_id != ' ' && (*it)->alt_loc_id != 'A') continue;
                
                if (((*it)->res_number - chain->protein->last_C->res_number) != 1) { //!NEU: 16.04.08
                    if (get_square_distance(chain->protein->last_C->coord,(*it)->coord) > 2.5) continue; //!NEU: 16.04.08
                }
                if (((*it)->bond_ind == 1) || ((*it)->bond_ind == 60)) {
                    (*it)->bonded_atoms.push_back(chain->protein->last_C);
                    chain->protein->last_C->bonded_atoms.push_back(*it);
                    dt = "BACKBONE | DICT | INTERRES";
                    create_bond(chain->protein->last_C,*it,"am",dt);
                    break;
                }
            }
    //    }
    }
    
    for (atoms_vec it=atoms.begin(); it!=atoms.end();++it) {
        if ((*it)->bond_ind == 2) {
            for (atoms_vec jt=atoms.begin(); jt!=atoms.end();++jt) {
                if ((*it)->alt_loc_id != (*jt)->alt_loc_id) continue;
                if (((*jt)->bond_ind == 1) || ((*jt)->bond_ind == 3) || ((*jt)->bond_ind == 60)) {
                    (*it)->bonded_atoms.push_back(*jt);
                    (*jt)->bonded_atoms.push_back(*it);
                    dt = "BACKBONE | DICT";
                    create_bond(*it,*jt,"1",dt);
                }
                if ((*jt)->bond_ind == 5) {
                    (*it)->bonded_atoms.push_back(*jt);
                    (*jt)->bonded_atoms.push_back(*it);
                    dt = "DICT";
                    create_bond(*it,*jt,"1",dt);
                }
            }
        }
        
        if ((*it)->bond_ind == 3) {
            chain->protein->last_C = *it; //!rein: 16.04.08
            for (atoms_vec jt=atoms.begin(); jt!=atoms.end();++jt) {
                if ((*it)->alt_loc_id != (*jt)->alt_loc_id) continue;
                if ((*jt)->bond_ind == 4) {
                    (*it)->bonded_atoms.push_back(*jt);
                    (*jt)->bonded_atoms.push_back(*it);
                    dt = "BACKBONE | DICT";
                    create_bond(*it,*jt,"2",dt);
                //    chain->protein->last_C = *it; //! raus: 16.04.08
                }
                
                //!NEU: 15.04.2008 : Bindung zwischen C und OXT =========
                if ((*jt)->bond_ind == 61) {
                    (*it)->bonded_atoms.push_back(*jt);
                    (*jt)->bonded_atoms.push_back(*it);
                    dt = "BACKBONE | DICT";
                    create_bond(*it,*jt,"ar",dt);
                }
                //!======================================================
                /*raus: 16.04.08
                if ((*jt)->is_ter) {
                    chain->protein->last_C = NULL; //!Es gibt auch TER-Eintraege in einer Kette
                }
                */
            }
        }
        
        dt = "DICT";
        if ((*it)->bond_ind == 5) {
            for (atoms_vec jt=atoms.begin(); jt!=atoms.end();++jt) {
                if ((*it)->alt_loc_id != (*jt)->alt_loc_id) continue;
                if (((*jt)->bond_ind == 6) || ((*jt)->bond_ind == 7)) {
                    (*it)->bonded_atoms.push_back(*jt);
                    (*jt)->bonded_atoms.push_back(*it);
                    create_bond(*it,*jt,"1",dt);
                }
            }
        }
        
        if ((*it)->bond_ind == 6) {
            for (atoms_vec jt=atoms.begin(); jt!=atoms.end();++jt) {
                if ((*it)->alt_loc_id != (*jt)->alt_loc_id) continue;
                if (((*jt)->bond_ind == 8) || ((*jt)->bond_ind == 9)) {
                    (*it)->bonded_atoms.push_back(*jt);
                    (*jt)->bonded_atoms.push_back(*it);
                    create_bond(*it,*jt,"1",dt);
                }
                if (((*jt)->bond_ind == 18) || ((*jt)->bond_ind == 19)) {
                    (*it)->bonded_atoms.push_back(*jt);
                    (*jt)->bonded_atoms.push_back(*it);
                    create_bond(*it,*jt,"ar",dt);
                }
                if ((*jt)->bond_ind == 38) {
                    (*it)->bonded_atoms.push_back(*jt);
                    (*jt)->bonded_atoms.push_back(*it);
                    create_bond(*it,*jt,"2",dt);
                }
                if ((*jt)->bond_ind == 49) {
                    (*it)->bonded_atoms.push_back(*jt);
                    (*jt)->bonded_atoms.push_back(*it);
                    create_bond(*it,*jt,"am",dt);
                }
            }
        }
        
        if ((*it)->bond_ind == 8) {
            for (atoms_vec jt=atoms.begin(); jt!=atoms.end();++jt) {
                if ((*it)->alt_loc_id != (*jt)->alt_loc_id) continue;
                if (((*jt)->bond_ind == 20) || ((*jt)->bond_ind == 21)) {
                    (*it)->bonded_atoms.push_back(*jt);
                    (*jt)->bonded_atoms.push_back(*it);
                    create_bond(*it,*jt,"ar",dt);
                }
                if ((*jt)->bond_ind == 40) {
                    (*it)->bonded_atoms.push_back(*jt);
                    (*jt)->bonded_atoms.push_back(*it);
                    create_bond(*it,*jt,"2",dt);
                }
                if ((*jt)->bond_ind == 51) {
                    (*it)->bonded_atoms.push_back(*jt);
                    (*jt)->bonded_atoms.push_back(*it);
                    create_bond(*it,*jt,"am",dt);
                }
                if ((*jt)->bond_ind == 60) {
                    (*it)->bonded_atoms.push_back(*jt);
                    (*jt)->bonded_atoms.push_back(*it);
                    create_bond(*it,*jt,"1",dt);
                }
            }
        }
        
        if ((*it)->bond_ind == 11) {
            for (atoms_vec jt=atoms.begin(); jt!=atoms.end();++jt) {
                if ((*it)->alt_loc_id != (*jt)->alt_loc_id) continue;
                if (((*jt)->bond_ind == 8) || ((*jt)->bond_ind == 13)) {
                    (*it)->bonded_atoms.push_back(*jt);
                    (*jt)->bonded_atoms.push_back(*it);
                    create_bond(*it,*jt,"1",dt);
                }
                if ((*jt)->bond_ind == 12) {
                    (*it)->bonded_atoms.push_back(*jt);
                    (*jt)->bonded_atoms.push_back(*it);
                    create_bond(*it,*jt,"2",dt);
                }
                if ((*jt)->bond_ind == 25) {
                    (*it)->bonded_atoms.push_back(*jt);
                    (*jt)->bonded_atoms.push_back(*it);
                    create_bond(*it,*jt,"ar",dt);
                }
            }
        }
        
        if ((*it)->bond_ind == 13) {
            for (atoms_vec jt=atoms.begin(); jt!=atoms.end();++jt) {
                if ((*it)->alt_loc_id != (*jt)->alt_loc_id) continue;
                if ((*jt)->bond_ind == 22) {
                    (*it)->bonded_atoms.push_back(*jt);
                    (*jt)->bonded_atoms.push_back(*it);
                    create_bond(*it,*jt,"ar",dt);
                }
            }
        }
        
        if ((*it)->bond_ind == 20) {
            for (atoms_vec jt=atoms.begin(); jt!=atoms.end();++jt) {
                if ((*it)->alt_loc_id != (*jt)->alt_loc_id) continue;
                if (((*jt)->bond_ind == 18) || ((*jt)->bond_ind == 22) || ((*jt)->bond_ind == 26)) {
                    (*it)->bonded_atoms.push_back(*jt);
                    (*jt)->bonded_atoms.push_back(*it);
                    create_bond(*it,*jt,"ar",dt);
                }
            }
        }
        
        if ((*it)->bond_ind == 21) {
            for (atoms_vec jt=atoms.begin(); jt!=atoms.end();++jt) {
                if ((*it)->alt_loc_id != (*jt)->alt_loc_id) continue;
                if (((*jt)->bond_ind == 19) || ((*jt)->bond_ind == 22) || ((*jt)->bond_ind == 9)) {
                    (*it)->bonded_atoms.push_back(*jt);
                    (*jt)->bonded_atoms.push_back(*it);
                    create_bond(*it,*jt,"ar",dt);
                }
                if ((*jt)->bond_ind == 10) {
                    (*it)->bonded_atoms.push_back(*jt);
                    (*jt)->bonded_atoms.push_back(*it);
                    create_bond(*it,*jt,"1",dt);
                }
            }
        }
        
        if ((*it)->bond_ind == 25) {
            for (atoms_vec jt=atoms.begin(); jt!=atoms.end();++jt) {
                if ((*it)->alt_loc_id != (*jt)->alt_loc_id) continue;
                if (((*jt)->bond_ind == 9) || ((*jt)->bond_ind == 26) || ((*jt)->bond_ind == 27)) {
                    (*it)->bonded_atoms.push_back(*jt);
                    (*jt)->bonded_atoms.push_back(*it);
                    create_bond(*it,*jt,"ar",dt);
                }
            }
        }
        
        if ((*it)->bond_ind == 38) {
            for (atoms_vec jt=atoms.begin(); jt!=atoms.end();++jt) {
                if ((*it)->alt_loc_id != (*jt)->alt_loc_id) continue;
                if (((*jt)->bond_ind == 10) || ((*jt)->bond_ind == 12)) {
                    (*it)->bonded_atoms.push_back(*jt);
                    (*jt)->bonded_atoms.push_back(*it);
                    create_bond(*it,*jt,"1",dt);
                }
                
            }
        }
        
        switch ((*it)->bond_ind) {
            case 80: for (atoms_vec jt=atoms.begin(); jt!=atoms.end();++jt) {
                    if ((*it)->alt_loc_id != (*jt)->alt_loc_id) continue;
                    if ((*jt)->bond_ind == 1) {
                        (*it)->bonded_atoms.push_back(*jt);
                        (*jt)->bonded_atoms.push_back(*it);
                        dt = "BACKBONE | DICT";
                        create_bond(*jt,*it,"1",dt);
                    }
                } break;
            case 81: for (atoms_vec jt=atoms.begin(); jt!=atoms.end();++jt) {
                    if ((*it)->alt_loc_id != (*jt)->alt_loc_id) continue;
                    if ((*jt)->bond_ind == 2) {
                        (*it)->bonded_atoms.push_back(*jt);
                        (*jt)->bonded_atoms.push_back(*it);
                        dt = "BACKBONE | DICT";
                        create_bond(*jt,*it,"1",dt);
                    }
                } break;
            case 82: for (atoms_vec jt=atoms.begin(); jt!=atoms.end();++jt) {
                    if ((*it)->alt_loc_id != (*jt)->alt_loc_id) continue;
                    if ((*jt)->bond_ind == 5) {
                        (*it)->bonded_atoms.push_back(*jt);
                        (*jt)->bonded_atoms.push_back(*it);
                        dt = "DICT";
                        create_bond(*jt,*it,"1",dt);
                    }
                } break;
            case 83: for (atoms_vec jt=atoms.begin(); jt!=atoms.end();++jt) {
                    if ((*it)->alt_loc_id != (*jt)->alt_loc_id) continue;
                    if ((*jt)->bond_ind == 6) {
                        (*it)->bonded_atoms.push_back(*jt);
                        (*jt)->bonded_atoms.push_back(*it);
                        dt = "DICT";
                        create_bond(*jt,*it,"1",dt);
                    }
                } break;
            case 84: for (atoms_vec jt=atoms.begin(); jt!=atoms.end();++jt) {
                    if ((*it)->alt_loc_id != (*jt)->alt_loc_id) continue;
                    if ((*jt)->bond_ind == 7) {
                        (*it)->bonded_atoms.push_back(*jt);
                        (*jt)->bonded_atoms.push_back(*it);
                        dt = "DICT";
                        create_bond(*jt,*it,"1",dt);
                    }
                } break;
            case 85: for (atoms_vec jt=atoms.begin(); jt!=atoms.end();++jt) {
                    if ((*it)->alt_loc_id != (*jt)->alt_loc_id) continue;
                    if ((*jt)->bond_ind == 8) {
                        (*it)->bonded_atoms.push_back(*jt);
                        (*jt)->bonded_atoms.push_back(*it);
                        dt = "DICT";
                        create_bond(*jt,*it,"1",dt);
                    }
                } break;
            case 86: for (atoms_vec jt=atoms.begin(); jt!=atoms.end();++jt) {
                    if ((*it)->alt_loc_id != (*jt)->alt_loc_id) continue;
                    if ((*jt)->bond_ind == 9) {
                        (*it)->bonded_atoms.push_back(*jt);
                        (*jt)->bonded_atoms.push_back(*it);
                        dt = "DICT";
                        create_bond(*jt,*it,"1",dt);
                    }
                } break;
            case 87: for (atoms_vec jt=atoms.begin(); jt!=atoms.end();++jt) {
                    if ((*it)->alt_loc_id != (*jt)->alt_loc_id) continue;
                    if ((*jt)->bond_ind == 18) {
                        (*it)->bonded_atoms.push_back(*jt);
                        (*jt)->bonded_atoms.push_back(*it);
                        dt = "DICT";
                        create_bond(*jt,*it,"1",dt);
                    }
                } break;
            case 88: for (atoms_vec jt=atoms.begin(); jt!=atoms.end();++jt) {
                    if ((*it)->alt_loc_id != (*jt)->alt_loc_id) continue;
                    if ((*jt)->bond_ind == 19) {
                        (*it)->bonded_atoms.push_back(*jt);
                        (*jt)->bonded_atoms.push_back(*it);
                        dt = "DICT";
                        create_bond(*jt,*it,"1",dt);
                    }
                } break;
            case 89: for (atoms_vec jt=atoms.begin(); jt!=atoms.end();++jt) {
                    if ((*it)->alt_loc_id != (*jt)->alt_loc_id) continue;
                    if ((*jt)->bond_ind == 20) {
                        (*it)->bonded_atoms.push_back(*jt);
                        (*jt)->bonded_atoms.push_back(*it);
                        dt = "DICT";
                        create_bond(*jt,*it,"1",dt);
                    }
                } break;
            case 90: for (atoms_vec jt=atoms.begin(); jt!=atoms.end();++jt) {
                    if ((*it)->alt_loc_id != (*jt)->alt_loc_id) continue;
                    if ((*jt)->bond_ind == 21) {
                        (*it)->bonded_atoms.push_back(*jt);
                        (*jt)->bonded_atoms.push_back(*it);
                        dt = "DICT";
                        create_bond(*jt,*it,"1",dt);
                    }
                } break;
            case 91: for (atoms_vec jt=atoms.begin(); jt!=atoms.end();++jt) {
                    if ((*it)->alt_loc_id != (*jt)->alt_loc_id) continue;
                    if ((*jt)->bond_ind == 22) {
                        (*it)->bonded_atoms.push_back(*jt);
                        (*jt)->bonded_atoms.push_back(*it);
                        dt = "DICT";
                        create_bond(*jt,*it,"1",dt);
                    }
                } break;
            case 92: for (atoms_vec jt=atoms.begin(); jt!=atoms.end();++jt) {
                    if ((*it)->alt_loc_id != (*jt)->alt_loc_id) continue;
                    if ((*jt)->bond_ind == 13) {
                        (*it)->bonded_atoms.push_back(*jt);
                        (*jt)->bonded_atoms.push_back(*it);
                        dt = "DICT";
                        create_bond(*jt,*it,"1",dt);
                    }
                } break;
            case 93: for (atoms_vec jt=atoms.begin(); jt!=atoms.end();++jt) {
                    if ((*it)->alt_loc_id != (*jt)->alt_loc_id) continue;
                    if ((*jt)->bond_ind == 49) {
                        (*it)->bonded_atoms.push_back(*jt);
                        (*jt)->bonded_atoms.push_back(*it);
                        dt = "DICT";
                        create_bond(*jt,*it,"1",dt);
                    }
                } break;
            case 94: for (atoms_vec jt=atoms.begin(); jt!=atoms.end();++jt) {
                    if ((*it)->alt_loc_id != (*jt)->alt_loc_id) continue;
                    if ((*jt)->bond_ind == 51) {
                        (*it)->bonded_atoms.push_back(*jt);
                        (*jt)->bonded_atoms.push_back(*it);
                        dt = "DICT";
                        create_bond(*jt,*it,"1",dt);
                    }
                } break;
            case 95: for (atoms_vec jt=atoms.begin(); jt!=atoms.end();++jt) {
                    if ((*it)->alt_loc_id != (*jt)->alt_loc_id) continue;
                    if ((*jt)->bond_ind == 38) {
                        (*it)->bonded_atoms.push_back(*jt);
                        (*jt)->bonded_atoms.push_back(*it);
                        dt = "DICT";
                        create_bond(*jt,*it,"1",dt);
                    }
                } break;
            case 96: for (atoms_vec jt=atoms.begin(); jt!=atoms.end();++jt) {
                    if ((*it)->alt_loc_id != (*jt)->alt_loc_id) continue;
                    if ((*jt)->bond_ind == 10) {
                        (*it)->bonded_atoms.push_back(*jt);
                        (*jt)->bonded_atoms.push_back(*it);
                        dt = "DICT";
                        create_bond(*jt,*it,"1",dt);
                    }
                } break;
            case 97: for (atoms_vec jt=atoms.begin(); jt!=atoms.end();++jt) {
                    if ((*it)->alt_loc_id != (*jt)->alt_loc_id) continue;
                    if ((*jt)->bond_ind == 26) {
                        (*it)->bonded_atoms.push_back(*jt);
                        (*jt)->bonded_atoms.push_back(*it);
                        dt = "DICT";
                        create_bond(*jt,*it,"1",dt);
                    }
                } break;
            case 98: for (atoms_vec jt=atoms.begin(); jt!=atoms.end();++jt) {
                    if ((*it)->alt_loc_id != (*jt)->alt_loc_id) continue;
                    if ((*jt)->bond_ind == 25) {
                        (*it)->bonded_atoms.push_back(*jt);
                        (*jt)->bonded_atoms.push_back(*it);
                        dt = "DICT";
                        create_bond(*jt,*it,"1",dt);
                    }
                } break;
            case 99: for (atoms_vec jt=atoms.begin(); jt!=atoms.end();++jt) {
                    if ((*it)->alt_loc_id != (*jt)->alt_loc_id) continue;
                    if ((*jt)->bond_ind == 11) {
                        (*it)->bonded_atoms.push_back(*jt);
                        (*jt)->bonded_atoms.push_back(*it);
                        dt = "DICT";
                        create_bond(*jt,*it,"1",dt);
                    }
                } break;
            case 100: for (atoms_vec jt=atoms.begin(); jt!=atoms.end();++jt) {
                    if ((*it)->alt_loc_id != (*jt)->alt_loc_id) continue;
                    if ((*jt)->bond_ind == 12) {
                        (*it)->bonded_atoms.push_back(*jt);
                        (*jt)->bonded_atoms.push_back(*it);
                        dt = "DICT";
                        create_bond(*jt,*it,"1",dt);
                    }
                } break;
            case 101: for (atoms_vec jt=atoms.begin(); jt!=atoms.end();++jt) {
                    if ((*it)->alt_loc_id != (*jt)->alt_loc_id) continue;
                    if ((*jt)->bond_ind == 61) {
                        (*it)->bonded_atoms.push_back(*jt);
                        (*jt)->bonded_atoms.push_back(*it);
                        dt = "CAP | DICT";
                        create_bond(*jt,*it,"1",dt);
                    }
                } break;
            case 102: for (atoms_vec jt=atoms.begin(); jt!=atoms.end();++jt) {
                    if ((*it)->alt_loc_id != (*jt)->alt_loc_id) continue;
                    if ((*jt)->bond_ind == 27) {
                        (*it)->bonded_atoms.push_back(*jt);
                        (*jt)->bonded_atoms.push_back(*it);
                        dt = "DICT";
                        create_bond(*jt,*it,"1",dt);
                    }
                } break;
        }
    }
}


//==============================================================================================
//Definitionen fuer CHAIN:
//==============================================================================================

CHAIN::CHAIN() :MOLECULE() {n_aacids = 0;}

CHAIN::~CHAIN() {
    for (aacids_map it=aacids.begin(); it!=aacids.end();++it) {
        it->second.kill();
    }
    for (reli_centers_vec it=reli_centers.begin(); it!=reli_centers.end();++it) {
        it->kill();
    }
}

void CHAIN::get_atom_types(int mode,const char *def_file,bool fill_X) {
    for (aacids_map jt=aacids.begin(); jt!=aacids.end();++jt) {
        jt->second->get_intern_types();
    }
    
    bool prot_acids = prot_acids_default;
    bool prot_guanidin = prot_guanidin_default;
    bool prot_amidin = prot_amidin_default;
    bool prot_amin = prot_amin_default;
    bool prot_phosphate = prot_phosphate_default;
    bool prot_sulfate = prot_sulfate_default;
    bool get_bonds = get_bonds_default;
    bool kekulize_aromatics = kekulize_aromatics_default;
    bool kekulize_charged = kekulize_charged_default;
    bool allow_charged_aromatics = allow_charged_aromatics_default;
    int max_ring_members = max_ring_members_default;
    try_for_def_file(mode,def_file,"X",prot_acids,prot_guanidin,prot_amidin,
                     prot_amin,prot_phosphate,prot_sulfate,get_bonds,max_ring_members,
                     kekulize_aromatics,kekulize_charged,allow_charged_aromatics);

    intern2mode(fill_X);
}

void CHAIN::get_aacids(int const& preset_n_aacids) { //!dringend noch kommentieren
    n_aacids = preset_n_aacids;
    last_name = "XXX";
    
    for (atoms_vec jt=atoms.begin(); jt!=atoms.end();++jt) {
        if ((*jt)->full_res_name != last_name) {
            ++n_aacids;
            stl_ptr<AACID> aac(new AACID());
            aac->chain = this;
            aac->res_name = (*jt)->res_name;
            aac->res_number = (*jt)->res_number;
            aacids[n_aacids] = aac;
            last_name = (*jt)->full_res_name;
        }
        
        aacids[n_aacids]->atoms.push_back(*jt);
    }
}

void CHAIN::build_bonds() { //nicht schoen, aber laeuft
    bool b_flag;
    protein->last_C = NULL;
    //Fuer alle aacids:
    for (aacids_map ut=aacids.begin(); ut!=aacids.end();++ut) {
        ut->second->build_bonds();
    }
    //Jetzt die Bindungen zu den CONECT-Eintraegen generieren
    //!Aber nur fuer die Proteinatome  --   bisher hoch ineffizient -> Problem dass die maps zum Parser gehoeren
    //Fuer alle CONECT-Eintraege:
    for (conects_vec it=protein->main_structure->conects.begin(); it!=protein->main_structure->conects.end();++it) {
        //Fuer alle Atome der aktuellen Chain:
        for (atoms_vec kt=atoms.begin(); kt!=atoms.end();++kt) {
            //Wenn das aktuelle Atom einen CONECT-Eintrag hat:
            if ((*it)->from_id == (*kt)->intern_id) {
                //Fuer alle kovalenten Bindungen:
                for (int i=0; i<4; ++i) {
                    //Eine Bindung existiert:
                    if ((*it)->cov_bonds[i] != 0) {
                        b_flag = false;
                        //Fuer alle Chains:
                        for (chains_vec lt=protein->chains.begin(); lt!=protein->chains.end();++lt) {
                            //Fuer alle Atome:
                            for (atoms_vec mt=(*lt)->atoms.begin(); mt!=(*lt)->atoms.end();++mt) {
                                //Wenn das aktuelle Atom ein Bindungspartner ist:
                                if ((*it)->cov_bonds[i] == (*mt)->intern_id) {
                                    //nur wenn das mt nicht schon in den bonded atm von kt ist:
                                    for (atoms_vec wt=(*kt)->bonded_atoms.begin();
                                            wt!=(*kt)->bonded_atoms.end();++wt) {
                                        if ((*wt)->intern_id == (*mt)->intern_id) {
                                            b_flag = true;
                                            break;
                                        }
                                    }
                                    if (!b_flag) {
                                        (*kt)->bonded_atoms.push_back(*mt);
                                        (*mt)->bonded_atoms.push_back(*kt);
                                        stl_ptr<BOND> bnd(new BOND());
                                        bonds.push_back(bnd);
                                        bnd->from = *kt;
                                        bnd->to = *mt;
                                        bnd->type = "1"; //!obwohl nicht sicher
                                        bnd->id = protein->bnd_id;
                                        protein->bnd_id += 1;
                                        bnd->dict_type = " ";
                                        //es muss nicht weiter gesucht werden:
                                        b_flag = true;
                                        break;
                                    } else break;
                                }
                                if (b_flag) break;
                            }
                        }
                    }
                }
            }
        }
    }
}

void CHAIN::get_reli_centers() {
    //! obsolet, weil eine veraltete Definition der Relibase Pseudozentren
    for (aacids_map jt=aacids.begin(); jt!=aacids.end();++jt) {
        //!zunaechst Backbone
        string name;
        istringstream is;
        for (atoms_vec it=atoms.begin(); it!=atoms.end();++it) {
            is.clear();
            is.str((*it)->name);
            is >> name;
            if (name == "O") {
                stl_ptr<RELIBASE_CENTER> rbc(new RELIBASE_CENTER());
                rbc->aacid = jt->second;
                rbc->type = "Acceptor";
                rbc->coord = (*it)->coord;
                reli_centers.push_back(rbc);
            } else if (name == "N") {
                stl_ptr<RELIBASE_CENTER> rbc(new RELIBASE_CENTER());
                rbc->aacid = jt->second;
                rbc->type = "Donor";
                rbc->coord = (*it)->coord;
                reli_centers.push_back(rbc);
            } else if (name == "C") {
                stl_ptr<RELIBASE_CENTER> rbc(new RELIBASE_CENTER());
                rbc->aacid = jt->second;
                rbc->type = "PI";
                rbc->coord = (*it)->coord;
                reli_centers.push_back(rbc);
            }
        }
        if (jt->second->res_name == "ALA") {
            for (atoms_vec it=jt->second->atoms.begin(); it!=jt->second->atoms.end();++it) {
                is.clear();
                is.str((*it)->name);
                is >> name;
                if (name == "CB") {
                    stl_ptr<RELIBASE_CENTER> rbc(new RELIBASE_CENTER());
                    rbc->aacid = jt->second;
                    rbc->type = "Aliphatic";
                    rbc->coord = (*it)->coord;
                    reli_centers.push_back(rbc);
                    break;
                }
            }
        } else if (jt->second->res_name == "ARG") {
            vec3d<float> buff_vec(0.,0.,0.);
            for (atoms_vec it=jt->second->atoms.begin(); it!=jt->second->atoms.end();++it) {
                is.clear();
                is.str((*it)->name);
                is >> name;
                if (name == "NE") {
                    stl_ptr<RELIBASE_CENTER> rbc(new RELIBASE_CENTER());
                    rbc->aacid = jt->second;
                    rbc->type = "DONOR";
                    rbc->coord = (*it)->coord;
                    reli_centers.push_back(rbc);
                } else if (name == "NH1") {
                    stl_ptr<RELIBASE_CENTER> rbc(new RELIBASE_CENTER());
                    rbc->aacid = jt->second;
                    rbc->type = "Donor";
                    rbc->coord = (*it)->coord;
                    reli_centers.push_back(rbc);
                } else if (name == "NH2") {
                    stl_ptr<RELIBASE_CENTER> rbc(new RELIBASE_CENTER());
                    rbc->aacid = jt->second;
                    rbc->type = "Donor";
                    rbc->coord = (*it)->coord;
                    reli_centers.push_back(rbc);
                } else if (name == "CB" || name == "CG" || name == "CD") {
                    buff_vec += (*it)->coord;
                }
            }
            buff_vec /= 3.;
            stl_ptr<RELIBASE_CENTER> rbc(new RELIBASE_CENTER());
            rbc->aacid = jt->second;
            rbc->type = "Aliphatic";
            rbc->coord = buff_vec;
            reli_centers.push_back(rbc);
        } else if (jt->second->res_name == "ASN") {
            for (atoms_vec it=jt->second->atoms.begin(); it!=jt->second->atoms.end();++it) {
                is.clear();
                is.str((*it)->name);
                is >> name;
                if (name == "OD1") {
                    stl_ptr<RELIBASE_CENTER> rbc(new RELIBASE_CENTER());
                    rbc->aacid = jt->second;
                    rbc->type = "Acceptor";
                    rbc->coord = (*it)->coord;
                    reli_centers.push_back(rbc);
                } else if (name == "ND2") {
                    stl_ptr<RELIBASE_CENTER> rbc(new RELIBASE_CENTER());
                    rbc->aacid = jt->second;
                    rbc->type = "Donor";
                    rbc->coord = (*it)->coord;
                    reli_centers.push_back(rbc);
                }
            }
        } else if (jt->second->res_name == "ASP" || jt->second->res_name == "ASH") {
            for (atoms_vec it=jt->second->atoms.begin(); it!=jt->second->atoms.end();++it) {
                is.clear();
                is.str((*it)->name);
                is >> name;
                if (name == "OD1") {
                    stl_ptr<RELIBASE_CENTER> rbc(new RELIBASE_CENTER());
                    rbc->aacid = jt->second;
                    rbc->type = "Acceptor";
                    rbc->coord = (*it)->coord;
                    reli_centers.push_back(rbc);
                } else if (name == "OD2") {
                    stl_ptr<RELIBASE_CENTER> rbc(new RELIBASE_CENTER());
                    rbc->aacid = jt->second;
                    rbc->type = "Acceptor";
                    rbc->coord = (*it)->coord;
                    reli_centers.push_back(rbc);
                }
            }
        } else if (jt->second->res_name == "ASX") {
            if (protein->main_structure->verb > 0) {
                cerr << c_message<cWARNING>("CHAIN::get_reli_centers() --> no centres will be set for ASX") << endl;
            }
        } else if (jt->second->res_name == "CYS" || jt->second->res_name == "CYX" || jt->second->res_name == "CYM") {
            vec3d<float> buff_vec(0.,0.,0.);
            for (atoms_vec it=jt->second->atoms.begin(); it!=jt->second->atoms.end();++it) {
                is.clear();
                is.str((*it)->name);
                is >> name;
                if (name == "CB" || name == "SG") {
                    buff_vec += (*it)->coord;
                }
            }
            buff_vec /= 2.;
            stl_ptr<RELIBASE_CENTER> rbc(new RELIBASE_CENTER());
            rbc->aacid = jt->second;
            rbc->type = "Aliphatic";
            rbc->coord = buff_vec;
            reli_centers.push_back(rbc);
        } else if (jt->second->res_name == "GLN") {
            for (atoms_vec it=jt->second->atoms.begin(); it!=jt->second->atoms.end();++it) {
                is.clear();
                is.str((*it)->name);
                is >> name;
                if (name == "OE1") {
                    stl_ptr<RELIBASE_CENTER> rbc(new RELIBASE_CENTER());
                    rbc->aacid = jt->second;
                    rbc->type = "Acceptor";
                    rbc->coord = (*it)->coord;
                    reli_centers.push_back(rbc);
                } else if (name == "NE2") {
                    stl_ptr<RELIBASE_CENTER> rbc(new RELIBASE_CENTER());
                    rbc->aacid = jt->second;
                    rbc->type = "Donor";
                    rbc->coord = (*it)->coord;
                    reli_centers.push_back(rbc);
                }
            }
        } else if (jt->second->res_name == "GLU" || jt->second->res_name == "GLH") {
            for (atoms_vec it=jt->second->atoms.begin(); it!=jt->second->atoms.end();++it) {
                is.clear();
                is.str((*it)->name);
                is >> name;
                if (name == "OE1") {
                    stl_ptr<RELIBASE_CENTER> rbc(new RELIBASE_CENTER());
                    rbc->aacid = jt->second;
                    rbc->type = "Acceptor";
                    rbc->coord = (*it)->coord;
                    reli_centers.push_back(rbc);
                } else if (name == "OE2") {
                    stl_ptr<RELIBASE_CENTER> rbc(new RELIBASE_CENTER());
                    rbc->aacid = jt->second;
                    rbc->type = "Acceptor";
                    rbc->coord = (*it)->coord;
                    reli_centers.push_back(rbc);
                }
            }
        } else if (jt->second->res_name == "GLX") {
            if (protein->main_structure->verb > 0) {
                cerr << c_message<cWARNING>("CHAIN::get_reli_centers() --> no centres will be set for GLX") << endl;
            }
        } else if (jt->second->res_name == "ILE") {
            vec3d<float> buff_vec(0.,0.,0.);
            for (atoms_vec it=jt->second->atoms.begin(); it!=jt->second->atoms.end();++it) {
                is.clear();
                is.str((*it)->name);
                is >> name;
                if (name == "CB" || name == "CG1" || name == "CG2" || name == "CD1") {
                    buff_vec += (*it)->coord;
                }
            }
            buff_vec /= 4.;
            stl_ptr<RELIBASE_CENTER> rbc(new RELIBASE_CENTER());
            rbc->aacid = jt->second;
            rbc->type = "Aliphatic";
            rbc->coord = buff_vec;
            reli_centers.push_back(rbc);
        } else if (jt->second->res_name == "LEU") {
            vec3d<float> buff_vec(0.,0.,0.);
            for (atoms_vec it=jt->second->atoms.begin(); it!=jt->second->atoms.end();++it) {
                is.clear();
                is.str((*it)->name);
                is >> name;
                if (name == "CB" || name == "CG" || name == "CD1" || name == "CD2") {
                    buff_vec += (*it)->coord;
                }
            }
            buff_vec /= 4.;
            stl_ptr<RELIBASE_CENTER> rbc(new RELIBASE_CENTER());
            rbc->aacid = jt->second;
            rbc->type = "Aliphatic";
            rbc->coord = buff_vec;
            reli_centers.push_back(rbc);
        } else if (jt->second->res_name == "LYS" || jt->second->res_name == "LYN") {
            vec3d<float> buff_vec(0.,0.,0.);
            for (atoms_vec it=jt->second->atoms.begin(); it!=jt->second->atoms.end();++it) {
                is.clear();
                is.str((*it)->name);
                is >> name;
                if (name == "CB" || name == "CG" || name == "CD" || name == "CE") {
                    buff_vec += (*it)->coord;
                } else if (name == "NZ") {
                    stl_ptr<RELIBASE_CENTER> rbc(new RELIBASE_CENTER());
                    rbc->aacid = jt->second;
                    rbc->type = "Donor";
                    rbc->coord = (*it)->coord;
                    reli_centers.push_back(rbc);
                }
            }
            buff_vec /= 4.;
            stl_ptr<RELIBASE_CENTER> rbc(new RELIBASE_CENTER());
            rbc->aacid = jt->second;
            rbc->type = "Aliphatic";
            rbc->coord = buff_vec;
            reli_centers.push_back(rbc);
        } else if (jt->second->res_name == "MET") {
            vec3d<float> buff_vec(0.,0.,0.);
            for (atoms_vec it=jt->second->atoms.begin(); it!=jt->second->atoms.end();++it) {
                is.clear();
                is.str((*it)->name);
                is >> name;
                if (name == "CB" || name == "CG" || name == "SD" || name == "CE") {
                    buff_vec += (*it)->coord;
                }
            }
            buff_vec /= 4.;
            stl_ptr<RELIBASE_CENTER> rbc(new RELIBASE_CENTER());
            rbc->aacid = jt->second;
            rbc->type = "Aliphatic";
            rbc->coord = buff_vec;
            reli_centers.push_back(rbc);
        } else if (jt->second->res_name == "PHE") {
            vec3d<float> buff_vec(0.,0.,0.);
            for (atoms_vec it=jt->second->atoms.begin(); it!=jt->second->atoms.end();++it) {
                is.clear();
                is.str((*it)->name);
                is >> name;
                if (name == "CG" || name == "CD1" || name == "CD2" ||
                    name == "CE1" || name == "CE2" || name == "CZ") {
                    buff_vec += (*it)->coord;
                }
            }
            buff_vec /= 6.;
            stl_ptr<RELIBASE_CENTER> rbc(new RELIBASE_CENTER());
            rbc->aacid = jt->second;
            rbc->type = "PI";
            rbc->coord = buff_vec;
            reli_centers.push_back(rbc);
        } else if (jt->second->res_name == "PRO") {
            vec3d<float> buff_vec(0.,0.,0.);
            for (atoms_vec it=jt->second->atoms.begin(); it!=jt->second->atoms.end();++it) {
                is.clear();
                is.str((*it)->name);
                is >> name;
                if (name == "CB" || name == "CG" || name == "CD") {
                    buff_vec += (*it)->coord;
                }
            }
            buff_vec /= 3.;
            stl_ptr<RELIBASE_CENTER> rbc(new RELIBASE_CENTER());
            rbc->aacid = jt->second;
            rbc->type = "Aliphatic";
            rbc->coord = buff_vec;
            reli_centers.push_back(rbc);
        } else if (jt->second->res_name == "SER") {
            for (atoms_vec it=jt->second->atoms.begin(); it!=jt->second->atoms.end();++it) {
                is.clear();
                is.str((*it)->name);
                is >> name;
                if (name == "OG") {
                    stl_ptr<RELIBASE_CENTER> rbc(new RELIBASE_CENTER());
                    rbc->aacid = jt->second;
                    rbc->type = "DON_ACC";
                    rbc->coord = (*it)->coord;
                    reli_centers.push_back(rbc);
                }
            }
        } else if (jt->second->res_name == "THR") {
            for (atoms_vec it=jt->second->atoms.begin(); it!=jt->second->atoms.end();++it) {
                is.clear();
                is.str((*it)->name);
                is >> name;
                if (name == "CD2") {
                    stl_ptr<RELIBASE_CENTER> rbc(new RELIBASE_CENTER());
                    rbc->aacid = jt->second;
                    rbc->type = "Aliphatic";
                    rbc->coord = (*it)->coord;
                    reli_centers.push_back(rbc);
                } else if (name == "OD1") {
                    stl_ptr<RELIBASE_CENTER> rbc(new RELIBASE_CENTER());
                    rbc->aacid = jt->second;
                    rbc->type = "DON_ACC";
                    rbc->coord = (*it)->coord;
                    reli_centers.push_back(rbc);
                }
            }
        } else if (jt->second->res_name == "TRP") {
            vec3d<float> buff_vec(0.,0.,0.);
            for (atoms_vec it=jt->second->atoms.begin(); it!=jt->second->atoms.end();++it) {
                is.clear();
                is.str((*it)->name);
                is >> name;
                if (name == "CG" || name == "CD1" || name == "CD2" || name == "NE1" ||
                    name == "CE2" || name == "CE3" || name == "CZ1" || name == "CZ3" || name == "CH2") {
                    buff_vec += (*it)->coord;
                } else if (name == "NE1") {
                    stl_ptr<RELIBASE_CENTER> rbc(new RELIBASE_CENTER());
                    rbc->aacid = jt->second;
                    rbc->type = "Donor";
                    rbc->coord = (*it)->coord;
                    reli_centers.push_back(rbc);
                }
            }
            buff_vec /= 9.;
            stl_ptr<RELIBASE_CENTER> rbc(new RELIBASE_CENTER());
            rbc->aacid = jt->second;
            rbc->type = "PI";
            rbc->coord = buff_vec;
            reli_centers.push_back(rbc);
        } else if (jt->second->res_name == "TYR" || jt->second->res_name == "TYM") {
            vec3d<float> buff_vec(0.,0.,0.);
            for (atoms_vec it=jt->second->atoms.begin(); it!=jt->second->atoms.end();++it) {
                is.clear();
                is.str((*it)->name);
                is >> name;
                if (name == "CG" || name == "CD1" || name == "CD2" ||
                    name == "CE1" || name == "CE2" || name == "CZ") {
                    buff_vec += (*it)->coord;
                } else if (name == "OH") {
                    stl_ptr<RELIBASE_CENTER> rbc(new RELIBASE_CENTER());
                    rbc->aacid = jt->second;
                    rbc->type = "DON_ACC";
                    rbc->coord = (*it)->coord;
                    reli_centers.push_back(rbc);
                }
            }
            buff_vec /= 6.;
            stl_ptr<RELIBASE_CENTER> rbc(new RELIBASE_CENTER());
            rbc->aacid = jt->second;
            rbc->type = "PI";
            rbc->coord = buff_vec;
            reli_centers.push_back(rbc);
        } else if (jt->second->res_name == "VAL") {
            vec3d<float> buff_vec(0.,0.,0.);
            for (atoms_vec it=jt->second->atoms.begin(); it!=jt->second->atoms.end();++it) {
                is.clear();
                is.str((*it)->name);
                is >> name;
                if (name == "CB" || name == "CG1" || name == "CG2") {
                    buff_vec += (*it)->coord;
                }
            }
            buff_vec /= 3.;
            stl_ptr<RELIBASE_CENTER> rbc(new RELIBASE_CENTER());
            rbc->aacid = jt->second;
            rbc->type = "Aliphatic";
            rbc->coord = buff_vec;
            reli_centers.push_back(rbc);
        } else if (jt->second->res_name == "HIS" || jt->second->res_name == "HID" ||
                   jt->second->res_name == "HIE" || jt->second->res_name == "HIP") {
            vec3d<float> buff_vec(0.,0.,0.);
            for (atoms_vec it=jt->second->atoms.begin(); it!=jt->second->atoms.end();++it) {
                is.clear();
                is.str((*it)->name);
                is >> name;
                if (name == "CG" || name == "ND1" || name == "CD2" ||
                    name == "CE1" || name == "CE2") {
                    buff_vec += (*it)->coord;
                } else if (name == "NE1") {
                    stl_ptr<RELIBASE_CENTER> rbc(new RELIBASE_CENTER());
                    rbc->aacid = jt->second;
                    rbc->type = "DON_ACC";
                    rbc->coord = (*it)->coord;
                    reli_centers.push_back(rbc);
                } else if (name == "NE2") {
                    stl_ptr<RELIBASE_CENTER> rbc(new RELIBASE_CENTER());
                    rbc->aacid = jt->second;
                    rbc->type = "DON_ACC";
                    rbc->coord = (*it)->coord;
                    reli_centers.push_back(rbc);
                }
            }
            buff_vec /= 5.;
            stl_ptr<RELIBASE_CENTER> rbc(new RELIBASE_CENTER());
            rbc->aacid = jt->second;
            rbc->type = "PI";
            rbc->coord = buff_vec;
            reli_centers.push_back(rbc);
        }
    }
}


//==============================================================================================
//Definitionen fuer CAVITY:
//==============================================================================================

CAVITY::CAVITY() {}

CAVITY::~CAVITY() {} //die aacids werden von protein geloescht!

bool CAVITY::operator<(CAVITY &rechts) {
    if (total_volume < rechts.total_volume) return true;
    return false;
}


void CAVITY::align(stl_ptr<CAVITY> ref) {
    //! Listen fuer das Sequenzalignment erstellen:
    vector<string> Aseq;
    vector<string> Bseq;
    string gap_val = "XXX";
    
    int a_gaps = 0;
    int b_gaps = 0;
    
    int last_res_num = aacids[0]->res_number;
    for (aacids_vec aa=aacids.begin(); aa!=aacids.end(); ++aa) {
        int n_gaps = abs((*aa)->res_number - last_res_num);
        if (n_gaps > 9) n_gaps = 10;
        for (int i=1; i<n_gaps; ++i) {
            Bseq.push_back(gap_val);
            a_gaps++;
        }
        last_res_num = (*aa)->res_number;
        Bseq.push_back((*aa)->res_name);
    }
    last_res_num = ref->aacids[0]->res_number;
    for (aacids_vec aa=ref->aacids.begin(); aa!=ref->aacids.end(); ++aa) {
        int n_gaps = abs((*aa)->res_number - last_res_num);
        if (n_gaps > 9) n_gaps = 10;
        for (int i=1; i<n_gaps; ++i) {
            Aseq.push_back(gap_val);
            b_gaps++;
        }
        last_res_num = (*aa)->res_number;
        Aseq.push_back((*aa)->res_name);
    }
    
    
    if (aacids[0]->chain->protein->main_structure->verb > 0) {
        cout << " -> " << Aseq.size() - a_gaps << " residues in reference structure" << "\n";
        cout << " -> " << Bseq.size() - b_gaps << " residues in alignment structure" << "\n";
    }
    
    SEQALIGN<string> my_sa(Aseq,Bseq,gap_val);
    my_sa.needleman_wunsch1(); //!Arbeitet hier eher schlecht => lieber Smith-Waterman nehmen
                               //!Dann allerdings nur je einen kuenstlichen Gap fuer unterbrochene
                               //!Ketten einsetzen!
    //my_sa.smith_waterman();
    
    //! Listen fuer das geometrische Alignment erstellen:
    vector<vec3d<float> > A_list;
    vector<vec3d<float> > B_list;
    
    unsigned int matched = 0;
    unsigned int Anum = 0;
    unsigned int Bnum = 0;
    tr1::unordered_set<int> Apos;
    tr1::unordered_set<int> Bpos;
    for (unsigned int i=0; i<my_sa.ires.size(); ++i) {
        if (i >= my_sa.jres.size()) break;
        
        //if (!(my_sa.ires[i] == "XXX" || my_sa.jres[i] == "XXX")) {
        if (my_sa.ires[i] == my_sa.jres[i] && !(my_sa.ires[i] == "XXX" && my_sa.jres[i] == "XXX")) {
            matched++;
            
            Apos.insert(Anum);
            Bpos.insert(Bnum);
            ++Anum;
            ++Bnum;
        } else {
            if (my_sa.ires[i] != "XXX") Anum++;
            if (my_sa.jres[i] != "XXX") Bnum++;
        }
    }
    
    if (aacids[0]->chain->protein->main_structure->verb > 0) cout << " -> " << matched << " residues matched in sequence alignment" << "\n";
    
    Anum = 0;
    Bnum = 0;
    for (aacids_vec aa=aacids.begin(); aa!=aacids.end(); ++aa) {
        if (Bpos.find(Bnum) != Bpos.end()) {
            istringstream is;
            string name;
            for (atoms_vec at=(*aa)->atoms.begin(); at!=(*aa)->atoms.end(); ++at) {
                is.clear(); //streamobjekt zurcksetzen
                is.str((*at)->name);
                is >> name;
                if (name == "CA") B_list.push_back(vec3d<float>((*at)->coord));
                //A_list.push_back(vec3d<float>((*at)->coord));
            }
        }
        Bnum++;
    }
    
    for (aacids_vec aa=ref->aacids.begin(); aa!=ref->aacids.end(); ++aa) {
        if (Apos.find(Anum) != Apos.end()) {
            istringstream is;
            string name;
            for (atoms_vec at=(*aa)->atoms.begin(); at!=(*aa)->atoms.end(); ++at) {
                is.clear(); //streamobjekt zurcksetzen
                is.str((*at)->name);
                is >> name;
                if (name == "CA") A_list.push_back(vec3d<float>((*at)->coord));
                //B_list.push_back(vec3d<float>((*at)->coord));
            }
        }
        Anum++;
    }
    
    //!Jetzt die Symmetrieelemente fuer das geometrische Alignment berechnen:
    get_align_matrix(A_list,B_list,optalign_rotm,optalign_trans[0],optalign_trans[1]);
}


void write_def_file(int const& mode,const char *name,vector<string> *alt_mode,
                    bool prot_acids,
                    bool prot_guanidin,
                    bool prot_amidin,
                    bool prot_amin,
                    bool prot_phosphate,
                    bool prot_sulfate,
                    bool get_bonds,
                    int max_ring_members,
                    bool kekulize_aromatics,
                    bool kekulize_charged,
                    bool allow_charged_aromatics) {
    ofstream f_out;
    f_out.open(name);
    
    //!zunaechst den Kommentar schreiben:
    f_out << "This is an atom type definitions file from 'molecule_GN.cpp'" << "\n" << "\n";
    f_out << "<COMMENT>" << "\n";
    f_out << "--> Everything in the COMMENT sections will NOT be read by the program (it's just comments)." << "\n";
    f_out << "--> In the HEADER section the 'internal_set_type' and the 'total' value is mandatory" << "\n";
    f_out << "    and should NEVER be changed by the user!" << "\n";
    f_out << "--> The FLAGS section and the DEF section are the places where you can make your changes." << "\n";
    f_out << "--> Each line of the DEF sections that starts with a '*' is one definition. The string following the '*'" << "\n";
    f_out << "    is the internal atom type. The next string seperated by a whitespace (or more) represents" << "\n";
    f_out << "    the type the internal type will be set to for your purpose. The rest of the line will" << "\n";
    f_out << "    be ignored (so you can make additional comments)." << "\n" << "\n";
    f_out << "--> A short explanation of the flags from the FLAGS section:\n";
    f_out << "    (Changes in the following explanations have no effect. Make your changes in the FLAGS section!)\n";
    f_out << "---------------------------------------------------" << "\n";
    f_out << "Please use '0' for false and '1' for true. If you do not want to override values given\n"
          << "in a function call to MOLECULE::get_atom_typing, please use '-1'\n\n";
    f_out << "protonate_acids 0         = acids with unknown protonation state will be typed charged\n"
          << "                            (if they are protonated in the input molecule, they keep uncharged)\n";
    f_out << "protonate_guanidin 1      = guanidino groups will be typed protonated (charged)\n";
    f_out << "protonate_amidin 1        = amidino groups will be protonated (charged)\n";
    f_out << "protonate_amine -1        = protonation state for amines is used as given by the\n"
          << "                            funtion call to get_atom_typing (e.g. if you specified\n"
          << "                            protonation states in fconv using '--p' the specified\n"
          << "                            value for amines will be used (while 0 or 1 would\n"
          << "                            override your fconv settings).";
    f_out << "protonate_phosphate 0     = HxPOy will be typed charged\n";
    f_out << "protonate_sulfate 0       = HxSOy and SO2NH will be typed charged\n";
    f_out << "set_bonds 1               = bond types will be assigned\n";
    f_out << "kekulize_aromatics 0      = use bond type 'ar' instead of alternating double bonds in aromatic systems\n";
//    f_out << "kekulize_charged 0        = use bond type 'ar' for charged systems like COO- or SO2NH-\n";
    f_out << "allow_charged_aromatics 1 = if possible, non aromatic systems will be made aromatic charging a nitrogen\n";
    f_out << "max_ring_members 10       = A ring with more than this number of members will not be considered a ring\n";
    f_out << "\n";
    f_out << "--> A short explanation of the internal atom types: " << "\n";
    f_out << "---------------------------------------------------" << "\n";
    f_out << "--> All elements having no type listed in the DEF section will be typed by the element name." << "\n";
    f_out << "--> The order of atom types also represents the priority if there is more than one" << "\n";
    f_out << "    type possible (e.g. a guanidino N that is also bonded to a carbonyl group will be" << "\n";
    f_out << "    set to a N.am type, because this type has a higher priority than N.gu types)." << "\n";
    f_out << "--> Heteroaromatics with unknown protonation state will be typed in their neutral form." << "\n";
    f_out << "--> Enoles without explicitly set hydrogens will be typed as ketones." << "\n";
    f_out << "--> Currently there are: (here comes just the explanation, NOT the definition)" << "\n";
    f_out << "H.ac   = acidic H (bonded to O.3ac, N.im, N.sam or N.ohac)" << "\n";
    f_out << "H.onh  = amide NH" << "\n";
    f_out << "H.n    = bonded to other nitrogens" << "\n";
    f_out << "H.o    = bonded to other oxygens" << "\n";
    f_out << "H.0    = all other hydrogens" << "\n" << "\n";
    f_out << "C.ar6p = sp2 carbon with a positive charged resonance structure in a protonated 6-membered heteroaromatic ring" << "\n";
    f_out << "C.ar6x = sp2 carbon in a 6-membered heteroaromatic ring" << "\n";
    f_out << "C.ar6  = sp2 carbon in a benzene ring" << "\n";
    f_out << "C.arp  = sp2 carbon with a positive charged resonance structure in other protonated heteroaromatic rings" << "\n";
    f_out << "C.arx  = sp2 carbon in other heteroaromatics" << "\n";
    f_out << "C.ar   = sp2 carbon in other aromatics" << "\n";
    f_out << "C.2r3o = carbonyl carbon in cyclopropanone or cyclopropenone" << "\n";
    f_out << "C.2r3x = sp2 carbon in heterocyclic 3-membered rings" << "\n";
    f_out << "C.2r3  = sp2 carbon in 3-membered rings" << "\n";
    f_out << "C.3r3x = sp3 carbon in heterocyclic 3-membered rings" << "\n";
    f_out << "C.3r3  = sp3 carbon in 3-membered rings" << "\n";
    f_out << "C.1n   = sp carbon in cyano groups" << "\n";
    f_out << "C.1p   = sp carbon with one heavy atom bonded" << "\n";
    f_out << "C.1s   = sp carbon with two heavy atoms bonded" << "\n";
    f_out << "C.co2h = sp2 carbon in explicitly protonated COOH groups" << "\n";
    f_out << "C.co2  = sp2 carbon in COO-  groups (also set if protonation state is unknown)" << "\n";
    f_out << "C.es   = carbonyl carbon in ester groups or anhydrides" << "\n";
    f_out << "C.hal  = carbonyl carbon in acidhalogenides" << "\n";
    f_out << "C.am   = carbonyl carbon in amides" << "\n";
    f_out << "C.o    = other carbonyl carbon" << "\n";
    f_out << "C.s    = thionyl carbon" << "\n";
    f_out << "C.gu   = sp2 carbon in unprotonated guanidino groups" << "\n";
    f_out << "C.guh  = sp2 carbon in protonated guanidino groups (also set if protonation state is unknown)" << "\n";
    f_out << "C.mi   = sp2 carbon in unprotonated amidino groups" << "\n";
    f_out << "C.mih  = sp2 carbon in protonated amidino groups (also set if protonation state is unknown)" << "\n";
    f_out << "C.n    = sp2 carbon in imines" << "\n";
    f_out << "C.2p   = other sp2 carbon with one heavy atom bonded" << "\n";
    f_out << "C.2s   = other sp2 carbon with two heavy atoms bonded" << "\n";
    f_out << "C.2t   = other sp2 carbon with 3 heavy atoms bonded" << "\n";
    f_out << "C.et   = sp3 carbon in ethers" << "\n";
    f_out << "C.ohp  = sp3 carbon in primary alcoholes" << "\n";
    f_out << "C.ohs  = sp3 carbon in secondary alcoholes" << "\n";
    f_out << "C.oht  = sp3 carbon in tertiary alcoholes" << "\n";
    f_out << "C.3n   = other sp3 carbon bonded to nitrogen" << "\n";
    f_out << "C.3p   = other sp3 carbon with one heavy atom bonded" << "\n";
    f_out << "C.3s   = other sp3 carbon with two heavy atoms bonded" << "\n";
    f_out << "C.3t   = other sp3 carbon with 3 heavy atoms bonded" << "\n";
    f_out << "C.3q   = other sp3 carbon with 4 heavy atoms bonded" << "\n" << "\n";
    f_out << "N.ar6p = positive charged nitrogen in 6-membered aromatics (e.g. pyridinium or NAD+)" << "\n";
    f_out << "N.ar6  = sp2 nitrogen in 6-membered aromatics" << "\n";
    f_out << "N.arp  = sp2 nitrogen in protonated aromatics (e.g both nitrogens in protonated imidazole" << "\n";
    f_out << "N.ar2  = sp2 nitrogen in aromatics with two bonded atoms (corresponding to sybyl type N.2)" << "\n";
    f_out << "N.ar3  = sp2 nitrogen in aromatics with 3 heavy atoms (corresponding to sybyl type N.pl3)" << "\n";
    f_out << "N.ar3h = sp2 nitrogen in aromatics with 2 heavy atoms and one hydrogen (corresponding to sybyl type N.pl3)" << "\n";
    f_out << "N.r3   = sp3 in aziridine or azirene rings" << "\n";
    f_out << "N.az   = middle nitrogen in azides" << "\n";
    f_out << "N.1    = other sp nitrogen" << "\n";
    f_out << "N.o2   = in nitro groups" << "\n";
    f_out << "N.ohac = in hydroxamic acids" << "\n";
    f_out << "N.oh   = in hydroxylamines" << "\n";
    f_out << "N.ims  = imide nitrogen with two heavy atoms bonded" << "\n";
    f_out << "N.imt  = imide nitrogen with 3 heavy atoms bonded" << "\n";
    f_out << "N.amp  = carbon- or thionamide with one heavy atom bonded" << "\n";
    f_out << "N.ams  = carbon- or thionamide with two heavy atoms bonded" << "\n";
    f_out << "N.amt  = carbon- or thionamide with 3 heavy atoms bonded" << "\n";
    f_out << "N.samp = sulfonamide with one heavy atom bonded" << "\n";
    f_out << "N.sams = sulfonamide with two heavy atoms bonded" << "\n";
    f_out << "N.samt = sulfonamide with 3 heavy atoms bonded" << "\n";
    f_out << "N.gu1  = NH in unprotonated guanidino group (only if explicitly protonated)" << "\n";
    f_out << "N.gu2  = NH2 in unprotonated guanidino group (only if explicitly protonated)" << "\n";
    f_out << "N.guh  = nitrogen in protonated guanidino group (also set if protonation state is unknown)" << "\n";
    f_out << "N.mi1  = NH in unprotonated amidino group (only if explicitly protonated)" << "\n";
    f_out << "N.mi2  = NH2 in unprotonated amidino group (only if explicitly protonated)" << "\n";
    f_out << "N.mih  = nitrogen in protonated amidino group (also set if protonation state is unknown)" << "\n";
    f_out << "N.aap  = primary aromatic amine (hybridization can't be determined exactly)" << "\n";
    f_out << "N.aas2 = sp2 hybridized secondary aromatic amine" << "\n";
    f_out << "N.aas3 = sp3 hybridized secondary aromatic amine" << "\n";
    f_out << "N.aat2 = sp2 hybridized tertiary aromatic amine" << "\n";
    f_out << "N.aat3 = sp3 hybridized tertiary aromatic amine" << "\n";
    f_out << "N.2n   = sp2 nitrogen bonded to another nitrogen" << "\n";
    f_out << "N.2p   = other sp2 nitrogen with one heavy atom" << "\n";
    f_out << "N.2s   = other sp2 nitrogen with two heavy atoms" << "\n";
    f_out << "N.2t   = other sp2 nitrogen with three heavy atoms" << "\n";
    f_out << "N.3n   = sp3 nitrogen bonded to another nitrogen" << "\n";
    f_out << "N.3p   = sp3 nitrogen with one heavy atom bonded" << "\n";
    f_out << "N.3s   = sp3 nitrogen with two heavy atoms bonded" << "\n";
    f_out << "N.3t   = sp3 nitrogen with 3 heavy atoms bonded" << "\n";
    f_out << "N.4q   = sp3 nitrogen with 4 bonded heavy atoms" << "\n";
    f_out << "N.4h   = sp3 nitrogen with 4 bonded atoms (at least 1 hydrogen)\n" << "\n";
    f_out << "O.ar   = aromatic oxygen" << "\n";
    f_out << "O.r3   = in oxiran ring" << "\n";
    f_out << "O.h2o  = water oxygen" << "\n";
    f_out << "O.n    = oxygen in nitro groups" << "\n";
    f_out << "O.noh  = sp3 oxygen in hydroxylamine or hydroxamic acid" << "\n";
    f_out << "O.2co2 = sp2 oxygen in COOH (sp2 bonded to C.co2h)" << "\n";
    f_out << "O.2es  = sp2 oxygen in esters or anhydrids" << "\n";
    f_out << "O.2hal = sp2 oxygen in acidhalogenides" << "\n";
    f_out << "O.am   = in carbonamides" << "\n";
    f_out << "O.co2  = in COO-  or CSO-" << "\n";
    f_out << "O.2po  = sp2 oxygen in P=O (non deprotonated groups)" << "\n";
    f_out << "O.2so  = sp2 oxygen in S=O (non deprotonated groups)" << "\n";
    f_out << "O.2p   = sp2 oxygen in OPO3H- or PO3H- or POO-" << "\n";
    f_out << "O.2s   = sp2 oxygen in OSO3- or SO3- or POO- or deprotonated sulfonamides" << "\n";
    f_out << "O.3po  = sp3 oxygen with 2 heavy atoms bonded to at least one phosphor" << "\n";
    f_out << "O.3so  = sp3 oxygen with 2 heavy atoms bonded to at least one sulfur" << "\n";
    f_out << "O.carb = in other carbonyl groups" << "\n";
    f_out << "O.o    = in peroxo groups" << "\n";
    f_out << "O.3ac  = OH in COOH, CSOH, POOHOH, POOH or SOOOH" << "\n";
    f_out << "O.ph   = phenolic hydroxyl group" << "\n";
    f_out << "O.3oh  = hydroxyl group" << "\n";
    f_out << "O.3es  = sp3 oxygen in esters or anhydrids" << "\n";
    f_out << "O.3eta = aromatic ether" << "\n";
    f_out << "O.3et  = aliphatic ether" << "\n" << "\n";
    f_out << "S.ar   = aromatic sulfur" << "\n";
    f_out << "S.r3   = in thiiran ring" << "\n";
    f_out << "S.thi  = thionyl group" << "\n";
    f_out << "S.o    = in SO" << "\n";
    f_out << "S.o2h  = in protonated sulfonamide or other SO2" << "\n";
    f_out << "S.o3h  = in SO3" << "\n";
    f_out << "S.o4h  = in OSO3" << "\n";
    f_out << "S.o2   = in SO2 or deprotonated sulfonamides (or unknown protonation state)" << "\n";
    f_out << "S.o3   = in SO3- (or unknown protonation state)" << "\n";
    f_out << "S.o4   = in OSO3- (or unknown protonation state)" << "\n";
    f_out << "S.2    = in CSO-  COS-  or other sp2" << "\n";
    f_out << "S.sh   = in SH groups" << "\n";
    f_out << "S.s    = in S-S bonds" << "\n";
    f_out << "S.3    = other sp3 sulfur" << "\n" << "\n";
    f_out << "P.r3   = in phosphiran rings" << "\n";
    f_out << "P.o    = in PO" << "\n";
    f_out << "P.o2h  = in not deprotonated PO2 groups" << "\n";
    f_out << "P.o3h  = in not deprotonated PO3 groups" << "\n";
    f_out << "P.o4h  = in not deprotonated PO4 groups" << "\n";
    f_out << "P.o2   = in deprotonated PO2 groups (or unknown protonation state)" << "\n";
    f_out << "P.o3   = in deprotonated PO3 groups (or unknown protonation state)" << "\n";
    f_out << "P.o4   = in deprotonated PO4 groups (or unknown protonation state)" << "\n";
    f_out << "P.3    = other sp3" << "\n";
    f_out << "F.0    = bonded fluor\n";
    f_out << "F.i    = fluor ion\n";
    f_out << "Cl.0   = bonded chlorine\n";
    f_out << "Cl.i   = chlorine ion\n";
    f_out << "Br.0   = bonded bromine\n";
    f_out << "Br.i   = bromine ion\n";
    f_out << "I.0    = bonded iod\n";
    f_out << "I.i    = iod ion\n";
    

    f_out << "------------------------------------------------------------------------------------------------------------- \n";
    f_out << "If you need more differentiation, another priority order or if you have just any other" << "\n";
    f_out << "suggestions, please contact me:" << "\n";
    f_out << " Gerd Neudert" << "\n";
    f_out << " neudert@staff.uni-marburg.de" << "\n";

    f_out << "\n";
    
    //!Jetzt der Header mit Versionsstring und Anzahl der internen Typen:
    f_out << "<HEADER>" << "\n";
    f_out << "internal_set_type " << a_t_version << "\n";
    f_out << "total " << n_intern_types << "\n" << "\n";

    f_out << "<FLAGS>\n";
    f_out << "protonate_acids " << prot_acids << "\n";
    f_out << "protonate_guanidin " << prot_guanidin << "\n";
    f_out << "protonate_amidin " << prot_amidin << "\n";
    f_out << "protonate_amine " << prot_amin << "\n";
    f_out << "protonate_phosphate " << prot_phosphate << "\n";
    f_out << "protonate_sulfate " << prot_sulfate << "\n";
    f_out << "kekulize_aromatics " << kekulize_aromatics << "\n";
//    f_out << "kekulize_charged " << kekulize_charged << "\n";
    f_out << "allow_charged_aromatics " << allow_charged_aromatics << "\n";
    f_out << "set_bonds " << get_bonds << "\n";
    f_out << "max_ring_members " << max_ring_members << "\n";
    
    //!Jetzt die eigentlichen Definitionen:
    f_out << "<DEF>" << "\n";
    for (int i=0; i<n_intern_types; ++i) {
        f_out << "* ";
        f_out.width(8); f_out << left << i_t[i];
        if (alt_mode) f_out << (*alt_mode)[i] << " \n";
        else if (mode == 0) f_out << i_t[i] << " \n";
        else if (mode == 1) f_out << mode_1[i] << " \n";
        else if (mode == 2) f_out << mode_2[i] << " \n";
    }
    f_out.close();
    f_out.clear();
}


void write_def_file(int const& mode,const char *name,vector<string> *alt_mode,
                    tr1::unordered_map<string,string> new_flags) {
    ofstream f_out;
    f_out.open(name);

    //!zunaechst den Kommentar schreiben:
    f_out << "This is an atom type definitions file from 'molecule_GN.cpp'" << "\n" << "\n";
    f_out << "<COMMENT>" << "\n";
    f_out << "--> Everything in the COMMENT sections will NOT be read by the program (it's just comments)." << "\n";
    f_out << "--> In the HEADER section the 'internal_set_type' and the 'total' value is mandatory" << "\n";
    f_out << "    and should NEVER be changed by the user!" << "\n";
    f_out << "--> The FLAGS section and the DEF section are the places where you can make your changes." << "\n";
    f_out << "--> Each line of the DEF sections that starts with a '*' is one definition. The string following the '*'" << "\n";
    f_out << "    is the internal atom type. The next string seperated by a whitespace (or more) represents" << "\n";
    f_out << "    the type the internal type will be set to for your purpose. The rest of the line will" << "\n";
    f_out << "    be ignored (so you can make additional comments)." << "\n" << "\n";
    f_out << "--> A short explanation of the flags from the FLAGS section:\n";
    f_out << "    (Changes in the following explanations have no effect. Make your changes in the FLAGS section!)\n";
    f_out << "---------------------------------------------------" << "\n";
    f_out << "Please use '0' for false and '1' for true. If you do not want to override values given\n"
          << "in a function call to MOLECULE::get_atom_typing, please use '-1'\n\n";
    f_out << "protonate_acids 0         = acids with unknown protonation state will be typed charged\n"
          << "                            (if they are protonated in the input molecule, they keep uncharged)\n";
    f_out << "protonate_guanidin 1      = guanidino groups will be typed protonated (charged)\n";
    f_out << "protonate_amidin 1        = amidino groups will be protonated (charged)\n";
    f_out << "protonate_amine -1        = protonation state for amines is used as given by the\n"
          << "                            funtion call to get_atom_typing (e.g. if you specified\n"
          << "                            protonation states in fconv using '--p' the specified\n"
          << "                            value for amines will be used (while 0 or 1 would\n"
          << "                            override your fconv settings).";
    f_out << "protonate_phosphate 0     = HxPOy will be typed charged\n";
    f_out << "protonate_sulfate 0       = HxSOy and SO2NH will be typed charged\n";
    f_out << "set_bonds 1               = bond types will be assigned\n";
    f_out << "kekulize_aromatics 0      = use bond type 'ar' instead of alternating double bonds in aromatic systems\n";
//    f_out << "kekulize_charged 0        = use bond type 'ar' for charged systems like COO- or SO2NH-\n";
    f_out << "allow_charged_aromatics 1 = if possible, non aromatic systems will be made aromatic charging a nitrogen\n";
    f_out << "max_ring_members 10       = A ring with more than this number of members will not be considered a ring\n";
    f_out << "\n";
    f_out << "--> A short explanation of the internal atom types: " << "\n";
    f_out << "---------------------------------------------------" << "\n";
    f_out << "--> All elements having no type listed in the DEF section will be typed by the element name." << "\n";
    f_out << "--> The order of atom types also represents the priority if there is more than one" << "\n";
    f_out << "    type possible (e.g. a guanidino N that is also bonded to a carbonyl group will be" << "\n";
    f_out << "    set to a N.am type, because this type has a higher priority than N.gu types)." << "\n";
    f_out << "--> Heteroaromatics with unknown protonation state will be typed in their neutral form." << "\n";
    f_out << "--> Enoles without explicitly set hydrogens will be typed as ketones." << "\n";
    f_out << "--> Currently there are: (here comes just the explanation, NOT the definition)" << "\n";
    f_out << "H.ac   = acidic H (bonded to O.3ac, N.im, N.sam or N.ohac)" << "\n";
    f_out << "H.onh  = amide NH" << "\n";
    f_out << "H.n    = bonded to other nitrogens" << "\n";
    f_out << "H.o    = bonded to other oxygens" << "\n";
    f_out << "H.0    = all other hydrogens" << "\n" << "\n";
    f_out << "C.ar6p = sp2 carbon with a positive charged resonance structure in a protonated 6-membered heteroaromatic ring" << "\n";
    f_out << "C.ar6x = sp2 carbon in a 6-membered heteroaromatic ring" << "\n";
    f_out << "C.ar6  = sp2 carbon in a benzene ring" << "\n";
    f_out << "C.arp  = sp2 carbon with a positive charged resonance structure in other protonated heteroaromatic rings" << "\n";
    f_out << "C.arx  = sp2 carbon in other heteroaromatics" << "\n";
    f_out << "C.ar   = sp2 carbon in other aromatics" << "\n";
    f_out << "C.2r3o = carbonyl carbon in cyclopropanone or cyclopropenone" << "\n";
    f_out << "C.2r3x = sp2 carbon in heterocyclic 3-membered rings" << "\n";
    f_out << "C.2r3  = sp2 carbon in 3-membered rings" << "\n";
    f_out << "C.3r3x = sp3 carbon in heterocyclic 3-membered rings" << "\n";
    f_out << "C.3r3  = sp3 carbon in 3-membered rings" << "\n";
    f_out << "C.1n   = sp carbon in cyano groups" << "\n";
    f_out << "C.1p   = sp carbon with one heavy atom bonded" << "\n";
    f_out << "C.1s   = sp carbon with two heavy atoms bonded" << "\n";
    f_out << "C.co2h = sp2 carbon in explicitly protonated COOH groups" << "\n";
    f_out << "C.co2  = sp2 carbon in COO-  groups (also set if protonation state is unknown)" << "\n";
    f_out << "C.es   = carbonyl carbon in ester groups or anhydrides" << "\n";
    f_out << "C.hal  = carbonyl carbon in acidhalogenides" << "\n";
    f_out << "C.am   = carbonyl carbon in amides" << "\n";
    f_out << "C.o    = other carbonyl carbon" << "\n";
    f_out << "C.s    = thionyl carbon" << "\n";
    f_out << "C.gu   = sp2 carbon in unprotonated guanidino groups" << "\n";
    f_out << "C.guh  = sp2 carbon in protonated guanidino groups (also set if protonation state is unknown)" << "\n";
    f_out << "C.mi   = sp2 carbon in unprotonated amidino groups" << "\n";
    f_out << "C.mih  = sp2 carbon in protonated amidino groups (also set if protonation state is unknown)" << "\n";
    f_out << "C.n    = sp2 carbon in imines" << "\n";
    f_out << "C.2p   = other sp2 carbon with one heavy atom bonded" << "\n";
    f_out << "C.2s   = other sp2 carbon with two heavy atoms bonded" << "\n";
    f_out << "C.2t   = other sp2 carbon with 3 heavy atoms bonded" << "\n";
    f_out << "C.et   = sp3 carbon in ethers" << "\n";
    f_out << "C.ohp  = sp3 carbon in primary alcoholes" << "\n";
    f_out << "C.ohs  = sp3 carbon in secondary alcoholes" << "\n";
    f_out << "C.oht  = sp3 carbon in tertiary alcoholes" << "\n";
    f_out << "C.3n   = other sp3 carbon bonded to nitrogen" << "\n";
    f_out << "C.3p   = other sp3 carbon with one heavy atom bonded" << "\n";
    f_out << "C.3s   = other sp3 carbon with two heavy atoms bonded" << "\n";
    f_out << "C.3t   = other sp3 carbon with 3 heavy atoms bonded" << "\n";
    f_out << "C.3q   = other sp3 carbon with 4 heavy atoms bonded" << "\n" << "\n";
    f_out << "N.ar6p = positive charged nitrogen in 6-membered aromatics (e.g. pyridinium or NAD+)" << "\n";
    f_out << "N.ar6  = sp2 nitrogen in 6-membered aromatics" << "\n";
    f_out << "N.arp  = sp2 nitrogen in protonated aromatics (e.g both nitrogens in protonated imidazole" << "\n";
    f_out << "N.ar2  = sp2 nitrogen in aromatics with two bonded atoms (corresponding to sybyl type N.2)" << "\n";
    f_out << "N.ar3  = sp2 nitrogen in aromatics with 3 heavy atoms (corresponding to sybyl type N.pl3)" << "\n";
    f_out << "N.ar3h = sp2 nitrogen in aromatics with 2 heavy atoms and one hydrogen (corresponding to sybyl type N.pl3)" << "\n";
    f_out << "N.r3   = sp3 in aziridine or azirene rings" << "\n";
    f_out << "N.az   = middle nitrogen in azides" << "\n";
    f_out << "N.1    = other sp nitrogen" << "\n";
    f_out << "N.o2   = in nitro groups" << "\n";
    f_out << "N.ohac = in hydroxamic acids" << "\n";
    f_out << "N.oh   = in hydroxylamines" << "\n";
    f_out << "N.ims  = imide nitrogen with two heavy atoms bonded" << "\n";
    f_out << "N.imt  = imide nitrogen with 3 heavy atoms bonded" << "\n";
    f_out << "N.amp  = carbon- or thionamide with one heavy atom bonded" << "\n";
    f_out << "N.ams  = carbon- or thionamide with two heavy atoms bonded" << "\n";
    f_out << "N.amt  = carbon- or thionamide with 3 heavy atoms bonded" << "\n";
    f_out << "N.samp = sulfonamide with one heavy atom bonded" << "\n";
    f_out << "N.sams = sulfonamide with two heavy atoms bonded" << "\n";
    f_out << "N.samt = sulfonamide with 3 heavy atoms bonded" << "\n";
    f_out << "N.gu1  = NH in unprotonated guanidino group (only if explicitly protonated)" << "\n";
    f_out << "N.gu2  = NH2 in unprotonated guanidino group (only if explicitly protonated)" << "\n";
    f_out << "N.guh  = nitrogen in protonated guanidino group (also set if protonation state is unknown)" << "\n";
    f_out << "N.mi1  = NH in unprotonated amidino group (only if explicitly protonated)" << "\n";
    f_out << "N.mi2  = NH2 in unprotonated amidino group (only if explicitly protonated)" << "\n";
    f_out << "N.mih  = nitrogen in protonated amidino group (also set if protonation state is unknown)" << "\n";
    f_out << "N.aap  = primary aromatic amine (hybridization can't be determined exactly)" << "\n";
    f_out << "N.aas2 = sp2 hybridized secondary aromatic amine" << "\n";
    f_out << "N.aas3 = sp3 hybridized secondary aromatic amine" << "\n";
    f_out << "N.aat2 = sp2 hybridized tertiary aromatic amine" << "\n";
    f_out << "N.aat3 = sp3 hybridized tertiary aromatic amine" << "\n";
    f_out << "N.2n   = sp2 nitrogen bonded to another nitrogen" << "\n";
    f_out << "N.2p   = other sp2 nitrogen with one heavy atom" << "\n";
    f_out << "N.2s   = other sp2 nitrogen with two heavy atoms" << "\n";
    f_out << "N.2t   = other sp2 nitrogen with three heavy atoms" << "\n";
    f_out << "N.3n   = sp3 nitrogen bonded to another nitrogen" << "\n";
    f_out << "N.3p   = sp3 nitrogen with one heavy atom bonded" << "\n";
    f_out << "N.3s   = sp3 nitrogen with two heavy atoms bonded" << "\n";
    f_out << "N.3t   = sp3 nitrogen with 3 heavy atoms bonded" << "\n";
    f_out << "N.4q   = sp3 nitrogen with 4 bonded heavy atoms" << "\n";
    f_out << "N.4h   = sp3 nitrogen with 4 bonded atoms (at least 1 hydrogen)\n" << "\n";
    f_out << "O.ar   = aromatic oxygen" << "\n";
    f_out << "O.r3   = in oxiran ring" << "\n";
    f_out << "O.h2o  = water oxygen" << "\n";
    f_out << "O.n    = oxygen in nitro groups" << "\n";
    f_out << "O.noh  = sp3 oxygen in hydroxylamine or hydroxamic acid" << "\n";
    f_out << "O.2co2 = sp2 oxygen in COOH (sp2 bonded to C.co2h)" << "\n";
    f_out << "O.2es  = sp2 oxygen in esters or anhydrids" << "\n";
    f_out << "O.2hal = sp2 oxygen in acidhalogenides" << "\n";
    f_out << "O.am   = in carbonamides" << "\n";
    f_out << "O.co2  = in COO-  or CSO-" << "\n";
    f_out << "O.2po  = sp2 oxygen in P=O (non deprotonated groups)" << "\n";
    f_out << "O.2so  = sp2 oxygen in S=O (non deprotonated groups)" << "\n";
    f_out << "O.2p   = sp2 oxygen in OPO3H- or PO3H- or POO-" << "\n";
    f_out << "O.2s   = sp2 oxygen in OSO3- or SO3- or POO- or deprotonated sulfonamides" << "\n";
    f_out << "O.3po  = sp3 oxygen with 2 heavy atoms bonded to at least one phosphor" << "\n";
    f_out << "O.3so  = sp3 oxygen with 2 heavy atoms bonded to at least one sulfur" << "\n";
    f_out << "O.carb = in other carbonyl groups" << "\n";
    f_out << "O.o    = in peroxo groups" << "\n";
    f_out << "O.3ac  = OH in COOH, CSOH, POOHOH, POOH or SOOOH" << "\n";
    f_out << "O.ph   = phenolic hydroxyl group" << "\n";
    f_out << "O.3oh  = hydroxyl group" << "\n";
    f_out << "O.3es  = sp3 oxygen in esters or anhydrids" << "\n";
    f_out << "O.3eta = aromatic ether" << "\n";
    f_out << "O.3et  = aliphatic ether" << "\n" << "\n";
    f_out << "S.ar   = aromatic sulfur" << "\n";
    f_out << "S.r3   = in thiiran ring" << "\n";
    f_out << "S.thi  = thionyl group" << "\n";
    f_out << "S.o    = in SO" << "\n";
    f_out << "S.o2h  = in protonated sulfonamide or other SO2" << "\n";
    f_out << "S.o3h  = in SO3" << "\n";
    f_out << "S.o4h  = in OSO3" << "\n";
    f_out << "S.o2   = in SO2 or deprotonated sulfonamides (or unknown protonation state)" << "\n";
    f_out << "S.o3   = in SO3- (or unknown protonation state)" << "\n";
    f_out << "S.o4   = in OSO3- (or unknown protonation state)" << "\n";
    f_out << "S.2    = in CSO-  COS-  or other sp2" << "\n";
    f_out << "S.sh   = in SH groups" << "\n";
    f_out << "S.s    = in S-S bonds" << "\n";
    f_out << "S.3    = other sp3 sulfur" << "\n" << "\n";
    f_out << "P.r3   = in phosphiran rings" << "\n";
    f_out << "P.o    = in PO" << "\n";
    f_out << "P.o2h  = in not deprotonated PO2 groups" << "\n";
    f_out << "P.o3h  = in not deprotonated PO3 groups" << "\n";
    f_out << "P.o4h  = in not deprotonated PO4 groups" << "\n";
    f_out << "P.o2   = in deprotonated PO2 groups (or unknown protonation state)" << "\n";
    f_out << "P.o3   = in deprotonated PO3 groups (or unknown protonation state)" << "\n";
    f_out << "P.o4   = in deprotonated PO4 groups (or unknown protonation state)" << "\n";
    f_out << "P.3    = other sp3" << "\n";
    f_out << "F.0    = bonded fluor\n";
    f_out << "F.i    = fluor ion\n";
    f_out << "Cl.0   = bonded chlorine\n";
    f_out << "Cl.i   = chlorine ion\n";
    f_out << "Br.0   = bonded bromine\n";
    f_out << "Br.i   = bromine ion\n";
    f_out << "I.0    = bonded iod\n";
    f_out << "I.i    = iod ion\n";


    f_out << "------------------------------------------------------------------------------------------------------------- \n";
    f_out << "If you need more differentiation, another priority order or if you have just any other" << "\n";
    f_out << "suggestions, please contact me:" << "\n";
    f_out << " Gerd Neudert" << "\n";
    f_out << " neudert@staff.uni-marburg.de" << "\n";

    f_out << "\n";

    //!Jetzt der Header mit Versionsstring und Anzahl der internen Typen:
    f_out << "<HEADER>" << "\n";
    f_out << "internal_set_type " << a_t_version << "\n";
    f_out << "total " << n_intern_types << "\n" << "\n";

    f_out << "<FLAGS>\n";

    for (tr1::unordered_map<string,string>::iterator it=new_flags.begin(); it!=new_flags.end(); ++it) {
        f_out << it->first << " " << it->second << "\n";
    }

    //!Jetzt die eigentlichen Definitionen:
    f_out << "<DEF>" << "\n";
    for (int i=0; i<n_intern_types; ++i) {
        f_out << "* ";
        f_out.width(8); f_out << left << i_t[i];
        if (alt_mode) f_out << (*alt_mode)[i] << " \n";
        else if (mode == 0) f_out << i_t[i] << " \n";
        else if (mode == 1) f_out << mode_1[i] << " \n";
        else if (mode == 2) f_out << mode_2[i] << " \n";
    }
    f_out.close();
    f_out.clear();
}


void load_def_file_global(char const* name,
                          bool &prot_acids,
                          bool &prot_guanidin,
                          bool &prot_amidin,
                          bool &prot_amin,
                          bool &prot_phosphate,
                          bool &prot_sulfate,
                          bool &get_bonds,
                          int &max_ring_members,
                          bool kekulize_aromatics,
                          bool kekulize_charged,
                          bool allow_charged_aromatics) {
    //! ACHTUNG: laedt immer in a_t_map
    istringstream is;
    string row,key,vers,tot;
    ifstream f_in;
    f_in.open(name);
    if (!f_in) {
        cerr << c_message<cERROR>("could not open ") << name << endl;
        exit(1);
    }
    while (!f_in.eof()) {
        getline(f_in,row);
        is.clear(); is.str(row);
        is >> key;
        if (key == "<HEADER>") break;
    }
    while (!f_in.eof()) {
        getline(f_in,row);
        is.clear(); is.str(row);
        is >> key;
        if (key == "internal_set_type") is >> vers;
        else if (key == "total") {is >> tot; break;}
    }
    if (vers != a_t_version) {
        f_in.close();
        cerr << c_message<cERROR>("load_def_file_global --> your definition file has a wrong version number") << endl;
        exit(1);
    }
    int n_types = 0;
    bool in_flags = false;
    while (!f_in.eof()) {
        getline(f_in,row);
        is.clear(); is.str(row);
        is >> key;
        if (key == "*") {
            is >> vers;
            is >> tot;
            def_map::a_t_map[vers] = tot;
            ++n_types;
            vers = "UNK";
            tot = "UNK";
            key = "X";
        } else if (key == "<FLAGS>") in_flags = true;
        else if (in_flags) {
            if (key == "protonate_acids") {
                int ti = -1;
                is >> ti;
                if (ti == 0) prot_acids = false;
                else if (ti > 0) prot_acids = true;
                else continue;
            } else if (key == "protonate_amidin") {
                int ti = -1;
                is >> ti;
                if (ti == 0) prot_amidin = false;
                else if (ti > 0) prot_amidin = true;
                else continue;
            } else if (key == "protonate_guanidin") {
                int ti = -1;
                is >> ti;
                if (ti == 0) prot_guanidin = false;
                else if (ti > 0) prot_guanidin = true;
                else continue;
            } else if (key == "protonate_amine") {
                int ti = -1;
                is >> ti;
                if (ti == 0) prot_amin = false;
                else if (ti > 0) prot_amin = true;
                else continue;
            } else if (key == "protonate_phosphate") {
                int ti = -1;
                is >> ti;
                if (ti == 0) prot_phosphate = false;
                else if (ti > 0) prot_phosphate = true;
                else continue;
            } else if (key == "protonate_sulfate") {
                int ti = -1;
                is >> ti;
                if (ti == 0) prot_sulfate = false;
                else if (ti > 0) prot_sulfate = true;
                else continue;
            } else if (key == "set_bonds") {
                int ti = -1;
                is >> ti;
                if (ti == 0) get_bonds = false;
                else if (ti > 0) get_bonds = true;
                else continue;
            } else if (key == "kekulize_aromatics") {
                int ti = -1;
                is >> ti;
                if (ti == 0) kekulize_aromatics = false;
                else if (ti > 0) kekulize_aromatics = true;
                else continue;
//            } else if (key == "kekulize_charged") {
//                int ti = -1;
//                is >> ti;
//                if (ti == 0) kekulize_charged = false;
//                else if (ti > 0) kekulize_charged = true;
//                else continue;
            } else if (key == "allow_charged_aromatics") {
                int ti = -1;
                is >> ti;
                if (ti == 0) allow_charged_aromatics = false;
                else if (ti > 0) allow_charged_aromatics = true;
                else continue;
            } else if (key == "max_ring_members") {
                int ti = -1;
                is >> ti;
                if (ti < 1) continue;
                else max_ring_members = ti;
            }
        }
    }
    def_map::last_a_t_map = 1;
    if (n_intern_types != n_types) {
        cerr << c_message<cWARNING>("load_def_file_global --> read ") << n_types << " atom type definitions from " << name << " , but " << n_intern_types << " were expected!" << endl;
    }
    f_in.close();
}

