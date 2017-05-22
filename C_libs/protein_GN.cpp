
//============================================================================
// protein_GN.cpp -*- C++ -*-; representation of molecular protein data
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
// This library implements protein objects, which are used
// as high-level abstraction for molecular data
//============================================================================


#include"protein_GN.h"
#include"structure_GN.h"

//==============================================================================================
//Definitionen fuer PROTEIN:
//==============================================================================================

PROTEIN::PROTEIN() : last_name("XXX"),/*n_aacids(0),*/bnd_id(1) {}

PROTEIN::~PROTEIN() {
    for (chains_vec it=chains.begin(); it!=chains.end();++it) {
        it->kill();
    }
}

void PROTEIN::get_atom_types(int mode,const char *def_file,bool fill_X) {
    for (chains_vec it=chains.begin(); it!=chains.end();++it) (*it)->get_atom_types(mode,def_file,fill_X);
}

void PROTEIN::get_aacids() { //!dringend noch kommentieren
    last_name = "XXX";
//    n_aacids = 0;
    for (chains_vec it=chains.begin(); it!=chains.end();++it) (*it)->get_aacids();
}

void PROTEIN::build_bonds() { //nicht schoen, aber laeuft
    bnd_id = 1;
    for (chains_vec jt=chains.begin(); jt!=chains.end();++jt) (*jt)->build_bonds();
}

void PROTEIN::ter_correct() { //!wenn diese Methode aufgerufen wurde kann nicht mehr als pdb gespeichert werden!!!
    //Die internen id's muessen umgesetzt werden, wenn vom pdb zum mol2 uebergegangen wird, weil die TER-Eintraege wegfallen
    //Entsprechend muessen natrlich auch die CONECT-Eintraege geaendert werden
    //!h_log kann nicht geaendert werden, weil es zu parser gehoert => keine gluetigkeit mehr als pdb
    int dec = 0;
    
    for (chains_vec it=chains.begin(); it!=chains.end();++it) {
        //Fr alle Atome der aktuellen Kette:
        for (atoms_vec jt=(*it)->atoms.begin(); jt!=(*it)->atoms.end();++jt) {
            //Korrektur der CONECT-Eintraege:
            for (conects_vec kt=main_structure->conects.begin(); kt!=main_structure->conects.end();++kt) {
                if ((*kt)->from_id == (*jt)->intern_id) {
                    (*kt)->from_id -= dec;
                }
                //Fr alle kovalenten Bindungen:
                for (int i=0; i<4; ++i) {
                    if ((*kt)->cov_bonds[i] == (*jt)->intern_id) {
                        (*kt)->cov_bonds[i] -= dec;
                    }
                }
            }
            (*jt)->intern_id -= dec;
            if ((*jt)->is_ter) dec++; //ab jetzt alles um 1 mehr erniedrigen
        }
    }
}

void PROTEIN::get_reli_centers() {
    for (chains_vec it=chains.begin(); it!=chains.end();++it) (*it)->get_reli_centers();
}

void PROTEIN::visualize_reli_centers() {
    ofstream f_out;
    f_out.open("pseudocenters_vis.py");

    f_out << "# This visualization file was automatically created by 'structure_GN.cpp'" << "\n";
    f_out << "# Start pymol and use 'run this_files_name.py' to start the visualization" << "\n" << "\n";
    
    float rad = 0.3;
    
    
    bool first = true;
    f_out << "PI = [";
    for (chains_vec it=chains.begin(); it!=chains.end();++it) {
        for (reli_centers_vec jt=(*it)->reli_centers.begin(); jt!=(*it)->reli_centers.end(); ++jt) {
            if ((*jt)->type != "PI") continue;
            if (!first) f_out << ",";
            else first = false;
            f_out << "7.0," << (*jt)->coord[0] << "," << (*jt)->coord[1] << "," << (*jt)->coord[2] << "," << rad;
        }
    }
    f_out << "]" << endl;
    f_out << "cmd.load_cgo(PI,'PI', 1)" << endl;
    f_out << "cmd.color('yellow','PI')" << "\n";
    
    
    first = true;
    f_out << "Aliphatic = [";
    for (chains_vec it=chains.begin(); it!=chains.end();++it) {
        for (reli_centers_vec jt=(*it)->reli_centers.begin(); jt!=(*it)->reli_centers.end(); ++jt) {
            if ((*jt)->type != "Aliphatic") continue;
            if (!first) f_out << ",";
            else first = false;
            f_out << "7.0," << (*jt)->coord[0] << "," << (*jt)->coord[1] << "," << (*jt)->coord[2] << "," << rad;
        }
    }
    f_out << "]" << endl;
    f_out << "cmd.load_cgo(Aliphatic,'Aliphatic', 1)" << endl;
    f_out << "cmd.color('white','Aliphatic')" << "\n";
    
    
    first = true;
    f_out << "DON_ACC = [";
    for (chains_vec it=chains.begin(); it!=chains.end();++it) {
        for (reli_centers_vec jt=(*it)->reli_centers.begin(); jt!=(*it)->reli_centers.end(); ++jt) {
            if ((*jt)->type != "DON_ACC") continue;
            if (!first) f_out << ",";
            else first = false;
            f_out << "7.0," << (*jt)->coord[0] << "," << (*jt)->coord[1] << "," << (*jt)->coord[2] << "," << rad;
        }
    }
    f_out << "]" << endl;
    f_out << "cmd.load_cgo(DON_ACC,'DON_ACC', 1)" << endl;
    f_out << "cmd.color('green','DON_ACC')" << "\n";
    
    
    first = true;
    f_out << "Donor = [";
    for (chains_vec it=chains.begin(); it!=chains.end();++it) {
        for (reli_centers_vec jt=(*it)->reli_centers.begin(); jt!=(*it)->reli_centers.end(); ++jt) {
            if ((*jt)->type != "Donor") continue;
            if (!first) f_out << ",";
            else first = false;
            f_out << "7.0," << (*jt)->coord[0] << "," << (*jt)->coord[1] << "," << (*jt)->coord[2] << "," << rad;
        }
    }
    f_out << "]" << endl;
    f_out << "cmd.load_cgo(Donor,'Donor', 1)" << endl;
    f_out << "cmd.color('red','Donor')" << "\n";
    
    
    first = true;
    f_out << "Acceptor = [";
    for (chains_vec it=chains.begin(); it!=chains.end();++it) {
        for (reli_centers_vec jt=(*it)->reli_centers.begin(); jt!=(*it)->reli_centers.end(); ++jt) {
            if ((*jt)->type != "Acceptor") continue;
            if (!first) f_out << ",";
            else first = false;
            f_out << "7.0," << (*jt)->coord[0] << "," << (*jt)->coord[1] << "," << (*jt)->coord[2] << "," << rad;
        }
    }
    f_out << "]" << endl;
    f_out << "cmd.load_cgo(Acceptor,'Acceptor', 1)" << endl;
    f_out << "cmd.color('blue','Acceptor')" << "\n";
    
    
    f_out.close();
    
    cout << "pymol debug visualization file written to  'pseudocenters_vis.py'" << endl;
}

void PROTEIN::align(PROTEIN *ref,bool debug,bool only_C_alpha) {
    //2 Proteinstrukturen ueberlagern
    //=> zunaechst Sequenzalignment, dann entsprechendes 3D-alignment
    //! Listen fuer das Sequenzalignment erstellen:
    vector<SA_CONTAINER<stl_ptr<AACID>,string> > Aseq;
    vector<SA_CONTAINER<stl_ptr<AACID>,string> > Bseq;
    
    SA_CONTAINER<stl_ptr<AACID>,string>  gap_val;
    AACID gap_at;
    gap_at.res_name = "XXX";
    stl_ptr<AACID> tatmp;
    tatmp = &gap_at;
    gap_val.at = tatmp;
    gap_val.comparer = "XXX";
    
    get_aacids();
    ref->get_aacids();
    
    int last_res_num = 1;
    int naA = 0;
    int naB = 0;
    for (chains_vec ct=chains.begin(); ct!=chains.end(); ++ct) {
        for (aacids_map aa=(*ct)->aacids.begin(); aa!=(*ct)->aacids.end(); ++aa) {
            if (abs(aa->second->res_number - last_res_num) > 1) Bseq.push_back(gap_val);
            last_res_num = aa->second->res_number;
            SA_CONTAINER<stl_ptr<AACID>,string> csa;
            csa.at = aa->second;
            csa.comparer = aa->second->res_name;
            ++naB;
            Bseq.push_back(csa);
        }
    }
    last_res_num = 1;
    for (chains_vec ct=ref->chains.begin(); ct!=ref->chains.end(); ++ct) {
        for (aacids_map aa=(*ct)->aacids.begin(); aa!=(*ct)->aacids.end(); ++aa) {
            if (abs(aa->second->res_number - last_res_num) > 1) Aseq.push_back(gap_val);
            last_res_num = aa->second->res_number;
            SA_CONTAINER<stl_ptr<AACID>,string> csa;
            csa.at = aa->second;
            csa.comparer = aa->second->res_name;
            ++naA;
            Aseq.push_back(csa);
        }
    }
    
    if (main_structure->verb > 0) {
        cout << " -> " << naA << " residues in reference structure" << "\n";
        cout << " -> " << naB << " residues in alignment structure" << "\n";
    }
    
    SEQALIGN<SA_CONTAINER<stl_ptr<AACID>,string> > my_sa(Aseq,Bseq,gap_val);
    my_sa.needleman_wunsch1();
    
    //! Listen fuer das geometrische Alignment erstellen:
    vector<vec3d<float> > A_list;
    vector<vec3d<float> > B_list;
    unsigned int matched = 0;
    for (unsigned int i=0; i<my_sa.ires.size(); ++i) {
        if (i >= my_sa.jres.size()) break;
        
        //!=== DEBUG =====
        if (debug) {
            cerr << my_sa.ires[i].at->res_name;
            if (my_sa.ires[i].at->res_name != "XXX") cerr << my_sa.ires[i].at->res_number;
            cerr << " <--> " << my_sa.jres[i].at->res_name;
            if (my_sa.jres[i].at->res_name != "XXX") cerr << my_sa.jres[i].at->res_number;
            cerr << endl;
        }
        //!===============
        
        if (my_sa.ires[i].at->res_name == my_sa.jres[i].at->res_name && 
            !(my_sa.ires[i].at->res_name == "XXX" && my_sa.jres[i].at->res_name == "XXX")) {
            ++matched;
            unsigned max_atm = my_sa.ires[i].at->atoms.size();
            if (my_sa.jres[i].at->atoms.size() < max_atm) max_atm = my_sa.jres[i].at->atoms.size();
            for (unsigned int ii=0; ii<max_atm; ++ii) {
                istringstream is;
                string name;
                is.str(my_sa.ires[i].at->atoms[ii]->name);
                is >> name;
                if (only_C_alpha && name != "CA") continue;
                
                A_list.push_back(my_sa.ires[i].at->atoms[ii]->coord);
                B_list.push_back(my_sa.jres[i].at->atoms[ii]->coord);
            }
        } else {
            
        }
    }
    
    if (main_structure->verb > 0) cout << " -> " << matched << " residues matched in sequence alignment" << "\n";
    
    get_align_matrix(A_list,B_list,optalign_rotm,optalign_trans[0],optalign_trans[1]);
    
    aligned_rmsd = 0.;
    vector<vec3d<float> > nA_list;
    vector<vec3d<float> > nB_list;
    for (unsigned int i=0; i<B_list.size(); ++i) {
        vec3d<float> tvb(B_list[i]);
        B_list[i] -= optalign_trans[1];
        B_list[i] *= optalign_rotm;
        B_list[i] += optalign_trans[0];
        float to_add = get_square_distance(B_list[i],A_list[i]);
        //! Neue Liste mit gematchten C-alphas machen, in die nur die jenigen Paare kommen,
        //! die nach diesem ersten Alignment schon nahe beieinander liegen:
        if (to_add > 6.) continue;
        nA_list.push_back(A_list[i]);
        nB_list.push_back(tvb);
        aligned_rmsd += to_add;
    }
    
    if (nB_list.size() > (matched / 2)) { //! Zweites Alignment nur, wenn nicht zuviele Paare rausgeworfen wurden!
        get_align_matrix(nA_list,nB_list,optalign_rotm,optalign_trans[0],optalign_trans[1]);
        aligned_rmsd = 0.;
        for (unsigned int i=0; i<nB_list.size(); ++i) {
            nB_list[i] -= optalign_trans[1];
            nB_list[i] *= optalign_rotm;
            nB_list[i] += optalign_trans[0];
            float to_add = get_square_distance(nB_list[i],nA_list[i]);
            aligned_rmsd += to_add;
        }
        aligned_rmsd /= nB_list.size();
        aligned_rmsd = sqrt(aligned_rmsd);
        matched = nB_list.size();
    } else {
        aligned_rmsd /= B_list.size();
        aligned_rmsd = sqrt(aligned_rmsd);
        matched = nB_list.size();
    }
    
    if (main_structure->verb > 0) cout << " -> rmsd = " << aligned_rmsd << "  (" << matched 
                                       << " C-alpha used in spatial alignment)" << "\n";
}


//##############################################################################
//##############################################################################
class C_AS {
public:
    stl_ptr<ATOM> ca;
    float dihed;
    vector<vec3d<float> > ref_system;
    C_AS(stl_ptr<ATOM> const& pN,stl_ptr<ATOM> const& pCA,stl_ptr<ATOM> const& pC,stl_ptr<ATOM> const& pO):ca(pCA) {
        dihed = dihedral(pN->coord,pCA->coord,pC->coord,pO->coord);
        vec3d<float> v1(pCA->coord); v1 -= pC->coord; v1.norm(); v1 *= 10.;
        ref_system.push_back(pC->coord+v1);
        vec3d<float> v2(pO->coord); v2 -= pC->coord; v2.norm(); v2 *= 10.;
        ref_system.push_back(pC->coord+v2);
        vec3d<float> v3(v1); v3 *= v2; v3.norm(); v3 *= 10.;
        ref_system.push_back(pC->coord+v3);
        ref_system.push_back(pC->coord);
    }
    ~C_AS() {}
};

class C_PAIR {
public:
    C_AS* aa;
    C_AS* ref_aa;
    vector<C_PAIR*> neighbours;
    C_PAIR(C_AS* a,C_AS* ref_a):aa(a),ref_aa(ref_a) {}
    ~C_PAIR() {}
};

class C_CLUSTER {
public:
    C_PAIR* p1;
    tr1::unordered_set<C_PAIR*> pairs;
    vector<C_PAIR*> new_pairs;
    vector<C_PAIR*> ref_points;
    C_CLUSTER(C_PAIR* cp):p1(cp) {pairs.insert(cp); new_pairs.push_back(cp); ref_points.push_back(cp);}
    ~C_CLUSTER() {}

    bool add_pair(C_PAIR* cp,float const& max_displace) {
        if (pairs.find(cp) != pairs.end()) return false;

        for (unsigned int i=0; i<4; ++i) {
            float dd = get_distance(p1->aa->ref_system[i],cp->aa->ref_system[i]);
            dd = fabs(dd - get_distance(p1->ref_aa->ref_system[i],cp->ref_aa->ref_system[i]));
            if (dd > max_displace) return false;
        }

        for (tr1::unordered_set<C_PAIR*>::iterator it=pairs.begin(); it!=pairs.end(); ++it) {
            if ((*it)->aa->ca->intern_id == cp->aa->ca->intern_id ||
                (*it)->ref_aa->ca->intern_id == cp->ref_aa->ca->intern_id) return false;
        }

        pairs.insert(cp);
        return true;
    }

    bool join(C_CLUSTER* const clq,float const& max_displace) {
        //! Pairs von clq zufuegen, wenn kein gemeinsames Pair und geometrisch kompatibel:
        for (tr1::unordered_set<C_PAIR*>::iterator it=pairs.begin(); it!=pairs.end(); ++it) {
            if (clq->pairs.find(*it) != clq->pairs.end()) return false;
        }

        for (unsigned int i=0; i<4; ++i) {
            float dd = get_distance(p1->aa->ref_system[i],clq->p1->aa->ref_system[i]);
            dd = fabs(dd - get_distance(p1->ref_aa->ref_system[i],clq->p1->ref_aa->ref_system[i]));
            if (dd > max_displace) return false;
        }
        for (tr1::unordered_set<C_PAIR*>::iterator it=pairs.begin(); it!=pairs.end(); ++it) {
            for (tr1::unordered_set<C_PAIR*>::iterator jt=clq->pairs.begin(); jt!=clq->pairs.end(); ++jt) {
                if ((*it)->aa->ca->intern_id == (*jt)->aa->ca->intern_id ||
                    (*it)->ref_aa->ca->intern_id == (*jt)->ref_aa->ca->intern_id) {
                    return false;
                }
            }
        }

        for (tr1::unordered_set<C_PAIR*>::iterator it=clq->pairs.begin(); it!=clq->pairs.end(); ++it) pairs.insert(*it);
        return true;
    }
};

float PROTEIN::align2(PROTEIN *ref,float const& max_displace,bool const& sim_only) {
    //! 1.) Relevante AS-Atome sortieren:
    vector<C_AS*> aa;
    vector<C_AS*> ref_aa;
    istringstream is;
    string name;
    for (chains_vec ct=chains.begin(); ct!=chains.end(); ++ct) {
        stl_ptr<ATOM> pN(0);
        stl_ptr<ATOM> pCA(0);
        stl_ptr<ATOM> pC(0);
        int last_num = -1;
        for (atoms_vec at=(*ct)->atoms.begin(); at!=(*ct)->atoms.end(); ++at) {
            is.clear();
            is.str((*at)->name);
            is >> name;
            if (name == "N") {
                pN = *at;
                last_num = (*at)->res_number;
            } else if (name == "CA") {
                if ((*at)->res_number != last_num) {
                    pN = 0; pCA = 0; pC = 0;
                    last_num = -1;
                    continue;
                }
                pCA = *at;
            } else if (name == "C") {
                if ((*at)->res_number != last_num) {
                    pN = 0; pCA = 0; pC = 0;
                    last_num = -1;
                    continue;
                }
                pC = *at;
            } else if (name == "O") {
                if ((*at)->res_number != last_num) {
                    pN = 0; pCA = 0; pC = 0;
                    last_num = -1;
                    continue;
                }
                if ((!pN.zero()) && (!pCA.zero()) && (!pC.zero())) {
                    aa.push_back(new C_AS(pN,pCA,pC,*at));
                    pN = 0; pCA = 0; pC = 0;
                }
            }
        }
    }
    for (chains_vec ct=ref->chains.begin(); ct!=ref->chains.end(); ++ct) {
        stl_ptr<ATOM> pN(0);
        stl_ptr<ATOM> pCA(0);
        stl_ptr<ATOM> pC(0);
        int last_num = -1;
        for (atoms_vec at=(*ct)->atoms.begin(); at!=(*ct)->atoms.end(); ++at) {
            is.clear();
            is.str((*at)->name);
            is >> name;
            if (name == "N") {
                pN = *at;
                last_num = (*at)->res_number;
            } else if (name == "CA") {
                if ((*at)->res_number != last_num) {
                    pN = 0; pCA = 0; pC = 0;
                    last_num = -1;
                    continue;
                }
                pCA = *at;
            } else if (name == "C") {
                if ((*at)->res_number != last_num) {
                    pN = 0; pCA = 0; pC = 0;
                    last_num = -1;
                    continue;
                }
                pC = *at;
            } else if (name == "O") {
                if ((*at)->res_number != last_num) {
                    pN = 0; pCA = 0; pC = 0;
                    last_num = -1;
                    continue;
                }
                if ((!pN.zero()) && (!pCA.zero()) && (!pC.zero())) {
                    ref_aa.push_back(new C_AS(pN,pCA,pC,*at));
                    pN = 0; pCA = 0; pC = 0;
                }
            }
        }
    }

    float dihedral_thresh = 0.2;
    unsigned int theo = aa.size() * ref_aa.size();
    if (theo > 1000000) dihedral_thresh = 0.1;
    if (theo > 5000000) dihedral_thresh = 0.05;
    dihedral_thresh *= dihedral_thresh;

//    cerr << "Pairing " << aa.size() << " * " << ref_aa.size() << endl;
    
    //! 2.) Moegliche Paarungen erstellen:
    tr1::unordered_map<int,vector<C_PAIR*> > pair_map;
    vector<C_PAIR*> to_kill;
    int last_id = -1;
    for (vector<C_AS*>::iterator at=aa.begin(); at!=aa.end(); ++at) {
        pair_map[(*at)->ca->intern_id] = vector<C_PAIR*>();
        for (vector<C_AS*>::iterator ref_at=ref_aa.begin(); ref_at!=ref_aa.end(); ++ref_at) {
            float err = ((*at)->dihed - (*ref_at)->dihed);
            err *= err;
            if (err < dihedral_thresh) {
                C_PAIR* cp = new C_PAIR(*at,*ref_at);
                to_kill.push_back(cp);
                pair_map[(*at)->ca->intern_id].push_back(cp);

                if (last_id != -1) {
                    for (vector<C_PAIR*>::iterator lt=pair_map[last_id].begin(); lt!=pair_map[last_id].end(); ++lt) {
                        if ((*lt)->ref_aa->ca->intern_id == cp->ref_aa->ca->intern_id) continue;
                        if (abs((*lt)->ref_aa->ca->res_number - cp->ref_aa->ca->res_number) == 1) {
                            (*lt)->neighbours.push_back(cp);
                            cp->neighbours.push_back(*lt);
                        }
                    }
                }
            }
        }
        last_id = (*at)->ca->intern_id;
    }
    list<C_CLUSTER*> clusts;
    for (vector<C_PAIR*>::iterator it=to_kill.begin(); it!=to_kill.end(); ++it) clusts.push_back(new C_CLUSTER(*it));


//    cerr << "Clustering for " << clusts.size() << endl;

    //! 3.) Heuristisch clustern:
    bool changed = true;
    while (changed) {
        changed = false;
        //! 1.) Cluster erweitern:
        for (list<C_CLUSTER*>::iterator it=clusts.begin(); it!=clusts.end(); ++it) {
            vector<C_PAIR*> ptp;
            for (vector<C_PAIR*>::iterator pt=(*it)->new_pairs.begin(); pt!=(*it)->new_pairs.end(); ++pt) ptp.push_back(*pt);
            (*it)->new_pairs.clear();
            for (vector<C_PAIR*>::iterator pt=ptp.begin(); pt!=ptp.end(); ++pt) {
                for (vector<C_PAIR*>::iterator nt=(*pt)->neighbours.begin(); nt!=(*pt)->neighbours.end(); ++nt) {
                    if ((*it)->add_pair(*nt,max_displace)) {
                        changed = true;
                        (*it)->new_pairs.push_back(*nt);
                    }
                }
            }
        }
    }

//    cerr << "Merging" << endl;

    //! 4.) Cluster mergen, wenn moeglich:
    unsigned int biggest = 3;
    for (list<C_CLUSTER*>::iterator it=clusts.begin(); it!=clusts.end(); ++it) {
        if ((*it)->pairs.size() > biggest) biggest = (*it)->pairs.size();
    }

//    cerr << "Biggest = " << biggest << endl;

    unsigned int min_frag = biggest / 3;
    if (min_frag < 3) min_frag = 3;
    multimap<float,C_CLUSTER*> sorted_clq;
    for (list<C_CLUSTER*>::iterator it=clusts.begin(); it!=clusts.end(); ++it) {
        if ((*it)->pairs.size() < min_frag) continue;

        vector<vec3d<float> > A_comp_list;
        vector<vec3d<float> > B_comp_list;
        matrix<float> toptalign_rotm;
        vec3d<float> toptalign_trans[2];
        for (tr1::unordered_set<C_PAIR*>::iterator at=(*it)->pairs.begin(); at!=(*it)->pairs.end(); ++at) {
            A_comp_list.push_back((*at)->aa->ca->coord);
            B_comp_list.push_back((*at)->ref_aa->ca->coord);
        }
        get_align_matrix(B_comp_list,A_comp_list,toptalign_rotm,toptalign_trans[0],toptalign_trans[1]);
        for (vector<vec3d<float> >::iterator at=A_comp_list.begin(); at!=A_comp_list.end(); ++at) {
            *at -= toptalign_trans[1];
            *at *= toptalign_rotm;
            *at += toptalign_trans[0];
        }
        float score = 0.;
        for (unsigned int i=0; i<A_comp_list.size(); ++i) score += get_square_distance(A_comp_list[i],B_comp_list[i]);
        score /= A_comp_list.size();
        

        sorted_clq.insert(pair<float,C_CLUSTER*>(score,*it));
    }

//    cerr << "scored" << endl;

    tr1::unordered_set<C_CLUSTER*> to_del;
    int count1 = 0;
    for (multimap<float,C_CLUSTER*>::iterator it=sorted_clq.begin(); it!=sorted_clq.end(); ++it) {
        if (it->second->pairs.size() < min_frag) continue;
        if (to_del.find(it->second) != to_del.end()) continue;
        ++count1;
        multimap<float,C_CLUSTER*>::iterator jt=it; ++jt;
        int count2 = count1;
        for (; jt!=sorted_clq.end(); ++jt) {
            if (jt->second->pairs.size() < min_frag) continue;
            if (it->second->join(jt->second,max_displace)) {
                to_del.insert(jt->second);
            }
            ++count2;
            if (count2 > 1000) break;
        }
        if (count1 > 1000) break;
    }


    //! 5.) Besten cluster ermitteln:
    biggest = 3;
    for (list<C_CLUSTER*>::iterator it=clusts.begin(); it!=clusts.end(); ++it) {
        if ((*it)->pairs.size() > biggest) biggest = (*it)->pairs.size();
    }
    matrix<float> toptalign_rotm;
    vec3d<float> toptalign_trans[2];
    float best_rmsd = 99999999.;
    int matched = biggest;
    for (list<C_CLUSTER*>::iterator it=clusts.begin(); it!=clusts.end(); ++it) {
        if ((*it)->pairs.size() < biggest) continue;
        vector<vec3d<float> > A_comp_list;
        vector<vec3d<float> > B_comp_list;
        for (tr1::unordered_set<C_PAIR*>::iterator at=(*it)->pairs.begin(); at!=(*it)->pairs.end(); ++at) {
            A_comp_list.push_back((*at)->aa->ca->coord);
            B_comp_list.push_back((*at)->ref_aa->ca->coord);
        }
        get_align_matrix(B_comp_list,A_comp_list,toptalign_rotm,toptalign_trans[0],toptalign_trans[1]);
        for (vector<vec3d<float> >::iterator at=A_comp_list.begin(); at!=A_comp_list.end(); ++at) {
            *at -= toptalign_trans[1];
            *at *= toptalign_rotm;
            *at += toptalign_trans[0];
        }
        float test_rmsd = 0.;
        for (unsigned int i=0; i<A_comp_list.size(); ++i) test_rmsd += get_square_distance(A_comp_list[i],B_comp_list[i]);
        test_rmsd /= A_comp_list.size();
//        cerr << "test_rmsd = " << test_rmsd << endl;
        if (test_rmsd < best_rmsd) {
            best_rmsd = test_rmsd;
            optalign_rotm = toptalign_rotm;
            optalign_trans[0] = toptalign_trans[0];
            optalign_trans[1] = toptalign_trans[1];
        }
    }

    if (best_rmsd > 99999998.) best_rmsd = -1.;
    else best_rmsd = sqrt(best_rmsd);
    if (main_structure->verb > 0 && sim_only == false) {
        cout << " -> rmsd = " << best_rmsd << "  (" << matched
             << " C-alpha used in spatial alignment  /  similarity score = ";
        cout << ((10./matched) + (2.*best_rmsd/matched)) << ")";
        cout << "\n";
    }

    for (vector<C_AS*>::iterator ref_at=aa.begin(); ref_at!=aa.end(); ++ref_at) delete *ref_at;
    for (vector<C_AS*>::iterator ref_at=ref_aa.begin(); ref_at!=ref_aa.end(); ++ref_at) delete *ref_at;
    for (vector<C_PAIR*>::iterator it=to_kill.begin(); it!=to_kill.end(); ++it) delete *it;
    for (list<C_CLUSTER*>::iterator it=clusts.begin(); it!=clusts.end(); ++it) delete *it;

    if (sim_only) {
        //return ((10./matched) + (best_rmsd/matched)); // GUT!
        return ((10./matched) + (2.*best_rmsd/matched)); // GUT!
    } else return best_rmsd;
//    if (best_rmsd < -0.01) return false;
//    else return true;
}
//##############################################################################
//##############################################################################


//class C_NODE {
//public:
//    vector<stl_ptr<ATOM> > cat1;
//    vector<stl_ptr<ATOM> > cat2;
//    vec3d<float> geo_mean1;
//    vec3d<float> geo_mean2;
//    vec3d<float> v_sum1;
//    vec3d<float> v_sum2;
//    tr1::unordered_set<int> id_hashs1;
//    tr1::unordered_set<int> id_hashs2;
//    C_NODE(vector<stl_ptr<ATOM> >& ca1,vector<stl_ptr<ATOM> >& ca2) {
//        for (atoms_vec at=ca1.begin(); at!=ca1.end(); ++at) cat1.push_back(*at);
//        for (atoms_vec at=ca2.begin(); at!=ca2.end(); ++at) cat2.push_back(*at);
//        for (int i=0; i<3; ++i) {
//            v_sum1[i] = 0.;
//            v_sum2[i] = 0.;
//        }
//        for (atoms_vec at=cat1.begin(); at!=cat1.end(); ++at) {
//            id_hashs1.insert((*at)->intern_id);
//            v_sum1 += (*at)->coord;
//        }
//        for (atoms_vec at=cat2.begin(); at!=cat2.end(); ++at) {
//            id_hashs2.insert((*at)->intern_id);
//            v_sum2 += (*at)->coord;
//        }
//        geo_mean1 = v_sum1;
//        geo_mean2 = v_sum2;
//        geo_mean1 /= float(cat1.size());
//        geo_mean2 /= float(cat1.size());
//    }
//    void update() {
//        id_hashs1.clear();
//        id_hashs2.clear();
//        for (int i=0; i<3; ++i) {
//            v_sum1[i] = 0.;
//            v_sum2[i] = 0.;
//        }
//        for (atoms_vec at=cat1.begin(); at!=cat1.end(); ++at) {
//            id_hashs1.insert((*at)->intern_id);
//            v_sum1 += (*at)->coord;
//        }
//        for (atoms_vec at=cat2.begin(); at!=cat2.end(); ++at) {
//            id_hashs2.insert((*at)->intern_id);
//            v_sum2 += (*at)->coord;
//        }
//        geo_mean1 = v_sum1;
//        geo_mean2 = v_sum2;
//        geo_mean1 /= float(cat1.size());
//        geo_mean2 /= float(cat1.size());
//    }
//    bool connectable(C_NODE const& cn) {
//        float dst = fabs(get_distance(geo_mean1,cn.geo_mean1) - get_distance(geo_mean2,cn.geo_mean2));
//        if (dst > 0.4) return false;
//        for (atoms_vec at=cat1.begin(); at!=cat1.end(); ++at) {
//            if (cn.id_hashs1.find((*at)->intern_id) != cn.id_hashs1.end()) return false;
//        }
//        for (atoms_vec at=cat2.begin(); at!=cat2.end(); ++at) {
//            if (cn.id_hashs2.find((*at)->intern_id) != cn.id_hashs2.end()) return false;
//        }
//        return true;
//    }
//    void merge(C_NODE const& cn) {
//        for (const_atoms_vec at=cn.cat1.begin(); at!=cn.cat1.end(); ++at) {
//            cat1.push_back(*at);
//            v_sum1 += (*at)->coord;
//            id_hashs1.insert((*at)->intern_id);
//        }
//        for (const_atoms_vec at=cn.cat2.begin(); at!=cn.cat2.end(); ++at) {
//            cat2.push_back(*at);
//            v_sum2 += (*at)->coord;
//            id_hashs2.insert((*at)->intern_id);
//        }
//        geo_mean1 = v_sum1;
//        geo_mean2 = v_sum2;
//        geo_mean1 /= float(cat1.size());
//        geo_mean2 /= float(cat1.size());
//    }
//};
//
//
//bool PROTEIN::align2(PROTEIN *ref,unsigned int const& early_term) {
//    //! Diese Funktion ist fuer Faelle gedacht in denen das Sequenzalignment versagt
//    //! Das Alignment beruht hier auf einem C-alpha Graphmatching:
//
//    //! PARAMETER ***********************
//    const float gamma_thresh = 0.140;
//    const int res_num_div = 3;
//    const float dist_thresh = 0.4;
//    const float gdiv_thresh = 0.297;
//    //! *********************************
//
//
//    //! 1.) Koordinaten der C-alphas sammeln:
//    vector<stl_ptr<ATOM> > ca;
//    vector<stl_ptr<ATOM> > ref_ca;
//    tr1::unordered_map<int,stl_ptr<ATOM> > cc;
//    tr1::unordered_map<int,stl_ptr<ATOM> > ref_cc;
//    istringstream is;
//    string name;
//    stl_ptr<ATOM> last_at = 0;
//    for (chains_vec ct=chains.begin(); ct!=chains.end(); ++ct) {
//        for (atoms_vec at=(*ct)->atoms.begin(); at!=(*ct)->atoms.end(); ++at) {
//            is.clear();
//            is.str((*at)->name);
//            is >> name;
//            if (name == "CA") {
//                ca.push_back(*at);
//                last_at = *at;
//            } else if (name == "C") {
//                if (!last_at.zero() && (*at)->res_number == last_at->res_number) cc[last_at->intern_id] = *at;
//            }
//        }
//    }
//    last_at = 0;
//    for (chains_vec ct=ref->chains.begin(); ct!=ref->chains.end(); ++ct) {
//        for (atoms_vec at=(*ct)->atoms.begin(); at!=(*ct)->atoms.end(); ++at) {
//            is.clear();
//            is.str((*at)->name);
//            is >> name;
//            if (name == "CA") {
//                ref_ca.push_back(*at);
//                last_at = *at;
//            } else if (name == "C") {
//                if (!last_at.zero() && (*at)->res_number == last_at->res_number) ref_cc[last_at->intern_id] = *at;
//            }
//        }
//    }
//
//    //! 2.) Jetzt von jedem C-alpha aus folgende Vektoren bestimmen und als Ursprungsvectoren ablegen:
//    //!     v1 = Vektor zum uebernaechsten C-alpha
//    //!     v3 = Vektor von CA zu C
//    //!     Gibt es kein "*naechstes" C-alpha oder entspricht die Differenz in den res_numbers nicht dem
//    //!     gewuenschten C-alpha-Abstand (Kettenteile fehlen), so wird ein Nullvector abgelegt.
//    tr1::unordered_map<int,vec3d<float> > v1;
//    tr1::unordered_map<int,bool> has_v1;
//    tr1::unordered_map<int,vec3d<float> > v3;
//    tr1::unordered_map<int,bool> has_v3;
//    tr1::unordered_map<int,float> gamma;
//    tr1::unordered_map<int,vec3d<float> > ref_v1;
//    tr1::unordered_map<int,bool> ref_has_v1;
//    tr1::unordered_map<int,vec3d<float> > ref_v3;
//    tr1::unordered_map<int,bool> ref_has_v3;
//    tr1::unordered_map<int,float> ref_gamma;
//    float t_gam;
//    for (atoms_vec at=ca.begin(); at!=ca.end(); ++at) {
//        has_v1[(*at)->intern_id] = false;
//        has_v3[(*at)->intern_id] = false;
//        if (cc.find((*at)->intern_id) != cc.end()) {
//            vec3d<float> tvec(cc[(*at)->intern_id]->coord); tvec -= (*at)->coord; tvec.norm();
//            v3[(*at)->intern_id] = tvec;
//            has_v3[(*at)->intern_id] = true;
//        }
//    }
//    for (unsigned int i=0; i<(ca.size()-3); ++i) {
//        t_gam = 0.;
//        if (abs(ca[i+2]->res_number - ca[i]->res_number) == 2) {
//            vec3d<float> tvec(ca[i+2]->coord); tvec -= ca[i]->coord; tvec.norm();
//            v1[ca[i]->intern_id] = tvec;
//            has_v1[ca[i]->intern_id] = true;
//            if (has_v3[ca[i]->intern_id]) t_gam = angle_for_normed(v1[ca[i]->intern_id],v3[ca[i]->intern_id]);
//        }
//        gamma[ca[i]->intern_id] = t_gam;
//    }
//
//    for (atoms_vec at=ref_ca.begin(); at!=ref_ca.end(); ++at) {
//        ref_has_v1[(*at)->intern_id] = false;
//        ref_has_v3[(*at)->intern_id] = false;
//        if (ref_cc.find((*at)->intern_id) != ref_cc.end()) {
//            vec3d<float> tvec(ref_cc[(*at)->intern_id]->coord); tvec -= (*at)->coord; tvec.norm();
//            ref_v3[(*at)->intern_id] = tvec;
//            ref_has_v3[(*at)->intern_id] = true;
//        }
//    }
//    for (unsigned int i=0; i<(ref_ca.size()-3); ++i) {
//        t_gam = 0.;
//        if (abs(ref_ca[i+2]->res_number - ref_ca[i]->res_number) == 2) {
//            vec3d<float> tvec(ref_ca[i+2]->coord); tvec -= ref_ca[i]->coord; tvec.norm();
//            ref_v1[ref_ca[i]->intern_id] = tvec;
//            ref_has_v1[ref_ca[i]->intern_id] = true;
//            if (ref_has_v3[ref_ca[i]->intern_id]) t_gam = angle_for_normed(ref_v1[ref_ca[i]->intern_id],ref_v3[ref_ca[i]->intern_id]);
//        }
//        ref_gamma[ref_ca[i]->intern_id] = t_gam;
//    }
//
//
//    //! 3.) Teilalignments von ca auf ca_ref bestimmen:
//    int chunk_size = 20;
//    int ref_chunk_size = (chunk_size * 2) - 1;
//    vector<C_NODE*> curr_cliques;
//    vector<C_NODE*> next_cliques;
//    int next_ii = 0;
//    int ii = 0;
//    for (int chunk=0; ii<int(ca.size()); ++chunk) {
//        int next_ref_ii = 0;
//        int ref_ii = 0;
//        for (int ref_chunk=0; ref_ii<int(ref_ca.size()); ++ref_chunk) {
//            bkn_container f_list;
//            for (ii=next_ii; ii<int(ca.size()) && ii<(next_ii+chunk_size); ++ii) {
//                for (ref_ii=next_ref_ii; ref_ii<int(ref_ca.size()) && ref_ii<(next_ref_ii+ref_chunk_size); ++ref_ii) {
//                    if (has_v1[ca[ii]->intern_id] && has_v3[ca[ii]->intern_id]) {
//                        if (ref_has_v1[ref_ca[ref_ii]->intern_id] && ref_has_v3[ref_ca[ref_ii]->intern_id]) {
//                            if (fabs(gamma[ca[ii]->intern_id]-ref_gamma[ref_ca[ref_ii]->intern_id]) > gamma_thresh) continue;
//                        }
//                    }
//                    bkn_pointer new_node = new bkn_node(ca[ii],ref_ca[ref_ii]);
//                    f_list.push_back(new_node);
//                }
//            }
//
//            next_ref_ii += chunk_size;
//
//            for (bkn_vec it=f_list.begin(); it!=f_list.end(); ++it) {
//                for (bkn_vec jt=it+1; jt!=f_list.end(); ++jt) {
//                    if ((*it)->obj1.equal_addr((*jt)->obj1) || (*it)->obj2.equal_addr((*jt)->obj2)) continue;
//
//                    int rd = abs((*it)->obj1->res_number - (*jt)->obj1->res_number) - abs((*it)->obj2->res_number - (*jt)->obj2->res_number);
//                    if (abs(rd) > res_num_div) continue;
//
//                    float sd = fabs(get_distance((*it)->obj1->coord,(*jt)->obj1->coord) - get_distance((*it)->obj2->coord,(*jt)->obj2->coord));
//                    if (sd > dist_thresh) continue;
//
//                    if (has_v3[(*it)->obj1->intern_id] && has_v3[(*jt)->obj1->intern_id] &&
//                        ref_has_v3[(*it)->obj2->intern_id] && ref_has_v3[(*jt)->obj2->intern_id]) {
//                        float sg = fabs(angle_for_normed(v3[(*it)->obj1->intern_id],v3[(*jt)->obj1->intern_id])
//                                        -angle_for_normed(ref_v3[(*it)->obj2->intern_id],ref_v3[(*jt)->obj2->intern_id]));
//                        if (sg > gdiv_thresh) continue;
//                    }
//
//                    if (has_v1[(*it)->obj1->intern_id] && has_v1[(*jt)->obj1->intern_id] &&
//                        ref_has_v1[(*it)->obj2->intern_id] && ref_has_v1[(*jt)->obj2->intern_id]) {
//                        float sg = fabs(angle_for_normed(v1[(*it)->obj1->intern_id],v1[(*jt)->obj1->intern_id])
//                                        -angle_for_normed(ref_v1[(*it)->obj2->intern_id],ref_v1[(*jt)->obj2->intern_id]));
//                        if (sg > gdiv_thresh) continue;
//                    }
//
//                    (*it)->edges.insert(*jt);
//                    (*jt)->edges.insert(*it);
//                }
//            }
//
//            BK_SOLVER<stl_ptr<ATOM> > msolver(f_list);
//            msolver.solve(3,false,10,10);
//
//            for (clique_map it=msolver.cliques.begin(); it!=msolver.cliques.end(); ++it) {
//                if (it->first < msolver.max_clique) continue;
//                vector<stl_ptr<ATOM> > nv1;
//                vector<stl_ptr<ATOM> > nv2;
//                for (bkn_vec jt=it->second.begin(); jt!=it->second.end(); ++jt) {
//                    nv1.push_back((*jt)->obj1);
//                    nv2.push_back((*jt)->obj2);
//                }
//                C_NODE* curr_cn = new C_NODE(nv1,nv2);
//                curr_cliques.push_back(curr_cn);
//            }
//
//            for (bkn_vec it=f_list.begin(); it!=f_list.end(); ++it) delete *it;
//        }
//        next_ii += chunk_size;
//    }
//
//    //cerr << curr_cliques.size() << "  initial cliques" << endl;
//
//    unsigned int biggest = 0;
//    tr1::unordered_set<unsigned int> was_merged;
//    for (unsigned int i=0; i<curr_cliques.size(); ++i) {
//        if (was_merged.find(i) != was_merged.end()) continue;
//        for (unsigned int j=i+1; j<curr_cliques.size(); ++j) {
//            if (curr_cliques[i]->connectable(*(curr_cliques[j]))) {
//                bkn_container f_list;
//                for (unsigned int ii=0; ii<curr_cliques[i]->cat1.size(); ++ii) {
//                    bkn_pointer new_node = new bkn_node(curr_cliques[i]->cat1[ii],curr_cliques[i]->cat2[ii]);
//                    f_list.push_back(new_node);
//                }
//                for (unsigned int ii=0; ii<curr_cliques[j]->cat1.size(); ++ii) {
//                    bkn_pointer new_node = new bkn_node(curr_cliques[j]->cat1[ii],curr_cliques[j]->cat2[ii]);
//                    f_list.push_back(new_node);
//                }
//                for (bkn_vec it=f_list.begin(); it!=f_list.end(); ++it) {
//                    for (bkn_vec jt=it+1; jt!=f_list.end(); ++jt) {
//                        if ((*it)->obj1.equal_addr((*jt)->obj1) || (*it)->obj2.equal_addr((*jt)->obj2)) continue;
//
//                        int rd = abs((*it)->obj1->res_number - (*jt)->obj1->res_number) - abs((*it)->obj2->res_number - (*jt)->obj2->res_number);
//                        if (abs(rd) > res_num_div) continue;
//
//                        float sd = fabs(get_distance((*it)->obj1->coord,(*jt)->obj1->coord) - get_distance((*it)->obj2->coord,(*jt)->obj2->coord));
//                        if (sd > dist_thresh) continue;
//
//                        if (has_v3[(*it)->obj1->intern_id] && has_v3[(*jt)->obj1->intern_id] &&
//                            ref_has_v3[(*it)->obj2->intern_id] && ref_has_v3[(*jt)->obj2->intern_id]) {
//                            float sg = fabs(angle_for_normed(v3[(*it)->obj1->intern_id],v3[(*jt)->obj1->intern_id])
//                                            -angle_for_normed(ref_v3[(*it)->obj2->intern_id],ref_v3[(*jt)->obj2->intern_id]));
//                            if (sg > gdiv_thresh) continue;
//                        }
//
//                        if (has_v1[(*it)->obj1->intern_id] && has_v1[(*jt)->obj1->intern_id] &&
//                            ref_has_v1[(*it)->obj2->intern_id] && ref_has_v1[(*jt)->obj2->intern_id]) {
//                            float sg = fabs(angle_for_normed(v1[(*it)->obj1->intern_id],v1[(*jt)->obj1->intern_id])
//                                            -angle_for_normed(ref_v1[(*it)->obj2->intern_id],ref_v1[(*jt)->obj2->intern_id]));
//                            if (sg > gdiv_thresh) continue;
//                        }
//                        (*it)->edges.insert(*jt);
//                        (*jt)->edges.insert(*it);
//                    }
//                }
//                BK_SOLVER<stl_ptr<ATOM> > msolver(f_list);
//                msolver.solve(3,false,10,10);
//                for (clique_map it=msolver.cliques.begin(); it!=msolver.cliques.end(); ++it) {
//                    if (it->first < msolver.max_clique) continue;
//                    curr_cliques[i]->cat1.clear();
//                    curr_cliques[i]->cat2.clear();
//                    for (bkn_vec jt=it->second.begin(); jt!=it->second.end(); ++jt) {
//                        curr_cliques[i]->cat1.push_back((*jt)->obj1);
//                        curr_cliques[i]->cat2.push_back((*jt)->obj2);
//                    }
//                    curr_cliques[i]->update();
//                    was_merged.insert(j);
//                    break; // Nur die erstbeste groesste Clique nehmen
//                }
//                for (bkn_vec it=f_list.begin(); it!=f_list.end(); ++it) delete *it;
//
//                if (curr_cliques[i]->cat1.size() > biggest) {
//                    biggest = curr_cliques[i]->cat1.size();
//                }
//            }
//        }
//
//        //cerr << "i = " << i << "   biggest = " << biggest << endl;
//
//        if (biggest > early_term) break;
//    }
//
//
//    //! 6.) C_NODE (nur die groessten) mit dem besten RMSD ermitteln:
//    matrix<float> toptalign_rotm;
//    vec3d<float> toptalign_trans[2];
//    float best_rmsd = 99999999.;
//    int matched = biggest;
//    for (unsigned int i=0; i<curr_cliques.size(); ++i) {
//        if (curr_cliques[i]->cat1.size() < biggest) continue;
//        vector<vec3d<float> > A_comp_list;
//        vector<vec3d<float> > B_comp_list;
//        for (atoms_vec at=curr_cliques[i]->cat1.begin(); at!=curr_cliques[i]->cat1.end(); ++at) {
//            A_comp_list.push_back((*at)->coord);
//        }
//        for (atoms_vec at=curr_cliques[i]->cat2.begin(); at!=curr_cliques[i]->cat2.end(); ++at) {
//            B_comp_list.push_back((*at)->coord);
//        }
//        get_align_matrix(B_comp_list,A_comp_list,toptalign_rotm,toptalign_trans[0],toptalign_trans[1]);
//        for (vector<vec3d<float> >::iterator at=A_comp_list.begin(); at!=A_comp_list.end(); ++at) {
//            *at -= toptalign_trans[1];
//            *at *= toptalign_rotm;
//            *at += toptalign_trans[0];
//        }
//        float test_rmsd = 0.;
//        for (unsigned int i=0; i<A_comp_list.size(); ++i) test_rmsd += get_square_distance(A_comp_list[i],B_comp_list[i]);
//        test_rmsd /= A_comp_list.size();
//        if (test_rmsd < best_rmsd) {
//            best_rmsd = test_rmsd;
//            optalign_rotm = toptalign_rotm;
//            optalign_trans[0] = toptalign_trans[0];
//            optalign_trans[1] = toptalign_trans[1];
//        }
//    }
//
//    if (best_rmsd > 99999998.) best_rmsd = -1.;
//    else best_rmsd = sqrt(best_rmsd);
//    if (main_structure->verb > 0) {
//        cout << " -> rmsd = " << best_rmsd << "  (" << matched
//             << " C-alpha used in spatial alignment)";
//        if (biggest > early_term) cout << " (early break after >" << early_term << "C-alpha)" << "\n";
//        else cout << "\n";
//    }
//
//    for (unsigned int i=0; i<curr_cliques.size(); ++i) delete curr_cliques[i];
//
//    if (best_rmsd < -0.01) return false;
//    else return true;
//}


void PROTEIN::align2b(PROTEIN *ref,vector<stl_ptr<ATOM> >&atms_a,vector<stl_ptr<ATOM> >&atms_b) {
    //! align2 Version fuer cavsimX
    
    //! PARAMETER ***********************
    //! grosse Proteine moeglich => harte Thresholds verwenden
    const int chunk_size = 10;//20;
    const float gamma_thresh = 0.140; //8 grad //0.105; //ca. 6 Grad
    const int res_num_div = 3;
    const float dist_thresh = 0.4;
    const float gdiv_thresh = 0.297; //17 grad //0.175; //ca. 10 Grad
    //! *********************************

    const int ref_chunk_size = (chunk_size * 2) - 1;

    //! 1.) Koordinaten der C-alphas sammeln:
    vector<stl_ptr<ATOM> > ca;
    vector<stl_ptr<ATOM> > ref_ca;
    tr1::unordered_map<int,stl_ptr<ATOM> > cc;
    tr1::unordered_map<int,stl_ptr<ATOM> > ref_cc;
    istringstream is;
    string name;
    stl_ptr<ATOM> last_at = 0;
    for (chains_vec ct=chains.begin(); ct!=chains.end(); ++ct) {
        for (atoms_vec at=(*ct)->atoms.begin(); at!=(*ct)->atoms.end(); ++at) {
            is.clear();
            is.str((*at)->name);
            is >> name;
            if (name == "CA") {
                ca.push_back(*at);
                last_at = *at;
            } else if (name == "C") {
                if (!last_at.zero() && (*at)->res_number == last_at->res_number) cc[last_at->intern_id] = *at;
            }
        }
    }
    last_at = 0;
    for (chains_vec ct=ref->chains.begin(); ct!=ref->chains.end(); ++ct) {
        for (atoms_vec at=(*ct)->atoms.begin(); at!=(*ct)->atoms.end(); ++at) {
            is.clear();
            is.str((*at)->name);
            is >> name;
            if (name == "CA") {
                ref_ca.push_back(*at);
                last_at = *at;
            } else if (name == "C") {
                if (!last_at.zero() && (*at)->res_number == last_at->res_number) ref_cc[last_at->intern_id] = *at;
            }
        }
    }
    
    //! 2.) Jetzt von jedem C-alpha aus folgende Vektoren bestimmen und als Ursprungsvectoren ablegen:
    //!     v1 = Vektor zum uebernaechsten C-alpha
    //!     v3 = Vektor von CA zu C
    //!     Gibt es kein "*naechstes" C-alpha oder entspricht die Differenz in den res_numbers nicht dem
    //!     gewuenschten C-alpha-Abstand (Kettenteile fehlen), so wird ein Nullvector abgelegt.
    tr1::unordered_map<int,vec3d<float> > v1;
    tr1::unordered_map<int,bool> has_v1;
    tr1::unordered_map<int,vec3d<float> > v3;
    tr1::unordered_map<int,bool> has_v3;
    tr1::unordered_map<int,float> gamma;
    tr1::unordered_map<int,vec3d<float> > ref_v1;
    tr1::unordered_map<int,bool> ref_has_v1;
    tr1::unordered_map<int,vec3d<float> > ref_v3;
    tr1::unordered_map<int,bool> ref_has_v3;
    tr1::unordered_map<int,float> ref_gamma;
    float t_gam;
    for (atoms_vec at=ca.begin(); at!=ca.end(); ++at) {
        has_v1[(*at)->intern_id] = false;
        has_v3[(*at)->intern_id] = false;
        if (cc.find((*at)->intern_id) != cc.end()) {
            vec3d<float> tvec(cc[(*at)->intern_id]->coord); tvec -= (*at)->coord; tvec.norm();
            v3[(*at)->intern_id] = tvec;
            has_v3[(*at)->intern_id] = true;
        }
    }
    for (unsigned int i=0; i<(ca.size()-3); ++i) {
        t_gam = -5.;
        if (ca[i+2]->res_number - ca[i]->res_number == 2) {
            vec3d<float> tvec(ca[i+2]->coord); tvec -= ca[i]->coord; tvec.norm();
            v1[ca[i]->intern_id] = tvec;
            has_v1[ca[i]->intern_id] = true;
            if (has_v3[ca[i]->intern_id]) t_gam = angle_for_normed(v1[ca[i]->intern_id],v3[ca[i]->intern_id]);
        }
        gamma[ca[i]->intern_id] = t_gam;
    }
    
    for (atoms_vec at=ref_ca.begin(); at!=ref_ca.end(); ++at) {
        ref_has_v1[(*at)->intern_id] = false;
        ref_has_v3[(*at)->intern_id] = false;
        if (ref_cc.find((*at)->intern_id) != ref_cc.end()) {
            vec3d<float> tvec(ref_cc[(*at)->intern_id]->coord); tvec -= (*at)->coord; tvec.norm();
            ref_v3[(*at)->intern_id] = tvec;
            ref_has_v3[(*at)->intern_id] = true;
        }
    }
    for (unsigned int i=0; i<(ref_ca.size()-3); ++i) {
        t_gam = -5.;
        if (ref_ca[i+2]->res_number - ref_ca[i]->res_number == 2) {
            vec3d<float> tvec(ref_ca[i+2]->coord); tvec -= ref_ca[i]->coord; tvec.norm();
            ref_v1[ref_ca[i]->intern_id] = tvec;
            ref_has_v1[ref_ca[i]->intern_id] = true;
            if (ref_has_v3[ref_ca[i]->intern_id]) t_gam = angle_for_normed(ref_v1[ref_ca[i]->intern_id],ref_v3[ref_ca[i]->intern_id]);
        }
        ref_gamma[ref_ca[i]->intern_id] = t_gam;
    }

    //! 3.) Teilloesungen bestimmen:
    //!     Teilalignments von ca auf ca_ref bestimmen und die entsprechenden Knoten
    //!     in die Liste fuer das globale Alignment uebernehmen.
    bkn_container test_list;
    unsigned int max_local = 3;
    int next_ii = 0;
    int ii = 0;
    tr1::unordered_set<uint64_t> is_in_test_list;
    for (int chunk=0; ii<int(ca.size()); ++chunk) {
        int next_ref_ii = 0;
        int ref_ii = 0;
        for (int ref_chunk=0; ref_ii<int(ref_ca.size()); ++ref_chunk) {
            bkn_container t_list;
            for (ii=next_ii; ii<int(ca.size()) && ii<((chunk+1)*chunk_size); ++ii) {
                //for (atoms_vec bt=ref_ca.begin(); bt!=ref_ca.end(); ++bt) {
                for (ref_ii=next_ref_ii; ref_ii<int(ref_ca.size()) && ref_ii<(next_ref_ii+ref_chunk_size); ++ref_ii) {
                    if (has_v1[ca[ii]->intern_id] && has_v3[ca[ii]->intern_id]) {
                        if (ref_has_v1[ref_ca[ref_ii]->intern_id] && ref_has_v3[ref_ca[ref_ii]->intern_id]) {
                            if (fabs(gamma[ca[ii]->intern_id]-ref_gamma[ref_ca[ref_ii]->intern_id]) > gamma_thresh) continue;
                        }
                    }
                    bkn_pointer new_node = new bkn_node(ca[ii],ref_ca[ref_ii]);
                    t_list.push_back(new_node);
                }
            }

            next_ref_ii += chunk_size;

            for (bkn_vec it=t_list.begin(); it!=t_list.end(); ++it) {
                for (bkn_vec jt=it+1; jt!=t_list.end(); ++jt) {
                    if ((*it)->obj1.equal_addr((*jt)->obj1) || (*it)->obj2.equal_addr((*jt)->obj2)) continue;

                    //int rd = abs(((*it)->obj1->res_number - (*jt)->obj1->res_number) - ((*it)->obj2->res_number - (*jt)->obj2->res_number));
                    //if (rd > res_num_div) continue;
                    int rd = abs((*it)->obj1->res_number - (*jt)->obj1->res_number) - abs((*it)->obj2->res_number - (*jt)->obj2->res_number);
                    if (abs(rd) > res_num_div) continue;

                    float sd = fabs(get_distance((*it)->obj1->coord,(*jt)->obj1->coord) - get_distance((*it)->obj2->coord,(*jt)->obj2->coord));
                    if (sd > dist_thresh) continue;

                    if (has_v3[(*it)->obj1->intern_id] && has_v3[(*jt)->obj1->intern_id] &&
                        ref_has_v3[(*it)->obj2->intern_id] && ref_has_v3[(*jt)->obj2->intern_id]) {
                        float sg = fabs(angle_for_normed(v3[(*it)->obj1->intern_id],v3[(*jt)->obj1->intern_id])
                                        -angle_for_normed(ref_v3[(*it)->obj2->intern_id],ref_v3[(*jt)->obj2->intern_id]));
                        if (sg > gdiv_thresh) continue;
                    }

                    if (has_v1[(*it)->obj1->intern_id] && has_v1[(*jt)->obj1->intern_id] &&
                        ref_has_v1[(*it)->obj2->intern_id] && ref_has_v1[(*jt)->obj2->intern_id]) {
                        float sg = fabs(angle_for_normed(v1[(*it)->obj1->intern_id],v1[(*jt)->obj1->intern_id])
                                        -angle_for_normed(ref_v1[(*it)->obj2->intern_id],ref_v1[(*jt)->obj2->intern_id]));
                        if (sg > gdiv_thresh) continue;
                    }

                    (*it)->edges.insert(*jt);
                    (*jt)->edges.insert(*it);
                }
            }

            bkn_container f_list;
            for (bkn_vec it=t_list.begin(); it!=t_list.end(); ++it) {
                if ((*it)->edges.size() > 0) f_list.push_back(*it);
                else delete *it; //# 0802
            }

            BK_SOLVER<stl_ptr<ATOM> > msolver(f_list);
            //msolver.solve(3,false,10,10);
            msolver.solve(max_local,false,10,10);

            if (msolver.max_clique > max_local) max_local = msolver.max_clique;

            for (clique_map it=msolver.cliques.begin(); it!=msolver.cliques.end(); ++it) {
                if (it->first < msolver.max_clique) continue;
                for (bkn_vec jt=it->second.begin(); jt!=it->second.end(); ++jt) (*jt)->edges.clear();
            }

            for (bkn_vec it=f_list.begin(); it!=f_list.end(); ++it) {
                if ((*it)->edges.size() == 0) {
                    uint64_t node_hash = (*it)->obj1->intern_id * 1000000 + (*it)->obj2->intern_id;
                    if (is_in_test_list.find(node_hash) == is_in_test_list.end()) {
                        test_list.push_back(*it);
                        is_in_test_list.insert(node_hash);
                    } else delete *it;
                } else delete *it;
            }

            //! DEBUG:
            //cout << "next_ref_ii = " << next_ref_ii
            //     << "   ref_ii = " << ref_ii << "   n_nodes = " << test_list.size() << endl;
            //char dummy; cin >> dummy;
        }
        next_ii += chunk_size;
        //! DEBUG:
        //cout << "n_nodes = " << test_list.size() << endl;
        //cout << "next_ii = " << next_ii
        //     << "   ii = " << ii << "   n_nodes = " << test_list.size()
        //     << "   is_in_test_list.size() = " << is_in_test_list.size() << endl;
        //char dummy; cin >> dummy;
    }

    //! jetzt das globale Alignment (Kanten muessen neu bestimmt werden):
    //!DEBUG:
    //int n_nodes = test_list.size();
    //cout << "-> " << n_nodes << " nodes for graph matching" << endl;

    //! 4.) Kanten erzeugen basierend auf den Abstaenden und Winkeln:
    //int n_edges = 0;
    for (bkn_vec it=test_list.begin(); it!=test_list.end(); ++it) {
        for (bkn_vec jt=it+1; jt!=test_list.end(); ++jt) {
            if ((*it)->obj1.equal_addr((*jt)->obj1) || (*it)->obj2.equal_addr((*jt)->obj2)) continue;

            //! Auch die res_number-Differenz koennte hier als Kriterium genommen werden:
            int rd = abs(((*it)->obj1->res_number - (*jt)->obj1->res_number) - ((*it)->obj2->res_number - (*jt)->obj2->res_number));
            if (rd > res_num_div) continue;

            float sd = abs(get_distance((*it)->obj1->coord,(*jt)->obj1->coord) - get_distance((*it)->obj2->coord,(*jt)->obj2->coord));
            if (sd > dist_thresh) continue;

            if (has_v3[(*it)->obj1->intern_id] && has_v3[(*jt)->obj1->intern_id] &&
                ref_has_v3[(*it)->obj2->intern_id] && ref_has_v3[(*jt)->obj2->intern_id]) {
                float sg = abs(angle_for_normed(v3[(*it)->obj1->intern_id],v3[(*jt)->obj1->intern_id])
                              -angle_for_normed(ref_v3[(*it)->obj2->intern_id],ref_v3[(*jt)->obj2->intern_id]));
                if (sg > gdiv_thresh) continue;
            }

            if (has_v1[(*it)->obj1->intern_id] && has_v1[(*jt)->obj1->intern_id] &&
                ref_has_v1[(*it)->obj2->intern_id] && ref_has_v1[(*jt)->obj2->intern_id]) {
                float sg = abs(angle_for_normed(v1[(*it)->obj1->intern_id],v1[(*jt)->obj1->intern_id])
                              -angle_for_normed(ref_v1[(*it)->obj2->intern_id],ref_v1[(*jt)->obj2->intern_id]));
                if (sg > gdiv_thresh) continue;
            }

            //++n_edges;
            (*it)->edges.insert(*jt);
            (*jt)->edges.insert(*it);
        }
    }
    //!DEBUG:
    //cout << "-> " << n_edges << " edges for graph matching" << endl;

    //! 5.) Jetzt die Cliquen suchen:
    BK_SOLVER<stl_ptr<ATOM> > mysolver(test_list);
    mysolver.solve(max_local,false,100,10);

    //! 6.) Es koennen mehrere groesste Cliquen gefunden worden sein => nun muss die Clique ermittelt werden, die
    //!     den besten RMSD liefert:
    //!DEBUG:
//    cout << "-> searching for best clique" << endl;
    
    matrix<float> toptalign_rotm;
    vec3d<float> toptalign_trans[2];
    float best_rmsd = 99999999.;
    int matched = 0;
    vector<BK_NODE<stl_ptr<ATOM> > * >* bc_ptr = 0;
    for (clique_map it=mysolver.cliques.begin(); it!=mysolver.cliques.end(); ++it) {
        if (it->first < mysolver.max_clique) continue; //! nur die groessten Cliquen nehmen
        vector<vec3d<float> > A_comp_list;
        vector<vec3d<float> > B_comp_list;
        for (bkn_vec jt=it->second.begin(); jt!=it->second.end(); ++jt) {
            A_comp_list.push_back((*jt)->obj1->coord);
            B_comp_list.push_back((*jt)->obj2->coord);
        }
        get_align_matrix(B_comp_list,A_comp_list,toptalign_rotm,toptalign_trans[0],toptalign_trans[1]);
        for (vector<vec3d<float> >::iterator at=A_comp_list.begin(); at!=A_comp_list.end(); ++at) {
            *at -= toptalign_trans[1];
            *at *= toptalign_rotm;
            *at += toptalign_trans[0];
        }
        float test_rmsd = 0.;
        for (unsigned int i=0; i<A_comp_list.size(); ++i) test_rmsd += get_square_distance(A_comp_list[i],B_comp_list[i]);
        test_rmsd /= A_comp_list.size();
        if (test_rmsd < best_rmsd) {
            bc_ptr = &(it->second);
            matched = A_comp_list.size();
            best_rmsd = test_rmsd;
            optalign_rotm = toptalign_rotm;
            optalign_trans[0] = toptalign_trans[0];
            optalign_trans[1] = toptalign_trans[1];
        }
    }
    if (bc_ptr) {
        for (bkn_vec jt=bc_ptr->begin(); jt!=bc_ptr->end(); ++jt) {
            atms_a.push_back((*jt)->obj1);
            atms_b.push_back((*jt)->obj2);
        }
    }
    
    for (bkn_vec it=test_list.begin(); it!=test_list.end(); ++it) delete *it;
}


void PROTEIN::get_CA_poc_atoms(vector< stl_ptr<ATOM> > &cv) {
    istringstream is;
    string name;
    for (chains_vec ct=chains.begin(); ct!=chains.end(); ++ct) {
        for (atoms_vec at=(*ct)->atoms.begin(); at!=(*ct)->atoms.end(); ++at) {
            is.clear();
            is.str((*at)->name);
            is >> name;
            if (name == "CA") {
                bool breaker = false;
                for (ligands_vec lig=main_structure->ligands.begin(); lig!=main_structure->ligands.end(); ++lig) {
                    if ((*lig)->atoms.size() < 8) continue;
                    for (atoms_vec lat=(*lig)->atoms.begin(); lat!=(*lig)->atoms.end(); ++lat) {
                        if (get_square_distance((*at)->coord,(*lat)->coord) < 65.) {
                            cv.push_back(*at);
                            breaker = true; break;
                        }
                    }
                    if (breaker) break;
                }
            }
        }
    }
}

pair<int,float> PROTEIN::get_similarity(vector< stl_ptr<ATOM> > &ref) {
    //! 1.) Koordinaten der Taschen-C-alphas sammeln:
    vector< stl_ptr<ATOM> > com;
    get_CA_poc_atoms(com);
    
    //! 2.) Jetzt alle moeglichen Paare erzeugen:
    bkn_container test_list;
    for (atoms_vec at=com.begin(); at!=com.end(); ++at) {
        for (atoms_vec bt=ref.begin(); bt!=ref.end(); ++bt) {
            bkn_pointer new_node = new bkn_node(*at,*bt);
            test_list.push_back(new_node);
        }
    }
    
    if (ref.size() < 5 || com.size() < 5) {
        for (bkn_vec it=test_list.begin(); it!=test_list.end(); ++it) delete *it;
        return (pair<int,float>(0,-1.));
    }
    
    //! 3.) Kanten erzeugen basierend auf den Abstaenden:
    for (bkn_vec it=test_list.begin(); it!=test_list.end(); ++it) {
        for (bkn_vec jt=it+1; jt!=test_list.end(); ++jt) {
            if ((*it)->obj1.equal_addr((*jt)->obj1) || (*it)->obj2.equal_addr((*jt)->obj2)) continue;
            float sd = abs(get_distance((*it)->obj1->coord,(*jt)->obj1->coord) - get_distance((*it)->obj2->coord,(*jt)->obj2->coord));
            
            if (sd < 0.8) {
                (*it)->edges.insert(*jt);
                (*jt)->edges.insert(*it);
            }
        }
    }
    
    //! 4.) Jetzt die Cliquen suchen:
    BK_SOLVER<stl_ptr<ATOM> > mysolver(test_list);
    mysolver.solve();
    
    //! 5.) Es koennen mehrere groesste Cliquen gefunden worden sein => nun muss die Clique ermittelt werden, die
    //!     den besten RMSD liefert:
    float rmsd = 999999999999.;
    bkn_container * best_clique = 0;
    matrix<float> optalign_rotm;
    vec3d<float> optalign_trans[2];
    for (clique_map it=mysolver.cliques.begin(); it!=mysolver.cliques.end(); ++it) {
        if (it->first < mysolver.max_clique) continue;
        float test_rmsd = 0.;
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
        if (test_rmsd < rmsd) {
            rmsd = test_rmsd;
            best_clique = &(it->second);
        }
    }
    
    if (mysolver.max_clique < 3) {
        for (bkn_vec it=test_list.begin(); it!=test_list.end(); ++it) delete *it;
        return (pair<int,float>(0,-1.));
    }
    
    rmsd /= mysolver.max_clique;
    
    for (bkn_vec it=test_list.begin(); it!=test_list.end(); ++it) delete *it;
    
    return (pair<int,float>(mysolver.max_clique,sqrt(rmsd)));
}

void PROTEIN::search_spattern(PROTEIN *ref,multimap<float,vector<stl_ptr<ATOM> > > &mmap,
                              vector<vec3d<float> > &t1,vector<vec3d<float> > &t2,vector<matrix<float> > &rm) {
    //! gucken, ob das Sekundaerstrukturelement ref im Protein enthalten ist
    //! -> Graph-matching der C-alphas:
    //! 1.) Koordinaten der C-alphas sammeln:
    vector<stl_ptr<ATOM> > ca;
    vector<stl_ptr<ATOM> > ref_ca;
    istringstream is;
    string name;
    for (chains_vec ct=chains.begin(); ct!=chains.end(); ++ct) {
        for (atoms_vec at=(*ct)->atoms.begin(); at!=(*ct)->atoms.end(); ++at) {
            is.clear();
            is.str((*at)->name);
            is >> name;
            if (name == "CA") ca.push_back(*at);
        }
    }
    for (chains_vec ct=ref->chains.begin(); ct!=ref->chains.end(); ++ct) {
        for (atoms_vec at=(*ct)->atoms.begin(); at!=(*ct)->atoms.end(); ++at) {
            is.clear();
            is.str((*at)->name);
            is >> name;
            if (name == "CA") ref_ca.push_back(*at);
        }
    }
    
    //! 2.) Jetzt alle moeglichen Paare erzeugen:
    int n_nodes = 0;
    bkn_container test_list;
    for (atoms_vec at=ca.begin(); at!=ca.end(); ++at) {
        for (atoms_vec bt=ref_ca.begin(); bt!=ref_ca.end(); ++bt) {
            
            ++n_nodes;
            bkn_pointer new_node = new bkn_node(*at,*bt);
            test_list.push_back(new_node);
        }
    }
    //!DEBUG:
    cout << "-> " << n_nodes << " Knoten erzeugt" << endl;
    
    //! 3.) Kanten erzeugen basierend auf den Abstaenden:
    int n_edges = 0;
    for (bkn_vec it=test_list.begin(); it!=test_list.end(); ++it) {
        for (bkn_vec jt=it+1; jt!=test_list.end(); ++jt) {
            if ((*it)->obj1.equal_addr((*jt)->obj1) || (*it)->obj2.equal_addr((*jt)->obj2)) continue;
            float sd = abs(get_distance((*it)->obj1->coord,(*jt)->obj1->coord) - get_distance((*it)->obj2->coord,(*jt)->obj2->coord));
            
            if (sd < 0.8) {
        //    if (sd < 0.9) {
                
                ++n_edges;
                (*it)->edges.insert(*jt);
                (*jt)->edges.insert(*it);
            }
        }
    }
    //!DEBUG:
    cout << "-> " << n_edges << " Kanten erzeugt" << endl;
    
    //! 4.) Jetzt die Cliquen suchen:
    BK_SOLVER<stl_ptr<ATOM> > mysolver(test_list);
    mysolver.solve();
    
    //! 5.) Es koennen mehrere groesste Cliquen gefunden worden sein => nun muss die Clique ermittelt werden, die
    //!     den besten RMSD liefert:
    matrix<float> toptalign_rotm;
    vec3d<float> toptalign_trans[2];
    for (clique_map it=mysolver.cliques.begin(); it!=mysolver.cliques.end(); ++it) {
        if (it->first < mysolver.max_clique) continue; //! nur die groessten Cliquen nehmen
        float test_rmsd = 0.;
        vector<vec3d<float> > A_comp_list;
        vector<vec3d<float> > B_comp_list;
        for (bkn_vec jt=it->second.begin(); jt!=it->second.end(); ++jt) {
            A_comp_list.push_back((*jt)->obj1->coord);
            B_comp_list.push_back((*jt)->obj2->coord);
        }
        get_align_matrix(B_comp_list,A_comp_list,toptalign_rotm,toptalign_trans[0],toptalign_trans[1]);
        for (vector<vec3d<float> >::iterator at=A_comp_list.begin(); at!=A_comp_list.end(); ++at) {
            *at -= toptalign_trans[1];
            *at *= toptalign_rotm;
            *at += toptalign_trans[0];
        }
        t1.push_back(toptalign_trans[1]); rm.push_back(toptalign_rotm); t2.push_back(toptalign_trans[0]);
        for (unsigned int i=0; i<A_comp_list.size(); ++i) test_rmsd += get_square_distance(A_comp_list[i],B_comp_list[i]);
        test_rmsd = sqrt(test_rmsd/it->first);
        vector<stl_ptr<ATOM> > cv;
        int last_res = (*(it->second.begin()))->obj1->res_number;
        for (bkn_vec jt=it->second.begin(); jt!=it->second.end(); ++jt) {
            //! Wenn zwischen der aktuellen und der vorhergehenden res_number ein kleiner Abstand ist => auffuellen!
            int res_diff = (*jt)->obj1->res_number - last_res;
            if (res_diff > 1 && res_diff < 6) { //auffuellen:
                for (atoms_vec tl=ca.begin(); tl!=ca.end(); ++tl) {
                    if ((*tl)->res_number > last_res && (*tl)->res_number < (*jt)->obj1->res_number) {
                        if ((*tl)->chain_id == (*jt)->obj1->chain_id) cv.push_back(*tl);
                    }
                }
            }
            cv.push_back((*jt)->obj1);
            last_res = (*jt)->obj1->res_number;
        }
        mmap.insert(pair<float,vector<stl_ptr<ATOM> > >(test_rmsd,cv));
    }
    for (bkn_vec it=test_list.begin(); it!=test_list.end(); ++it) delete *it;
}

void PROTEIN::search_spattern2(PROTEIN *ref,multimap<float,vector<stl_ptr<ATOM> > > &mmap,
                               vector<vec3d<float> > &t1,vector<vec3d<float> > &t2,vector<matrix<float> > &rm) {
    //! In dieser Variante werden zusatzlich Vektoren zwischen den C-alphas als Information genutzt
    //! => weniger Knoten und Kanten im Graphen => deutlich schnellere Clique-Suche
    
    //! PARAMETER ***********************
    const float gamma_thresh = 0.175; //ca. 10�
    const float dist_thresh = 0.7;
    const float gdiv_thresh = 0.26; //ca. 15�
    const float dist_thresh2 = 0.8;
    const float gdiv_thresh2 = 0.785; //ca. 45�
    //! *********************************
    
    //! 1.) Koordinaten der C-alphas sammeln:
    vector<stl_ptr<ATOM> > ca;
    vector<stl_ptr<ATOM> > ref_ca;
    tr1::unordered_map<int,stl_ptr<ATOM> > cc;
    tr1::unordered_map<int,stl_ptr<ATOM> > ref_cc;
    istringstream is;
    string name;
    stl_ptr<ATOM> last_at = 0;
    for (chains_vec ct=chains.begin(); ct!=chains.end(); ++ct) {
        for (atoms_vec at=(*ct)->atoms.begin(); at!=(*ct)->atoms.end(); ++at) {
            is.clear();
            is.str((*at)->name);
            is >> name;
            if (name == "CA") {
                ca.push_back(*at);
                last_at = *at;
            } else if (name == "C") {
                if (!last_at.zero() && (*at)->res_number == last_at->res_number) cc[last_at->intern_id] = *at;
            }
        }
    }
    last_at = 0;
    for (chains_vec ct=ref->chains.begin(); ct!=ref->chains.end(); ++ct) {
        for (atoms_vec at=(*ct)->atoms.begin(); at!=(*ct)->atoms.end(); ++at) {
            is.clear();
            is.str((*at)->name);
            is >> name;
            if (name == "CA") {
                ref_ca.push_back(*at);
                last_at = *at;
            } else if (name == "C") {
                if (!last_at.zero() && (*at)->res_number == last_at->res_number) ref_cc[last_at->intern_id] = *at;
            }
        }
    }
    
    //! 2.) Jetzt von jedem C-alpha aus folgende Vektoren bestimmen und als Ursprungsvectoren ablegen:
    //!     v1 = Vektor zum uebernaechsten C-alpha
    //!     v3 = Vektor von CA zu C
    //!     Gibt es kein "*naechstes" C-alpha oder entspricht die Differenz in den res_numbers nicht dem
    //!     gewuenschten C-alpha-Abstand (Kettenteile fehlen), so wird ein Nullvector abgelegt.
    tr1::unordered_map<int,vec3d<float> > v1;
    tr1::unordered_map<int,bool> has_v1;
    tr1::unordered_map<int,vec3d<float> > v3;
    tr1::unordered_map<int,bool> has_v3;
    tr1::unordered_map<int,float> gamma;
    tr1::unordered_map<int,vec3d<float> > ref_v1;
    tr1::unordered_map<int,bool> ref_has_v1;
    tr1::unordered_map<int,vec3d<float> > ref_v3;
    tr1::unordered_map<int,bool> ref_has_v3;
    tr1::unordered_map<int,float> ref_gamma;
    float t_gam;
    for (atoms_vec at=ca.begin(); at!=ca.end(); ++at) {
        has_v1[(*at)->intern_id] = false;
        has_v3[(*at)->intern_id] = false;
        if (cc.find((*at)->intern_id) != cc.end()) {
            vec3d<float> tvec(cc[(*at)->intern_id]->coord); tvec -= (*at)->coord; tvec.norm();
            v3[(*at)->intern_id] = tvec;
            has_v3[(*at)->intern_id] = true;
        }
    }
    for (unsigned int i=0; i<(ca.size()-3); ++i) {
        t_gam = -5.;
        if (ca[i+2]->res_number - ca[i]->res_number == 2) {
            vec3d<float> tvec(ca[i+2]->coord); tvec -= ca[i]->coord; tvec.norm();
            v1[ca[i]->intern_id] = tvec;
            has_v1[ca[i]->intern_id] = true;
            if (has_v3[ca[i]->intern_id]) t_gam = angle_for_normed(v1[ca[i]->intern_id],v3[ca[i]->intern_id]);
        }
        gamma[ca[i]->intern_id] = t_gam;
    }
    
    for (atoms_vec at=ref_ca.begin(); at!=ref_ca.end(); ++at) {
        ref_has_v1[(*at)->intern_id] = false;
        ref_has_v3[(*at)->intern_id] = false;
        if (ref_cc.find((*at)->intern_id) != ref_cc.end()) {
            vec3d<float> tvec(ref_cc[(*at)->intern_id]->coord); tvec -= (*at)->coord; tvec.norm();
            ref_v3[(*at)->intern_id] = tvec;
            ref_has_v3[(*at)->intern_id] = true;
        }
    }
    for (unsigned int i=0; i<(ref_ca.size()-3); ++i) {
        t_gam = -5.;
        if (ref_ca[i+2]->res_number - ref_ca[i]->res_number == 2) {
            vec3d<float> tvec(ref_ca[i+2]->coord); tvec -= ref_ca[i]->coord; tvec.norm();
            ref_v1[ref_ca[i]->intern_id] = tvec;
            ref_has_v1[ref_ca[i]->intern_id] = true;
            if (ref_has_v3[ref_ca[i]->intern_id]) t_gam = angle_for_normed(ref_v1[ref_ca[i]->intern_id],ref_v3[ref_ca[i]->intern_id]);
        }
        ref_gamma[ref_ca[i]->intern_id] = t_gam;
    }
    
    //! 3.) Jetzt alle moeglichen Paare erzeugen:
    int n_nodes = 0;
    bkn_container test_list;
    for (atoms_vec at=ca.begin(); at!=ca.end(); ++at) {
        for (atoms_vec bt=ref_ca.begin(); bt!=ref_ca.end(); ++bt) {
            
            if (has_v1[(*at)->intern_id] && has_v3[(*at)->intern_id]) {
                if (ref_has_v1[(*bt)->intern_id] && ref_has_v3[(*bt)->intern_id]) {
                    if (abs(gamma[(*at)->intern_id]-ref_gamma[(*bt)->intern_id]) > gamma_thresh) continue;
                }
            }
            
            ++n_nodes;
            bkn_pointer new_node = new bkn_node(*at,*bt);
            test_list.push_back(new_node);
        }
    }
    //!DEBUG:
    cout << "-> " << n_nodes << " Knoten erzeugt" << endl;
    
    //! 4.) Kanten erzeugen basierend auf den Abstaenden und Winkeln:
    int n_edges = 0;
    for (bkn_vec it=test_list.begin(); it!=test_list.end(); ++it) {
        for (bkn_vec jt=it+1; jt!=test_list.end(); ++jt) {
            if ((*it)->obj1.equal_addr((*jt)->obj1) || (*it)->obj2.equal_addr((*jt)->obj2)) continue;
            
            float sd = abs(get_distance((*it)->obj1->coord,(*jt)->obj1->coord) - get_distance((*it)->obj2->coord,(*jt)->obj2->coord));
            if (sd > dist_thresh) continue;
            
            if (has_v3[(*it)->obj1->intern_id] && has_v3[(*jt)->obj1->intern_id] &&
                ref_has_v3[(*it)->obj2->intern_id] && ref_has_v3[(*jt)->obj2->intern_id]) {
                float sg = abs(angle_for_normed(v3[(*it)->obj1->intern_id],v3[(*jt)->obj1->intern_id])
                              -angle_for_normed(ref_v3[(*it)->obj2->intern_id],ref_v3[(*jt)->obj2->intern_id]));
                if (sg > gdiv_thresh) continue;
            }
            
            if (has_v1[(*it)->obj1->intern_id] && has_v1[(*jt)->obj1->intern_id] &&
                ref_has_v1[(*it)->obj2->intern_id] && ref_has_v1[(*jt)->obj2->intern_id]) {
                float sg = abs(angle_for_normed(v1[(*it)->obj1->intern_id],v1[(*jt)->obj1->intern_id])
                              -angle_for_normed(ref_v1[(*it)->obj2->intern_id],ref_v1[(*jt)->obj2->intern_id]));
                if (sg > gdiv_thresh) continue;
            }
            
            ++n_edges;
            (*it)->edges.insert(*jt);
            (*jt)->edges.insert(*it);
        }
    }
    //!DEBUG:
    cout << "-> " << n_edges << " Kanten erzeugt" << endl;
    
    //! 5.) Jetzt die Cliquen suchen:
    BK_SOLVER<stl_ptr<ATOM> > mysolver(test_list);
    mysolver.solve();
    
    //! 6.) Es koennen mehrere groesste Cliquen gefunden worden sein => nun muss die Clique ermittelt werden, die
    //!     den besten RMSD liefert:
    
    tr1::unordered_set<int> next_round;
    
    unsigned int min_clique_size = mysolver.max_clique - (mysolver.max_clique / 2);
    if (min_clique_size < 6) min_clique_size = mysolver.max_clique;
    
    for (clique_map it=mysolver.cliques.begin(); it!=mysolver.cliques.end(); ++it) {
        if (it->first < min_clique_size) continue;
        for (bkn_vec jt=it->second.begin(); jt!=it->second.end(); ++jt) {
            if (next_round.find((*jt)->obj1->intern_id) == next_round.end()) {
                next_round.insert((*jt)->obj1->intern_id);
            }
        }
    }
    
    bkn_container new_test_list;
    for (bkn_vec it=test_list.begin(); it!=test_list.end(); ++it) {
        if (next_round.find((*it)->obj1->intern_id) != next_round.end()) {
            bkn_pointer new_node = new bkn_node((*it)->obj1,(*it)->obj2);
            new_test_list.push_back(new_node);
        }
        delete *it;
    }
    
    cout << "-> " << new_test_list.size() << " Knoten im 2. Durchlauf" << endl;
    
    n_edges = 0;
    for (bkn_vec it=new_test_list.begin(); it!=new_test_list.end(); ++it) {
        for (bkn_vec jt=it+1; jt!=new_test_list.end(); ++jt) {
            if ((*it)->obj1.equal_addr((*jt)->obj1) || (*it)->obj2.equal_addr((*jt)->obj2)) continue;
            
            float sd = abs(get_distance((*it)->obj1->coord,(*jt)->obj1->coord) - get_distance((*it)->obj2->coord,(*jt)->obj2->coord));
            if (sd > dist_thresh2) continue;
            
            if (has_v3[(*it)->obj1->intern_id] && has_v3[(*jt)->obj1->intern_id] &&
                ref_has_v3[(*it)->obj2->intern_id] && ref_has_v3[(*jt)->obj2->intern_id]) {
                float sg = abs(angle_for_normed(v3[(*it)->obj1->intern_id],v3[(*jt)->obj1->intern_id])
                              -angle_for_normed(ref_v3[(*it)->obj2->intern_id],ref_v3[(*jt)->obj2->intern_id]));
                if (sg > gdiv_thresh2) continue;
            }
            
            ++n_edges;
            (*it)->edges.insert(*jt);
            (*jt)->edges.insert(*it);
        }
    }
    //!DEBUG:
    cout << "-> " << n_edges << " Kanten im 2. Durchlauf" << endl;
    
    BK_SOLVER<stl_ptr<ATOM> > new_mysolver(new_test_list);
    new_mysolver.solve();
    
    matrix<float> toptalign_rotm;
    vec3d<float> toptalign_trans[2];
    
    for (clique_map it=new_mysolver.cliques.begin(); it!=new_mysolver.cliques.end(); ++it) {
        if (it->first < new_mysolver.max_clique) continue;
        float test_rmsd = 0.;
        vector<vec3d<float> > A_comp_list;
        vector<vec3d<float> > B_comp_list;
        for (bkn_vec jt=it->second.begin(); jt!=it->second.end(); ++jt) {
            A_comp_list.push_back((*jt)->obj1->coord);
            B_comp_list.push_back((*jt)->obj2->coord);
        }
        get_align_matrix(B_comp_list,A_comp_list,toptalign_rotm,toptalign_trans[0],toptalign_trans[1]);
        for (vector<vec3d<float> >::iterator at=A_comp_list.begin(); at!=A_comp_list.end(); ++at) {
            *at -= toptalign_trans[1];
            *at *= toptalign_rotm;
            *at += toptalign_trans[0];
        }
        t1.push_back(toptalign_trans[1]); rm.push_back(toptalign_rotm); t2.push_back(toptalign_trans[0]);
        for (unsigned int i=0; i<A_comp_list.size(); ++i) test_rmsd += get_square_distance(A_comp_list[i],B_comp_list[i]);
        test_rmsd = sqrt(test_rmsd/it->first);
        int last_res = (*(it->second.begin()))->obj1->res_number;
        vector<stl_ptr<ATOM> > cv;
        for (bkn_vec jt=it->second.begin(); jt!=it->second.end(); ++jt) {
            int res_diff = (*jt)->obj1->res_number - last_res;
            if (res_diff > 1 && res_diff < 6) { //auffuellen:
                for (atoms_vec tl=ca.begin(); tl!=ca.end(); ++tl) {
                    if ((*tl)->res_number > last_res && (*tl)->res_number < (*jt)->obj1->res_number) {
                        if ((*tl)->chain_id == (*jt)->obj1->chain_id) cv.push_back(*tl);
                    }
                }
            }
            
            cv.push_back((*jt)->obj1);
            last_res = (*jt)->obj1->res_number;
        }
        mmap.insert(pair<float,vector<stl_ptr<ATOM> > >(test_rmsd,cv));
    }
    
    for (bkn_vec it=new_test_list.begin(); it!=new_test_list.end(); ++it) delete *it;
}

void PROTEIN::search_spattern3(PROTEIN *ref,tr1::unordered_set<int> &set1,tr1::unordered_set<int> &set2) {
    //! PARAMETER ***********************
    const float gamma_thresh = 0.175; //ca. 10�
    const float dist_thresh = 0.8;
    const float gdiv_thresh = 0.26; //ca. 15�
    //! *********************************
    
    //! 1.) Koordinaten der C-alphas sammeln:
    vector<stl_ptr<ATOM> > ca;
    vector<stl_ptr<ATOM> > ref_ca;
    tr1::unordered_map<int,stl_ptr<ATOM> > cc;
    tr1::unordered_map<int,stl_ptr<ATOM> > ref_cc;
    istringstream is;
    string name;
    stl_ptr<ATOM> last_at = 0;
    for (chains_vec ct=chains.begin(); ct!=chains.end(); ++ct) {
        for (atoms_vec at=(*ct)->atoms.begin(); at!=(*ct)->atoms.end(); ++at) {
            is.clear();
            is.str((*at)->name);
            is >> name;
            if (name == "CA") {
                ca.push_back(*at);
                last_at = *at;
            } else if (name == "C") {
                if (!last_at.zero() && (*at)->res_number == last_at->res_number) cc[last_at->intern_id] = *at;
            }
        }
    }
    last_at = 0;
    for (chains_vec ct=ref->chains.begin(); ct!=ref->chains.end(); ++ct) {
        for (atoms_vec at=(*ct)->atoms.begin(); at!=(*ct)->atoms.end(); ++at) {
            is.clear();
            is.str((*at)->name);
            is >> name;
            if (name == "CA") {
                ref_ca.push_back(*at);
                last_at = *at;
            } else if (name == "C") {
                if (!last_at.zero() && (*at)->res_number == last_at->res_number) ref_cc[last_at->intern_id] = *at;
            }
        }
    }
    
    //! 2.) Jetzt von jedem C-alpha aus folgende Vektoren bestimmen und als Ursprungsvectoren ablegen:
    //!     v1 = Vektor zum uebernaechsten C-alpha
    //!     v3 = Vektor von CA zu C
    //!     Gibt es kein "*naechstes" C-alpha oder entspricht die Differenz in den res_numbers nicht dem
    //!     gewuenschten C-alpha-Abstand (Kettenteile fehlen), so wird ein Nullvector abgelegt.
    tr1::unordered_map<int,vec3d<float> > v1;
    tr1::unordered_map<int,bool> has_v1;
    tr1::unordered_map<int,vec3d<float> > v3;
    tr1::unordered_map<int,bool> has_v3;
    tr1::unordered_map<int,float> gamma;
    tr1::unordered_map<int,vec3d<float> > ref_v1;
    tr1::unordered_map<int,bool> ref_has_v1;
    tr1::unordered_map<int,vec3d<float> > ref_v3;
    tr1::unordered_map<int,bool> ref_has_v3;
    tr1::unordered_map<int,float> ref_gamma;
    float t_gam;
    for (atoms_vec at=ca.begin(); at!=ca.end(); ++at) {
        has_v1[(*at)->intern_id] = false;
        has_v3[(*at)->intern_id] = false;
        if (cc.find((*at)->intern_id) != cc.end()) {
            vec3d<float> tvec(cc[(*at)->intern_id]->coord); tvec -= (*at)->coord; tvec.norm();
            v3[(*at)->intern_id] = tvec;
            has_v3[(*at)->intern_id] = true;
        }
    }
    for (unsigned int i=0; i<(ca.size()-3); ++i) {
        t_gam = -5.;
        if (ca[i+2]->res_number - ca[i]->res_number == 2) {
            vec3d<float> tvec(ca[i+2]->coord); tvec -= ca[i]->coord; tvec.norm();
            v1[ca[i]->intern_id] = tvec;
            has_v1[ca[i]->intern_id] = true;
            if (has_v3[ca[i]->intern_id]) t_gam = angle_for_normed(v1[ca[i]->intern_id],v3[ca[i]->intern_id]);
        }
        gamma[ca[i]->intern_id] = t_gam;
    }
    
    for (atoms_vec at=ref_ca.begin(); at!=ref_ca.end(); ++at) {
        ref_has_v1[(*at)->intern_id] = false;
        ref_has_v3[(*at)->intern_id] = false;
        if (ref_cc.find((*at)->intern_id) != ref_cc.end()) {
            vec3d<float> tvec(ref_cc[(*at)->intern_id]->coord); tvec -= (*at)->coord; tvec.norm();
            ref_v3[(*at)->intern_id] = tvec;
            ref_has_v3[(*at)->intern_id] = true;
        }
    }
    for (unsigned int i=0; i<(ref_ca.size()-3); ++i) {
        t_gam = -5.;
        if (ref_ca[i+2]->res_number - ref_ca[i]->res_number == 2) {
            vec3d<float> tvec(ref_ca[i+2]->coord); tvec -= ref_ca[i]->coord; tvec.norm();
            ref_v1[ref_ca[i]->intern_id] = tvec;
            ref_has_v1[ref_ca[i]->intern_id] = true;
            if (ref_has_v3[ref_ca[i]->intern_id]) t_gam = angle_for_normed(ref_v1[ref_ca[i]->intern_id],ref_v3[ref_ca[i]->intern_id]);
        }
        ref_gamma[ref_ca[i]->intern_id] = t_gam;
    }
    
    //! 3.) Jetzt alle moeglichen Paare erzeugen:
    int n_nodes = 0;
    bkn_container test_list;
    for (atoms_vec at=ca.begin(); at!=ca.end(); ++at) {
        for (atoms_vec bt=ref_ca.begin(); bt!=ref_ca.end(); ++bt) {
            
            if (has_v1[(*at)->intern_id] && has_v3[(*at)->intern_id]) {
                if (ref_has_v1[(*bt)->intern_id] && ref_has_v3[(*bt)->intern_id]) {
                    if (abs(gamma[(*at)->intern_id]-ref_gamma[(*bt)->intern_id]) > gamma_thresh) continue;
                }
            }
            
            ++n_nodes;
            bkn_pointer new_node = new bkn_node(*at,*bt);
            test_list.push_back(new_node);
        }
    }
    //!DEBUG:
    cout << "     -> " << n_nodes << " Knoten" << endl;
    
    //! 4.) Kanten erzeugen basierend auf den Abstaenden und Winkeln:
    int n_edges = 0;
    for (bkn_vec it=test_list.begin(); it!=test_list.end(); ++it) {
        for (bkn_vec jt=it+1; jt!=test_list.end(); ++jt) {
            if ((*it)->obj1.equal_addr((*jt)->obj1) || (*it)->obj2.equal_addr((*jt)->obj2)) continue;
            
            float sd = abs(get_distance((*it)->obj1->coord,(*jt)->obj1->coord) - get_distance((*it)->obj2->coord,(*jt)->obj2->coord));
            if (sd > dist_thresh) continue;
            
            if (has_v3[(*it)->obj1->intern_id] && has_v3[(*jt)->obj1->intern_id] &&
                ref_has_v3[(*it)->obj2->intern_id] && ref_has_v3[(*jt)->obj2->intern_id]) {
                float sg = abs(angle_for_normed(v3[(*it)->obj1->intern_id],v3[(*jt)->obj1->intern_id])
                              -angle_for_normed(ref_v3[(*it)->obj2->intern_id],ref_v3[(*jt)->obj2->intern_id]));
                if (sg > gdiv_thresh) continue;
            }
            
            if (has_v1[(*it)->obj1->intern_id] && has_v1[(*jt)->obj1->intern_id] &&
                ref_has_v1[(*it)->obj2->intern_id] && ref_has_v1[(*jt)->obj2->intern_id]) {
                float sg = abs(angle_for_normed(v1[(*it)->obj1->intern_id],v1[(*jt)->obj1->intern_id])
                              -angle_for_normed(ref_v1[(*it)->obj2->intern_id],ref_v1[(*jt)->obj2->intern_id]));
                if (sg > gdiv_thresh) continue;
            }
            
            ++n_edges;
            (*it)->edges.insert(*jt);
            (*jt)->edges.insert(*it);
        }
    }
    //!DEBUG:
    cout << "     -> " << n_edges << " Kanten" << endl;
    
    //! 5.) Jetzt die Cliquen suchen:
    BK_SOLVER<stl_ptr<ATOM> > mysolver(test_list);
    mysolver.solve();
    
    //! 6.) ...:
    unsigned int min_clique_size = mysolver.max_clique - (mysolver.max_clique / 3);
    if (min_clique_size < 7) min_clique_size = mysolver.max_clique;
    
    for (clique_map it=mysolver.cliques.begin(); it!=mysolver.cliques.end(); ++it) {
        if (it->first < min_clique_size) continue;
        int last_res1 = (*(it->second.begin()))->obj1->res_number;
        int last_res2 = (*(it->second.begin()))->obj2->res_number;
        for (bkn_vec jt=it->second.begin(); jt!=it->second.end(); ++jt) {
            
            int res_diff = (*jt)->obj1->res_number - last_res1;
            if (res_diff > 1 && res_diff < 6) { //auffuellen:
                for (atoms_vec tl=ca.begin(); tl!=ca.end(); ++tl) {
                    if ((*tl)->res_number > last_res1 && (*tl)->res_number < (*jt)->obj1->res_number) {
                        if ((*tl)->chain_id == (*jt)->obj1->chain_id) {
                            if (set1.find((*tl)->intern_id) == set1.end()) set1.insert((*tl)->intern_id);
                        }
                    }
                }
            }
            res_diff = (*jt)->obj2->res_number - last_res2;
            if (res_diff > 1 && res_diff < 6) { //auffuellen:
                for (atoms_vec tl=ref_ca.begin(); tl!=ref_ca.end(); ++tl) {
                    if ((*tl)->res_number > last_res2 && (*tl)->res_number < (*jt)->obj2->res_number) {
                        if ((*tl)->chain_id == (*jt)->obj2->chain_id) {
                            if (set2.find((*tl)->intern_id) == set2.end()) set2.insert((*tl)->intern_id);
                        }
                    }
                }
            }
            
            if (set1.find((*jt)->obj1->intern_id) == set1.end()) set1.insert((*jt)->obj1->intern_id);
            if (set2.find((*jt)->obj2->intern_id) == set2.end()) set2.insert((*jt)->obj2->intern_id);
            last_res1 = (*jt)->obj1->res_number;
            last_res2 = (*jt)->obj2->res_number;
        }
    }
    
    for (bkn_vec it=test_list.begin(); it!=test_list.end(); ++it) delete *it;
}

