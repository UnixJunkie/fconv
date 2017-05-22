
//============================================================================
// structure_GN.cpp -*- C++ -*-; high level representation of molecular data
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
// This library implements structure objects as high level abstraction of
// molecular data.
//============================================================================



#include"structure_GN.h"


#if defined (_SPEEDTEST)
#include"timer_GN.hpp"
TIMER_GN<double> t_tot;
TIMER_GN<double> t_del;
TIMER_GN<double> t_alpha;
TIMER_GN<double> t_diff;
TIMER_GN<double> t_clust;
TIMER_GN<double> t_refine;
#endif

//==============================================================================================
//Definitionen fuer HOLE:
//==============================================================================================

HOLE::HOLE() {}

HOLE::~HOLE() {}


//==============================================================================================
//Definitionen fuer CAVDUM:
//==============================================================================================

CAVDUM::CAVDUM() {}

CAVDUM::~CAVDUM() {}

bool CAVDUM::operator<(CAVDUM &rechts) {
    if (total_volume < rechts.total_volume) return true;
    return false;
}

void CAVDUM::get_volume() {
    total_volume = 0.;
    for (holes_map hit=holes.begin(); hit!=holes.end(); ++hit) {
        for (holes_vec hiv=hit->second.begin(); hiv!=hit->second.end(); ++hiv) {
            total_volume += (*hiv)->cells.size();
        }
    }
}

void CAVDUM::pymol_vis(int j,float spacing) {
    ofstream f_out;
    ostringstream os;
    os << "poc_" << j << "_vis.py";
    string filename = os.str();
    f_out.open(filename.c_str());

    f_out << "# This visualization file was automatically created by 'structure_GN.cpp'" << "\n";
    f_out << "# Start pymol and use 'run this_files_name.py' to start the visualization" << "\n" << "\n";
    
    float rad = 0.4;
    
    bool first = true;
    
    f_out << "cavpoints" << j << " = [";
    for (holes_map hit=holes.begin(); hit!=holes.end(); ++hit) {
        for (holes_vec hiv=hit->second.begin(); hiv!=hit->second.end(); ++hiv) {
            for (cells_vec cc=(*hiv)->cells.begin(); cc!=(*hiv)->cells.end(); ++cc) {
                if (!first) f_out << ",";
                else first = false;
                f_out << "7.0," << (*cc)->middle[0] << "," << (*cc)->middle[1] << "," << (*cc)->middle[2] << "," << rad;
            }
        }
    }
    f_out << "]" << endl;
    
    f_out << "cmd.load_cgo(cavpoints" << j << ", 'cavpoints" << j << "', 1)" << endl;
    f_out << "cmd.color('red', 'cavpoints" << j << "')" << "\n";
    
    f_out.close();
    
    cout << "pymol debug visualization file written to " << os.str() << endl;
}


//==============================================================================================
//Definitionen fuer STRUCTURE:
//==============================================================================================

STRUCTURE::STRUCTURE(int verbosity) : name("X"), has_peptides(false), protein(NULL), verb(verbosity), model_number(0) {}

STRUCTURE::~STRUCTURE() {    
    for (ligands_vec it=ligands.begin(); it!=ligands.end();++it) {
        it->kill();
    }
    for (ligands_vec it=splitted_ligands.begin(); it!=splitted_ligands.end();++it) {
        it->kill();
    }
    for (ligands_vec it=mol2_protein.begin(); it!=mol2_protein.end();++it) {
        it->kill();
    }    
    for (het_entrys_vec it=het_entrys.begin(); it!=het_entrys.end();++it) {
        it->kill();
    }    
    for (metals_vec it=metals.begin(); it!=metals.end();++it) {
        it->kill();
    }    
    for (cavities_vec it=cavities.begin(); it!=cavities.end();++it) {
        it->kill();
    }
    for (waters_vec it=waters.begin(); it!=waters.end();++it) {
        it->kill();
    }
    for (conects_vec it=conects.begin(); it!=conects.end();++it) {
        it->kill();
    }
    for (alt_struct_vec it=alt_struct.begin(); it!=alt_struct.end();++it) {
        it->kill();
    }
    if (protein) delete protein;
    if (!water_ligand.zero()) water_ligand.zero_kill();
    if (!metal_ligand.zero()) metal_ligand.zero_kill();
}

void STRUCTURE::clear() {
    if (ligands.size() > 0) {
        for (ligands_vec it=ligands.begin(); it!=ligands.end();++it) {
            it->kill();
        }
        ligands.clear();
    }
    
    if (splitted_ligands.size() > 0) {
        for (ligands_vec it=splitted_ligands.begin(); it!=splitted_ligands.end();++it) {
            it->kill();
        }
        splitted_ligands.clear();
    }
    
    if (mol2_protein.size() > 0) {
        for (ligands_vec it=mol2_protein.begin(); it!=mol2_protein.end();++it) {
            it->kill();
        }
        mol2_protein.clear();
    }
    
    if (het_entrys.size() > 0) {
        for (het_entrys_vec it=het_entrys.begin(); it!=het_entrys.end();++it) {
            it->kill();
        }
        het_entrys.clear();
    }
    
    if (metals.size() > 0) {
        for (metals_vec it=metals.begin(); it!=metals.end();++it) {
            it->kill();
        }
        metals.clear();
    }
    
    if (cavities.size() > 0) {
        for (cavities_vec it=cavities.begin(); it!=cavities.end();++it) {
            it->kill();
        }
        cavities.clear();
    }
    
    if (waters.size() > 0) {
        for (waters_vec it=waters.begin(); it!=waters.end();++it) {
            it->kill();
        }
        waters.clear();
    }
    
    if (conects.size() > 0) {
        for (conects_vec it=conects.begin(); it!=conects.end();++it) {
            it->kill();
        }
        conects.clear();
    }
    
    if (alt_struct.size() > 0) {
        for (alt_struct_vec it=alt_struct.begin(); it!=alt_struct.end();++it) {
            it->kill();
        }
        alt_struct.clear();
    }
    
    if (protein) delete protein;
    
    if (!water_ligand.zero()) water_ligand.zero_kill();
    if (!metal_ligand.zero()) metal_ligand.zero_kill();

    name = "X";
    protein = NULL;
    model_number = 0;
}

void STRUCTURE::clear_ligands() {
    if (ligands.size() > 0) {
        for (ligands_vec it=ligands.begin(); it!=ligands.end();++it) {
            it->kill();
        }
        ligands.clear();
    }
    if (metals.size() > 0) {
        for (metals_vec it=metals.begin(); it!=metals.end();++it) {
            it->kill();
        }
        metals.clear();
    }
    if (waters.size() > 0) {
        for (waters_vec it=waters.begin(); it!=waters.end();++it) {
            it->kill();
        }
        waters.clear();
    }
}

bool STRUCTURE::pdb2mol2(int mode,char const* def_file,bool fill_X,bool bond_building) {
    if (!protein) {
        cerr << c_message<cWARNING>("STRUCTURE::pdb2mol2 --> found no protein in input structure") << endl;
        return false;
    }
    protein->ter_correct(); //jetzt koennen ja wieder einfach die id's genommen werden
//    if (protein->n_aacids == 0) protein->get_aacids();
    for (chains_vec cv=protein->chains.begin(); cv!=protein->chains.end(); ++cv) if ((*cv)->n_aacids == 0) (*cv)->get_aacids();
    protein->get_atom_types(mode,def_file,fill_X);
    if (bond_building) protein->build_bonds();
    return true;
}

void STRUCTURE::full_pdb2mol2(int mode,bool get_ele,const char *def_file,bool get_bonds,int verb,
                           bool kill_ext,bool fill_X,const char *alt_def_file,int max_ring_members,
                           bool no_free_rot_planar,LIGAND* at_ref_mol,LIGAND* at_ref_mol2,bool prot_acids,
                           bool prot_guanidin,bool prot_amidin,bool prot_amin,bool prot_phosphate,bool prot_sulfate) {
    //zunaechst protein als Ligand in mol2_protein anlegen:
    stl_ptr<LIGAND> lig(new LIGAND());
    lig->name = "protein";
    lig->type = "PROTEIN";
    stl_ptr<STRUCTURE> this_struct = this;
    lig->main_structure = this_struct;
    int cur_id = 1;
    bool kill = true;
    
    if (protein != NULL) for (chains_vec it=protein->chains.begin(); it!=protein->chains.end();++it) {
        for (atoms_vec jt=(*it)->atoms.begin(); jt!=(*it)->atoms.end();++jt) {
            stl_ptr<ATOM> atm(new ATOM(**jt));
            atm->intern_id = cur_id; //interne id's neu vergeben um probleme mit ausgelassenen
                                     //atomen zu vermeiden
            atm->res_name = atm->full_res_name;
            cur_id++;
            lig->atoms.push_back(atm);
            kill = false;
        }
    }
    
    if (!kill) mol2_protein.push_back(lig);
    else lig.kill();
    
    //jetzt die Liganden anhaengen:
    for (ligands_vec lit=ligands.begin(); lit!=ligands.end(); ++lit) {
        stl_ptr<LIGAND> nlig(new LIGAND(**lit));
        nlig->main_structure = this_struct;
        mol2_protein.push_back(nlig);
    }
    
    //!Metalle auch mit nehmen: 20.08.2008
    for (metals_vec met=metals.begin(); met!=metals.end(); ++met) {
        stl_ptr<LIGAND> nlig = new LIGAND();
        stl_ptr<ATOM> natm(new ATOM(*((*met)->atom)));
        nlig->atoms.push_back(natm);
        nlig->main_structure = this_struct;
        mol2_protein.push_back(nlig);
    }
    //!=====================================
    
    //jetzt die alternativen Locations aufteilen:
    get_alt_loc_ligs(true);
    
    //jetzt die Atomtypzuweisung:
    for (ligands_vec lit=mol2_protein.begin(); lit!=mol2_protein.end(); ++lit) {
        (*lit)->get_atom_typing(mode,get_ele,def_file,get_bonds,verb,kill_ext,fill_X,alt_def_file,max_ring_members,no_free_rot_planar,at_ref_mol,
                                at_ref_mol2,prot_acids,prot_guanidin,prot_amidin,prot_amin,prot_phosphate,prot_sulfate);
    }
}

void STRUCTURE::get_lig_cavity(float dist, LIGAND* lig,bool full_res,bool no_hyd) {
    if (!protein) {
        cerr << c_message<cWARNING>("STRUCTURE::get_lig_cavity --> found no protein atoms in input structure") << endl;
        exit(1);
    }
    bool breaker;
    stl_ptr<CAVITY> cavity = new CAVITY();
    if (no_hyd) {
        lig->get_elements(false);
    }
    for (chains_vec cv=protein->chains.begin(); cv!=protein->chains.end(); ++cv) if ((*cv)->n_aacids == 0) (*cv)->get_aacids();
    dist *= dist;
    for (chains_vec jt=protein->chains.begin(); jt!=protein->chains.end();++jt) {
        for (aacids_map ut=(*jt)->aacids.begin(); ut!=(*jt)->aacids.end();++ut) {
            for (atoms_vec at=ut->second->atoms.begin(); at!=ut->second->atoms.end(); ++at) {
                (*at)->bond_ind = 0;
                breaker = false;
                for (atoms_vec it=lig->atoms.begin(); it!=lig->atoms.end(); ++it) {
                    if (no_hyd && (*at)->element == "H") continue;
                    if (get_square_distance((*it)->coord,(*at)->coord) <= dist) {
                        if (full_res) {
                            cavity->aacids.push_back(ut->second);
                            breaker = true;
                            break;
                        } else {
                            (*at)->bond_ind = 1;
                        }
                    }
                }
                if (breaker) break;
            }
        }
    }
    cavities.push_back(cavity);
}

void STRUCTURE::get_buriedness() {
    //! 1.) Elementtypen fuer alle Proteinatome bestimmen
    for (chains_vec ct=protein->chains.begin(); ct != protein->chains.end(); ++ct) {
        for (atoms_vec at=(*ct)->atoms.begin(); at!=(*ct)->atoms.end(); ++at) {
            //!neues ATOM_EXT Objekt:
            if ((*at)->ext.zero()) (*at)->ext = new ATOM_EXT();
        }
        (*ct)->get_elements(false);
    }
    
    //! 2.) Jetzt das Protein in ein Grid legen:
    float spacing = 1.; // erster Parameter zum spielen
    int wsteps = int((5. / spacing) + 0.5) + 5;
    stl_ptr<STRUCTURE_GRID> grid = new STRUCTURE_GRID(protein,spacing,wsteps);
    
    //! 3.) Leere Zellen bestimmen:
    grid->check_empty();
    
    //! 4.) Surface-Zellen bestimmen:
    grid->get_next_to_surf();
    
    //! 5.) Konvexe Surface-Zellen bestimmen, die mindestens in 8er Clustern auftreten:
    grid->get_real_surf(18);
    typedef GRID<3,CELL>::iterator g_iterator;
    for (g_iterator it=grid->grid->begin(); it!=grid->grid->end(); ++it) {
        if (it->next_to_surf) {
            it->friends = it->burial;
        }
    }
    grid->get_buriedness();
    
    
    //! 6.) Jetzt die alte Variante mit dem Schalenweisen Clustern, nur mit real_surf-cells statt empty:
    map<int,vector<stl_ptr<HOLE> > > holes; //! schale:holes
    vector<stl_ptr<CAVDUM> > cavdums;
    
    float con_dist = (2. * spacing * spacing) + 0.2;
    
    //Ueber alle Wuerfelschalen:
    for (int i=grid->grid->max_shell; i >=0; i--) {
        //Ueber alle Zellen der Schale:
        for (GRID<3,CELL>::shell_iterator it=grid->grid->sbegin(i); it!=grid->grid->send(i); ++it) {
            if (it->real_surf) {
                int to_con = -1;
                int hole_count = 0;
                for (holes_vec hvi=holes[i].begin(); hvi!=holes[i].end(); ++hvi) {
                    //pruefen, ob die Zelle zu einem bereits vorhandenen HOLE dieser Schale gehoert:
                    for (cells_vec ct=(*hvi)->cells.begin(); ct!=(*hvi)->cells.end(); ++ct) {
                        float helper = get_square_distance(it->middle,(*ct)->middle);
                        if ((helper < con_dist) && (helper > 0.2)) {
                            if (to_con != -1) { //!aktuelles HOLE mit to_con HOLE verbinden
                                //die Zelle *it gehoert also bereits zum to_con HOLE
                                for (cells_vec cbt=(*hvi)->cells.begin(); cbt!=(*hvi)->cells.end(); ++cbt) {
                                    holes[i][to_con]->cells.push_back(*cbt);
                                }
                                //! jetzt das aktuelle HOLE loeschen und sowohl hvi, 
                                //! als auch hole_count zuruecksetzen:
                                hvi->kill();
                                holes[i].erase(hvi); hvi--;
                                hole_count--;
                                break;
                            } else {
                                stl_ptr<CELL> ccc = &(*it);
                                (*hvi)->cells.push_back(ccc);
                                to_con = hole_count;
                                break;
                            }
                        }
                    }
                    hole_count++;
                }
                if (to_con == -1) { //Zelle gehoert noch zu keinem HOLE dieser Schale
                    stl_ptr<HOLE> hol(new HOLE());
                    stl_ptr<CELL> ccc = &(*it);
                    hol->cells.push_back(ccc);
                    holes[i].push_back(hol);
                }
            }
        }
    }
    
    //HOLEs zu CAVDUMs vereinigen
    for (int i=grid->grid->max_shell; i >=0; --i) {
        //Ueber alle HOLEs der schale:
        for (holes_vec hvi=holes[i].begin(); hvi!=holes[i].end(); ++hvi) {
            //Ueber alle CAVDUM Objekte:
            bool new_dum = true;
            int to_con = -1;
            int cav_count = 0;
            for (cavdums_vec cvi=cavdums.begin(); cvi!=cavdums.end(); ++cvi) {
                if ((*cvi)->holes.find(i+1) != (*cvi)->holes.end()) { //das CAVDUM hat HOLEs in der Schale eins tiefer
                    bool bigbreak = false;
                    //Ueber alle HOLEs unter der aktuellen Schale:
                    for (holes_vec dvi=(*cvi)->holes[i+1].begin(); dvi!=(*cvi)->holes[i+1].end(); ++dvi) {
                    for (cells_vec bt=(*hvi)->cells.begin(); bt!=(*hvi)->cells.end(); ++bt) {
                    for (cells_vec ct=(*dvi)->cells.begin(); ct!=(*dvi)->cells.end(); ++ct) {
                        float helper = get_square_distance((*bt)->middle,(*ct)->middle);
                        if ((helper < con_dist) && (helper > 0.2)) { //!von vert_con_dist auf con_dist umgestellt!!!
                            if (to_con != -1) {
                                //aktuelles CAVDUM mit to_con CAVDUM vereinen:
                                for (holes_map hx=(*cvi)->holes.begin(); hx!=(*cvi)->holes.end(); ++hx) {
                                    for (holes_vec hvx=(*cvi)->holes[hx->first].begin();
                                                   hvx!=(*cvi)->holes[hx->first].end(); ++hvx) {
                                        cavdums[to_con]->holes[hx->first].push_back(*hvx);
                                    }
                                }
                                cvi->kill();
                                cavdums.erase(cvi);
                                cvi--;
                                bigbreak = true;
                                break;
                            } else {
                                //das HOLE hvi an das CAVDUM cvi anhaengen:
                                (*cvi)->holes[i].push_back(*hvi);
                                new_dum = false;
                                bigbreak = true;
                                to_con = cav_count;
                                break;
                            }
                        }
                        if (bigbreak) break;
                    }
                    if (bigbreak) break;
                    }
                    if (bigbreak) break;
                    }
                    if (bigbreak) continue;
                }
                cav_count++;
            }
            if (new_dum) { //das HOLE hvi gehoert zu keinem bereits vorhandenen CAVDUM
                stl_ptr<CAVDUM> cav(new CAVDUM());
                cav->holes[i].push_back(*hvi);
                cavdums.push_back(cav);
            }
        }
    }
    
    //nach total_volume sortieren:
    //sort(cavdums.begin(),cavdums.end());
    //reverse(cavdums.begin(),cavdums.end());
    
    //! Hier jetzt Vergrabenheiten fuer die Atome zuweisen
    //! Jedes cavdum Objekt hat eine map holes, die als Schluessel integer hat:
    //! Je hoeher dieser integer, desto weiter aussen liegt die entsprechende Schale.
    //! Der Wert ist jeweils ein vector mit zeigern auf hole Objekte dieser Schale.
    //! Jedes hole Objekt hat einen vector mit den entsprechenden cell-Objekten.
    //!
    //! Die Vergrabenheit sollte Beeinflusst werden durch:
    //!  - den Wert des cell Objekts
    //!  - die Schale
    for (chains_vec ct=protein->chains.begin(); ct != protein->chains.end(); ++ct) {
        for (atoms_vec at=(*ct)->atoms.begin(); at!=(*ct)->atoms.end(); ++at) {
            (*at)->ext->sp = 999999.;
        }
    }
    
    
    for (cavdums_vec cvi=cavdums.begin(); cvi!=cavdums.end(); ++cvi) {
        int min_shell = 0;
        int max_shell = 999999;
        for (holes_map hit=(*cvi)->holes.begin(); hit!=(*cvi)->holes.end(); ++hit) {
            if (hit->first > max_shell) max_shell = hit->first;
            if (hit->first < min_shell) min_shell = hit->first;
        }
        //float scaler = (max_shell - min_shell) / 20.;
        for (holes_map hit=(*cvi)->holes.begin(); hit!=(*cvi)->holes.end(); ++hit) {
            float multiplier = 1. * (hit->first - min_shell + 1.);
            //float multiplier = 1.;
            for (holes_vec hiv=hit->second.begin(); hiv!=hit->second.end(); ++hiv) {
                for (cells_vec ci=(*hiv)->cells.begin(); ci!=(*hiv)->cells.end(); ++ci) {
                    if ((*ci)->next_to_surf == false) continue;
                    float curr_bury = (*ci)->burial * /*(*ci)->friends * */multiplier;
                    for (chains_vec ct=protein->chains.begin(); ct != protein->chains.end(); ++ct) {
                        for (atoms_vec at=(*ct)->atoms.begin(); at!=(*ct)->atoms.end(); ++at) {
                            float curr_sdst = get_square_distance((*ci)->middle,(*at)->coord);
                            if (curr_sdst < 16.) {
                            //    (*at)->ext->buriedness += curr_bury;
                                if ((*at)->ext->sp > curr_sdst) {
                                    (*at)->ext->buriedness = curr_bury;
                                    (*at)->ext->sp = curr_sdst;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    
    //Aufraeumen:
    grid.kill();
    
    for (holes_map hit=holes.begin(); hit!=holes.end(); ++hit) {
        for (holes_vec hiv=hit->second.begin(); hiv!=hit->second.end(); ++hiv) {
            hiv->kill();
        }
    }
    for (cavdums_vec cvi=cavdums.begin(); cvi!=cavdums.end(); ++cvi) {
        cvi->kill();
    }
    
}

void STRUCTURE::get_cavity_auto(float radius,int min_buried,bool debug_mode,unsigned int vis_number) {
    //! 1.) Elementtypen fuer alle Proteinatome bestimmen
    for (chains_vec ct=protein->chains.begin(); ct != protein->chains.end(); ++ct) {
        for (atoms_vec at=(*ct)->atoms.begin(); at!=(*ct)->atoms.end(); ++at) {
            //!neues ATOM_EXT Objekt:
            if ((*at)->ext.zero()) (*at)->ext = new ATOM_EXT();
        }
        (*ct)->get_elements(false);
    }
    
    //! 2.) Jetzt das Protein in ein Grid legen:
    float spacing = 1.0; // erster Parameter zum spielen
    stl_ptr<STRUCTURE_GRID> grid = new STRUCTURE_GRID(protein,spacing);
    
    //! 3.) Leere Zellen bestimmen:
    grid->check_empty();
    
    //! 4.) Surface-Zellen bestimmen:
    grid->get_next_to_surf();
    
    //! 5.) Konvexe Surface-Zellen bestimmen, die mindestens in 8er Clustern auftreten:
    grid->get_real_surf(min_buried,8);
    
    
    //! 6.) Jetzt die alte Variante mit dem Schalenweisen Clustern, nur mit real_surf-cells statt empty:
    map<int,vector<stl_ptr<HOLE> > > holes; //! schale:holes
    vector<stl_ptr<CAVDUM> > cavdums;
    
    float con_dist = (2. * spacing * spacing) + 0.2;
    
    //Ueber alle Wuerfelschalen:
    for (int i=grid->grid->max_shell; i >=0; i--) {
        //Ueber alle Zellen der Schale:
        for (GRID<3,CELL>::shell_iterator it=grid->grid->sbegin(i); it!=grid->grid->send(i); ++it) {
            if (it->real_surf) {
                int to_con = -1;
                int hole_count = 0;
                for (holes_vec hvi=holes[i].begin(); hvi!=holes[i].end(); ++hvi) {
                    //pruefen, ob die Zelle zu einem bereits vorhandenen HOLE dieser Schale gehoert:
                    for (cells_vec ct=(*hvi)->cells.begin(); ct!=(*hvi)->cells.end(); ++ct) {
                        float helper = get_square_distance(it->middle,(*ct)->middle);
                        if ((helper < con_dist) && (helper > 0.2)) {
                            if (to_con != -1) { //!aktuelles HOLE mit to_con HOLE verbinden
                                //die Zelle *it gehoert also bereits zum to_con HOLE
                                for (cells_vec cbt=(*hvi)->cells.begin(); cbt!=(*hvi)->cells.end(); ++cbt) {
                                    holes[i][to_con]->cells.push_back(*cbt);
                                }
                                //! jetzt das aktuelle HOLE loeschen und sowohl hvi, 
                                //! als auch hole_count zuruecksetzen:
                                hvi->kill();
                                holes[i].erase(hvi); hvi--;
                                hole_count--;
                                break;
                            } else {
                                stl_ptr<CELL> ccc = &(*it);
                                (*hvi)->cells.push_back(ccc);
                                to_con = hole_count;
                                break;
                            }
                        }
                    }
                    hole_count++;
                }
                if (to_con == -1) { //Zelle gehoert noch zu keinem HOLE dieser Schale
                    stl_ptr<HOLE> hol(new HOLE());
                    stl_ptr<CELL> ccc = &(*it);
                    hol->cells.push_back(ccc);
                    holes[i].push_back(hol);
                }
            }
        }
    }
    
    //HOLEs zu CAVDUMs vereinigen
    for (int i=grid->grid->max_shell; i >=0; i--) {
        //Ueber alle HOLEs der schale:
        for (holes_vec hvi=holes[i].begin(); hvi!=holes[i].end(); ++hvi) {
            //Ueber alle CAVDUM Objekte:
            bool new_dum = true;
            int to_con = -1;
            int cav_count = 0;
            for (cavdums_vec cvi=cavdums.begin(); cvi!=cavdums.end(); ++cvi) {
                if ((*cvi)->holes.find(i+1) != (*cvi)->holes.end()) { //das CAVDUM hat HOLEs in der Schale eins tiefer
                    bool bigbreak = false;
                    //Ueber alle HOLEs unter der aktuellen Schale:
                    for (holes_vec dvi=(*cvi)->holes[i+1].begin(); dvi!=(*cvi)->holes[i+1].end(); ++dvi) {
                    for (cells_vec bt=(*hvi)->cells.begin(); bt!=(*hvi)->cells.end(); ++bt) {
                    for (cells_vec ct=(*dvi)->cells.begin(); ct!=(*dvi)->cells.end(); ++ct) {
                        float helper = get_square_distance((*bt)->middle,(*ct)->middle);
                        if ((helper < con_dist) && (helper > 0.2)) { //!von vert_con_dist auf con_dist umgestellt!!!
                            if (to_con != -1) {
                                //aktuelles CAVDUM mit to_con CAVDUM vereinen:
                                for (holes_map hx=(*cvi)->holes.begin(); hx!=(*cvi)->holes.end(); ++hx) {
                                    for (holes_vec hvx=(*cvi)->holes[hx->first].begin();
                                                   hvx!=(*cvi)->holes[hx->first].end(); ++hvx) {
                                        cavdums[to_con]->holes[hx->first].push_back(*hvx);
                                    }
                                }
                                cvi->kill();
                                cavdums.erase(cvi);
                                cvi--;
                                bigbreak = true;
                                break;
                            } else {
                                //das HOLE hvi an das CAVDUM cvi anhaengen:
                                (*cvi)->holes[i].push_back(*hvi);
                                new_dum = false;
                                bigbreak = true;
                                to_con = cav_count;
                                break;
                            }
                        }
                        if (bigbreak) break;
                    }
                    if (bigbreak) break;
                    }
                    if (bigbreak) break;
                    }
                    if (bigbreak) continue;
                }
                cav_count++;
            }
            if (new_dum) { //das HOLE hvi gehoert zu keinem bereits vorhandenen CAVDUM
                stl_ptr<CAVDUM> cav(new CAVDUM());
                cav->holes[i].push_back(*hvi);
                cavdums.push_back(cav);
            }
        }
    }
    
    for (cavdums_vec cvi=cavdums.begin(); cvi!=cavdums.end(); ++cvi) {
        (*cvi)->get_volume();
        (*cvi)->total_volume *= (spacing * spacing * spacing);
        if ((*cvi)->total_volume < 100.) {
            cvi->kill();
            cavdums.erase(cvi);
            cvi--;
        }
    }
    
    //nach total_volume sortieren:
    sort(cavdums.begin(),cavdums.end());
    reverse(cavdums.begin(),cavdums.end());
    
    //anhand der CAVDUMs CAVITY Obj. erzeugen
    for (cavdums_vec cvi=cavdums.begin(); cvi!=cavdums.end(); ++cvi) {
        get_cavdum_cavity(radius,*cvi);
    }
    
    if (debug_mode) {
        for (unsigned int j=0; j<vis_number; ++j) {
            if (cavdums.size() > j) {
                cavdums[j]->pymol_vis(j+1,spacing);
            }
        }
    }
    
    //Aufraeumen:
    grid.kill();
    for (holes_map hit=holes.begin(); hit!=holes.end(); ++hit) {
        for (holes_vec hiv=hit->second.begin(); hiv!=hit->second.end(); ++hiv) {
            hiv->kill();
        }
    }
    for (cavdums_vec cvi=cavdums.begin(); cvi!=cavdums.end(); ++cvi) {
        cvi->kill();
    }
}


void STRUCTURE::get_cavdum_cavity(float dist,stl_ptr<CAVDUM> &cavdum) {
    vec3d<float> vv;
    bool breaker;
    stl_ptr<CAVITY> cavity = new CAVITY();
//    if (protein->n_aacids == 0) protein->get_aacids();
    for (chains_vec cv=protein->chains.begin(); cv!=protein->chains.end(); ++cv) if ((*cv)->n_aacids == 0) (*cv)->get_aacids();
    
    for (chains_vec jt=protein->chains.begin(); jt!=protein->chains.end();++jt) {
        for (aacids_map ut=(*jt)->aacids.begin(); ut!=(*jt)->aacids.end();++ut) {
            for (atoms_vec at=ut->second->atoms.begin(); at!=ut->second->atoms.end(); ++at) {
                breaker = false;        
                for (holes_map hit=cavdum->holes.begin(); hit!=cavdum->holes.end(); ++hit) {
                    for (holes_vec hiv=hit->second.begin(); hiv!=hit->second.end(); ++hiv) {
                    for (cells_vec ct=(*hiv)->cells.begin(); ct!=(*hiv)->cells.end(); ++ct) {
                        vv = (*ct)->middle;
                        vv -= (*at)->coord;
                        if (vv.value() <= dist) {
                            cavity->aacids.push_back(ut->second);
                            breaker = true;
                            break;
                        }
                    }
                    if (breaker) break;
                    }
                    if (breaker) break;
                }
                if (breaker) break;
            }
        }
    }
    cavity->total_volume = cavdum->total_volume;
    cavities.push_back(cavity);
}


void STRUCTURE::get_alpha_cavity(vector<stl_ptr<ATOM> >& pv,map<int,vector<stl_ptr<ATOM> > >& p_clust,map<int,double> &output_vol,
                         vector<stl_ptr<ATOM> >& lig_atoms,double alpha1,double alpha2,bool refine,double min_vol,
                         bool visualize,double def_rad,double surr_radius,bool get_cavity_obj,bool profile) {
    //! Diese Funktion ermittelt Cavitaeten als Differenzvolumen zwischen 2 Alphashapes.
    //! ACHTUNG: Zumindest die Elementtypen muessen fuer alle Atome in pv bereits gesetzt sein!!!
    #if defined (_SPEEDTEST)
    t_tot.continue_t();
    #endif
    
    if (alpha2 == 8.0 && lig_atoms.size() > 0) alpha2 = 15.0;
    
    //! 1.) DT_POINTs erzeugen:
    atom_properties::initialize();
    vector<DT_POINT> dtpoints;
    int id = 1;
    map<int,stl_ptr<ATOM> > points2atoms;
    
//    cerr << "ACHTUNG: Gewichte zu Debugzwecken auf NULL gesetzt!!!" << endl;
    for (atoms_vec at=pv.begin(); at!=pv.end(); ++at) {
        if ((*at)->element == "H") continue;
        double w = atom_properties::square_vdW_map[(*at)->element];
        double r = atom_properties::vdW_map[(*at)->element];
        DT_POINT tmp(id,(*at)->coord[0],(*at)->coord[1],(*at)->coord[2],w,r,(ATOM*)&(**at));
        points2atoms[id] = *at;
        ++id;
        dtpoints.push_back(tmp);
    }
    
    
    //! 2.) Triangulieren:
    if (verb > 0) cout << "-> calculating triangulation for " << id-1 << " weighted points" << endl;
    DT_SOLVER mysolver(dtpoints);
    mysolver.solve();
    
    if (visualize) cerr << "    => " << mysolver.delaunay_tetrahedrons.size() << " tetrahedrons in final triangulation" << endl;
    
    vec3d<double> shift(0.,0.,0.);
    vector<stl_ptr<L_FACET> > alpha_facets1;
    vector<stl_ptr<L_FACET> > alpha_facets2;
    if (visualize) {
        //!===== Fuer das Rausschreiben von wrl-Files ====================================
        float vals = 0.;
        for (vector<DT_POINT>::iterator it=dtpoints.begin(); it!=dtpoints.end(); ++it) {
            shift += it->coord;
            vals += 1.;
        }
        if (vals > 0.) shift /= vals;
        mysolver.get_alpha_shape(alpha_facets1,alpha1);
        mysolver.facets2wrl("facets1.wrl",alpha_facets1,shift);
        alpha_facets1.clear();
        mysolver.get_alpha_shape(alpha_facets2,alpha2);
        mysolver.facets2wrl("facets2.wrl",alpha_facets2,shift);
        alpha_facets2.clear();
        //!==============================================================================
        
        
        cerr << "    => calculating " << alpha1 << "A alpha shape" << endl;
        vector<stl_ptr<L_FACET> > alpha_facets1;
        mysolver.get_alpha_shape(alpha_facets1,alpha1);
        cerr << "    => writing " << alpha_facets1.size() << " " << alpha1 << "A facets as facets1.off" << endl;
        mysolver.facets2off("facets1.off",alpha_facets1);
        cerr << "    => calculating " << alpha2 << "A alpha shape" << endl;
    }
    
    mysolver.get_alpha_shape(alpha_facets2,alpha2);
    if (visualize) {
        cerr << "    => writing "  << alpha_facets2.size() << " " << alpha2 << "A facets as facets2.off" << endl;
        mysolver.facets2off("facets2.off",alpha_facets2);
    }
    
    //! 3.) Ein Set mit IDs aller Punkte anlegen, die Teil vom alpha2-Shape sind:
    tr1::unordered_set<int> fset;
    for (lf_vec lt=alpha_facets2.begin(); lt!=alpha_facets2.end(); ++lt) {
        fset.insert((*lt)->p0->id);
        fset.insert((*lt)->p1->id);
        fset.insert((*lt)->p2->id);
    }
    
    
    //! 4.) Differenz berechnen und Clustern:
    if (verb > 0) cout << "-> clustering difference volume" << endl;
    vector<stl_ptr<TETRAHEDRON> > pt;
    tr1::unordered_map<int64_t,int> id2idx;
    int cindex = 0;
    double srad2 = alpha2 * alpha2;
    double srad1 = alpha1 * alpha1;
    aforbidden.clear();
    for (th_vec it=mysolver.delaunay_tetrahedrons.begin(); it!=mysolver.delaunay_tetrahedrons.end(); ++it) {
        id2idx[(*it)->id] = cindex; ++cindex;
        if ((*it)->sphere_square_radius > srad2) continue; // gehoert nicht zum alpha2-Komplex
        if ((*it)->sphere_square_radius <= srad1) continue; // gehoert zum alpha1-Komplex
        
        //! ACHTUNG: Das gewichtete Volumen der Tetraeder ist nur ein approximierter Wert
        //!          Die exakte Berechnung erfolgt spaeter fuer eine definierte Menge von
        //!          Tetraedern ueber mysolver.calc_exact_volume
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
    
    if (visualize) {
        cerr << "    => " << pt.size() << " tetrahedrons in difference volume" << endl;
        string pname = "tets";
        string fname = pname + ".py";
        mysolver.visualize_tetrahedrons(pt,fname,pname);
    }
    
    avisited.clear();
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
    
    if (visualize) cerr << "    => " << clusters.size() << " clusters generated" << endl;
    
    //! 5.) Jetzt die aus aforbidden dazunehmen:
    if (refine) {
        if (verb > 0) cout << "-> refine clusters" << endl;
        for (list<ALPHACLUSTER*>::iterator it=clusters.begin(); it!=clusters.end(); ++it) {
            tr1::unordered_set<int64_t> add;
            for (tr1::unordered_set<int64_t>::iterator jt=(*it)->elements.begin(); jt!=(*it)->elements.end(); ++jt) {
                if (aforbidden.find(mysolver.delaunay_tetrahedrons[id2idx[*jt]]->p0n->id) != aforbidden.end()) {
                    add.insert(mysolver.delaunay_tetrahedrons[id2idx[*jt]]->p0n->id);
                }
                if (aforbidden.find(mysolver.delaunay_tetrahedrons[id2idx[*jt]]->p1n->id) != aforbidden.end()) {
                    add.insert(mysolver.delaunay_tetrahedrons[id2idx[*jt]]->p1n->id);
                }
                if (aforbidden.find(mysolver.delaunay_tetrahedrons[id2idx[*jt]]->p2n->id) != aforbidden.end()) {
                    add.insert(mysolver.delaunay_tetrahedrons[id2idx[*jt]]->p2n->id);
                }
                if (aforbidden.find(mysolver.delaunay_tetrahedrons[id2idx[*jt]]->p3n->id) != aforbidden.end()) {
                    add.insert(mysolver.delaunay_tetrahedrons[id2idx[*jt]]->p3n->id);
                }
            }
            for (tr1::unordered_set<int64_t>::iterator jt=add.begin(); jt!=add.end(); ++jt) {
                (*it)->elements.insert(*jt);
            }
        }
    }
    
    //! 6.) Cluster nach Volumen filtern:
    map<int,vector<stl_ptr<TETRAHEDRON> > > final_tets;
    map<int,vector<stl_ptr<TETRAHEDRON> > > output_tets;
    map<int,double> final_vol;
    int n_clust = 0;
    for (list<ALPHACLUSTER*>::iterator it=clusters.begin(); it!=clusters.end(); ++it) {
        double vol = 0.;
        for (tr1::unordered_set<int64_t>::iterator jt=(*it)->elements.begin(); jt!=(*it)->elements.end(); ++jt) {
            vol += mysolver.delaunay_tetrahedrons[id2idx[*jt]]->volume;
        }
        if (vol < min_vol) continue;
        final_tets[n_clust] = vector<stl_ptr<TETRAHEDRON> >();
        for (tr1::unordered_set<int64_t>::iterator jt = (*it)->elements.begin(); jt!=(*it)->elements.end(); ++jt) {
            final_tets[n_clust].push_back(mysolver.delaunay_tetrahedrons[id2idx[*jt]]);
        }
        ///    final_vol[n_clust] = vol;
        final_vol[n_clust] = mysolver.get_exact_void_volume(final_tets[n_clust]);
        
        if (visualize) {
            int n_sig = 0;
            for (th_vec tt=final_tets[n_clust].begin(); tt!=final_tets[n_clust].end(); ++tt) {
                if ((*tt)->sphere_square_radius > 1.) ++n_sig;
            }
            cerr << "       -> cluster " << n_clust << " with " << final_tets[n_clust].size() 
                 << " tetrahedrons (" << n_sig << " significant)" << endl;
            
        }
        
        ++n_clust;
    }
    for (list<ALPHACLUSTER*>::iterator it=clusters.begin(); it!=clusters.end(); ++it) delete *it;
    
    
    //! 7.) Ergebnis aufbereiten:
    if (lig_atoms.size() > 0) {
        float sradius = def_rad * def_rad;
        output_tets[0] = vector<stl_ptr<TETRAHEDRON> >();
        for (map<int,vector<stl_ptr<TETRAHEDRON> > >::iterator it=final_tets.begin(); it!=final_tets.end(); ++it)
        for (th_vec jt=it->second.begin(); jt!=it->second.end(); ++jt) {
            for (atoms_vec at=lig_atoms.begin(); at!=lig_atoms.end(); ++at) {
                vec3d<float> sm((*jt)->sphere_middle[0],(*jt)->sphere_middle[1],(*jt)->sphere_middle[2]);
                if (get_square_distance(sm,(*at)->coord) <= sradius) {
                    output_tets[0].push_back(*jt);
                    break;
                }
            }
        }
        output_vol[0] = 0.;
        for (th_vec jt=output_tets[0].begin(); jt!=output_tets[0].end(); ++jt) output_vol[0] += (*jt)->volume;
    } else {
        multimap<double,int,greater<double> > inv_vol;
        for (map<int,double>::iterator it=final_vol.begin(); it!=final_vol.end(); ++it) {
            inv_vol.insert(pair<double,int>(it->second,it->first));
        }
        int n_c = 0;
        for (multimap<double,int,greater<double> >::iterator it=inv_vol.begin(); it!=inv_vol.end(); ++it) {
            output_vol[n_c] = final_vol[it->second];
            output_tets[n_c] = final_tets[it->second];
            ++n_c;
        }
    }
    
    
    if (visualize) {
        cerr << "    => " << final_tets.size() << " clusters after volume filter" << endl;
        //! Shape der gefundenen Cavities rausschreiben:
        for (map<int,vector<stl_ptr<TETRAHEDRON> > >::iterator it=output_tets.begin(); it!=output_tets.end(); ++it) {
            vector<stl_ptr<L_FACET> > afc;
            mysolver.get_tet_shape(afc,it->second);
            ostringstream os; os << it->first;
            string sname = "cav_" + os.str() + ".off";
            cerr << "    => writing cavity " << it->first << " shape as " << sname << endl;
            mysolver.facets2off(sname.c_str(),afc);
            
            
            sname = "cav_" + os.str() + ".wrl";
            mysolver.facets2wrl(sname.c_str(),afc,shift);
        }
    }
    
    
    if (profile) {
        for (map<int,vector<stl_ptr<TETRAHEDRON> > >::iterator it=output_tets.begin(); it!=output_tets.end(); ++it) {
            ostringstream os; os << it->first;
            string sname = "cav_" + os.str() + "_profile.txt";
            get_cavity_profile2(mysolver,id2idx,it->second,alpha1,alpha2,sname);
    ///        break;///
        }
    }
    
    
    for (map<int,vector<stl_ptr<TETRAHEDRON> > >::iterator it=output_tets.begin(); it!=output_tets.end(); ++it) {
        //if (visualize) {
        //    string pname = "alpha_tets"; string_fu::add2s(pname,it->first);
        //    string fname = pname + ".py";
        //    mysolver.visualize_tetrahedrons(it->second,fname,pname);
        //    cout << " -> written cavity " << it->first << " tetrahedrons as " << fname << endl;
        //}
        tr1::unordered_set<int> pids;
        p_clust[it->first] = vector<stl_ptr<ATOM> >();
        for (th_vec jt=it->second.begin(); jt!=it->second.end(); ++jt) {
            if (pids.find((*jt)->p0->id) == pids.end()) {
                p_clust[it->first].push_back(points2atoms[(*jt)->p0->id]);
                pids.insert((*jt)->p0->id);
            }
            if (pids.find((*jt)->p1->id) == pids.end()) {
                p_clust[it->first].push_back(points2atoms[(*jt)->p1->id]);
                pids.insert((*jt)->p1->id);
            }
            if (pids.find((*jt)->p2->id) == pids.end()) {
                p_clust[it->first].push_back(points2atoms[(*jt)->p2->id]);
                pids.insert((*jt)->p2->id);
            }
            if (pids.find((*jt)->p3->id) == pids.end()) {
                p_clust[it->first].push_back(points2atoms[(*jt)->p3->id]);
                pids.insert((*jt)->p3->id);
            }
        }
    }
    
    if (surr_radius > 0.5) {
        float sura = surr_radius * surr_radius;
        for (map<int,vector<stl_ptr<ATOM> > >::iterator it=p_clust.begin(); it!=p_clust.end(); ++it) {
            tr1::unordered_set<int> at_ids;
            for (atoms_vec at=it->second.begin(); at!=it->second.end(); ++at) at_ids.insert((*at)->intern_id);
            vector<stl_ptr<ATOM> > adder;
            for (atoms_vec at=pv.begin(); at!=pv.end(); ++at) {
                if (at_ids.find((*at)->intern_id) != at_ids.end()) continue;
                for (atoms_vec bt=it->second.begin(); bt!=it->second.end(); ++bt) {
                    if (get_square_distance((*at)->coord,(*bt)->coord) < sura) {
                        adder.push_back(*at);
                        at_ids.insert((*at)->intern_id);
                        break;
                    }
                }
            }
            for (atoms_vec at=adder.begin(); at!=adder.end(); ++at) it->second.push_back(*at);
        }
    }
    
    //! eventuell noch CAVITY Objekte erzeugen, falls ueber write_pdb_cav() rausgeschrieben werden soll:
    if (get_cavity_obj) {
        for (chains_vec cv=protein->chains.begin(); cv!=protein->chains.end(); ++cv) if ((*cv)->n_aacids == 0) (*cv)->get_aacids();
        for (map<int,vector<stl_ptr<ATOM> > >::iterator it=p_clust.begin(); it!=p_clust.end(); ++it) {
            stl_ptr<CAVITY> cavity = new CAVITY();
            cavity->total_volume = output_vol[it->first];
            for (chains_vec jt=protein->chains.begin(); jt!=protein->chains.end();++jt) {
                bool breaker;
                for (aacids_map ut=(*jt)->aacids.begin(); ut!=(*jt)->aacids.end();++ut) {
                    for (atoms_vec at=ut->second->atoms.begin(); at!=ut->second->atoms.end(); ++at) {
                        breaker = false;
                        for (atoms_vec bt=it->second.begin(); bt!=it->second.end(); ++bt) {
                            if (((*at)->res_name == (*bt)->res_name) && ((*at)->res_number == (*bt)->res_number)) {
                                cavity->aacids.push_back(ut->second);
                                breaker = true;
                                break;
                            }
                        }
                        if (breaker) break;
                    }
                }
            }
            cavities.push_back(cavity);
        }
    }
    
    #if defined (_SPEEDTEST)
    t_tot.stop_t();
    cout << "total    = " << t_tot.get_time() << endl;
    #endif
}

void STRUCTURE::get_cavity_profile(DT_SOLVER& mysolver,vector<stl_ptr<TETRAHEDRON> >& pt,double alpha1,double alpha2,string& fname) {
    //! pt ist ein vector mit Tetraedern einer Cavity, welche von alpha1 und alpha2 begrenzt wird
    //! Wenn jetzt alpha1 stueckweise angehoben wird, sollte bestimmte Tetraeder in den alpha1-Komplex fallen
    //! und somit das Volumen der Cavity abnehmen. Eine sprunghafte Abnahme deutet auf eine abgeschnittene
    //! Subcavity hin:
    if (verb > 0) cout << "-> writing " << fname << endl;
    ofstream fout;
    fout.open(fname.c_str());
    double srad2 = alpha2 * alpha2;
    double vol = 0.;
    for (th_vec it=pt.begin(); it!=pt.end(); ++it) vol += (*it)->volume;
    fout << alpha1 << " " << vol << "\n";
    for (double la=alpha1+0.1; la<alpha2; la += 0.1) {
        double sla = la * la;
        vol = 0.;
        
        for (th_vec it=pt.begin(); it!=pt.end(); ++it) {
            if ((*it)->sphere_square_radius <= srad2 && (*it)->sphere_square_radius > sla) {
                vol += (*it)->volume;
            }
        }
        fout << la << " " << vol << "\n";
    }
    fout.close();
}

void STRUCTURE::get_cavity_profile2(DT_SOLVER& mysolver,tr1::unordered_map<int64_t,int>& id2idx,vector<stl_ptr<TETRAHEDRON> >& pt,
                                    double alpha1,double alpha2,string& fname) {
    //! alpha2 schrittweise absenken und jeweils clustern
    //! Punkte, an denen aus einem Cluster 2 werden markieren Subcavities:
    if (verb > 0) cout << "-> writing " << fname << endl;
    ofstream fout;
    fout.open(fname.c_str());
    double srad1 = alpha1 * alpha1;
    tr1::unordered_set<int64_t> allowed;
    for (th_vec it=pt.begin(); it!=pt.end(); ++it) allowed.insert((*it)->id);
    for (double la=alpha2; la>alpha1; la -= 0.1) {
        double sla = la * la;
        
        tr1::unordered_set<int64_t> forbidden;
        for (th_vec it=pt.begin(); it!=pt.end(); ++it) {
            if ((*it)->sphere_square_radius <= sla) {
                if ((*it)->p0n->sphere_square_radius > sla ||
                    (*it)->p1n->sphere_square_radius > sla ||
                    (*it)->p2n->sphere_square_radius > sla ||
                    (*it)->p3n->sphere_square_radius > sla) {
                    if ((*it)->p0n->sphere_square_radius <= srad1 ||
                        (*it)->p1n->sphere_square_radius <= srad1 ||
                        (*it)->p2n->sphere_square_radius <= srad1 ||
                        (*it)->p3n->sphere_square_radius <= srad1) {
                        forbidden.insert((*it)->id);
                    }
                }
            }
        }
        
        tr1::unordered_set<int64_t> visited;
        list<ALPHACLUSTER*> clusters;
        for (th_vec it=pt.begin(); it!=pt.end(); ++it) {
            if (visited.find((*it)->id) == visited.end() && (*it)->sphere_square_radius <= sla &&
                forbidden.find((*it)->id) == forbidden.end()) {
                visited.insert((*it)->id);
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
                    if (visited.find(t->id) != visited.end()) {
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
                        if (t->sphere_square_radius <= sla && t->sphere_square_radius > srad1 &&
                            allowed.find(t->id) != allowed.end() && forbidden.find(t->id) == forbidden.end()) {
                            clust->elements.insert(t->id);
                            visited.insert(t->id);
                            tstack.push(t->p0n);
                            tstack.push(t->p1n);
                            tstack.push(t->p2n);
                            tstack.push(t->p3n);
                        }
                    }
                }
            }
        }
        
        for (list<ALPHACLUSTER*>::iterator it=clusters.begin(); it!=clusters.end(); ++it) {
            tr1::unordered_set<int64_t> add;
            for (tr1::unordered_set<int64_t>::iterator jt=(*it)->elements.begin(); jt!=(*it)->elements.end(); ++jt) {
                if (forbidden.find(mysolver.delaunay_tetrahedrons[id2idx[*jt]]->p0n->id) != forbidden.end()) {
                    add.insert(mysolver.delaunay_tetrahedrons[id2idx[*jt]]->p0n->id);
                }
                if (forbidden.find(mysolver.delaunay_tetrahedrons[id2idx[*jt]]->p1n->id) != forbidden.end()) {
                    add.insert(mysolver.delaunay_tetrahedrons[id2idx[*jt]]->p1n->id);
                }
                if (forbidden.find(mysolver.delaunay_tetrahedrons[id2idx[*jt]]->p2n->id) != forbidden.end()) {
                    add.insert(mysolver.delaunay_tetrahedrons[id2idx[*jt]]->p2n->id);
                }
                if (forbidden.find(mysolver.delaunay_tetrahedrons[id2idx[*jt]]->p3n->id) != forbidden.end()) {
                    add.insert(mysolver.delaunay_tetrahedrons[id2idx[*jt]]->p3n->id);
                }
            }
            for (tr1::unordered_set<int64_t>::iterator jt=add.begin(); jt!=add.end(); ++jt) {
                (*it)->elements.insert(*jt);
            }
            
            //! Noch weitere forbidden dazunehmen!!! (die, die gar nicht an einen Cluster grenzen)
            
        }
        
        fout << la << " ";
        int n_sc = 0;
        for (list<ALPHACLUSTER*>::iterator it=clusters.begin(); it!=clusters.end(); ++it) {
            double vol = 0.;
            for (tr1::unordered_set<int64_t>::iterator jt=(*it)->elements.begin(); jt!=(*it)->elements.end(); ++jt) {
                vol += mysolver.delaunay_tetrahedrons[id2idx[*jt]]->volume;
            }
            if (vol < 20.) continue;
            ++n_sc;
            fout << vol << " ";
        }
        if (n_sc == 0) fout << 0.0;
        fout << "\n";
        for (list<ALPHACLUSTER*>::iterator it=clusters.begin(); it!=clusters.end(); ++it) delete *it;
    }
    fout.close();
}

void STRUCTURE::shorten_lig_names() {
    //neu 17.01.07: um ewig lange Namen bei langen Peptidketten zu vermeiden
    for (ligands_vec lig=ligands.begin(); lig!=ligands.end(); ++lig) {
        if ((*lig)->name.size() > 100) {
            if ((*lig)->atoms.size() > 0) {
                if ((*lig)->atoms[0]->type != 1) {
                    string name_buf;
                    string orig_name = (*lig)->name;
                    name_buf.assign((*lig)->name,0,3);
                    (*lig)->name = name_buf;
                    bool peptide = false;
                    for (atoms_vec atm=(*lig)->atoms.begin(); atm!=(*lig)->atoms.end(); ++atm) {
                        if ((*atm)->type == 1) {peptide = true; break;}
                    }
                    if (peptide) {
                        (*lig)->name += "_peptide";
                    } else (*lig)->name = "big_molecule";
                    if (peptide && (*lig)->atoms[(*lig)->atoms.size()-1]->type != 1) {
                        name_buf.assign(orig_name,orig_name.size()-4,3);
                        (*lig)->name += "_";
                        (*lig)->name += name_buf;
                    } 
                    if (orig_name[orig_name.size()-2] == '_') {
                        (*lig)->name += "_";
                        (*lig)->name += orig_name[orig_name.size()-1];
                    }
                } else {
                    string name_buf;
                    string orig_name = (*lig)->name;
                    (*lig)->name = "peptide";
                    if ((*lig)->atoms[(*lig)->atoms.size()-1]->type != 1) {
                        name_buf.assign(orig_name,orig_name.size()-4,3);
                        (*lig)->name += "_";
                        (*lig)->name += name_buf;
                    }
                    if (orig_name[orig_name.size()-2] == '_') {
                        (*lig)->name += "_";
                        (*lig)->name += orig_name[orig_name.size()-1];
                    }
                }
            }
        }
    }
}

void STRUCTURE::merge_ligands() {
    //!diese Funktion nur aufrufen, wenn kovalente Bindungen zwischen Liganden wahrscheinlich sind!!!
    //vector<int> con_vec[ligands.size()];
    EDnew::initialize();

    for (ligands_vec alig=ligands.begin(); alig!=ligands.end(); ++alig) {
        for (atoms_vec at=(*alig)->atoms.begin(); at!=(*alig)->atoms.end(); ++at) {
            //!neues ATOM_EXT Objekt:
            if ((*at)->ext.zero()) (*at)->ext = new ATOM_EXT();
            else (*at)->ext->clear();
            //!eventuell vorhandene Bindungen loeschen:
            (*at)->bonded_atoms.clear();
        }
        (*alig)->get_elements();
    }
    float dst;
    bool breaker;
    for (ligands_vec alig=ligands.begin(); alig!=ligands.end(); ++alig) {
        bool no_merge = false;
        for (ligands_vec blig=(alig+1); blig!=ligands.end(); ++blig) {
            breaker = false;
            for (atoms_vec at=(*alig)->atoms.begin(); at!=(*alig)->atoms.end(); ++at) {
                if ((*at)->element == "H") continue;
                for (atoms_vec bt=(*blig)->atoms.begin(); bt!=(*blig)->atoms.end(); ++bt) {
                    if ((*bt)->element == "H") continue;
                    if ((*at)->alt_loc_id != (*bt)->alt_loc_id) continue;

                    if ((*at)->res_number == (*bt)->res_number) { //!neu: 05.08.2010
                        no_merge = true;
                        breaker = false;
                        break;
                    }

                    dst = get_distance((*at)->coord,(*bt)->coord);
                    if (dst > EDnew::ED_max_dist) continue;
                    if (!EDnew::EDAccess::connection_is_possible((*at)->element,(*bt)->element)) continue;

                    if (dst < 1.0) { //! 28.07.2010: auf Clash checken (siehe 7ABP):
                        no_merge = true;
                        breaker = false;
                        break;
                    }

                    if (breaker) continue;

                    if (dst < 1.5 || EDnew::EDAccess::get_max_probability((*at)->element,(*bt)->element,dst) > 0.0005) {
                        breaker = true;
                    }
                    
                }
                if (breaker || no_merge) break;
            }
            if (breaker) {
                for (atoms_vec bt=(*blig)->atoms.begin(); bt!=(*blig)->atoms.end(); ++bt) {
                    stl_ptr<ATOM> atm(new ATOM(**bt));
                    (*alig)->atoms.push_back(atm);
                }
                (*alig)->name += "_" + (*blig)->name;
                blig->kill();
                ligands.erase(blig);
                alig--;
                break;
            }
        }
    }
}


void STRUCTURE::merge_peptides() {
    //!diese Funktion nur aufrufen, wenn kovalente Bindungen zu Peptiden wahrscheinlich sind und diese mit
    //!dem Liganden zusammengefasst werden sollen
    if (!protein) return;
    //vector<int> con_vec[ligands.size()];
    EDnew::initialize();

    for (ligands_vec alig=ligands.begin(); alig!=ligands.end(); ++alig) {
        for (atoms_vec at=(*alig)->atoms.begin(); at!=(*alig)->atoms.end(); ++at) {
            //!neues ATOM_EXT Objekt:
            if ((*at)->ext.zero()) (*at)->ext = new ATOM_EXT();
            else (*at)->ext->clear();
            //!eventuell vorhandene Bindungen loeschen:
            (*at)->bonded_atoms.clear();
        }
        (*alig)->get_elements();
        if ((*alig)->next_chain.zero()) continue; //!nur fuer in Frage kommende Peptidketten
    //    if ((*alig)->next_chain->atoms.size() > 140) { //!definitiv kein Peptid mehr
        if ((*alig)->next_chain->atoms.size() > 100) { //! 27.04.2010:  runtergesetzt
            (*alig)->next_chain->prev_lig = 0;
            (*alig)->next_chain = 0;
            continue;
        }
        for (atoms_vec at=(*alig)->next_chain->atoms.begin(); at!=(*alig)->next_chain->atoms.end(); ++at) {
            //!neues ATOM_EXT Objekt:
            if ((*at)->ext.zero()) (*at)->ext = new ATOM_EXT();
            else (*at)->ext->clear();
            //!eventuell vorhandene Bindungen loeschen:
            (*at)->bonded_atoms.clear();
        }
        (*alig)->next_chain->get_elements();
    }

    float dst;
    bool breaker;
    for (ligands_vec alig=ligands.begin(); alig!=ligands.end(); ++alig) {
        for (chains_vec blig=protein->chains.begin(); blig!=protein->chains.end(); ++blig) {
            if ((*blig)->prev_lig.zero()) continue;
            breaker = false;
            for (atoms_vec at=(*alig)->atoms.begin(); at!=(*alig)->atoms.end(); ++at) {
                for (atoms_vec bt=(*blig)->atoms.begin(); bt!=(*blig)->atoms.end(); ++bt) {
                    dst = get_distance((*at)->coord,(*bt)->coord);
                    if (dst < EDnew::ED_min_dist || dst > EDnew::ED_max_dist) continue;
                    if (!EDnew::EDAccess::connection_is_possible((*at)->element,(*bt)->element)) continue;

                    if (dst < 1.5 || EDnew::EDAccess::get_max_probability((*at)->element,(*bt)->element,dst) > 0.0001) {
                        breaker = true;
                        break;
                    }
                }
                if (breaker) break;
            }
            if (breaker) {
                for (atoms_vec bt=(*blig)->atoms.begin(); bt!=(*blig)->atoms.end(); ++bt) {
                    stl_ptr<ATOM> atm(new ATOM(**bt));
                    (*alig)->atoms.push_back(atm);
                }
                (*alig)->name += "_peptide";
                (*alig)->is_peptide = true;
                (*blig)->prev_lig = 0; //nicht nochmal diese Kette pruefen!!!
                alig--;
                break;
            }
        }
    }
    for (ligands_vec alig=ligands.begin(); alig!=ligands.end(); ++alig) {
        if (!((*alig)->is_peptide)) continue;
        for (ligands_vec blig=(alig+1); blig!=ligands.end(); ++blig) {
            breaker = false;
            for (atoms_vec at=(*alig)->atoms.begin(); at!=(*alig)->atoms.end(); ++at) {
                for (atoms_vec bt=(*blig)->atoms.begin(); bt!=(*blig)->atoms.end(); ++bt) {
                    dst = get_distance((*at)->coord,(*bt)->coord);
                    if (dst < EDnew::ED_min_dist || dst > EDnew::ED_max_dist) continue;
                    if (!EDnew::EDAccess::connection_is_possible((*at)->element,(*bt)->element)) continue;

                    if (dst < 1.5 || EDnew::EDAccess::get_max_probability((*at)->element,(*bt)->element,dst) > 0.0001) {
                        breaker = true;
                        break;
                    }
                }
                if (breaker) break;
            }
            if (breaker) {
                for (atoms_vec bt=(*blig)->atoms.begin(); bt!=(*blig)->atoms.end(); ++bt) {
                    stl_ptr<ATOM> atm(new ATOM(**bt));
                    (*alig)->atoms.push_back(atm);
                }
                (*alig)->name += "_" + (*blig)->name;
                blig->kill();
                ligands.erase(blig);
                alig--;
                break;
            }
        }
    }
}

void STRUCTURE::get_alt_loc_ligs(bool prot_mode) {
    //Fuer jeden Liganden mit alt_loc Eintrag einen entsprechenden alternativen Liganden erzeugen:
    vector< stl_ptr<LIGAND> > alt_ligs;
    vector<char> used_alt_locs;
    //!es muessen neue intern_id's vergeben werden:
    int cur_id = 1;
    if (prot_mode) { //!auf mol2_protein anwenden:
        if (mol2_protein.size() > 0) {
            cur_id = mol2_protein.back()->atoms.back()->intern_id + 1;
        }
        for (ligands_vec lig=mol2_protein.begin(); lig!=mol2_protein.end(); ++lig) {
            for (atoms_vec at=(*lig)->atoms.begin(); at!=(*lig)->atoms.end(); ++at) {
                if ((*at)->alt_loc_id != ' ') {
                    bool new_alt = true;
                    for (vector<char>::iterator cit=used_alt_locs.begin(); cit!=used_alt_locs.end(); ++cit) {
                        if ((*cit) == (*at)->alt_loc_id) new_alt = false;
                    }
                    if (new_alt) used_alt_locs.push_back((*at)->alt_loc_id);
                }
            }
            if (used_alt_locs.size() == 0) continue;
            for (vector<char>::iterator cit=used_alt_locs.begin(); cit!=used_alt_locs.end(); ++cit) {
                stl_ptr<LIGAND> new_lig(new LIGAND());
                stl_ptr<STRUCTURE> this_struct = this;
                new_lig->main_structure = this_struct;
                new_lig->name = (*lig)->name + "_" + *cit;
                new_lig->type = (*lig)->type;
                new_lig->charge_type = (*lig)->charge_type;
                for (atoms_vec at=(*lig)->atoms.begin(); at!=(*lig)->atoms.end(); ++at) {
                    if ((*at)->alt_loc_id == ' ' || (*at)->alt_loc_id == *cit) {
                        stl_ptr<ATOM> atm(new ATOM(**at));
                        atm->intern_id = cur_id; ++cur_id;
                        new_lig->atoms.push_back(atm);
                    }
                }
                alt_ligs.push_back(new_lig);
            }
            //jetzt den urspruenglichen Liganden aus mol2_protein loeschen:
            mol2_protein.erase(lig);
            lig--;
            used_alt_locs.clear();
        }
        //und jetzt die neuen Liganden an die alten anfuegen:
        for (ligands_vec lig=alt_ligs.begin(); lig!=alt_ligs.end(); ++lig) {
            mol2_protein.push_back(*lig);
        }
    } else { //!nur fuer ligands:
        if (ligands.size() > 0) {
            cur_id = ligands.back()->atoms.back()->intern_id + 1;
        }
        for (ligands_vec lig=ligands.begin(); lig!=ligands.end(); ++lig) {
            for (atoms_vec at=(*lig)->atoms.begin(); at!=(*lig)->atoms.end(); ++at) {
                if ((*at)->alt_loc_id != ' ') {
                    bool new_alt = true;
                    for (vector<char>::iterator cit=used_alt_locs.begin(); cit!=used_alt_locs.end(); ++cit) {
                        if ((*cit) == (*at)->alt_loc_id) new_alt = false;
                    }
                    if (new_alt) used_alt_locs.push_back((*at)->alt_loc_id);
                }
            }
            if (used_alt_locs.size() == 0) continue;
            for (vector<char>::iterator cit=used_alt_locs.begin(); cit!=used_alt_locs.end(); ++cit) {
                stl_ptr<LIGAND> new_lig(new LIGAND());
                stl_ptr<STRUCTURE> this_struct = this;
                new_lig->main_structure = this_struct;
                new_lig->name = (*lig)->name + "_" + *cit;
                new_lig->type = (*lig)->type;
                new_lig->charge_type = (*lig)->charge_type;
                for (atoms_vec at=(*lig)->atoms.begin(); at!=(*lig)->atoms.end(); ++at) {
                    if ((*at)->alt_loc_id == ' ' || (*at)->alt_loc_id == *cit) {
                        stl_ptr<ATOM> atm(new ATOM(**at));
                        atm->intern_id = cur_id; cur_id++;
                        new_lig->atoms.push_back(atm);
                    }
                }
                alt_ligs.push_back(new_lig);
            }
            //jetzt den urspruenglichen Liganden aus ligands loeschen:
            ligands.erase(lig);
            lig--;
            used_alt_locs.clear();
        }
        //und jetzt die neuen Liganden an die alten anfuegen:
        for (ligands_vec lig=alt_ligs.begin(); lig!=alt_ligs.end(); ++lig) {
            ligands.push_back(*lig);
        }
    }
}

void STRUCTURE::water2lig() {
    water_ligand = new LIGAND();
    for (waters_vec wt=waters.begin(); wt!=waters.end(); ++wt) {
        water_ligand->atoms.push_back((*wt)->atom);
    }
    waters.clear();
}

void STRUCTURE::metal2lig() {
    metal_ligand = new LIGAND();
    for (metals_vec wt=metals.begin(); wt!=metals.end(); ++wt) {
        metal_ligand->atoms.push_back((*wt)->atom);
    }
    metals.clear();
}

void STRUCTURE::align_by_protein(PROTEIN *ref,bool spatial,bool debug,bool only_C_alpha) {
    if (!protein) return;
    protein->align(ref,debug,only_C_alpha);

    if (spatial) {
        for (chains_vec it=protein->chains.begin(); it!=protein->chains.end();++it) {
            for (atoms_vec at=(*it)->atoms.begin(); at!=(*it)->atoms.end();++at) {
                (*at)->coord -= protein->optalign_trans[1];
                (*at)->coord *= protein->optalign_rotm;
                (*at)->coord += protein->optalign_trans[0];
            }
        }
        
        for (ligands_vec it=ligands.begin(); it!=ligands.end(); ++it) {
            for (atoms_vec at=(*it)->atoms.begin(); at!=(*it)->atoms.end();++at) {
                (*at)->coord -= protein->optalign_trans[1];
                (*at)->coord *= protein->optalign_rotm;
                (*at)->coord += protein->optalign_trans[0];
            }
        }
        for (waters_vec it=waters.begin(); it!=waters.end(); ++it) {
            (*it)->atom->coord -= protein->optalign_trans[1];
            (*it)->atom->coord *= protein->optalign_rotm;
            (*it)->atom->coord += protein->optalign_trans[0];
        }
        for (metals_vec it=metals.begin(); it!=metals.end(); ++it) {
            (*it)->atom->coord -= protein->optalign_trans[1];
            (*it)->atom->coord *= protein->optalign_rotm;
            (*it)->atom->coord += protein->optalign_trans[0];
        }
    }
}

float STRUCTURE::get_shape_sim(stl_ptr<STRUCTURE> const& ref) {
    if (!protein) return -1.;
    if (!ref->protein) return -1.;
    return protein->align2(ref->protein,1.4,true);
}

bool STRUCTURE::align_by_protein2(PROTEIN *ref,bool spatial,float const& max_displacement) {
    if (!protein) return false;
    
    if (protein->align2(ref,max_displacement) < -0.01) return false;
    
    if (spatial) {
        for (chains_vec it=protein->chains.begin(); it!=protein->chains.end();++it) {
            for (atoms_vec at=(*it)->atoms.begin(); at!=(*it)->atoms.end();++at) {
                (*at)->coord -= protein->optalign_trans[1];
                (*at)->coord *= protein->optalign_rotm;
                (*at)->coord += protein->optalign_trans[0];
            }
        }
        
        for (ligands_vec it=ligands.begin(); it!=ligands.end(); ++it) {
            for (atoms_vec at=(*it)->atoms.begin(); at!=(*it)->atoms.end();++at) {
                (*at)->coord -= protein->optalign_trans[1];
                (*at)->coord *= protein->optalign_rotm;
                (*at)->coord += protein->optalign_trans[0];
            }
        }
        for (waters_vec it=waters.begin(); it!=waters.end(); ++it) {
            (*it)->atom->coord -= protein->optalign_trans[1];
            (*it)->atom->coord *= protein->optalign_rotm;
            (*it)->atom->coord += protein->optalign_trans[0];
        }
        for (metals_vec it=metals.begin(); it!=metals.end(); ++it) {
            (*it)->atom->coord -= protein->optalign_trans[1];
            (*it)->atom->coord *= protein->optalign_rotm;
            (*it)->atom->coord += protein->optalign_trans[0];
        }
    }
    return true;
}

void STRUCTURE::align_by_cavity(STRUCTURE *ref) {
    if (ref->cavities.size() < 1 || cavities.size() < 1) return;
    cavities[0]->align(ref->cavities[0]);
    
    for (chains_vec it=protein->chains.begin(); it!=protein->chains.end();++it) {
        for (atoms_vec at=(*it)->atoms.begin(); at!=(*it)->atoms.end();++at) {
            (*at)->coord -= cavities[0]->optalign_trans[1];
            (*at)->coord *= cavities[0]->optalign_rotm;
            (*at)->coord += cavities[0]->optalign_trans[0];
        }
    }
    
    for (ligands_vec it=ligands.begin(); it!=ligands.end(); ++it) {
        for (atoms_vec at=(*it)->atoms.begin(); at!=(*it)->atoms.end();++at) {
            (*at)->coord -= cavities[0]->optalign_trans[1];
            (*at)->coord *= cavities[0]->optalign_rotm;
            (*at)->coord += cavities[0]->optalign_trans[0];
        }
    }
    for (waters_vec it=waters.begin(); it!=waters.end(); ++it) {
        (*it)->atom->coord -= cavities[0]->optalign_trans[1];
        (*it)->atom->coord *= cavities[0]->optalign_rotm;
        (*it)->atom->coord += cavities[0]->optalign_trans[0];
    }
    for (metals_vec it=metals.begin(); it!=metals.end(); ++it) {
        (*it)->atom->coord -= cavities[0]->optalign_trans[1];
        (*it)->atom->coord *= cavities[0]->optalign_rotm;
        (*it)->atom->coord += cavities[0]->optalign_trans[0];
    }
}


vec3d<float> STRUCTURE::get_ligand_diversity() {
    //! Sollte die Diversitaet aller Liganden in Bezug auf deren rmsd
    //! untereinander wiedergeben. => Es geht also um gleiche Liganden
    
    //! 1.) rmsd-werte fuer alle Kombinationen bestimmen:
    float total_val = 0.;
    int n_rmsds = 0;
    vector<float> s_vals;
    for (ligands_vec lig=ligands.begin(); lig!=ligands.end(); ++lig) {
        for (ligands_vec blig=(lig+1); blig!=ligands.end(); ++blig) {
            float curr_val = (*lig)->get_rmsd(*blig,false,false);
            if (curr_val < 0.) (*lig)->get_bk_rmsd(*blig,false,false);
            ++n_rmsds;
            total_val += curr_val;
            s_vals.push_back(curr_val);
        }
    }
    
    //! 2.) Mittelwert berechnen:
    float tval = n_rmsds + 0.0001;
    if (n_rmsds < 2) n_rmsds = 2;
    float m_val = total_val / n_rmsds;
    
    //! 3.) Standardabweichung berechnen:
    float vbuf = 0.;
    float bb;
    for (vector<float>::iterator it=s_vals.begin(); it!=s_vals.end(); ++it) {
        bb = m_val - *it;
        vbuf += bb * bb;
    }
    
    --n_rmsds;
    vbuf /= n_rmsds;
    vbuf = sqrt(vbuf);
    return vec3d<float>(m_val,vbuf,tval);
}


//==============================================================================================
//Definitionen fuer CELL:
//==============================================================================================
CELL::CELL() : empty(true), surf(false), real_surf(false),
               next_to_surf(false), burial(0), friends(0) {}

CELL::~CELL() {}



//==============================================================================================
//Definitionen fuer BURIAL_CELL:
//==============================================================================================
BURIAL_CELL::BURIAL_CELL() : empty(true), surf(false), next_to_surf(false),
                             real_surf(false), burial(0), friends(0) {}

BURIAL_CELL::~BURIAL_CELL() {}


//==============================================================================================
//Definitionen fuer STRUCTURE_GRID:
//==============================================================================================
STRUCTURE_GRID::STRUCTURE_GRID(vector< stl_ptr<ATOM> > &atms, float c_size, int add): cell_size(c_size),overadd(add) {
    for (atoms_vec at=atms.begin(); at!=atms.end(); ++at) {
        atoms.push_back(*at);
    }
    calc_properties();
    generate_grid();
    get_cell_middle();
}

STRUCTURE_GRID::STRUCTURE_GRID(PROTEIN *protein, float c_size, int add): cell_size(c_size),overadd(add) {
    for (chains_vec cit=protein->chains.begin(); cit!=protein->chains.end(); ++cit) {
        for (atoms_vec at=(*cit)->atoms.begin(); at!=(*cit)->atoms.end(); ++at) {
            atoms.push_back(*at);
        }
    }
    calc_properties();
    generate_grid();
    get_cell_middle();
}

STRUCTURE_GRID::~STRUCTURE_GRID() {
    grid.kill();
}

void STRUCTURE_GRID::calc_properties() {
    half_size = cell_size / 2.;
    min = atoms[0]->coord;
    max = atoms[0]->coord;
    for (atoms_vec at=atoms.begin(); at!=atoms.end(); ++at) { //Eckpunkte bestimmen
        if ((*at)->coord[0] < min[0]) min[0] = (*at)->coord[0];
        else if ((*at)->coord[0] > max[0]) max[0] = (*at)->coord[0];
        if ((*at)->coord[1] < min[1]) min[1] = (*at)->coord[1];
        else if ((*at)->coord[1] > max[1]) max[1] = (*at)->coord[1];
        if ((*at)->coord[2] < min[2]) min[2] = (*at)->coord[2];
        else if ((*at)->coord[2] > max[2]) max[2] = (*at)->coord[2];
    }
    float addi = half_size + (overadd * cell_size) + 0.0001;
    min[0] -= addi; max[0] += addi;
    min[1] -= addi; max[1] += addi;
    min[2] -= addi; max[2] += addi;
    max -= min;
    sizes[0] = int((max[0] / cell_size) + 0.5);
    sizes[1] = int((max[1] / cell_size) + 0.5);
    sizes[2] = int((max[2] / cell_size) + 0.5);
    max[0] = min[0] + (sizes[0] * cell_size);
    max[1] = min[1] + (sizes[1] * cell_size);
    max[2] = min[2] + (sizes[2] * cell_size);
    center = max; center -= min;
    center *= 0.5; center += min;
}

void STRUCTURE_GRID::generate_grid() {
    int bs[3] = {sizes[0],sizes[1],sizes[2]};
    grid = new GRID<3,CELL>(bs);
    for (g_iterator it=grid->begin(); it!=grid->end(); ++it) {
        it->index = it.get_index();
    }
}

void STRUCTURE_GRID::check_empty(float add_to_vdw) {
    
    atom_properties::initialize();
    int g_ind[3];
    CELL *cav;
    
    int curr_vec[3];
    float tot_vdw;
    float add_cell_rad = sqrt(3. * cell_size * cell_size) / 2.;
    float check_rad;
    
    //!ACHTUNG: Wasserstoffe werden dadurch beruecksichtigt, dass vdW-Radien eines united
    //!         atom Modells benutzt werden, bei denen moegliche Hs mit eingerechnet sind
    for (atoms_vec at=atoms.begin(); at!=atoms.end(); ++at) {
        get_index(g_ind,(*at)->coord);
        grid->value(g_ind).empty = false;
        tot_vdw = atom_properties::vdW_map[(*at)->element] + add_to_vdw + add_cell_rad;
        check_rad = tot_vdw * tot_vdw;
        
        int n_steps = int((tot_vdw / cell_size) + 0.5);
        
        for (int x=g_ind[0]-n_steps; x<=g_ind[0]+n_steps; ++x) {
            if (x<0 || x>=sizes[0]) continue;
            curr_vec[0] = x;
            for (int y=g_ind[1]-n_steps; y<=g_ind[1]+n_steps; ++y) {
                if (y<0 || y>=sizes[1]) continue;
                curr_vec[1] = y;
                for (int z=g_ind[2]-n_steps; z<=g_ind[2]+n_steps; ++z) {
                    if (z<0 || z>=sizes[2]) continue;
                    curr_vec[2] = z;
                    cav = &(grid->value(curr_vec));
                    if (get_square_distance(cav->middle,(*at)->coord) < check_rad) {
                        cav->empty = false;
                    }
                }
            }
        }
    }
    //!=> Fuer alle Zellen im vdW-Radius des Proteins ist 'empty' jetzt false
}

void STRUCTURE_GRID::get_next_to_surf() {
    //!darf erst nach check_empty aufgerufen werden (zumindest wenn die Ergebnisse
    //!sinvoll sein sollen)
    int g_ind[3];
    CELL *cav;
    int curr_vec[3];
    for (g_iterator it=grid->begin(); it!=grid->end(); ++it) {
        if (it->empty) continue;
        grid->get_coords(g_ind,it.get_index());
        for (int x=g_ind[0]-1; x<=g_ind[0]+1; ++x) {
            if (x<0 || x>=sizes[0]) continue;
            curr_vec[0] = x;
            for (int y=g_ind[1]-1; y<=g_ind[1]+1; ++y) {
                if (y<0 || y>=sizes[1]) continue;
                curr_vec[1] = y;
                for (int z=g_ind[2]-1; z<=g_ind[2]+1; ++z) {
                    if (z<0 || z>=sizes[2]) continue;
                    curr_vec[2] = z;
                    cav = &(grid->value(curr_vec));
                    if (cav->empty) {
                        it->surf = true;
                        cav->next_to_surf = true;
                    }
                }
            }
        }
    }
    //!=> Fuer alle non-empty Zellen in Nachbarschaft zu leeren Zellen ist 'surf' jetzt true
    //!=> Fuer alle empty-cells in Nachbarschaft zu einer surf-cell ist 'next_to_surf' jetzt true
}

void STRUCTURE_GRID::get_real_surf(int min_buried,int min_clust) {
    
    vec3d<int> v[26];
    
    v[0] = vec3d<int>(1,0,0); v[9] = vec3d<int>(-1,0,0);
    v[1] = vec3d<int>(1,1,0); v[10] = vec3d<int>(-1,1,0); v[18] = vec3d<int>(0,1,0);
    v[2] = vec3d<int>(1,-1,0); v[11] = vec3d<int>(-1,-1,0); v[19] = vec3d<int>(0,-1,0);
    
    v[3] = vec3d<int>(1,0,1); v[12] = vec3d<int>(-1,0,1); v[20] = vec3d<int>(0,0,1);
    v[4] = vec3d<int>(1,1,1); v[13] = vec3d<int>(-1,1,1); v[21] = vec3d<int>(0,1,1);
    v[5] = vec3d<int>(1,-1,1); v[14] = vec3d<int>(-1,-1,1); v[22] = vec3d<int>(0,-1,1);
    
    v[6] = vec3d<int>(1,0,-1); v[15] = vec3d<int>(-1,0,-1); v[23] = vec3d<int>(0,0,-1);
    v[7] = vec3d<int>(1,1,-1); v[16] = vec3d<int>(-1,1,-1); v[24] = vec3d<int>(0,1,-1);
    v[8] = vec3d<int>(1,-1,-1); v[17] = vec3d<int>(-1,-1,-1); v[25] = vec3d<int>(0,-1,-1);
    
    CELL *cav;
    int g_ind[3];
    
    int n_steps = int(16. / cell_size);
    
    for (g_iterator it=grid->begin(); it!=grid->end(); ++it) {
        if (it->next_to_surf) {
            int freeline = 0;
            for (int i=0; i<26; ++i) {
                get_index(g_ind,it->middle);
                vec3d<int> curr_vec(g_ind);
                
                bool freeway = true;
                for (int j=0; j<n_steps; ++j) {
                    curr_vec += v[i];
                    if (curr_vec[0] < 0 || curr_vec[1] < 0 || curr_vec[2] < 0 ||
                        curr_vec[0] >= sizes[0] || curr_vec[1] >= sizes[1] ||
                        curr_vec[2] >= sizes[2]) break;
                    cav = &(grid->value(curr_vec.get_vec()));
                    if (!cav->empty) {
                        freeway = false;
                        //!wenn Leerzellen ueber mehr als 5A, dann entsprechenden counter erhoehen:
                        if (get_square_distance(it->middle,cav->middle) > 16.) ++freeline;
                        break;
                    }
                }
                if (!freeway) {
                    it->burial += 1;
                }
            }
            if (freeline < 1) it->burial = 27;
        }
    }
    
    for (g_iterator it=grid->begin(); it!=grid->end(); ++it) {
        //if (it->burial > 16 && it->burial < 27) it->real_surf = true;
        if (it->burial > min_buried && it->burial < 27) it->real_surf = true; //min_buried=18 scheint besser zu performen
    }
    //!=> Alle surf-cells die im konvexen Teil der Proteinoberfl�che liegen und in mindestens eine
    //!   freie Linie ueber 4A haben sind jetzt 'real_surf' cells
    //!   (burial hat einen Wert zwischen 0 und 27
    //!    27 bedeutet dabei, dass in keine Richtung mehr als 4A Platz ist)
    
    //!Jetzt die real_surf_cells clustern:
    for (g_iterator it=grid->begin(); it!=grid->end(); ++it) {
        if (it->real_surf) {
            for (int i=0; i<26; ++i) {
                get_index(g_ind,it->middle);
                vec3d<int> curr_vec(g_ind);
                curr_vec += v[i];
                cav = &(grid->value(curr_vec.get_vec()));
                if (cav->real_surf) it->friends += 1;
            }
        }
    }
    for (g_iterator it=grid->begin(); it!=grid->end(); ++it) {
        if (it->real_surf) {
            if (it->friends < min_clust) it->real_surf = false;
        }
    }
    //! => Jetzt sind nur noch die Zellen real_surf, die an mindestens 8 andere real_surf-cells angrenzen
}

void STRUCTURE_GRID::get_buriedness() {
    //! Fuer alle next_to_surf:
    //! Die umgebende Schale durchgehen und alle empty's in einen vector schreiben.
    //! Dann fuer jede Zelle aus dem vector das gleiche tun (in einen neuen vector schreiben),
    //! allerdings nur fuer Zellen, die nicht schon im alten vector stehen.
    //! Insgesamt duerfen maximal  (8.879 / spacing)^3  Leerzellen vorhanden sein, um von lokaler
    //! Konkavitaet sprechen zu koennen.
    
    int g_ind[3];
    
    stl_ptr<CELL> cav;
    
    int steps = int((5. / cell_size) + 0.5);
    for (g_iterator it=grid->begin(); it!=grid->end(); ++it) {
        it->burial = 0;
        //it->friends = 0;
        if (!(it->next_to_surf)) continue;
        get_index(g_ind,it->middle);
        vec3d<int> curr_vec(g_ind);
        for (int x=-steps; x<=steps; ++x) {
            for (int y=-steps; y<=steps; ++y) {
                for (int z=-steps; z<=steps; ++z) {
                    curr_vec[0] = g_ind[0] + x;
                    curr_vec[1] = g_ind[1] + y;
                    curr_vec[2] = g_ind[2] + z;
                    cav = &(grid->value(curr_vec.get_vec()));
                    if (!(cav->empty)) it->burial += 1;
                }
            }
        }
    }
}

vec3d<int> STRUCTURE_GRID::get_sizes() {
    return vec3d<int>(sizes);
}

void STRUCTURE_GRID::get_index(int *index,vec3d<float> &coord) {
    index[0] = int((coord[0] - min[0]) / cell_size);
    index[1] = int((coord[1] - min[1]) / cell_size);
    index[2] = int((coord[2] - min[2]) / cell_size);
}

void STRUCTURE_GRID::get_cell_middle() {
    for (g_iterator it=grid->begin(); it!=grid->end(); ++it) { //Mittelpunkte aller Zellen bestimmen
        grid->get_coords(grid_pointer,it.get_index()); //!schreibt die Koordinaten der aktuellen Zelle in grid_pointer
        it->middle[0] = (grid_pointer[0] * cell_size) + half_size;
        it->middle[1] = (grid_pointer[1] * cell_size) + half_size;
        it->middle[2] = (grid_pointer[2] * cell_size) + half_size;
        it->middle += min;
    }
}



//==============================================================================================
//Definitionen fuer BURIAL_STRUCTURE_GRID:
//==============================================================================================
BURIAL_STRUCTURE_GRID::BURIAL_STRUCTURE_GRID(vector< stl_ptr<ATOM> > &atms, float c_size, int add): cell_size(c_size),overadd(add) {
    for (atoms_vec at=atms.begin(); at!=atms.end(); ++at) {
        atoms.push_back(*at);
    }
    calc_properties();
    generate_grid();
    get_cell_middle();
}

BURIAL_STRUCTURE_GRID::BURIAL_STRUCTURE_GRID(PROTEIN *protein, float c_size, int add): cell_size(c_size),overadd(add) {
    for (chains_vec cit=protein->chains.begin(); cit!=protein->chains.end(); ++cit) {
        for (atoms_vec at=(*cit)->atoms.begin(); at!=(*cit)->atoms.end(); ++at) {
            atoms.push_back(*at);
        }
    }
    calc_properties();
    generate_grid();
    get_cell_middle();
}

BURIAL_STRUCTURE_GRID::~BURIAL_STRUCTURE_GRID() {
    grid.kill();
}

void BURIAL_STRUCTURE_GRID::calc_properties() {
    half_size = cell_size / 2.;
    min = atoms[0]->coord;
    max = atoms[0]->coord;
    for (atoms_vec at=atoms.begin(); at!=atoms.end(); ++at) { //Eckpunkte bestimmen
        if ((*at)->coord[0] < min[0]) min[0] = (*at)->coord[0];
        else if ((*at)->coord[0] > max[0]) max[0] = (*at)->coord[0];
        if ((*at)->coord[1] < min[1]) min[1] = (*at)->coord[1];
        else if ((*at)->coord[1] > max[1]) max[1] = (*at)->coord[1];
        if ((*at)->coord[2] < min[2]) min[2] = (*at)->coord[2];
        else if ((*at)->coord[2] > max[2]) max[2] = (*at)->coord[2];
    }
    float addi = half_size + (overadd * cell_size) + 0.0001;
    min[0] -= addi; max[0] += addi;
    min[1] -= addi; max[1] += addi;
    min[2] -= addi; max[2] += addi;
    max -= min;
    sizes[0] = int((max[0] / cell_size) + 0.5);
    sizes[1] = int((max[1] / cell_size) + 0.5);
    sizes[2] = int((max[2] / cell_size) + 0.5);
    max[0] = min[0] + (sizes[0] * cell_size);
    max[1] = min[1] + (sizes[1] * cell_size);
    max[2] = min[2] + (sizes[2] * cell_size);
    center = max; center -= min;
    center *= 0.5; center += min;
}

void BURIAL_STRUCTURE_GRID::generate_grid() {
    int bs[3] = {sizes[0],sizes[1],sizes[2]};
    grid = new GRID<3,BURIAL_CELL>(bs);
}

void BURIAL_STRUCTURE_GRID::check_empty(float add_to_vdw) {
    
    atom_properties::initialize();
    int g_ind[3];
    BURIAL_CELL *cav;
    
    int curr_vec[3];
    float tot_vdw;
    float add_cell_rad = sqrt(3. * cell_size * cell_size) / 2.;
    float check_rad;
    
    //!ACHTUNG: Wasserstoffe werden dadurch beruecksichtigt, dass vdW-Radien eines united
    //!         atom Modells benutzt werden, bei denen moegliche Hs mit eingerechnet sind
    for (atoms_vec at=atoms.begin(); at!=atoms.end(); ++at) {
        get_index(g_ind,(*at)->coord);
        grid->value(g_ind).empty = false;
        tot_vdw = atom_properties::vdW_map[(*at)->element] + add_to_vdw + add_cell_rad;
        check_rad = tot_vdw * tot_vdw;
        
        int n_steps = int((tot_vdw / cell_size) + 0.5);
        
        for (int x=g_ind[0]-n_steps; x<=g_ind[0]+n_steps; ++x) {
            if (x<0 || x>=sizes[0]) continue;
            curr_vec[0] = x;
            for (int y=g_ind[1]-n_steps; y<=g_ind[1]+n_steps; ++y) {
                if (y<0 || y>=sizes[1]) continue;
                curr_vec[1] = y;
                for (int z=g_ind[2]-n_steps; z<=g_ind[2]+n_steps; ++z) {
                    if (z<0 || z>=sizes[2]) continue;
                    curr_vec[2] = z;
                    cav = &(grid->value(curr_vec));
                    if (get_square_distance(cav->middle,(*at)->coord) < check_rad) {
                        cav->empty = false;
                        cav->atoms.push_back(*at);
                    }
                }
            }
        }
    }
    //!=> Fuer alle Zellen im vdW-Radius des Proteins ist 'empty' jetzt false
}

void BURIAL_STRUCTURE_GRID::get_next_to_surf() {
    //!darf erst nach check_empty aufgerufen werden (zumindest wenn die Ergebnisse
    //!sinvoll sein sollen)
    int g_ind[3];
    BURIAL_CELL *cav;
    int curr_vec[3];
    for (g_iterator it=grid->begin(); it!=grid->end(); ++it) {
        if (it->empty) continue;
        grid->get_coords(g_ind,it.get_index());
        for (int x=g_ind[0]-1; x<=g_ind[0]+1; ++x) {
            if (x<0 || x>=sizes[0]) continue;
            curr_vec[0] = x;
            for (int y=g_ind[1]-1; y<=g_ind[1]+1; ++y) {
                if (y<0 || y>=sizes[1]) continue;
                curr_vec[1] = y;
                for (int z=g_ind[2]-1; z<=g_ind[2]+1; ++z) {
                    if (z<0 || z>=sizes[2]) continue;
                    curr_vec[2] = z;
                    cav = &(grid->value(curr_vec));
                    if (cav->empty) {
                        it->surf = true;
                        cav->next_to_surf = true;
                    }
                }
            }
        }
    }
    //!=> Fuer alle non-empty Zellen in Nachbarschaft zu leeren Zellen ist 'surf' jetzt true
    //!=> Fuer alle empty-cells in Nachbarschaft zu einer surf-cell ist 'next_to_surf' jetzt true
}

void BURIAL_STRUCTURE_GRID::get_buriedness() {
    vec3d<int> v[26];
    
    v[0] = vec3d<int>(1,0,0); v[9] = vec3d<int>(-1,0,0);
    v[1] = vec3d<int>(1,1,0); v[10] = vec3d<int>(-1,1,0); v[18] = vec3d<int>(0,1,0);
    v[2] = vec3d<int>(1,-1,0); v[11] = vec3d<int>(-1,-1,0); v[19] = vec3d<int>(0,-1,0);
    
    v[3] = vec3d<int>(1,0,1); v[12] = vec3d<int>(-1,0,1); v[20] = vec3d<int>(0,0,1);
    v[4] = vec3d<int>(1,1,1); v[13] = vec3d<int>(-1,1,1); v[21] = vec3d<int>(0,1,1);
    v[5] = vec3d<int>(1,-1,1); v[14] = vec3d<int>(-1,-1,1); v[22] = vec3d<int>(0,-1,1);
    
    v[6] = vec3d<int>(1,0,-1); v[15] = vec3d<int>(-1,0,-1); v[23] = vec3d<int>(0,0,-1);
    v[7] = vec3d<int>(1,1,-1); v[16] = vec3d<int>(-1,1,-1); v[24] = vec3d<int>(0,1,-1);
    v[8] = vec3d<int>(1,-1,-1); v[17] = vec3d<int>(-1,-1,-1); v[25] = vec3d<int>(0,-1,-1);
    
    BURIAL_CELL *cav;
    int g_ind[3];
    
    int n_steps = int(16. / cell_size);
    
    for (g_iterator it=grid->begin(); it!=grid->end(); ++it) {
        if (it->next_to_surf) {
            int freeline = 0;
            for (int i=0; i<26; ++i) {
                get_index(g_ind,it->middle);
                vec3d<int> curr_vec(g_ind);
                
                bool freeway = true;
                for (int j=0; j<n_steps; ++j) {
                    curr_vec += v[i];
                    if (curr_vec[0] < 0 || curr_vec[1] < 0 || curr_vec[2] < 0 ||
                        curr_vec[0] >= sizes[0] || curr_vec[1] >= sizes[1] ||
                        curr_vec[2] >= sizes[2]) break;
                    cav = &(grid->value(curr_vec.get_vec()));
                    if (!cav->empty) {
                        freeway = false;
                        //!wenn Leerzellen ueber mehr als 5A, dann entsprechenden counter erhoehen:
                        if (get_square_distance(it->middle,cav->middle) > 16.) ++freeline;
                        break;
                    }
                }
                if (!freeway) it->burial += 1;
            }
            if (freeline < 1) it->burial = 27;
        }
    }
    
    for (g_iterator it=grid->begin(); it!=grid->end(); ++it) {
        //if (it->burial > 16 && it->burial < 27) it->real_surf = true;
        if (it->burial > 18 && it->burial < 27) it->real_surf = true; //min_buried=18 scheint besser zu performen
    }
    
    for (g_iterator it=grid->begin(); it!=grid->end(); ++it) {
        if (it->real_surf) {
            for (int i=0; i<26; ++i) {
                get_index(g_ind,it->middle);
                vec3d<int> curr_vec(g_ind);
                curr_vec += v[i];
                cav = &(grid->value(curr_vec.get_vec()));
                if (cav->real_surf) it->friends += 1;
            }
        }
    }
    
    
    for (g_iterator it=grid->begin(); it!=grid->end(); ++it) {
        //if (!it->empty) {
        if (it->surf) {
            for (int i=0; i<26; ++i) {
                get_index(g_ind,it->middle);
                vec3d<int> curr_vec(g_ind);
                curr_vec += v[i];
                if (curr_vec[0] < 0 || curr_vec[1] < 0 || curr_vec[2] < 0 ||
                    curr_vec[0] >= sizes[0] || curr_vec[1] >= sizes[1] ||
                    curr_vec[2] >= sizes[2]) continue;
                cav = &(grid->value(curr_vec.get_vec()));
                if (!cav->empty) {
                    it->burial += 38;
                } else {
                    it->burial += cav->burial;
                    if (cav->friends > 7) it->burial += 2;
                }
            }
            //Jetzt den enthaltenen atomen zuweisen:
            for (atoms_vec at=it->atoms.begin(); at!=it->atoms.end(); ++at) {
                (*at)->ext->buriedness += it->burial;
            }
        }
    }
    //! Alle Atome haben jetzt einen buriedness-Wert
}

vec3d<int> BURIAL_STRUCTURE_GRID::get_sizes() {
    return vec3d<int>(sizes);
}

void BURIAL_STRUCTURE_GRID::get_index(int *index,vec3d<float> &coord) {
    index[0] = int((coord[0] - min[0]) / cell_size);
    index[1] = int((coord[1] - min[1]) / cell_size);
    index[2] = int((coord[2] - min[2]) / cell_size);
}

void BURIAL_STRUCTURE_GRID::get_cell_middle() {
    for (g_iterator it=grid->begin(); it!=grid->end(); ++it) { //Mittelpunkte aller Zellen bestimmen
        grid->get_coords(grid_pointer,it.get_index()); //!schreibt die Koordinaten der aktuellen Zelle in grid_pointer
        it->middle[0] = (grid_pointer[0] * cell_size) + half_size;
        it->middle[1] = (grid_pointer[1] * cell_size) + half_size;
        it->middle[2] = (grid_pointer[2] * cell_size) + half_size;
        it->middle += min;
    }
}


//==============================================================================================
//Definitionen fuer COORD_GRID:
//==============================================================================================

COORD_GRID::COORD_GRID(vector< stl_ptr<ATOM> > &atoms, float spa, float add_size) {
    spacing = spa;
    min = atoms[0]->coord;
    max = atoms[0]->coord;
    
    for (atoms_vec at=atoms.begin(); at!=atoms.end(); ++at) { //zunaechst die Eckpunkte bestimmen
        if ((*at)->coord[0] < min[0]) min[0] = (*at)->coord[0];
        else if ((*at)->coord[0] > max[0]) max[0] = (*at)->coord[0];
        if ((*at)->coord[1] < min[1]) min[1] = (*at)->coord[1];
        else if ((*at)->coord[1] > max[1]) max[1] = (*at)->coord[1];
        if ((*at)->coord[2] < min[2]) min[2] = (*at)->coord[2];
        else if ((*at)->coord[2] > max[2]) max[2] = (*at)->coord[2];
    }
    
    min[0] -= add_size; max[0] += add_size;
    min[1] -= add_size; max[1] += add_size;
    min[2] -= add_size; max[2] += add_size;
    max -= min;
    sizes[0] = int((max[0] / spacing)+0.5);
    sizes[1] = int((max[1] / spacing)+0.5);
    sizes[2] = int((max[2] / spacing)+0.5);
    max[0] = min[0] + (sizes[0] * spacing);
    max[1] = min[1] + (sizes[1] * spacing);
    max[2] = min[2] + (sizes[2] * spacing);
    vec3d<float> z_vec(0.,0.,0.);
    int bs[3] = {1+sizes[0],1+sizes[1],1+sizes[2]};
    grid.resize(bs,z_vec);
    vec3d<float> hvec;
    for (int z=0; z<=sizes[2]; ++z) {
        bs[2] = z;
        hvec[2] = min[2] + z*spacing;
        for (int y=0; y<=sizes[1]; ++y) {
            bs[1] = y;
            hvec[1] = min[1] + y*spacing;
            for (int x=0; x<=sizes[0]; ++x) {
                bs[0] = x;
                hvec[0] = min[0] + x*spacing;
                grid.value(bs) = hvec;
            }
        }
    }
}

COORD_GRID::COORD_GRID(vec3d<float> &mi,vec3d<float> &ma,float spa,float add_size) {
    spacing = spa;
    min = mi;
    max = ma;
    min[0] -= add_size; max[0] += add_size;
    min[1] -= add_size; max[1] += add_size;
    min[2] -= add_size; max[2] += add_size;
    max -= min;
    sizes[0] = int((max[0] / spacing)+0.5);
    sizes[1] = int((max[1] / spacing)+0.5);
    sizes[2] = int((max[2] / spacing)+0.5);
    max[0] = min[0] + (sizes[0] * spacing);
    max[1] = min[1] + (sizes[1] * spacing);
    max[2] = min[2] + (sizes[2] * spacing);
    vec3d<float> z_vec(0.,0.,0.);
    int bs[3] = {1+sizes[0],1+sizes[1],1+sizes[2]};
    grid.resize(bs,z_vec);
    vec3d<float> hvec;
    for (int z=0; z<=sizes[2]; ++z) {
        bs[2] = z;
        hvec[2] = min[2] + z*spacing;
        for (int y=0; y<=sizes[1]; ++y) {
            bs[1] = y;
            hvec[1] = min[1] + y*spacing;
            for (int x=0; x<=sizes[0]; ++x) {
                bs[0] = x;
                hvec[0] = min[0] + x*spacing;
                grid.value(bs) = hvec;
            }
        }
    }
}

COORD_GRID::~COORD_GRID() {}

