
//============================================================================
// atom_GN.cpp -*- C++ -*-; objects representing atoms and related stuff
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


#include"atom_GN.h"


//==============================================================================================
//Definitionen fr ATOM:
//==============================================================================================

ATOM::ATOM() : type(5),sybyl_type("XX"),intern_type("X"),alt_loc_id(' '),res_name("X"),chain_id(' '),res_number(0),element("X"),
               charge(" 0"),is_ter(false),is_model(false),is_endmdl(false),dict_type(" "),bond_ind(0) {}

ATOM::ATOM(ATOM const& atom) :
    id(atom.id), intern_id(atom.intern_id), name(atom.name),
    type(atom.type), sybyl_type(atom.sybyl_type), intern_type(atom.intern_type), alt_loc_id(atom.alt_loc_id),
    res_name(atom.res_name), chain_id(atom.chain_id), res_number(atom.res_number),
    alt_res(atom.alt_res), coord(atom.coord), occupancy(atom.occupancy),
    b_factor(atom.b_factor), element(atom.element), charge(atom.charge),
    is_ter(atom.is_ter), is_model(atom.is_model), is_endmdl(atom.is_endmdl),
    model_number(atom.model_number), dict_type(atom.dict_type), bond_ind(atom.bond_ind),
    full_res_name(atom.full_res_name) {
    if (!atom.ext.zero()) {
        ext = new ATOM_EXT(*(atom.ext));
    }
}
    //!bonded_atoms werden nicht bernommen

ATOM& ATOM::operator=(ATOM const& atom) {
    id = atom.id;
    intern_id = atom.intern_id;
    name = atom.name;
    type = atom.type;
    sybyl_type = atom.sybyl_type;
    intern_type = atom.intern_type;
    alt_loc_id = atom.alt_loc_id;
    res_name = atom.res_name; 
    chain_id = atom.chain_id;
    res_number = atom.res_number;
    alt_res = atom.alt_res;
    coord = atom.coord;
    occupancy = atom.occupancy;
    b_factor = atom.b_factor;
    element = atom.element;
    charge = atom.charge;
    is_ter = atom.is_ter;
    is_model = atom.is_model;
    is_endmdl = atom.is_endmdl;
    model_number = atom.model_number;
    dict_type = atom.dict_type;
    bond_ind = atom.bond_ind;
    full_res_name = atom.full_res_name;
    if (!atom.ext.zero()) {
        ext = new ATOM_EXT(*(atom.ext));
    }
    //!bonded_atoms werden nicht bernommen
    return *this;
}

ATOM::~ATOM() {
    if (!ext.zero()) ext.kill();
}

void ATOM::remove_ext() {
    if (!ext.zero()) ext.zero_kill();
}

bool ATOM::operator>(ATOM const& rechts) {
    if (intern_id > rechts.intern_id) return true;
    return false;
}

bool ATOM::operator<(ATOM const& rechts) {
    if (intern_id < rechts.intern_id) return true;
    return false;
}

bool ATOM::operator==(ATOM const& rechts) {
    if (intern_id == rechts.intern_id) return true;
    return false;
}

bool ATOM::operator!=(ATOM const& rechts) {
    if (intern_id != rechts.intern_id) return true;
    return false;
}

void ATOM::get_compare_number() {
    if (ext.zero()) {
        ext = new ATOM_EXT(); // GN: 07.01.2011  warum nicht einfach anlegen statt error
        //cerr << c_message<cERROR>("ATOM::get_compare_number --> atoms have no EXT object") << endl;
        //exit(1);
    }
    ext->compare_number = 0;
    for (atoms_vec bt=bonded_atoms.begin(); bt!=bonded_atoms.end(); ++bt) {
        if ((*bt)->element == "C") ext->compare_number += 1;
        else if ((*bt)->element == "H") continue;
        else if ((*bt)->element == "N") ext->compare_number += 5;
        else if ((*bt)->element == "O") ext->compare_number += 21;
        else if ((*bt)->element == "S") ext->compare_number += 85;
        else if ((*bt)->element == "P") ext->compare_number += 341;
        else if ((*bt)->element == "F") ext->compare_number += 1365;
        else if ((*bt)->element == "Cl") ext->compare_number += 5461;
        else if ((*bt)->element == "Br") ext->compare_number += 21845;
        else if ((*bt)->element == "I") ext->compare_number += 87381;
    }
}

void ATOM::get_compare_number2() {
    if (ext.zero()) {
        ext = new ATOM_EXT(); // GN: 07.01.2011  warum nicht einfach anlegen statt error
        //cerr << c_message<cERROR>("ATOM::get_compare_number --> atoms have no EXT object") << endl;
        //exit(1);
    }
    ext->compare_number2 = 0;
    for (atoms_vec bt=bonded_atoms.begin(); bt!=bonded_atoms.end(); ++bt) {
        if ((*bt)->element == "H") continue; //!wichtig!: sonst wird bei ungleicher H-Anzahl der rmsd -1
        ext->compare_number2 += (*bt)->ext->compare_number;
    }
}

bool ATOM::is_equal(ATOM &ref,vector<stl_ptr<ATOM> > &a_roots,vector<stl_ptr<ATOM> > &b_roots,
                    tr1::unordered_map<int,tr1::unordered_set<int> > &known,
                    tr1::unordered_map<int,tr1::unordered_set<int> > &negknown) {
    if (ext->compare_number != ref.ext->compare_number) return false;
    if (ext->compare_number2 != ref.ext->compare_number2) return false;
    if (ref.intern_type != "X" && ext->is_ring != ref.ext->is_ring) return false;

    if (known.find(intern_id) != known.end()) {
        if (known[intern_id].find(ref.intern_id) != known[intern_id].end()) return true;
    }

    if (negknown.find(intern_id) != negknown.end()) {
        if (negknown[intern_id].find(ref.intern_id) != negknown[intern_id].end()) return false;
    }
    
    a_roots.push_back(stl_ptr<ATOM>(*this));
    b_roots.push_back(stl_ptr<ATOM>(ref));

    //cerr << endl << "here for " << name << " and " << ref.name << endl;
    
    //----------------------------------------------------------------------------
    //! GN 07.01.2011: Folgender Mehraufwand ist noetig fuer AJUYIN / ACIHIE,
    //! da sonst Atome als aequivalent erkannt werden, die nicht aequivalent sind:
    int needed = 0;
    int got = 0;
    for (atoms_vec at=bonded_atoms.begin(); at!=bonded_atoms.end(); ++at) {
        if ((*at)->element == "H") continue;

        bool conti = false;
        for (atoms_vec pt=a_roots.begin(); pt!=a_roots.end(); ++pt) {
            if (*at == *pt) { //!dieses Atom wurde in diesem Durchlauf bereits verglichen
                conti = true;
                break;
            }
        }
        if (conti) continue;
        ++needed;
    }
    for (atoms_vec bt=ref.bonded_atoms.begin(); bt!=ref.bonded_atoms.end(); ++bt) {
        if ((*bt)->element == "H") continue;
        bool conti = false;
        for (atoms_vec pt2=b_roots.begin(); pt2!=b_roots.end(); ++pt2) {
            if (*bt == *pt2) {
                conti = true;
                break;
            }
        }
        if (conti) continue;
        ++got;
    }
    if (needed != got) {
        negknown[intern_id].insert(ref.intern_id);
        a_roots.pop_back();
        b_roots.pop_back();
        return false;
    }
    //----------------------------------------------------------------------------

    for (atoms_vec at=bonded_atoms.begin(); at!=bonded_atoms.end(); ++at) {
        if ((*at)->element == "H") continue;
        
        bool conti = false;
        for (atoms_vec pt=a_roots.begin(); pt!=a_roots.end(); ++pt) {
            if (*at == *pt) { //!dieses Atom wurde in diesem Durchlauf bereits verglichen
                conti = true;
                break;
            }
        }
        if (conti) continue;

        //cerr << "searching match for " << (*at)->name << endl;

        bool has_eq = false;
        for (atoms_vec bt=ref.bonded_atoms.begin(); bt!=ref.bonded_atoms.end(); ++bt) {
            if ((*bt)->element == "H") continue;
            conti = false;
            for (atoms_vec pt2=b_roots.begin(); pt2!=b_roots.end(); ++pt2) {
                if (*bt == *pt2) {
                    conti = true;
                    break;
                }
            }
            if (conti) continue;
            if ((*at)->is_equal(**bt,a_roots,b_roots,known,negknown)) {

                //cerr << (*at)->name << " equals " << (*bt)->name << endl;

                has_eq = true;
                break;
            }
        }
        if (!has_eq) {
            a_roots.pop_back();
            b_roots.pop_back();
            negknown[intern_id].insert(ref.intern_id);
            return false;
        }
    }

    known[intern_id].insert(ref.intern_id);
    return true;
}

void ATOM::get_element(bool reset_names) {
    string ele_help = element;
    ostringstream os2;
    os2 << name;
    element = os2.str();
    if (!(sybyl_type == "XX")) {
        if (sybyl_type ==  "H") {element = "H"; return;}
        else if (sybyl_type.size() > 1) {
            if (sybyl_type[0] ==  'C' && sybyl_type[1] == '.') {element = "C"; return;}
            else if (sybyl_type[0] ==  'N' && sybyl_type[1] == '.') {element = "N"; return;}
            else if (sybyl_type[0] ==  'O' && sybyl_type[1] == '.') {element = "O"; return;}
            else if (sybyl_type[0] ==  'S' && sybyl_type[1] == '.') {element = "S"; return;}
            else if (sybyl_type[0] ==  'P' && sybyl_type[1] == '.') {element = "P"; return;}
            else if (sybyl_type[0] ==  'N' && sybyl_type[1] == '.') {element = "N"; return;}
            else if (sybyl_type[1] == '.') {element = sybyl_type[0]; return;}
            else {
                if (sybyl_type[1] < 91) sybyl_type[1] += 32;
                if (sybyl_type ==  "Lp") {
                    sybyl_type = "LP";
                    element = sybyl_type;
                    return;
                }
                element.assign(sybyl_type,0,2);
                if (element == "Cl" || element == "Br" || element == "Du") return;
                else type = 3;
                return;
            }
        } else {
            if (sybyl_type ==  "X") { //den hat fconv dann schonmal auf X gesetzt
                if (name.size() > 0) {
                    if (name[0] == 'C') {element = "C"; return;}
                    if (name[0] == 'N') {element = "N"; return;}
                    if (name[0] == 'O') {element = "O"; return;}
                    if (name[0] == 'S') {element = "S"; return;}
                    if (name[0] == 'P') {element = "P"; return;}
                }
            } else {
                element = sybyl_type;
                if (!(element == "C" || element == "N" || element == "O" || element == "S" ||
                      element == "P" || element == "F" || element == "I" || element == "B")) type = 3;
                return;
            }
        }
        
        /*
        else if (sybyl_type ==  "C.1") {element = "C"; return;}
        else if (sybyl_type ==  "C.2") {element = "C"; return;}
        else if (sybyl_type ==  "C.3") {element = "C"; return;}
        else if (sybyl_type ==  "C.cat") {element = "C"; return;}
        else if (sybyl_type ==  "C.ar") {element = "C"; return;}
        else if (sybyl_type ==  "N.1") {element = "N"; return;}
        else if (sybyl_type ==  "N.2") {element = "N"; return;}
        else if (sybyl_type ==  "N.3") {element = "N"; return;}
        else if (sybyl_type ==  "N.pl3") {element = "N"; return;}
        else if (sybyl_type ==  "N.ar") {element = "N"; return;}
        else if (sybyl_type ==  "N.am") {element = "N"; return;}
        else if (sybyl_type ==  "N.4") {element = "N"; return;}
        else if (sybyl_type[0] == 'O') {element = "O"; return;}
        else if (sybyl_type ==  "S.2") {element = "S"; return;}
        else if (sybyl_type ==  "S.3") {element = "S"; return;}
        else if (sybyl_type ==  "S.o") {element = "S"; return;}
        else if (sybyl_type ==  "S.o2") {element = "S"; return;}
        else if (sybyl_type ==  "P.3") {element = "P"; return;}
        else if (sybyl_type ==  "Cl") {element = "Cl"; return;}
        else if (sybyl_type ==  "Cl.0") {element = "Cl"; return;}
        else if (sybyl_type ==  "Cl.i") {element = "Cl"; return;}
        else if (sybyl_type ==  "CL") {element = "Cl"; return;}
        else if (sybyl_type ==  "Br") {element = "Br"; return;}
        else if (sybyl_type ==  "Br.0") {element = "Br"; return;}
        else if (sybyl_type ==  "Br.i") {element = "Br"; return;}
        else if (sybyl_type ==  "BR") {element = "Br"; return;}
        else if (sybyl_type ==  "I") {element = "I"; return;}
        else if (sybyl_type ==  "I.0") {element = "I"; return;}
        else if (sybyl_type ==  "I.i") {element = "I"; return;}
        else if (sybyl_type ==  "F") {element = "F"; return;}
        else if (sybyl_type ==  "F.0") {element = "F"; return;}
        else if (sybyl_type ==  "F.i") {element = "F"; return;}
        else if (sybyl_type ==  "Hal") {element = "Hal"; return;}
        else if (sybyl_type ==  "Het") {element = "Het"; return;}
        else if (sybyl_type ==  "Hev") {element = "Hev"; return;}
        else if (sybyl_type ==  "LP") {element = "LP"; return;} //!wichtig
        else if (sybyl_type ==  "Du") {element = "Du"; return;}
        else if (sybyl_type ==  "DU") {element = "Du"; return;}
        
        else if (sybyl_type ==  "Cu") {element = "Cu"; type = 3; return;}
        else if (sybyl_type ==  "Fe") {element = "Fe"; type = 3; return;}
        else if (sybyl_type ==  "Zn") {element = "Zn"; type = 3; return;}
        else if (sybyl_type ==  "Mg") {element = "Mg"; type = 3; return;}
        else if (sybyl_type ==  "Ni") {element = "Ni"; type = 3; return;}
        else if (sybyl_type ==  "Ca") {element = "Ca"; type = 3; return;}
        
        else if (sybyl_type ==  "CU") {element = "Cu"; type = 3; return;}
        else if (sybyl_type ==  "FE") {element = "Fe"; type = 3; return;}
        else if (sybyl_type ==  "ZN") {element = "Zn"; type = 3; return;}
        else if (sybyl_type ==  "MG") {element = "Mg"; type = 3; return;}
        else if (sybyl_type ==  "NI") {element = "Ni"; type = 3; return;}
        else if (sybyl_type ==  "CA") {element = "Ca"; type = 3; return;}
        
        else if (sybyl_type ==  "Co") {element = "Co"; type = 3; return;}
        else if (sybyl_type ==  "Co.oh") {element = "Co"; type = 3; return;}
        else if (sybyl_type ==  "Hg") {element = "Hg"; type = 3; return;}
        else if (sybyl_type ==  "Mn") {element = "Mn"; type = 3; return;}
        else if (sybyl_type ==  "Na") {element = "Na"; type = 3; return;}
        
        else if (sybyl_type ==  "NA") {element = "Na"; type = 3; return;}
        
        else if (sybyl_type ==  "K") {element = "K"; type = 3; return;}
        else if (sybyl_type ==  "Si") {element = "Si"; type = 3; return;}
        else if (sybyl_type ==  "Se") {element = "Se"; type = 3; return;}
        else if (sybyl_type ==  "B") {element = "B"; return;}
        else if (sybyl_type ==  "Al") {element = "Al"; type = 3; return;}
        
        else if (sybyl_type ==  "AL") {element = "Al"; type = 3; return;}
        
        else if (sybyl_type ==  "Li") {element = "Li"; type = 3; return;}
        else if (sybyl_type ==  "Cr") {element = "Cr"; type = 3; return;}
        else if (sybyl_type ==  "Cr.th") {element = "Cr"; type = 3; return;}
        else if (sybyl_type ==  "Cr.oh") {element = "Cr"; type = 3; return;}
        else if (sybyl_type ==  "Mo") {element = "Mo"; type = 3; return;}
        else if (sybyl_type ==  "Sn") {element = "Sn"; type = 3; return;}
        else if (sybyl_type[0] == 'C') {element = "C"; return;}
        else if (sybyl_type[0] == 'N') {element = "N"; return;}
        else if (sybyl_type[0] == 'S') {element = "S"; return;}
        else if (sybyl_type[0] == 'P') {element = "P"; return;}
        else if (sybyl_type[0] == 'H') {element = "H"; return;}
        else if (sybyl_type ==  "X") { //den hat fconv dann schonmal auf X gesetzt
            if (name.size() > 0) {
                if (name[0] == 'C') {element = "C"; return;}
                if (name[0] == 'N') {element = "N"; return;}
                if (name[0] == 'O') {element = "O"; return;}
                if (name[0] == 'S') {element = "S"; return;}
                if (name[0] == 'P') {element = "P"; return;}
            }
        }
        */
    } else {
        if (res_name == "HOH" || res_name == "H2O") { //! neu: 10.10.2008
            element = "O"; type = 4; return;
        }
        
        if (ele_help != "X") {
            if (ele_help.size() == 2) {
                if (ele_help == "C " || ele_help == "O " || ele_help == "N " ||
                    ele_help == "P " || ele_help == "S " || ele_help == "H " ||
                    ele_help == "F " || ele_help == "I " || ele_help == "B " ||
                    ele_help == "K ") {element.assign(1,ele_help[0]); return;}
                else if (ele_help == " C" || ele_help == " O" || ele_help == " N" ||
                        ele_help == " P" || ele_help == " S" || ele_help == " H" ||
                        ele_help == " F" || ele_help == " I" || ele_help == " B" ||
                        ele_help == " K") {element.assign(1,ele_help[1]); return;}
                else if (ele_help == "Cl" || ele_help == "Br") {element = ele_help; return;}
                else if (ele_help == "Na" || ele_help == "Mg" || ele_help == "Ca" ||
                        ele_help == "Al" || ele_help == "Si" || ele_help == "Cr" ||
                        ele_help == "Mn" || ele_help == "Fe" || ele_help == "Co" ||
                        ele_help == "Ni" || ele_help == "Cu" || ele_help == "Zn" ||
                        ele_help == "Se" || ele_help == "Ag" || ele_help == "Au" ||
                        ele_help == "Sn" || ele_help == "Pt" || ele_help == "Hg" ||
                        ele_help == "Pb" || ele_help == "Bi") {
                    element = ele_help;
                    type = 3;
                    return;
                }
                else if (ele_help == "CL") {element = "Cl"; return;}
                else if (ele_help == "BR") {element = "Br"; return;}
                else if (ele_help == "MG") {element = "Mg"; type = 3; return;}
                else if (ele_help == "CA") {element = "Ca"; type = 3; return;}
                else if (ele_help == "MN") {element = "Mn"; type = 3; return;}
                else if (ele_help == "CR") {element = "Cr"; type = 3; return;}
                else if (ele_help == "AL") {element = "Al"; type = 3; return;}
                else if (ele_help == "SI") {element = "Si"; type = 3; return;}
                else if (ele_help == "NA") {element = "Na"; type = 3; return;}
                else if (ele_help == "FE") {element = "Fe"; type = 3; return;}
                else if (ele_help == "CO") {element = "Co"; type = 3; return;}
                else if (ele_help == "NI") {element = "Ni"; type = 3; return;}
                else if (ele_help == "CU") {element = "Cu"; type = 3; return;}
                else if (ele_help == "ZN") {element = "Zn"; type = 3; return;}
                else if (ele_help == "AG") {element = "Ag"; type = 3; return;}
                else if (ele_help == "AU") {element = "Au"; type = 3; return;}
                else if (ele_help == "SN") {element = "Sn"; type = 3; return;}
                else if (ele_help == "PT") {element = "Pt"; type = 3; return;}
                else if (ele_help == "HG") {element = "Hg"; type = 3; return;}
                else if (ele_help == "PB") {element = "Pb"; type = 3; return;}
                else if (ele_help == "BI") {element = "Bi"; type = 3; return;}
                else if (ele_help == "XE") {element = "Xe"; return;}
            } else if (ele_help.size() == 1) { //eigentlich doch immer 2 ???
                if (ele_help == "C" || ele_help == "O" || ele_help == "N" ||
                    ele_help == "P" || ele_help == "S" || ele_help == "H" ||
                    ele_help == "F" || ele_help == "I" || ele_help == "B" ||
                    ele_help == "K") {element = ele_help; return;}
            }
        }
    
        ostringstream os;
        
        if (name == " UNK") {element = "UNK"; return;}
        if (res_name == "NAD" || res_name == "NDP" || res_name == "GPC") {
            if (name[2] == 'H' || name[2] == 'C' || name[2] == 'N' ||
                name[2] == 'O' || name[2] == 'S' || name[2] == 'P') {
                if (name[2] == 'P') {
                    if (name[1] == 'O') { //!Ausnahme bei NDP
                        element.assign(1,name[1]);
                        if (reset_names) {
                            os << name[1] << name[2] << name[3];
                            name = os.str();
                        }
                        return;
                    }
                }
                element.assign(1,name[2]);
                if (reset_names) {
                    os << (name[2]) << (name[3]);
                    name = os.str();
                }
            } else {
                element.assign(1,name[1]);
                if (reset_names) {
                    os << name[1] << name[2] << name[3];
                    name = os.str();
                }
            }
            return;
        }

        unsigned int ta = 0;
        
        if (name[0] == ' ' || name[0] == '1' || name[0] == '2' ||
            name[0] == '3' || name[0] == '4' || name[0] == '\'' ||
            name[0] == '"' || name[0] == '*') ta = 1;
        if (name[0+ta] ==  'H') {
            if (name.size() > 1+ta) {
                if (name[1+ta] ==  'G' || name[1+ta] ==  'g') {
                    if (ext.zero()) ext = new ATOM_EXT();
                    ext->pos_metal = true;
                }
                if (reset_names) {
                    if (name.size() > 2+ta) os << name[0+ta] << name[1+ta] << name[2+ta];
                    else os << name[0+ta] << name[1+ta];
                    name = os.str();
                }
            }
            element = "H"; 
            return;
        }
        if (name[0+ta] ==  'C') {
            if (name.size() > 1+ta) {
                if (name[1+ta] ==  'L' || name[1+ta] ==  'l') {
                    element = "Cl"; 
                    if (reset_names) {
                        if (name.size() > 2+ta) os << name[0+ta] 
                            << name[1+ta] << name[2+ta];
                        else os << name[0+ta] << name[1+ta];
                        name = os.str();
                    }
                    return;
                }
                if (name[1+ta] ==  'A' || name[1+ta] ==  'a' ||
                    name[1+ta] ==  'U' || name[1+ta] ==  'u') {
                        if (ext.zero()) ext = new ATOM_EXT();
                        ext->pos_metal = true;
                    }
                if (reset_names) {
                    if (name.size() > 2+ta) os << name[0+ta] << name[1+ta] << name[2+ta];
                    else os << name[0+ta] << name[1+ta];
                    name = os.str();
                }
            }
            element = "C";
            return;
        }
        if (name[0+ta] ==  'N') {
            if (name.size() > 1+ta) {
                if (name[1+ta] ==  'I' || name[1+ta] ==  'i') {
                    if (ext.zero()) ext = new ATOM_EXT();
                    ext->pos_metal = true;
                }
                if (name[1+ta] ==  'A' || name[1+ta] ==  'a') {
                    if (ext.zero()) ext = new ATOM_EXT();
                    ext->pos_metal = true;
                }
                if (reset_names) {
                    if (name.size() > 2+ta) os << name[0+ta] << name[1+ta] << name[2+ta];
                    else os << name[0+ta] << name[1+ta];
                    name = os.str();
                }
            }
            element = "N"; return;
        }
        if (name[0+ta] ==  'O') {
            element = "O"; 
            if (reset_names) {
                if (name.size() > 2+ta) os << name[0+ta] << name[1+ta] << name[2+ta];
                else os << name[0+ta] << name[1+ta];
                name = os.str();
            }
            return;
        }
        if (name[0+ta] ==  'P') {
            element = "P";
            if (reset_names) {
                if (name.size() > 2+ta) os << name[0+ta] << name[1+ta] << name[2+ta];
                else os << name[0+ta] << name[1+ta];
                name = os.str();
            }
            return;
        }
        if (name[0+ta] ==  'I') {
            element = "I"; 
            if (reset_names) {
                if (name.size() > 2+ta) os << name[0+ta] << name[1+ta] << name[2+ta];
                else os << name[0+ta] << name[1+ta];
                name = os.str();
            }
            return;
        }
        if (name[0+ta] ==  'B') {
            element = "B";
            if (name.size() > 1+ta) {
                if (name[1+ta] ==  'R' || name[1+ta] ==  'r') {
                    element = "Br";
                    if (reset_names) {
                        if (name.size() > 2+ta) os << name[0+ta] 
                            << name[1+ta] << name[2+ta];
                        else os << name[0+ta] << name[1+ta];
                        name = os.str();
                    }
                    return;
                }
                if (reset_names) {
                    if (name.size() > 2+ta) os << name[0+ta] << name[1+ta] << name[2+ta];
                    else os << name[0+ta] << name[1+ta];
                    name = os.str();
                }
            }
            return;
        }
        if (name[0+ta] ==  'S') {
            element = "S";
            if (name.size() > 1+ta) {
                if (name[1+ta] ==  'I' || name[1+ta] ==  'i') {
                    if (ext.zero()) ext = new ATOM_EXT();
                    ext->pos_metal = true;
                }
                if (name[1+ta] ==  'E' || name[1+ta] ==  'e') {
                    if (ext.zero()) ext = new ATOM_EXT();
                    ext->pos_metal = true;
                }
                if (reset_names) {
                    if (name.size() > 2+ta) os << name[0+ta] << name[1+ta] << name[2+ta];
                    else os << name[0+ta] << name[1+ta];
                    name = os.str();
                }
            }
            return;
        }
        if (name[0+ta] ==  'F') {
            element = "F";
            if (name.size() > 1+ta) {
                if (name[1+ta] ==  'E' || name[1+ta] ==  'e') {
                    if (ext.zero()) ext = new ATOM_EXT();
                    ext->pos_metal = true;
                }
                if (reset_names) {
                    if (name.size() > 2+ta) os << name[0+ta] << name[1+ta] << name[2+ta];
                    else os << name[0+ta] << name[1+ta];
                    name = os.str();
                }
            }
            return;
        }
        
        element.assign(1,name[0+ta]);
        if (reset_names) {
            if (name.size() > 2+ta) os << name[0+ta] << name[1+ta] << name[2+ta];
            else os << name[0+ta] << name[1+ta];
            name = os.str();
        }
        if (ext.zero()) ext = new ATOM_EXT();
        ext->pos_metal = true; //!wenn sonst nicht zuzuordnen ist ein metall recht wahrscheinlich
    }
}

bool ATOM::is_in_same_ring(stl_ptr<ATOM>& ref) {
    if (ext.zero() || ref->ext.zero()) {
        cerr << c_message<cERROR>("ATOM::is_in_same_ring --> no EXT object") << endl;
        exit(1);
    }
    if (ext->ring_ptrs.size() == 0 || ref->ext->ring_ptrs.size() == 0) return false;
    for (tr1::unordered_set<RING*>::iterator it=ext->ring_ptrs.begin(); it!=ext->ring_ptrs.end(); ++it) {
        if (ref->ext->ring_ptrs.find(*it) != ref->ext->ring_ptrs.end()) return true;
    }
    return false;
}

float const& ATOM::get_smallest_bond_angle() {
    if (ext->smallest_ba > 0. || bonded_atoms.size() < 2) return ext->smallest_ba;
    vec3d<float> v1 = bonded_atoms[0]->coord;
    v1 -= coord;
    vec3d<float> v2 = bonded_atoms[1]->coord;
    v2 -= coord;
    ext->smallest_ba = angle(v1,v2);
    if (bonded_atoms.size() > 2) {
        vec3d<float> v3 = bonded_atoms[2]->coord;
        v3 -= coord;
        if (angle(v1,v3) < ext->smallest_ba) ext->smallest_ba = angle(v1,v3);
        if (angle(v2,v3) < ext->smallest_ba) ext->smallest_ba = angle(v2,v3);
        if (bonded_atoms.size() > 3) {
            vec3d<float> v4 = bonded_atoms[3]->coord;
            v4 -= coord;
            if (angle(v1,v4) < ext->smallest_ba) ext->smallest_ba = angle(v1,v4);
            if (angle(v2,v4) < ext->smallest_ba) ext->smallest_ba = angle(v2,v4);
            if (angle(v3,v4) < ext->smallest_ba) ext->smallest_ba = angle(v3,v4);
        }
    }
    return ext->smallest_ba;
}

float const& ATOM::get_trigo_angle() {
    if (ext->trigo_angle > 0. || bonded_atoms.size() < 3) return ext->trigo_angle;
    vec3d<float> v1 = bonded_atoms[0]->coord; v1 -= coord;
    vec3d<float> v2 = bonded_atoms[1]->coord; v2 -= coord;
    vec3d<float> v3 = bonded_atoms[2]->coord; v3 -= coord;
    v1 *= v2;
    ext->trigo_angle = fabs(v1.skalar_product(v3));
    return ext->trigo_angle;
}

//==============================================================================================
//Definitionen fr ATOM_EXT:
//==============================================================================================

ATOM_EXT::ATOM_EXT() : compare_number(0), compare_number2(0), hybridization(0), n_heavy_bonded(0),
                       sp(0.), sp2(0.), sp3(0.), buriedness(-1.),smallest_ba (-1.),trigo_angle(-1.),
                       check_hyb(false), sec_check(false), is_aromatic(false), is_ring(false),
                       is_planar_ring(false), not_full_sp2_ring(false), pos_metal(false) {}

ATOM_EXT::ATOM_EXT(ATOM_EXT const& ext) : compare_number(ext.compare_number),compare_number2(ext.compare_number2),
                                    hybridization(ext.hybridization), n_heavy_bonded(ext.n_heavy_bonded),
                                    sp(ext.sp), sp2(ext.sp2), sp3(ext.sp3), buriedness(ext.buriedness),smallest_ba(ext.smallest_ba),
                                    trigo_angle(ext.trigo_angle),check_hyb(ext.check_hyb), sec_check(ext.sec_check), is_aromatic(ext.is_aromatic),
                                    is_ring(ext.is_ring), is_planar_ring(ext.is_planar_ring),
                                    not_full_sp2_ring(ext.not_full_sp2_ring), pos_metal(ext.pos_metal) {}

ATOM_EXT::~ATOM_EXT() {} //!rings wird von molecule gekillt

void ATOM_EXT::clear() {
    compare_number = 0;
    compare_number2 = 0;
    hybridization = 0; n_heavy_bonded = 0; pos_metal = false; check_hyb = false; sec_check = false;
    is_aromatic = false; is_ring = false; is_planar_ring = false; sp = 0.; sp2 = 0.; sp3 = 0.;
    buriedness = -1.;
    smallest_ba = -1.;
    trigo_angle = -1.;
    existing_bonds.clear();
    ring_ptrs.clear();
    //best_bonds.clear(); sec_bonds.clear();
    //rings.clear();
}


//==============================================================================================
//Definitionen fuer COMPARENODE:
//==============================================================================================

COMPARENODE::COMPARENODE(stl_ptr<ATOM> &atom,vector<stl_ptr<ATOM> > &prev) : root_atom(atom) {
    for (atoms_vec at=prev.begin(); at!=prev.end(); ++at) {
        prev_roots.push_back(*at);
    }
    if (atom->element != "H") prev_roots.push_back(atom);
}

COMPARENODE::~COMPARENODE() {
    for (compares_vec it=sub_nodes.begin(); it!=sub_nodes.end(); ++it) {
        it->kill();
    }
}

void COMPARENODE::make_subnodes() {
    //!Hier werden auch immer wieder neu subnodes angelegt
    //!Lieber vorher einmal fuer jedes ATOM-Object die Comparenodes (mit allen Subnodes) anlegen
    //!Dann werden zwar einige Subnodes umsonst angelegt, aber dafuer passiert das ganze nicht fuer ein Atom 
    //!mehrfach
    
    for (atoms_vec at=root_atom->bonded_atoms.begin(); at!=root_atom->bonded_atoms.end(); ++at) {
        if ((*at)->element == "H") continue;
        bool conti = false;
        for (atoms_vec bt=prev_roots.begin(); bt!=prev_roots.end(); ++bt) {
            if (*at == *bt) { //keinen neuen subnode machen, weil dieses atom schonmal root war
                                          //verhindert auch das endloslaufen in ringen
                conti = true;
                break;
            }
        }
        if (conti) continue;
        stl_ptr<COMPARENODE> cn(new COMPARENODE(*at,prev_roots));
        sub_nodes.push_back(cn);
    }
}

bool COMPARENODE::is_equal(COMPARENODE& cnode) {
    if (root_atom->ext->compare_number != cnode.root_atom->ext->compare_number) return false;
    make_subnodes();
    cnode.make_subnodes();
    for (compares_vec it=sub_nodes.begin(); it!=sub_nodes.end(); ++it) {
        bool has_eq = false;
        for (compares_vec jt=cnode.sub_nodes.begin(); jt!=cnode.sub_nodes.end();++jt) {
            if ((*it)->is_equal(*(*jt))) {
                has_eq = true;
                break;
            }
        }
    
        if (!has_eq) return false;
        continue;
    }
    return true;
}


//==============================================================================================
//Definitionen fr RING:
//==============================================================================================

RING::RING() : n_hetero(0), n_members(0), n_pi_ele(0), is_planar(false), is_aromatic(false), 
               is_pos(false), is_prot(false) {}

RING::~RING() {}

bool RING::fused(stl_ptr<RING> const& rg) {
    if (atom_set.size() == 0) {
        for (atoms_vec at=ring_atoms.begin(); at!=ring_atoms.end(); ++at) atom_set.insert((*at)->intern_id);
    }
    int n_shared = 0;
    for (const_atoms_vec at=rg->ring_atoms.begin(); at!=rg->ring_atoms.end(); ++at) {
        if (atom_set.find((*at)->intern_id) != atom_set.end()) ++n_shared;
        if (n_shared > 1) return true;
    }
    return false;
}
//==============================================================================================
//Definitionen fuer NODE:
//==============================================================================================

NODE::NODE(stl_ptr<ATOM> &at) : is_end(false), atm(at), prev(NULL)  {}

NODE::~NODE() {}

//==============================================================================================
//Definitionen fr BOND:
//==============================================================================================

BOND::BOND() : free_rot(false), dict_type(" "), type("un") {}

BOND::~BOND() {}

bool BOND::operator<(BOND const& rechts) {
    if (from < rechts.from) return true;
    if (from == rechts.from && to < rechts.to) return true;
    return false;
}


//==============================================================================================
//Definitionen fr PDBCONECT:
//==============================================================================================

PDBCONECT::PDBCONECT() {
    for (int i=0; i<4 ; ++i) {cov_bonds[i] = 0; h_bonds[i] = 0;}
    salt_bonds[0] = 0; salt_bonds[1] = 0;
}

PDBCONECT::~PDBCONECT() {}

//==============================================================================================
//Definitionen fr HET:
//==============================================================================================

HET::HET() {}

HET::~HET() {}


//==============================================================================================
//Definitionen fr COMMENT:
//==============================================================================================

COMMENT::COMMENT() {
    text = "";
}

COMMENT::~COMMENT() {}


//==============================================================================================
//Definitionen fr WATER:
//==============================================================================================

WATER::WATER(bool calc) : calculated(calc){}

WATER::~WATER() {
    atom.kill();
}


//==============================================================================================
//Definitionen fr METAL:
//==============================================================================================

METAL::METAL(bool cov) : is_covalent(cov) {}

METAL::~METAL() {
    atom.kill();
}


ostream &operator<<(ostream &os,ATOM const& ref) {
    os << "[ atomname  intern_id  element  intern_type  sybyl_type             coords              has_ext ]" << "\n";
    os << "[   "; os.width(10); os << left << ref.name; os.width(11); os << left << ref.intern_id; os.width(9); os << left << ref.element;
    os.width(13); os << left << ref.intern_type; os.width(10); os << left << ref.sybyl_type; os << ref.coord;
    if (!(ref.ext.zero())) {
        os << "    yes    ]" << "\n";
        os << "[  hybrid.   n_heavy   bonded_atoms   is_ring  planar_ring   is_aro.                            ]" << "\n";
        os << "[   "; os.width(12); os << left << ref.ext->hybridization; os.width(11); os << left << ref.ext->n_heavy_bonded;
        os.width(15); os << left << ref.bonded_atoms.size();
        os.width(9); os << left << ref.ext->is_ring; os.width(13); os << left << ref.ext->is_planar_ring;
        os.width(10); os << left << ref.ext->is_aromatic;
        os << "                      ]";
    } else os << "    no     ]";
    return os;
}
