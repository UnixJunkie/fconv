
//============================================================================
// files_GN.cpp -*- C++ -*-; reading and writing of molecular data files
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
// This library is for parsing and writing of molecular data files, such
// as PDB(QT), MOL2, DLG, CIF, SDF
//============================================================================


#include"files_GN.h"

//==============================================================================================
//Definitionen fr PARSER:
//==============================================================================================

PARSER::PARSER(STRUCTURE *st, int vb, bool parse_hyd, bool parse_wat, bool parse_ligs, bool parse_stuff) :
    main_structure(st), structure(st) {
    still_open = "XXX";
    if ((vb >= 0)&&(vb <= 2)) verbosity = vb;
    else { //das default-Level fr die verbosity ist 1
        cout << "possible verbosity-levels are 0, 1, or 2  --  setting to default (1)" << endl;
        verbosity = 1;
    }
    parse_hydrogens = parse_hyd; //Wasserstoffe mit parsen? (default = false) 
    parse_water = parse_wat; //Wasser mit parsen? (default = false)
    parse_ligands = parse_ligs;
    parse_comments = parse_stuff;
    n_conects = 0;
    pdbqt = false;
    atom_properties::initialize();
}

PARSER::~PARSER() {
    if (f_in) f_in.close();
    if (f_out) f_out.close();
}

void PARSER::clear() {
    structure = main_structure;
    still_open ="XXX";
    last_atm = NULL;
    h_log.clear();
    rh_log.clear();
    cav_log.clear();
    n_conects = 0;
}


void PARSER::set_verbosity(int vb) {
    //! level 0: nur error-Meldungen
    //! level 1: auch warnings
    //! level 2: auch status-informationen
    if ((vb >= 0)&&(vb <= 2)) verbosity = vb;
    else {
        cout << "possible verbosity-levels are 0, 1, or 2  --  setting to default (1)" << endl;
        verbosity = 1;
    }
}


//eine ATOM-Line einer pdb-Datei auswerten
void PARSER::parse_atom_line(string const& row) { //Die weitere Auswertung mu�das STRUCTURE-Obj. bernehmen
    //! wenn kein H gelesen wird muessen noch alle anderen id's angepasst werden
    //! => noch Variable anlegen, die die H's zaehlt und den Wert von folgenden id's abzieht
    //!    (auch im CONECT Teil nicht vergessen!)
    //! (erst beim Schreiben bercksichtigen?)
    //! gleiches muss eventuell auch noch fuer HOH implementiert werden!
    
    //!11.04.07:  Ketten nach einem HETATM trotzdem weiterparsen, obwohl dann auch Peptide mit geparsed werden
    //!if (hetatoms && !(parse_ligands)) return;
    
    if (!parse_hydrogens) { //kein Wasserstoff mit einlesen
        //!GN 03.05.06:  auch row[12]
        if ((row[13] == 'H') || (row[12] == 'H')) { //kann in diesem Fall nur Wasserstoff sein
            return;
        }
    }
    
    //!18.08.06:
    string s2(row,12,4);
    if ((s2 == " DUM") || (s2[0] == ' ' && s2[1] == 'Q')) return;
    
    //!alle STRUCTURE-Objekte (ausser protein) werden ueber stl_ptr verwaltet (siehe 'stl_ptr_GN.hpp')
    stl_ptr<ATOM> atm(new ATOM()); //neues ATOM-Objekt erzeugen
    atm->type = 1;
    atm->intern_id = cur_id; //fortlaufende id fuer interne Verwaltung
    cur_id += 1; //fortlaufende Nummer hochzaehlen
    atm->model_number = structure->model_number;
    //line in einzelne strings zerlegen
    string s1(row,6,5), s4(row,17,3), s6(row,22,4), //id / name / res_name / res_number
           s7(row,30,8), s8(row,38,8), s9(row,46,8); //Koordinaten
    
    //! 07.10.09: Fuer den Fall das die Kristallographen wieder vergessen haben HETATMs auch
    //! als solche zu kennzeichnen
    if (atom_properties::known_aacids.find(s4) == atom_properties::known_aacids.end()) {
        if (atom_properties::known_mod_aacids.find(s4) == atom_properties::known_mod_aacids.end()) {
            if (atom_properties::known_nucleic_acids.find(s4) == atom_properties::known_nucleic_acids.end()) {
                atm.kill();
                cur_id -= 1;
                parse_hetatm_line(row);
                return;
            }
        }
    }
    
    string s10, s11, s12, s13;
    
    try {
        if (pdbqt) {
            s10.assign(row,54,6); s11.assign(row,60,6); s12.assign(row,77,2); s13.assign(row,70,6); //occup. / b_factor / element / charge
        } else {
            s10.assign(row,54,6); s11.assign(row,60,6); s12.assign(row,76,2); s13.assign(row,78,2); //occup. / b_factor / element / charge
        }
    } catch(...) {};
                    
    is.clear();
    is.str(s1);
    
    is >> atm->id;    //Die "Seriennummer" des Atoms
    h_log[atm->id] = atm->intern_id; //um spaeter die BONDs auswerten zu koennen
    atm->name = s2;  //Atom-name  --  z.B.  CA  fuer C-alpha
    atm->alt_loc_id = row[16]; //alternative location fuer das Atom?
    atm->res_name = s4; //residue-name
    atm->chain_id = row[21]; //chain_id
    atm->alt_res = row[26]; //neu: fuer alternative Seitenketten
    is.clear(); is.str(s6); is >> atm->res_number; //integer
    
    if (atm->alt_res != ' ') { //!Neu: 17.12.2009
        if ((atm->res_number != last_res_add_num) || (atm->alt_res != last_alt_res_num)) {
            ++res_add;
            last_res_add_num = atm->res_number;
            last_alt_res_num = atm->alt_res;
        }
        atm->alt_res = ' ';
    }
    atm->res_number += res_add;  //!Neu: 17.12.2009
    
    
    //!TER einfuegen, wenn AS in der Kette fehlen (z.B. fuer ausgeschnittene Tasche)
    if (!(last_atm.zero())) {
        if (((last_atm->res_number - atm->res_number) > 1) ||
            ((atm->res_number - last_atm->res_number) > 1)) {
            if (last_atm->chain_id == atm->chain_id) {
                if (!(last_atm->is_ter)) {
                    last_atm->is_ter = true; //auf das zuletzt geparste Atom folgt ein TER-Entry
                    atm->intern_id = cur_id;
                    cur_id += 1;
                    h_log[atm->id] = atm->intern_id;
                }
            }
        }
    }
    //!Kommentar: Es gibt Probleme, wenn sich im pdb 2 chains (z.B. A und B) staendig abwechseln und keine TER-Eintraege
    //!           vorhanden sind und in der Kette AS fehlen
    
    is.clear(); is.str(s7); is >> atm->coord[0]; //x
    is.clear(); is.str(s8); is >> atm->coord[1]; //y
    is.clear(); is.str(s9); is >> atm->coord[2]; //z
    try { //siehe oben
        is.clear(); is.str(s10); 
        if (is.str() == "      ") atm->occupancy = 1.; 
        else {
            is >> ftmp; atm->occupancy = ftmp; //occupancy
        }
        is.clear(); is.str(s11);
        if (is.str() == "      ") atm->b_factor = 0.; 
        else {
            is >> ftmp; atm->b_factor = ftmp; //temp_factor
        }

        if (pdbqt) {
            if (s12 == "Cl" || s12 == "Br") atm->element = s12;
            else if (s12[0] == 'A' || s12[0] == 'C') atm->element = "C";
            else if (s12[0] == 'O') atm->element = "O";
            else if (s12[0] == 'N') atm->element = "N";
            else if (s12[0] == 'S') atm->element = "S";
            else if (s12[0] == 'P') atm->element = "P";
            else if (s12[0] == 'H') atm->element = "H";
            else if (s12[0] == 'F') atm->element = "F";
            else if (s12[0] == 'I') atm->element = "I";
            else {
                atm->element.assign(1,s12[0]);
                if (atm->element.size() > 1) {
                    atm->element = s12;
                    if (s12[1] < 97) atm->element[1] += 32;
                }
            }
        } else atm->element = s12; //!element -- stimmt oft nicht, obwohl gefordert

        if (s13.size() > 0 && s13 != "  ") atm->charge = s13; //!so gut wie nie vorhanden
    } catch(...) {};
                    
    ostringstream hos;
    hos << atm->res_name << atm->res_number << atm->alt_res;
    atm->full_res_name = hos.str();
    
    
    if (structure->protein == NULL) {
        structure->protein = new PROTEIN(); //neues PROTEIN-Obj.
            
        structure->protein->main_structure = main_structure;
    }
    
    if (last_atm.zero() || last_atm->type > 1 || last_atm->is_endmdl ||
        (atm->chain_id != last_atm->chain_id)) { //neue Kette ?
        stl_ptr<CHAIN> chn(new CHAIN());
        if (pdbqt) chn->change_ele_set(true); // Die Elementtypen fuer PDBQTs werden schon hier bestimmt
        structure->protein->chains.push_back(chn);
        chn->protein = structure->protein; //neu: jede Chain bekommt einen Zeiger auf ihr Protein
        if (hetatoms && structure->ligands.size() > 0) {
            structure->protein->chains.back()->prev_lig = structure->ligands.back();
            structure->ligands.back()->next_chain = structure->protein->chains.back();
            hetatoms = false; //nur die erste kette nach den hetatoms
            main_structure->has_peptides = true;
        }
    }
    
    structure->protein->chains.back()->atoms.push_back(atm); //Atom an die aktuelle Chain anhaengen
        
    last_atm = &(*atm);
}

//eine HETATM-Line einer pdb-Datei auswerten
void PARSER::parse_hetatm_line(string const& row) {
    hetatoms = true;
    
//    if (!parse_ligands) return; //! 07.10.09: modifizierte AS sind auch HETATM!!! => siehe unten
    
    if (!parse_water) { //kein Wasser mit einlesen
        if (row[13] == 'O') {
            if ((row[17] == 'H')&&(row[18] == 'O')&&(row[19] == 'H')) return;
            if ((row[17] == 'H')&&(row[18] == '2')&&(row[19] == 'O')) return;
        }
    }
    if (!parse_hydrogens) { //kein Wasserstoff mit einlesen
        if ((row[13] == 'H') || (row[12] == 'H')) {
            return;
        }
    }
    
    //!18.08.06:
    string s2(row,12,4);
    if (s2 == " DUM" || (s2[0] == ' ' && s2[1] == 'Q')) return;
    
    stl_ptr<ATOM> atm(new ATOM()); //neues ATOM-Objekt erzeugen
    atm->type = 2;
    atm->intern_id = cur_id; //fortlaufende id fr interne Verwaltung
    cur_id += 1; //fortlaufende Nummer hochzaehlen
    atm->model_number = structure->model_number;
    string s1(row,6,5), s4(row,17,3), s6(row,22,4), //id / name / res_name / res_number
           s7(row,30,8), s8(row,38,8), s9(row,46,8); //Koordinaten
    
    //! 07.10.09: weil modifizierte AS als HETATM abgelegt werden: =======
    if (!parse_ligands) {
        if (atom_properties::known_mod_aacids.find(s4) == atom_properties::known_mod_aacids.end()) {
            atm.kill();
            cur_id -= 1;
            return;
        } else {
            atm.kill();
            cur_id -= 1;
            parse_atom_line(row);
            return;
        }
    } else {
        if (atom_properties::known_mod_aacids.find(s4) != atom_properties::known_mod_aacids.end()) {
            atm.kill();
            cur_id -= 1;
            parse_atom_line(row);
            return;
        } else if (atom_properties::known_nucleic_acids.find(s4) != atom_properties::known_nucleic_acids.end()) {
            atm.kill();
            cur_id -= 1;
            parse_atom_line(row);
            return;
        }
    }
    //!===================================================================
    
    string s10, s11, s12, s13;
    try { //siehe ATOM
        if (pdbqt) {
            s10.assign(row,54,6); s11.assign(row,60,6); s12.assign(row,77,2); s13.assign(row,70,6); //occup. / b_factor / element / charge
        } else {
            s10.assign(row,54,6); s11.assign(row,60,6); s12.assign(row,76,2); s13.assign(row,78,2); //occup. / b_factor / element / charge
        }
    } catch(...) {};
    is.clear(); //streamobjekt zuruecksetzen
    is.str(s1); //string dem streamobjekt zuweisen
    is >> atm->id;     //Die "Seriennummer" des Atoms
    
    atm->name = s2;  //Atom-name  --  z.B.  CA  fr C-alpha
    atm->alt_loc_id = row[16]; //alternative location fr das Atom?
    atm->res_name = s4; //residue-name
    atm->chain_id = row[21]; //chain_id
    atm->alt_res = ' ';
    is.clear(); is.str(s6); is >> atm->res_number; //integer
    is.clear(); is.str(s7); is >> atm->coord[0]; //x
    is.clear(); is.str(s8); is >> atm->coord[1]; //y
    is.clear(); is.str(s9); is >> atm->coord[2]; //z
    try {
        is.clear(); is.str(s10); 
        if (is.str() == "      ") atm->occupancy = 1.; 
        else {
            is >> ftmp; atm->occupancy = ftmp; //occupancy
        }
        is.clear(); is.str(s11);
        if (is.str() == "      ") atm->b_factor = 0.; 
        else {
            is >> ftmp; atm->b_factor = ftmp; //temp_factor
        }

        if (pdbqt) {
            if (s12 == "Cl" || s12 == "Br") atm->element = s12;
            else if (s12[0] == 'A' || s12[0] == 'C') atm->element = "C";
            else if (s12[0] == 'O') atm->element = "O";
            else if (s12[0] == 'N') atm->element = "N";
            else if (s12[0] == 'S') atm->element = "S";
            else if (s12[0] == 'P') atm->element = "P";
            else if (s12[0] == 'H') atm->element = "H";
            else if (s12[0] == 'F') atm->element = "F";
            else if (s12[0] == 'I') atm->element = "I";
            else {
                atm->element.assign(1,s12[0]);
                if (atm->element.size() > 1) {
                    atm->element = s12;
                    if (s12[1] < 97) atm->element[1] += 32;
                }
            }
        } else atm->element = s12; //element
        
        if (s13.size() > 0 && s13 != "  ") atm->charge = s13;
    } catch(...) {};
    
    atm->full_res_name = " ";
    //!Als Ligand werden nur Molekuele mit HET-Entry erkannt
    //!NEIN - Das ist jetzt nicht mehr so!!!
    for (het_entrys_vec it=structure->het_entrys.begin(); it!=structure->het_entrys.end();++it) { //ueber alle HET-Obj.
        if (atm->res_name == (*it)->res_name) { //wenn das HETATM zu einem Liganden gehoert
            //am 15.08.07 rausgenommen, damit alle Metalle erstmal in ligands sind (und wieder rein :) )
            for (int i=0; i<n_metal_list; ++i) {
                if (atm->res_name == metal_list_ra[i]) { //handelt es sich um ein Metall, dass in elements.hpp bercksichtigt ist?
                    atm->type = 3;
                    stl_ptr<METAL> met(new METAL());
                    met->atom = atm;
                    structure->metals.push_back(met);
                    last_atm = atm;
                    h_log[atm->id] = atm->intern_id;
                    return;
                }
            }
            
            if ((structure->ligands.size() == 0) || (structure->ligands.back()->atoms.size() == 0) ||
                (atm->res_name != structure->ligands.back()->atoms.back()->res_name) ||
                (atm->chain_id != structure->ligands.back()->atoms.back()->chain_id) ||
                (atm->res_number != structure->ligands.back()->atoms.back()->res_number)) { //neuer Ligand
            //    //!16.08.06: siehe parse_atom_line:  wenn peptidverknuepfung, dann mit vorherigem verbinden:
            //    if (merge_ligands) {
            //        merge_ligands = false;
            //        structure->ligands.back()->name += "_" + atm->res_name;
            //    } else {
                    stl_ptr<LIGAND> lig(new LIGAND());
                    if (pdbqt) lig->change_ele_set(true); // Die Elementtypen fuer PDBQTs werden schon hier bestimmt
                    lig->main_structure = main_structure;
                    lig->name = atm->res_name + "_";
                    string_fu::add2s(lig->name,atm->res_number);
                    if ((structure->ligands.size() != 0) &&
                        (structure->ligands.back()->atoms.size() != 0) &&
                        (atm->res_name == structure->ligands.back()->atoms.back()->res_name)) {
                        if (atm->chain_id != structure->ligands.back()->atoms.back()->chain_id) {
                            structure->ligands.back()->name = structure->ligands.back()->name + "_" +
                                                              structure->ligands.back()->atoms.back()->chain_id;
                            lig->name = lig->name + "_" + atm->chain_id;
                        }
                    }
                    
                    structure->ligands.push_back(lig);
            //    }
            }
            structure->ligands.back()->atoms.push_back(atm);
            last_atm = atm;
            h_log[atm->id] = atm->intern_id;
            return;
        }
    }
    
    // jetzt kommt eigentlich nur noch kleines Ion oder Wasser in frage
    if (row[13] == 'O' || row[13] == 'H') {
        if (((row[17] == 'H')&&(row[18] == 'O')&&(row[19] == 'H'))||
            ((row[17] == 'H')&&(row[18] == '2')&&(row[19] == 'O'))) {
            atm->type = 4;
            stl_ptr<WATER> wat(new WATER());
            wat->atom = atm;
            structure->waters.push_back(wat);
            last_atm = atm;
            h_log[atm->id] = atm->intern_id;
            return;
        }
    }
    
    if (verbosity > 1) {
        cerr << "note :   no HET-Entry for HETATM line:" << endl << row << endl << " in file: " << filename << endl;
        cerr << "         regardless it's taken as ligand atom in this version" << endl;
        cerr << "         (if it is a metal this can be corrected automatically later by get_atom_typing(bool mode) )" << endl << endl;
    }
    
    if ((structure->ligands.size() == 0) || (structure->ligands.back()->atoms.size() == 0) ||
        (atm->res_name != structure->ligands.back()->atoms.back()->res_name) ||
        (atm->chain_id != structure->ligands.back()->atoms.back()->chain_id) ||
        (atm->res_number != structure->ligands.back()->atoms.back()->res_number)) { //neuer Ligand
    //    //!16.08.06: siehe parse_atom_line:  wenn peptidverknuepfung, dann mit vorherigem verbinden:
    //    if (merge_ligands) {
    //        merge_ligands = false;
    //        structure->ligands.back()->name += "_" + atm->res_name;
    //    } else {
            stl_ptr<LIGAND> lig(new LIGAND());
            if (pdbqt) lig->change_ele_set(true); // Die Elementtypen fuer PDBQTs werden schon hier bestimmt
            lig->main_structure = main_structure;
            lig->name = atm->res_name + "_";
            string_fu::add2s(lig->name,atm->res_number);
            if ((structure->ligands.size() != 0) &&
                (structure->ligands.back()->atoms.size() != 0) &&
                (atm->res_name == structure->ligands.back()->atoms.back()->res_name)) {
                if (atm->chain_id != structure->ligands.back()->atoms.back()->chain_id) {
                    structure->ligands.back()->name = structure->ligands.back()->name + "_" +
                                        structure->ligands.back()->atoms.back()->chain_id;
                    lig->name = lig->name + "_" + atm->chain_id;
                }
            }
            
            structure->ligands.push_back(lig);
    //    }
    }
    
    structure->ligands.back()->atoms.push_back(atm);
    last_atm = atm;
    h_log[atm->id] = atm->intern_id;
    return;

//alte version, die keine HETATMs ohne HET-Entry zulies:
/*    //HETATM konnte nicht vernnftig ausgewertet werden
    if (verbosity > 0) cerr << "warning: could not process HETATM line:" << endl << row << endl <<
                               " in file: " << filename << endl;
    atm.kill(); //das ATOM-Obj. wieder loeschen, wenn es nicht zugeordnet werden konnte!
    cur_id -= 1;
*/
}

//eine HET-Line einer pdb-Datei auswerten
void PARSER::parse_het_line(string const& row) {
    stl_ptr<HET> het(new HET());
    string s1(row,7,3); //res_name(het_id) und Anzahl der korrespondierenden HETATMs (inkl. H's)
    het->res_name = s1;
    het->chain_id = row[12];
    main_structure->het_entrys.push_back(het);
}

//eine TER-Line einer pdb-Datei auswerten
void PARSER::parse_ter_line(string const& row) {
    last_atm->is_ter = true; //auf das zuletzt geparste Atom folgt ein TER-Entry
    cur_id += 1;
}

//eine MODEL-Line einer pdb-Datei verarbeiten (z.B. bei NMR-Strukturen)
void PARSER::parse_model_line(string const& row) {
    int mn;
    if (!main_structure->alt_struct.empty()) {
        mn = main_structure->alt_struct.back()->model_number;
    } else mn = main_structure->model_number;
    if (&(*last_atm) != NULL) { //nicht erster MODEL-Eintrag
        last_atm->is_model = true;
        stl_ptr<STRUCTURE> new_struct(new STRUCTURE());
        main_structure->alt_struct.push_back(new_struct);
    //    structure = &(*new_struct); //aktuelle STRUCTURE auf new_struct setzen
        structure = new_struct;
    }
    structure->model_number = mn+1; //!sind keine MODEL-Entrys vorhanden hat die STRUCTURE die
                                    //!model_number 0  (sonst 1..n)
    cur_id = 1; //Nummerierung von 1 an fuer jedes MODEL erzwingen (auch wenn im Widerspruch zu pdb-Standard)
}

//eine ENDMDL-Line einer pdb-Datei verarbeiten
void PARSER::parse_endmdl_line(string const& row) {
    last_atm->is_endmdl = true;
    structure = main_structure; //wieder die Hauptstruktur zur aktuellen Struktur machen
}

//Verknuepfungen der pdb-Atome auswerten
void PARSER::parse_conect_line(string const& row) {
    //!noch keine BOND-Objekte erzeugen (erst wenn diese gebraucht werden)
    stl_ptr<PDBCONECT> con(new PDBCONECT());
    string from(row,6,5); //from_atom_id
    
    string c1, c2, c3, c4, h1, h2, h3, h4, s1, s2;
    try {//siehe ATOM und HETATM
        c1.assign(row,11,5); c2.assign(row,16,5); c3.assign(row,21,5); c4.assign(row,26,5); //covalent bonds
        h1.assign(row,31,5); h2.assign(row,36,5); h3.assign(row,46,5); h4.assign(row,51,5); //h-bonds
        s1.assign(row,41,5); s2.assign(row,56,5);                             //salt-bridges
    } catch(...) {}; //!Achtung: nicht geklaert, ob dann unten Leerzeichen zugewiesen werden!!!
           
    int help; //erst in help lesen und dann pruefen, ob ein H beteiligt ist
    is.clear(); is.str(from); is >> help;
    
    con->from_id = help;
    is.clear(); is.str(c1);help = 0; is >> help; con->cov_bonds[0] = help;
    is.clear(); is.str(c2);help = 0; is >> help; con->cov_bonds[1] = help;
    is.clear(); is.str(c3);help = 0; is >> help; con->cov_bonds[2] = help;
    is.clear(); is.str(c4);help = 0; is >> help; con->cov_bonds[3] = help;
    
    is.clear(); is.str(h1);help = 0; is >> help; con->h_bonds[0] = help;
    is.clear(); is.str(h2);help = 0; is >> help; con->h_bonds[1] = help;
    is.clear(); is.str(h3);help = 0; is >> help; con->h_bonds[2] = help;
    is.clear(); is.str(h4);help = 0; is >> help; con->h_bonds[3] = help;
    
    is.clear(); is.str(s1);help = 0; is >> help; con->salt_bonds[0] = help;
    is.clear(); is.str(s2);help = 0; is >> help; con->salt_bonds[1] = help;
    
    main_structure->conects.push_back(con);
    
    //!CONECT-Eintraege fuer NMR-Strukturen sollten nur einmal (nach dem letzten MODEL) auftauchen!
    
}

//pdb-Line parsen
void PARSER::compute_pdb_line(string const& row) {
    if (row.compare(0,2,"AT") == 0) parse_atom_line(row);
    else if (row.compare(0,4,"HETA") == 0) parse_hetatm_line(row);
    else if (row.compare(0,3,"CON") == 0 && parse_ligands) parse_conect_line(row);
    else if (row.compare(0,4,"HET ") == 0 && parse_ligands) parse_het_line(row);
    else if (row.compare(0,4,"MODE") == 0) parse_model_line(row);
    else if (row.compare(0,4,"ENDM") == 0) parse_endmdl_line(row);
    else if (row.compare(0,3,"TER") == 0) parse_ter_line(row); 
    else if(parse_comments && ((row.compare(0,3,"HEA") == 0) || //auf unwichtigen Kram pruefen
            (row.compare(0,3,"COM") == 0) || //und diesen in string-vector ablegen
                (row.compare(0,3,"AUT") == 0) ||
            (row.compare(0,4,"TITL") == 0) ||
            (row.compare(0,4,"JRNL") == 0) ||
            (row.compare(0,4,"REMA") == 0) ||
            (row.compare(0,4,"KEYW") == 0) ||
            (row.compare(0,4,"REVD") == 0) ||
            (row.compare(0,4,"EXPD") == 0) ||
            (row.compare(0,4,"SOUR") == 0) ||
            (row.compare(0,4,"DBRE") == 0) ||
            (row.compare(0,4,"SEQR") == 0) ||
            (row.compare(0,4,"HELI") == 0) ||
            (row.compare(0,4,"FORM") == 0) ||
            (row.compare(0,4,"HETN") == 0) ||
            (row.compare(0,4,"HETS") == 0) ||
            (row.compare(0,4,"SHEE") == 0) ||
            (row.compare(0,4,"SITE") == 0) ||
            (row.compare(0,4,"CRYS") == 0) ||
            (row.compare(0,4,"ORIG") == 0) ||
            (row.compare(0,4,"SCAL") == 0) ||
            (row.compare(0,4,"REMA") == 0))) pdb_headers.push_back(string(row));
//    else if (row.compare(0,6,"SSBOND") == 0) parse_ssbond_line(row);
//    else if (row.compare(0,4,"END ") == 0) parse_end_line(row);
//    else if (row.compare(0,3,"ANISOU") == 0) parse_anisou_line(row);
}

//Haupt-Methode, die zum Lesen einer pdb-Datei aufgerufen wird
bool PARSER::read_pdb(char const* file,bool pdbqt_mode) {
    if (pdbqt_mode) pdbqt = true;
    else pdbqt = false;
    last_res_name = "@@@";
    res_add = 0;
    last_res_add_num = -1;
    last_alt_res_num = '#';
    hetatoms = false;
    h_log.clear();
    pdb_headers.clear();
    cur_id = 1;
    filename = file;
    if (verbosity > 1) cout << "reading " << filename << " ..." << endl;
    
    f_in.open(filename);
    
    if (!f_in) {cerr << c_message<cERROR>("PARSER::read_pdb --> could not open ") << filename << " !" << endl; f_in.clear(); return 0;}
    
    vector<string> readbuf;
    while (!f_in.eof()) {
        string row; //aktuelle Zeile
        getline(f_in,row); //Zeile einlesen
        readbuf.push_back(row);
    }
    
    merge_ligands = false;
    
    for (vector<string>::iterator it=readbuf.begin(); it!=readbuf.end(); ++it) {
        compute_pdb_line(*it);
    }
    
    readbuf.clear();
    
    f_in.close();
    f_in.clear();
    
//    if (merge_ligands) main_structure->merge_ligands();
    
    if (verbosity > 1) cout << "finished reading " << filename << endl;

    pdbqt = false;

    return 1;
}

void PARSER::write_header_lines() {
    for (vector<string>::iterator it=pdb_headers.begin(); it!=pdb_headers.end();++it) {
        f_out << *it << "\n"; //!kein endl, damit der Puffer nicht geleert wird!!!
    }
}

void PARSER::write_het_lines() {
    int n_hetatoms = 0;
    if (main_structure->het_entrys.size() > 0) { //GN: 04.07.06
        for (het_entrys_vec it=main_structure->het_entrys.begin(); it!=main_structure->het_entrys.end();++it) {
            f_out << "HET    ";
            f_out << (*it)->res_name; //res_name schreiben
            f_out << "  " << (*it)->chain_id << "       ";
            for (ligands_vec iat=main_structure->ligands.begin(); iat!=main_structure->ligands.end();++iat) {
                if ((*iat)->name == (*it)->res_name) n_hetatoms = (*iat)->atoms.size();
            }
            f_out.width(5);
            f_out << right << n_hetatoms << "\n";
        }
    }
}

void PARSER::write_first_model_line() { //Wenn mehrere Models, dann die erste MODEL-Line vor dem ersten ATOM schreiben
    if (!(main_structure->protein)) return;
    if (main_structure->protein->chains.size() == 0) return;
    if (main_structure->protein->chains[0]->atoms.size() == 0) return;
    if (main_structure->protein->chains[0]->atoms.front()->model_number > 0) {
        f_out << "MODEL        1" << "\n";
    }
}

void PARSER::write_pdb_atom(stl_ptr<ATOM> const& at,bool cav_mode) {
    if (at->type == 1) {
        //! 07.10.09: modifizierte AS muessen HETATM sein!!!
        if (atom_properties::known_mod_aacids.find(at->res_name) == atom_properties::known_mod_aacids.end()) f_out << "ATOM  ";
        else f_out << "HETATM";
    } else f_out << "HETATM";

    f_out.width(5);
    f_out << right << curr_id;
    
    rh_log[at->intern_id] = curr_id;

    f_out << " " << at->name; //name schreiben

    f_out << at->alt_loc_id; //alt_loc_id schreiben

    f_out << at->res_name; //res_name schreiben

    f_out << " " << at->chain_id; //chain_id schreiben

    f_out.width(4);
    f_out << right << at->res_number; //res_number schreiben
    
    f_out << at->alt_res << "   "; //neu: alt_res
    
    f_out.width(8); f_out.setf(ios::fixed); f_out.precision(3); //Koordinaten schreiben
    f_out << right << at->coord[0];
    f_out.width(8); f_out.precision(3);
    f_out << right << at->coord[1];
    f_out.width(8); f_out.precision(3);
    f_out << right << at->coord[2];
    
    f_out.width(6); f_out.precision(2); //Occupancy schreiben
    f_out << right << at->occupancy;
    
    f_out.width(6); f_out.precision(2); //temp_factor schreiben
    f_out << right << at->b_factor << "          ";
    
    f_out.unsetf(ios::fixed);

    f_out << at->element; //element schreiben
    
    f_out << at->charge; //Ladung schreiben
    
    f_out << "\n"; //Zeilenende
    
    if (at->is_ter) { //nach dem aktuellen ATOM folgt ein TER-Entry
        f_out << "TER   ";
        f_out.width(5);
        ++curr_id;
        f_out << right << curr_id << "      "; //id schreiben
        f_out << at->res_name << " "; //res_name schreiben
        f_out << at->chain_id; //chain_id schreiben
        f_out.width(4);
        f_out << right << at->res_number; //res_number schreiben
        f_out << "\n";
    }

    if (at->is_endmdl) { //nach dem aktuellen ATOM folgt ein ENDMDL-Entry
        f_out << "ENDMDL" << "\n";
    }
    
    if (at->is_model) { //nach dem aktuellen ATOM folgt ein MODEL-Entry
        f_out << "MODEL     ";
        f_out.width(4);
        f_out << right << (at->model_number)+1;
        f_out << "\n";
    }

    ++curr_id;
    ++total_written_atoms;
    if (cav_mode) {
        ++n_cav_atoms;
    //    cav_log[at->intern_id] = curr_id;
        last_atm = at;
    }
}

void PARSER::write_atoms() {
    rh_log.clear();
    curr_id = 1;
    total_written_atoms = 0;

    if (main_structure->protein) { //! auch noch fuer die alt_structs
        for (chains_vec cvt=main_structure->protein->chains.begin(); cvt!=main_structure->protein->chains.end(); ++cvt) {
            if (!((*cvt)->prev_lig.zero())) for (atoms_vec at=(*cvt)->prev_lig->atoms.begin(); at!=(*cvt)->prev_lig->atoms.end(); ++at) {
                write_pdb_atom(*at);
            }
            for (atoms_vec at=(*cvt)->atoms.begin(); at!=(*cvt)->atoms.end(); ++at) {
                write_pdb_atom(*at);
            }
        }
    }

    for (ligands_vec lig=main_structure->ligands.begin(); lig!=main_structure->ligands.end(); ++lig) {
        if (!(*lig)->next_chain.zero()) continue; //die wurden schon geschrieben
        for (atoms_vec at=(*lig)->atoms.begin(); at!=(*lig)->atoms.end(); ++at) {
            write_pdb_atom(*at);
        }
    }

    for (metals_vec met=main_structure->metals.begin(); met!=main_structure->metals.end(); ++met) {
        write_pdb_atom((*met)->atom);
    }
    for (waters_vec wat=main_structure->waters.begin(); wat!=main_structure->waters.end(); ++wat) {
        write_pdb_atom((*wat)->atom);
    }

    for (alt_struct_vec ast=main_structure->alt_struct.begin(); ast!=main_structure->alt_struct.end(); ++ast) {
        curr_id = 1;
        if ((*ast)->protein) { //! auch noch fuer die alt_structs
            for (chains_vec cvt=(*ast)->protein->chains.begin(); cvt!=(*ast)->protein->chains.end(); ++cvt) {
                if (!((*cvt)->prev_lig.zero())) for (atoms_vec at=(*cvt)->prev_lig->atoms.begin(); at!=(*cvt)->prev_lig->atoms.end(); ++at) {
                    write_pdb_atom(*at);
                }
                for (atoms_vec at=(*cvt)->atoms.begin(); at!=(*cvt)->atoms.end(); ++at) {
                    write_pdb_atom(*at);
                }
            }
        }
        for (ligands_vec lig=(*ast)->ligands.begin(); lig!=(*ast)->ligands.end(); ++lig) {
            if (!(*lig)->next_chain.zero()) continue; //die wurden schon geschrieben
            for (atoms_vec at=(*lig)->atoms.begin(); at!=(*lig)->atoms.end(); ++at) {
                write_pdb_atom(*at);
            }
        }
        for (metals_vec met=(*ast)->metals.begin(); met!=(*ast)->metals.end(); ++met) {
            write_pdb_atom((*met)->atom);
        }
        for (waters_vec wat=(*ast)->waters.begin(); wat!=(*ast)->waters.end(); ++wat) {
            write_pdb_atom((*wat)->atom);
        }
    }
}

void PARSER::write_cav_atoms(stl_ptr<CAVITY> const& cav) {
    rh_log.clear();
    curr_id = 1;
    total_written_atoms = 0;
    last_atm = 0;
    for (aacids_vec at=cav->aacids.begin(); at!=cav->aacids.end(); ++at) {
        for (atoms_vec it=(*at)->atoms.begin(); it!=(*at)->atoms.end(); ++it) {
            if (!(last_atm.zero())) {
                if (!(last_atm->is_ter)) {
                    if (((last_atm->res_number - (*it)->res_number) > 1) ||
                    (((*it)->res_number - last_atm->res_number) > 1)) {
                        f_out << "TER   ";
                        f_out.width(5);
                        f_out << right << curr_id << "      "; //id schreiben
                        f_out << last_atm->res_name << " "; //res_name schreiben
                        f_out << last_atm->chain_id; //chain_id schreiben
                        f_out.width(4);
                        f_out << right << last_atm->res_number; //res_number schreiben
                        f_out << "\n";
                        curr_id++;
                    }
                }
            }
            write_pdb_atom(*it,true);
        }
    }
}

void PARSER::write_atoms_special() {
    rh_log.clear();
    curr_id = 1;
    total_written_atoms = 0;
    int last_resn = -1;
    stl_ptr<ATOM> latm;
    if (main_structure->protein) { //! auch noch fuer die alt_structs
        for (chains_vec cvt=main_structure->protein->chains.begin(); cvt!=main_structure->protein->chains.end(); ++cvt) {
            if (!((*cvt)->prev_lig.zero())) for (atoms_vec at=(*cvt)->prev_lig->atoms.begin(); at!=(*cvt)->prev_lig->atoms.end(); ++at) {
                if ((*at)->bond_ind == 0) continue;
                if ((last_resn != -1) && ((*at)->res_number - last_resn) > 1) latm->is_ter = true;
                last_resn = (*at)->res_number;
                latm = *at;
            }
            for (atoms_vec at=(*cvt)->atoms.begin(); at!=(*cvt)->atoms.end(); ++at) {
                if ((*at)->bond_ind == 0) continue;
                if ((last_resn != -1) && ((*at)->res_number - last_resn) > 1) latm->is_ter = true;
                last_resn = (*at)->res_number;
                latm = *at;
            }
        }
    }
    
    for (ligands_vec lig=main_structure->ligands.begin(); lig!=main_structure->ligands.end(); ++lig) {
        if (!(*lig)->next_chain.zero()) continue; //die wurden schon geschrieben
        for (atoms_vec at=(*lig)->atoms.begin(); at!=(*lig)->atoms.end(); ++at) {
            if ((*at)->bond_ind == 0) continue;
            if ((last_resn != -1) && ((*at)->res_number - last_resn) > 1) latm->is_ter = true;
            last_resn = (*at)->res_number;
            latm = *at;
        }
    }
    
    for (alt_struct_vec ast=main_structure->alt_struct.begin(); ast!=main_structure->alt_struct.end(); ++ast) {
        curr_id = 1;
        if ((*ast)->protein) { //! auch noch fuer die alt_structs
            for (chains_vec cvt=(*ast)->protein->chains.begin(); cvt!=(*ast)->protein->chains.end(); ++cvt) {
                if (!((*cvt)->prev_lig.zero())) for (atoms_vec at=(*cvt)->prev_lig->atoms.begin(); at!=(*cvt)->prev_lig->atoms.end(); ++at) {
                    if ((*at)->bond_ind == 0) continue;
                    if ((last_resn != -1) && ((*at)->res_number - last_resn) > 1) latm->is_ter = true;
                    last_resn = (*at)->res_number;
                    latm = *at;
                }
                for (atoms_vec at=(*cvt)->atoms.begin(); at!=(*cvt)->atoms.end(); ++at) {
                    if ((*at)->bond_ind == 0) continue;
                    if ((last_resn != -1) && ((*at)->res_number - last_resn) > 1) latm->is_ter = true;
                    last_resn = (*at)->res_number;
                    latm = *at;
                }
            }
        }
        for (ligands_vec lig=(*ast)->ligands.begin(); lig!=(*ast)->ligands.end(); ++lig) {
            if (!(*lig)->next_chain.zero()) continue; //die wurden schon geschrieben
            for (atoms_vec at=(*lig)->atoms.begin(); at!=(*lig)->atoms.end(); ++at) {
                if ((*at)->bond_ind == 0) continue;
                if ((last_resn != -1) && ((*at)->res_number - last_resn) > 1) latm->is_ter = true;
                last_resn = (*at)->res_number;
                latm = *at;
            }
        }
    }
    
    if (!(latm.zero())) latm->is_ter = true;
    
    if (main_structure->protein) { //! auch noch fuer die alt_structs
        for (chains_vec cvt=main_structure->protein->chains.begin(); cvt!=main_structure->protein->chains.end(); ++cvt) {
            if (!((*cvt)->prev_lig.zero())) for (atoms_vec at=(*cvt)->prev_lig->atoms.begin(); at!=(*cvt)->prev_lig->atoms.end(); ++at) {
                if ((*at)->bond_ind == 0) continue;
                write_pdb_atom(*at);
            }
            for (atoms_vec at=(*cvt)->atoms.begin(); at!=(*cvt)->atoms.end(); ++at) {
                if ((*at)->bond_ind == 0) continue;
                write_pdb_atom(*at);
            }
        }
    }
    
    for (ligands_vec lig=main_structure->ligands.begin(); lig!=main_structure->ligands.end(); ++lig) {
        if (!(*lig)->next_chain.zero()) continue; //die wurden schon geschrieben
        for (atoms_vec at=(*lig)->atoms.begin(); at!=(*lig)->atoms.end(); ++at) {
            if ((*at)->bond_ind == 0) continue;
            write_pdb_atom(*at);
        }
    }
    for (metals_vec met=main_structure->metals.begin(); met!=main_structure->metals.end(); ++met) {
        if ((*met)->atom->bond_ind == 0) continue;
        write_pdb_atom((*met)->atom);
    }
    for (waters_vec wat=main_structure->waters.begin(); wat!=main_structure->waters.end(); ++wat) {
        if ((*wat)->atom->bond_ind == 0) continue;
        write_pdb_atom((*wat)->atom);
    }
    for (alt_struct_vec ast=main_structure->alt_struct.begin(); ast!=main_structure->alt_struct.end(); ++ast) {
        curr_id = 1;
        if ((*ast)->protein) { //! auch noch fuer die alt_structs
            for (chains_vec cvt=(*ast)->protein->chains.begin(); cvt!=(*ast)->protein->chains.end(); ++cvt) {
                if (!((*cvt)->prev_lig.zero())) for (atoms_vec at=(*cvt)->prev_lig->atoms.begin(); at!=(*cvt)->prev_lig->atoms.end(); ++at) {
                    if ((*at)->bond_ind == 0) continue;
                    write_pdb_atom(*at);
                }
                for (atoms_vec at=(*cvt)->atoms.begin(); at!=(*cvt)->atoms.end(); ++at) {
                    if ((*at)->bond_ind == 0) continue;
                    write_pdb_atom(*at);
                }
            }
        }
        for (ligands_vec lig=(*ast)->ligands.begin(); lig!=(*ast)->ligands.end(); ++lig) {
            if (!(*lig)->next_chain.zero()) continue; //die wurden schon geschrieben
            for (atoms_vec at=(*lig)->atoms.begin(); at!=(*lig)->atoms.end(); ++at) {
                if ((*at)->bond_ind == 0) continue;
                write_pdb_atom(*at);
            }
        }
        for (metals_vec met=(*ast)->metals.begin(); met!=(*ast)->metals.end(); ++met) {
            if ((*met)->atom->bond_ind == 0) continue;
            write_pdb_atom((*met)->atom);
        }
        for (waters_vec wat=(*ast)->waters.begin(); wat!=(*ast)->waters.end(); ++wat) {
            if ((*wat)->atom->bond_ind == 0) continue;
            write_pdb_atom((*wat)->atom);
        }
    }
}

void PARSER::write_conect_lines() {
    if (main_structure->conects.size() > 0) { //GN: 04.07.06
        n_conects = 0;
        for (conects_vec it=main_structure->conects.begin(); it!=main_structure->conects.end();++it) {
            //gucken, ob es die Atome zu diesem conect eintrag ueberhaupt gibt:
            if (h_log.find((*it)->from_id) != h_log.end()) {
                if (rh_log.find(h_log[(*it)->from_id]) == rh_log.end()) continue;
                
                bool hat_partner = false;
                for (int i=0; i<4; ++i) {
                    if (h_log.find((*it)->h_bonds[i]) != h_log.end()) {
                        if (rh_log.find(h_log[(*it)->h_bonds[i]]) != rh_log.end()
                            && (*it)->h_bonds[i] != 0) {
                            hat_partner = true;
                            break;
                        }
                    }
                    if (h_log.find((*it)->cov_bonds[i]) != h_log.end()) {
                        if (rh_log.find(h_log[(*it)->cov_bonds[i]]) != rh_log.end()
                            && (*it)->cov_bonds[i] != 0) {
                            hat_partner = true;
                            break;
                        }
                    }
                }
                if (!hat_partner) for (int i=0; i<4; ++i) {
                    if (h_log.find((*it)->salt_bonds[i]) != h_log.end()) {
                        if (rh_log.find(h_log[(*it)->salt_bonds[i]]) != rh_log.end()
                            && (*it)->salt_bonds[i] != 0) {
                            hat_partner = true;
                            break;
                        }
                    }
                }
                if (!hat_partner) continue; //diese connection gibt es nicht mehr
            } else continue;
            
            f_out << "CONECT";
            f_out.width(5);
            f_out << right << rh_log[h_log[(*it)->from_id]];
            
            for (int i=0; i<4; ++i) {
                if ((*it)->cov_bonds[i] == 0) f_out << "     ";
                else {f_out.width(5); f_out << right << rh_log[h_log[(*it)->cov_bonds[i]]];}
            }
            
            if ((*it)->h_bonds[0] == 0) f_out << "     ";
            else {f_out.width(5); f_out << right << rh_log[h_log[(*it)->h_bonds[0]]];}
            if ((*it)->h_bonds[1] == 0) f_out << "     ";
            else {f_out.width(5); f_out << right << rh_log[h_log[(*it)->h_bonds[1]]];}
            
            if ((*it)->salt_bonds[0] == 0) f_out << "     ";
            else {f_out.width(5); f_out << right << rh_log[h_log[(*it)->salt_bonds[0]]];}
            
            if ((*it)->h_bonds[2] == 0) f_out << "     ";
            else {f_out.width(5); f_out << right << rh_log[h_log[(*it)->h_bonds[2]]];}
            if ((*it)->h_bonds[3] == 0) f_out << "     ";
            else {f_out.width(5); f_out << right << rh_log[h_log[(*it)->h_bonds[3]]];}
            
            if ((*it)->salt_bonds[1] == 0) f_out << "     ";
            else {f_out.width(5); f_out << right << rh_log[h_log[(*it)->salt_bonds[1]]];}
            
            ++n_conects;
            
            f_out << "\n";
        }
    }
}

void PARSER::write_cav_conect_lines() {
    write_conect_lines();
    /*
    for (conects_vec it=main_structure->conects.begin(); it!=main_structure->conects.end();++it) {
        bool dflag = false;
        if (!(cav_log.find((*it)->from_id) == cav_log.end())) {
            for (int i=0; i<4; i++) {
                if (!(cav_log.find((*it)->cov_bonds[i]) == cav_log.end())) dflag = true;  
            }
        }
        if (dflag) {
            f_out << "CONECT";
            n_cav_conects++;
            f_out.width(5);
            f_out << right << cav_log[(*it)->from_id];
            for (int i=0; i<4; i++) {
                if (!(cav_log.find((*it)->cov_bonds[i]) == cav_log.end())) {
                    f_out.width(5); f_out << right << cav_log[(*it)->cov_bonds[i]];
                }
            }
            f_out << "\n";
        }
    }
    */
}

void PARSER::write_master() {
    f_out << "MASTER              ";
    f_out.width(5); f_out << right << main_structure->het_entrys.size() << "                         ";
    f_out.width(5); f_out << right << total_written_atoms << "     ";
    f_out.width(5); f_out << right << n_conects << "\n";
}

void PARSER::write_cav_master() {
    f_out << "MASTER              ";
    f_out.width(5); f_out << right << 0 << "                         ";
    f_out.width(5); f_out << right << n_cav_atoms << "     ";
    f_out.width(5); f_out << right << /*n_cav_conects*/ n_conects << "\n";
}

bool PARSER::write_pdb(char const* file) {
    filename = file;
    if (verbosity > 1) cout << "writing " << filename << " ..." << endl;
    f_out.open(filename);
    if (!f_out) {cerr << c_message<cERROR>("PARSER::write_pdb --> could not open ") << filename << " for writing!" << endl; f_out.clear(); return 0;}

    write_header_lines();
    f_out << "REMARK    This file was generated by PARSER (from files_GN.cpp)" << "\n";
    write_het_lines();
//    write_ssbond_lines();
    write_first_model_line();

    write_atoms();

    write_conect_lines();
    
    write_master();
    f_out << "END" << endl;
    f_out.close();
    f_out.clear();
    if (verbosity > 1) cout << "finished writing " << filename << endl;
    
    return 1;
}

bool PARSER::write_pdb_cav(char const* file,stl_ptr<CAVITY> const& cav) {
    filename = file;
    if (verbosity > 1) cout << "writing " << filename << " ..." << endl;
    f_out.open(filename);
    if (!f_out) {cerr << c_message<cERROR>("PARSER::write_pdb_cav --> could not open ") << filename << " for writing!" << endl; f_out.clear(); return 0;}
    f_out << "REMARK    This cavity-file was generated by PARSER (from files_GN.cpp)" << "\n";
    last_atm = NULL;
    n_cav_atoms = 0;
    n_cav_conects = 0;
    cav_log.clear();
    write_cav_atoms(cav);
    write_cav_conect_lines();
    write_cav_master();
    f_out << "END" << endl;
    f_out.close();
    f_out.clear();
    if (verbosity > 1) cout << "finished writing " << filename << endl;
    return 1;
}

bool PARSER::write_pdb_cav_atoms_only(char const* file) {
    filename = file;
    if (verbosity > 1) cout << "writing " << filename << " ..." << endl;
    f_out.open(filename);
    if (!f_out) {cerr << c_message<cERROR>("PARSER::write_pdb_cav --> could not open ") << filename << " for writing!" << endl; f_out.clear(); return 0;}
    f_out << "REMARK    This cavity-file was generated by PARSER (from files_GN.cpp)" << "\n";
    
    write_atoms_special();
    
    f_out << "END" << endl;
    f_out.close();
    f_out.clear();
    if (verbosity > 1) cout << "finished writing " << filename << endl;
    return 1;
}

bool PARSER::write_pdb_cav_all(map<int,map<int,vec3d<float> > > &flex_map,char const* file) {
    filename = file;
    if (verbosity > 1) cout << "writing " << filename << " ..." << endl;
    f_out.open(filename);
    if (!f_out) {
        cerr << c_message<cERROR>("PARSER::write_pdb_cav_all --> could not open ") << filename << " for writing!" << endl;
        f_out.clear(); return 0;
    }
    f_out << "REMARK    This cavity-file was generated by PARSER (from files_GN.cpp)" << "\n";
    last_atm = NULL;
    n_cav_atoms = 0;
    n_cav_conects = 0;
    cav_log.clear();
    
    map<int,vec3d<float> > orig_pos;
    if (!(main_structure->protein)) return false;
    for (chains_vec cvt=main_structure->protein->chains.begin(); cvt!=main_structure->protein->chains.end(); ++cvt) {
        for (atoms_vec at=(*cvt)->atoms.begin(); at!=(*cvt)->atoms.end(); ++at) {
            orig_pos[(*at)->intern_id] = (*at)->coord;
        }
    }
            
    int modelnumber = 1;
    for (cavities_vec cav=main_structure->cavities.begin(); cav!=main_structure->cavities.end(); ++cav) {
        for (aacids_vec at=(*cav)->aacids.begin(); at!=(*cav)->aacids.end(); ++at) {
            for (atoms_vec it=(*at)->atoms.begin(); it!=(*at)->atoms.end(); ++it) {
                (*it)->coord = orig_pos[(*it)->intern_id];
                if (flex_map[modelnumber-1].find((*it)->intern_id) == flex_map[modelnumber-1].end()) continue;
                (*it)->coord = flex_map[modelnumber-1][(*it)->intern_id];
            }
        }
        f_out << "MODEL        " << modelnumber << "\n";
        write_cav_atoms(*cav);
        f_out << "ENDMDL" << "\n";
        ++modelnumber;
    }
    write_cav_conect_lines();
    write_cav_master();
    f_out << "END" << endl;
    f_out.close();
    f_out.clear();
    if (verbosity > 1) cout << "finished writing " << filename << endl;
    return 1;
}


//------------------------------------------------------------------------------------------------------
//mol2-Methoden fuer PARSER:
//------------------------------------------------------------------------------------------------------

void PARSER::parse_mol2_atoms(bool block_mode) {
    while (!f_in.eof()) {
        string row; //aktuelle Zeile
        ios::pos_type old_pos = f_in.tellg();
        getline(f_in,row); //Zeile einlesen
        if (row == "" || row[0] == '#') continue;
        if (row[0] == '@') {if (block_mode) {compute_mol2_line_nm(row,false,old_pos,true); break;} compute_mol2_line(row); break;}
        is.clear(); //streamobjekt zurcksetzen
        is.str(row);
        stl_ptr<ATOM> atm(new ATOM());//neues ATOM-Objekt erzeugen
        int adder = 0;
        try {
            atm->type = 5;
            is >> atm->id >> atm->name;
            
            atm->intern_id = cur_id; //fortlaufende id fuer interne Verwaltung
            cur_id += 1; //fortlaufende Nummer hochzaehlen
            
            is >> atm->coord[0] >> atm->coord[1] >> atm->coord[2];
            is >> atm->sybyl_type;
            
            //!NEU: 11.03.2010:  Lonepairs NICHT mit einlesen
            if (atm->sybyl_type == "LP" || atm->sybyl_type == "Lp") {
                atm.kill();
                cur_id -= 1;
                continue;
            }
            
            if (!is.eof()) {
            //    is >> atm->res_number; //als int gelesen
                //! 14.8.09: Leider steht da nicht immer ein int => kleiner Umstand
                is >> dummy;
                is_tmp.clear(); is_tmp.str(dummy);
                is_tmp >> atm->res_number;
                
                if (atm->res_number == 0) adder = 1;
                atm->res_number += adder;
            }
            
            if (!is.eof()) {
                is >> atm->res_name;
                
                //! neu: 08.04.09: Laenge des resname auf 3 beschraenken (fuer Sequenzalignment):
                //! --> NEIN! : lieber nur dort (im Alignment) beschneiden !!!
                //    if (atm->res_name.size() > 3) atm->res_name = atm->res_name.substr(0,3);
            }
            if (!is.eof()) is >> atm->charge; //als string gelesen
            if (!parse_hydrogens) { //kein Wasserstoff mit einlesen
                if (atm->sybyl_type == "H") {
                    atm.kill();
                    cur_id -= 1;
                    continue;
                } else if (atm->sybyl_type[0] == 'H') { //! 08.08.2010 fuer H. Typen
                    if (atm->sybyl_type[1] == '.') {
                        atm.kill();
                        cur_id -= 1;
                        continue;
                    }
                }
            }
            
            if (atm->res_name != last_res_name) {
                if (atm->name[0] == 'C') { //!neues root-atom
                    main_structure->ligands.back()->sub_map[atm->intern_id] = atm;
                    last_res_name = atm->res_name;
                }
            }
            
            h_log[atm->id] = atm->intern_id;
            main_structure->ligands.back()->atoms.push_back(atm);
        //    main_structure->atoms.push_back(atm);

        } catch(...) {
            cerr << c_message<cERROR>("PARSER::parse_mol2_atoms --> could not process mol2-line:") << endl << row << endl << "in file: " << filename << endl;
            atm.kill();
            cur_id -= 1;
        }
    }
}

void PARSER::parse_mol2_bonds(bool block_mode) {
    int help1, help2;
    cur_id = 1;
    if (!parse_bonds_flag) return;
    while (!f_in.eof()) {
        string row; //aktuelle Zeile
        ios::pos_type old_pos = f_in.tellg();
        getline(f_in,row); //Zeile einlesen
        if (row == "" || row[0] == '#') continue;
        if (row[0] == '@') {if (block_mode) {compute_mol2_line_nm(row,false,old_pos,true); break;} compute_mol2_line(row); break;}
        is.clear(); //streamobjekt zurcksetzen
        is.str(row);
        stl_ptr<BOND> bnd(new BOND()); //neues BOND-Objekt erzeugen
        try { //!Achtung: evtl. Ressourcenleck, wenn nach Erzeugung eines Obj. ein Fehler auftrat
            is >> bnd->id >> help1 >> help2 >> bnd->type;
            
            if (help1 == help2) { //! weil 2x gleiche ID in Schuwi-mol2 => segfault, weil unten kein ->to zugeordnet wird
                cerr << c_message<cWARNING>("PARSER::parse_mol2_bonds --> bullshit: ")  << row << endl << "in file: " << filename << endl;
                bnd.kill();
                continue;
            }
            
            bnd->id = cur_id;
            
            //jetzt die BONDs zuordnen - ist leider etwas zeitaufwendig, spart aber den Einsatz von maps fuer alle atome
            if ((h_log.find(help1) != h_log.end())&&(h_log.find(help2) != h_log.end())) {
                bool got1 = false; bool got2 = false;
                for (atoms_vec it=main_structure->ligands.back()->atoms.begin();
                                it!=main_structure->ligands.back()->atoms.end();++it) {
                    if ((*it)->intern_id == h_log[help1]) {
                        bnd->from = (*it);
                        if (got2) break;
                        else got1 = true;
                    } else if ((*it)->intern_id == h_log[help2]) {
                        bnd->to = (*it);
                        if (got1) break;
                        else got2 = true;
                    }
                }
                main_structure->ligands.back()->bonds.push_back(bnd);
                ++cur_id;
            } else bnd.kill();
            
        } catch(...) {
            cerr << c_message<cERROR>("PARSER::parse_mol2_bonds --> could not process mol2-line:") << endl << row << endl << "in file: " << filename << endl;
            bnd.kill();
        }
    }
}

void PARSER::parse_mol2_substructure() {
    //!wird anders ausgewertet
}

void PARSER::parse_mol2_comment(bool block_mode) {
    while (!f_in.eof()) {
        string row; //aktuelle Zeile
        ios::pos_type old_pos = f_in.tellg();
        getline(f_in,row); //Zeile einlesen
        if (row[0] == '@') {if (block_mode) {compute_mol2_line_nm(row,false,old_pos,true); break;} compute_mol2_line(row); break;}
        stl_ptr<COMMENT> new_com(new COMMENT());
        new_com->text = row;
        main_structure->ligands.back()->comments.push_back(new_com);
    }
}

void PARSER::parse_mol2_molecule() {
    if ((!multi)&&(!first_mol2)) {f_in.seekg(0,ios::end); return;}//zum Dateiende springen
    first_mol2 = false;
    last_res_name = "@@@";
    h_log.clear();
    cur_id = 1; //fuer jeden Liganden wieder auf 1 setzen
    string row;
    stl_ptr<LIGAND> lig(new LIGAND());
    lig->main_structure = main_structure;
    try { 
        getline(f_in,lig->name); //Name des Liganden
        getline(f_in,row);
        
        getline(f_in,row);
        is.clear(); is.str(row);
        is >> lig->type; //Typ des mol2-Molekuels
        
        getline(f_in,row);
        is.clear(); is.str(row);
        is >> lig->charge_type;
        main_structure->ligands.push_back(lig);
    } catch(...) {
        cerr << c_message<cERROR>("PARSER::parse_mol2_molecule --> could not process mol2-line:") << endl << row << endl << "in file: " << filename << endl;
        lig.kill();
    }
}

void PARSER::parse_mol2_crysin() {
    string row;
    const double to_rad = 0.0174532925199;
    double factor,a,b,c,alpha,beta,gamma;
    double cosa, cosb, cosg, sing;
    int sg,setting;
    getline(f_in,row);
    if (row == "") return;
    is.clear();
    is.str(row);
    is >> a >> b >> c >> alpha >> beta >> gamma >> sg >> setting;
    
    stl_ptr<CRYSIN> crys(new CRYSIN());
    
    crys->a = a; crys->b = b; crys->c = c;
    crys->wa = alpha; crys->wb = beta; crys->wc = gamma;
    
    alpha *= to_rad; beta *= to_rad; gamma *= to_rad; //ins Bogenmass
    if ((abs(alpha - beta) < 0.00001) && (abs(alpha - gamma) < 0.00001)) crys->rhombo = true;
    else crys->rhombo = false;
    
    cosa = cos(alpha); cosb = cos(beta); cosg = cos(gamma); sing = sin(gamma);
    factor = sqrt(1. - (cosa*cosa) - (cosb*cosb) - (cosg*cosg) + (2.*cosa*cosb*cosg) ); //Faktor zum berechnen der Transformationsmatrizen
    crys->group = sg;
    crys->setting = setting;
    
    crys->cryst2cart[0][0] = a; crys->cryst2cart[0][1] = b*cosg; crys->cryst2cart[0][2] = c*cosb;
    crys->cryst2cart[1][0] = 0.; crys->cryst2cart[1][1] = b*sing; crys->cryst2cart[1][2] = c*(cosa-(cosb*cosg)) / sing;
    crys->cryst2cart[2][0] = 0.; crys->cryst2cart[2][1] = 0.; crys->cryst2cart[2][2] = c * factor / sing;
    
    
    crys->cart2cryst[0][0] = 1./a; crys->cart2cryst[0][1] = -1.*cosg/(a*sing); crys->cart2cryst[0][2] = ((cosa*cosg)-cosb)/(a*sing*factor);
    crys->cart2cryst[1][0] = 0.; crys->cart2cryst[1][1] = 1./(b*sing); crys->cart2cryst[1][2] = ((cosb*cosg)-cosa)/(b*sing*factor); //!b oder beta?
//    crys->cart2cryst[1][0] = 0.; crys->cart2cryst[1][1] = 1./(b*sing); crys->cart2cryst[1][2] = ((cosb*cosg)-cosa)/(beta*sing*factor);
    //! b ist richtig!!!   //!08.08.07:  Tripos Dokumentation sagt aber beta!!!
    //! 11.09.07: Definitiv getestet (z.B.: UCOJOL): ES IST b und NICHT beta  => Fehler in Tripos Dokumentation

    crys->cart2cryst[2][0] = 0.; crys->cart2cryst[2][1] = 0.; crys->cart2cryst[2][2] = sing/(c*factor);
    
    main_structure->ligands.back()->crysin = crys;
}

void PARSER::compute_mol2_line(string const& row) {
    if (row.compare(9,4,"ATOM") == 0) parse_mol2_atoms();
    else if (row.compare(9,4,"BOND") == 0) parse_mol2_bonds();
    else if (row.compare(9,8,"MOLECULE") == 0) parse_mol2_molecule();
    else if (row.compare(9,6,"CRYSIN") == 0) parse_mol2_crysin();
    else if (parse_comments) {
        if (row.compare(9,7,"COMMENT") == 0) parse_mol2_comment();
    }
}

bool PARSER::read_mol2(char const* file,bool first_only,bool atoms_only) {
    if (first_only) multi = false; else multi = true;
    first_mol2 = true;
    filename = file;
    if (verbosity > 1) cout << "reading " << filename << " ..." << endl;
    f_in.open(filename);
    if (!f_in) {cerr << c_message<cERROR>("PARSER::read_mol2 --> could not open ") << filename << " !" << endl; f_in.clear(); return 0;}
    
    if (atoms_only) parse_bonds_flag = false;
    else parse_bonds_flag = true;
    
    while (!f_in.eof()) {
        string row; //aktuelle Zeile
        getline(f_in,row); //Zeile einlesen
        
        if (row[0] == '@') compute_mol2_line(row);
    }
    
    //!neu 08.01.07:
    if (!parse_water) { //!nachtraeglich entfernen ist effektiver als bei jeder Atomline zu pruefen
        for (ligands_vec lt=main_structure->ligands.begin(); lt!=main_structure->ligands.end(); ++lt) {
            (*lt)->kick_water();
        }
    }
    
    f_in.close();
    f_in.clear();
    if (verbosity > 1) cout << "finished reading " << filename << endl;
    return 1;
}

int PARSER::compute_mol2_line_nm(string const& row,bool next_allowed,ios::pos_type old_pos,bool returner) {
    if (row.compare(9,4,"ATOM") == 0) {parse_mol2_atoms(true); return 0;}
    else if (row.compare(9,4,"BOND") == 0) {parse_mol2_bonds(true); return 0;}
    else if (row.compare(9,8,"MOLECULE") == 0) {
        if (next_allowed) {parse_mol2_molecule(); return 1;}
        else {
            f_in.seekg(old_pos);
            if (returner) return 0; else return 1;
        }
    }
    else if (row.compare(9,6,"CRYSIN") == 0) {parse_mol2_crysin(); return 0;}
    else if (parse_comments) {
        if (row.compare(9,7,"COMMENT") == 0) {parse_mol2_comment(); return 0;}
    }
    return 0;
}

int PARSER::read_next_mol2(char const* file,int n_mols,bool atoms_only) {
    multi = true;
    filename = file;
    string scompi = file;
    
    if (still_open != scompi) {
        still_open = scompi;
        f_in.open(filename);
        if (!f_in) {
            cerr << c_message<cERROR>("PARSER::read_next_mol2 --> could not open ") << filename << " !" << endl;
            f_in.clear();
            return -1; // return 0;      //! changed: 30.09.08
        }
    }
    
    if (atoms_only) parse_bonds_flag = false;
    else parse_bonds_flag = true;
    
    int n_parsed = 0;
    bool clearer = true;
    bool next_allowed = true;
    while (!f_in.eof()) {
        string row; //aktuelle Zeile
        ios::pos_type old_pos = f_in.tellg();
        getline(f_in,row); //Zeile einlesen
        
        if (n_parsed == n_mols) next_allowed = false;
        
        if (row[0] == '@') n_parsed += compute_mol2_line_nm(row,next_allowed,old_pos);
        
        if (n_parsed > n_mols) {
            clearer = false;
            break;
        }
    }
        
    if (!parse_water) { //!nachtraeglich entfernen ist effektiver als bei jeder Atomline zu pruefen
        for (ligands_vec lt=main_structure->ligands.begin(); lt!=main_structure->ligands.end(); ++lt) {
            (*lt)->kick_water();
        }
    }
    
    if (clearer) {
        f_in.close();
        f_in.clear();
        return n_parsed;
    }
        
    return n_mols;
}

void PARSER::write_mol2_file(stl_ptr<LIGAND> const& ligand) {
    f_out << "@<TRIPOS>MOLECULE" << "\n";
    f_out << ligand->name << "\n";
    
    f_out.width(6); f_out << right << ligand->atoms.size();
    f_out.width(6); f_out << right << ligand->bonds.size() << "\n";
    f_out << ligand->type << "\n";
    f_out << ligand->charge_type << "\n";
    f_out << "****" << "\n";
    f_out << "This file was generated by PARSER (from files_GN.cpp)" << "\n" << "\n";
    
    f_out << "@<TRIPOS>ATOM" << "\n";
    
    rh_log.clear();
    curr_id = 1;
    
    int last_at_resnum = -1;
    last_res_number = 0;
    
    for (atoms_vec it=ligand->atoms.begin(); it!=ligand->atoms.end();++it) {
        f_out.setf(ios::fixed);
        f_out.width(7); f_out << right << curr_id << " "; 
        
        unsigned int ni=0; while ((*it)->name[ni++] == ' ') {if (ni >= (*it)->name.size()) break;}
        for (; ni<(*it)->name.size(); ++ni) if ((*it)->name[ni] == ' ') (*it)->name[ni] = '_'; //! kein Whitespace erlaubt!
        
        f_out.width(8); f_out << left << (*it)->name << " ";
        f_out.width(9); f_out.precision(4); f_out << right << (*it)->coord[0] << " ";
        f_out.width(9); f_out.precision(4); f_out << right << (*it)->coord[1] << " ";
        f_out.width(9); f_out.precision(4); f_out << right << (*it)->coord[2] << " ";
        f_out.width(7); f_out << left << (*it)->sybyl_type << " ";
        
        //! 29.04.09 -----------------
    //    f_out.width(5); f_out << right << (*it)->res_number << " ";
        if ((*it)->res_number < last_res_number) {
            if ((*it)->res_number != last_at_resnum) {
                last_at_resnum = (*it)->res_number;
                (*it)->res_number = last_res_number + 1;
            } else (*it)->res_number = last_res_number;
        }
        
        f_out.width(5); f_out << right << (*it)->res_number << " ";
        last_res_number = (*it)->res_number;
        //! --------------------------
        
        f_out.width(6); f_out << left << (*it)->res_name << " ";
        f_out.width(9); f_out << right << (*it)->charge << "\n";
        f_out.unsetf(ios::fixed);
        rh_log[(*it)->intern_id] = curr_id;
        ++curr_id;
    }
    
    curr_id = 1;
    if (ligand->bonds.size() > 0) f_out << "@<TRIPOS>BOND" << "\n";
    for (bonds_vec it=ligand->bonds.begin(); it!=ligand->bonds.end(); ++it) {
        if (rh_log.find((*it)->from->intern_id) == rh_log.end()) continue;
        if (rh_log.find((*it)->to->intern_id) == rh_log.end()) continue;
        f_out.width(6); f_out << right << curr_id; 
        f_out.width(6); f_out << right << rh_log[(*it)->from->intern_id];
        f_out.width(6); f_out << right << rh_log[(*it)->to->intern_id];
        f_out.width(4); f_out << right << (*it)->type << "\n";
        ++curr_id;
    }
    
    if (ligand->sub_map.size() > 0) {
        f_out << "@<TRIPOS>SUBSTRUCTURE" << "\n";
    //    int sn = 1;
        for (map<int,stl_ptr<ATOM> >::iterator it=ligand->sub_map.begin(); it!=ligand->sub_map.end();++it) {
    //        f_out.width(6); f_out << right << sn << " "; sn++;
            f_out.width(6); f_out << it->second->res_number << " "; //!muss laut Olli die res_number sein
        
            f_out.width(4); f_out << right << it->second->res_name << " ";
            f_out.width(6); f_out << right << it->first << "\n";
    //        f_out.width(6); f_out << right << rh_log[it->second->intern_id] << "\n";
        }
    }
    
    if (!(ligand->crysin.zero())) {
        f_out << "@<TRIPOS>CRYSIN" << "\n";
        f_out.setf(ios::fixed);
        f_out.width(9); f_out.precision(4); f_out << right << ligand->crysin->a << " ";
        f_out.width(9); f_out.precision(4); f_out << right << ligand->crysin->b << " ";
        f_out.width(9); f_out.precision(4); f_out << right << ligand->crysin->c << " ";
        f_out.width(9); f_out.precision(4); f_out << right << ligand->crysin->wa << " ";
        f_out.width(9); f_out.precision(4); f_out << right << ligand->crysin->wb << " ";
        f_out.width(9); f_out.precision(4); f_out << right << ligand->crysin->wc << " ";
        f_out << ligand->crysin->group << " " << ligand->crysin->setting << "\n";
        f_out.unsetf(ios::fixed);
    }
    if (ligand->comments.size() > 0) {    
        f_out << "@<TRIPOS>COMMENT" << "\n";
        for (comments_vec ct=ligand->comments.begin(); ct!=ligand->comments.end(); ++ct) {
            f_out << (*ct)->text << "\n";
        }
    }
}

void PARSER::write_mol2_file(PROTEIN *protein) {
    int cur_num = 0;
    
    string last_name = "XXX";
    f_out << "@<TRIPOS>MOLECULE" << "\n";
    f_out << "Protein file" << "\n";
    int n_atoms = 0;
    int n_bonds = 0;
    for (chains_vec it=protein->chains.begin(); it!=protein->chains.end();++it) {
        n_atoms += (*it)->atoms.size();
        n_bonds += (*it)->bonds.size();
    }
    f_out.width(6); f_out << right << n_atoms; //!falsch
    f_out.width(6); f_out << right << n_bonds << "\n"; //!falsch
    f_out << "PROTEIN" << "\n";
    f_out << "NO_CHARGES" << "\n";
    f_out << "****" << "\n";
    f_out << "This file was generated by PARSER (from files_GN.cpp)" << "\n" << "\n";
    f_out << "@<TRIPOS>DICT" << "\n";
    f_out << "PROTEIN PROTEIN" << "\n";
    f_out << "@<TRIPOS>ATOM" << "\n";
    
    //int test_size = 0;
    rh_log.clear();
    curr_id = 1;
    
    for (chains_vec jt=protein->chains.begin(); jt!=protein->chains.end();++jt) {
        for (atoms_vec it=(*jt)->atoms.begin(); it!=(*jt)->atoms.end();++it) {
            
            if (last_name != (*it)->full_res_name) {
                //++cur_num; //22.04.09: durch nachfolgende Zeile ersetzt
                cur_num = (*it)->res_number;
                last_name = (*it)->full_res_name;
            }
            f_out.setf(ios::fixed);
            
            f_out.width(7); f_out << right << curr_id << " "; 
            
            unsigned int ni=0; while ((*it)->name[ni++] == ' ') {if (ni >= (*it)->name.size()) break;}
            for (; ni<(*it)->name.size(); ++ni) if ((*it)->name[ni] == ' ') (*it)->name[ni] = '_'; //! kein Whitespace erlaubt!
            
            f_out.width(8); f_out << left << (*it)->name << " ";
            f_out.width(9); f_out.precision(4); f_out << right << (*it)->coord[0] << " ";
            f_out.width(9); f_out.precision(4); f_out << right << (*it)->coord[1] << " ";
            f_out.width(9); f_out.precision(4); f_out << right << (*it)->coord[2] << " ";
            f_out.width(7); f_out << left << (*it)->sybyl_type << " ";
            f_out.width(5); f_out << right << cur_num << " ";
            f_out.width(9); f_out << left << (*it)->full_res_name << " ";
            f_out.width(9); f_out << right << "0.0000 "; //bisher keine Ladungen
            f_out << (*it)->dict_type << "\n";
            f_out.unsetf(ios::fixed);
            rh_log[(*it)->intern_id] = curr_id;
            ++curr_id;
        }
    }
    curr_id = 1;
    f_out << "@<TRIPOS>BOND" << "\n";
    for (chains_vec jt=protein->chains.begin(); jt!=protein->chains.end();++jt) {
        for (bonds_vec it=(*jt)->bonds.begin(); it!=(*jt)->bonds.end();++it) {
            if (rh_log.find((*it)->from->intern_id) == rh_log.end()) continue;
            if (rh_log.find((*it)->to->intern_id) == rh_log.end()) continue;
            f_out.width(6); f_out << right << curr_id;
            f_out.width(6); f_out << right << rh_log[(*it)->from->intern_id];
            f_out.width(6); f_out << right << rh_log[(*it)->to->intern_id];
            f_out.width(4); f_out << right << (*it)->type << "   ";
            f_out << (*it)->dict_type << "\n";
            ++curr_id;
        }
    }
    f_out << "@<TRIPOS>SUBSTRUCTURE" << "\n";
    //int sn = 1;
    bool new_chain;
    for (chains_vec jt=protein->chains.begin(); jt!=protein->chains.end();++jt) {
        new_chain = true;
        for (aacids_map it=(*jt)->aacids.begin(); it!=(*jt)->aacids.end();++it) {
            if (rh_log.find(it->second->atoms[0]->intern_id) == rh_log.end()) continue;
            //f_out.width(6); f_out << right << sn << " "; ++sn; //22.04.09: durch nachfolgende Zeile ersetzt
            f_out.width(6); f_out << right << it->second->atoms[0]->res_number << " ";
            f_out.width(9); f_out << left << it->second->atoms[0]->full_res_name << " ";
            f_out.width(6); f_out << right << rh_log[it->second->atoms[0]->intern_id] << " "; //!fehlt noch check, ob es atoms[0] gibt!!!
            f_out << " RESIDUE   1 "; f_out << it->second->atoms[0]->chain_id << "  ";//!       -"-
            f_out << it->second->res_name;
            if (new_chain) {f_out << " 1  root" << "\n"; new_chain = false;}
            else {
                if (it->second->last_res) f_out << " 1" << "\n";
                else f_out << " 2" << "\n";
            }
        }
    }
}

bool PARSER::write_mol2(char const* file) {
    filename = file;
    if (verbosity > 1) cout << "writing " << filename << " ..." << endl;
    f_out.open(filename);
    if (!f_out) {cerr << c_message<cERROR>("PARSER::write_mol2 --> could not open ") << filename << " for writing!" << endl; f_out.clear(); return 0;}
    
    for (ligands_vec it=main_structure->ligands.begin(); it!=main_structure->ligands.end();++it) {
        write_mol2_file(*it);
        f_out << "\n";
    }
    
    f_out.close();
    f_out.clear();
    if (verbosity > 1) cout << "finished writing " << filename << endl;
    return 1;
}

bool PARSER::write_mol2(vector<stl_ptr<LIGAND> > &ligands,char const* file) {
    filename = file;
    if (verbosity > 1) cout << "writing " << filename << " ..." << endl;
    f_out.open(filename);
    if (!f_out) {cerr << c_message<cERROR>("PARSER::write_mol2 --> could not open ") << filename << " for writing!" << endl; f_out.clear(); return 0;}
    
    for (ligands_vec it=ligands.begin(); it!=ligands.end();++it) {
        write_mol2_file(*it);
        f_out << "\n";
    }
    
    f_out.close();
    f_out.clear();
    if (verbosity > 1) cout << "finished writing " << filename << endl;
    return 1;
}

bool PARSER::write_next_mol2(vector<stl_ptr<LIGAND> > &ligands,char const* file) {
    //filename = file;
    alt_filename = file;
    string scompi = file;
    //if (still_open != scompi) {
    if (alt_still_open != scompi) {
        //still_open = scompi;
        alt_still_open = scompi;
        //f_out.open(filename);
        f_out.open(alt_filename);
        
        if (!f_out) {
            //cerr << c_message<cERROR>("PARSER::write_next_mol2 --> could not open ") << filename << " for writing!" << endl;
            cerr << c_message<cERROR>("PARSER::write_next_mol2 --> could not open ") << alt_filename << " for writing!" << endl; 
            f_out.clear();
            return 0;
        }
    }
    
    for (ligands_vec it=ligands.begin(); it!=ligands.end(); ++it) {
        write_mol2_file(*it);
        f_out << endl;
    }
    return 1;
}

void PARSER::write_next_mol2_end() {
    f_out.close();
    f_out.clear();
}

bool PARSER::write_mol2_ligand(stl_ptr<LIGAND> const& ligand,char const* file) {
    filename = file;
    if (verbosity > 1) cout << "writing " << filename << " ..." << endl;
    f_out.open(filename);
    if (!f_out) {cerr << c_message<cERROR>("PARSER::write_mol2_ligand --> could not open ") << filename << " for writing!" << endl; f_out.clear(); return 0;}

    write_mol2_file(ligand);

    f_out.close();
    f_out.clear();
    if (verbosity > 1) cout << "finished writing " << filename << endl;
    return 1;
}

bool PARSER::write_mol2_protein(PROTEIN *protein,char const* file) {
    filename = file;
    if (verbosity > 1) cout << "writing " << filename << " ..." << endl;
    f_out.open(filename);
    if (!f_out) {cerr << c_message<cERROR>("PARSER::write_mol2_protein --> could not open ") << filename << " for writing!" << endl; f_out.clear(); return 0;}

    write_mol2_file(protein);

    f_out.close();
    f_out.clear();
    if (verbosity > 1) cout << "finished writing " << filename << endl;
    return 1;
}


//------------------------------------------------------------------------------------------------------
//dlg-Methoden fuer PARSER:
//------------------------------------------------------------------------------------------------------

bool PARSER::read_dlg(char const* file,int high, int low,int read_energy,bool get_flex_res,
                      float cut_radius,char const* pre_cav_name) {
    if (verbosity > 1) cout << "reading " << file << " ..." << endl;
    f_in.open(file);
    filename = file;
    if (!f_in) {
        cerr << c_message<cERROR>("PARSER::read_dlg --> could not load ") << filename << endl;
        return false;
    }
    int curr_id = 1;
    string row,key,key2,d1,d2,d3,last_res,at_name_buf;
    int last_model = 0;
    int new_model = 0;
    int mol_number = 0;
    float curr_energy = 0.;
    float curr_energy2 = 0.;
    //float x,y,z;
    bool is_actual = false;
    bool first = false;
    bool first2 = false;
    while (!f_in.eof()) {
        getline(f_in,row);
        is.clear(); is.str(row);
        key = "X";
        is >> key;
        if (key != "DOCKED:") continue;
        is >> key;
        if (key == "MODEL") {
            is >> new_model;
            ++mol_number;
            if (mol_number >= low && mol_number <= high) {
                if (main_structure->ligands.size() > 0) {
                    //Dem letzten Liganden einen Namen geben:
                    ostringstream hos;
                    hos << "_model_" << last_model;
                    if (read_energy == 1) hos << " " << curr_energy;
                    else if (read_energy == 2) hos << " " << curr_energy2;
                    else if (read_energy == 3) hos << " " << curr_energy << " " << curr_energy2;
                    main_structure->ligands.back()->name = hos.str(); //!spaeter den eigentlichen Namen voranstellen
                }
                is_actual = true;
                stl_ptr<LIGAND> lig(new LIGAND());
                lig->main_structure = main_structure; //!neu: 11.10.07
                first = true;
                first2 = true;
                main_structure->ligands.push_back(lig);
                curr_energy = 0.;
                curr_energy2 = 0.;
                last_model = new_model;
            } else is_actual = false;
        } else if (key == "USER") {
            d1 = "X";
            is >> d1;
            if (d1 == "Final") {
                d2 = "X";
                is >> d2;
                if (d2 == "Docked") {
                    d3 = "X";
                    is >> d3;
                    if (d3 == "Energy") {
                        is >> d1 >> curr_energy;
                    }
                }
            } else if (d1 == "Estimated") {
                d2 = "X";
                is >> d2;
                if (d2 == "Free") {
                    d3 = "X";
                    is >> d3;
                    if (d3 == "Energy") {
                        d3 = "X";
                        is >> d3;
                        if (d3 == "of") {
                            d3 = "X";
                            is >> d3;
                            if (d3 == "Binding") {
                                is >> d1 >> curr_energy2;
                            }
                        }
                    }
                }
            }
        } else if (key == "ATOM" && is_actual) {
            string s_id(row,14,5), s_name(row,20,4), s_res(row,25,3), s_resnum(row,30,4),
                   s_x(row,38,8), s_y(row,46,8), s_z(row,54,8), s_charge(row,76,8);
            
            string s12 = "X";
            try { //siehe ATOM
                s12.assign(row,85,2);
            } catch(...) {};
            
            stl_ptr<ATOM> atm(new ATOM());
            is.clear(); is.str(s_id); is >> atm->id;
            atm->intern_id = curr_id;
            ++curr_id;
            is.clear(); is.str(s_name); is >> atm->name;
            atm->chain_id = row[29];
            
            atm->res_name = s_res;
            if (first) {
                first = false;
                last_res = atm->res_name;
            }
            is.clear(); is.str(s_resnum); is >> atm->res_number;
            is.clear(); is.str(s_x); is >> atm->coord[0];
            is.clear(); is.str(s_y); is >> atm->coord[1];
            is.clear(); is.str(s_z); is >> atm->coord[2];
            is.clear(); is.str(s_charge); is >> atm->charge;
            
            if (s12 == "X") {
                d2 = atm->name;
            } else {
                if (s12[0] > 64) { //! neu 03.08.2010
                //    cerr << "using element column: " << s12 << endl;
                    d2 = s12;
                } else d2 = atm->name;
            }
            
            if (d2[0] == 'H') {
                atm->element = "H";
            }
            else if (d2[0] == 'A') {
                atm->element = "C";
                d2[0] = 'C';
            }
            else if (d2[0] == 'C') {
                atm->element = "C";
                if (d2.size() > 1) {
                    if (d2[1] == 'l') atm->element = "Cl";
                }
            }
            else if (d2[0] == 'b') {
                atm->element = "Br";
            }
            else if (d2[0] == 'B') {
                atm->element = "B";
                if (d2.size() > 1) {
                    if (d2[1] == 'r') atm->element = "Br";
                }
            }
            else if (d2[0] == 'c') {
                atm->element = "Cl";
            }
            else if (d2[0] == 'j' || d2[0] == 'I') {
                atm->element = "I";
            }
            else if (d2[0] == 'f' || d2[0] == 'F') {
                atm->element = "F";
            } else {
                //!nochmal pruefen wie S, P, usw. in den dlg files heissen
                atm->element.assign(1,d2[0]);
                if ((d2.size() > 1) && !(d2[0] == 'O' || d2[0] == 'N' || d2[0] == 'S' ||
                     d2[0] == 'P' || d2[0] == 'F' || d2[0] == 'I')) atm->element += d2[1];
            }
            if (get_flex_res) {
                if (atm->res_name != last_res /*&& atm->chain_id != ' '*/) {
                    if (first2) {
                        if (main_structure->cavities.size() > 0) {
                            string cav_name = pre_cav_name;
                            cav_name = cav_name + 
                            main_structure->ligands[main_structure->ligands.size()-2]->name + "_poc.pdb";
                            write_pdb_cav(cav_name.c_str(),main_structure->cavities.back());
                        }
                        main_structure->cavities.clear();
                        main_structure->get_lig_cavity(cut_radius,&(*main_structure->ligands.back()));
                        first2 = false;
                    }
                    for (aacids_vec ut=main_structure->cavities.back()->aacids.begin();
                                    ut!=main_structure->cavities.back()->aacids.end();++ut) {
                        if ((*ut)->res_name == atm->res_name && (*ut)->res_number == atm->res_number) {
                            for (atoms_vec at=(*ut)->atoms.begin(); at!=(*ut)->atoms.end(); ++at) {
                                is.clear(); is.str((*at)->name); is >> at_name_buf;
                                if (at_name_buf == atm->name) {
                                    (*at)->coord = atm->coord;
                                    break;
                                }
                            }
                            break;
                        }
                    }
                } else {
                    main_structure->ligands.back()->atoms.push_back(atm);
                    last_res = atm->res_name;
                }
            } else {
                main_structure->ligands.back()->atoms.push_back(atm);
                last_res = atm->res_name;
            }
        }
    }
    //Dem letzten Liganden einen Namen geben:
    if (main_structure->ligands.size() > 0) {
        ostringstream hos;
        hos << "_model_" << last_model;
        if (read_energy == 1) hos << " " << curr_energy;
        else if (read_energy == 2) hos << " " << curr_energy2;
        else if (read_energy == 3) hos << " " << curr_energy << " " << curr_energy2;
        main_structure->ligands.back()->name = hos.str(); //!spaeter den eigentlichen Namen voranstellen
        if (main_structure->cavities.size() > 0) {
            string cav_name = pre_cav_name;
            cav_name = cav_name + 
            main_structure->ligands[main_structure->ligands.size()-1]->name + "_poc.pdb";
            write_pdb_cav(cav_name.c_str(),main_structure->cavities.back());
        }
    }
    f_in.close();
    f_in.clear();
    if (verbosity > 1) cout << "finished reading " << file << endl;
    return true;
}

//Das ganze nochmal fuer dsX --- sehr unschoen => spaeter mal EINE generische Form draus machen
bool PARSER::read_dlg(map<int,map<int,vec3d<float> > > &flex_map,PROTEIN *protein,
                      char const* file,bool get_flex_res) {
    if (verbosity > 1) cout << "reading " << file << " ..." << endl;
    f_in.open(file);
    filename = file;
    if (!f_in) {
        cerr << c_message<cERROR>("PARSER::read_dlg --> could not load ") << filename << endl;
        return false;
    }
    int curr_id = 1;
    string row,key,key2,d1,d2,d3,last_res,at_name_buf;
    int last_model = 0;
    int new_model = 0;
    int mol_number = 0;
    bool is_actual = false;
    bool first = false;
    bool first2 = false;
    while (!f_in.eof()) {
        getline(f_in,row);
        is.clear(); is.str(row);
        key = "X";
        is >> key;
        if (key != "DOCKED:") continue;
        is >> key;
        if (key == "MODEL") {
            is >> new_model;
            ++mol_number;
            is_actual = true;
            stl_ptr<LIGAND> lig(new LIGAND());
            lig->main_structure = main_structure; //!neu: 11.10.07
            first = true;
            first2 = true;
            main_structure->ligands.push_back(lig);
            last_model = new_model;
        } else if (key == "ATOM" && is_actual) {
            string s_id(row,14,5), s_name(row,20,4), s_res(row,25,3), s_resnum(row,30,4),
                   s_x(row,38,8), s_y(row,46,8), s_z(row,54,8), s_charge(row,76,8);
            
            string s12 = "X";
            try { //siehe ATOM
                s12.assign(row,85,2);
            } catch(...) {};
            
            stl_ptr<ATOM> atm(new ATOM());
            is.clear(); is.str(s_id); is >> atm->id;
            atm->intern_id = curr_id;
            ++curr_id;
            is.clear(); is.str(s_name); is >> atm->name;
            atm->chain_id = row[29];
            atm->res_name = s_res;
            if (first) {
                first = false;
                last_res = atm->res_name;
            }
            is.clear(); is.str(s_resnum); is >> atm->res_number;
            is.clear(); is.str(s_x); is >> atm->coord[0];
            is.clear(); is.str(s_y); is >> atm->coord[1];
            is.clear(); is.str(s_z); is >> atm->coord[2];
            is.clear(); is.str(s_charge); is >> atm->charge;
            
            if (s12 == "X") {
                d2 = atm->name;
            } else {
                if (s12[0] > 64) { //! neu 03.08.2010
                    d2 = s12;
                } else d2 = atm->name;
            }
            
            if (d2[0] == 'H') {
                atm->element = "H";
            }
            else if (d2[0] == 'A') {
                atm->element = "C";
                d2[0] = 'C';
            }
            else if (d2[0] == 'C') {
                atm->element = "C";
                if (d2.size() > 1) {
                    if (d2[1] == 'l') atm->element = "Cl";
                }
            }
            else if (d2[0] == 'b') {
                atm->element = "Br";
            }
            else if (d2[0] == 'B') {
                atm->element = "B";
                if (d2.size() > 1) {
                    if (d2[1] == 'r') atm->element = "Br";
                }
            }
            else if (d2[0] == 'c') {
                atm->element = "Cl";
            }
            else if (d2[0] == 'j' || d2[0] == 'I') {
                atm->element = "I";
            }
            else if (d2[0] == 'f' || d2[0] == 'F') {
                atm->element = "F";
            } else {
                //!nochmal pruefen wie S, P, usw. in den dlg files heissen
                atm->element.assign(1,d2[0]);
                if ((d2.size() > 1) && !(d2[0] == 'O' || d2[0] == 'N' || d2[0] == 'S' ||
                     d2[0] == 'P' || d2[0] == 'F' || d2[0] == 'I')) atm->element += d2[1];
            }
            if (get_flex_res) {
                if (atm->res_name != last_res /*&& atm->chain_id != ' '*/) {
                    
                    for (chains_vec cv=protein->chains.begin(); cv!=protein->chains.end(); ++cv) {
                    for (aacids_map ut=(*cv)->aacids.begin(); ut!=(*cv)->aacids.end(); ++ut) {
                        if (ut->second->res_name == atm->res_name && ut->second->res_number == atm->res_number) {
                            for (atoms_vec at=ut->second->atoms.begin(); at!=ut->second->atoms.end(); ++at) {
                                is.clear(); is.str((*at)->name); is >> at_name_buf;
                                if (at_name_buf == atm->name) {
                                    flex_map[mol_number-1].insert(pair<int,vec3d<float> >((*at)->intern_id,
                                                                                        atm->coord));    
                                    break;
                                }
                            }
                            break;
                        }
                    }
                    }
                } else {
                    main_structure->ligands.back()->atoms.push_back(atm);
                    last_res = atm->res_name;
                }
            } else {
                main_structure->ligands.back()->atoms.push_back(atm);
                last_res = atm->res_name;
            }
        }
        
    }
    
    f_in.close();
    f_in.clear();
    if (verbosity > 1) cout << "finished reading " << file << endl;
    return true;
}

bool PARSER::read_dlg_for_full_pdb(char const* file,int high,int low,char const* pre_cav_name) {
    if (verbosity > 1) cout << "reading " << file << " ..." << endl;
    f_in.open(file);
    filename = file;
    if (!f_in) {
        cerr << c_message<cERROR>("PARSER::read_dlg --> could not load ") << filename << endl;
        return false;
    }
    string row,key,key2,d1,last_res,at_name_buf;
    int last_model = 0;
    int new_model = 0;
    int mol_number = 0;
    bool is_actual = false;
    bool first = true;
    
    map<stl_ptr<ATOM>,vec3d<float> > save_coord;
    
    while (!f_in.eof()) {
        getline(f_in,row);
        is.clear(); is.str(row);
        key = "X";
        is >> key;
        if (key != "DOCKED:") continue;
        is >> key;
        if (key == "MODEL") {
            is >> new_model;
            ++mol_number;
            if (mol_number >= low && mol_number <= high) is_actual = true;
            else is_actual = false;
            if (!first && is_actual) {
                ostringstream os;
                os << pre_cav_name << "_model_" << last_model << ".pdb";
                string pname = os.str();
                write_pdb(pname.c_str());
                for (map<stl_ptr<ATOM>,vec3d<float> >::iterator it=save_coord.begin(); it!=save_coord.end(); ++it) {
                    it->first->coord = it->second;
                }
            } else if (is_actual) first = false;
            if (is_actual) last_model = new_model;
        } else if (key == "ATOM" && is_actual) {
            string s_id(row,14,5), s_name(row,20,4), s_res(row,25,3), s_resnum(row,30,4),
                   s_x(row,38,8), s_y(row,46,8), s_z(row,54,8), s_charge(row,76,8);
            
            string s12 = "X";
            try { //siehe ATOM
                s12.assign(row,85,2);
            } catch(...) {};
            
            stl_ptr<ATOM> atm(new ATOM());
            is.clear(); is.str(s_id); is >> atm->id;
            is.clear(); is.str(s_name); is >> atm->name;
            atm->chain_id = row[29];
            atm->res_name = s_res;
            is.clear(); is.str(s_resnum); is >> atm->res_number;
            is.clear(); is.str(s_x); is >> atm->coord[0];
            is.clear(); is.str(s_y); is >> atm->coord[1];
            is.clear(); is.str(s_z); is >> atm->coord[2];
            if (s12 != "X") atm->element = s12;
            
            bool breaker = false;
            for (ligands_vec lig=main_structure->ligands.begin(); lig!=main_structure->ligands.end();++lig) {
                for (atoms_vec at=(*lig)->atoms.begin(); at!=(*lig)->atoms.end(); ++at) {
                    if ((*at)->res_name != atm->res_name) continue;
                    is.clear(); is.str((*at)->name); is >> at_name_buf;
                    if (at_name_buf == atm->name && (*at)->res_number == atm->res_number) {
                        save_coord[*at] = vec3d<float>((*at)->coord);
                        (*at)->coord = atm->coord;
                        breaker = true;
                        break;
                    }
                }
                if (breaker) break;
            }
            if (!breaker) for (chains_vec ch=main_structure->protein->chains.begin(); ch!=main_structure->protein->chains.end(); ++ch) {
                for (atoms_vec at=(*ch)->atoms.begin(); at!=(*ch)->atoms.end(); ++at) {
                    if ((*at)->res_name != atm->res_name) continue;
                    is.clear(); is.str((*at)->name); is >> at_name_buf;
                    if (at_name_buf == atm->name && (*at)->res_number == atm->res_number) {
                        save_coord[*at] = vec3d<float>((*at)->coord);
                        (*at)->coord = atm->coord;
                        breaker = true;
                        break;
                    }
                }
                if (breaker) break;
            }
            atm.kill();
            
        }
        
    }
    
    ostringstream os;
    os << pre_cav_name << "_model_" << last_model << ".pdb";
    string pname = os.str();
    write_pdb(pname.c_str());
    
    f_in.close();
    f_in.clear();
    if (verbosity > 1) cout << "finished reading " << file << endl;
    return true;
}

//------------------------------------------------------------------------------------------------------
//sdf-Methoden fuer PARSER:
//------------------------------------------------------------------------------------------------------
bool PARSER::read_sdf(char const* file,bool read_energy) {
    if (verbosity > 1) cout << "reading " << file << " ..." << endl;
    f_in.open(file);
    filename = file;
    if (!f_in) {
        cerr << c_message<cERROR>("PARSER::read_sdf --> could not load ") << filename << endl;
        return false;
    }
    
    string row,key,name_buf;
    int n_atoms = 0;
    //bool first;
    while (!f_in.eof()) {
        getline(f_in,row);
        if (f_in.eof()) break;

        name_buf = "X";
        if (row.size() > 0) {
            string_fu::replace_char(row,' ','_');
            name_buf = row;
        }
        
        getline(f_in,row);
        if (f_in.eof()) break;
        getline(f_in,row);
        if (f_in.eof()) break;
        
        stl_ptr<LIGAND> lig(new LIGAND());
        lig->name = name_buf;
        lig->main_structure = main_structure; //!neu: 11.10.07
        lig->ele_already_set = true; //! GN: 20.01.11
        int a_intern_id = 1;
        getline(f_in,row);
        if (f_in.eof()) {
            cerr << c_message<cERROR>("PARSER::read_sdf --> there is a problem with ") << filename << endl;
            lig.kill();
            return false;
        }
        try {
            string s1(row,0,3);
            is.clear(); is.str(s1);
            n_atoms = 0;
            is >> n_atoms;
            
        } catch(...) {
            cerr << c_message<cERROR>("PARSER::read_sdf --> could not read atom number in ") << filename << endl;
            lig.kill();
            return false;
        }
        //jetzt die Atome einlesen:
        for (int i=0; i<n_atoms; ++i) {
            getline(f_in,row);
            if (f_in.eof()) {
                cerr << c_message<cERROR>("PARSER::read_sdf --> could not find ") << n_atoms << " in " << filename << endl;
                lig.kill();
                return false;
            }
            stl_ptr<ATOM> atm(new ATOM());
            string sx(row,0,10), sy(row,10,10), sz(row,20,10), se(row,31,3);
            is.clear(); is.str(sx);
            is >> atm->coord[0];
            is.clear(); is.str(sy);
            is >> atm->coord[1];
            is.clear(); is.str(sz);
            is >> atm->coord[2];
            is.clear(); is.str(se);
            is >> atm->element;

            if (atm->element.size() > 1) {
                if (atm->element[1] < 97) atm->element[1] += 32;
            }

            atm->name = atm->element;
            atm->type = 5;
            atm->intern_id = a_intern_id; a_intern_id++;
            atm->id = atm->intern_id;
            atm->res_number = 0;
            lig->atoms.push_back(atm);
        }
        //jetzt den molekuelnamen und eventuell die energy lesen:
        if (lig->name == "X") {
            ostringstream cnam; cnam << "mol " << glob_sdf_read; ++glob_sdf_read;
            lig->name = cnam.str();
        }
        while (!f_in.eof()) {
            getline(f_in,row);
            is.clear(); is.str(row);
            key = "*";
            name_buf = "*";
            is >> key;
            if (key == ">") is >> name_buf;
            else {
                if (key == "$$$$") break;
                continue;
            }
            
            if (name_buf.size() == 0) continue;
            if (name_buf.find("name") != name_buf.npos ||
                name_buf.find("NAME") != name_buf.npos ||
                name_buf.find("Name") != name_buf.npos) {
                getline(f_in,row);
                if (row.size() > 0) {
                    string_fu::replace_char(row,' ','_');
                    lig->name = row;
                }
            }
            try {
                if (name_buf.find("score") != name_buf.npos ||
                name_buf.find("Score") != name_buf.npos ||
                name_buf.find("SCORE") != name_buf.npos ||
                name_buf.find("energy") != name_buf.npos ||
                name_buf.find("Energy") != name_buf.npos ||
                name_buf.find("ENERGY") != name_buf.npos) {
                    getline(f_in,row);
                    is.clear(); is.str(row);
                    float energy = 0.;
                    is >> energy;
                    ostringstream hos;
                    hos << lig->name << " " << energy;
                    lig->name = hos.str();
                }
            } catch(...) {}
        }
        main_structure->ligands.push_back(lig);
    }
    
    f_in.close();
    f_in.clear();
    if (verbosity > 1) cout << "finished reading " << file << endl;
    return true;
}


int PARSER::read_next_sdf(char const* file,int n_mols) {
    multi = true;
    filename = file;
    string scompi = file;

    if (still_open != scompi) {
        still_open = scompi;
        f_in.clear();
        f_in.open(filename);
        glob_sdf_read = 0;
        if (!f_in) {
            cerr << c_message<cERROR>("PARSER::read_next_sdf --> could not open ") << filename << " !" << endl;
            f_in.clear();
            return -1;
        }
    }

    int n_parsed = 0;
    string row,key,name_buf;
    int n_atoms = 0;
    while (!f_in.eof()) {
        getline(f_in,row);
        if (f_in.eof()) break;

        name_buf = "X";
        if (row.size() > 0) {
            string_fu::replace_char(row,' ','_');
            name_buf = row;
        }

        getline(f_in,row);
        if (f_in.eof()) break;
        getline(f_in,row);
        if (f_in.eof()) break;

        stl_ptr<LIGAND> lig(new LIGAND());
        lig->name = name_buf;
        lig->main_structure = main_structure;
        lig->ele_already_set = true; //! GN: 20.01.11
        int a_intern_id = 1;
        getline(f_in,row);
        if (f_in.eof()) {
            cerr << c_message<cERROR>("PARSER::read_next_sdf --> there is a problem with ") << filename << endl;
            lig.kill();
            return n_parsed;
        }
        try {
            string s1(row,0,3);
            is.clear(); is.str(s1);
            n_atoms = 0;
            is >> n_atoms;

        } catch(...) {
            cerr << c_message<cERROR>("PARSER::read_next_sdf --> could not read atom number in ") << filename << endl;
            lig.kill();
            return n_parsed;
        }
        //jetzt die Atome einlesen:
        for (int i=0; i<n_atoms; ++i) {
            getline(f_in,row);
            if (f_in.eof()) {
                cerr << c_message<cERROR>("PARSER::read_next_sdf --> could not find ") << n_atoms << " in " << filename << endl;
                lig.kill();
                return n_parsed;
            }
            stl_ptr<ATOM> atm(new ATOM());
            string sx(row,0,10), sy(row,10,10), sz(row,20,10), se(row,31,3);
            is.clear(); is.str(sx);
            is >> atm->coord[0];
            is.clear(); is.str(sy);
            is >> atm->coord[1];
            is.clear(); is.str(sz);
            is >> atm->coord[2];
            is.clear(); is.str(se);
            is >> atm->element;

            if (atm->element.size() > 1) {
                if (atm->element[1] < 97) atm->element[1] += 32;
            }

            atm->name = atm->element;
            atm->type = 5;
            atm->intern_id = a_intern_id; a_intern_id++;
            atm->id = atm->intern_id;
            atm->res_number = 0;
            lig->atoms.push_back(atm);
        }
        //jetzt den molekuelnamen lesen (sonst default-Namen nehmen):
        if (lig->name == "X") {
            ostringstream cnam; cnam << "mol " << glob_sdf_read; ++glob_sdf_read;
            lig->name = cnam.str();
        }
        while (!f_in.eof()) {
            getline(f_in,row);
            is.clear(); is.str(row);
            key = "*";
            name_buf = "*";
            is >> key;
            if (key == ">") is >> name_buf;
            else {
                if (key == "$$$$") break;
                continue;
            }

            if (name_buf.find("name") != name_buf.npos ||
                name_buf.find("NAME") != name_buf.npos ||
                name_buf.find("Name") != name_buf.npos) {
                getline(f_in,row);
                if (row.size() > 0) {
                    string_fu::replace_char(row,' ','_');
                    lig->name = row;
                }
            }
        }
        main_structure->ligands.push_back(lig);
        ++n_parsed;

        if (n_parsed == n_mols) return n_parsed;
    }

    f_in.close();
    f_in.clear();
    return n_parsed;
}


bool PARSER::s2cryst(string const& s,stl_ptr<CRYSIN_POSITION> const& ccp) {
    //!erstmal splitten:
    vector<string> sv;
    if (string_fu::mysplit(s,sv,',') < 3) {
        cerr << c_message<cERROR>("PARSER::s2cryst --> could not parse ") << s << endl;
        return false;
    }
    for (int i=0; i<3; ++i) {
        unsigned int ii=0;
        float mod = 1.;
        while (ii<sv[i].size()) {
            unsigned int last_ii = ii;
            if (sv[i][ii] == '-') {mod *= -1.; ++ii;}
            
            if (sv[i][ii] == '1') {++ii;}
            else if (sv[i][ii] == '2') {mod *= 2.; ++ii;}
            else if (sv[i][ii] == '3') {mod *= 3.; ++ii;}
            else if (sv[i][ii] == '4') {mod *= 4.; ++ii;}
            else if (sv[i][ii] == '5') {mod *= 5.; ++ii;}
            else if (sv[i][ii] == '6') {mod *= 6.; ++ii;}
            else if (sv[i][ii] == '7') {mod *= 7.; ++ii;}
            else if (sv[i][ii] == '8') {mod *= 8.; ++ii;}
            else if (sv[i][ii] == '9') {mod *= 9.; ++ii;}
            
            if (sv[i][ii] == '/') {
                ++ii;
                if (sv[i][ii] == '2') {mod /= 2.; ++ii;}
                else if (sv[i][ii] == '3') {mod /= 3.; ++ii;}
                else if (sv[i][ii] == '4') {mod /= 4.; ++ii;}
                else if (sv[i][ii] == '5') {mod /= 5.; ++ii;}
                else if (sv[i][ii] == '6') {mod /= 6.; ++ii;}
                else if (sv[i][ii] == '7') {mod /= 7.; ++ii;}
                else if (sv[i][ii] == '8') {mod /= 8.; ++ii;}
                else if (sv[i][ii] == '9') {mod /= 9.; ++ii;}
            }
            
            if (sv[i][ii] == '+') {ccp->const_add[i] = mod; mod = 1.;++ii;}
            else if (sv[i][ii] == '-') {ccp->const_add[i] = mod; mod = 1.;} //!ACHTUNG: hier ii NICHT hochzaehlen
            else if (sv[i][ii] == 'x') {ccp->x_add[i] = mod; mod = 1.; ++ii;}
            else if (sv[i][ii] == 'y') {ccp->y_add[i] = mod; mod = 1.; ++ii;}
            else if (sv[i][ii] == 'z') {ccp->z_add[i] = mod; mod = 1.; ++ii;}
            
            if (ii == last_ii) {
                cerr << c_message<cERROR>("PARSER::s2cryst --> could not parse ") << s << endl;
                return false;
            }
        }
    }
    return true;
}

bool PARSER::read_cif_molecule(ios::pos_type &old_pos) {
    bool get_line = true;
    string row;
    bool finished = false;
    while (!f_in.eof()) {
        if (get_line) row = "";
        string key;
        if (get_line) {
            old_pos = f_in.tellg();
            getline(f_in,row);
        }
        get_line = true;
        if (row == "" || row[0] == '#') continue;
        is.clear(); is.str(row);
        is >> key;
        if (key == "_database_code_CSD") {
            if (finished) return true;
            else finished = true;
            stl_ptr<LIGAND> lig(new LIGAND());
            stl_ptr<CRYSIN> crys(new CRYSIN());
            lig->crysin = crys;
            lig->main_structure = main_structure;
            is >> lig->name;
            main_structure->ligands.push_back(lig);
        }
        else if (key == "_symmetry_Int_Tables_number") {
            is >> main_structure->ligands.back()->crysin->group;
        }
        else if (key == "_cell_length_a") {
            is >> main_structure->ligands.back()->crysin->a; //! pruefen, ob der komplette Wert bis zur Klammer gelesen wird
        }
        else if (key == "_cell_length_b") {
            is >> main_structure->ligands.back()->crysin->b;
        }
        else if (key == "_cell_length_c") {
            is >> main_structure->ligands.back()->crysin->c;
        }
        else if (key == "_cell_length_a_pm") {
            is >> main_structure->ligands.back()->crysin->a;
            main_structure->ligands.back()->crysin->a /= 100.;
        }
        else if (key == "_cell_length_b_pm") {
            is >> main_structure->ligands.back()->crysin->b;
            main_structure->ligands.back()->crysin->a /= 100.;
        }
        else if (key == "_cell_length_c_pm") {
            is >> main_structure->ligands.back()->crysin->c;
            main_structure->ligands.back()->crysin->a /= 100.;
        }
        else if (key == "_cell_length_a_nm") {
            is >> main_structure->ligands.back()->crysin->a;
            main_structure->ligands.back()->crysin->a *= 10.;
        }
        else if (key == "_cell_length_b_nm") {
            is >> main_structure->ligands.back()->crysin->b;
            main_structure->ligands.back()->crysin->a *= 10.;
        }
        else if (key == "_cell_length_c_nm") {
            is >> main_structure->ligands.back()->crysin->c;
            main_structure->ligands.back()->crysin->a *= 10.;
        }
        else if (key == "_cell_angle_alpha") {
            float atmp = 450.;
            is >> atmp;
            if (atmp > 400.) main_structure->ligands.back()->crysin->wa = 90.;
            else main_structure->ligands.back()->crysin->wa = atmp;
        }
        else if (key == "_cell_angle_beta") {
            float atmp = 450.;
            is >> atmp;
            if (atmp > 400.) main_structure->ligands.back()->crysin->wb = 90.;
            else main_structure->ligands.back()->crysin->wb = atmp;
        }
        else if (key == "_cell_angle_gamma") {
            float atmp = 450.;
            is >> atmp;
            if (atmp > 400.) main_structure->ligands.back()->crysin->wc = 90.;
            else main_structure->ligands.back()->crysin->wc = atmp;
        }
        else if (key == "_symmetry_equiv_pos_as_xyz") {
            //! Positions parsen:
            bool breaker = false;
            while (!breaker) {
                string row2;
                getline(f_in,row2);
                if (row2 == "" || row2[0] == '_') {
                    row = row2;
                    get_line = false;
                    break;
                }
                is.clear(); is.str(row2);
                int num = 0;
                is >> num; if (num < 2) continue;
                string sposi = "";
                is >> sposi;
                stl_ptr<CRYSIN_POSITION> ccp = new CRYSIN_POSITION;
                if (!s2cryst(sposi,ccp)) return false;
                //! neuen posi_vec in CRYSIN einhaengen:
                main_structure->ligands.back()->crysin->positions.push_back(ccp);
            }
        }
        else if (key == "_atom_site_label") {
            getline(f_in,row);
            if (row != "_atom_site_type_symbol") {
                cerr << c_message<cERROR>("PARSER::read_cif --> atom lines are not in format 'name ele x y z") << endl;
                return false;
            }
            getline(f_in,row);
            if (row != "_atom_site_fract_x") {
                cerr << c_message<cERROR>("PARSER::read_cif --> atom lines are not in format 'name ele x y z") << endl;
                return false;
            }
            getline(f_in,row);
            if (row != "_atom_site_fract_y") {
                cerr << c_message<cERROR>("PARSER::read_cif --> atom lines are not in format 'name ele x y z") << endl;
                return false;
            }
            getline(f_in,row);
            if (row != "_atom_site_fract_z") {
                cerr << c_message<cERROR>("PARSER::read_cif --> atom lines are not in format 'name ele x y z") << endl;
                return false;
            }
            //! spaetestens jetzt muessten die crystallographischen Daten fuer den Liganden komplett eingeparsed sein
            //! => jetzt die Umwandlungsmatrizen berechnen:
            const double to_rad = 0.0174532925199;
            double factor,a,b,c,alpha,beta,gamma;
            double cosa, cosb, cosg, sing;
            a = main_structure->ligands.back()->crysin->a; b = main_structure->ligands.back()->crysin->b;
            c = main_structure->ligands.back()->crysin->c;
            alpha = main_structure->ligands.back()->crysin->wa; beta = main_structure->ligands.back()->crysin->wb;
            gamma = main_structure->ligands.back()->crysin->wc;
            
            alpha *= to_rad; beta *= to_rad; gamma *= to_rad;
            if ((abs(alpha - beta) < 0.00001) && (abs(alpha - gamma) < 0.00001)) main_structure->ligands.back()->crysin->rhombo = true;
            else main_structure->ligands.back()->crysin->rhombo = false;
            
            cosa = cos(alpha); cosb = cos(beta); cosg = cos(gamma); sing = sin(gamma);
            factor = sqrt(1. - (cosa*cosa) - (cosb*cosb) - (cosg*cosg) + (2.*cosa*cosb*cosg) );
            
            main_structure->ligands.back()->crysin->cryst2cart[0][0] = a;
            main_structure->ligands.back()->crysin->cryst2cart[0][1] = b*cosg;
            main_structure->ligands.back()->crysin->cryst2cart[0][2] = c*cosb;
            main_structure->ligands.back()->crysin->cryst2cart[1][0] = 0.;
            main_structure->ligands.back()->crysin->cryst2cart[1][1] = b*sing;
            main_structure->ligands.back()->crysin->cryst2cart[1][2] = c*(cosa-(cosb*cosg)) / sing;
            main_structure->ligands.back()->crysin->cryst2cart[2][0] = 0.;
            main_structure->ligands.back()->crysin->cryst2cart[2][1] = 0.;
            main_structure->ligands.back()->crysin->cryst2cart[2][2] = c * factor / sing;
            
            main_structure->ligands.back()->crysin->cart2cryst[0][0] = 1./a;
            main_structure->ligands.back()->crysin->cart2cryst[0][1] = -1.*cosg/(a*sing);
            main_structure->ligands.back()->crysin->cart2cryst[0][2] = ((cosa*cosg)-cosb)/(a*sing*factor);
            main_structure->ligands.back()->crysin->cart2cryst[1][0] = 0.;
            main_structure->ligands.back()->crysin->cart2cryst[1][1] = 1./(b*sing);
            main_structure->ligands.back()->crysin->cart2cryst[1][2] = ((cosb*cosg)-cosa)/(b*sing*factor);
            main_structure->ligands.back()->crysin->cart2cryst[2][0] = 0.;
            main_structure->ligands.back()->crysin->cart2cryst[2][1] = 0.;
            main_structure->ligands.back()->crysin->cart2cryst[2][2] = sing/(c*factor);
    
            bool breaker = false;
            int cid = 1;
            while (!breaker) {
                string row2;
                getline(f_in,row2);
                if (row2 == "" || row2[0] == '_' || row2 == "#END" || row2 == "loop_") break;
                is.clear(); is.str(row2);
                //! Atom Objekt erzeugen:
                stl_ptr<ATOM> atm(new ATOM());
                atm->type = 5;
                atm->id = cid; atm->intern_id = cid; ++cid;
                is >> atm->name;
                is >> atm->element;
                
                if (!parse_hydrogens) { //kein Wasserstoff mit einlesen
                    if (atm->element == "H") {
                        atm.kill();
                        --cid;
                        continue;
                    }
                }
                
                string ca, cb, cc;
                is >> ca >> cb >> cc;
                is.clear(); is.str(ca);
                is >> atm->coord[0]; is.clear(); is.str(cb);
                is >> atm->coord[1]; is.clear(); is.str(cc);
                is >> atm->coord[2];
                
                atm->coord *= main_structure->ligands.back()->crysin->cryst2cart;
                
                atm->sybyl_type = atm->element;
                
                h_log[atm->id] = atm->intern_id;
                main_structure->ligands.back()->atoms.push_back(atm);
            }
        }
    }
    if (finished) return true;
    else return false;
}

bool PARSER::read_cif(char const* file) {
    if (verbosity > 1) cout << "reading " << file << " ..." << endl;
    f_in.open(file);
    filename = file;
    if (!f_in) {
        cerr << c_message<cERROR>("PARSER::read_cif --> could not load ") << filename << endl;
        return false;
    }
    
    ios::pos_type old_pos = 0;
    while (read_cif_molecule(old_pos)) {
        if (f_in) f_in.seekg(old_pos);
    }
    
    f_in.close();
    f_in.clear();
    if (verbosity > 1) cout << "finished reading " << file << endl;
    return true;
}

int PARSER::read_next_cif(char const* file,int n_mols) {
    multi = true;
    filename = file;
    string scompi = file;
    
    if (still_open != scompi) {
        still_open = scompi;
        f_in.open(filename);
        if (!f_in) {
            cerr << c_message<cERROR>("PARSER::read_next_cif --> could not open ") << filename << " !" << endl;
            f_in.clear();
            return -1;
        }
    }
    
    int n_parsed = 0;
    bool clearer = true;
    ios::pos_type old_pos = 0;
    
    while (read_cif_molecule(old_pos)) {
        ++n_parsed;
        if (f_in) f_in.seekg(old_pos);
        if (n_parsed == n_mols) {
            clearer = false;
            break;
        }
    }
    
    if (clearer) {
        f_in.close();
        f_in.clear();
    }
    
    return n_parsed;
}

