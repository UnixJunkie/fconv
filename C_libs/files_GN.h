
//============================================================================
// files_GN.h -*- C++ -*-; reading and writing of molecular data files
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


#ifndef __FILESGN
#define __FILESGN


#include<fstream>
#include<string>
#include<sstream>
#include<map>
#include<tr1/unordered_map>
#include"linalg_GN.hpp"
#include"structure_GN.h"
#include"stl_ptr_GN.hpp"
#include"message_GN.hpp"
#include"string_fu_GN.hpp"

#include"atom_GN.h"
#include"molecule_GN.h"
#include"protein_GN.h"
#include"structure_additional_GN.h"
#include"atom_properties_GN.h"

using namespace std;
using namespace TXT;


//==================================================================================================
//Forward-Deklarationen:
//==================================================================================================

class PARSER;


//==================================================================================================
//Deklaration der Klasse Parser:
//==================================================================================================

class PARSER {
    private:
        const char *filename;
        const char *alt_filename;
        string still_open;
        string alt_still_open;
        string dummy;
        istringstream is;
        istringstream is_tmp;
        ostringstream os;
        ifstream f_in;
        ofstream f_out;
        
        stl_ptr<STRUCTURE> main_structure; //STRUCTURE-Obj. in das geparsed wird
        stl_ptr<STRUCTURE> structure; //aktuelles STRUCTURE-Obj. (wichtig bei MODEL-Entrys)
        
        int verbosity; //Level 0 bis 2
        int n_cav_atoms;
        int n_cav_conects;
        int curr_id;
        int total_written_atoms;
        int n_conects;
        int res_add;
        int last_res_add_num;
        int glob_sdf_read;
        char last_alt_res_num;
        bool parse_hydrogens; //sollen H-Atome mit eingelesen werden
        bool parse_water; //soll Wasser mit eingelesen werden
        bool parse_ligands;
        bool parse_comments;
        
        bool multi; //soll mol2 als multimol gelesen werden
        bool first_mol2; //erstes Molekuel im multimol2?
        bool hetatoms;
        bool parse_bonds_flag;
        bool merge_ligands;

        bool pdbqt;
        
        vector<string> pdb_headers; //HEADER, COMPND, AUTHOR, JRNL
        string last_res_name; //letzter res_name (zur Bestimmung, ob neues root-atom fuer mol2)
        stl_ptr<ATOM> last_atm; //letztes bearbeitetes ATOM-Obj.
        int cur_id; //aktuelle (fortlaufende id)
        int last_res_number;
        float ftmp; //Zwischenspeichern von Fliessommazahlen
        tr1::unordered_map<int,int> h_log; // speichert zu den originalen id's die intern_id's
        tr1::unordered_map<int,int> rh_log; // zu den intern_id's die geschriebenen id's
        tr1::unordered_map<int,int> cav_log; //speichert zu den originalen id's die cavity-id's
        
        inline void compute_pdb_line(string const& row); //Zeile einer pdb-Datei auswerten
        inline void parse_atom_line(string const& row); //eine ATOM-Zeile einer pdb-Datei in ein ATOM-Objekt verwandeln
        inline void parse_hetatm_line(string const& row); //           -"-
        inline void parse_het_line(string const& row);
        inline void parse_ter_line(string const& row);
        inline void parse_model_line(string const& row);
        inline void parse_endmdl_line(string const& row);
        inline void parse_conect_line(string const& row);
    //    void parse_ssbond_line(string const& row);
    //    void parse_anisou_line(string const& row);
        
        inline void write_header_lines();
        inline void write_het_lines();
    //    void write_ssbond_lines();
        inline void write_first_model_line();
        inline void write_pdb_atom(stl_ptr<ATOM> const& at,bool cav_mode = false);
        inline void write_atoms();
        inline void write_cav_atoms(stl_ptr<CAVITY> const& cav);
        inline void write_atoms_special();
        inline void write_conect_lines();
        inline void write_cav_conect_lines();
        inline void write_master();
        inline void write_cav_master();
        
        inline void compute_mol2_line(string const& row);
        inline int compute_mol2_line_nm(string const& row,bool next_allowed,ios::pos_type old_pos,bool returner = false);
        inline void parse_mol2_atoms(bool block_mode = false);
        inline void parse_mol2_bonds(bool block_mode = false);
        inline void parse_mol2_substructure();
        inline void parse_mol2_comment(bool block_mode = false);
        inline void parse_mol2_molecule();
        inline void parse_mol2_crysin();
        
        inline void write_mol2_file(stl_ptr<LIGAND> const& ligand);
        inline void write_mol2_file(PROTEIN *protein);
        inline bool s2cryst(string const& s,stl_ptr<CRYSIN_POSITION> const& ccp);
        inline bool read_cif_molecule(ios::pos_type &old_pos);
    public:
        PARSER(STRUCTURE *st, int vb = 1, bool parse_hyd = true, bool parse_wat = true,
               bool parse_ligs = true, bool parse_stuff = false);
        ~PARSER();
        
        void clear();
        
        void set_verbosity(int vb); //Ausgabelevel aendern
        
        bool read_pdb(char const* file,bool pdbqt_mode = false); //liest pdb in das STRUCTURE-Objekt
        
        bool write_pdb(char const* file); //schreibt pdb-file fr das STRUCTURE-Objekt
        bool write_pdb_cav(char const* file,stl_ptr<CAVITY> const& cav); //schreibt pdb-file fuer das CAVITY-Objekt
        bool write_pdb_cav_all(map<int,map<int,vec3d<float> > > &flex_map,char const* file); //schreibt alle cavities in ein pdb mit
                                                                                             //verschiedenen MODEL Eintraegen
        bool write_pdb_cav_atoms_only(char const* file);
        
        bool read_mol2(char const* file,bool first_only = false,bool atoms_only = false); //liest mol2 in das STRUCTURE-Objekt
        int read_next_mol2(char const* file,int n_mols = 1,bool atoms_only = false); //liest die naechsten n_mols Molecules
        bool write_mol2(char const* file); //schreibt mol2-file fuer das STRUCTURE-Objekt
        bool write_mol2(vector<stl_ptr<LIGAND> > &ligands,char const* file);
        bool write_next_mol2(vector<stl_ptr<LIGAND> > &ligands,char const* file);
        void write_next_mol2_end(); //!File schliessen!!!
        
        bool write_mol2_ligand(stl_ptr<LIGAND> const& ligand,char const* file); //schreibt einzelnen Liganden in mol2-file
        
        bool write_mol2_protein(PROTEIN *protein,char const* file);

        bool read_dlg(char const* file,int high = 99999999, int low = 1,int read_energy = 0,
                      bool get_flex_res = false,float cut_radius = 8.,char const* pre_cav_name = "protein");
        
        bool read_dlg_for_full_pdb(char const* file,int high = 99999999,int low = 1,char const* pre_cav_name = "protein");
        
        bool read_dlg(map<int,map<int,vec3d<float> > > &flex_map,PROTEIN *protein,
                      char const* file,bool get_flex_res);
        
        bool read_sdf(char const* file,bool read_energy = false);
        int read_next_sdf(char const* file,int n_mols = 1);
        
        bool read_cif(char const* file);
        int read_next_cif(char const* file,int n_mols = 1);
};

#endif //__FILESGN

