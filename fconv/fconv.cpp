
//============================================================================
// fconv.cpp -*- C++ -*-; main file of the program fconv
//
// Copyright (C) 2007, 2008, 2009, 2010, 2011, 2012 Gerd Neudert
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
// We acknowledge the Cambridge Crystallographic Data Centre (CCDC) for
// the permission to use CSD-derived data for bond lengths and bond angles,
// as found in 'ELEMENT_DATA_GN.hpp'.
// 
// fconv is a program, intended for parsing and manipulating multiple aspects
// and properties of molecular data. Typical tasks are: Conversion and error
// correction of multiple formats such as PDB, MOL2, SDF, DLG and CIF,
// extracting ligands from PDB as MOL2, automatic or ligand-based cavity
// detection, rmsd calculation and clustering, substructure searches,
// functional alignment and structural superposition, building of crystal
// packings, adding hydrogens, calculation of various properties like the
// number of rotatable bonds, molecular weights or vdW-volumes. The atom type
// classification is based on a consistent assignment of internal atom types,
// which are more differentiated compared to e.g. Sybyl atom types. Apart
// from the predefined mapping of these types onto Sybyl types, the user is
// able to assign own mappings by providing modified template files.
//----------------------------------------------------------------------------
//
// The program is developed since 2007, but I decided to make it open source
// just in 2010 as part of the publication. I appreciate the idea behind free
// software, especially in the field of scientific research. Nevertheless,
// at the moment of publication the documentation inside the code demands
// significant improvement. Furthermore, this program is just a side product
// of my PhD thesis. While having a robust understanding of good
// programming practise nowadays, at the beginning of my PhD I was completely
// new to C++. Thus, there are many parts in the code showing an ugly style
// and not regarding preferable programming rules.
// It's also up to you, to help on improving fconv. I would appreciate, if
// people contact me and share their opinions and changes.
// I am currently using gcc 4.5.x on Ubuntu and gcc 4.5.x on Windows7 (MinGW)
// and Mac. The code compiles without any errors or warnings with options
// '-pedantic -Wall -O2' on all 3 systems. Use '-D_LINUX_OS' if compiling
// on a linux system. I also prefer to link statically.
// 
// Gerd Neudert, Marburg 27.09.2010
//----------------------------------------------------------------------------
//
//============================================================================


#include<stdlib.h>
#include<iostream>
#include<sstream>
#include<fstream>
#include<string>
#include<vector>
#include<map>
#include<tr1/unordered_map>
#include<dirent.h>
#include<sys/stat.h>
#include<unistd.h>
#include<time.h>

#include"structure_GN.h"
#include"files_GN.h"
#include"message_GN.hpp"
#include"string_fu_GN.hpp"
#include"linalg_GN.hpp"
#include"cluster_GN.hpp"
#include"atom_GN.h"
#include"molecule_GN.h"
#include"protein_GN.h"
#include"structure_additional_GN.h"
#include"atom_properties_GN.h"
#include"crystallizer_GN.h"

using namespace std;
using namespace TXT;
using namespace string_fu;

string const version_string = "1.24  (30.01.2012)";

float rmsd_dist_function(stl_ptr<LIGAND>& l1,stl_ptr<LIGAND>& l2) {
    float val = l1->get_rmsd(l2,false);
    if (val < 0.) val = l1->get_bk_rmsd(l2,false,false,false,false);
    return val;
}

float protein_sim_function(stl_ptr<STRUCTURE>& l1,stl_ptr<STRUCTURE>& l2) {
    return l1->get_shape_sim(l2);
}

float rmsd_dist_ordered(stl_ptr<LIGAND>& l1,stl_ptr<LIGAND>& l2) {
    if (l1->atoms.size() != l2->atoms.size()) {
        cerr << c_message<cERROR>("molecules with different number of atoms -> abort") << endl;
        exit(0);
    }
    float val = 0.;
    for (unsigned int i=0; i<l1->atoms.size(); ++i) {
        val += get_square_distance(l1->atoms[i]->coord,l2->atoms[i]->coord);
    }
    val /= l1->atoms.size();
    return sqrt(val);
}

string get_clust_label(stl_ptr<LIGAND> const& lig) {
    return lig->name;
}

string get_pro_label(stl_ptr<STRUCTURE> const& pro) {
    return pro->name;
}

#if defined (_LINUX_OS)
int nutzlos(const struct dirent *useless) { // filter for distinct names here
    return 1;
}

void rename_file(string &dname,string &nname) {
    string command = "mv " + dname + " " + nname;
    
    cout << command << endl;
    
    int serr = system(command.c_str());
    if (serr) cerr << "fconv::rename_file returned with an error" << endl;
}

void rename_by_dir(vector<string> &sfiles,char const* dir_name,string &cur_ext) {
    struct dirent **tree_names;
    struct stat attrib; //!neu: fuer Systeme bei denen d_type unknown ist
    int n_names;
    n_names = scandir(dir_name, &tree_names, nutzlos, alphasort);
    for (int i=0; i<n_names; ++i) {
        string ws = tree_names[i]->d_name;
        if ((ws == ".") or (ws == "..")) continue;
        //!es folgt der neue Teil mit stat:
        stat(tree_names[i]->d_name,&attrib);
        if (S_ISDIR(attrib.st_mode)) {
            string dname = dir_name;
            string new_ext;
            dname = dname + ws + "/";
            //!hier new_ext bestimmen
            bool ws_is_int = true;
            try {
                int test;
                istringstream is;
                is.str(ws);
                is >> test;
            } catch(...) {ws_is_int = false;}
            if (ws_is_int) {
                new_ext = cur_ext + ws;
            } else new_ext = cur_ext + "_" + ws;
            
            rename_by_dir(sfiles,dname.c_str(),new_ext);
        } else {
            string dname = dir_name;
            string buf = dname + ws;
            string buf1 = buf;
            string buf2 = cur_ext + get_ext(ws);
            remove_ext(buf1);
            string nname = buf1 + buf2;
            for (vector<string>::iterator fit=sfiles.begin(); fit!=sfiles.end(); ++fit) {
                if (*fit == ws) {
                    rename_file(buf,nname);
                }
            }
        }
    }
}
#endif

bool yes_no() {
    string ui = "";
    bool ask = true;
    while (ask) {
        ask = false;
        cin >> ui;
        if (ui.size() == 0) {
            cout << " Please type 'y' or 'n': ";
            ask = true;
        } else if (!(ui[0] == 'y' || ui[0] == 'n' ||
                     ui[0] == 'Y' || ui[0] == 'N')) {
            cout << " Please type 'y' or 'n': ";
            ask = true;
        }
    }
    if (ui[0] == 'y' || ui[0] == 'Y') return true;
    if (ui[0] == 'n' || ui[0] == 'N') return false;
    return false;
}

void not_supported() {
    cout << "\n";
    cout << "So you are looking for a functionality that is currently not supported by the" << "\n";
    cout << "wizard. You should have a look on fconv -h  to see if your desired function" << "\n";
    cout << "is available. There are also some combinations that are not covered by the" << "\n";
    cout << "wizard (e.g. '-pl' and '--t'). Please contact me if you have additional tasks" << "\n";
    cout << "for 'fconv'. If I think it's useful I will implement it" << "\n" << endl;
    exit(0);
}

void pdb_choice1(string &ui) {
    cout << "\n";
    cout << "What do you want to do with your pdb file(s)?" << "\n" << "\n";
    cout << "  (1)  Just rewrite it to fix errors" << "\n";
    cout << "  (2)  Remove water" << "\n";
    cout << "  (3)  Remove hydrogens" << "\n";
    cout << "  (4)  Remove all HETATMs (ligands, water and metals)" << "\n";
    cout << "  (5)  Convert to a mol2 file" << "\n";
    cout << "  (6)  Extract ligands to mol2 files" << "\n";
    cout << "  (7)  Extract water to a mol2 file" << "\n";
    cout << "  (8)  Extract binding pocket" << "\n";
    cout << "  (9)  Give out the number of free rotatable bonds for all molecules" << "\n";
    cout << "  (q)  Do something else" << "\n" << "\n";
    cout << "  Your choice: ";
    ui = "";
    bool ask = true;
    while (ask) {
        ask = false;
        cin >> ui;
        if (ui.size() == 0) {
            cout << " Please give a number from '1' to '9' or 'q': ";
            ask = true;
        } else if (!(ui[0] == '1' || ui[0] == '2' || ui[0] == '3' ||
                     ui[0] == '4' || ui[0] == '5' || ui[0] == '6' ||
                     ui[0] == '7' || ui[0] == '8' || ui[0] == '9' ||
                     ui[0] == 'q')) {
            cout << " Please give a number from '1' to '9' or 'q': ";
            ask = true;
        }
    }
}

void mol2_choice1(string &ui) {
    cout << "\n";
    cout << "What do you want to do with your mol2 file(s)?" << "\n" << "\n";
    cout << "  (1)  Just rewrite it to fix errors" << "\n";
    cout << "  (2)  Remove water" << "\n";
    cout << "  (3)  Remove hydrogens" << "\n";
    cout << "  (4)  Recalculate all atom and bond types" << "\n";
    cout << "  (5)  Rename all atoms (using serially numbered elements)" << "\n";
    cout << "  (6)  Rename molecules (using serially numbered filename without extension)" << "\n";
    cout << "  (7)  Merge mol2 files" << "\n";
    cout << "  (8)  Add molecule to a mol2 file" << "\n";
    cout << "  (9)  Extract a molecule from a mol2 file" << "\n";
    cout << "  (a)  Split a multimol2 file into single mol2 files" << "\n";
    cout << "  (b)  Build the crystal package from a mol2 file" << "\n";
    cout << "  (c)  Give out the number of free rotatable bonds for all molecules" << "\n";
    cout << "  (d)  Calculate rmsd values for all molecules compared to a given reference" << "\n";
    cout << "  (q)  Do something else" << "\n" << "\n";
    cout << "  Your choice: ";
    ui = "";
    bool ask = true;
    while (ask) {
        ask = false;
        cin >> ui;
        if (ui.size() == 0) {
            cout << " Please give a number from '1' to '9' or\n"
                             << " letter 'a' to 'd' or 'q': ";
            ask = true;
        } else if (!(ui[0] == '1' || ui[0] == '2' || ui[0] == '3' ||
                     ui[0] == '4' || ui[0] == '5' || ui[0] == '6' ||
                     ui[0] == '7' || ui[0] == '8' || ui[0] == '9' ||
                     ui[0] == 'a' || ui[0] == 'b' || ui[0] == 'c' ||
                     ui[0] == 'd' || ui[0] == 'q')) {
            cout << " Please give a number from '1' to '9' or\n"
                             << " letter 'a' to 'd' or 'q': ";
            ask = true;
        }
    }
}

void dlg_choice1(string &ui) {
    cout << "\n";
    cout << "What do you want to do with your dlg file(s)?" << "\n" << "\n";
    cout << "  (1)  Convert them to mol2 files using a reference mol2" << "\n";
    cout << "  (2)  Convert them automatically without a reference" << "\n";
    cout << "  (q)  Do something else" << "\n" << "\n";
    cout << "  Your choice: ";
    ui = "";
    bool ask = true;
    while (ask) {
        ask = false;
        cin >> ui;
        if (ui.size() == 0) {
            cout << " Please give a number from '1' to '2' or 'q': ";
            ask = true;
        } else if (!(ui[0] == '1' || ui[0] == '2' || ui[0] == 'q')) {
            cout << " Please give a number from '1' to '2' or 'q': ";
            ask = true;
        }
    }
}

void sdf_choice1(string &ui) {
    cout << "\n";
    cout << "What do you want to do with your sdf file(s)?" << "\n" << "\n";
    cout << "  (1)  Convert them to mol2 files" << "\n";
    cout << "  (q)  Do something else" << "\n" << "\n";
    cout << "  Your choice: ";
    ui = "";
    bool ask = true;
    while (ask) {
        ask = false;
        cin >> ui;
        if (ui.size() == 0) {
            cout << " Please give a number from '1' or 'q': ";
            ask = true;
        } else if (!(ui[0] == '1' || ui[0] == 'q')) {
            cout << " Please give a number from '1' or 'q': ";
            ask = true;
        }
    }
}

void inform_patrick() {
    cout << "\n";
    cout << "This special function was exclusively introduced for Dr. Patrick Pfeffer." << "\n";
    cout << "It renames all given files in all subdirectories using the names of these" << "\n";
    cout << "directories. If you are in a directory with the following subdirs:" << "\n";
    cout << "   bla/blub/0/1/0/name.mol2" << "\n";
    cout << "   bzzz/1/2/3/zzzb/name.mol2" << "\n";
    cout << "and you call this function with '*.mol2' the two example mol2 files would be" << "\n";
    cout << "renamed to  'name_bla_blub_010.mol2'  and  'name_bzzz_123_zzzb.mol2'" << "\n";
    cout << "It is possible that this function is not working on all filesystems and with\n"
         << "all kernels!" << endl;
}

bool get_multi() {
    cout << "\n";
    cout << "Do you want to process more than one file (y/n): ";
    return yes_no();
}

bool more_opt() {
    cout << "\n";
    cout << "Do you want to do more things with your file(s) (y/n): ";
    return yes_no();
}

bool to_m_type() {
    cout << "\n";
    cout << "Writing pdb files is the default for your task." << "\n";
    cout << "Do you want to write out a mol2 file instead (y/n): ";
    return yes_no();
}

void get_multi_filenames(string &ui) {
    ui = "";
    string buf;
    cout << "\n";
    cout << "Please type in the filenames seperated by a whitespace." << "\n";
    cout << "You can also use regular expressions (e.g. '*.pdb' to process all pdb files): " << endl;
    cin >> ui;
    getline(cin,buf);
    ui += buf;
}

void get_filename(string &ui) {
    ui = "";
    string buf;
    cout << "\n";
    cout << "Please type in the name of the file you want to process:" << endl;
    cin >> ui;
    getline(cin,buf);
    ui += buf;
}

void get_target_ext(string &ui) {
    ui = "";
    cout << "\n";
    cout << "You want to process multiple files." << "\n";
    cout << "Please type in the name extension for your targetfiles (e.g. if you type 'new'\n"
             << "a file with the name 'myfile.pdb' would be written as 'myfile_new.pdb'):" << "\n";
    cin >> ui;
}

void get_special_ext(string &ui) {
    ui = "";
    cout << "\n";
    cout << "Please type in the base name for your files: " << "\n";
    cin >> ui;
}

void get_target(string &ui) {
    ui = "";
    cout << "\n";
    cout << "Please type in the name of your targetfile:" << endl;
    cin >> ui;
}

void get_def(string &ui) {
    ui = "";
    cout << "\n";
    cout << "Please type in the name of your definition file:" << endl;
    cin >> ui;
}


void get_add_source(string &ui) {
    ui = "";
    cout << "\n";
    cout << "Please type in the name of the additional sourcefile that is neccesary\n"
             << "for your task:" << endl;
    cin >> ui;
}

bool full_convert() {
    cout << "\n";
    cout << "Want to convert your pdb file to a mol2 file. There are two ways to do this." << "\n";
    cout << "1.) The atom types and the connectivity can be set using a dictionary." << "\n";
    cout << "    -> Advantages    :  Very fast / always same types for same aminoacids" << "\n";
    cout << "    -> Disadvantages :  Currently no dict for RNA or DNA (so in that case you\n"
         << "                         have to use the second mode) /" << "\n";
    cout << "                        More vulnerable to errors in the pdb file / no HETATMs" << "\n";
    cout << "2.) Full conversion of the pdb using the atom typing routines." << "\n";
    cout << "    -> Advantages    :  More flexible considering modified aminoacids (e.g. \n"
         << "                         different protonation states) /" << "\n";
    cout << "                        Also the HETATMs (ligands, metals, water) will be \n"
         << "                         converted" << "\n";
    cout << "    -> Disadvantages :  Slow compared to the first mode / Same aminoacids may\n"
         << "                         be typed different, e.g. HIS (note that this can also\n"
         << "                         be an advantage) / you are not one of the cool people\n"
         << "                         if you don't use my routines :)" << "\n" << "\n";
    cout << "Do you want to use mode 2 for a full conversion (y/n): ";
    return yes_no();
}

bool overwrite() {
    cout << "\n";
    cout << "With your current options the original file(s) would be overwritten!" << "\n";
    cout << "Do you want this (y/n): ";
    return yes_no();
}

bool spec_target() {
    cout << "\n";
    cout << "Do you want to specify a name for your targetfile (y/n): ";
    return yes_no();
}

bool get_energy() {
    cout << "\n";
    cout << "Do you want to extend the molecules names by the final docked energy (y/n): ";
    return yes_no();
}

bool get_free_energy() {
    cout << "\n";
    cout << "Do you want to extend the molecules names by the estimated free energy\n"
         << "of binding (y/n): ";
    return yes_no();
}

bool get_energy2() {
    cout << "\n";
    cout << "The sdf files can be scanned for energy or score entries." << "\n";
    cout << "Do you want to extend the molecules names by the corresponding values (y/n): ";
    return yes_no();
}

bool special_names() {
    cout << "\n";
    cout << "Splitting a multimol2 'my.mol2' would normally result in 'my_0.mol2',\n"
         <<" 'my_1.mol2',..." << "\n";
    cout << "Do you want to specify another base name than the original filename (y/n): ";
    return yes_no();
}

bool special_names2() {
    cout << "\n";
    cout << "Normally 'my.mol2' would become 'my_cryst.mol2'" << "\n";
    cout << "Do you want to specify another base name than the original filename (y/n): ";
    return yes_no();
}

bool special_names3() {
    cout << "\n";
    cout << "Do you want to use a special extension for your filenames (y/n): ";
    return yes_no();
}

void get_at_mode(string &ui) {
    cout << "\n";
    cout << "Now you have to make a choice between the differnt atom typing modes." << "\n";
    cout << "Regardless of the mode all atoms will be set to the internal types initially.\n"
         << "The different modes only specify to which types the internal types will be set\n"
         << "in the output file(s). You can also use your own definition file with" << "\n";
    cout << "rules how to combine the internal types to your own. This wizard can also\n"
         << "teach you how to get an example def file." << "\n";
    cout << "  (1)  Use the internal atom types" << "\n";
    cout << "  (2)  Use the sybyl atom types (note that this won't result in an option as\n"
         << "        it is default)" << "\n";
    cout << "  (3)  Use the modified sybyl atom types" << "\n";
    cout << "  (4)  Use your own mode (so you will need a def file)" << "\n";
    cout << "  (q)  Do something else" << "\n" << "\n";
    cout << "  Your choice: ";
    ui = "";
    bool ask = true;
    while (ask) {
        ask = false;
        cin >> ui;
        if (ui.size() == 0) {
            cout << " Please give a number from '1' to '4' or 'q': ";
            ask = true;
        } else if (!(ui[0] == '1' || ui[0] == '2' || ui[0] == '3' || ui[0] == '4' || ui[0] == 'q')) {
            cout << " Please give a number from '1' to '4' or 'q': ";
            ask = true;
        }
    }
}

void inform_at_mode(string &ui) {
    cout << "\n";
    cout << "You want to get informations about the atom typing modes or you want to get a\n"
         << "template for your own def files. Regardless of the mode all atoms will be set\n"
         << "to the internal types initially. The different modes only specify to which\n"
         << "types the internal types will be set in the output file(s). You can also use\n"
         << "your own definition file with rules how to combine the internal types to your\n"
         << "own types. What do you want:" << "\n";
    cout << "  (1)  Write a definition file for the internal atom types" << "\n";
    cout << "  (2)  Write a definition file for the sybyl atom types" << "\n";
    cout << "  (3)  Write a definition file for the modified sybyl atom types" << "\n";
    cout << "  (q)  Do something else" << "\n" << "\n";
    cout << "  Your choice: ";
    ui = "";
    bool ask = true;
    while (ask) {
        ask = false;
        cin >> ui;
        if (ui.size() == 0) {
            cout << " Please give a number from '1' to '3' or 'q': ";
            ask = true;
        } else if (!(ui[0] == '1' || ui[0] == '2' || ui[0] == '3' || ui[0] == 'q')) {
            cout << " Please give a number from '1' to '3' or 'q': ";
            ask = true;
        }
    }
}

void get_extract_mode(string &ui) {
    cout << "\n";
    cout << "You want to cut out the binding pocket of a protein. Well there are different\n"
         << "modes to do this and there are also some additional options. First you have\n"
         << "the following choice:" << "\n";
    cout << "  (1)  Extract a pocket for each ligand in your sourcefile(s)" << "\n";
    cout << "  (2)  Extract a pocket for each molecule from a seperate file (mol2 or pdb)" << "\n";
    cout << "  (3)  Try to detect possible binding pockets automatically" << "\n";
    cout << "  (q)  Do something else" << "\n" << "\n";
    cout << "  Your choice: ";
    ui = "";
    bool ask = true;
    while (ask) {
        ask = false;
        cin >> ui;
        if (ui.size() == 0) {
            cout << " Please give a number from '1' to '3' or 'q': ";
            ask = true;
        } else if (!(ui[0] == '1' || ui[0] == '2' || ui[0] == '3' || ui[0] == 'q')) {
            cout << " Please give a number from '1' to '3' or 'q': ";
            ask = true;
        }
    }
}

void get_build_mode(string &ui) {
    cout << "\n";
    cout << "So you are intersted in a crystal package and you hopefully have sourcefiles\n"
         << "with a CRYSIN entry. You have three different options:" << "\n";
    cout << "  (1)  Only create the unit cell" << "\n";
    cout << "  (2)  Create the unit cell and duplicate it n times in all directions" << "\n";
    cout << "  (3)  Create the unit cell and duplicate it within a given cutoff radius" << "\n";
    cout << "  (q)  Do something else" << "\n" << "\n";
    cout << "  Your choice: ";
    ui = "";
    bool ask = true;
    while (ask) {
        ask = false;
        cin >> ui;
        if (ui.size() == 0) {
            cout << " Please give a number from '1' to '3' or 'q': ";
            ask = true;
        } else if (!(ui[0] == '1' || ui[0] == '2' || ui[0] == '3' || ui[0] == 'q')) {
            cout << " Please give a number from '1' to '3' or 'q': ";
            ask = true;
        }
    }
}

void get_cut_radius(string &ui) {
    cout << "\n";
    cout << "Now you can specify a radius (in Angstrom) for your task. If you want to \n"
         << "extract binding pockets this radius is used as cutoff. So all aminoacids, \n"
         << "having a distance to your ligand equal or shorter than this value, are\n"
         << "considered as binding pocket. In automatic pocket detection mode it's the \n"
         << "cutoff distance to the determined cavity grid points." << "\n";
    cout << "If you are building crystals you can specify here the radius for duplicating \n"
         << "the unit cell. However, the default value here is 8.0 . Please type in another\n"
         << "float or type 'd' for default: ";
    ui = "";
    cin >> ui;
}

void how_many_times(string &ui) {
    cout << "\n";
    cout << "If you want to extract binding pockets automatically you can specify here the" << "\n";
    cout << "number of pockets to write out. If you type in '2' for example the two 'best'\n"
         << "pockets will be written. The default value here is 1. So only the best pocket\n"
         << "is written by default." << "\n";
    cout << "If you are building crystals you can specify here how often the unit cell \n"
         << "should be duplicated in all directions. The default here is 0." << "\n";
    cout << "Please type in another integer or type 'd' for default: ";
    ui = "";
    cin >> ui;
}

void get_min_buried(string &ui) {
    cout << "\n";
    cout << "If the automatic detection does not perform very well for your target it may\n"
         << "help to change the minimum buriedness for all cavity surface points. \n"
         << "Increasing this number results in smaller pockets." << "\n";
    cout << "Note that the result of the automatic detection is also affected by the \n"
         << "initial orientation of the protein as the algorithm is grid based. This will\n"
         << "be changed in one of the next versions." << "\n";
    cout << "Things may come clear if you are using the debug mode to visualize the \n"
         << "pocket surface points. Therefor you have to type '---d' at the end of your\n"
         << "commandline. This results in a python file for each pocket that you can\n"
         << "run from pymol." << "\n";
    cout << "The default value here is 18" << "\n";
    cout << "Please type in another float or type 'd' for default: ";
    ui = "";
    cin >> ui;
}

void inform_about_ra() {
    cout << "\n";
    cout << "OK - your atoms will be renamed by their element type. The first O will be\n"
         << "named 'O1' the second 'O2'..." << "\n";
    cout << endl;
}

void inform_about_rm() {
    cout << "\n";
    cout << "OK - your molecules will be renamed using the filename. If you have a file\n"
         << "'my.mol2' with two maolecules, these molecules will be named 'my_0' and 'my_1'" << "\n";
    cout << endl;
}

void get_extract_number(string &ui) {
    cout << "\n";
    cout << "Please type in the number of the molecule you want to extract (the first\n"
         << "molecule is 0): ";
    ui = "";
    cin >> ui;
}

void start_wizard() {
    bool ask;
    bool multi = false;
    string user_input;
    ostringstream command;
    ostringstream add_options2;
    command << "fconv ";
    cout << "                                               ?" << "\n";
    cout << "                                             _/_\\_" << "\n";
    cout << "                                             (o o)" << "\n";
    cout << "+----------------------------------------ooO--(_)--Ooo---------+" << "\n";
    cout << "|          Welcome to the fconv commandline wizard !!!         |" << "\n";
    cout << "|                                                              |" << "\n";
    cout << "| -> This tool does NOT perform any operations on your files!  |" << "\n";
    cout << "|    It's just an interactive help that tries to find out      |" << "\n";
    cout << "|    what you like to do with fconv and in the end it will     |" << "\n";
    cout << "|    give you the corresponding commandline for calling fconv. |" << "\n";
    cout << "| -> There are still some functions you won't get using this   |" << "\n";
    cout << "|    tool (The wizard is based on vers. 0.34), so it's just    |" << "\n";
    cout << "|    meant as an introduction to get a feeling for the syntax. |" << "\n";
    cout << "| -> You can leave this tool anytime by pressing  <STRG>+C     |" << "\n";
    cout << "+--------------------------------------------------------------+" << "\n";
    cout << "\n" << "\n";
    cout << "What do you want to do?" << "\n" << "\n";
    cout << "  (1)  Do something with pdb files" << "\n";
    cout << "  (2)  Do something with mol2 files" << "\n";
    cout << "  (3)  Do something with dlg files" << "\n";
    cout << "  (4)  Do something with sdf files" << "\n";
    cout << "  (5)  Get an atom types definition file" << "\n";
    cout << "  (6)  Rename files by extending the names with the superior directory names" << "\n";
    cout << "  (q)  Do something else" << "\n" << "\n";
    cout << "  Your choice: ";
    user_input = "";
    ask = true;
    while (ask) {
        ask = false;
        cin >> user_input;
        if (user_input.size() == 0) {
            cout << " Please give a number from '1' to '6' or 'q': ";
            ask = true;
        } else if (!(user_input[0] == '1' || user_input[0] == '2' || user_input[0] == '3' ||
                     user_input[0] == '4' || user_input[0] == '5' || user_input[0] == '6' ||
                     user_input[0] == 'q')) {
            cout << " Please give a number from '1' to '6' or 'q': ";
            ask = true;
        }
    }
    cout << "\n";
    
    //!==============================================================================================
    if (user_input[0] == '1') { //!pdb wizard
        ask = true;
        bool ask_for_target = false;
        bool special_target = false;
        while (ask) {
            ask = false;
            pdb_choice1(user_input);
            if (user_input[0] == '1') { //pdb file neu schreiben
                command << "-R "; ask_for_target = true;
            } else if (user_input[0] == '2') { //wasser entfernen
                command << "-W "; ask_for_target = true;
                if (more_opt()) ask = true;
            } else if (user_input[0] == '3') { //wasserstoff entfernen
                command << "-H "; ask_for_target = true;
                if (more_opt()) ask = true;
            } else if (user_input[0] == '4') { //alle HETATMs entfernen
                command << "-L "; ask_for_target = true;
                if (more_opt()) ask = true;
            } else if (user_input[0] == '5') { //convert to mol2
                if (full_convert()) {
                    command << "-f ";
                    get_at_mode(user_input);
                    if (user_input[0] == '1') command << "--m=0 ";
                    else if (user_input[0] == '3') command << "--m=2 ";
                    else if (user_input[0] == '4') {
                        get_def(user_input);
                        add_options2 << "--d=" << user_input << " ";
                    } else if (user_input[0] == 'q') not_supported();
                } else {
                    command << "-c ";
                    get_at_mode(user_input);
                    if (user_input[0] == '1') command << "--m=0 ";
                    else if (user_input[0] == '3') command << "--m=2 ";
                    else if (user_input[0] == '4') {
                        get_def(user_input);
                        add_options2 << "--d=" << user_input << " ";
                    } else if (user_input[0] == 'q') not_supported();
                }
                if (more_opt()) ask = true;
            } else if (user_input[0] == '6') { //Liganden nach mol2 extrahieren
                command << "-l ";
                get_at_mode(user_input);
                if (user_input[0] == '1') command << "--m=0 ";
                else if (user_input[0] == '3') command << "--m=2 ";
                else if (user_input[0] == '4') {
                    get_def(user_input);
                    add_options2 << "--d=" << user_input << " ";
                } else if (user_input[0] == 'q') not_supported();
            } else if (user_input[0] == '7') { //Wasser in mol2 file schreiben
                command << "-w ";
                get_at_mode(user_input);
                if (user_input[0] == '1') command << "--m=0 ";
                else if (user_input[0] == '3') command << "--m=2 ";
                else if (user_input[0] == '4') {
                    get_def(user_input);
                    add_options2 << "--d=" << user_input << " ";
                } else if (user_input[0] == 'q') not_supported();
            } else if (user_input[0] == '8') { //Bindetasche extrahieren
                get_extract_mode(user_input);
                if (user_input[0] == '1') {
                    command << "-pl ";
                } else if (user_input[0] == '2') {
                    command << "-pl ";
                    get_add_source(user_input);
                    add_options2 << "--s=" << user_input << " ";
                } else if (user_input[0] == '3') {
                    command << "-pa ";
                    how_many_times(user_input);
                    if (user_input[0] == '0' || user_input[0] == '1' || user_input[0] == '2' ||
                        user_input[0] == '3' || user_input[0] == '4' || user_input[0] == '5' ||
                        user_input[0] == '6' || user_input[0] == '7' || user_input[0] == '8' ||
                        user_input[0] == '9') add_options2 << "--n=" << user_input << " ";
                    get_min_buried(user_input);
                    if (user_input[0] == '0' || user_input[0] == '1' || user_input[0] == '2' ||
                        user_input[0] == '3' || user_input[0] == '4' || user_input[0] == '5' ||
                        user_input[0] == '6' || user_input[0] == '7' || user_input[0] == '8' ||
                        user_input[0] == '9') add_options2 << "--v=" << user_input << " ";
                } else if (user_input[0] == 'q') not_supported();
                special_target = true;
                get_cut_radius(user_input);
                if (user_input[0] == '0' || user_input[0] == '1' || user_input[0] == '2' ||
                    user_input[0] == '3' || user_input[0] == '4' || user_input[0] == '5' ||
                    user_input[0] == '6' || user_input[0] == '7' || user_input[0] == '8' ||
                    user_input[0] == '9') add_options2 << "--r=" << user_input << " ";
                if (to_m_type()) {
                    add_options2 << "--f=m ";
                    get_at_mode(user_input);
                    if (user_input[0] == '1') command << "--m=0 ";
                    else if (user_input[0] == '3') command << "--m=2 ";
                    else if (user_input[0] == '4') {
                        get_def(user_input);
                        add_options2 << "--d=" << user_input << " ";
                    } else if (user_input[0] == 'q') not_supported();
                }
                if (more_opt()) ask = true;
            } else if (user_input[0] == '9') { //rotierbare Bindungen ausgeben
                command << "-gfr ";
                special_target = true;
            } else if (user_input[0] == 'q') not_supported();
        }
        multi = get_multi();
        if (multi) get_multi_filenames(user_input);
        else get_filename(user_input);
        command << user_input << " ";
        if (!special_target) {
            if (ask_for_target) {
                if (overwrite()) ask_for_target = false;
            } else if (spec_target()) ask_for_target = true;
            if (ask_for_target) {
                if (multi) get_target_ext(user_input);
                else get_target(user_input);
                command << "--t=" << user_input << " ";
            }
        }
    
    //!==============================================================================================
    } else if (user_input[0] == '2') { //!mol2 wizard
        ask = true;
        bool ask_for_target = false;
        bool special_target = false;
        bool ask_auto = false;
        bool cry = false;
        bool no_target = false;
        while (ask) {
            ask = false;
            mol2_choice1(user_input);
            if (user_input[0] == '1') { //pdb file neu schreiben
                command << "-R "; ask_for_target = true;
            } else if (user_input[0] == '2') { //wasser entfernen
                command << "-W "; ask_for_target = true;
                if (more_opt()) ask = true;
            } else if (user_input[0] == '3') { //wasserstoff entfernen
                command << "-H "; ask_for_target = true;
                if (more_opt()) ask = true;
            } else if (user_input[0] == '4') { //atom typing
                command << "-T "; ask_for_target = true;
                get_at_mode(user_input);
                if (user_input[0] == '1') command << "--m=0 ";
                else if (user_input[0] == '3') command << "--m=2 ";
                else if (user_input[0] == '4') {
                    get_def(user_input);
                    add_options2 << "--d=" << user_input << " ";
                } else if (user_input[0] == 'q') not_supported();
                if (more_opt()) ask = true;
            } else if (user_input[0] == '5') { //Atome umbenennen
                command << "-ra "; ask_for_target = true;
                inform_about_ra();
                if (more_opt()) ask = true;
            } else if (user_input[0] == '6') { //Molekuele umbenennen
                command << "-rm "; ask_for_target = true;
                inform_about_rm();
                if (more_opt()) ask = true;
            } else if (user_input[0] == '7') { //merge mol2 files 
                command << "-m "; special_target = true;
            } else if (user_input[0] == '8') { //add mol2 file
                command << "-a "; special_target = true;
                if (more_opt()) ask = true;
            } else if (user_input[0] == '9') { //extract mol2 molecule 
                get_extract_number(user_input);
                command << "--e=" << user_input << " ";
                if (more_opt()) ask = true;
            } else if (user_input[0] == 'a') { //split multimol2 
                command << "-s ";
                ask_auto = true;
                if (more_opt()) ask = true;
            } else if (user_input[0] == 'b') { //Kristalle bauen
                command << "-C ";
                ask_auto = true;
                cry = true;
                get_build_mode(user_input);
                if (user_input[0] == '2') {
                    how_many_times(user_input);
                    if (user_input[0] == '0' || user_input[0] == '1' || user_input[0] == '2' ||
                        user_input[0] == '3' || user_input[0] == '4' || user_input[0] == '5' ||
                        user_input[0] == '6' || user_input[0] == '7' || user_input[0] == '8' ||
                        user_input[0] == '9') add_options2 << "--n=" << user_input << " ";
                } else if (user_input[0] == '3') {
                    get_cut_radius(user_input);
                    if (user_input[0] == '0' || user_input[0] == '1' || user_input[0] == '2' ||
                        user_input[0] == '3' || user_input[0] == '4' || user_input[0] == '5' ||
                        user_input[0] == '6' || user_input[0] == '7' || user_input[0] == '8' ||
                        user_input[0] == '9') add_options2 << "--r=" << user_input << " ";
                } else if (user_input[0] == 'q') not_supported();
            } else if (user_input[0] == 'c') { //freerot bonds 
                command << "-gfr ";
                special_target = true;
                no_target = true;
            } else if (user_input[0] == 'd') { //rmsd berechnen
                command << "-rmsd ";
                get_add_source(user_input);
                add_options2 << "--s=" << user_input << " ";
                special_target = true;
                no_target = true;
            } else if (user_input[0] == 'q') not_supported();
        }
        if (!special_target) {
            multi = get_multi();
            if (multi) get_multi_filenames(user_input);
            else get_filename(user_input);
        } else get_multi_filenames(user_input);
        command << user_input << " ";
        if (!special_target) {
            if (!ask_auto) {
                if (ask_for_target) {
                    if (overwrite()) ask_for_target = false;
                } else if (spec_target()) ask_for_target = true;
                if (ask_for_target) {
                    if (multi) get_target_ext(user_input);
                    else get_target(user_input);
                    command << "--t=" << user_input << " ";
                }
            } else {
                if (!cry) {
                    if (special_names()) get_special_ext(user_input);
                    command << "--t=" << user_input << " ";
                } else {
                    if (special_names2()) get_special_ext(user_input);
                    command << "--t=" << user_input << " ";
                }
            }
        } else {
            if (!no_target) {
                get_target(user_input);
                command << "--t=" << user_input << " ";
            }
        }
    
    //!==============================================================================================
    } else if (user_input[0] == '3') { //!dlg wizard
        dlg_choice1(user_input);
        if (user_input[0] == '1') {
            command << "-c ";
            get_add_source(user_input);
            add_options2 << "--s=" << user_input << " ";
        } else if (user_input[0] == '2') {
            command << "-a ";
            get_at_mode(user_input);
            if (user_input[0] == '1') command << "--m=0 ";
            else if (user_input[0] == '3') command << "--m=2 ";
            else if (user_input[0] == '4') {
                get_def(user_input);
                add_options2 << "--d=" << user_input << " ";
            } else if (user_input[0] == 'q') not_supported();
        } else if (user_input[0] == 'q') not_supported();
        if (get_energy()) command << "-E ";
        if (get_free_energy()) command << "-EF ";
        multi = get_multi();
        if (multi) get_multi_filenames(user_input);
        else get_filename(user_input);
        if (multi) {if (special_names3()) get_special_ext(user_input);}
        else get_target(user_input);
        command << "--t=" << user_input << " ";
    
    //!==============================================================================================
    } else if (user_input[0] == '4') { //!sdf wizard
        sdf_choice1(user_input);
        if (user_input[0] == '1') {
            command << "-c ";
            get_at_mode(user_input);
            if (user_input[0] == '1') command << "--m=0 ";
            else if (user_input[0] == '3') command << "--m=2 ";
            else if (user_input[0] == '4') {
                get_def(user_input);
                add_options2 << "--d=" << user_input << " ";
            } else if (user_input[0] == 'q') not_supported();
        } else if (user_input[0] == 'q') not_supported();
        if (get_energy2()) command << "-E ";
        multi = get_multi();
        if (multi) get_multi_filenames(user_input);
        else get_filename(user_input);
        if (special_names3()) get_special_ext(user_input);
        command << "--t=" << user_input << " ";
    
    //!==============================================================================================
    } else if (user_input[0] == '5') { //!atom typing
        inform_at_mode(user_input);
        if (user_input[0] == '1') {
            command << "--M=0";
        } else if (user_input[0] == '2') {
            command << "--M=1";
        } else if (user_input[0] == '3') {
            command << "--M=2";
        } else if (user_input[0] == 'q') not_supported();
    
    //!==============================================================================================
    } else if (user_input[0] == '6') { //!patrick :)
        inform_patrick();
        command << "-rf ";
        get_multi_filenames(user_input);
        command << user_input;
    } else if (user_input[0] == 'q') not_supported();
    
    command << add_options2.str();
    cout << "\n";
    cout << "                                               !" << "\n";
    cout << "                                             _/_\\_" << "\n";
    cout << "                                             (o o)" << "\n";
    cout << "-----------------------------------------ooO--(_)--Ooo------------------------" << "\n";
    cout << "Your commandline is:\n" << command.str() << endl << endl;
}

void show_example() {
    cout << "\n";
    cout << "Some examples using 'fconv' :" << "\n" << "\n";
    cout << "Examples for pdb files:" << "\n";
    cout << "=======================" << "\n";
    cout << "'fconv -L *.pdb' :" << "\n";
    cout << "    This command converts all pdb files in your current directory to clean\n"
         << "    pdb files without any HETATMs (or peptides). The original files will\n"
         << "    be overwritten!!!" << "\n" << "\n";
    cout << "'fconv -L my.pdb --t=myclean.pdb' :" << "\n";
    cout << "    Here the clean pdb is written to a new file specified by '--t='" << "\n" << "\n";
    cout << "'fconv -c *.pdb --t=blablub' :" << "\n";
    cout << "    Converts all pdb files from the current directory to mol2 files using\n"
         << "    dictionaries for the atom typing. The targetfiles will be called\n"
         << "    'filename_blablub.mol2'. (without the --t option it would be \n"
         << "    'filename.mol2' --- so there is no need for --t)" << "\n";
    cout << "'fconv -l *.pdb' :" << "\n";
    cout << "    Extracts all ligands of all pdb files in your current directory to mol2\n"
         << "    files using sybyl atom types." << "\n" << "\n";
    cout << "'fconv --m=0 -l my.pdb' :" << "\n";
    cout << "    Extracts all ligands from 'my.pdb' to seperate mol2 files using\n"
         << "    internal atom types." << "\n" << "\n";
    cout << "'fconv --pl=2 my.pdb' :" << "\n";
    cout << "    Gives you the 8A pocket around ligand number 3 (the ligand count\n"
         << "    starts with 0). The target files name will be 'my_ligname_poc.pdb'." << "\n" << "\n";
    cout << "'fconv --pl=0 my.pdb --t=mypocket.pdb --r=8' :" << "\n";
    cout << "    Gives you the 8A pocket around the first ligand and writes it\n"
         << "    to 'mypocket.pdb'." << "\n" << "\n";
    cout << "'fconv -pl *.pdb --r=8 --f=m' :" << "\n";
    cout << "    For all ligands of all pdb files in your current directory a 8A pocket\n"
         << "    is written as mol2 file" << "\n" << "\n";
    cout << "'fconv -pa *.pdb --n=2' :" << "\n";
    cout << "    For all pdb files in your current directory the two 'best' 6A pockets\n"
         << "    are calculated and written as pdb file" << "\n" << "\n";
    cout << "'fconv -W -f -ha my.pdb' :" << "\n";
    cout << "    All water molecules will be removed (-W), the pdb will be converted to\n"
         << "    a mol2 using the full atom typing" << "\n";
    cout << "    routines and hydrogens will be added to the molecules." << "\n" << "\n";
    
    cout << "Examples for mol2 files:" << "\n";
    cout << "========================" << "\n";
    cout << "'fconv -T *.mol2 --t=new' :" << "\n";
    cout << "    Recalculates all atom types for all mol2 files in the current directory" << "\n";
    cout << "    and writes them to 'filename_new.mol2'" << "\n" << "\n";
    cout << "'fconv -T my.mol2 --t=trallala.mol2 --d=my_own_types.def' :" << "\n";
    cout << "    Recalculates all atom types for 'my.mol2' using the atom types defined\n"
         << "    in 'my_own_types.def' and writes them to 'trallala.mol2'" << "\n" << "\n";
    cout << "'fconv -ha my.mol2' :" << "\n";
    cout << "    This will add missing hydrogens to the molecule(s)" << "\n" << "\n";
    cout << "'fconv -rmsd my.mol2 --s=reference.mol2' :" << "\n";
    cout << "    Rmsd values for all molecules in my.mol2 will be calculated with respect\n"
         << "    to the given reference." << "\n" << "\n";
    cout << "'fconv -m *.mol2 --t=all_mols.mol2' :" << "\n";
    cout << "    All mol2 files in your current directory are merged to 'all_mols.mol2'" << "\n" << "\n";
    cout << "'fconv -s multi.mol2' :" << "\n";
    cout << "    The multimol is split to single mol2 files named 'multi_0.mol2' to\n"
         << "    'multi_n.mol2'" << "\n" << "\n";
    cout << "Examples for dlg files:" << "\n";
    cout << "=======================" << "\n";
    cout << "'fconv -c --n=1-2 mydock.dlg --s=origin.mol2 --t=results.mol2' :" << "\n";
    cout << "    Converts the first two solutions from your autodock dlg file to a\n"
         << "    multimol2 with the name 'results.mol2'. For this conversion you need\n"
         << "    the original mol2 file of your docked molecule." << "\n" << "\n";
    cout << "'fconv -a mydock.dlg' :" << "\n";
    cout << "    Converts mydock.dlg to mydock.mol2 using atom typing routines" << "\n" << "\n";
    cout << "'fconv -f mydock.dlg --s2=used_for_docking.pdb' :" << "\n";
    cout << "    This is for dlg files from a docking run with flexible residues.\n"
         << "    Additionally to each docking solution the corresponding binding pocket\n"
         << "    will be written." << "\n" << "\n";
    cout << "Examples for sdf files:" << "\n";
    cout << "=======================" << "\n";
    cout << "'fconv -c my.sdf' : " << "\n";
    cout << "    Converts the file my.sdf to my.mol2" << "\n" << "\n";
    cout << "'fconv -c --m==2 my.sdf --t=my_new.mol2' :" << "\n";
    cout << "    Converts the file my.sdf to my_new.mol2 using modified sybyl types" << "\n" << endl;
}

void show_header() {
    cout << "\n";
    cout << " +---------------------------------------------------------------------------+" << "\n";
    cout << " | 'fconv'         conversion and manipulation of some standard file formats |" << "\n";
    cout << " |  author     :   Gerd Neudert                                              |" << "\n";
    cout << " |  supervisor :   Prof. Dr. G. Klebe                          " << c_message<cGREY>("___    _ _") << "    |" << "\n";
    cout << " |  mailto     :   neudert@staff.uni-marburg.de                " << c_message<cGREY>("))_    )\\`)") << "   |" << "\n";
    cout << " |  version    :   " << version_string << "                         " << c_message<cGREY>("((_( o ((\\( o") << "  |" << "\n";
    //cout << " |                 development version                                       |" << "\n";
    cout << " |                 public version (check for new versions on www.agklebe.de) |" << "\n";
    cout << " +---------------------------------------------------------------------------+" << "\n" << "\n";
    cout << "  Copyright (C) 2007, 2008, 2009, 2010, 2011, 2012 Gerd Neudert\n";
    cout << "  This is free software; see the source for copying conditions. There is NO\n"
         << "  warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n\n";
    cout << "=> Usage:    fconv  Options1  Files  [Options2]" << "\n";
    cout << "=================================================" << "\n" << "\n";
}

void show_help() {
    show_header();
    cout << " Type 'fconv -h' to get detailed help  or  'fconv -eg' to get some examples" << "\n";
    cout << " Use '-hp', '-hm', '-hd', '-hs' or '-hc' to get help only on options\n";
    cout << " for pdb (pdbqt), mol2, dlg, sdf or cif files\n";
    cout << " If you still don't know what to do try out the wizard with 'fconv -wiz'" << "\n" << endl;
}

void show_gen_opt() {
    cout << "  Files    :  Valid file types are pdb, pdbqt, mol2, dlg, sdf and cif\n";
    cout << "------------" << "\n";
    cout << "  Options1 :" << "\n";
    cout << "------------" << "\n";
    cout << "    -h[s]  :  Show help (e.g. s='p', for help on pdb only)" << "\n";
    cout << "    -eg    :  Gives some examples" << "\n";
    cout << "    -wiz   :  Start an interactive mode to determine the commandline \n"
         << "              for your specific task" << "\n";
    cout << "   --m=n   :  Specify the atom typing mode. Valid modes are:" << "\n";
    cout << "              0: internal atom types /   1: sybyl atom types (DEFAULT) /  \n"
         << "              2: modified sybyl types" << "\n";
    cout << "              You can also use your own definitions for the atom typing.\n"
         << "              Therefor you have to supply a definition file (see --d under\n"
         << "              options2)" << "\n";
    cout << "   --M=n   :  Write out a definition file for the atom typing mode specified\n"
         << "              by n (use this as a template for your own definition files)" << "\n";
    cout << "   --p=b   :  Specify standard protonation state for some groups by\n"
         << "              a 6-digit binary b. The default is 111110\n";
    cout << "              '1' stands for charged (e.g. deprotonated acid or protonated\n"
         << "              guanidine) and '0' for uncharged.\n";
    cout << "               digit 1: acids\n";
    cout << "               digit 2: HxSOy and SO2NH\n";
    cout << "               digit 3: HxPOy\n";
    cout << "               digit 4: guanidines\n";
    cout << "               digit 5: amidines\n";
    cout << "               digit 6: amines\n";
    cout << "    -KA    :  Kekulize aromatics (use alternating bond types '1' and '2'\n"
         << "              instead of type 'ar')\n";
//    cout << "    -kc    :  Kekulize charged: also use '2' and '1' in case of charged\n"
//         << "              groups like COO- or SO2NH- (instead of 'ar')\n";
    cout << "    -NCA   :  Allow no charged aromatics\n";
    cout << "   --v=n   :  Set verbosity (default: n=1)  |  'n' must be 0 (just errors),\n"
         << "              1 (warnings and workflow) or 2 (info)" << "\n";
    cout << "              (Modes 0 and 2 will be improved in one of the next versions.)" << "\n";
    cout << "    -rf    :  Renames all given files from all subdirectories using the names\n"
         << "              of those subdirs. '-rf name.ext' renames the file 'name.ext' in\n"
         << "              the subdir '/type/0/1/' to 'name_type_01.ext'" << "\n";
    cout << "\n";
}

void show_pdb_opt() {
    cout << "-> Valid options for pdb (or pdbqt) input files:" << "\n";
    cout << "    -c     :  Convert protein to mol2 (without ligands)" << "\n";
    cout << "              Currently there is no dict for DNA and RNA (try '-f' until it is)" << "\n";
    cout << "    -f     :  Full conversion of the whole pdb file to a mol2 file" << "\n";
    cout << "              With that option the atom typing routines are also used for the\n"
         << "              protein (this will take some time)" << "\n";
    cout << "    -R     :  Just rewrite the pdb (this may fix minor errors)\n";
    cout << "              pdbqt files will be rewritten as pdb files\n";
    //cout << "    -F     :  As -R but this time trying to fix also major errors and to\n"
    //     << "              write out a complete pdb format conform file." << "\n";
    cout << "    -W     :  Remove water" << "\n";
    cout << "    -H     :  Remove hydrogens" << "\n";
    cout << "    -V     :  Parse REMARKS and other header stuff" << "\n";
    cout << "    -ha    :  Add hydrogens (using standard protonation states)" << "\n";
    cout << "              This is currently only possible in  -f  mode!" << "\n";
    cout << "    -L     :  Remove ligands and other HETATMs (also water)" << "\n";
    cout << "    -LL    :  As -L, but also removes peptidic ligands" << "\n";
    cout << "    -l     :  Extract ligands to seperate mol2 files" << "\n";
    cout << "              Use --s to specify a ligand by its res_name\n";
    cout << "              Only ligands with 6 or more heavy atoms will be considered.\n"
         << "              You can change this number with --n from options2.\n";
//    cout << "    -lm    :  Extract ligands to multimol2 file" << "\n";
    cout << "    -oa    :  Aligns all input structures on a reference structure specified\n"
         << "              by --s=ref.pdb. Please note that the alignment is based on\n"
         << "              C-alpha atoms only. You can also specify a mol2 or PDB with\n"
         << "              --s2=to_align_file in order to align all structures from" << "\n";
    cout << "              this file instead of the pdb. (this also works for -oa2)" << "\n";
    cout << "              Use '-fr' as additional option if the spatial alignment should\n"
         << "              be based on all atoms!" << "\n";
    cout << "    -oa2   :  As '-oa', but this time the spatial alignment is NOT based on a\n"
         << "              sequence alignment. Instead of this, a backbone-shape-based matching\n"
         << "              is used. You can change the tolerance (max_displacement) with --r.\n"
         << "              The default max_displacement is 1.4A and the value should be always\n"
         << "              in a range between 0.5 and 2.0A.\n";
    cout << "    -clust :  Hierarchical complete linkage clustering with respect to backbone.\n"
         << "              shape similarity (see -oa2).\n";
    cout << "              You will get the corresponding dendrogram and similarity matrix\n"
         << "              as svg files.\n";
//    cout << "    -poa   :  Aligns two binding pockets." << "\n";
    cout << "    -pl    :  Extract binding pockets around all ligands (gives a mol2 or pdb\n"
         << "              file for each ligand). You can also supply a seperate ligand\n"
         << "              file with --s=ligandfile (pdb or mol2). In that case only the\n"
         << "              pocket for that ligand (or the ligands if you supply a multimol)\n"
         << "              will be extracted." << "\n";
    cout << "    -pp    :  The same as -pl, but only atoms within the given radius are\n"
         << "              written instead of whole residues." << "\n";
//    cout << "    -plu   :  Extract one binding pocket around all ligands" << "\n";
//    cout << "   --pl=s  :  Extract binding pocket around ligand, where 's' is the name of that ligand" << "\n";
    cout << "   --pl=n  :  Extract binding pocket around ligand, where 'n' is the number\n"
         << "              of that ligand (starting with 0)" << "\n";
    cout << "    -pa    :  Extract possible binding pockets without ligand knowledge\n"
         << "              Please also see '-pa2'!. The number of best pockets can be\n"
         << "              specified in Options2 (default: --n=1 => only the best pocket)." << "\n";
    cout << "              As with -pl the cut radius can be specified in Options2\n"
         << "              (default: 6.0A). You can also try to play around with the\n"
         << "              min_buried factor (see Options2). See ---d mode for\n"
         << "              visualization of the pocket-surface-points." << "\n";
    cout << "    -pa2   :  As -pa an automatic pocket detection. While -pa uses a grid\n"
         << "              based algorithm, -pa2 uses a delaunay triangulation with\n"
         << "              weighted points and determines the possible cavities as the\n"
         << "              difference volume of two alpha shapes. It is slower than -pa,\n"
         << "              but more stable with respect to litte changes in the binding\n"
         << "              pocket and more precise with respect to the determined volume!\n"
         << "              Additionally you can supply a ligand by --s to restrict the\n"
         << "              possible binding pocket to an area 4A around that ligand.\n"
         << "              By --s2 you can specify a multimol2 containing cofactor,\n"
         << "              metals or other atoms that should be considered as protein.\n";
    cout << "              Use --r to specify a radius around the detected border atoms\n"
         << "              (default: 0A). Use this to get not only the bordering atoms,\n"
         << "              but also their neighbours. By --r2 you can specify the lower\n"
         << "              alpha level (default: 2.5A). A lower value will result in finer\n"
         << "              pockets, but you should never go below 2.2A. By --r3 you can\n"
         << "              specify the upper alpha level (default: 8.0A / 15.0A if you\n"
         << "              supply a ligand by --s). A higher value leads to larger pockets.\n";
    cout << "    -pa3   :  As -pa2, but writing pdb-cavities with full residues.\n";
    cout << "    -paf   :  Get alpha cavity profile.\n";
    cout << "    -rmsd  :  Calculates the rmsd values after optimal alignment compared to\n"
         << "              the reference structure specified by --s=ref.pdb" << "\n";
    cout << "              Please note that only CA carbons are used for calculation (as\n"
         << "              for the alignment too)" << "\n";
    cout << "    -sim   :  For an arbitrary number of pdb-files an 8A binding pocket around\n"
         << "              all ligands is determined and the shape similarity compared to a\n"
         << "              reference given by '--s=ref.pdb' is calculated based on an\n"
         << "              optimal superposition of C-alphas. In contrast to '-oa'\n"
         << "              this is not based on a preliminary sequence alignment,\n"
         << "              but only on the spatial arrangement of C-alphas. The output" << "\n";
    cout << "              consists of a rmsd after superposition and the number of\n"
         << "              C-alphas used in superposition. (This function is dated and)\n"
         << "              will be replaced by a faster variant soon." << "\n";
    cout << "    -spat  :  Searches within pdb-files for a secondary structure motif given\n"
         << "              by '--s=ref.pdb'. As in '-sim' mode only the spatial arrangement\n"
         << "              of C-alphas is used for comparison. (This function is dated and)\n"
         << "              will be replaced by a faster variant soon." << "\n";
    cout << "    -ss    :  Search for PDBs with ligands which contain a substructure\n"
         << "              specified by --s=ref.mol2. If you want to have certain atoms\n"
         << "              with same number of bonded heavy atoms as in the reference you" << "\n";
    cout << "              can specify this in your reference file. Therefore you have to\n"
         << "              replace the charge of those atoms by the string 'exact_bonded'.\n"
         << "              e.g. if your reference is an acid you have to specify \n"
         << "              'exact_bonded' for both oxygens, because otherwise you will\n"
         << "              also find esters. The number following the matched structures\n"
         << "              specifies the number of possible matches." << "\n";
    cout << "    -ssh   :  As -ss, but this time requireing the same hybridization states\n"
         << "              for matched atoms. As the substructure search does not care\n"
         << "              about hydrogens, with '-ss' you would also find secondary\n"
         << "              alcohols if searching for ketones. With '-ssh' matching atoms\n"
         << "              are required to have the same hybridization state to avoid this.\n"
         << "              Nevertheless you still want to have atoms with a variable state.\n"
         << "              Please mark those atoms in your reference by replacing their\n"
         << "              charge by the string 'free_hyb'. If you specify 'exact_bonded'\n"
         << "              this also implies 'free_hyb'. If you don't want this please use\n"
         << "              'exact_match' instead of 'exact_bonded'!" << "\n";
    cout << "    -w     :  Extract all water molecules to a mol2 file" << "\n";
    cout << "    -me    :  Extract all metal atoms to a mol2 file" << "\n";
    cout << "    -gfr   :  Give out the number of free rotatable bonds for each molecule" << "\n";
    cout << "    -bval  :  Returns information about the b-values around all ligands. You\n"
         << "              can specify the radius by --r (default = 6A)\n"
         << "              You could also supply a ligand seperately as pdb or mol2 by --s\n";
    cout << "    -bval2 :  As -bval, but this time considering complete residues within\n"
         << "              the given radius\n";
    cout << "    -tc    :  Returns the occurrence of each atom type per molecule. This may\n"
         << "              be useful if you supply a definition file with special types\n"
         << "              (e.g. donor/acceptor/doneptor) and you want to know how many\n"
         << "              atoms of such a specific type your molecules have. The atom\n"
         << "              type 'X' is never regarded here." << "\n";
    cout << "\n";
}

void show_mol2_opt() {
    cout << "-> Valid options for mol2 input files:" << "\n";
    cout << "    -W     :  Remove water" << "\n";
    cout << "    -H     :  Remove hydrogens" << "\n";
    cout << "    -ha    :  Add hydrogens (using standard protonation states)" << "\n";
    cout << "              Note that this will also reset all atom types!" << "\n";
    cout << "    -R     :  Just rewrite the mol2 (this may fix possible errors of the\n"
         << "              original file)" << "\n";
    cout << "    -T     :  Recalculate all atom and bond types (use --m to specify the\n"
         << "              mode or give a definition file by --d). You can also supply a\n"
         << "              reference by --s2 to ensure the same atom types as in the\n"
         << "              reference. For --s2 the atoms must have the same order in\n"
         << "              all files\n";
    cout << "    -ra    :  Rename atoms" << "\n";
    cout << "    -rm    :  Rename molecules using the input file name\n"
         << "              (my.mol2  ==> my_0 ... my_n)" << "\n";
    cout << "    -rd    :  Remove duplicate molecules from a multimol2 (keeping only the\n"
         << "              first occurrence). Stereo chemistry is not considered (use\n"
         << "              -rd2 if you want this).\n";
    cout << "              Currently also resets all atom types!\n";
    cout << "    -rd2   :  As -rd, but also considering stereo chemistry.\n";
    cout << "    -rda   :  Remove duplicate atoms (atoms within one MOLECULE and equal\n"
         << "              coordinates (equal with respect to a threshold of 0.15A dist.)).\n"
         << "              The first occurrence will be kept.\n";
    cout << "    -ro    :  Rearrange atom entries in a way that hydrogen atoms are always\n"
         << "              behind all other atoms." << "\n";
    cout << "    -rs    :  Keep biggest fragment. Keeps only the biggest structure within\n"
         << "              each MOLECULE entry (e.g. removes chlorines etc.)." << "\n";
    cout << "    -m     :  Merge mol2 (or multimol2) files to one multimol2" << "\n";
    cout << "    -m2    :  As '-m', but merging into a single MOLECULE entry\n";
    cout << "    -a     :  Add mol2 file to multimol2 (specified by --t=yourmulti.mol2)" << "\n";
    cout << "    -clash :  Detect steric clashes with a reference given by '--s=ref.mol2'\n"
         << "              (outputs the vdW-overlap). Please note, that this is intended\n"
         << "              for protein-ligand-clashes only. In other cases the clash-volume\n"
         << "              will not be very precise (changed in a future version)." << "\n";
    cout << "   --es=s  :  Extract mol2 file from multimol2, where 's' is the name\n"
         << "              of the molecule" << "\n";
    cout << "   --e=n   :  Extract mol2 file from multimol2, where 'n' is the number\n"
         << "              of the molecule (starting with 0)" << "\n";
    cout << "    -pl    :  As '-pl' for pdb input." << "\n";
    cout << "    -pp    :  As '-pp' for pdb input." << "\n";
    cout << "    -s     :  Split multimol2 to individual files" << "\n";
    cout << "    -sl    :  If there is more than one molecule within one MOLECULE entry\n"
         << "              it will be splitted. Note that also metals will become\n"
         << "              individual MOLECULEs." << "\n";
    cout << "    -ss    :  Search in an arbitrary number of mol2 files for a substructure\n"
         << "              specified by --s=ref.mol2. If you want to have certain atoms\n"
         << "              with same number of bonded heavy atoms as in the reference you\n"
         << "              can specify this in your reference file. Therefore you have to\n"
         << "              replace the charge of those atoms by the string 'exact_bonded'.\n"
         << "              e.g. if your reference is an acid you have to specify" << "\n";
    cout << "              'exact_bonded' for both oxygens, because otherwise you will also\n"
         << "              find esters. If you specify a target (by --t) all matching\n"
         << "              molecules will be written to that file. Otherwise the output\n"
         << "              (stdout) will be the names of matched structures. The number" << "\n";
    cout << "              following the matched structures specifies the number of\n"
         << "              possible matches." << "\n";
    cout << "    -ssh   :  As -ss, but this time requireing the same hybridization states\n"
         << "              for matched atoms. As the substructure search does not care\n"
         << "              about hydrogens, with '-ss' you would also find secondary\n"
         << "              alcohols if searching for ketones. With '-ssh' matching atoms\n"
         << "              are required to have the same hybridization state to avoid this.\n"
         << "              Nevertheless you still want to have atoms with a variable state.\n"
         << "              Please mark those atoms in your reference by replacing their\n"
         << "              charge by the string 'free_hyb'. If you specify 'exact_bonded'\n"
         << "              this also implies 'free_hyb'. If you don't want this please use\n"
         << "              'exact_match' instead of 'exact_bonded'!" << "\n";
    cout << "    -sso   :  As -ss, but this time requireing the same atom types for matched\n"
         << "              atoms, based on the types as they are specified in the input\n"
         << "              files (So NOT after own atom typing!). You may also supply a\n"
         << "              definition file (--d) to set the corresponding types." << "\n";
    cout << "   --bs=n  :  Blocksplit: Split multimol2 into mol2 files\n"
         << "              containing n molecules" << "\n";
    cout << "    -C     :  Builds the unitcell for the given input structure resulting in\n"
         << "              a multimol2 file. This function requires a CRYSIN-entry in your\n"
         << "              mol2 file to get the cell parameters. If you specify a radius\n"
         << "              (see Options2) the unit cell will be duplicated within that\n"
         << "              radius. Alternatively can can specify the number of duplications\n"
         << "              (see Options2). If no target is specified (see --t in Options2)\n"
         << "              the resulting file will be called 'origname_cryst.mol2'" << "\n";
    cout << "    -V     :  Parse COMMENT-entries" << "\n";
    cout << "    -gfr   :  Gives out the number of free rotatable bonds for each molecule" << "\n";
    cout << "              (use ---d to get the atom id's of all rotatable bonds," << "\n";
    cout << "               additionally use --t=file.txt to redirect the output into a\n"
         << "               textfile). You can specify a protein file (PDB or MOL2) by\n"
         << "              '--s=pro_file'. In this case fconv will check for covalent bonds\n"
         << "              to the protein, consider them as rotatable and report them.\n";
    cout << "              You can specify a reference by '--s2=ref.mol2' to ensure that\n"
         << "              all ligands have the same atom typing as in the reference.\n"
         << "              All ligands MUST be the same molecules as the reference, with\n"
         << "              the same order of atoms!!! (If the reference is a multimol2,\n"
         << "              only the first molecule in this file will be taken as reference)\n";
    cout << "    -gfr2  :  As -gfr but here also single bonds within a planar conjugated\n"
         << "              system are counted\n";
    cout << "    -rmsd  :  Gives out the rmsd values compared to the reference specified\n"
         << "              by --s=ref.mol2. You will get 3 values. The rmsd without\n"
         << "              hydrogens, with hydrogens and the rmsd after optimal alignment\n"
         << "              (without hydrogens)." << "\n";
    cout << "              (See ---d mode to get information which atoms are compared to\n"
         << "               each other!)" << "\n";
    cout << "              A result of -2 indicates different molecules with respect to\n"
         << "              their stereo chemistry, while -1 indicates different molecules\n"
         << "              due to any other reasons.\n";
    cout << "              You may also supply a reference2 by --s2 which ensures the same\n"
         << "              atom types for all molecules as in the reference2. For --s2, the\n"
         << "              atoms must have the same order in all files\n";
    cout << "    -rmsd2 :  While '-rmsd' only works for identical molecules this option\n"
         << "              also works on matching substructures (determining the MCIS\n"
         << "              resulting in the best rmsd value). Due to the harder problem\n"
         << "              and the reason rmsd values have to be calculated for possible\n"
         << "              MCISs 3 times, this option works much slower than '-rmsd'.\n"
         << "              It is also possible to get different optimal MCIS for RMSD\n"
         << "              and RMSD_opt. See ---d for this." << "\n";
    cout << "    -rmsd3 :  Calculates per-residue-rmsd values and outputs to shell.\n"
         << "              Total-rmsd, backbone-rmsd, residue-rmsd and also the\n"
         << "              corresponding values after spatial alignment will be given." << "\n";
    cout << "              (Here a sequence alignment is used for atom-atom matching)" << "\n";
    cout << "    -clust :  Hierarchical complete linkage clustering by rmsd values.\n"
         << "              The default rmsd threshold is 2.0A, but you can specify your\n"
         << "              own threshold by --r= (see options2). The output will be a\n"
         << "              multimol2 file with one representative for each cluster." << "\n";
    cout << "              As for all other cluster modes you can get the corresponding\n"
         << "              dendrogram and similarity matrix as svg files if you specifiy a\n"
         << "              name by --s=name.svg (options 2). You can also get a textfile\n"
         << "              with all cluster information by --s2=name.dat\n";
    cout << "    -clust2:  As -clust, but with a multimol2 as output for each cluster." << "\n";
    cout << "    -clust3:  As -clust, but without writing mol2 files (meant for usage\n"
         << "              with --s and/or --s2).\n";
    cout << "    -clusto:  ADDITIONAL flag for clust or clust2. The rmsd values will be\n"
         << "              calculated based on the order of atoms, thus all molecules must\n"
         << "              have the same number of atoms!\n";
    cout << "    -oa    :  Aligns all input structures on a reference structure specified\n"
         << "              by --s=ref.mol2. Note that this works only for equal structures\n"
         << "              (use -oa2 for substructure alignment)" << "\n";
    cout << "              You can also specify a mol2 or PDB with --s2=to_align_file\n"
         << "              in order to align all structures from this file instead of" << "\n";
    cout << "              the sourcefile. (this also works for -oa2)" << "\n";
    cout << "    -oa2   :  As '-oa', but is also able to align non-identical structures\n"
         << "              based on a subgraph match. (watch your memory if you are\n"
         << "              working with huge molecules!!!)" << "\n";
    cout << "    -coa   :  Calculates the optimal alignment on a reference given\n"
         << "              by --s=ref.mol2 and gives out the corresponding rotation matrix\n"
         << "              and translation vectors (stdout)" << "\n";
    cout << "    -pp    :  Extract binding pockets (from a mol2 protein) around a reference\n"
         << "              ligand (given by --s=ligandfile). The default cut radius may be\n"
         << "              changed by --r=float" << "\n";
    cout << "    -pl    :  As -pp but this time cutting whole residues (Therefore you need\n"
         << "              residue numbers in your mol2)." << "\n";
    cout << "    -sdiv  :  Calculate the geometrical diversity between similar molecules." << "\n";
    cout << "              (Slow and outdated => replaced soon)" << "\n";
    cout << "    -nr    :  Returns the number of of smallest rings for all molecules\n"
         << "              (e.g. 3 for an anthracene)." << "\n";
    cout << "              Use debug_mode (---d) for detailed ring information." << "\n";
    cout << "    -nnr   :  Returns 2 numbers: How many heavy atoms are part of a ring\n"
         << "              and how many are not." << "\n";
    cout << "    -mw    :  Returns the molecular weight for all molecules. Please note that\n"
         << "              the function uses the theoretical number of hydrogens based on\n"
         << "              hybridization and bonding. There may be cases where it fails due\n"
         << "              to strange protonation states." << "\n";
    cout << "    -vol   :  Returns the vdW-volume for each molecule\n";
    cout << "    -tc    :  Returns the occurrence of each atom type per molecule. This may\n"
         << "              be useful if you supply a definition file with special types\n"
         << "              (e.g. donor/acceptor/doneptor) and you want to know how many\n"
         << "              atoms of such a specific type your molecules have. The atom type\n"
         << "              'X' is never regarded here." << "\n";
    cout << "\n";
}

void show_dlg_opt() {
    cout << "-> valid options for dlg input files:" << "\n";
    cout << "    -c     :  Convert to multimol2 using the original mol2\n"
         << "              (specified with --s=origin.mol2)" << "\n";
    cout << "    -a     :  Convert to a multimol2 without the knowledge of an original mol2" << "\n";
    cout << "    -f     :  For dlg files from a docking run with flexible residues" << "\n";
    cout << "              fconv will write out the binding pocket for each ligand,\n"
         << "              considering changed coordinates. In this mode you have to\n"
         << "              supply the corresponding pdb file by --s2=my.pdb" << "\n";
    cout << "              You also have to specify -c or -a (e.g. NOT '-fa', BUT '-f -a'\n"
         << "              for the conversion mode of your ligands." << "\n";
    cout << "    -g     :  As -f, but this time writing a complete pdb file for each\n"
         << "              solution and NO ligands (convert them in a second fconv-run\n"
         << "              by -a or -c\n";
    cout << "   --n=a-b :  Specify which solutions you want to write out (all by default)" << "\n";
    cout << "              --n=10-20  => write out solution 10 to 20 from the dlg file" << "\n";
    cout << "              --n=2-2    => write out only the second solution from the dlg" << "\n";
    cout << "              (this time the counting starts with 1)" << "\n";
//    cout << "    -ra    :  Rename atoms" << "\n";
    cout << "    -E     :  Write out the final docked energy to the multimol2" << "\n";
    cout << "    -EF    :  Write out estimated free energy of binding to the multimol2" << "\n";
    cout << "\n";
}

void show_sdf_opt() {
    cout << "-> valid options for sdf input files:" << "\n";
    cout << "    -c     :  Convert to multimol2" << "\n";
    cout << "   --es=s  :  Extract mol2 file from sdf, where 's' is the name\n"
         << "              of the molecule" << "\n";
    cout << "   --e=n   :  Extract mol2 file from sdf, where 'n' is the number\n"
         << "              of the molecule (starting with 0)" << "\n";
    cout << "    -s     :  Split into individual mol2 files" << "\n";
    cout << "    -ss    :  As -ss for mol2. The substructure must be supplied as mol2.\n";
    cout << "    -ssh   :  As -ssh for mol2. The substructure must be supplied as mol2.\n";
    cout << "    -sso   :  As -ss, but this time requireing the same atom types for matched\n"
         << "              atoms. In contrast to -sso for mol2, here the atom types will\n"
         << "              be resetted.\n";
    cout << "   --bs=n  :  Blocksplit: Split sdf into mol2 files containing n molecules\n";
    cout << "    -gfr   :  Gives out the number of free rotatable bonds for each molecule\n";
    cout << "              (use ---d to get the atom id's of all rotatable bonds" << "\n";
    cout << "               additionally use --t=file.txt to redirect the output into a\n"
         << "               textfile). You can specify a protein file (PDB or MOL2) by\n"
         << "              '--s=pro_file'. In this case fconv will check for covalent bonds\n"
         << "              to the protein, consider them as rotatable and report them.\n";
    cout << "              You can specify a reference by '--s2=ref.mol2' to ensure that\n"
         << "              all ligands have the same atom typing as in the reference.\n"
         << "              All ligands MUST be the same Molecules as the reference, with\n"
         << "              the same order of atoms!!! (If the reference is a multimol2,\n"
         << "              only the first molecule in this file will be taken as reference)\n";
    cout << "    -gfr2  :  As -gfr but here also single bonds within a planar conjugated\n"
         << "              system are counted\n";
    cout << "    -rmsd  :  Gives out the rmsd values compared to the reference specified\n"
         << "              by --s=ref.mol2. You will get 3 values. The rmsd without\n"
         << "              hydrogens, with hydrogens and the rmsd after optimal alignment\n"
         << "              (without hydrogens)." << "\n";
    cout << "              (See ---d mode to get information which atoms are compared to\n"
         << "               each other!)" << "\n";
    cout << "              A result of -2 indicates different molecules with respect to\n"
         << "              their stereo chemistry, while -1 indicates different molecules\n"
         << "              due to any other reasons.\n";
    cout << "              You may also supply a reference2 by --s2 which ensures the same\n"
         << "              atom types for all molecules as in the reference2. For --s2, the\n"
         << "              atoms must have the same order in all files\n";
    cout << "    -rmsd2 :  While '-rmsd' only works for identical molecules this option\n"
         << "              also works on matching substructures (determining the MCIS\n"
         << "              resulting in the best rmsd value). Due to the harder problem\n"
         << "              and the reason rmsd values have to be calculated for possible\n"
         << "              MCISs 3 times, this option works much slower than '-rmsd'.\n"
         << "              It is also possible to get different optimal MCIS for RMSD\n"
         << "              and RMSD_opt. See ---d for this." << "\n";
    cout << "    -clust :  Hierarchical complete linkage clustering by rmsd values.\n"
         << "              The default rmsd threshold is 2.0A, but you can specify your\n"
         << "              own threshold by --r= (see options2). The output will be a\n"
         << "              multimol2 file with one representative for each cluster." << "\n";
    cout << "              As for all other cluster modes you can get the corresponding\n"
         << "              dendrogram and similarity matrix as svg files if you specifiy a\n"
         << "              name by --s=name.svg (options 2). You can also get a textfile\n"
         << "              with all cluster information by --s2=name.dat\n";
    cout << "    -clust2:  As -clust, but with a multimol2 as output for each cluster." << "\n";
    cout << "    -clust3:  As -clust, but without writing mol2 files (meant for usage\n"
         << "              with --s and/or --s2).\n";
    cout << "    -clusto:  ADDITIONAL Flag for clust or clust2. The rmsd values will be\n"
         << "              calculated based on the order of atoms, thus all molecules must\n"
         << "              have the same number of atoms!\n";
    cout << "    -oa    :  Aligns all input structures on a reference structure specified\n"
         << "              by --s=ref.mol2. Note that this works only for equal structures\n"
         << "              (use -oa2 for substructure alignment)" << "\n";
    cout << "    -oa2   :  As '-oa', but is also able to align non-identical structures\n"
         << "              based on a subgraph match. (watch your memory if you are\n"
         << "              working with huge molecules!!!)" << "\n";
    cout << "    -coa   :  Calculates the optimal alignment on a reference given\n"
         << "              by --s=ref.mol2 and gives out the corresponding rotation matrix\n"
         << "              and translation vectors (stdout)" << "\n";
    cout << "    -nr    :  Returns the number of of smallest rings for all molecules\n"
         << "              (e.g. 3 for an anthracene)." << "\n";
    cout << "              Use debug_mode (---d) for detailed ring information." << "\n";
    cout << "    -nnr   :  Returns 2 numbers: How many heavy atoms are part of a ring\n"
         << "              and how many are not." << "\n";
    cout << "    -mw    :  Returns the molecular weight for all molecules. Please note that\n"
         << "              the function uses the theoretical number of hydrogens based on\n"
         << "              hybridization and bonding. There may be cases where it fails due\n"
         << "              to strange protonation states." << "\n";
    cout << "    -vol   :  Returns the vdW-volume for each molecule\n";
    cout << "    -tc    :  Returns the occurrence of each atom type per molecule. This may\n"
         << "              be useful if you supply a definition file with special types\n"
         << "              (e.g. donor/acceptor/doneptor) and you want to know how many\n"
         << "              atoms of such a specific type your molecules have. The atom type\n"
         << "              'X' is never regarded here." << "\n";
    cout << "\n";
}

void show_cif_opt() {
    cout << "-> valid options for cif input files:" << "\n";
    cout << "    -c     :  Convert to multimol2" << "\n";
    cout << "    -C     :  Build the crystal package (see '-C' for mol2 files). In contrast\n"
         << "              to the same option for mol2 files, here it is not possible to\n"
         << "              use the debug mode in the moment" << "\n";
    cout << "\n";
}

void show_opt2() {
    cout << "------------" << "\n";    
    cout << "  Options2 :" << "\n";
    cout << "------------" << "\n";
    cout << "    -h     :  Do not use ligand hydrogens in pocket extraction.\n";
    cout << "    -s     :  Output single mol2 files instead of a multimol2 using the\n"
         << "              molecules names as filenames (dlg only)" << "\n";
    cout << "   --t=s   :  Use 's' as name for the target file." << "\n";
    cout << "   --d=s   :  Specify a definition file for the atom typing." << "\n";
    cout << "              To get an example file see --M in options1." << "\n"; 
    cout << "   --s=s   :  Specify a source file if necessary" << "\n";
    cout << "              Specify a name for svg files in cluster modes\n";
    cout << "   --s2=s  :  Specify an additional source if necessary" << "\n";
    cout << "              Specify a name for a cluster data file\n";
    cout << "   --f=t   :  Target file type (t = 'p' for pdb  /  t = 'm' for mol2)\n"
         << "              (default is input file type)" << "\n";
    cout << "   --r=f   :  Specify the cut radius for the binding pocket (default: 6.0)" << "\n";
    cout << "              Also used as radius for unit cell duplication (see -C mode for\n"
         << "              mol2 files)" << "\n";
    cout << "              Also used as threshold for rmsd clustering (see -clust)" << "\n";
    cout << "              Also used as 'surrounding radius' in -pa2 mode.\n";
    cout << "              Also specifies the maximum displacement for -oa2 (default: 1.4A)\n";
    cout << "   --r2=f  :  Specify the lower alpha value for -pa2 mode.\n";
    cout << "   --r3=f  :  Specify the upper alpha value for -pa2 mode.\n";
    cout << "   --n=n   :  Specify the number of pockets for '-pa' (default: 1)" << "\n";
    cout << "              Also used to specify how often a unit cell should be duplicated\n"
         << "              in all directions in the -C mode for mol2 files (default is 0)" << "\n";
    cout << "              Also gives the threshold for free rotatable bonds in -gfr[2]\n"
         << "              mode (only molecules with FR <= n will be written).\n";
    cout << "   --m=n   :  Specify the maximum number of ring-members for atom typing\n"
         << "              (default = 10). This number only refers to the smallest rings,\n"
         << "              not fused systems (e.g. in anthracene the max number of\n"
         << "              ringmembers is 6)" << "\n";
    cout << "   --v=i   :  Specify the minimum buriedness for calculation of binding\n"
         << "              pockets (default: 18)" << "\n";
    cout << "              The buriedness is an integer between 0 and 26, but only values\n"
         << "              between 15 and 22 are usefull here. A higher value would\n"
         << "              result in smaller pockets." << "\n";
    cout << "  ---d     :  This turns on the unofficial debug mode. For users this may be\n"
         << "              interesting in '-pa' mode as it writes out pymol visualization\n"
         << "              files for the cavity surface points (start them from pymol with\n"
         << "              the command 'run'). In '-C' mode it writes out a pymol\n"
         << "              visualisation file for each symmetry element. In '-nr' mode\n"
         << "              it gives detailed ring information." << "\n";
    cout << "              It also writes out which atoms are compared to each other for\n"
         << "              rmsd caculation." << "\n";
    cout << endl;
}

void show_pdb_opt2() {
    cout << "------------" << "\n";    
    cout << "  Options2 :" << "\n";
    cout << "------------" << "\n";
    cout << "    -h     :  Do not use ligand hydrogens in pocket extraction.\n";
    cout << "   --t=s   :  Use 's' as name for the target file." << "\n";
    cout << "   --d=s   :  Specify a definition file for the atom typing." << "\n";
    cout << "              To get an example file see --M in options1." << "\n"; 
    cout << "   --s=s   :  Specify a source file if necessary" << "\n";
    cout << "   --s2=s  :  Specify an additional source if necessary" << "\n";
    cout << "   --f=t   :  Target file type (t = 'p' for pdb  /  t = 'm' for mol2)\n"
         << "              (default is input file type)" << "\n";
    cout << "   --r=f   :  Specify the cut radius for the binding pocket (default: 6.0)" << "\n";
    cout << "              Also used as 'surrounding radius' in -pa2 mode.\n";
    cout << "              Also specifies the maximum displacement for -oa2 (default: 1.4A)\n";
    cout << "   --r2=f  :  Specify the lower alpha value for -pa2 mode.\n";
    cout << "   --r3=f  :  Specify the upper alpha value for -pa2 mode.\n";
    cout << "   --n=n   :  Specify the number of pockets for '-pa' (default: 1)" << "\n";
    cout << "   --m=n   :  Specify the maximum number of ring-members for atom typing\n"
         << "              (default = 10). This number only refers to the smallest rings,\n"
         << "              not fused systems (e.g. in anthracene the max number of\n"
         << "              ringmembers is 6)" << "\n";
    cout << "   --v=i   :  Specify the minimum buriedness for calculation of binding\n"
         << "              pockets (default: 18)" << "\n";
    cout << "              The buriedness is an integer between 0 and 26, but only values\n"
         << "              between 15 and 22 are usefull here. A higher value would\n"
         << "              result in smaller pockets." << "\n";
    cout << "  ---d     :  This turns on the unofficial debug mode. For users this may be\n"
         << "              interesting in '-pa' mode as it writes out python visualization\n"
         << "              files for the cavity surface points." << "\n";
    cout << endl;
}

void show_mol2_opt2() {
    cout << "------------" << "\n";
    cout << "  Options2 :" << "\n";
    cout << "------------" << "\n";
    cout << "    -h     :  Do not use ligand hydrogens in pocket extraction.\n";
    cout << "   --t=s   :  Use 's' as name for the target file." << "\n";
    cout << "   --d=s   :  Specify a definition file for the atom typing." << "\n";
    cout << "              To get an example file see --M in options1." << "\n"; 
    cout << "   --s=s   :  Specify a source file if necessary" << "\n";
    cout << "              Specify a name for svg files in cluster modes\n";
    cout << "   --s2=s     Specify a name for a cluster data file\n";
    cout << "   --f=t   :  Target file type (t = 'p' for pdb  /  t = 'm' for mol2)\n"
         << "              (default is input file type)" << "\n";
    cout << "   --r=f   :  Specify the radius for unit cell duplication (see -C mode)\n"
         << "              or the radius for -pp mode" << "\n";
    cout << "   --n=n   :  Specify how often a unit cell should be duplicated in\n"
         << "              all directions in the -C mode (default is 0). Also gives the\n"
         << "              threshold for free rotatable bonds in -gfr[2] mode\n"
         << "              (only molecules with FR <= n will be written).\n";
    cout << "   --m=n   :  Specify the maximum number of ring-members for atom typing\n"
         << "              (default = 10). This number only refers to the smallest rings,\n"
         << "              not fused systems (e.g. in anthracene the max number of\n"
         << "              ringmembers is 6)" << "\n";
    cout << "  ---d     :  This turns on the unofficial debug mode." << "\n";
    cout << "              It writes out which atoms are compared to each other for\n"
         << "              rmsd caculation." << "\n";
    cout << endl;
}

void show_dlg_opt2() {
    cout << "------------" << "\n";
    cout << "  Options2 :" << "\n";
    cout << "------------" << "\n";
    cout << "    -s     :  Output single mol2 files instead of a multimol2 using the\n"
         << "              molecules names as filenames" << "\n";
    cout << "   --t=s   :  Use 's' as name for the target file." << "\n";
    cout << "   --d=s   :  Specify a definition file for the atom typing." << "\n";
    cout << "              To get an example file see --M in options1." << "\n";
    cout << "   --m=n   :  Specify the maximum number of ring-members for atom typing\n"
         << "              (default = 10). This number only refers to the smallest rings,\n"
         << "              not fused systems (e.g. in anthracene the max number of\n"
         << "              ringmembers is 6)" << "\n";
    cout << "   --s=s   :  Specify a source file if necessary" << "\n";
    cout << "   --s2=s  :  Specify an additional source if necessary" << "\n";
    cout << endl;
}

void show_cif_opt2() {
    cout << "------------" << "\n";
    cout << "  Options2 :" << "\n";
    cout << "------------" << "\n";
    cout << "   --t=s   :  Use 's' as name for the target file." << "\n";
    cout << "   --d=s   :  Specify a definition file for the atom typing." << "\n";
    cout << "              To get an example file see --M in options1." << "\n"; 
    cout << "   --r=f   :  Specify the radius for unit cell duplication (see -C mode)" << "\n";
    cout << "   --n=n   :  Specify how often a unit cell should be duplicated\n"
         << "              in all directions" << "\n";
    cout << "              in the -C mode (default is 0)" << "\n";
    cout << "   --m=n   :  Specify the maximum number of ring-members for atom typing\n"
         << "              (default = 10). This number only refers to the smallest rings,\n"
         << "              not fused systems (e.g. in anthracene the max number of\n"
         << "              ringmembers is 6)" << "\n";
    cout << endl;
}

void show_notes() {
    cout << "----------------" << "\n";
    cout << "Important notes:" << "\n";
    cout << "----------------" << "\n";
    cout << "  If your target is of the same file type as your source and you don't specify\n"
         << "  a name for it, your original file will be !!! OVERWRITTEN !!! in modes\n"
         << "  '-c', '-f', '-oa', '-rd', '-R', '-W', '-H' and '-T' !!!" << "\n";
    cout << "  In other cases a new filename is generated automatically if you don't\n"
         << "  specify it. Also if you specify a target file name it may be modified by\n"
         << "  fconv (e.g. if you split a multimol2 with the option '--t=splitted.mol2'\n"
         << "  your files will be called 'splitted_0.mol2' to 'splitted_n.mol2')." << "\n" << "\n";
}

void show_big_help() {
    show_header();
    
    show_gen_opt();
    
    show_pdb_opt();
    
    show_mol2_opt();
    
    show_dlg_opt();
    
    show_sdf_opt();
    
    show_cif_opt();
    
    show_opt2();
    
    show_notes();
    
    cout << endl;
}

int main(int argc,char *argv[]) {
    vector<string> s;
    
    int const max_mol2_block = 5000; //! how many molecules may be processed in one block
    
    vector<string> sourcefiles;
    string sourcefile = "X";
    string pdb_source = "X";
    string targetfile = "X";
    string targettype = "X";
    string def_file ="X";
    string tfile;
    string mode;
    string dummy;
    bool gethyd = true;
    bool getwat = true;
    bool convert = false;
    bool full_convert = false;
    bool rewrite = false;
    bool fix_pdb = false;
    int atom_typing = 1;
    bool extract_ligs_seperate = false;
    bool retype = false;
    bool parse_ligands = true;
    bool parse_comments = false;
    bool extract_all_ligs_cav = false;
    int verb = 1;
    float cut_radius = 6.;
    float max_displacement = 1.4;
    float unit_cell_exp_rad = 0.;
    unsigned int unit_cell_multi = 0;
    bool merge_mol2 = false;
    bool merge_mol2_in1mol = false;
    bool add_mol2 = false;
    bool extract_mol2 = false;
    bool extract_mol2_by_name = false;
    bool extract_water = false;
    bool extract_metals = false;
    unsigned int extract_number = 0;
    string extract_name = "X";
    bool split_multimol = false;
    bool blocksplit = false;
    int blocksplit_size = 0;
    bool cut_by_lig = false;
    unsigned int cut_lig_number = 0;
    bool rename_atoms = false;
    bool rename_molecules = false;
    bool remove_duplicates = false;
    bool remove_duplicate_atoms = false;
    bool rename_files = false;
    bool dlg2mol2 = false;
    bool dlg2singlemol2 = false;
    unsigned int high = 99999999;
    unsigned int low = 1;
    bool write_energy = false;
    bool write_free_energy = false;
    bool auto_dlg = false;
    bool auto_cav = false;
    bool auto_cav2 = false;
    bool auto_cav3 = false;
    bool cut_full_res = true;
    bool no_hyd_for_poc = false;
    int min_buried = 18;
    unsigned int max_poc = 1;
    bool debug_mode = false;
    bool build_crystal = false;
    bool get_rot_bonds = false;
    bool rmsd = false;
    bool rmsd2 = false;
    bool rmsd3 = false;
    bool optalign = false;
    bool optalign2 = false;
    bool pock_align = false;
    bool sdiv = false;
    bool flex_res_dlg = false;
    bool full_pdb_from_dlg = false;
    bool get_n_rings = false;
    bool get_n_rings2 = false;
    bool get_mol_weight = false;
    bool type_count = false;
    bool add_hyd = false;
    bool calc_optalign_matrix = false;
    bool search_substruct = false;
    bool ss_with_exact_hybrid = false;
    bool ss_with_exact_type = false;
    bool shape_sim = false;
    bool s_pattern = false;
    bool clashes = false;
    bool intern_split = false;
    bool reorder_atoms = false;
    bool remove_peptidic = false;
    bool keep_biggest_fragment = false;
    bool cluster = false;
    bool write_full_clusters = false;
    bool cluster_without_mol2 = false;
    bool only_C_alpha = true;
    bool alpha_profile = false;
    float cluster_threshold = 2.0;
    bool cluster_rmsds_ordered = false;
    float alpha1 = 2.5;
    float alpha2 = 8.0;
    float surr_radius = 0.;
    bool get_vol = false;
    
    int max_ring_members = 10;
    
    bool no_planar_free_rot = true;
    
    bool prot_acids = false;
    bool prot_guanidin = true;
    bool prot_amidin = true;
    bool prot_amin = false;
    bool prot_phosphate = false;
    bool prot_sulfate = false;

    bool kekulize_aromatics = false;
    bool kekulize_charged = false;
    bool allow_charged_aromatics = true;

    bool bvals = false;

    bool pdbqt_mode = false;

    bool check_stereo = false;

    s.resize(argc);
    for (int i=1; i<argc; ++i) {
        s[i].assign(argv[i]);
    }
    
    bool opt2 = false;

    if (argc < 2) {show_help(); exit(0);}
    
    for (int i=1; i<argc; ++i) {
        if (s[i][0] == '-') {
            if (s[i].size() < 2) continue;
            if (s[i][1] == 'H') gethyd = false;
            else if (s[i][1] == 'W') getwat = false;
            else if (s[i][1] == 'C') build_crystal = true;
            else if (s[i][1] == 'f') {
                if (s[i].size() > 2) {
                    if (s[i] == "-fr") only_C_alpha = false;
                } else full_convert = true; flex_res_dlg = true;
            }
            else if (s[i][1] == 'R') rewrite = true;
            else if (s[i][1] == 'F') fix_pdb = true;
            else if (s[i][1] == 'T') retype = true;
            else if (s[i][1] == 'L') {
                parse_ligands = false;
                if (s[i].size() > 2) {
                    if (s[i] == "-LL") remove_peptidic = true;
                }
            } else if (s[i][1] == 'a') {add_mol2 = true; auto_dlg = true;}
            else if (s[i][1] == 'l') extract_ligs_seperate = true;
            else if (s[i][1] == 's' && !opt2) {
                if (s[i].size() > 2) {
                    if (s[i] == "-ss") {
                        search_substruct = true;
                    } else if (s[i] == "-sl") {
                        intern_split = true;
                    } else if (s[i] == "-ssh") {
                        ss_with_exact_hybrid = true;
                        search_substruct = true;
                    } else if (s[i] == "-sso") {
                        ss_with_exact_type = true;
                    } else if (s[i] == "-sdiv") {
                        sdiv = true;
                    } else if (s[i] == "-sim") {
                        shape_sim = true;
                    } else if (s[i] == "-spat") {
                        s_pattern = true;
                    } else {
                        cerr << c_message<cERROR>(s[i]) << " is no valid option (see 'fconv -h')" << endl;
                        exit(1);
                    }
                } else split_multimol = true;
            }
            else if (s[i][1] == 'c') {
                if (s[i].size() > 2) {
                    if (s[i] == "-coa") {
                        calc_optalign_matrix = true;
                    } else if (s[i] == "-clash") {
                        clashes = true;
                    } else if (s[i] == "-clust") {
                        cluster = true;
                    } else if (s[i] == "-clust2") {
                        cluster = true;
                        write_full_clusters = true;
                    } else if (s[i] == "-clust3") {
                        cluster = true;
                        cluster_without_mol2 = true;
                    } else if (s[i] == "-clusto") {
                        gethyd = false;
                        cluster = true;
                        cluster_rmsds_ordered = true;
                    } else {
                        cerr << c_message<cERROR>(s[i]) << " is no valid option (see 'fconv -h')" << endl;
                        exit(1);
                    }
                } else {
                    convert = true;
                    dlg2mol2 = true;
                }
            }
            else if (s[i][1] == 's' && opt2) dlg2singlemol2 = true;
            else if (s[i][1] == 'h' && opt2) no_hyd_for_poc = true;
            else if (s[i][1] == 'V') parse_comments = true;
            else if (s[i] == "-bval") {bvals = true; cut_full_res = false;}
            else if (s[i] == "-bval2") {bvals = true; cut_full_res = true;}
            else if (s[i][1] == 'e') {show_example(); exit(0);}
            else if (s[i][1] == 'w') {
                if (s[i].size() > 3) {
                    if (s[i][2] == 'i') {
                        start_wizard();
                        exit(0);
                    } else {
                        cerr << c_message<cERROR>(s[i]) << " is no valid option (see 'fconv -h')" << endl;
                        exit(1);
                    }
                } else extract_water = true;
            }
            else if (s[i][1] == 'o') {
                if (s[i].size() > 2) {
                    if (s[i][2] == 'a') optalign = true;
                    else {
                        cerr << c_message<cERROR>(s[i]) << " is no valid option (see 'fconv -h')" << endl;
                        exit(1);
                    }
                    if (s[i].size() > 3) {
                        if (s[i][3] == '2') {optalign = false; optalign2 = true;}
                        else {
                            cerr << c_message<cERROR>(s[i]) << " is no valid option (see 'fconv -h')" << endl;
                            exit(1);
                        }
                    }
                } else {
                    cerr << c_message<cERROR>(s[i]) << " is no valid option (see 'fconv -h')" << endl;
                    exit(1);
                }
            }
            else if (s[i][1] == 'p') {
                if (s[i].size() > 2) {
                    if (s[i][2] == 'l') extract_all_ligs_cav = true;
                    else if (s[i][2] == 'p') {extract_all_ligs_cav = true; cut_full_res = false;}
                    else if (s[i][2] == 'a') {
                        if (s[i].size() > 3) {
                            if (s[i][3] == '2') auto_cav2 = true;
                            else if (s[i][3] == '3') auto_cav3 = true;
                            else if (s[i][3] == 'f') alpha_profile = true;
                            else {
                                cerr << c_message<cERROR>(s[i]) << " is no valid option (see 'fconv -h')" << endl;
                                exit(1);
                            }
                        } else auto_cav = true;
                    } else if (s[i][2] == 'o') pock_align = true;
                    else {
                        cerr << c_message<cERROR>(s[i]) << " is no valid option (see 'fconv -h')" << endl;
                        exit(1);
                    }
                } else {
                    cerr << c_message<cERROR>(s[i]) << " is no valid option (see 'fconv -h')" << endl;
                    exit(1);
                }
            }
            else if (s[i][1] == 'E') {
                if (s[i].size() > 2) {
                    if (s[i][2] == 'F') write_free_energy = true;
                    else {
                        cerr << c_message<cERROR>(s[i]) << " is no valid option (see 'fconv -h')" << endl;
                        exit(1);
                    }
                } else {
                    write_energy = true;
                }
            }
            else if (s[i][1] == 'h') {
                if (s[i].size() > 2) {
                    if (s[i][2] == 'p') {
                        show_header(); show_gen_opt(); show_pdb_opt(); show_pdb_opt2();
                        exit(0);
                    }
                    if (s[i][2] == 'm') {
                        show_header(); show_gen_opt(); show_mol2_opt(); show_mol2_opt2();
                        exit(0);
                    }
                    if (s[i][2] == 'd') {
                        show_header(); show_gen_opt(); show_dlg_opt(); show_dlg_opt2();
                        exit(0);
                    }
                    if (s[i][2] == 's') {
                        show_header(); show_gen_opt(); show_sdf_opt(); show_mol2_opt2();
                        exit(0);
                    }
                    if (s[i][2] == 'c') {
                        show_header(); show_gen_opt(); show_cif_opt(); show_cif_opt2();
                        exit(0);
                    }
                    if (s[i][2] == 'a') {
                        add_hyd = true;
                    }
                } else {
                    show_big_help(); exit(0);
                }
            }
            else if (s[i][1] == 'm') {
                if (s[i].size() == 3) {
                    if (s[i][2] == 'w') {
                        get_mol_weight = true;
                    } else if (s[i][2] == 'e') {
                        extract_metals = true;
                    } else if (s[i][2] == '2') {
                        merge_mol2_in1mol = true;
                    } else {
                        cerr << c_message<cERROR>(s[i]) << " is no valid option (see 'fconv -h')" << endl;
                        exit(1);
                    }
                } else if (s[i].size() == 2 ) {
                    merge_mol2 = true;
                } else {
                    cerr << c_message<cERROR>(s[i]) << " is no valid option (see 'fconv -h')" << endl;
                    exit(1);
                }
            }
            else if (s[i][1] == 'n') {
                if (s[i].size() == 3) {
                    if (s[i][2] == 'r') get_n_rings = true;
                } else if (s[i] == "-nnr") {
                    get_n_rings2 = true;
                } else {
                    cerr << c_message<cERROR>(s[i]) << " is no valid option (see 'fconv -h')" << endl;
                    exit(1);
                }
            }
            else if (s[i][1] == 't') {
                if (s[i].size() == 3) {
                    if (s[i][2] == 'c') type_count = true;
                } else {
                    cerr << c_message<cERROR>(s[i]) << " is no valid option (see 'fconv -h')" << endl;
                    exit(1);
                }
            }
            else if (s[i][1] == 'r') {
                if (s[i].size() > 2) {
                    if (s[i][2] == 'a') rename_atoms = true;
                    else if (s[i][2] == 'o') reorder_atoms = true;
                    else if (s[i][2] == 's') keep_biggest_fragment = true;
                    else if (s[i][2] == 'm') {
                        if (s[i].size() > 4) {
                            if (s[i].size() > 5) {
                                if (s[i][3] == 's' && s[i][4] == 'd' && s[i][5] == '2') rmsd2 = true;
                                if (s[i][3] == 's' && s[i][4] == 'd' && s[i][5] == '3') rmsd3 = true;
                            } else if (s[i][3] == 's' && s[i][4] == 'd') rmsd = true;
                            else rename_molecules = true;
                        } else rename_molecules = true;
                    }
                    else if (s[i][2] == 'f') rename_files = true;
                    else if (s[i][2] == 'd') {
                        if (s[i].size() > 3) {
                            if (s[i][3] == 'a') remove_duplicate_atoms = true;
                            else if (s[i][3] == '2') {
                                remove_duplicates = true;
                                check_stereo = true;
                            } else {
                                cerr << c_message<cERROR>(s[i]) << " is no valid option (see 'fconv -h')" << endl;
                                exit(1);
                            }
                        } else remove_duplicates = true;
                    }
                } else {
                    cerr << c_message<cERROR>(s[i]) << " is no valid option (see 'fconv -h')" << endl;
                    exit(1);
                }
            }
            else if (s[i][1] == 'g') {
                if (s[i].size() > 3) {
                    if (s[i][2] == 'f' && s[i][3] == 'r') get_rot_bonds = true;
                    if (s[i].size() > 4) {
                        if (s[i][4] == '2') no_planar_free_rot = false;
                    }
                } else {
                    full_pdb_from_dlg = true;
                }
            }
            else if (s[i][1] == 'v') {
                if (s[i] == "-vol") get_vol = true;
            }
            else if (s[i] == "-KA") {
                kekulize_aromatics = true;
            }
            else if (s[i] == "-KC") {
                kekulize_charged = true;
            }
            else if (s[i] == "-NCA") {
                allow_charged_aromatics = false;
            }
            else if (s[i][1] == '-') {
                if (s[i].size() < 3) {
                    cerr << c_message<cERROR>(s[i]) << " is no valid option (see 'fconv -h')" << endl;
                    exit(1);
                }
                if (s[i][2] == 'm' && !opt2) {
                    if (s[i].size() == 5) {
                        if (s[i][3] == '=') {
                            if (s[i][4] == '0') atom_typing = 0;
                            if (s[i][4] == '1') atom_typing = 1;
                            if (s[i][4] == '2') atom_typing = 2;
                        }
                    } else {
                        cerr << c_message<cERROR>(s[i]) << " is no valid option (see 'fconv -h')" << endl;
                        exit(1);
                    }
                }
                else if (s[i][2] == 'M' && !opt2) {
                    if (s[i].size() == 5) {
                        if (s[i][3] == '=') {
                            if (s[i][4] == '0') {
                                write_def_file(0,"Atom_types_mode0.def");
                                cerr << "Definition file written to Atom_types_mode0.def" << endl;
                                exit(0);
                            } else if (s[i][4] == '1') {
                                write_def_file(1,"Atom_types_mode1.def");
                                cerr << "Definition file written to Atom_types_mode1.def" << endl;
                                exit(0);
                            } else if (s[i][4] == '2') {
                                write_def_file(2,"Atom_types_mode2.def");
                                cerr << "Definition file written to Atom_types_mode2.def" << endl;
                                exit(0);
                            }
                        }
                    } else {
                        cerr << c_message<cERROR>(s[i]) << " is no valid option (see 'fconv -h')" << endl;
                        exit(1);
                    }
                }
                else if (s[i][2] == 'p' && !opt2) {
                    if (s[i].size() == 10) {
                        if (s[i][4] == '0') prot_acids = true;
                        if (s[i][5] == '0') prot_sulfate = true;
                        if (s[i][6] == '0') prot_phosphate = true;
                        if (s[i][7] == '0') prot_guanidin = false;
                        if (s[i][8] == '0') prot_amidin = false;
                        if (s[i][9] == '1') prot_amin = true;
                    } else {
                        cerr << c_message<cERROR>(s[i]) << " is no valid option (see 'fconv -h')" << endl;
                        exit(1);
                    }
                }
                else if (s[i][2] == 'v' && !opt2) {
                    if (s[i].size() == 5) {
                        if (s[i][3] == '=') {
                            if (s[i][4] == '0') verb = 0;
                            if (s[i][4] == '1') verb = 1;
                            if (s[i][4] == '2') verb = 2;
                        }
                    } else {
                        cerr << c_message<cERROR>(s[i]) << " is no valid option (see 'fconv -h')" << endl;
                        exit(1);
                    }
                }
                else if (s[i][2] == 'p') {
                    if (s[i].size() > 5) {
                        dummy.assign(s[i],5,s[i].size()-5);
                        istringstream is;
                        is.str(dummy);
                        if (s[i][3] == 'l') {
                            try {
                                is >> cut_lig_number;
                            } catch(...) {
                                cerr << c_message<cERROR>(dummy) << " is no valid number" << endl;
                                exit(1);
                            }
                        }
                        cut_by_lig = true;
                    } else {
                        cerr << c_message<cERROR>(s[i]) << " is no valid option (see 'fconv -h')" << endl;
                        exit(1);
                    }
                }
                else if (s[i][2] == 'e') {
                    if (s[i].size() > 4) {
                        if (s[i][3] == '=') {
                            dummy.assign(s[i],4,s[i].size()-4);
                            istringstream is;
                            is.str(dummy);
                            try {
                                is >> extract_number;
                            } catch(...) {
                                cerr << c_message<cERROR>(dummy) << " is no valid number" << endl;
                                exit(1);
                            }
                            extract_mol2 = true;
                        } else if (s[i][3] == 's') {
                            dummy.assign(s[i],5,s[i].size()-5);
                            istringstream is;
                            is.str(dummy);
                            try {
                                is >> extract_name;
                            } catch(...) {
                                cerr << c_message<cERROR>(dummy) << " is no valid name" << endl;
                                exit(1);
                            }
                            extract_mol2_by_name = true;
                        }
                    } else {
                        cerr << c_message<cERROR>(s[i]) << " is no valid option (see 'fconv -h')" << endl;
                        exit(1);
                    }
                }
                else if (s[i][2] == 'b') {
                    if (s[i].size() > 4) {
                        if (s[i][3] == 's') {
                            dummy.assign(s[i],5,s[i].size()-5);
                            istringstream is;
                            is.str(dummy);
                            try {
                                is >> blocksplit_size;
                            } catch(...) {
                                cerr << c_message<cERROR>(dummy) << " is no valid number" << endl;
                                exit(1);
                            }
                            blocksplit = true;
                        }
                    } else {
                        cerr << c_message<cERROR>(s[i]) << " is no valid option (see 'fconv -h')" << endl;
                        exit(1);
                    }
                }
                else if (s[i][2] == 'n' && !opt2) {
                    if (s[i].size() > 4) {
                        bool breaker = false;
                        int counter = 4;
                        while (!breaker) {
                            counter++;
                            if (s[i][counter] == '-') breaker = true;
                        }
                        string low_num(s[i],4,counter-2);
                        istringstream is;
                        is.clear(); is.str(low_num); is >> low;
                        string high_num(s[i],counter+1,s[i].size()-counter-1);
                        is.clear(); is.str(high_num); is >> high;
                    } else {
                        cerr << c_message<cERROR>(s[i]) << " is no valid option (see 'fconv -h')" << endl;
                        exit(1);
                    }
                }
                
                else if (s[i][2] == 't' && opt2) {
                    if (s[i].size() > 4) {
                        targetfile.assign(s[i],4,s[i].size()-4);
                    } else {
                        cerr << c_message<cERROR>(s[i]) << " is no valid option (see 'fconv -h')" << endl;
                        exit(1);
                    }
                }
                else if (s[i][2] == 'd' && opt2) {
                    if (s[i].size() > 4) {
                        def_file.assign(s[i],4,s[i].size()-4);
                    } else {
                        cerr << c_message<cERROR>(s[i]) << " is no valid option (see 'fconv -h')" << endl;
                        exit(1);
                    }
                }
                else if (s[i][2] == 's' && opt2) {
                    if (s[i].size() > 4) {
                        if (s[i][3] == '=') sourcefile.assign(s[i],4,s[i].size()-4);
                        else if (s[i][3] == '2') pdb_source.assign(s[i],5,s[i].size()-5);
                        else {
                            cerr << c_message<cERROR>(s[i]) << " is no valid option (see 'fconv -h')" << endl;
                            exit(1);
                        }
                    } else {
                        cerr << c_message<cERROR>(s[i]) << " is no valid option (see 'fconv -h')" << endl;
                        exit(1);
                    }
                }
                else if (s[i][2] == 'f' && opt2) {
                    if (s[i].size() > 4) {
                        if (s[i][4] == 'p') targettype = "pdb";
                        else if (s[i][4] == 'm') targettype = "mol";
                    } else {
                        cerr << c_message<cERROR>(s[i]) << " is no valid option (see 'fconv -h')" << endl;
                        exit(1);
                    }
                }
                else if (s[i][2] == 'r' && opt2) {
                    if (s[i].size() > 4) {
                        if (s[i][3] == '2' ) {
                            dummy.assign(s[i],5,s[i].size()-5);
                            istringstream is;
                            is.str(dummy);
                            is >> alpha1;
                        } else if (s[i][3] == '3' ) {
                            dummy.assign(s[i],5,s[i].size()-5);
                            istringstream is;
                            is.str(dummy);
                            is >> alpha2;
                        } else {
                            dummy.assign(s[i],4,s[i].size()-4);
                            istringstream is;
                            is.str(dummy);
                            is >> cut_radius;
                            unit_cell_exp_rad = cut_radius;
                            cluster_threshold = cut_radius;
                            surr_radius = cut_radius;
                            max_displacement = cut_radius;
                        }
                    } else {
                        cerr << c_message<cERROR>(s[i]) << " is no valid option (see 'fconv -h')" << endl;
                        exit(1);
                    }
                }
                else if (s[i][2] == 'v' && opt2) {
                    if (s[i].size() > 4) {
                        dummy.assign(s[i],4,s[i].size()-4);
                        istringstream is;
                        is.str(dummy);
                        is >> min_buried;
                    } else {
                        cerr << c_message<cERROR>(s[i]) << " is no valid option (see 'fconv -h')" << endl;
                        exit(1);
                    }
                }
                else if (s[i][2] == 'n' && opt2) {
                    if (s[i].size() > 4) {
                        dummy.assign(s[i],4,s[i].size()-4);
                        istringstream is;
                        is.str(dummy);
                        is >> max_poc;
                        unit_cell_multi = max_poc;
                    } else {
                        cerr << c_message<cERROR>(s[i]) << " is no valid option (see 'fconv -h')" << endl;
                        exit(1);
                    }
                }
                else if (s[i][2] == 'm' && opt2) {
                    if (s[i].size() > 4) {
                        dummy.assign(s[i],4,s[i].size()-4);
                        istringstream is;
                        is.str(dummy);
                        is >> max_ring_members;
                    } else {
                        cerr << c_message<cERROR>(s[i]) << " is no valid option (see 'fconv -h')" << endl;
                        exit(1);
                    }
                }
                else if (s[i][2] == '-' && opt2) {
                    if (s[i].size() > 3) {
                        if (s[i][3] == 'd') debug_mode = true;
                        cerr << "debug mode switched on" << endl;
                    } else {
                        cerr << c_message<cERROR>(s[i]) << " is no valid option (see 'fconv -h')" << endl;
                        exit(1);
                    }
                }
                else {
                    cerr << c_message<cERROR>(s[i]) << " is no valid option (see 'fconv -h')" << endl;
                    exit(1);
                }
            }
        } else {
            if (opt2 == false) {
                //!hier den mode anhand der Dateiendung ermitteln:
                int hs = s[i].size();
                if (s[i][hs-4] == 'm' && s[i][hs-3] == 'o' &&
                    s[i][hs-2] == 'l' && s[i][hs-1] == '2') mode = "mol";
                else if (s[i][hs-3] == 'p' && s[i][hs-2] == 'd' &&
                         s[i][hs-1] == 'b') mode = "pdb";
                else if (s[i][hs-3] == 'P' && s[i][hs-2] == 'D' &&
                         s[i][hs-1] == 'B') mode = "pdb";
                else if (s[i][hs-3] == 'd' && s[i][hs-2] == 'l' &&
                         s[i][hs-1] == 'g') mode = "dlg";
                else if (s[i][hs-3] == 's' && s[i][hs-2] == 'd' &&
                         s[i][hs-1] == 'f') mode = "sdf";
                else if (s[i][hs-3] == 'e' && s[i][hs-2] == 'n' &&
                         s[i][hs-1] == 't') mode = "pdb";
                else if (s[i][hs-3] == 'E' && s[i][hs-2] == 'N' &&
                         s[i][hs-1] == 'T') mode = "pdb";
                else if (s[i][hs-3] == 'c' && s[i][hs-2] == 'i' &&
                         s[i][hs-1] == 'f') mode = "cif";
                else if (s[i][hs-3] == 'C' && s[i][hs-2] == 'I' &&
                         s[i][hs-1] == 'F') mode = "cif";
                else if (s[i][hs-5] == 'p' && s[i][hs-4] == 'd' && s[i][hs-3] == 'b' && s[i][hs-2] == 'q' &&
                         s[i][hs-1] == 't') {mode = "pdb"; pdbqt_mode = true;}
                else if (s[i][hs-5] == 'P' && s[i][hs-4] == 'D' && s[i][hs-3] == 'B' && s[i][hs-2] == 'Q' &&
                         s[i][hs-1] == 'T') {mode = "pdb"; pdbqt_mode = true;}
                else if (!rename_files) {
                    cerr << c_message<cERROR>(s[i]) << "  has no supported file type extension" << endl;
                    exit(1);
                }
                opt2 = true;
            }
            sourcefiles.push_back(s[i]);
        }
    }

    if (sourcefiles.size() == 0) {
        cerr << c_message<cERROR>("could not detect any sourcefile (see 'fconv -h')") << endl;
        exit(1);
    }
    
    STRUCTURE *instruct;
    PARSER *inparser;
    
    instruct = new STRUCTURE();
    inparser = new PARSER(instruct,verb,gethyd,getwat,parse_ligands,parse_comments);
    
    if (rename_files) {
        #if defined (_LINUX_OS)
        char const* dir_name = "./";
        string start_ext = "_";
        rename_by_dir(sourcefiles,dir_name,start_ext);
        exit(0);
        #else
        cerr << "sorry, but this options is only available under LINUX/UNIX systems" << endl;
        exit(1);
        #endif
    }
    
    stl_ptr<LIGAND> m2_lig;
    
    if ((merge_mol2 || add_mol2 || merge_mol2_in1mol) && mode == "mol") {
        if (targetfile != "X") tfile = targetfile;
        else {
            cerr << c_message<cERROR>("you must specify a target file  (--t=targetfile)") << endl;
            delete instruct;
            delete inparser;
            exit(1);
        }
        if (merge_mol2_in1mol) {
            m2_lig = new LIGAND();
        }
        if (add_mol2) {
            if (!(inparser->read_mol2(targetfile.c_str()))) {
                cerr << c_message<cERROR>("could not read '") << targetfile << "' as mol2 file" << endl;
                delete instruct;
                delete inparser;
                exit(1);
            }
        }
    }

    if (cluster && mode=="pdb") {
        vector<stl_ptr<STRUCTURE> > strvec;
        for (vector<string>::iterator fit=sourcefiles.begin(); fit!=sourcefiles.end(); ++fit) {
            stl_ptr<STRUCTURE> inst;
            stl_ptr<PARSER> inpa;
            inst = new STRUCTURE();
            inpa = new PARSER(&(*inst),verb,false,false,false,false);
            if (!(inpa->read_pdb((*fit).c_str(),pdbqt_mode))) { //!pdb einlesen
                if (pdbqt_mode) cerr << c_message<cERROR>("could not read '") << (*fit).c_str() << "' as pdbqt file" << endl;
                else cerr << c_message<cERROR>("could not read '") << (*fit).c_str() << "' as pdb file" << endl;
                inst.kill();
                inpa.kill();
                continue;
            }
            inst->name = *fit;
            strvec.push_back(inst);
            inpa.kill();
        }

        cout << "clustering " << strvec.size() << " proteins" << endl;
        HIERA_CLUST<stl_ptr<STRUCTURE> > myclust(strvec,protein_sim_function);

        myclust.cluster_complete_link(99999999.);

        if (sourcefile != "X") {
            string dendro_name = "dendrogram_" + sourcefile;
            cout << "writing dendrogram as " << dendro_name << endl;
            myclust.write_dendrogram(dendro_name.c_str(),get_pro_label);

            string matrix_name = "matrix_" + sourcefile;
            cout << "writing similarity matrix as " << matrix_name << endl;
            myclust.write_matrix(matrix_name.c_str(),get_pro_label);
        } else {
            string dendro_name = "dendrogram.svg";
            cout << "writing dendrogram as " << dendro_name << endl;
            myclust.write_dendrogram(dendro_name.c_str(),get_pro_label);

            string matrix_name = "matrix.svg";
            cout << "writing similarity matrix as " << matrix_name << endl;
            myclust.write_matrix(matrix_name.c_str(),get_pro_label);
        }

        for (vector<stl_ptr<STRUCTURE> >::iterator it=strvec.begin(); it!=strvec.end(); ++it) it->kill();

        exit(0);
    }


    //int file_counter = 1;
    for (vector<string>::iterator fit=sourcefiles.begin(); fit!=sourcefiles.end(); ++fit) {
        inparser->clear();
        bool must_append = false;
        if (mode == "pdb") {
            if (add_hyd && !full_convert) {
                cerr << c_message<cWARNING>("adding hydrogens is currently only possible in -f mode!") << endl;
            }
            if (!(inparser->read_pdb((*fit).c_str(),pdbqt_mode))) { //!pdb einlesen
                if (pdbqt_mode) cerr << c_message<cERROR>("could not read '") << (*fit).c_str() << "' as pdbqt file" << endl;
                else cerr << c_message<cERROR>("could not read '") << (*fit).c_str() << "' as pdb file" << endl;
                instruct->clear();
                continue;
            }
            if (rmsd) {
                if (sourcefile == "X") {
                    cerr << c_message<cERROR>("you have to specify a reference file by --s=your_ref.pdb") << endl;
                    exit(1);
                }
                STRUCTURE *rstructure;
                PARSER *rparser;
                rstructure = new STRUCTURE();
                rparser = new PARSER(rstructure,verb,true,true,false,parse_comments);
                if (!(rparser->read_pdb(sourcefile.c_str()))) { //!pdb einlesen
                    cerr << c_message<cERROR>("could not read '") << sourcefile << "' as pdb file" << endl;
                    delete rstructure;
                    delete rparser;
                    exit(1);
                }
                if (rstructure->protein == 0) {
                    cerr << c_message<cERROR>("found no protein residues in '") << sourcefile << endl;
                    delete rstructure;
                    delete rparser;
                    exit(1);
                }
                if (instruct->protein == 0) {
                    cerr << c_message<cERROR>("found no protein residues in '") << *fit << endl;
                    instruct->clear();
                    continue;
                }
                cout << *fit << " :" << "\n";
                instruct->align_by_protein(rstructure->protein,true,debug_mode);
                instruct->clear();
                delete rstructure;
                delete rparser;
                continue;
            }
            if (remove_peptidic) {
                if (instruct->protein == 0) {
                    cerr << c_message<cERROR>("found no protein in ") << *fit << endl;
                    instruct->clear();
                    continue;
                }
                map<int,int> chain_n_res;
                for (int ci=0; ci<int(instruct->protein->chains.size()); ++ci) {
                    int old_num = 0;
                    chain_n_res[ci] = 0;
                    for (atoms_vec at=instruct->protein->chains[ci]->atoms.begin();
                                   at!=instruct->protein->chains[ci]->atoms.end(); ++at) {
                        if ((*at)->res_number != old_num) {
                            old_num = (*at)->res_number;
                            chain_n_res[ci] += 1;
                        }
                    }
                }
                
                vector<stl_ptr<CHAIN> > chain_buf;
                for (map<int,int>::iterator cit=chain_n_res.begin(); cit!=chain_n_res.end(); ++cit) {
                    if (cit->second < 4) instruct->protein->chains[cit->first].kill();
                    else chain_buf.push_back(instruct->protein->chains[cit->first]);
                }
                instruct->protein->chains.clear();
                for (vector<stl_ptr<CHAIN> >::iterator cit=chain_buf.begin(); cit!=chain_buf.end(); ++cit) {
                    instruct->protein->chains.push_back(*cit);
                }
                
                if (sourcefiles.size() > 1) {
                    tfile = *fit;
                    if (pdbqt_mode) replace_ext(tfile,".pdbqt",".pdb");
                    if (targetfile != "X") {
                        ostringstream os;
                        os << "_" << targetfile << ".pdb";
                        replace_ext(tfile,".pdb",os.str());
                    }
                } else {
                    if (targetfile != "X") tfile = targetfile;
                    else {
                        tfile = *fit;
                        if (pdbqt_mode) replace_ext(tfile,".pdbqt",".pdb");
                    }
                }
                string_fu::remove_char(tfile);
                inparser->write_pdb(tfile.c_str());
                instruct->clear();
                continue;
            }
            if (shape_sim) {
                if (sourcefile == "X") {
                    cerr << c_message<cERROR>("you have to specify a reference file by --s=your_ref.pdb") << endl;
                    exit(1);
                }
                STRUCTURE *rstructure;
                PARSER *rparser;
                rstructure = new STRUCTURE();
                rparser = new PARSER(rstructure,0,false,false,true,false);
                if (!(rparser->read_pdb(sourcefile.c_str()))) {
                    cerr << c_message<cERROR>("could not read '") << sourcefile << "' as pdb file" << endl;
                    delete rstructure;
                    delete rparser;
                    exit(1);
                }
                if (rstructure->ligands.size() == 0) {
                    cerr << c_message<cERROR>("found no ligand in ") << sourcefile << endl;
                    delete rstructure;
                    delete rparser;
                    exit(1);
                }
                if (instruct->ligands.size() == 0) {
                    cerr << c_message<cWARNING>("found no ligand in ") << *fit << "  => skipped" << endl;
                    delete rstructure;
                    delete rparser;
                    instruct->clear();
                    continue;
                }
                
                vector< stl_ptr<ATOM> > ref;
                rstructure->protein->get_CA_poc_atoms(ref);
                pair<int,float> result;
                result = instruct->protein->get_similarity(ref);
                
                cout << *fit << " :  matched C-alphas = " << result.first << "  rmsd = " << result.second << endl;
                
                delete rstructure;
                delete rparser;
                instruct->clear();
                continue;
            }
            if (s_pattern) {
                if (sourcefile == "X") {
                    cerr << c_message<cERROR>("you have to specify a reference file by --s=your_ref.pdb") << endl;
                    exit(1);
                }
                STRUCTURE *rstructure;
                PARSER *rparser;
                rstructure = new STRUCTURE();
                rparser = new PARSER(rstructure,0,false,false,true,false);
                if (!(rparser->read_pdb(sourcefile.c_str()))) {
                    cerr << c_message<cERROR>("could not read '") << sourcefile << "' as pdb file" << endl;
                    delete rstructure;
                    delete rparser;
                    exit(1);
                }
                multimap<float,vector<stl_ptr<ATOM> > > mmap;
                vector<vec3d<float> > t1; vector<vec3d<float> > t2; vector<matrix<float> > rm;
                instruct->protein->search_spattern2(rstructure->protein,mmap,t1,t2,rm);
                set<int> rnums;
                int nnc = -1;
                
                for (multimap<float,vector<stl_ptr<ATOM> > >::iterator it=mmap.begin(); it!=mmap.end(); ++it) {
                    ++nnc;
                    float ctotal = 0.; float cequal = 0.;
                    for (vector<stl_ptr<ATOM> >::iterator jt=it->second.begin(); jt!=it->second.end(); ++jt) {
                        ctotal += 1.;
                        if (rnums.find((*jt)->res_number) != rnums.end()) cequal += 1.;
                    }
                    if ((cequal/ctotal) > 0.8) continue;
                    rnums.clear();
                    cout << "\n@ " << *fit << ":\n";
                    
                    ostringstream os;
                    os.width(6); os.setf(ios::fixed); os.precision(3);
                    os << "     [" << right << rm[nnc][0][0] << "  ";
                    os.width(6); os << right << rm[nnc][0][1] << "  ";
                    os.width(6); os << right << rm[nnc][0][2] << "]   [    (";
                    os.width(6); os << right << t1[nnc][0] << ")]   (";
                    os.width(6); os << right << t2[nnc][0] << ")\n";
                    
                    os.width(6); os << "x' = [" << right << rm[nnc][1][0] << "  ";
                    os.width(6); os << right << rm[nnc][1][1] << "  ";
                    os.width(6); os << right << rm[nnc][1][2] << "] * [x - (";
                    os.width(6); os << right << t1[nnc][1] << ")] + (";
                    os.width(6); os << right << t2[nnc][1] << ")\n";
                    
                    os.width(6); os << "     [" << right << rm[nnc][2][0] << "  ";
                    os.width(6); os << right << rm[nnc][2][1] << "  ";
                    os.width(6); os << right << rm[nnc][2][2] << "]   [    (";
                    os.width(6); os << right << t1[nnc][2] << ")]   (";
                    os.width(6); os << right << t2[nnc][2] << ")" << endl;
                    os.unsetf(ios::fixed);
                    cout << os.str();
                    
                    cout << "rmsd = ";
                    cout << it->first << "   ";
                    cout << it->second.size() << " C-alpha matched:" << "\n";
                    for (vector<stl_ptr<ATOM> >::iterator jt=it->second.begin(); jt!=it->second.end(); ++jt) {
                        rnums.insert((*jt)->res_number);
                        cout << "   " << (*jt)->id << "  " << (*jt)->name << "  " << (*jt)->chain_id << "  " 
                             << (*jt)->res_name << "  " << (*jt)->res_number << "\n";
                    }
                }
                
                delete rstructure;
                delete rparser;
                instruct->clear();
                continue;
            }
            if (get_rot_bonds) {
            //    cout << "Calculating free rotatable bonds for '" << (*fit).c_str() << "':" << "\n";
                instruct->full_pdb2mol2(atom_typing,true,def_file.c_str(),true,
                                        verb,true,false,"X",max_ring_members,
                                        true,0,0,prot_acids,prot_guanidin,prot_amidin,prot_amin,
                                        prot_phosphate,prot_sulfate);
                int rotis;
                for (ligands_vec lit=instruct->mol2_protein.begin(); lit!=instruct->mol2_protein.end(); ++lit) {
                    rotis = (*lit)->get_freerot_bonds();
                    cout.width(30); cout << left << (*lit)->name;
                    cout << ":  "; cout.width(6); cout << left << rotis;
                    cout << "('" << (*fit).c_str() << "')" << "\n";
                }
                instruct->clear();
                continue;
            }
            if (type_count) {
                instruct->full_pdb2mol2(atom_typing,true,def_file.c_str(),true,
                                        verb,true,false,"X",max_ring_members,
                                        true,0,0,prot_acids,prot_guanidin,prot_amidin,prot_amin,
                                        prot_phosphate,prot_sulfate);
                for (ligands_vec lit=instruct->mol2_protein.begin(); lit!=instruct->mol2_protein.end(); ++lit) {
                    map<string,int> cmap;
                    for (atoms_vec at=(*lit)->atoms.begin(); at!=(*lit)->atoms.end(); ++at) {
                        if ((*at)->sybyl_type == "X") continue;
                        if (cmap.find((*at)->sybyl_type) == cmap.end()) cmap[(*at)->sybyl_type] = 1;
                        else cmap[(*at)->sybyl_type] += 1;
                    }
                    cout << (*lit)->name << ":" << "\n";
                    for (map<string,int>::iterator it=cmap.begin(); it!=cmap.end(); ++it) {
                        cout << "   -> ";
                        cout.width(8); cout << left << it->first << " : " << it->second << "\n";
                    }
                }
                instruct->clear();
                continue;
            }
            if (extract_ligs_seperate) {
                if (instruct->ligands.size() == 0) {
                    cerr << c_message<cERROR>("no ligands found in ") << *fit << endl;
                    instruct->clear();
                    continue;
                }
                
                instruct->get_alt_loc_ligs();
                
                if (instruct->has_peptides) instruct->merge_peptides(); //!Peptidketten mit zum Liganden
                
                instruct->merge_ligands(); //! 28.07.10: immer checken (ACHTUNG: muss NACH merge_peptides)

                int min_heavy_atoms_to_extract = 6;
                if (unit_cell_multi > 0) min_heavy_atoms_to_extract = unit_cell_multi;

                for (ligands_vec lig=instruct->ligands.begin(); lig!=instruct->ligands.end(); ++lig) {
                    if (sourcefile != "X") {
                        if ((*lig)->atoms[0]->res_name != sourcefile) continue;
                    }
                    
                ///    if ((*lig)->atoms.size() < 24) {
                        (*lig)->get_elements();
                        int n_heavy = 0;
                        bool has_carbon = false;
                        for (atoms_vec at=(*lig)->atoms.begin(); at!=(*lig)->atoms.end(); ++at) {
                            if ((*at)->element == "C") {
                                ++n_heavy;
                                has_carbon = true;
                            } else if ((*at)->element == "N" ||
                                       (*at)->element == "S" || (*at)->element == "O" ||
                                       (*at)->element == "P") ++n_heavy;
                        }
                        if (n_heavy < min_heavy_atoms_to_extract || !has_carbon) continue;
                ///    }

                    STRUCTURE* lpstruct = 0;
                    PARSER* lpparser = 0;
                    if (targettype == "pdb") {
                        lpstruct = new STRUCTURE();
                        lpparser = new PARSER(lpstruct,verb);
                        stl_ptr<LIGAND> lplig = new LIGAND(**lig);
                        lpstruct->ligands.push_back(lplig);
                    } else (*lig)->get_atom_typing(atom_typing,true,def_file.c_str(),true,verb,true,false,"X",max_ring_members,
                                                   true,0,0,prot_acids,prot_guanidin,prot_amidin,prot_amin,prot_phosphate,
                                                   prot_sulfate,kekulize_aromatics,kekulize_charged,allow_charged_aromatics);

                    string filename = "_" + (*lig)->name;
                    if ((*lig)->atoms[0]->chain_id != ' ') filename = filename + "_" + (*lig)->atoms[0]->chain_id; //! neu 07.07.2010
                    tfile = *fit;
                    if (pdbqt_mode) replace_ext(tfile,".pdbqt",".pdb");
                    replace_ext(tfile,".pdb",filename);
                    if (targetfile != "X") {
                        if (sourcefiles.size() == 1 && instruct->ligands.size() == 1) tfile = targetfile;
                        else {
                            tfile = targetfile;
                            replace_ext(tfile,"",filename);
                            if (targettype == "pdb") replace_ext(tfile,"",".pdb");
                            else replace_ext(tfile,"",".mol2");
                        }
                    } else {
                        if (targettype == "pdb") replace_ext(tfile,"",".pdb");
                        else replace_ext(tfile,"",".mol2");
                    }
                    
                    string_fu::remove_char(tfile);
                    if (targettype == "pdb") {
                        lpparser->write_pdb(tfile.c_str());
                        if (lpstruct) delete lpstruct;
                        if (lpparser) delete lpparser;
                    } else inparser->write_mol2_ligand(*lig,tfile.c_str());
                    if (verb > 0) cout << (*lig)->name << "  written to " << tfile << endl;
                }
                instruct->clear();
                continue;
            }
            if (bvals) {
                if (sourcefile != "X") { //Ligand extra mitgegeben
                    STRUCTURE *mstruct;
                    PARSER *mparser;
                    mstruct = new STRUCTURE();
                    mparser = new PARSER(mstruct,verb,true,true,true,parse_comments);
                    if (sourcefile[sourcefile.size()-1] == '2') {
                        if (!(mparser->read_mol2(sourcefile.c_str()))) {
                            cerr << c_message<cERROR>("could not read '") << sourcefile.c_str() << "' as mol2 file" << endl;
                            delete mstruct;
                            delete mparser;
                            break;
                        }
                    } else {
                        if (!(mparser->read_pdb(sourcefile.c_str()))) {
                            cerr << c_message<cERROR>("could not read '") << sourcefile.c_str() << "' as pdb file" << endl;
                            delete mstruct;
                            delete mparser;
                            break;
                        }
                        if (mstruct->ligands.size() == 0) {
                            cerr << c_message<cERROR>("found no ligands in '") << sourcefile.c_str() << endl;
                            delete mstruct;
                            delete mparser;
                            break;
                        }
                    }
                    for (ligands_vec lig=mstruct->ligands.begin(); lig!=mstruct->ligands.end(); ++lig) {
                        string ligname = (*lig)->name;
                    
                        float sdist = cut_radius * cut_radius;
                        vector<stl_ptr<ATOM> > bats;                    
                        for (chains_vec cv=instruct->protein->chains.begin(); cv!=instruct->protein->chains.end();++cv) {
                            tr1::unordered_set<int> taken_bats;
                            tr1::unordered_set<int> taken_res;
                            for (atoms_vec at=(*cv)->atoms.begin(); at!=(*cv)->atoms.end(); ++at) {
                                for (atoms_vec bt=(*lig)->atoms.begin(); bt!=(*lig)->atoms.end(); ++bt) {
                                    if (get_square_distance((*at)->coord,(*bt)->coord) <= sdist) {
                                        bats.push_back(*at);
                                        if (cut_full_res) {
                                            taken_bats.insert((*at)->intern_id);
                                            taken_res.insert((*at)->res_number);
                                        }
                                        break;
                                    }
                                }
                            }
                            if (cut_full_res) {
                                for (atoms_vec at=(*cv)->atoms.begin(); at!=(*cv)->atoms.end(); ++at) {
                                    if (taken_bats.find((*at)->intern_id) != taken_bats.end()) continue;
                                    if (taken_res.find((*at)->res_number) != taken_res.end()) {
                                        bats.push_back(*at);
                                    }
                                }
                            }
                        }
                        if (bats.size() < 2) continue;
                        float bval_sum = 0.;
                        float bval_mean = 0.;
                        float bval_sdev = 0.;
                        float bval_sumb = 0.; int back = 0;
                        float bval_meanb = 0.;
                        float bval_sdevb = 0.;
                        float bval_sumr = 0.; int res = 0;
                        float bval_meanr = 0.;
                        float bval_sdevr = 0.;
                        for (atoms_vec at=bats.begin(); at!=bats.end(); ++at) {
                            bval_sum += (*at)->b_factor;
                            
                            istringstream is;
                            string aname;
                            is.str((*at)->name);
                            is >> aname;
                            if (aname == "N" || aname == "CA" || aname == "H" || aname == "C" ||
                                aname == "O" || aname == "HA" || aname == "HO" || aname == "HN") {
                                   bval_sumb += (*at)->b_factor; ++back;                 
                            } else {bval_sumr += (*at)->b_factor; ++res;}
                        }
                        bval_mean = bval_sum / bats.size();
                        if (back > 0) bval_meanb = bval_sumb / back;
                        if (res > 0) bval_meanr = bval_sumr / res;                    
                       
                        for (atoms_vec at=bats.begin(); at!=bats.end(); ++at) {
                            float tval = (*at)->b_factor - bval_mean;
                            tval *= tval;                       
                            bval_sdev += tval;
                            istringstream is;
                            string aname;
                            is.str((*at)->name);
                            is >> aname;                        
                            if (aname == "N" || aname == "CA" || aname == "H" || aname == "C" ||
                                aname == "O" || aname == "HA" || aname == "HO" || aname == "HN") {
                                tval = (*at)->b_factor - bval_meanb;
                                tval *= tval;                       
                                bval_sdevb += tval;    
                            } else {
                                tval = (*at)->b_factor - bval_meanr;
                                tval *= tval;                       
                                bval_sdevr += tval; 
                            }
                        }
                        bval_sdev /= bats.size() - 1.;
                        if (back > 1) bval_sdevb /= back - 1.;
                        if (res > 1) bval_sdevr /= res - 1.;

                        cout << "-> ligname: " << ligname << "\n"
                             << "    protein atoms in " << cut_radius << "A distance = " << bats.size() << "\n"
                             << "    mean b-value                 = " << bval_mean << "\n"
                             << "    standard deviation           = " << sqrt(bval_sdev) << "\n"
                             << "    backbone atoms               = " << back << "\n"
                             << "    mean backbone b-value        = " << bval_meanb << "\n"
                             << "    backbone standard deviation  = " << sqrt(bval_sdevb) << "\n"
                             << "    residue atoms                = " << res << "\n"
                             << "    mean residue b-value         = " << bval_meanr << "\n"
                             << "    residue standard deviation   = " << sqrt(bval_sdevr) << "\n" << endl; 
                    }
                    delete mstruct;
                    delete mparser;
                } else {                
                    if (instruct->ligands.size() == 0) {
                        cerr << c_message<cERROR>("no ligands found in ") << *fit << endl;
                        instruct->clear();
                        continue;
                    }
                
                    instruct->get_alt_loc_ligs();
                    if (instruct->has_peptides) instruct->merge_peptides(); //!Peptidketten mit zum Liganden
                    instruct->merge_ligands(); //! 28.07.10: immer checken (ACHTUNG: muss NACH merge_peptides)
                
                    for (ligands_vec lig=instruct->ligands.begin(); lig!=instruct->ligands.end(); ++lig) {
                        (*lig)->get_elements();
                        int n_heavy = 0;
                        bool has_carbon = false;
                        for (atoms_vec at=(*lig)->atoms.begin(); at!=(*lig)->atoms.end(); ++at) {
                            if ((*at)->element == "C") {
                                ++n_heavy;
                                has_carbon = true;
                            } else if ((*at)->element == "N" ||
                                       (*at)->element == "S" || (*at)->element == "O" ||
                                       (*at)->element == "P") ++n_heavy;
                        }
                        if (n_heavy < 6 || !has_carbon) continue;
                    
                        string ligname = (*lig)->name + (*lig)->atoms[0]->chain_id;
                    
                        float sdist = cut_radius * cut_radius;
                        vector<stl_ptr<ATOM> > bats;                    
                        for (chains_vec cv=instruct->protein->chains.begin(); cv!=instruct->protein->chains.end();++cv) {
                            tr1::unordered_set<int> taken_bats;
                            tr1::unordered_set<int> taken_res;
                            for (atoms_vec at=(*cv)->atoms.begin(); at!=(*cv)->atoms.end(); ++at) {
                                for (atoms_vec bt=(*lig)->atoms.begin(); bt!=(*lig)->atoms.end(); ++bt) {
                                    if (get_square_distance((*at)->coord,(*bt)->coord) <= sdist) {
                                        bats.push_back(*at);
                                        if (cut_full_res) {
                                            taken_bats.insert((*at)->intern_id);
                                            taken_res.insert((*at)->res_number);
                                        }
                                        break;
                                    }
                                }
                            }
                            if (cut_full_res) {
                                for (atoms_vec at=(*cv)->atoms.begin(); at!=(*cv)->atoms.end(); ++at) {
                                    if (taken_bats.find((*at)->intern_id) != taken_bats.end()) continue;
                                    if (taken_res.find((*at)->res_number) != taken_res.end()) {
                                        bats.push_back(*at);
                                    }
                                }
                            }
                        }
                        if (bats.size() < 2) continue;
                        float bval_sum = 0.;
                        float bval_mean = 0.;
                        float bval_sdev = 0.;
                        float bval_sumb = 0.; int back = 0;
                        float bval_meanb = 0.;
                        float bval_sdevb = 0.;
                        float bval_sumr = 0.; int res = 0;
                        float bval_meanr = 0.;
                        float bval_sdevr = 0.;
                        for (atoms_vec at=bats.begin(); at!=bats.end(); ++at) {
                            bval_sum += (*at)->b_factor;
                            
                            istringstream is;
                            string aname;
                            is.str((*at)->name);
                            is >> aname;
                            if (aname == "N" || aname == "CA" || aname == "H" || aname == "C" ||
                                aname == "O" || aname == "HA" || aname == "HO" || aname == "HN") {
                                   bval_sumb += (*at)->b_factor; ++back;                 
                            } else {bval_sumr += (*at)->b_factor; ++res;}
                        }
                        bval_mean = bval_sum / bats.size();
                        if (back > 0) bval_meanb = bval_sumb / back;
                        if (res > 0) bval_meanr = bval_sumr / res;                    
                       
                        for (atoms_vec at=bats.begin(); at!=bats.end(); ++at) {
                            float tval = (*at)->b_factor - bval_mean;
                            tval *= tval;                       
                            bval_sdev += tval;
                            istringstream is;
                            string aname;
                            is.str((*at)->name);
                            is >> aname;                            
                            if (aname == "N" || aname == "CA" || aname == "H" || aname == "C" ||
                                aname == "O" || aname == "HA" || aname == "HO" || aname == "HN") {
                                tval = (*at)->b_factor - bval_meanb;
                                tval *= tval;                       
                                bval_sdevb += tval;    
                            } else {
                                tval = (*at)->b_factor - bval_meanr;
                                tval *= tval;                       
                                bval_sdevr += tval; 
                            }
                        }
                        bval_sdev /= bats.size() - 1.;
                        if (back > 1) bval_sdevb /= back - 1.;
                        if (res > 1) bval_sdevr /= res - 1.;

                        cout << "-> ligname: " << ligname << "\n"
                             << "    protein atoms in " << cut_radius << "A distance = " << bats.size() << "\n"
                             << "    mean b-value                 = " << bval_mean << "\n"
                             << "    standard deviation           = " << sqrt(bval_sdev) << "\n"
                             << "    backbone atoms               = " << back << "\n"
                             << "    mean backbone b-value        = " << bval_meanb << "\n"
                             << "    backbone standard deviation  = " << sqrt(bval_sdevb) << "\n"
                             << "    residue atoms                = " << res << "\n"
                             << "    mean residue b-value         = " << bval_meanr << "\n"
                             << "    residue standard deviation   = " << sqrt(bval_sdevr) << "\n" << endl;  
                    }
                }
                instruct->clear();
                continue;
            }
            if (extract_water) {
                instruct->water2lig();
                instruct->water_ligand->get_atom_typing(atom_typing,true,def_file.c_str(),true,verb,true,false,"X",max_ring_members,
                                                   true,0,0,prot_acids,prot_guanidin,prot_amidin,prot_amin,prot_phosphate,
                                                   prot_sulfate,kekulize_aromatics,kekulize_charged,allow_charged_aromatics);
                if (targetfile != "X") tfile = targetfile;
                else {
                    ostringstream os;
                    os << "_water" << ".mol2";
                    if (targetfile != "X") tfile = targetfile;
                    else {
                        tfile = *fit;
                        if (pdbqt_mode) replace_ext(tfile,".pdbqt",".pdb");
                        replace_ext(tfile,".pdb",os.str());
                    }
                }
                if (instruct->water_ligand->atoms.size() > 0) {
                    string_fu::remove_char(tfile);
                    inparser->write_mol2_ligand(instruct->water_ligand,tfile.c_str());
                    if (verb > 0) cout << instruct->water_ligand->atoms.size() 
                                       << " waters from " << *fit << " written to " << tfile << endl;
                } else if (verb > 0) cout << "no waters found in " << *fit << endl;
                instruct->clear();
                continue;
            }
            if (extract_metals) {
                
                instruct->metal2lig();
                instruct->metal_ligand->get_atom_typing(atom_typing,true,def_file.c_str(),true,verb,true,false,"X",max_ring_members,
                                                   true,0,0,prot_acids,prot_guanidin,prot_amidin,prot_amin,prot_phosphate,
                                                   prot_sulfate,kekulize_aromatics,kekulize_charged,allow_charged_aromatics);
                if (targetfile != "X") tfile = targetfile;
                else {
                    ostringstream os;
                    os << "_metals" << ".mol2";
                    if (targetfile != "X") tfile = targetfile;
                    else {
                        tfile = *fit;
                        if (pdbqt_mode) replace_ext(tfile,".pdbqt",".pdb");
                        replace_ext(tfile,".pdb",os.str());
                    }
                }
                // setzt Atomtypen der Liganden neu, geht Liganden und Atome der Liganden durch
                for (ligands_vec lig=instruct->ligands.begin(); lig!=instruct->ligands.end(); ++lig) {
                    (*lig)->get_atom_typing(atom_typing,true,def_file.c_str(),true,verb,true,false,"X",max_ring_members,
                                                   true,0,0,prot_acids,prot_guanidin,prot_amidin,prot_amin,prot_phosphate,
                                                   prot_sulfate,kekulize_aromatics,kekulize_charged,allow_charged_aromatics);
                    for (atoms_vec at=(*lig)->atoms.begin(); at!=(*lig)->atoms.end(); ++at) {
                        if ((*at)->type == 3) {
                            //instruct->metal_ligand->atoms.push_back(*at);
                            instruct->metal_ligand->atoms.push_back(new ATOM(**at));
                        }
                    }
                }
                if (instruct->metal_ligand->atoms.size() > 0) {
                    string_fu::remove_char(tfile);
                    inparser->write_mol2_ligand(instruct->metal_ligand,tfile.c_str());
                    if (verb > 0) cout << instruct->metal_ligand->atoms.size() 
                                       << " metals from " << *fit << " written to " << tfile << endl;
                } else if (verb > 0) cout << "no metals found in " << *fit << endl;
                
                instruct->clear();
                
                continue;
            }
            if (extract_all_ligs_cav) {
                if (sourcefile != "X") { //Ligand extra mitgegeben
                    STRUCTURE *mstruct;
                    PARSER *mparser;
                    mstruct = new STRUCTURE();
                    mparser = new PARSER(mstruct,verb,true,true,true,parse_comments);
                    if (sourcefile[sourcefile.size()-1] == '2') {
                        if (!(mparser->read_mol2(sourcefile.c_str()))) {
                            cerr << c_message<cERROR>("could not read '") << sourcefile.c_str() << "' as mol2 file" << endl;
                            delete mstruct;
                            delete mparser;
                            break;
                        }
                    } else {
                        if (!(mparser->read_pdb(sourcefile.c_str()))) {
                            cerr << c_message<cERROR>("could not read '") << sourcefile.c_str() << "' as pdb file" << endl;
                            delete mstruct;
                            delete mparser;
                            break;
                        }
                        if (mstruct->ligands.size() == 0) {
                            cerr << c_message<cERROR>("found no ligands in '") << sourcefile.c_str() << endl;
                            delete mstruct;
                            delete mparser;
                            break;
                        }
                    }
                    for (ligands_vec lig=mstruct->ligands.begin(); lig!=mstruct->ligands.end(); ++lig) {
                        /*
                        if ((*lig)->atoms.size() < 24) {
                            int n_hyd = 0;
                            for (atoms_vec at=(*lig)->atoms.begin(); at!=(*lig)->atoms.end(); ++at) {
                                if ((*at)->element == "H") n_hyd++;
                            }
                            if (((*lig)->atoms.size() - n_hyd) < 6) continue;
                        }
                        */
                        string filename;
                        if (targettype == "mol") filename = "_" + (*lig)->name + "_poc" + ".mol2";
                        else filename = "_" + (*lig)->name + "_poc" + ".pdb";
                        if (targetfile != "X") tfile = targetfile;
                        else {
                            tfile = *fit;
                            if (pdbqt_mode) replace_ext(tfile,".pdbqt",".pdb");
                            replace_ext(tfile,".pdb",filename);
                        }
                        instruct->get_lig_cavity(cut_radius,&(**lig),cut_full_res,no_hyd_for_poc);
                        string_fu::remove_char(tfile);
                        if (cut_full_res) inparser->write_pdb_cav(tfile.c_str(),instruct->cavities.back());
                        else inparser->write_pdb_cav_atoms_only(tfile.c_str());
                        if (verb > 0) cout << (*lig)->name << " pocket written to " << tfile << endl;
                        if (targettype == "mol") {
                            STRUCTURE *tstruct;
                            PARSER *tparser;
                            tstruct = new STRUCTURE();
                            tparser = new PARSER(tstruct,0,true,true,true,parse_comments);
                            tparser->read_pdb(tfile.c_str());
                            tstruct->pdb2mol2(atom_typing,def_file.c_str());
                            tparser->write_mol2_protein(tstruct->protein,tfile.c_str());
                            delete tstruct;
                            delete tparser;
                        }
                    }
                    delete mstruct;
                    delete mparser;
                    instruct->clear();
                    continue;
                }
                //Liganden aus pdb nehmen
                if (instruct->has_peptides) instruct->merge_peptides();
                for (ligands_vec lig=instruct->ligands.begin(); lig!=instruct->ligands.end(); ++lig) {
                    if ((*lig)->atoms.size() < 24) {
                        int n_hyd = 0;
                        for (atoms_vec at=(*lig)->atoms.begin(); at!=(*lig)->atoms.end(); ++at) {
                            if ((*at)->element == "H") n_hyd++;
                        }
                        if (((*lig)->atoms.size() - n_hyd) < 6) continue;
                    }
                    
                    string filename;
                    if (targettype == "mol") filename = "_" + (*lig)->name + "_poc" + ".mol2";
                    else filename = "_" + (*lig)->name + "_poc" + ".pdb";
                    if (targetfile != "X") tfile = targetfile;
                    else {
                        tfile = *fit;
                        if (pdbqt_mode) replace_ext(tfile,".pdbqt",".pdb");
                    }
                    replace_ext(tfile,".pdb",filename);
                    instruct->get_lig_cavity(cut_radius,&(**lig),cut_full_res,no_hyd_for_poc);
                    string_fu::remove_char(tfile);
                    inparser->write_pdb_cav(tfile.c_str(),instruct->cavities.back());
                    if (verb > 0) cout << (*lig)->name << " pocket written to " << tfile << endl;
                    if (targettype == "mol") {
                        STRUCTURE *mstruct;
                        PARSER *mparser;
                        mstruct = new STRUCTURE();
                        mparser = new PARSER(mstruct,0,true,true,true,parse_comments);
                        mparser->read_pdb(tfile.c_str());
                        mstruct->pdb2mol2(atom_typing,def_file.c_str());
                        mparser->write_mol2_protein(mstruct->protein,tfile.c_str());
                        delete mstruct;
                        delete mparser;
                    }
                }
                instruct->clear();
                continue;
            }
            else if (cut_by_lig) {
                if (instruct->ligands.size() < (cut_lig_number + 1)) {
                    cerr << c_message<cERROR>((*fit).c_str()) << " has only " << instruct->ligands.size() << " ligands" << endl;
                    instruct->clear();
                    continue;
                }
                if (targetfile != "X") tfile = targetfile;
                else {
                    string filename;
                    if (targettype == "mol") filename = "_" + instruct->ligands[cut_lig_number]->name + "_poc" + ".mol2";
                    else filename = "_" + instruct->ligands[cut_lig_number]->name + "_poc" + ".pdb";
                    if (targetfile != "X") tfile = targetfile;
                    else {
                        tfile = *fit;
                        if (pdbqt_mode) replace_ext(tfile,".pdbqt",".pdb");
                    }
                    replace_ext(tfile,".pdb",filename);
                }
                instruct->get_lig_cavity(cut_radius,&(*(instruct->ligands[cut_lig_number])),cut_full_res,no_hyd_for_poc);
                string_fu::remove_char(tfile);
                inparser->write_pdb_cav(tfile.c_str(),instruct->cavities.back());
                if (verb > 0) cout << instruct->ligands[cut_lig_number]->name << " pocket written to " << tfile << endl;
                if (targettype == "mol") {
                    STRUCTURE *mstruct;
                    PARSER *mparser;
                    mstruct = new STRUCTURE();
                    mparser = new PARSER(mstruct,0,true,true,true,parse_comments);
                    mparser->read_pdb(tfile.c_str());
                    mstruct->pdb2mol2(atom_typing,def_file.c_str());
                    mparser->write_mol2_protein(mstruct->protein,tfile.c_str());
                    delete mstruct;
                    delete mparser;
                }    
                instruct->clear();
                continue;
            }
            else if (auto_cav) {
                instruct->get_cavity_auto(cut_radius,min_buried,debug_mode,max_poc);
                unsigned int poccount = 1;
                for (cavities_vec cav=instruct->cavities.begin(); cav!=instruct->cavities.end(); ++cav) {
                    string filename;
                    ostringstream hos;
                    hos << poccount;
                    if (targettype == "mol") filename = "_poc_" + hos.str() + ".mol2";
                    else filename = "_poc" + hos.str() + ".pdb";
                    if (targetfile != "X") tfile = targetfile;
                    else {
                        tfile = *fit;
                        if (pdbqt_mode) replace_ext(tfile,".pdbqt",".pdb");
                    }
                    replace_ext(tfile,".pdb",filename);
                    string_fu::remove_char(tfile);
                    inparser->write_pdb_cav(tfile.c_str(),*cav);
                    if (verb > 0) cout << (*cav)->total_volume << " A^3 pocket written to " << tfile << endl;
                    if (targettype == "mol") {
                        STRUCTURE *mstruct;
                        PARSER *mparser;
                        mstruct = new STRUCTURE();
                        mparser = new PARSER(mstruct,0,true,true,true,parse_comments);
                        mparser->read_pdb(tfile.c_str());
                        mstruct->pdb2mol2(atom_typing,def_file.c_str());
                        mparser->write_mol2_protein(mstruct->protein,tfile.c_str());
                        delete mstruct;
                        delete mparser;
                    }
                    poccount++;
                    if (poccount > max_poc) break;
                }
                instruct->clear();
                continue;
            }
            else if (auto_cav2) {
                vector<stl_ptr<ATOM> > pv;
                for (chains_vec ct=instruct->protein->chains.begin(); ct!=instruct->protein->chains.end(); ++ct) {
                    (*ct)->get_atom_typing(1,true,"X",true,0,false,false,"X",max_ring_members,
                                           true,0,0,prot_acids,prot_guanidin,prot_amidin,prot_amin,prot_phosphate,
                                           prot_sulfate,kekulize_aromatics,kekulize_charged,allow_charged_aromatics);
                    for (atoms_vec at=(*ct)->atoms.begin(); at!=(*ct)->atoms.end(); ++at) {
                        if ((*at)->alt_loc_id != ' ' && (*at)->alt_loc_id != 'A') continue; //!sonst unter Umstaenden zweimal
                                                                                            //!gleiche Coords!!!
                        pv.push_back(*at);
                    }
                }
                STRUCTURE *rstructure = 0;
                PARSER *rparser = 0;
                if (pdb_source != "X") {
                    rstructure = new STRUCTURE();
                    rparser = new PARSER(rstructure,verb,false,true,true,false);
                    if (!(rparser->read_mol2(pdb_source.c_str(),false,true))) {
                        cerr << c_message<cERROR>("could not load ") << pdb_source << " as mol2 file" << endl;
                        delete rstructure;
                        delete rparser;
                        exit(1);
                    }
                    for (ligands_vec lig=rstructure->ligands.begin(); lig!=rstructure->ligands.end(); ++lig) {
                        (*lig)->get_atom_typing(1,true,"X",true,0,false,false,"X",max_ring_members,
                                           true,0,0,prot_acids,prot_guanidin,prot_amidin,prot_amin,prot_phosphate,
                                           prot_sulfate,kekulize_aromatics,kekulize_charged,allow_charged_aromatics);
                        for (atoms_vec at=(*lig)->atoms.begin(); at!=(*lig)->atoms.end(); ++at) pv.push_back(*at);
                    }
                }
                vector<stl_ptr<ATOM> > lig_atoms;
                STRUCTURE *lig_struct = 0;
                PARSER *lig_parser = 0;
                if (sourcefile != "X") {
                    lig_struct = new STRUCTURE();
                    lig_parser = new PARSER(lig_struct,1,false,false,true,false);
                    if (!(lig_parser->read_mol2(sourcefile.c_str(),false,true))) {
                        cerr << c_message<cERROR>("could not load ") << sourcefile << " as mol2 file" << endl;
                        delete lig_struct;
                        delete lig_parser;
                        exit(1);
                    }
                    for (ligands_vec lig=lig_struct->ligands.begin(); lig!=lig_struct->ligands.end(); ++lig) {
                        (*lig)->get_atom_typing(1,true,"X",true,0,false,false,"X",max_ring_members,
                                           true,0,0,prot_acids,prot_guanidin,prot_amidin,prot_amin,prot_phosphate,
                                           prot_sulfate,kekulize_aromatics,kekulize_charged,allow_charged_aromatics);
                        for (atoms_vec at=(*lig)->atoms.begin(); at!=(*lig)->atoms.end(); ++at) lig_atoms.push_back(*at);
                    }
                }
                map<int,vector<stl_ptr<ATOM> > > p_clust;
                map<int,double> output_vol;
                instruct->get_alpha_cavity(pv,p_clust,output_vol,lig_atoms,alpha1,alpha2,true,300.,
                                           debug_mode,4.0,surr_radius,false);
                
                if (p_clust.size() == 0) {
                    cout << "--> no cavities found => searching with lower minimum volume:" << endl;
                    instruct->get_alpha_cavity(pv,p_clust,output_vol,lig_atoms,alpha1,alpha2,true,200.,
                                               debug_mode,4.0,surr_radius,false);
                }
                if (p_clust.size() == 0) {
                    cout << "--> again no cavities found => searching with higher alpha2 level:" << endl;
                    alpha2 += 2.;
                    instruct->get_alpha_cavity(pv,p_clust,output_vol,lig_atoms,alpha1,alpha2,true,200.,
                                               debug_mode,4.0,surr_radius,false);
                }
                
                for (map<int,vector<stl_ptr<ATOM> > >::iterator it=p_clust.begin(); it!=p_clust.end(); ++it) {
                    cout << "--> volume of cavity " << it->first << " = " << output_vol[it->first] << " A^3" << endl;
                    stl_ptr<LIGAND> plig = new LIGAND();
                    for (atoms_vec at=it->second.begin(); at!=it->second.end(); ++at) plig->atoms.push_back(*at);
                    string n_ext = "_"; string_fu::add2s(n_ext,it->first); string_fu::add2s(n_ext,".mol2");
                    string lname = string_fu::get_text_after(*fit,'/');
                    string_fu::replace_ext(lname,".pdb",n_ext);
                    inparser->write_mol2_ligand(plig,lname.c_str());
                    plig->atoms.clear();
                    plig.kill();
                    cout << "    written cavity " << it->first << " as " << lname << endl;
                }
                if (sourcefile != "X") {
                    delete lig_struct;
                    delete lig_parser;
                }
                if (pdb_source != "X") {
                    delete rstructure;
                    delete rparser;
                }
                instruct->clear();
                continue;
            }
            else if (alpha_profile) {
                vector<stl_ptr<ATOM> > pv;
                for (chains_vec ct=instruct->protein->chains.begin(); ct!=instruct->protein->chains.end(); ++ct) {
                    (*ct)->get_atom_typing(1,true,"X",true,0,false,false,"X",max_ring_members,
                                           true,0,0,prot_acids,prot_guanidin,prot_amidin,prot_amin,prot_phosphate,
                                           prot_sulfate,kekulize_aromatics,kekulize_charged,allow_charged_aromatics);
                    for (atoms_vec at=(*ct)->atoms.begin(); at!=(*ct)->atoms.end(); ++at) {
                        if ((*at)->alt_loc_id != ' ' && (*at)->alt_loc_id != 'A') continue; //!sonst unter Umstaenden zweimal
                                                                                            //!gleiche Coords!!!
                        pv.push_back(*at);
                    }
                }
                STRUCTURE *rstructure = 0;
                PARSER *rparser = 0;
                if (pdb_source != "X") {
                    rstructure = new STRUCTURE();
                    rparser = new PARSER(rstructure,verb,false,true,true,false);
                    if (!(rparser->read_mol2(pdb_source.c_str(),false,true))) {
                        cerr << c_message<cERROR>("could not load ") << pdb_source << " as mol2 file" << endl;
                        delete rstructure;
                        delete rparser;
                        exit(1);
                    }
                    for (ligands_vec lig=rstructure->ligands.begin(); lig!=rstructure->ligands.end(); ++lig) {
                        (*lig)->get_atom_typing(1,true,"X",true,0,false,false,"X",max_ring_members,
                                           true,0,0,prot_acids,prot_guanidin,prot_amidin,prot_amin,prot_phosphate,
                                           prot_sulfate,kekulize_aromatics,kekulize_charged,allow_charged_aromatics);
                        for (atoms_vec at=(*lig)->atoms.begin(); at!=(*lig)->atoms.end(); ++at) pv.push_back(*at);
                    }
                }
                vector<stl_ptr<ATOM> > lig_atoms;
                STRUCTURE *lig_struct = 0;
                PARSER *lig_parser = 0;
                if (sourcefile != "X") {
                    lig_struct = new STRUCTURE();
                    lig_parser = new PARSER(lig_struct,1,false,false,true,false);
                    if (!(lig_parser->read_mol2(sourcefile.c_str(),false,true))) {
                        cerr << c_message<cERROR>("could not load ") << sourcefile << " as mol2 file" << endl;
                        delete lig_struct;
                        delete lig_parser;
                        exit(1);
                    }
                    for (ligands_vec lig=lig_struct->ligands.begin(); lig!=lig_struct->ligands.end(); ++lig) {
                        (*lig)->get_atom_typing(1,true,"X",true,0,false,false,"X",max_ring_members,
                                           true,0,0,prot_acids,prot_guanidin,prot_amidin,prot_amin,prot_phosphate,
                                           prot_sulfate,kekulize_aromatics,kekulize_charged,allow_charged_aromatics);
                        for (atoms_vec at=(*lig)->atoms.begin(); at!=(*lig)->atoms.end(); ++at) lig_atoms.push_back(*at);
                    }
                }
                float sr = surr_radius;
                if (sr < 0.001) sr = 0.1;
            //!    instruct->get_alpha_profile(pv,lig_atoms,debug_mode);
                cerr << "This option is currently disabled!" << endl;
                
                if (sourcefile != "X") {
                    delete lig_struct;
                    delete lig_parser;
                }
                if (pdb_source != "X") {
                    delete rstructure;
                    delete rparser;
                }
                instruct->clear();
                continue;
            }
            else if (auto_cav3) {
                vector<stl_ptr<ATOM> > pv;
                for (chains_vec ct=instruct->protein->chains.begin(); ct!=instruct->protein->chains.end(); ++ct) {
                    (*ct)->get_atom_typing(1,true,"X",true,0,false,false,"X",max_ring_members,
                                           true,0,0,prot_acids,prot_guanidin,prot_amidin,prot_amin,prot_phosphate,
                                           prot_sulfate,kekulize_aromatics,kekulize_charged,allow_charged_aromatics);
                    for (atoms_vec at=(*ct)->atoms.begin(); at!=(*ct)->atoms.end(); ++at) {
                        if ((*at)->alt_loc_id != ' ' && (*at)->alt_loc_id != 'A') continue; //!sonst unter Umstaenden zweimal
                                                                                            //!gleiche Coords!!!
                        pv.push_back(*at);
                    }
                }
                STRUCTURE *rstructure = 0;
                PARSER *rparser = 0;
                if (pdb_source != "X") {
                    rstructure = new STRUCTURE();
                    rparser = new PARSER(rstructure,verb,false,true,true,false);
                    if (!(rparser->read_mol2(pdb_source.c_str(),false,true))) {
                        cerr << c_message<cERROR>("could not load ") << pdb_source << " as mol2 file" << endl;
                        delete rstructure;
                        delete rparser;
                        exit(1);
                    }
                    for (ligands_vec lig=rstructure->ligands.begin(); lig!=rstructure->ligands.end(); ++lig) {
                        (*lig)->get_atom_typing(1,true,"X",true,0,false,false,"X",max_ring_members,
                                           true,0,0,prot_acids,prot_guanidin,prot_amidin,prot_amin,prot_phosphate,
                                           prot_sulfate,kekulize_aromatics,kekulize_charged,allow_charged_aromatics);
                        for (atoms_vec at=(*lig)->atoms.begin(); at!=(*lig)->atoms.end(); ++at) pv.push_back(*at);
                    }
                }
                vector<stl_ptr<ATOM> > lig_atoms;
                STRUCTURE *lig_struct = 0;
                PARSER *lig_parser = 0;
                if (sourcefile != "X") {
                    lig_struct = new STRUCTURE();
                    lig_parser = new PARSER(lig_struct,1,false,false,true,false);
                    if (!(lig_parser->read_mol2(sourcefile.c_str(),false,true))) {
                        cerr << c_message<cERROR>("could not load ") << sourcefile << " as mol2 file" << endl;
                        delete lig_struct;
                        delete lig_parser;
                        exit(1);
                    }
                    for (ligands_vec lig=lig_struct->ligands.begin(); lig!=lig_struct->ligands.end(); ++lig) {
                        (*lig)->get_atom_typing(1,true,"X",true,0,false,false,"X",max_ring_members,
                                           true,0,0,prot_acids,prot_guanidin,prot_amidin,prot_amin,prot_phosphate,
                                           prot_sulfate,kekulize_aromatics,kekulize_charged,allow_charged_aromatics);
                        for (atoms_vec at=(*lig)->atoms.begin(); at!=(*lig)->atoms.end(); ++at) lig_atoms.push_back(*at);
                    }
                }
                map<int,vector<stl_ptr<ATOM> > > p_clust;
                map<int,double> output_vol;
                instruct->get_alpha_cavity(pv,p_clust,output_vol,lig_atoms,alpha1,alpha2,true,300.,
                                           debug_mode,4.0,surr_radius,true);
                if (p_clust.size() == 0) {
                    cout << "--> no cavities found => searching with lower minimum volume:" << endl;
                    instruct->get_alpha_cavity(pv,p_clust,output_vol,lig_atoms,alpha1,alpha2,true,200.,
                                               debug_mode,4.0,surr_radius,false);
                }
                if (p_clust.size() == 0) {
                    cout << "--> again no cavities found => searching with higher alpha2 level:" << endl;
                    alpha2 += 2.;
                    instruct->get_alpha_cavity(pv,p_clust,output_vol,lig_atoms,alpha1,alpha2,true,200.,
                                               debug_mode,4.0,surr_radius,false);
                }
                
                int n_cavs = 0;
                for (cavities_vec cav=instruct->cavities.begin(); cav!=instruct->cavities.end(); ++cav) {
                    cout << "--> volume of cavity " << n_cavs << " = " << (*cav)->total_volume << " A^3" << endl;
                    string n_ext = "_"; string_fu::add2s(n_ext,n_cavs); string_fu::add2s(n_ext,".pdb");
                    string lname = string_fu::get_text_after(*fit,'/');
                    string_fu::replace_ext(lname,".pdb",n_ext);
                    inparser->write_pdb_cav(lname.c_str(),*cav);
                    cout << "    written cavity " << n_cavs << " as " << lname << endl;
                    ++n_cavs;
                }
                
                if (sourcefile != "X") {
                    delete lig_struct;
                    delete lig_parser;
                }
                if (pdb_source != "X") {
                    delete rstructure;
                    delete rparser;
                }
                instruct->clear();
                continue;
            }
            else if (search_substruct) {
                if (sourcefile == "X") {
                    cerr << c_message<cERROR>("you have to specify a reference file by --s=your_ref.mol2") << endl;
                    exit(1);
                }
                STRUCTURE *rstructure;
                PARSER *rparser;
                rstructure = new STRUCTURE();
                rparser = new PARSER(rstructure,verb,true,true,true,parse_comments);
                if (!(rparser->read_mol2(sourcefile.c_str()))) {
                    cerr << c_message<cERROR>("could not read '") << sourcefile << "' as mol2 file" << endl;
                    delete rstructure;
                    delete rparser;
                    exit(1);
                }
                if (rstructure->ligands.size() == 0) {
                    cerr << c_message<cERROR>("found no molecule in '") << sourcefile << endl;
                    delete rstructure;
                    delete rparser;
                    exit(1);
                }
                
                rstructure->ligands[0]->get_hybrid_only();
                for (ligands_vec lig=instruct->ligands.begin(); lig!=instruct->ligands.end(); ++lig) {
                    if ((*lig)->atoms.size() < rstructure->ligands[0]->atoms.size()) continue;
                    (*lig)->get_hybrid_only();
                    int n_sub = (*lig)->has_substructure(rstructure->ligands[0],ss_with_exact_hybrid);
                    if (n_sub > 0) cout << *fit << " :  " << (*lig)->name << "  " << n_sub << "\n";
                }
                cout.flush();
                
                instruct->clear();
                delete rstructure;
                delete rparser;
                continue;
            }
            else if ((!gethyd) || (!getwat) || rewrite || fix_pdb || convert || full_convert || 
                     !(parse_ligands) || optalign || optalign2 || pock_align) {
                if (pdbqt_mode && convert) {
                    cout << "warning: '-c' is not allowed for pdbqt => switching to '-f'" << endl;
                    convert = false;
                    full_convert = true;
                }
                if (convert) { //!protein als mol2 schreiben
                    if (instruct->pdb2mol2(atom_typing,def_file.c_str())) {
                        if (sourcefiles.size() > 1) {
                            tfile = *fit;
                            if (pdbqt_mode) replace_ext(tfile,".pdbqt",".pdb");
                            if (targetfile != "X") {
                                ostringstream os;
                                os << "_" << targetfile << ".mol2";
                                replace_ext(tfile,".pdb",os.str());
                            } else replace_ext(tfile,"pdb","mol2");
                        } else {
                            if (targetfile != "X") tfile = targetfile;
                            else {
                                tfile = *fit;
                                if (pdbqt_mode) replace_ext(tfile,".pdbqt",".pdb");
                                replace_ext(tfile,"pdb","mol2");
                            }
                        }
                        /*
                        if (add_hyd) {
                            for (chains_vec it=chains.begin(); it!=chains.end();++it)
                                (*it)->set_standard_protonation(verb);
                            }
                        }
                        */
                        string_fu::remove_char(tfile);
                        inparser->write_mol2_protein(instruct->protein,tfile.c_str());
                    }
                    instruct->clear();
                    continue;
                } else if (full_convert) { //!ganzes pdb als mol2 schreiben und atom_typing nutzen
                    instruct->full_pdb2mol2(atom_typing,true,def_file.c_str(),true,verb,true,
                                            false,"X",max_ring_members,no_planar_free_rot,0,
                                            0,prot_acids,prot_guanidin,prot_amidin,
                                            prot_amin,prot_phosphate,prot_sulfate);
                    if (sourcefiles.size() > 1) {
                        tfile = *fit;
                        if (pdbqt_mode) replace_ext(tfile,".pdbqt",".pdb");
                        if (targetfile != "X") {
                            ostringstream os;
                            os << "_" << targetfile << ".mol2";
                            replace_ext(tfile,".pdb",os.str());
                        } else replace_ext(tfile,"pdb","mol2");
                    } else {
                        if (targetfile != "X") tfile = targetfile;
                        else {
                            tfile = *fit;
                            if (pdbqt_mode) replace_ext(tfile,".pdbqt",".pdb");
                            replace_ext(tfile,"pdb","mol2");
                        }
                    }
                    if (add_hyd) {
                        for (ligands_vec lit=instruct->mol2_protein.begin(); lit!=instruct->mol2_protein.end(); ++lit) {
                            (*lit)->set_standard_protonation(verb);
                        }
                    }
                    string_fu::remove_char(tfile);
                    inparser->write_mol2(instruct->mol2_protein,tfile.c_str());
                    instruct->clear();
                    continue;
                } if (optalign || optalign2) {
                    if (sourcefile == "X") {
                        cerr << c_message<cERROR>("you have to specify a reference file by --s=your_ref.pdb") << endl;
                        exit(1);
                    }
                    STRUCTURE *rstructure;
                    PARSER *rparser;
                    rstructure = new STRUCTURE();
                    rparser = new PARSER(rstructure,verb,true,true,false,parse_comments);
                    if (!(rparser->read_pdb(sourcefile.c_str()))) { //!pdb einlesen
                        cerr << c_message<cERROR>("could not read '") << sourcefile << "' as pdb file" << endl;
                        delete rstructure;
                        delete rparser;
                        exit(1);
                    }
                    if (rstructure->protein == 0) {
                        cerr << c_message<cERROR>("found no protein residues in '") << sourcefile << endl;
                        delete rstructure;
                        delete rparser;
                        exit(1);
                    }
                    if (instruct->protein == 0) {
                        cerr << c_message<cERROR>("found no protein residues in '") << *fit << endl;
                        instruct->clear();
                        continue;
                    }
                    
                    bool spatial = true;
                    string pro_mode = "pdb";
                    if (pdb_source != "X") {
                        spatial = false;
                        if (pdb_source[pdb_source.size()-1] == '2') pro_mode = "mol";
                    }
                    cout << "aligning " << *fit << " :" << endl;
                    if (optalign2) {
                        if (!(instruct->align_by_protein2(rstructure->protein,spatial,max_displacement))) {
                            if (verb > 0) cerr << c_message<cWARNING>("found no reasonable alignment for ") << *fit << endl;
                        }
                    } else {
                        instruct->align_by_protein(rstructure->protein,spatial,debug_mode,only_C_alpha);
                    }
                    /*
                    cout << "1. translation: " << instruct->protein->optalign_trans[1] << "\n";
                    cout << "rotation matrix: " << instruct->protein->optalign_rotm << "\n";
                    cout << "2. translation: " << instruct->protein->optalign_trans[0] << "\n" << "\n";
                    */
                    if (sourcefiles.size() > 1) {
                        tfile = *fit;
                        if (pdbqt_mode) replace_ext(tfile,".pdbqt",".pdb");
                        if (targetfile != "X") {
                            ostringstream os;
                            os << "_" << targetfile << ".pdb";
                            replace_ext(tfile,".pdb",os.str());
                        }
                        if (!spatial) {
                            tfile = pdb_source;
                            if (targetfile != "X") {
                                ostringstream os;
                                os << "_" << targetfile << ".mol2";
                                replace_ext(tfile,".mol2",os.str());
                            }
                        }
                    } else {
                        if (targetfile != "X") tfile = targetfile;
                        else {
                            tfile = *fit;
                            if (pdbqt_mode) replace_ext(tfile,".pdbqt",".pdb");
                            if (!spatial) tfile = pdb_source;
                        }
                    }
                    string_fu::remove_char(tfile);
                    if (spatial) inparser->write_pdb(tfile.c_str());
                    else {
                        STRUCTURE *rstructure2;
                        PARSER *rparser2;
                        rstructure2 = new STRUCTURE();
                        rparser2 = new PARSER(rstructure2,0,true,true,true,parse_comments);
                        
                        if (pro_mode == "mol") {
                            if (!(rparser2->read_mol2(pdb_source.c_str()))) {
                                cerr << c_message<cERROR>("could not read '") << pdb_source.c_str() << "' as mol2 file" << endl;
                                delete rstructure2;
                                delete rparser2;
                                exit(1);
                            }
                            for (ligands_vec lig=rstructure2->ligands.begin(); lig!=rstructure2->ligands.end(); ++lig) {
                                for (atoms_vec at=(*lig)->atoms.begin(); at!=(*lig)->atoms.end();++at) {
                                    (*at)->coord -= instruct->protein->optalign_trans[1];
                                    (*at)->coord *= instruct->protein->optalign_rotm;
                                    (*at)->coord += instruct->protein->optalign_trans[0];
                                }
                            }
                            rparser2->write_mol2(tfile.c_str());
                        } else {
                            if (!(rparser2->read_pdb(pdb_source.c_str()))) {
                                cerr << c_message<cERROR>("could not read '") << pdb_source.c_str() << "' as pdb file" << endl;
                                delete rstructure2;
                                delete rparser2;
                                exit(1);
                            }
                            if (rstructure2->protein) for (chains_vec ct=rstructure2->protein->chains.begin(); ct!=rstructure2->protein->chains.end(); ++ct) {
                                for (atoms_vec at=(*ct)->atoms.begin(); at!=(*ct)->atoms.end();++at) {
                                    (*at)->coord -= instruct->protein->optalign_trans[1];
                                    (*at)->coord *= instruct->protein->optalign_rotm;
                                    (*at)->coord += instruct->protein->optalign_trans[0];
                                }
                            }
                            for (ligands_vec lig=rstructure2->ligands.begin(); lig!=rstructure2->ligands.end(); ++lig) {
                                for (atoms_vec at=(*lig)->atoms.begin(); at!=(*lig)->atoms.end();++at) {
                                    (*at)->coord -= instruct->protein->optalign_trans[1];
                                    (*at)->coord *= instruct->protein->optalign_rotm;
                                    (*at)->coord += instruct->protein->optalign_trans[0];
                                }
                            }
                            for (waters_vec wat=rstructure2->waters.begin(); wat!=rstructure2->waters.end(); ++wat) {
                                (*wat)->atom->coord -= instruct->protein->optalign_trans[1];
                                (*wat)->atom->coord *= instruct->protein->optalign_rotm;
                                (*wat)->atom->coord += instruct->protein->optalign_trans[0];
                            }
                            for (metals_vec met=rstructure2->metals.begin(); met!=rstructure2->metals.end(); ++met) {
                                (*met)->atom->coord -= instruct->protein->optalign_trans[1];
                                (*met)->atom->coord *= instruct->protein->optalign_rotm;
                                (*met)->atom->coord += instruct->protein->optalign_trans[0];
                            }
                            rparser2->write_pdb(tfile.c_str());
                        }
                        
                        delete rstructure2;
                        delete rparser2;
                    }
                    instruct->clear();
                    delete rstructure;
                    delete rparser;
                    continue;
                } else if (pock_align) {
                    /*
                    if (sourcefile == "X") {
                        cerr << c_message<cERROR>("you have to specify a reference file by --s=your_ref.mol2") << endl;
                        exit(1);
                    }
                    STRUCTURE *rstructure;
                    PARSER *rparser;
                    rstructure = new STRUCTURE();
                    rparser = new PARSER(rstructure,verb,true,true,false,parse_comments);
                    if (!(rparser->read_pdb(sourcefile.c_str()))) { //!pdb einlesen
                        cerr << c_message<cERROR>("could not read '") << sourcefile << "' as pdb file" << endl;
                        delete rstructure;
                        delete rparser;
                        exit(1);
                    }
                    if (rstructure->protein == 0) {
                        cerr << c_message<cERROR>("found no protein residues in '") << sourcefile << endl;
                        delete rstructure;
                        delete rparser;
                        exit(1);
                    }
                    if (instruct->protein == 0) {
                        cerr << c_message<cERROR>("found no protein residues in '") << *fit << endl;
                        instruct->clear();
                        continue;
                    }
                    instruct->get_cavity_auto(cut_radius,min_buried,debug_mode,1);
                    rstructure->get_cavity_auto(cut_radius,min_buried,debug_mode,1);
                    instruct->align_by_cavity(rstructure);
                    if (sourcefiles.size() > 1) {
                        tfile = *fit;
                        if (targetfile != "X") {
                            ostringstream os;
                            os << "_" << targetfile << ".pdb";
                            replace_ext(tfile,".pdb",os.str());
                        }
                    } else {
                        if (targetfile != "X") tfile = targetfile;
                        else tfile = *fit;
                    }
                    string_fu::remove_char(tfile);
                    inparser->write_pdb(tfile.c_str());
                    instruct->clear();
                    delete rstructure;
                    delete rparser;
                    */
                    continue;
                } else if (fix_pdb) {
                    if (sourcefiles.size() > 1) {
                        tfile = *fit;
                        if (pdbqt_mode) replace_ext(tfile,".pdbqt",".pdb");
                        if (targetfile != "X") {
                            ostringstream os;
                            os << "_" << targetfile << ".pdb";
                            replace_ext(tfile,".pdb",os.str());
                        }
                    } else {
                        if (targetfile != "X") tfile = targetfile;
                        else {
                            tfile = *fit;
                            if (pdbqt_mode) replace_ext(tfile,".pdbqt",".pdb");
                        }
                    }
                    //! Atomtyping, damit das richtige Element im pdb File steht
                    for (chains_vec cvt=instruct->protein->chains.begin(); cvt!=instruct->protein->chains.end(); ++cvt) {
                        for (atoms_vec at=(*cvt)->atoms.begin(); at!=(*cvt)->atoms.end(); ++at) {
                            (*at)->ext = new ATOM_EXT();
                            (*at)->charge = "";
                        }
                        (*cvt)->get_elements(false);
                    }
                    for (ligands_vec cvt=instruct->ligands.begin(); cvt!=instruct->ligands.end(); ++cvt) {
                        for (atoms_vec at=(*cvt)->atoms.begin(); at!=(*cvt)->atoms.end(); ++at) {
                            (*at)->ext = new ATOM_EXT();
                            (*at)->charge = "";
                        }
                        (*cvt)->get_elements(false);
                    }
                    for (metals_vec cvt=instruct->metals.begin(); cvt!=instruct->metals.end(); ++cvt) {
                        (*cvt)->atom->ext = new ATOM_EXT();
                        (*cvt)->atom->get_element(false);
                        (*cvt)->atom->charge = "";
                    }
                    for (waters_vec cvt=instruct->waters.begin(); cvt!=instruct->waters.end(); ++cvt) {
                        (*cvt)->atom->ext = new ATOM_EXT();
                        (*cvt)->atom->get_element(false);
                        (*cvt)->atom->charge = "";
                    }
                    
                    //! Eventuell fehlende HET-Eintraege erzeugen
                    
                    string_fu::remove_char(tfile);
                    inparser->write_pdb(tfile.c_str());
                    instruct->clear();
                    continue;
                } else { //!pdb neu schreiben (rewrite)
                    if (sourcefiles.size() > 1) {
                        tfile = *fit;
                        if (pdbqt_mode) replace_ext(tfile,".pdbqt",".pdb");
                        if (targetfile != "X") {
                            ostringstream os;
                            os << "_" << targetfile << ".pdb";
                            replace_ext(tfile,".pdb",os.str());
                        }
                    } else {
                        if (targetfile != "X") tfile = targetfile;
                        else {
                            tfile = *fit;
                            if (pdbqt_mode) replace_ext(tfile,".pdbqt",".pdb");
                        }
                    }
                    string_fu::remove_char(tfile);
                    inparser->write_pdb(tfile.c_str());
                    instruct->clear();
                    continue;
                }
            }
            
        } else if (mode == "mol" || mode == "sdf") {
            //! Neu: partiell einlesen *a*
            unsigned int total_mol_read = 0;
            int block_size = max_mol2_block;
            if (remove_duplicates) block_size = 999999999;
            int curr_bs_needed = 0; // fuer Blocksplit
            int fnmw = 1;
            int tnmw = 0;
            if (blocksplit) {
                if (block_size > blocksplit_size) block_size = blocksplit_size;
                curr_bs_needed = blocksplit_size;
            }
            bool go_on = true;
            must_append = false;
            
            if (extract_all_ligs_cav && mode == "mol") { //! am 08.01.09 fuer Andreas S. reingenommen
                if (sourcefile == "X") {
                    cerr << c_message<cERROR>("you have to specify a reference file by --s=your_lig.mol2") << endl;
                    exit(1);
                }
                
                if (!(inparser->read_mol2((*fit).c_str()))) {
                    cerr << c_message<cERROR>("could not read '") << (*fit).c_str() << "' as mol2 file" << endl;
                    continue;
                }
                
                STRUCTURE *mstruct;
                PARSER *mparser;
                mstruct = new STRUCTURE();
                bool phyd = !no_hyd_for_poc;
                mparser = new PARSER(mstruct,verb,phyd,true,false);
                if (!(mparser->read_mol2(sourcefile.c_str()))) {
                    cerr << c_message<cERROR>("could not read '") << sourcefile.c_str() << "' as mol2 file" << endl;
                    delete mstruct;
                    delete mparser;
                    exit(1);
                }
                
                if (mstruct->ligands.size() == 0) {
                    cerr << c_message<cERROR>("found no ligands in '") << sourcefile.c_str() << "'" << endl;
                    delete mstruct;
                    delete mparser;
                    exit(1);
                }
                
                float srad = cut_radius * cut_radius;
                vector<stl_ptr<LIGAND> > cav_ligs;
                for (ligands_vec lig=instruct->ligands.begin(); lig!=instruct->ligands.end(); ++lig) {
                    stl_ptr<LIGAND> nl = new LIGAND();
                    nl->name = (*lig)->name;
                    nl->type = (*lig)->type;
                    nl->charge_type = (*lig)->charge_type;
                    for (bonds_vec bv=(*lig)->bonds.begin(); bv!=(*lig)->bonds.end(); ++bv) {
                        nl->bonds.push_back(*bv);
                    }
                    bool added = false;
                    
                    tr1::unordered_set<int> poc_atoms;
                    tr1::unordered_set<int> poc_res;
                    
                    for (atoms_vec at=(*lig)->atoms.begin(); at!=(*lig)->atoms.end(); ++at) {
                        bool breaker = false;
                        for (ligands_vec mlig=mstruct->ligands.begin(); mlig!=mstruct->ligands.end(); ++mlig) {
                        for (atoms_vec bt=(*mlig)->atoms.begin(); bt!=(*mlig)->atoms.end(); ++bt) {
                            if (get_square_distance((*at)->coord,(*bt)->coord) <= srad) {
                                nl->atoms.push_back(*at);
                                if (cut_full_res) {
                                    poc_atoms.insert((*at)->intern_id);
                                    poc_res.insert((*at)->res_number);
                                }
                                if (!added) {
                                    added = true;
                                    cav_ligs.push_back(nl);
                                }
                                breaker =true;
                                break;
                            }
                        }
                        if (breaker) break;
                        }
                    }
                    
                    if (cut_full_res) {
                        for (atoms_vec at=(*lig)->atoms.begin(); at!=(*lig)->atoms.end(); ++at) {
                            if (poc_res.find((*at)->res_number) != poc_res.end()) {
                                if (poc_atoms.find((*at)->intern_id) == poc_atoms.end()) nl->atoms.push_back(*at);
                            }
                        }
                    }
                }
                
                string sf_name = sourcefile; remove_ext(sf_name);
                string filename = "_" + sf_name + "_poc" + ".mol2";
                if (targetfile != "X") tfile = targetfile;
                else {
                    tfile = *fit;
                    replace_ext(tfile,".mol2",filename);
                }
                string_fu::remove_char(tfile);
                inparser->write_mol2(cav_ligs,tfile.c_str());
                delete mstruct;
                delete mparser;
                
                STRUCTURE *tstruct;
                PARSER *tparser;
                tstruct = new STRUCTURE();
                tparser = new PARSER(tstruct,0,true,true,true,parse_comments);
                tparser->read_mol2(tfile.c_str());
                tparser->write_mol2(tfile.c_str());
                delete tstruct;
                delete tparser;
                
                instruct->clear();
                continue;
            } else if (extract_all_ligs_cav && mode == "sdf") {
                cerr << c_message<cERROR>("option not available for sdf files") << endl;
                exit(1);
            }
            
            while (go_on) {
                int n_mol_read = 0;
            
            if (mode == "mol" && (retype || add_hyd || rmsd || rmsd2 || rmsd3 || optalign || optalign2 || sdiv || get_n_rings || get_n_rings2 ||
                search_substruct || get_mol_weight || calc_optalign_matrix || get_rot_bonds || type_count ||
                remove_duplicates || clashes || intern_split || keep_biggest_fragment || cluster || get_vol)) {
                n_mol_read = inparser->read_next_mol2((*fit).c_str(),block_size,true);
                if (n_mol_read == -1) break;
                total_mol_read += n_mol_read;
            } else if (mode == "sdf" && (convert || rmsd || rmsd2 || optalign || optalign2 || get_rot_bonds ||
                                         type_count || get_vol || cluster || search_substruct || ss_with_exact_type ||
                                         calc_optalign_matrix || get_n_rings || get_n_rings2 || get_mol_weight ||
                                         split_multimol || extract_mol2 || extract_mol2_by_name)) {
                n_mol_read = inparser->read_next_sdf((*fit).c_str(),block_size);
                if (n_mol_read < 1) break;
                total_mol_read += n_mol_read;
            } else if (mode == "mol") {
                if (blocksplit && curr_bs_needed < block_size) n_mol_read = inparser->read_next_mol2((*fit).c_str(),curr_bs_needed);
                else n_mol_read = inparser->read_next_mol2((*fit).c_str(),block_size);
                tnmw += n_mol_read;
                if (n_mol_read == -1) break;
                total_mol_read += n_mol_read;
            } else if (mode == "sdf" && blocksplit) {
                if (curr_bs_needed < block_size) n_mol_read = inparser->read_next_sdf((*fit).c_str(),curr_bs_needed);
                else n_mol_read = inparser->read_next_sdf((*fit).c_str(),block_size);
                tnmw += n_mol_read;
                if (n_mol_read == -1) break;
                total_mol_read += n_mol_read;
            } else if (mode == "sdf") {
                cerr << c_message<cERROR>("option is currently not available for sdf files") << endl;
                exit(1);
            }
            if (n_mol_read < block_size) {
                if (blocksplit) {
                    if (n_mol_read < curr_bs_needed) go_on = false;
                } else go_on = false; //! ***
            }

            if (get_rot_bonds) {
            //    cout << "Calculating free rotatable bonds for '" << (*fit).c_str() << "':" << "\n";
                int rotis;
                STRUCTURE *rstructure = 0;
                PARSER *rparser = 0;
                LIGAND* at_ref_mol = 0;
                if (pdb_source != "X") {
                    rstructure = new STRUCTURE();
                    rparser = new PARSER(rstructure,0,true,true,true,parse_comments);
                    if (!(rparser->read_mol2(pdb_source.c_str(),true))) {
                        cerr << c_message<cERROR>("could not read '") << pdb_source << "' as mol2 file" << endl;
                        delete rstructure;
                        delete rparser;
                        exit(1);
                    }
                    if (rstructure->ligands.size() == 0) {
                        cerr << c_message<cERROR>("found no molecule in '") << pdb_source << endl;
                        delete rstructure;
                        delete rparser;
                        exit(1);
                    }
                    at_ref_mol = &(*(rstructure->ligands[0]));
                }
                
                STRUCTURE *mstruct = 0;
                PARSER *mparser = 0;
                string pro_mode = "pdb";
                vector<stl_ptr<ATOM> > pro_atoms;
                if (sourcefile != "X") { //protein => auf kovalente Bindungen checken
                    if (sourcefile[sourcefile.size()-1] == '2') pro_mode = "mol";
                    mstruct = new STRUCTURE();
                    mparser = new PARSER(mstruct,0,true,true,false);
                    if (pro_mode == "mol") {
                        if (!(mparser->read_mol2(sourcefile.c_str()))) {
                            cerr << c_message<cERROR>("could not read '") << sourcefile.c_str() 
                                 << "' as mol2 file" << endl;
                            delete mstruct;
                            delete mparser;
                            exit(1);
                        }
                        for (ligands_vec lig=mstruct->ligands.begin(); lig!=mstruct->ligands.end(); ++lig) {
                            if (debug_mode) {
                                (*lig)->get_atom_typing(atom_typing,true,def_file.c_str(),true,0,false,false,"X",
                                                            max_ring_members,no_planar_free_rot,at_ref_mol,0,prot_acids,prot_guanidin,prot_amidin,prot_amin,prot_phosphate,
                                                   prot_sulfate,kekulize_aromatics,kekulize_charged,allow_charged_aromatics);
                            }
                            for (atoms_vec at=(*lig)->atoms.begin(); at!=(*lig)->atoms.end(); ++at) {
                                pro_atoms.push_back(*at);
                            }
                        }
                    } else {
                        if (!(mparser->read_pdb(sourcefile.c_str()))) { //!pdb einlesen
                            cerr << c_message<cERROR>("could not read '") << sourcefile << "' as pdb file" << endl;
                            delete mstruct;
                            delete mparser;
                            exit(1);
                        }
                        if (mstruct->protein == 0) {
                            cerr << c_message<cERROR>("found no residues in '") << sourcefile << "'" << endl;
                            delete mstruct;
                            delete mparser;
                            exit(1);
                        }
                        if (debug_mode) {
                            //! Weil die Bindungen gebraucht werden, um die komplette Torsion der kovalenten Bindung auszugeben:
                            for (chains_vec ct=mstruct->protein->chains.begin(); ct!=mstruct->protein->chains.end(); ++ct) {
                                (*ct)->get_atom_typing(atom_typing,true,def_file.c_str(),true,0,false,false,"X",max_ring_members,
                                                       no_planar_free_rot,at_ref_mol,0,prot_acids,prot_guanidin,prot_amidin,prot_amin,prot_phosphate,
                                                   prot_sulfate,kekulize_aromatics,kekulize_charged,allow_charged_aromatics);
                            }
                        }
                        for (chains_vec ct=mstruct->protein->chains.begin(); ct!=mstruct->protein->chains.end(); ++ct) {
                            for (atoms_vec at=(*ct)->atoms.begin(); at!=(*ct)->atoms.end(); ++at) {
                                pro_atoms.push_back(*at);
                            }
                        }
                    }
                }
                
                vector<stl_ptr<LIGAND> > ligs2write;
                int l_idx = 0;
                ofstream rout;
                if (debug_mode) {
                    if (targetfile != "X") rout.open(targetfile.c_str());
                }
                for (ligands_vec lig=instruct->ligands.begin(); lig!=instruct->ligands.end(); ++lig) {
                    (*lig)->get_atom_typing(atom_typing,true,def_file.c_str(),true,verb,false,false,"X",max_ring_members,
                                            no_planar_free_rot,at_ref_mol,0,prot_acids,prot_guanidin,prot_amidin,prot_amin,prot_phosphate,
                                                   prot_sulfate,kekulize_aromatics,kekulize_charged,allow_charged_aromatics);
                    rotis = (*lig)->get_freerot_bonds();
                    
                    int dat = 0;
                    stl_ptr<ATOM> dat1;
                    stl_ptr<ATOM> dat2;
                    stl_ptr<ATOM> pdat1;
                    stl_ptr<ATOM> pdat2;
                    if (sourcefile != "X") {
                        for (atoms_vec at=(*lig)->atoms.begin(); at!=(*lig)->atoms.end(); ++at) {
                            for (atoms_vec bt=pro_atoms.begin(); bt!=pro_atoms.end(); ++bt) {
                                if (get_square_distance((*at)->coord,(*bt)->coord) < 0.07) {
                                    ++dat;
                                    if (dat == 1) {
                                        dat1 = *at;
                                        pdat1 = *bt;

//                                        cerr << "gemeinsames Atom 1: lig = " << dat1->id << "   pro = " << pdat1->id << endl;

                                    } else {
                                        dat2 = *at;
                                        pdat2 = *bt;

//                                        cerr << "gemeinsames Atom 2: lig = " << dat2->id << "   pro = " << pdat2->id << endl;
//                                        cerr << dat2->name << ", " << dat2->res_name << ", " << dat2->res_number << endl;
//                                        cerr << pdat2->name << ", " << pdat2->res_name << ", " << pdat2->res_number << endl;
//                                        cerr << "dist = " << get_distance(dat2->coord,pdat2->coord) << endl;
                                    }
                                    break;
                                }
                            }
                            if (dat == 2) break;
                        }
                    }
                    
                    if (dat > 0) ++rotis;
                    
                    if (dat == 1) {
                        for (atoms_vec xt=dat1->bonded_atoms.begin(); xt!=dat1->bonded_atoms.end(); ++xt) {
                            if ((*xt)->element == "H") continue;
                            dat2 = *xt;
                            ++dat;

//                            cerr << "dat2 = " << dat2->id << endl;

                        }
                        if (dat < 2) {
                            cerr << c_message<cWARNING>("detected shared atoms ") << dat1->name << "(protein) and "
                                 << pdat1->name << "(ligand), but no heavy atom is connected to " << dat1->name << endl;
                            dat = 0;
                        }
                        bool has_t1 = false;
                        for (atoms_vec at=pdat1->bonded_atoms.begin(); at!=pdat1->bonded_atoms.end(); ++at) {
                            if ((*at)->element == "H") continue;
                            has_t1 = true;
                            break;
                        }
                        if (!has_t1) {
                            cerr << c_message<cWARNING>("detected shared atoms ") << dat1->name << "(protein) and "
                                 << pdat1->name << "(ligand), but no heavy atom is connected to " << pdat1->name << endl;
                            dat = 0;
                        }
                    }
                    
                    if (targetfile == "X") {
                        cout << "# " << l_idx << " ";
                        cout.width(30); cout << left << (*lig)->name;
                        cout << ":  "; cout.width(6); cout << left << rotis;
                        cout << "('" << (*fit).c_str() << "')" << "\n";
                    } else if (!debug_mode) {
                        if (rotis <= int(unit_cell_multi)) ligs2write.push_back(*lig);
                    }
                    if (debug_mode) { //die entsprechenden bond id's ausgeben
                        if (targetfile != "X") {
                            cout << "# " << l_idx << " ";
                            cout.width(30); cout << left << (*lig)->name;
                            cout << ":  "; cout.width(6); cout << left << rotis;
                            cout << "('" << (*fit).c_str() << "')" << "\n";
                            if (sourcefile != "X" && dat > 0) {
                                rout << "   " << dat1->id << " " << dat2->id << " coval" << "\n";
                            }
                            for (bonds_vec tit=(*lig)->bonds.begin(); tit!=(*lig)->bonds.end();++tit) {
                                if ((*tit)->free_rot) {
                                    if (targetfile != "X") rout << (*tit)->from->id 
                                                    << " " << (*tit)->to->id << "\n";
                                }
                            }
                        } else {
                            
                            //!=========================================================================
                            //! Nicht nur IDs fuer die Bindungen, sondern komplett 4 IDs fuer Torsion ausgeben
                            stl_ptr<ATOM> t1 = 0;
                            stl_ptr<ATOM> t2 = 0;
                            stl_ptr<ATOM> t3 = 0;
                            stl_ptr<ATOM> t4 = 0;
                            
                            if (sourcefile != "X" && dat > 0) {
                                if (pdat2.zero()) { //! Nur ein gemeinsames Atom => dat1 ist auf Proteinseite der Torsion
                                    t2 = dat1;
                                    t3 = dat2;
                                    for (atoms_vec at=pdat1->bonded_atoms.begin(); at!=pdat1->bonded_atoms.end(); ++at) {
                                        if ((*at)->element == "H") continue;
                                        t1 = *at;
                                        //cerr << "got t1" << endl;
                                        break;
                                    }                                
                                    for (atoms_vec at=t3->bonded_atoms.begin(); at!=t3->bonded_atoms.end(); ++at) {
                                        if ((*at)->element == "H") continue;
                                        if (*at == t2) continue;
                                        t4 = *at;
                                        //cerr << "got t4" << endl;
                                        break;
                                    }
                                } else {    
                                    if (dat1->ext->n_heavy_bonded < dat2->ext->n_heavy_bonded) { //! dat1 ist auf Proteinseite
                                        t2 = dat1;
                                        t3 = dat2;
                                        for (atoms_vec at=pdat1->bonded_atoms.begin(); at!=pdat1->bonded_atoms.end(); ++at) {
                                            if ((*at)->element == "H") continue;
                                            t1 = *at;
                                            break;
                                        }
                                        for (atoms_vec at=t3->bonded_atoms.begin(); at!=t3->bonded_atoms.end(); ++at) {
                                            if ((*at)->element == "H") continue;
                                            if (*at == t2) continue;
                                            t4 = *at;
                                            break;
                                        }
                                    } else { //! dat2 ist auf Proteinseite:
                                        t2 = dat2;
                                        t3 = dat1;
                                        for (atoms_vec at=pdat2->bonded_atoms.begin(); at!=pdat2->bonded_atoms.end(); ++at) {
                                            if ((*at)->element == "H") continue;
                                            t1 = *at;
                                            break;
                                        }
                                        for (atoms_vec at=t3->bonded_atoms.begin(); at!=t3->bonded_atoms.end(); ++at) {
                                            if ((*at)->element == "H") continue;
                                            if (*at == t2) continue;
                                            t4 = *at;
                                            break;
                                        }
                                    }
                                }

                                cout << "   " << t1->id << " " << t2->id << " " << t3->id << " " << t4->id << " coval (first id is protein)\n";
                                //! Die weiteren Bindungen so ausgeben, dass t1 immer naeher an coval_t3 ist als t4:
                                if ((*lig)->SP_map == 0) (*lig)->calc_SP_map();
                                int ref_id = t3->intern_id;
                                for (bonds_vec tit=(*lig)->bonds.begin(); tit!=(*lig)->bonds.end();++tit) {
                                    if ((*tit)->free_rot) {
                                        t2 = (*tit)->from;
                                        t3 = (*tit)->to;
                                        for (atoms_vec at=t2->bonded_atoms.begin(); at!=t2->bonded_atoms.end(); ++at) {
                                            if ((*at)->element == "H") continue;
                                            if (*at == t3) continue;
                                            t1 = *at;
                                            break;
                                        }
                                        for (atoms_vec at=t3->bonded_atoms.begin(); at!=t3->bonded_atoms.end(); ++at) {
                                            if ((*at)->element == "H") continue;
                                            if (*at == t2) continue;
                                            t4 = *at;
                                            break;
                                        }
                                        
                                        int sp_1,sp_2;
                                        if (ref_id > t1->intern_id) sp_1 = (*lig)->SP_map[ref_id][t1->intern_id];
                                        else sp_1 = (*lig)->SP_map[t1->intern_id][ref_id];
                                        if (ref_id > t4->intern_id) sp_2 = (*lig)->SP_map[ref_id][t4->intern_id];
                                        else sp_2 = (*lig)->SP_map[t4->intern_id][ref_id];
                                        
                                        if (sp_1 < sp_2) {
                                            cout << "   " << t1->id << " " << t2->id << " " << t3->id << " " << t4->id << "\n";
                                        } else {
                                            cout << "   " << t4->id << " " << t3->id << " " << t2->id << " " << t1->id << "\n";
                                        }
                                    }
                                }
                            } else {
                                for (bonds_vec tit=(*lig)->bonds.begin(); tit!=(*lig)->bonds.end();++tit) {
                                    if ((*tit)->free_rot) {
                                        t2 = (*tit)->from;
                                        t3 = (*tit)->to;
                                        for (atoms_vec at=t2->bonded_atoms.begin(); at!=t2->bonded_atoms.end(); ++at) {
                                            if ((*at)->element == "H") continue;
                                            if (*at == t3) continue;
                                            t1 = *at;
                                            break;
                                        }
                                        for (atoms_vec at=t3->bonded_atoms.begin(); at!=t3->bonded_atoms.end(); ++at) {
                                            if ((*at)->element == "H") continue;
                                            if (*at == t2) continue;
                                            t4 = *at;
                                            break;
                                        }
                                        cout << "   " << t1->id << " " << t2->id << " " << t3->id << " " << t4->id << "\n";
                                    }
                                }
                            }
                            //!=========================================================================
                        }
                    }
                    ++l_idx;
                }
                if (debug_mode) {
                    if (targetfile != "X") rout.close();
                }
                if ((!debug_mode) && targetfile != "X") {
                    inparser->write_next_mol2(ligs2write,targetfile.c_str());
                    if (go_on) must_append = true;
                }
                if (pdb_source != "X") {
                    delete rstructure;
                    delete rparser;
                }
                if (sourcefile != "X") {
                    delete mstruct;
                    delete mparser;
                }
                instruct->clear();
                continue;
            }
            if (type_count) {
                for (ligands_vec lit=instruct->ligands.begin(); lit!=instruct->ligands.end(); ++lit) {
                    (*lit)->get_atom_typing(atom_typing,true,def_file.c_str(),false,verb,true,true,"X",max_ring_members,
                                            true,0,0,prot_acids,prot_guanidin,prot_amidin,prot_amin,prot_phosphate,
                                            prot_sulfate,kekulize_aromatics,kekulize_charged,allow_charged_aromatics);
                    map<string,int> cmap;
                    for (atoms_vec at=(*lit)->atoms.begin(); at!=(*lit)->atoms.end(); ++at) {
                        if ((*at)->sybyl_type == "X") continue;
                        if (cmap.find((*at)->sybyl_type) == cmap.end()) cmap[(*at)->sybyl_type] = 1;
                        else cmap[(*at)->sybyl_type] += 1;
                    }
                    cout << (*lit)->name << ":" << "\n";
                    for (map<string,int>::iterator it=cmap.begin(); it!=cmap.end(); ++it) {
                        cout << "   -> ";
                        cout.width(8); cout << left << it->first << " : " << it->second << "\n";
                    }
                }
                instruct->clear();
                continue;
            }
            if (rmsd) {
                if (sourcefile == "X") {
                    cerr << c_message<cERROR>("you have to specify a reference file by --s=your_ref.mol2") << endl;
                    exit(1);
                }
                STRUCTURE *rstructure;
                PARSER *rparser;
                rstructure = new STRUCTURE();
                rparser = new PARSER(rstructure,verb,true,true,true,parse_comments);
                if (!(rparser->read_mol2(sourcefile.c_str()))) {
                    cerr << c_message<cERROR>("could not read '") << sourcefile << "' as mol2 file" << endl;
                    delete rstructure;
                    delete rparser;
                    exit(1);
                }
                if (rstructure->ligands.size() == 0) {
                    cerr << c_message<cERROR>("found no molecule in '") << sourcefile << endl;
                    delete rstructure;
                    delete rparser;
                    exit(1);
                }


                STRUCTURE *refstructure = 0;
                PARSER *refparser = 0;
                LIGAND* at_ref_mol = 0;
                if (pdb_source != "X") {
                    refstructure = new STRUCTURE();
                    refparser = new PARSER(refstructure,0,true,true,true,parse_comments);
                    if (!(refparser->read_mol2(pdb_source.c_str(),true))) {
                        cerr << c_message<cERROR>("could not read '") << pdb_source << "' as mol2 file" << endl;
                        delete refstructure;
                        delete refparser;
                        exit(1);
                    }
                    if (refstructure->ligands.size() == 0) {
                        cerr << c_message<cERROR>("found no molecule in '") << pdb_source << endl;
                        delete refstructure;
                        delete refparser;
                        exit(1);
                    }
                    at_ref_mol = &(*(refstructure->ligands[0]));
                    at_ref_mol->get_elements(false);
                }

                rstructure->ligands[0]->get_atom_typing(0,true,"X",true,0,false,false,"X",max_ring_members,false,0,at_ref_mol,prot_acids,prot_guanidin,prot_amidin,prot_amin,prot_phosphate,
                                                   prot_sulfate,kekulize_aromatics,kekulize_charged,allow_charged_aromatics);

                if (debug_mode) {
                    cout << "please note, that different optimal atom-atom-matches are possible for" << "\n";
                    cout << "RMSD and RMSD_opt => in debug mode you will get two match lists" << "\n";
                }
                for (ligands_vec lig=instruct->ligands.begin(); lig!=instruct->ligands.end(); ++lig) {
                    (*lig)->get_atom_typing(0,true,"X",true,0,false,false,"X",max_ring_members,false,0,at_ref_mol,prot_acids,prot_guanidin,prot_amidin,prot_amin,prot_phosphate,
                                                   prot_sulfate,kekulize_aromatics,kekulize_charged,allow_charged_aromatics);
                    cout << "@ "; cout.width(24); cout << left << (*lig)->name;
                    cout << ":  ";
                    float trms = (*lig)->get_rmsd(rstructure->ligands[0],false,debug_mode);
                    cout << "RMSD = "; cout.width(12); cout << left << trms << "   ";
                    cout << "RMSD_with_H = "; cout.width(12); cout << left << 
                    (*lig)->get_rmsd(rstructure->ligands[0],true,false,false);
                    if (debug_mode) {cout << "\n@ "; cout << left << (*lig)->name; cout << ":  ";}
                    cout << "RMSD_opt = " << (*lig)->get_rmsd(rstructure->ligands[0],false,debug_mode,true) << "\n";
                }
                instruct->clear();
                delete rstructure;
                delete rparser;
                continue;
            }
            if (rmsd2) {
                if (sourcefile == "X") {
                    cerr << c_message<cERROR>("you have to specify a reference file by --s=your_ref.mol2") << endl;
                    exit(1);
                }
                STRUCTURE *rstructure;
                PARSER *rparser;
                rstructure = new STRUCTURE();
                rparser = new PARSER(rstructure,verb,true,true,true,parse_comments);
                if (!(rparser->read_mol2(sourcefile.c_str()))) {
                    cerr << c_message<cERROR>("could not read '") << sourcefile << "' as mol2 file" << endl;
                    delete rstructure;
                    delete rparser;
                    exit(1);
                }
                if (rstructure->ligands.size() == 0) {
                    cerr << c_message<cERROR>("found no molecule in '") << sourcefile << endl;
                    delete rstructure;
                    delete rparser;
                    exit(1);
                }
                rstructure->ligands[0]->get_hybrid_only();
                if (debug_mode) {
                    cout << "please note, that different optimal atom-atom-matches are possible for" << "\n";
                    cout << "RMSD and RMSD_opt => in debug mode you will get two match lists" << "\n";
                }
                for (ligands_vec lig=instruct->ligands.begin(); lig!=instruct->ligands.end(); ++lig) {
                    cout << "@ "; cout.width(24); cout << left << (*lig)->name;
                    cout << ":  RMSD = "; cout.width(12); cout << left <<
                    (*lig)->get_bk_rmsd(rstructure->ligands[0],false,debug_mode);
                    cout << "RMSD_with_H = "; cout.width(12); cout << left << 
                    (*lig)->get_bk_rmsd(rstructure->ligands[0],true,false,false,false);
                    if (debug_mode) {cout << "\n@ "; cout << left << (*lig)->name; cout << ":  ";}
                    cout << "RMSD_opt = " << (*lig)->get_bk_rmsd(rstructure->ligands[0],false,debug_mode,true,false) << "\n";
                }
                instruct->clear();
                delete rstructure;
                delete rparser;
                continue;
            }
            if (rmsd3) {
                if (sourcefile == "X") {
                    cerr << c_message<cERROR>("you have to specify a reference file by --s=your_ref.mol2") << endl;
                    exit(1);
                }
                STRUCTURE *rstructure;
                PARSER *rparser;
                rstructure = new STRUCTURE();
                rparser = new PARSER(rstructure,verb,true,true,true,parse_comments);
                if (!(rparser->read_mol2(sourcefile.c_str()))) {
                    cerr << c_message<cERROR>("could not read '") << sourcefile << "' as mol2 file" << endl;
                    delete rstructure;
                    delete rparser;
                    exit(1);
                }
                for (ligands_vec lig=instruct->ligands.begin(); lig!=instruct->ligands.end(); ++lig) {
                    vector<string> vrn;
                    vector<float> vtr; vector<float> otr;
                    vector<float> vbr; vector<float> obr;
                    vector<float> vrr; vector<float> orr;
                    (*lig)->rmsd_by_mol2_residues(vrn,vtr,vbr,vrr,otr,obr,orr,rstructure->ligands[0],debug_mode);
                    cout << "@ " << *fit << " : " << (*lig)->name << "\n";
                    for (unsigned int i=0; i<vrn.size(); ++i) {
                        cout << vrn[i] << " :  tot = " << vtr[i] << "  back = " << vbr[i] << "  res = " << vrr[i] 
                             << "  o_tot = " << otr[i] << "  o_back = " << obr[i] << "  o_res = " << orr[i] << "\n";
                    }
                }
                instruct->clear();
                delete rstructure;
                delete rparser;
                continue;
            }
            if (get_vol) {
                for (ligands_vec lig=instruct->ligands.begin(); lig!=instruct->ligands.end(); ++lig) {
                    cout.width(30); cout << left << (*lig)->name;
                    if (debug_mode) {
                        cout << ":  vol = "; cout.width(12); cout << left <<
                        (*lig)->get_vdw_volume_grid() << " A^3\n";
                    } else {
                        cout << ":  vol = "; cout.width(12); cout << left <<
                        (*lig)->get_vdw_volume() << " A^3\n";
                    }
                }
                instruct->clear();
                continue;
            }
            if (cluster) {
                if (go_on) {
                    cerr << c_message<cERROR>("too many molecules in ") << sourcefile << " for hierarchical clustering" << endl;
                    exit(1);
                }
                //! Hierarchisch mit complete linkage nach RMSD clustern:
                if (!cluster_rmsds_ordered) for (ligands_vec lig=instruct->ligands.begin(); lig!=instruct->ligands.end(); ++lig) {
                    (*lig)->get_atom_typing(atom_typing,true,def_file.c_str(),true,verb,false,false,"X",max_ring_members,
                                            true,0,0,prot_acids,prot_guanidin,prot_amidin,prot_amin,prot_phosphate,
                                            prot_sulfate,kekulize_aromatics,kekulize_charged,allow_charged_aromatics);
                }
                
                cout << "clustering " << instruct->ligands.size() << " ligands..." << endl;
                
                HIERA_CLUST<stl_ptr<LIGAND> > myclust(instruct->ligands,rmsd_dist_function);
                if (cluster_rmsds_ordered) myclust.get_elem_distance = rmsd_dist_ordered;
            
            //!    DEBUG:
            ///    HIERA_CLUST<stl_ptr<LIGAND> > myclust("test.dat");
                
                myclust.cluster_complete_link(cluster_threshold);
                
                if (sourcefile != "X") {
                    if (myclust.clusters.size() > 1) {
                        cout << "warning: due to the given low threshold of " << cluster_threshold << " angstrom the "
                             << "clustering and thus the dendrogram is incomplete" << endl;
                    }
                    
                    string dendro_name = "dendrogram_" + sourcefile;
                    cout << "writing dendrogram as " << dendro_name << endl;
                    myclust.write_dendrogram(dendro_name.c_str(),get_clust_label);
                    
                    string matrix_name = "matrix_" + sourcefile;
                    cout << "writing similarity matrix as " << matrix_name << endl;
                    myclust.write_matrix(matrix_name.c_str(),get_clust_label);
                    
                //!    DEBUG:
                ///    myclust.write_matrix_values("test.dat",get_clust_label);
                    
                }
                
                if (pdb_source != "X") {
                    cout << "writing cluster file..." << endl;
                    myclust.write_clusters(pdb_source.c_str(),get_clust_label);
                }
                
                if (write_full_clusters) {
                    int ncl = 0;
                    for (HIERA_CLUST<stl_ptr<LIGAND> >::cluster_iter it=myclust.clusters.begin();
                                                                     it!=myclust.clusters.end(); ++it) {
                        vector<stl_ptr<LIGAND> > ligs2write;
                        HIERA_CLUST<stl_ptr<LIGAND> >::ele_ptr cele = (*it)->elements;
                        while (cele) {
                            ligs2write.push_back(instruct->ligands[cele->element]);
                            cele = cele->next;
                        }
                        
                        ++ncl;
                        cout << " -> cluster " << ncl << " :  ligands = " << (*it)->n_elements << endl;
                        
                        if (sourcefiles.size() > 1) {
                            tfile = *fit;
                            if (targetfile != "X") {
                                ostringstream os;
                                os << "_" << targetfile << "_cluster" << ncl << ".mol2";
                                replace_ext(tfile,".mol2",os.str());
                            }
                        } else {
                            if (targetfile != "X") tfile = targetfile;
                            else {
                                tfile = *fit;
                                ostringstream os;
                                os << "_cluster" << ncl << ".mol2";
                                replace_ext(tfile,".mol2",os.str());
                            }
                        }
                        string_fu::remove_char(tfile);
                        inparser->write_mol2(ligs2write,tfile.c_str());
                    }
                    
                    cout << instruct->ligands.size() << " ligands merged to " << ncl << " clusters" << endl;
                } else if (cluster_without_mol2 == false) {
                    //! Jetzt in jedem Cluster die "mittlere" Struktur bestimmen und diese rausschreiben:
                    //! -> einfach die Struktur nehmen, die insgesamt die geringsten RMSDs zu allen anderen hat:
                    vector<stl_ptr<LIGAND> > ligs2write;
                    int ncl = 0;
                    for (HIERA_CLUST<stl_ptr<LIGAND> >::cluster_iter it=myclust.clusters.begin();
                                                                     it!=myclust.clusters.end(); ++it) {
                        int to_take = (*it)->elements->element;
                        float min_sum = 999999999.;
                        HIERA_CLUST<stl_ptr<LIGAND> >::ele_ptr cele = (*it)->elements;
                        while (cele) {
                            float rsum = 0.;
                            HIERA_CLUST<stl_ptr<LIGAND> >::ele_ptr cele2 = (*it)->elements;
                            while (cele2) {
                                if (cele == cele2) {
                                    cele2 = cele2->next;
                                    continue;
                                }
                                rsum += myclust.get_wert(cele->element,cele2->element);
                                cele2 = cele2->next;
                            }
                            if (rsum < min_sum) {
                                to_take = cele->element;
                                min_sum = rsum;
                            }
                            cele = cele->next;
                        }
                        
                        if ((*it)->n_elements > 1) min_sum /= (*it)->n_elements - 1;
                        else min_sum = 0.;
                        ++ncl;
                        cout << " -> cluster " << ncl << " :  ligands = " << (*it)->n_elements
                             << "   mean rmsd = " << min_sum << endl;
                        ligs2write.push_back(instruct->ligands[to_take]);
                    }
                    
                    if (sourcefiles.size() > 1) {
                        tfile = *fit;
                        if (targetfile != "X") {
                            ostringstream os;
                            os << "_" << targetfile << ".mol2";
                            replace_ext(tfile,".mol2",os.str());
                        }
                    } else {
                        if (targetfile != "X") tfile = targetfile;
                        else {
                            tfile = *fit;
                            replace_ext(tfile,".mol2","_clustered.mol2");
                        }
                    }
                    string_fu::remove_char(tfile);
                    inparser->write_mol2(ligs2write,tfile.c_str());
                    
                    cout << instruct->ligands.size() << " ligands merged to " << ligs2write.size() << " clusters" << endl;
                }
                
                instruct->clear();
                continue;
            }
            if (clashes) {
                if (sourcefile == "X") {
                    cerr << c_message<cERROR>("you have to specify a reference file by --s=your_ref.mol2") << endl;
                    exit(1);
                }
                STRUCTURE *rstructure;
                PARSER *rparser;
                rstructure = new STRUCTURE();
                rparser = new PARSER(rstructure,verb,true,true,true,parse_comments);
                if (!(rparser->read_mol2(sourcefile.c_str()))) {
                    cerr << c_message<cERROR>("could not read '") << sourcefile << "' as mol2 file" << endl;
                    delete rstructure;
                    delete rparser;
                    exit(1);
                }
                if (rstructure->ligands.size() == 0) {
                    cerr << c_message<cERROR>("found no molecule in '") << sourcefile << endl;
                    delete rstructure;
                    delete rparser;
                    exit(1);
                }
                
                
                atom_properties::initialize();
                rstructure->ligands[0]->get_elements(false);
                
            //    if (instruct->ligands.size() > 1) {
                    //!========================================
                    //! Raumaufteilung fuer Performance:
                    vector<stl_ptr<ATOM> > grid_atoms;
                    for (atoms_vec bt=rstructure->ligands[0]->atoms.begin();
                            bt!=rstructure->ligands[0]->atoms.end(); ++bt) {
                        if ((*bt)->element == "H") continue;
                        grid_atoms.push_back(*bt);
                    }
                    float spacing;
                    
                    if (instruct->ligands.size() < 10) spacing = 5.;
                    else if (instruct->ligands.size() < 50) spacing = 4.;
                    else if (instruct->ligands.size() < 300) spacing = 3.;
                    else spacing = 2.;
                    float hsp = sqrt(3. * spacing * spacing);
                    hsp = (hsp/2.) + 0.1;
                    COORD_GRID pro_grid(grid_atoms,spacing,4.4+hsp);
                    map<int,vector<int> > rel_pro;
                    vec3d<float> g_vec;
                    float max_lim = 4.4 + hsp;
                    max_lim *= max_lim;
                    for (int indi=0; indi<pro_grid.grid.size(); ++indi) {
                        g_vec = pro_grid.grid.dvalue(indi); //Gitterpunkt vom Proteingrid
                        int pro_count = 0;
                        for (atoms_vec apro=grid_atoms.begin(); apro!=grid_atoms.end(); ++apro) {
                            if (get_square_distance((*apro)->coord,g_vec) <= max_lim) {
                                rel_pro[indi].push_back(pro_count);
                            }
                            pro_count++;
                        }
                    }
                    int g_ind[3]; int indi;
                    for (ligands_vec lig=instruct->ligands.begin(); lig!=instruct->ligands.end(); ++lig) {
                        (*lig)->get_elements(false);
                        float overlap = 0.; int nat = 0;
                        for (atoms_vec at=(*lig)->atoms.begin(); at!=(*lig)->atoms.end(); ++at) {
                            
                            if ((*at)->element == "H") continue;
                            ++nat;
                            
                            g_ind[0] = int((((*at)->coord[0] - pro_grid.min[0]) / spacing) + 0.5);
                            g_ind[1] = int((((*at)->coord[1] - pro_grid.min[1]) / spacing) + 0.5);
                            g_ind[2] = int((((*at)->coord[2] - pro_grid.min[2]) / spacing) + 0.5);
                            indi = pro_grid.grid.get_direct_index(g_ind);
                            if (indi < 0 || indi > pro_grid.grid.size()) continue;
                            
                            for (vector<int>::iterator apro=rel_pro[indi].begin(); apro!=rel_pro[indi].end(); ++apro) {
                                overlap += get_sphere_overlap(atom_properties::vdW_map[(*at)->element],
                                    atom_properties::vdW_map[grid_atoms[*apro]->element],
                                    (*at)->coord,grid_atoms[*apro]->coord);
                            }
                        }
                        /*if (overlap > 0.1)*/ cout << "@ " << *fit << " : " << (*lig)->name << "  vdw_overlap = " 
                                    << overlap << " A^3  per_atom_overlap = " 
                                    << overlap/nat << " A^3" << "\n";
                    }
                    //!========================================
                /*} else for (ligands_vec lig=instruct->ligands.begin(); lig!=instruct->ligands.end(); ++lig) {
                    (*lig)->get_elements(false);
                    float overlap = 0.;
                    for (atoms_vec at=(*lig)->atoms.begin(); at!=(*lig)->atoms.end(); ++at) {
                        for (atoms_vec bt=rstructure->ligands[0]->atoms.begin();
                                       bt!=rstructure->ligands[0]->atoms.end(); ++bt) {
                            overlap += get_sphere_overlap(atom_properties::vdW_map[(*at)->element],
                                       atom_properties::vdW_map[(*bt)->element],(*at)->coord,(*bt)->coord);
                        }
                    }
                    if (overlap > 0.1) cout << "@ " << *fit << " : " << (*lig)->name << "  vdw_overlap = " 
                                            << overlap << " A^3" << "\n";
                }*/
                instruct->clear();
                delete rstructure;
                delete rparser;
                continue;
            }
            if (optalign) {
                if (sourcefile == "X") {
                    cerr << c_message<cERROR>("you have to specify a reference file by --s=your_ref.mol2") << endl;
                    exit(1);
                }
                if (targetfile == "X" && go_on) {
                    cerr << c_message<cERROR>("your input file has more than ") << max_mol2_block
                         << " molecules => thus you have to specify a target by --t=name" << endl;
                    exit(1);
                }
                STRUCTURE *rstructure;
                PARSER *rparser;
                rstructure = new STRUCTURE();
                rparser = new PARSER(rstructure,verb,true,true,true,parse_comments);
                if (!(rparser->read_mol2(sourcefile.c_str()))) {
                    cerr << c_message<cERROR>("could not read '") << sourcefile << "' as mol2 file" << endl;
                    delete rstructure;
                    delete rparser;
                    exit(1);
                }
                if (rstructure->ligands.size() == 0) {
                    cerr << c_message<cERROR>("found no molecule in '") << sourcefile << endl;
                    delete rstructure;
                    delete rparser;
                    exit(1);
                }
                rstructure->ligands[0]->get_atom_typing(1,true,"X",true,0,false,false,"X",max_ring_members,false,0,0,prot_acids,prot_guanidin,prot_amidin,prot_amin,prot_phosphate,
                                                   prot_sulfate,kekulize_aromatics,kekulize_charged,allow_charged_aromatics);

                bool spatial = true;
                string pro_mode = "pdb";
                if (pdb_source != "X") {
                    if (instruct->ligands.size() > 1) {
                        cerr << c_message<cERROR>("found more than one molecules in ")
                             << *fit << " => not clear which one to take for the "
                             << "alignment of " << pdb_source << endl;
                        exit(1);
                    }
                    spatial = false;
                    if (pdb_source[pdb_source.size()-1] == '2') pro_mode = "mol";
                }


                for (ligands_vec lig=instruct->ligands.begin(); lig!=instruct->ligands.end(); ++lig) {
                    (*lig)->get_atom_typing(1,true,"X",true,0,false,false,"X",max_ring_members,false,
                                           0,0,prot_acids,prot_guanidin,prot_amidin,prot_amin,prot_phosphate,
                                           prot_sulfate,kekulize_aromatics,kekulize_charged,allow_charged_aromatics);
                    if ((*lig)->get_rmsd(rstructure->ligands[0],false,false,true) < 0.) {
                        cerr << c_message<cWARNING>("found no alingnment for '") << (*lig)->name 
                             << " from " << *fit << endl;
                        continue;
                    }
                    if (spatial) (*lig)->opt_align();
                }
                
                if (sourcefiles.size() > 1) {
                    tfile = *fit;
                    if (targetfile != "X") {
                        ostringstream os;
                        os << "_" << targetfile << ".mol2";
                        replace_ext(tfile,".mol2",os.str());
                    }
                } else {
                    if (targetfile != "X") tfile = targetfile;
                    else {
                        tfile = *fit;
                        if (!spatial) tfile = pdb_source;
                    }
                }

                if (spatial) {
                    string_fu::remove_char(tfile);
                    if (go_on || must_append) {
                        if (int(instruct->ligands.size()) >= max_mol2_block && tfile == *fit) {
                            cerr << c_message<cERROR>("if there are more than ") << max_mol2_block-1
                                 << " molecules in your mol2 (the current blocksize of fconv), you have to specify a target with '--t'" << endl;
                            exit(0);
                        }
                        if (mode == "sdf") {
                            for (ligands_vec lig=instruct->ligands.begin(); lig!=instruct->ligands.end(); ++lig) {
                                (*lig)->rename_atoms(false);
                            }
                        }
                        inparser->write_next_mol2(instruct->ligands,tfile.c_str());
                    } else {
                        if (mode == "sdf") {
                            for (ligands_vec lig=instruct->ligands.begin(); lig!=instruct->ligands.end(); ++lig) {
                                (*lig)->rename_atoms(false);
                            }
                        }
                        inparser->write_mol2(tfile.c_str());
                    }
                    if (go_on) must_append = true;
                } else {
                    STRUCTURE *rstructure2;
                    PARSER *rparser2;
                    rstructure2 = new STRUCTURE();
                    rparser2 = new PARSER(rstructure2,0,true,true,true,parse_comments);

                    if (pro_mode == "mol") {
                        if (!(rparser2->read_mol2(pdb_source.c_str()))) {
                            cerr << c_message<cERROR>("could not read '") << pdb_source.c_str() << "' as mol2 file" << endl;
                            delete rstructure2;
                            delete rparser2;
                            exit(1);
                        }
                        for (ligands_vec lig=rstructure2->ligands.begin(); lig!=rstructure2->ligands.end(); ++lig) {
                            for (atoms_vec at=(*lig)->atoms.begin(); at!=(*lig)->atoms.end();++at) {
                                (*at)->coord -= instruct->ligands[0]->optalign_trans[1];
                                (*at)->coord *= instruct->ligands[0]->optalign_rotm;
                                (*at)->coord += instruct->ligands[0]->optalign_trans[0];
                            }
                        }
                        rparser2->write_mol2(tfile.c_str());
                    } else {
                        if (!(rparser2->read_pdb(pdb_source.c_str()))) {
                            cerr << c_message<cERROR>("could not read '") << pdb_source.c_str() << "' as pdb file" << endl;
                            delete rstructure2;
                            delete rparser2;
                            exit(1);
                        }
                        if (rstructure2->protein) for (chains_vec ct=rstructure2->protein->chains.begin(); ct!=rstructure2->protein->chains.end(); ++ct) {
                            for (atoms_vec at=(*ct)->atoms.begin(); at!=(*ct)->atoms.end();++at) {
                                (*at)->coord -= instruct->ligands[0]->optalign_trans[1];
                                (*at)->coord *= instruct->ligands[0]->optalign_rotm;
                                (*at)->coord += instruct->ligands[0]->optalign_trans[0];
                            }
                        }
                        for (ligands_vec lig=rstructure2->ligands.begin(); lig!=rstructure2->ligands.end(); ++lig) {
                            for (atoms_vec at=(*lig)->atoms.begin(); at!=(*lig)->atoms.end();++at) {
                                (*at)->coord -= instruct->ligands[0]->optalign_trans[1];
                                (*at)->coord *= instruct->ligands[0]->optalign_rotm;
                                (*at)->coord += instruct->ligands[0]->optalign_trans[0];
                            }
                        }
                        for (waters_vec wat=rstructure2->waters.begin(); wat!=rstructure2->waters.end(); ++wat) {
                            (*wat)->atom->coord -= instruct->ligands[0]->optalign_trans[1];
                            (*wat)->atom->coord *= instruct->ligands[0]->optalign_rotm;
                            (*wat)->atom->coord += instruct->ligands[0]->optalign_trans[0];
                        }
                        for (metals_vec met=rstructure2->metals.begin(); met!=rstructure2->metals.end(); ++met) {
                            (*met)->atom->coord -= instruct->ligands[0]->optalign_trans[1];
                            (*met)->atom->coord *= instruct->ligands[0]->optalign_rotm;
                            (*met)->atom->coord += instruct->ligands[0]->optalign_trans[0];
                        }
                        rparser2->write_pdb(tfile.c_str());
                    }

                    delete rstructure2;
                    delete rparser2;
                }


                instruct->clear();
                delete rstructure;
                delete rparser;
                continue;
            }
            if (optalign2) {
                if (sourcefile == "X") {
                    cerr << c_message<cERROR>("you have to specify a reference file by --s=your_ref.mol2") << endl;
                    exit(1);
                }
                if (targetfile == "X" && go_on) {
                    cerr << c_message<cERROR>("your input file has more than ") << max_mol2_block
                         << " molecules => thus you have to specify a target by --t=name" << endl;
                    exit(1);
                }
                STRUCTURE *rstructure;
                PARSER *rparser;
                rstructure = new STRUCTURE();
                rparser = new PARSER(rstructure,verb,true,true,true,parse_comments);
                if (!(rparser->read_mol2(sourcefile.c_str()))) {
                    cerr << c_message<cERROR>("could not read '") << sourcefile << "' as mol2 file" << endl;
                    delete rstructure;
                    delete rparser;
                    exit(1);
                }
                if (rstructure->ligands.size() == 0) {
                    cerr << c_message<cERROR>("found no molecule in '") << sourcefile << endl;
                    delete rstructure;
                    delete rparser;
                    exit(1);
                }
                
                rstructure->ligands[0]->get_hybrid_only();


                bool spatial = true;
                string pro_mode = "pdb";
                if (pdb_source != "X") {
                    if (instruct->ligands.size() > 1) {
                        cerr << c_message<cERROR>("found more than one molecules in ")
                             << *fit << " => not clear which one to take for the "
                             << "alignment of " << pdb_source << endl;
                        exit(1);
                    }
                    spatial = false;
                    if (pdb_source[pdb_source.size()-1] == '2') pro_mode = "mol";
                }


                for (ligands_vec lig=instruct->ligands.begin(); lig!=instruct->ligands.end(); ++lig) {
                    (*lig)->get_atom_typing(1,true,"X",true,0,false,false,"X",max_ring_members,
                                           true,0,0,prot_acids,prot_guanidin,prot_amidin,prot_amin,prot_phosphate,
                                           prot_sulfate,kekulize_aromatics,kekulize_charged,allow_charged_aromatics);
                    if((*lig)->get_bk_rmsd(rstructure->ligands[0],false,debug_mode,true,false) < 0.) {
                        cerr << c_message<cWARNING>("found no alingnment for '") << (*lig)->name 
                             << " from " << *fit << endl;
                        continue;
                    }
                    if (spatial) (*lig)->opt_align();
                }
                
                if (sourcefiles.size() > 1) {
                    tfile = *fit;
                    if (targetfile != "X") {
                        ostringstream os;
                        os << "_" << targetfile << ".mol2";
                        replace_ext(tfile,".mol2",os.str());
                    }
                } else {
                    if (targetfile != "X") tfile = targetfile;
                    else {
                        tfile = *fit;
                        if (!spatial) tfile = pdb_source;
                    }
                }

                if (spatial) {
                    string_fu::remove_char(tfile);
                    if (go_on || must_append) {
                        if (int(instruct->ligands.size()) >= max_mol2_block && tfile == *fit) {
                            cerr << c_message<cERROR>("if there are more than ") << max_mol2_block-1
                                 << " molecules in your mol2 (the current blocksize of fconv), you have to specify a target with '--t'" << endl;
                            exit(0);
                        }
                        if (mode == "sdf") {
                            for (ligands_vec lig=instruct->ligands.begin(); lig!=instruct->ligands.end(); ++lig) {
                                (*lig)->rename_atoms(false);
                            }
                        }
                        inparser->write_next_mol2(instruct->ligands,tfile.c_str());
                    } else {
                        if (mode == "sdf") {
                            for (ligands_vec lig=instruct->ligands.begin(); lig!=instruct->ligands.end(); ++lig) {
                                (*lig)->rename_atoms(false);
                            }
                        }
                        inparser->write_mol2(tfile.c_str());
                    }
                    if (go_on) must_append = true;
                } else {
                    STRUCTURE *rstructure2;
                    PARSER *rparser2;
                    rstructure2 = new STRUCTURE();
                    rparser2 = new PARSER(rstructure2,0,true,true,true,parse_comments);

                    if (pro_mode == "mol") {
                        if (!(rparser2->read_mol2(pdb_source.c_str()))) {
                            cerr << c_message<cERROR>("could not read '") << pdb_source.c_str() << "' as mol2 file" << endl;
                            delete rstructure2;
                            delete rparser2;
                            exit(1);
                        }
                        for (ligands_vec lig=rstructure2->ligands.begin(); lig!=rstructure2->ligands.end(); ++lig) {
                            for (atoms_vec at=(*lig)->atoms.begin(); at!=(*lig)->atoms.end();++at) {
                                (*at)->coord -= instruct->ligands[0]->optalign_trans[1];
                                (*at)->coord *= instruct->ligands[0]->optalign_rotm;
                                (*at)->coord += instruct->ligands[0]->optalign_trans[0];
                            }
                        }
                        rparser2->write_mol2(tfile.c_str());
                    } else {
                        if (!(rparser2->read_pdb(pdb_source.c_str()))) {
                            cerr << c_message<cERROR>("could not read '") << pdb_source.c_str() << "' as pdb file" << endl;
                            delete rstructure2;
                            delete rparser2;
                            exit(1);
                        }
                        if (rstructure2->protein) for (chains_vec ct=rstructure2->protein->chains.begin(); ct!=rstructure2->protein->chains.end(); ++ct) {
                            for (atoms_vec at=(*ct)->atoms.begin(); at!=(*ct)->atoms.end();++at) {
                                (*at)->coord -= instruct->ligands[0]->optalign_trans[1];
                                (*at)->coord *= instruct->ligands[0]->optalign_rotm;
                                (*at)->coord += instruct->ligands[0]->optalign_trans[0];
                            }
                        }
                        for (ligands_vec lig=rstructure2->ligands.begin(); lig!=rstructure2->ligands.end(); ++lig) {
                            for (atoms_vec at=(*lig)->atoms.begin(); at!=(*lig)->atoms.end();++at) {
                                (*at)->coord -= instruct->ligands[0]->optalign_trans[1];
                                (*at)->coord *= instruct->ligands[0]->optalign_rotm;
                                (*at)->coord += instruct->ligands[0]->optalign_trans[0];
                            }
                        }
                        for (waters_vec wat=rstructure2->waters.begin(); wat!=rstructure2->waters.end(); ++wat) {
                            (*wat)->atom->coord -= instruct->ligands[0]->optalign_trans[1];
                            (*wat)->atom->coord *= instruct->ligands[0]->optalign_rotm;
                            (*wat)->atom->coord += instruct->ligands[0]->optalign_trans[0];
                        }
                        for (metals_vec met=rstructure2->metals.begin(); met!=rstructure2->metals.end(); ++met) {
                            (*met)->atom->coord -= instruct->ligands[0]->optalign_trans[1];
                            (*met)->atom->coord *= instruct->ligands[0]->optalign_rotm;
                            (*met)->atom->coord += instruct->ligands[0]->optalign_trans[0];
                        }
                        rparser2->write_pdb(tfile.c_str());
                    }

                    delete rstructure2;
                    delete rparser2;
                }

                instruct->clear();
                delete rstructure;
                delete rparser;
                continue;
            }
            
            if (search_substruct) {
                if (sourcefile == "X") {
                    cerr << c_message<cERROR>("you have to specify a reference file by --s=your_ref.mol2") << endl;
                    exit(1);
                }
                STRUCTURE *rstructure;
                PARSER *rparser;
                rstructure = new STRUCTURE();
                rparser = new PARSER(rstructure,verb,true,true,true,parse_comments);
                if (!(rparser->read_mol2(sourcefile.c_str()))) {
                    cerr << c_message<cERROR>("could not read '") << sourcefile << "' as mol2 file" << endl;
                    delete rstructure;
                    delete rparser;
                    exit(1);
                }
                if (rstructure->ligands.size() == 0) {
                    cerr << c_message<cERROR>("found no molecule in '") << sourcefile << endl;
                    delete rstructure;
                    delete rparser;
                    exit(1);
                }
                
                rstructure->ligands[0]->get_hybrid_only();
                
                vector<stl_ptr<LIGAND> > ligs2write;
                for (ligands_vec lig=instruct->ligands.begin(); lig!=instruct->ligands.end(); ++lig) {
                    if ((*lig)->atoms.size() < rstructure->ligands[0]->atoms.size()) continue;
                    if (targetfile != "X") (*lig)->get_atom_typing(atom_typing,true,def_file.c_str(),true,
                                                                   0,false,false,"X",max_ring_members,
                                                                   true,0,0,prot_acids,prot_guanidin,prot_amidin,prot_amin,
                                                                   prot_phosphate,prot_sulfate,kekulize_aromatics,kekulize_charged,
                                                                   allow_charged_aromatics);
                    else (*lig)->get_hybrid_only();
                    int n_sub = (*lig)->has_substructure(rstructure->ligands[0],ss_with_exact_hybrid);
                    if (n_sub > 0) {
                        if (targetfile != "X") {
                            ligs2write.push_back(*lig);
                        } else {
                            cout << *fit << " :  " << (*lig)->name << "  " << n_sub << endl;
                        }
                    }
                }
                if (targetfile != "X") {
                    string_fu::remove_char(targetfile);
                    if (mode == "sdf") {
                        for (ligands_vec lig=ligs2write.begin(); lig!=ligs2write.end(); ++lig) {
                            (*lig)->rename_atoms(false);
                        }
                    }
                    inparser->write_next_mol2(ligs2write,targetfile.c_str());
                    if (go_on) must_append = true;
                }
                
                instruct->clear();
                delete rstructure;
                delete rparser;
                continue;
            }
            
            if (ss_with_exact_type) {
                if (sourcefile == "X") {
                    cerr << c_message<cERROR>("you have to specify a reference file by --s=your_ref.mol2") << endl;
                    exit(1);
                }
                STRUCTURE *rstructure;
                PARSER *rparser;
                rstructure = new STRUCTURE();
                rparser = new PARSER(rstructure,verb,true,true,true,parse_comments);
                if (!(rparser->read_mol2(sourcefile.c_str()))) {
                    cerr << c_message<cERROR>("could not read '") << sourcefile << "' as mol2 file" << endl;
                    delete rstructure;
                    delete rparser;
                    exit(1);
                }
                if (rstructure->ligands.size() == 0) {
                    cerr << c_message<cERROR>("found no molecule in '") << sourcefile << endl;
                    delete rstructure;
                    delete rparser;
                    exit(1);
                }
                
                if (def_file != "X") {
                    rstructure->ligands[0]->get_atom_typing(atom_typing,true,def_file.c_str(),true,
                                                            0,false,false,"X",max_ring_members,
                                                            true,0,0,prot_acids,prot_guanidin,prot_amidin,prot_amin,
                                                            prot_phosphate,prot_sulfate,
                                                            kekulize_aromatics,kekulize_charged,allow_charged_aromatics);
                } else rstructure->ligands[0]->get_hybrid_only();
                
                vector<stl_ptr<LIGAND> > ligs2write;
                for (ligands_vec lig=instruct->ligands.begin(); lig!=instruct->ligands.end(); ++lig) {
                    if ((*lig)->atoms.size() < rstructure->ligands[0]->atoms.size()) continue;
                    
                    if (targetfile != "X" || def_file != "X") {
                        if (def_file != "X") {
                            for (atoms_vec at=(*lig)->atoms.begin(); at!=(*lig)->atoms.end(); ++at) {
                                (*at)->full_res_name = (*at)->sybyl_type;
                            }
                        }
                        (*lig)->get_atom_typing(atom_typing,true,def_file.c_str(),true,
                                                    0,false,false,"X",max_ring_members,
                                                    true,0,0,prot_acids,prot_guanidin,prot_amidin,prot_amin,prot_phosphate,
                                                    prot_sulfate,kekulize_aromatics,kekulize_charged,allow_charged_aromatics);
                    }
                    else (*lig)->get_hybrid_only();
                    int n_sub;
                    if (def_file == "X") n_sub = (*lig)->has_substructure(rstructure->ligands[0],ss_with_exact_hybrid,true);
                    else n_sub = (*lig)->has_substructure(rstructure->ligands[0],false,true,true);
                    if (n_sub > 0) {
                        if (targetfile != "X") {
                            ligs2write.push_back(*lig);
                        } else {
                            cout << *fit << " :  " << (*lig)->name << "  " << n_sub << endl;
                        }
                    }
                    if (def_file != "X") {
                        for (atoms_vec at=(*lig)->atoms.begin(); at!=(*lig)->atoms.end(); ++at) {
                            (*at)->sybyl_type = (*at)->full_res_name;
                        }
                    }
                }
                if (targetfile != "X") {
                    string_fu::remove_char(targetfile);
                    if (mode == "sdf") {
                        for (ligands_vec lig=ligs2write.begin(); lig!=ligs2write.end(); ++lig) {
                            (*lig)->rename_atoms(false);
                        }
                    }
                    inparser->write_next_mol2(ligs2write,targetfile.c_str());
                    if (go_on) must_append = true;
                }
                
                instruct->clear();
                delete rstructure;
                delete rparser;
                continue;
            }
            
            if (calc_optalign_matrix) {
                if (sourcefile == "X") {
                    cerr << c_message<cERROR>("you have to specify a reference file by --s=your_ref.mol2") << endl;
                    exit(1);
                }
                STRUCTURE *rstructure;
                PARSER *rparser;
                rstructure = new STRUCTURE();
                rparser = new PARSER(rstructure,verb,true,true,true,parse_comments);
                if (!(rparser->read_mol2(sourcefile.c_str()))) {
                    cerr << c_message<cERROR>("could not read '") << sourcefile << "' as mol2 file" << endl;
                    delete rstructure;
                    delete rparser;
                    exit(1);
                }
                if (rstructure->ligands.size() == 0) {
                    cerr << c_message<cERROR>("found no molecule in '") << sourcefile << endl;
                    delete rstructure;
                    delete rparser;
                    exit(1);
                }
                rstructure->ligands[0]->get_atom_typing(1,true,"X",true,0,false,false,"X",max_ring_members,
                                                   true,0,0,prot_acids,prot_guanidin,prot_amidin,prot_amin,prot_phosphate,
                                                   prot_sulfate,kekulize_aromatics,kekulize_charged,allow_charged_aromatics);
                for (ligands_vec lig=instruct->ligands.begin(); lig!=instruct->ligands.end(); ++lig) {
                    (*lig)->get_atom_typing(1,true,"X",true,0,false,false,"X",max_ring_members,
                                                   true,0,0,prot_acids,prot_guanidin,prot_amidin,prot_amin,prot_phosphate,
                                                   prot_sulfate,kekulize_aromatics,kekulize_charged,allow_charged_aromatics);
                    if((*lig)->get_rmsd(rstructure->ligands[0],false,false,true) < 0.) continue;
                    cout << (*lig)->name << " :" << "\n";
                    cout << "1. translation: " << (*lig)->optalign_trans[1] << "\n";
                    cout << "rotation matrix: " << (*lig)->optalign_rotm << "\n";
                    cout << "2. translation: " << (*lig)->optalign_trans[0] << "\n";
                    cout << endl;
                }
                
                instruct->clear();
                delete rstructure;
                delete rparser;
                continue;
            }
            
            if (sdiv) {
                
                clock_t start,ende; //Variablen zum Zeitstoppen
                double zeit,t1;//
                
                start=clock();
                
                for (ligands_vec lig=instruct->ligands.begin(); lig!=instruct->ligands.end(); ++lig) {
                    (*lig)->get_atom_typing(0,true,"X",false,0,false,false,"X",max_ring_members,
                                                   true,0,0,prot_acids,prot_guanidin,prot_amidin,prot_amin,prot_phosphate,
                                                   prot_sulfate,kekulize_aromatics,kekulize_charged,allow_charged_aromatics);
                }
                vec3d<float> svbuf = instruct->get_ligand_diversity();
                cout << "\n" << *fit << " :" << "\n";
                cout << "   -> mean rmsd          : " << svbuf[0] << "\n";
                cout << "   -> standard deviation : " << svbuf[1] << "\n";
                cout << "   -> similar molecules  : " << int(svbuf[2]) << "\n";
                
                ende=clock();//
                zeit=(ende-start);//
                t1 = zeit/CLOCKS_PER_SEC;//
                cout << "   -> calculation time   : " << t1 << " s" << endl;
                
                instruct->clear();
                continue;
            }
            if (get_n_rings) {
                cout << "@ " << *fit << ":\n";
                for (ligands_vec lig=instruct->ligands.begin(); lig!=instruct->ligands.end(); ++lig) {
                    int n_r = (*lig)->get_number_of_rings(max_ring_members,debug_mode);


                    tr1::unordered_map<int,tr1::unordered_set<int> > fused_map;
                    for (unsigned int i=0; i<(*lig)->rings.size(); ++i) {
                        bool conti = false;
                        for (unsigned int k=0; k<i; ++k) {
                            if (fused_map.find(k) != fused_map.end()) {
                                if (fused_map[k].find(i) != fused_map[k].end()) {
                                    conti = true;
                                    break;
                                }
                            }
                        }
                        if (conti) continue;
                        fused_map[i] = tr1::unordered_set<int>();
                        fused_map[i].insert(i);
                        bool added = false;
                        for (unsigned int j=i+1; j<(*lig)->rings.size(); ++j) {
                            if ((*lig)->rings[i]->fused((*lig)->rings[j])) {
                                fused_map[i].insert(j);
                                added = true;
                            }
                        }
                        while (added) {
                            added = false;
                            for (unsigned int j=i+1; j<(*lig)->rings.size(); ++j) {
                                if (fused_map[i].find(j) != fused_map[i].end()) continue;
                                for (tr1::unordered_set<int>::iterator it=fused_map[i].begin(); it!=fused_map[i].end(); ++it) {
                                    if ((*lig)->rings[*it]->fused((*lig)->rings[j])) {
                                        added = true;
                                        fused_map[i].insert(j);
                                        break;
                                    }
                                }
                            }
                        }
                    }

                    cout.width(30); cout << left << (*lig)->name;
                    cout << ":  n_rings = "; cout.clear(); cout << n_r;
                    cout << "   n_fused_rings = " << fused_map.size() << "  ( ";
                    for (tr1::unordered_map<int,tr1::unordered_set<int> >::iterator mm=fused_map.begin();
                                                                                    mm!=fused_map.end(); ++mm) {
                        cout << mm->second.size() << " ";
                    }
                    cout << ")\n";
                }
                instruct->clear();
                continue;
            }
            if (get_n_rings2) {
                cout << "@ " << *fit << ":\n";
                for (ligands_vec lig=instruct->ligands.begin(); lig!=instruct->ligands.end(); ++lig) {
                    (*lig)->get_number_of_rings(max_ring_members,debug_mode);
                    cout.width(30); cout << left << (*lig)->name;
                    set<int> id_set;
                    for (rings_vec rt=(*lig)->rings.begin(); rt!=(*lig)->rings.end(); ++rt) {
                        for (atoms_vec at=(*rt)->ring_atoms.begin(); at!=(*rt)->ring_atoms.end(); ++at) {
                            id_set.insert((*at)->intern_id);
                        }
                    }
                    cout << ":  ring_atoms = " << id_set.size() 
                         << "  other_heavy_atoms = " << (*lig)->get_n_heavy()-id_set.size() 
                         << "  n_rings = " << (*lig)->rings.size() << "\n";
                }
                instruct->clear();
                continue;
            }
            if (get_mol_weight) {
                for (ligands_vec lig=instruct->ligands.begin(); lig!=instruct->ligands.end(); ++lig) {
                    cout.width(30); cout << left << (*lig)->name;
                    cout << ":  weight = "; cout.width(12); cout << left <<
                    (*lig)->get_mol_weight() << "\n";
                }
                instruct->clear();
                continue;
            }
            if (build_crystal) {
                //int lig_num = 0;
                int lig_num = 0 + total_mol_read - n_mol_read; //! ***
                
                for (ligands_vec lig=instruct->ligands.begin(); lig!=instruct->ligands.end(); ++lig) {
                    lig_num++;
                    CRYSTALIZER *mycrys;
                    mycrys = new CRYSTALIZER(*lig,debug_mode);
                    
                //    if (unit_cell_multi > 0) mycrys->build_crystal2(unit_cell_multi);
                    if (unit_cell_multi > 0) mycrys->build_crystal2_splitted(unit_cell_multi);
                    else mycrys->build_crystal_splitted(unit_cell_exp_rad);
                //    else mycrys->build_crystal(unit_cell_exp_rad);
                    
                    PARSER *pars2;
                    pars2 = new PARSER(mycrys->structure,1,true,true,true,parse_comments);
                
                    if (instruct->ligands.size() > 1) {
                        tfile = *fit;
                        if (targetfile == "X") {
                            ostringstream os;
                            os << "_" << lig_num << "_crys.mol2";
                            replace_ext(tfile,".mol2",os.str());
                        } else {
                            tfile = targetfile;
                            ostringstream os;
                            os << "_" << lig_num << ".mol2";
                            replace_ext(tfile,".mol2",os.str());
                        }
                    } else {
                        tfile = *fit;
                        if (targetfile == "X") {
                            ostringstream os;
                            os << "_crys.mol2";
                            replace_ext(tfile,".mol2",os.str());
                        } else {
                            tfile = targetfile;
                        }
                    }
                    
                    for (ligands_vec lig=mycrys->structure->splitted_ligands.begin(); 
                                     lig!=mycrys->structure->splitted_ligands.end(); ++lig) {
                        (*lig)->get_atom_typing(atom_typing,true,def_file.c_str(),true,0,true,false,"X",max_ring_members,
                                                true,0,0,prot_acids,prot_guanidin,prot_amidin,prot_amin,prot_phosphate,
                                                    prot_sulfate,kekulize_aromatics,kekulize_charged,allow_charged_aromatics);
                    }
                    string_fu::remove_char(tfile);
                    pars2->write_mol2(mycrys->structure->splitted_ligands,tfile.c_str());
                //    pars2->write_mol2(tfile.c_str());
                    delete pars2;
                    delete mycrys;
                }
                
                instruct->clear();
                continue;
            }
            if (split_multimol) {
                
            //    int lig_num = 0;
                int lig_num = 0 + total_mol_read - n_mol_read; //! ***
                
                for (ligands_vec lig=instruct->ligands.begin(); lig!=instruct->ligands.end(); ++lig) {
                    ostringstream os;
                    os << "_" << lig_num << ".mol2";
                    lig_num++;
                    if (targetfile != "X") tfile = targetfile;
                    else tfile = *fit;
                    if (mode == "sdf") replace_ext(tfile,".sdf",os.str());
                    else replace_ext(tfile,".mol2",os.str());
                    string_fu::remove_char(tfile);
                    if (mode == "sdf") {
                        (*lig)->get_atom_typing(atom_typing,false,def_file.c_str(),true,verb,true,false,"X",max_ring_members,no_planar_free_rot,0,0,
                                                    prot_acids,prot_guanidin,prot_amidin,prot_amin,prot_phosphate,prot_sulfate,
                                                    kekulize_aromatics,kekulize_charged,allow_charged_aromatics);
                        (*lig)->rename_atoms(false);
                    }
                    inparser->write_mol2_ligand(*lig,tfile.c_str());
                }
                instruct->clear();
                continue;
            }
            if (blocksplit) {
                if (curr_bs_needed == blocksplit_size) {
                    ostringstream os;
                    os << "_" << total_mol_read-n_mol_read+1 << ".mol2";
                    if (targetfile != "X") tfile = targetfile;
                    else tfile = *fit;
                    if (mode == "sdf") replace_ext(tfile,".sdf",os.str());
                    else replace_ext(tfile,".mol2",os.str());
                    string_fu::remove_char(tfile);
                }
                if (mode == "sdf") {
                    for (ligands_vec lig=instruct->ligands.begin(); lig!=instruct->ligands.end(); ++lig) {
                        (*lig)->get_atom_typing(atom_typing,false,def_file.c_str(),true,verb,true,false,"X",max_ring_members,no_planar_free_rot,0,0,
                                                    prot_acids,prot_guanidin,prot_amidin,prot_amin,prot_phosphate,prot_sulfate,
                                                    kekulize_aromatics,kekulize_charged,allow_charged_aromatics);
                        (*lig)->rename_atoms(false);
                    }
                }
                inparser->write_next_mol2(instruct->ligands,tfile.c_str());
                curr_bs_needed -= instruct->ligands.size();
                if (curr_bs_needed == 0 || go_on == false) {
                    inparser->write_next_mol2_end();
                    curr_bs_needed = blocksplit_size;
                    cout << " molecules " << fnmw << " - " << tnmw << " written as " << tfile << endl;
                    fnmw = tnmw + 1;
                }
                instruct->clear();
                continue;
            }
            if (extract_mol2) {
                if (total_mol_read < (extract_number + 1)) break; //! ***
                if (instruct->ligands.size() < (extract_number + 1)) {
                    cerr << c_message<cERROR>((*fit).c_str()) << " has only " << instruct->ligands.size() << " molecules" << endl;
                    instruct->clear();
                    continue;
                }
                ostringstream os;
                os << "_" << extract_number << ".mol2";
                if (targetfile != "X") tfile = targetfile;
                else {
                    tfile = *fit;
                    if (mode == "sdf") replace_ext(tfile,".sdf",os.str());
                    else replace_ext(tfile,".mol2",os.str());
                }
                string_fu::remove_char(tfile);
                if (mode == "sdf") {
                    instruct->ligands[extract_number]->get_atom_typing(atom_typing,false,def_file.c_str(),true,verb,true,false,"X",max_ring_members,no_planar_free_rot,0,0,
                                                prot_acids,prot_guanidin,prot_amidin,prot_amin,prot_phosphate,prot_sulfate,
                                                kekulize_aromatics,kekulize_charged,allow_charged_aromatics);
                    instruct->ligands[extract_number]->rename_atoms(false);
                }
                inparser->write_mol2_ligand(instruct->ligands[extract_number],tfile.c_str());
                instruct->clear();
                continue;
            }
            if (extract_mol2_by_name) {
                vector<stl_ptr<LIGAND> > tligs;
                for (ligands_vec lig=instruct->ligands.begin(); lig!=instruct->ligands.end(); ++lig) {
                ///    string::size_type i1 = (*lig)->name.find(extract_name);
                ///    if (i1 != string::npos) {
                ///        tligs.push_back(*lig);
                ///    }
                    if ((*lig)->name == extract_name) tligs.push_back(*lig);
                }
                ostringstream os;
                os << "_" << extract_name << ".mol2";
                if (targetfile != "X") tfile = targetfile;
                else {
                    tfile = *fit;
                    if (mode == "sdf") replace_ext(tfile,".sdf",os.str());
                    else replace_ext(tfile,".mol2",os.str());
                }
                string_fu::remove_char(tfile);
                if (mode == "sdf") {
                    for (ligands_vec lig=tligs.begin(); lig!=tligs.end(); ++lig) {
                        (*lig)->get_atom_typing(atom_typing,false,def_file.c_str(),true,verb,true,false,"X",max_ring_members,no_planar_free_rot,0,0,
                                                    prot_acids,prot_guanidin,prot_amidin,prot_amin,prot_phosphate,prot_sulfate,
                                                    kekulize_aromatics,kekulize_charged,allow_charged_aromatics);
                        (*lig)->rename_atoms(false);
                    }
                }
                inparser->write_mol2(tligs,tfile.c_str());
                instruct->clear();
                continue;
            }

            if ((!gethyd) || (!getwat) || convert || retype || add_hyd || rename_atoms || rename_molecules || rewrite ||
                  remove_duplicates || intern_split || reorder_atoms || keep_biggest_fragment || remove_duplicate_atoms) {
                vector<stl_ptr<LIGAND> > first_ligs;

                if (convert) {
                    for (ligands_vec lig=instruct->ligands.begin(); lig!=instruct->ligands.end(); ++lig) {
                        (*lig)->get_atom_typing(atom_typing,false,def_file.c_str(),true,verb,true,false,"X",max_ring_members,no_planar_free_rot,0,0,
                                                prot_acids,prot_guanidin,prot_amidin,prot_amin,prot_phosphate,prot_sulfate,
                                                kekulize_aromatics,kekulize_charged,allow_charged_aromatics);
                        (*lig)->rename_atoms(false);
                        
                    }
                }

                if (retype || add_hyd || remove_duplicates || intern_split || keep_biggest_fragment) {
                    if (targetfile == "X" && go_on && !remove_duplicates) {
                        cerr << c_message<cERROR>("your input file has more than ") << max_mol2_block
                        << " molecules => thus you have to specify a target by --t=name" << endl;
                        exit(1);
                    }
                    STRUCTURE *refstructure = 0;
                    PARSER *refparser = 0;
                    LIGAND* at_ref_mol = 0;
                    if (pdb_source != "X") {
                        refstructure = new STRUCTURE();
                        refparser = new PARSER(refstructure,0,true,true,true,parse_comments);
                        if (!(refparser->read_mol2(pdb_source.c_str(),true))) {
                            cerr << c_message<cERROR>("could not read '") << pdb_source << "' as mol2 file" << endl;
                            delete refstructure;
                            delete refparser;
                            exit(1);
                        }
                        if (refstructure->ligands.size() == 0) {
                            cerr << c_message<cERROR>("found no molecule in '") << pdb_source << endl;
                            delete refstructure;
                            delete refparser;
                            exit(1);
                        }
                        at_ref_mol = &(*(refstructure->ligands[0]));
                    }
                    for (ligands_vec lig=instruct->ligands.begin(); lig!=instruct->ligands.end(); ++lig) {
                        if (remove_duplicates || intern_split || keep_biggest_fragment) (*lig)->get_atom_typing(atom_typing,true,def_file.c_str(),true,verb,false,false,"X",max_ring_members,
                                                true,0,0,prot_acids,prot_guanidin,prot_amidin,prot_amin,prot_phosphate,
                                                    prot_sulfate);
                        else (*lig)->get_atom_typing(atom_typing,true,def_file.c_str(),true,verb,true,false,"X",max_ring_members,no_planar_free_rot,0,at_ref_mol,prot_acids,prot_guanidin,prot_amidin,prot_amin,prot_phosphate,
                                                     prot_sulfate,kekulize_aromatics,kekulize_charged,allow_charged_aromatics);
                    }
                    if (pdb_source != "X") {
                        delete refstructure;
                        delete refparser;
                    }
                    if (add_hyd) {
                        for (ligands_vec lig=instruct->ligands.begin(); lig!=instruct->ligands.end(); ++lig) {
                            
                            //!DEBUG:
                        //    cerr << (*lig)->name << endl;
                            
                            (*lig)->set_standard_protonation(verb);
                        }
                    }
                    if (intern_split || keep_biggest_fragment) {
                        for (ligands_vec lig=instruct->ligands.begin(); lig!=instruct->ligands.end(); ++lig) {
                            (*lig)->split(keep_biggest_fragment);
                        }
                    }
                    if (remove_duplicates) {
                        //! ToDo: blockweise einlesen und nur die uniquen im Speicher halten
                        int removed = 0;
                        if (instruct->ligands.size() > 0) {
                            instruct->ligands[0]->get_compare_numbers();
                            instruct->ligands[0]->has_stereo = true;
                            first_ligs.push_back(instruct->ligands[0]);
                            for (ligands_vec lig=instruct->ligands.begin()+1; lig!=instruct->ligands.end(); ++lig) {
                                (*lig)->get_compare_numbers();
                                bool conti2 = false;
                                (*lig)->has_stereo = true;
                                for (ligands_vec blig=first_ligs.begin(); blig!=first_ligs.end(); ++blig) {

                                    if ((*lig)->check_equality(*blig,false,check_stereo)) {
                                        
                                        ++removed;
                                        //cerr << (*lig)->name << " == " << (*blig)->name << endl;

                                        conti2 = true;
                                        break;
                                    }
                                }
                                if (conti2) continue;
                                first_ligs.push_back(*lig);
                            }
                        }
                        
                        cout << " -> " << removed << " duplicate structures removed!" << endl;
                    }
                }
                if (rename_atoms) {
                    for (ligands_vec lig=instruct->ligands.begin(); lig!=instruct->ligands.end(); ++lig) {
                        if (retype) (*lig)->rename_atoms(false);
                        else (*lig)->rename_atoms(true);
                    }
                }
                if (reorder_atoms) {
                    for (ligands_vec lig=instruct->ligands.begin(); lig!=instruct->ligands.end(); ++lig) {
                        if (retype) (*lig)->reorder_atoms(false);
                        else (*lig)->reorder_atoms(true);
                    }
                }
                if (rename_molecules) {
                    //int lig_num = 0;
                    int lig_num = 0 + total_mol_read - n_mol_read; //! ***
                    
                    for (ligands_vec lig=instruct->ligands.begin(); lig!=instruct->ligands.end(); ++lig) {
                        if (retype) (*lig)->get_atom_typing(atom_typing,true,def_file.c_str(),true,verb,true,false,"X",max_ring_members,
                                                true,0,0,prot_acids,prot_guanidin,prot_amidin,prot_amin,prot_phosphate,
                                                    prot_sulfate,kekulize_aromatics,kekulize_charged,allow_charged_aromatics);
                        ostringstream os;
                        os << "_" << lig_num;
                        lig_num++;
                        string newname = *fit;
                        replace_ext(newname,".mol2",os.str());
                        (*lig)->name = newname;
                    }
                }
                if (remove_duplicate_atoms) {
                    for (ligands_vec lig=instruct->ligands.begin(); lig!=instruct->ligands.end(); ++lig) {
                        vector<stl_ptr<ATOM> > new_at;
                        for (atoms_vec at=(*lig)->atoms.begin(); at!=(*lig)->atoms.end(); ++at) {
                            bool to_take = true;
                            for (atoms_vec bt=new_at.begin(); bt!=new_at.end(); ++bt) {
                                if (get_square_distance((*at)->coord,(*bt)->coord) < 0.0225) {
                                    to_take = false;
                                    
                                    //! Bindungen von at auf bt umlinken:
                                    for (bonds_vec bv=(*lig)->bonds.begin(); bv!=(*lig)->bonds.end(); ++bv) {
                                        if ((*bv)->from->intern_id == (*at)->intern_id) {
                                            bool change = true;
                                            for (bonds_vec bv2=(*lig)->bonds.begin();
                                                           bv2!=(*lig)->bonds.end(); ++bv2) {
                                                if ((*bv2)->from == *bt && (*bv2)->to == (*bv)->to) {
                                                    change = false;
                                                    break;
                                                }
                                            }
                                            if (change) (*bv)->from = *bt;
                                        } else if ((*bv)->to->intern_id == (*at)->intern_id) {
                                            bool change = true;
                                            for (bonds_vec bv2=(*lig)->bonds.begin();
                                                           bv2!=(*lig)->bonds.end(); ++bv2) {
                                                if ((*bv2)->to == *bt && (*bv2)->from == (*bv)->from) {
                                                    change = false;
                                                    break;
                                                }
                                            }
                                            if (change) (*bv)->to = *bt;
                                        }
                                    }
                                    
                                    break;
                                }
                            }
                            if (to_take) new_at.push_back(*at);
                            else at->kill();
                        }
                        (*lig)->atoms.clear();
                        for (atoms_vec bt=new_at.begin(); bt!=new_at.end(); ++bt) (*lig)->atoms.push_back(*bt);
                    }
                }
            
                if (sourcefiles.size() > 1) {
                    tfile = *fit;
                    if (targetfile != "X") {
                        ostringstream os;
                        os << "_" << targetfile << ".mol2";
                        if (mode == "sdf") replace_ext(tfile,".sdf",os.str());
                        else replace_ext(tfile,".mol2",os.str());
                    } else if (mode == "sdf") replace_ext(tfile,"sdf","mol2");
                } else {
                    if (targetfile != "X") tfile = targetfile;
                    else {
                        tfile = *fit;
                        if (mode == "sdf") replace_ext(tfile,"sdf","mol2");
                    }
                }
                string_fu::remove_char(tfile);
                if (remove_duplicates) {
                    inparser->write_mol2(first_ligs,tfile.c_str());
                } else if (intern_split || keep_biggest_fragment) {
                    if (go_on || must_append) {
                        if (int(instruct->splitted_ligands.size()) >= max_mol2_block && tfile == *fit) {
                            cerr << c_message<cERROR>("if there are more than ") << max_mol2_block-1 
                            << " molecules in your mol2 (the current blocksize of fconv), you have to specify a target with '--t'" << endl;
                            exit(0);
                        }
                        for (ligands_vec lig=instruct->splitted_ligands.begin();
                                         lig!=instruct->splitted_ligands.end(); ++lig) {
                            (*lig)->get_atom_typing(atom_typing,true,def_file.c_str(),true,verb,
                                                    true,false,"X",max_ring_members,
                                                    true,0,0,prot_acids,prot_guanidin,prot_amidin,prot_amin,
                                                    prot_phosphate,prot_sulfate,kekulize_aromatics,kekulize_charged,allow_charged_aromatics);
                        }
                        inparser->write_next_mol2(instruct->splitted_ligands,tfile.c_str());
                    } else {
                        for (ligands_vec lig=instruct->splitted_ligands.begin();
                                         lig!=instruct->splitted_ligands.end(); ++lig) {
                            (*lig)->get_atom_typing(atom_typing,true,def_file.c_str(),true,verb,
                                                    true,false,"X",max_ring_members,
                                                    true,0,0,prot_acids,prot_guanidin,prot_amidin,prot_amin,
                                                    prot_phosphate,prot_sulfate,kekulize_charged,allow_charged_aromatics);
                        }
                        inparser->write_mol2(instruct->splitted_ligands,tfile.c_str());
                    }
                    if (go_on) must_append = true;
                } else {
                    if (go_on || must_append) {
                        if (int(instruct->ligands.size()) >= max_mol2_block && tfile == *fit) {
                            cerr << c_message<cERROR>("if there are more than ") << max_mol2_block-1 
                            << " molecules in your mol2 (the current blocksize of fconv), you have to specify a target with '--t'" << endl;
                            exit(0);
                        }
                        inparser->write_next_mol2(instruct->ligands,tfile.c_str());
                    } else inparser->write_mol2(tfile.c_str());
                    if (go_on) must_append = true;
                }
                instruct->clear();
                continue;
            }
            
            //! *a*
            /*
            if (!(merge_mol2 || add_mol2)) {
                instruct->clear_ligands();
            }
            */
            
            if (merge_mol2 || add_mol2 || merge_mol2_in1mol) {
                if (int(instruct->ligands.size()) >= max_mol2_block && tfile == *fit) {
                    cerr << c_message<cERROR>("if there are more than ") << max_mol2_block-1 
                    << " molecules in your mol2 (the current blocksize of fconv), you have to specify a target with '--t'" << endl;
                    exit(0);
                }
                string_fu::remove_char(tfile);
                if (!merge_mol2_in1mol) {
                    inparser->write_next_mol2(instruct->ligands,tfile.c_str());
                    instruct->clear_ligands();
                }
            } else instruct->clear_ligands();
            
            }
            continue;
            //! *e*
            
        } else if (mode == "dlg") {
            STRUCTURE *mstructure;
            PARSER *mparser;
            mstructure = new STRUCTURE();
            mparser = new PARSER(mstructure,verb,true,true,true,parse_comments);
            int energy_mode = 0;
            if (write_free_energy && write_energy) energy_mode = 3;
            else if (write_free_energy) energy_mode = 2;
            else if (write_energy) energy_mode = 1;
            string pre_cav_name = "cav";
            if (flex_res_dlg) {
                if (pdb_source == "X") {
                    cerr << c_message<cERROR>("in -f mode you have to specify a pdb file by --s2=my.pdb") << endl;
                    delete mstructure;
                    delete mparser;
                    exit(1);
                }
                if (!(mparser->read_pdb(pdb_source.c_str()))) {
                    cerr << c_message<cERROR>("could not read '") << pdb_source.c_str() << "' as pdb file" << endl;
                    delete mstructure;
                    delete mparser;
                    exit(1);
                }
                string namebuf = pdb_source;
                remove_ext(namebuf);
                pre_cav_name = namebuf;
            } else if (full_pdb_from_dlg) {
                if (pdb_source == "X") {
                    cerr << c_message<cERROR>("in -g mode you have to specify a pdb file by --s2=my.pdb") << endl;
                    delete mstructure;
                    delete mparser;
                    exit(1);
                }
                if (!(mparser->read_pdb(pdb_source.c_str()))) {
                    cerr << c_message<cERROR>("could not read '") << pdb_source.c_str() << "' as pdb file" << endl;
                    delete mstructure;
                    delete mparser;
                    exit(1);
                }
                string namebuf = pdb_source;
                remove_ext(namebuf);
                pre_cav_name = namebuf;
                
                if (!(mparser->read_dlg_for_full_pdb((*fit).c_str(),high,low,pre_cav_name.c_str()))) {
                    delete mstructure;
                    delete mparser;
                    continue;
                }
                
                delete mstructure;
                delete mparser;
                instruct->clear();
                continue;
            }
            
            if (!(mparser->read_dlg((*fit).c_str(),high,low,energy_mode,flex_res_dlg,cut_radius,pre_cav_name.c_str()))) {
                delete mstructure;
                delete mparser;
                continue;
            }
            if (mstructure->ligands.size() == 0) {
                cerr << c_message<cERROR>("no molecules found in ") << *fit << endl;
                delete mstructure;
                delete mparser;
                continue;
            }
            if (auto_dlg) {
                for (ligands_vec lv=mstructure->ligands.begin(); lv!=mstructure->ligands.end(); ++lv) {
                    (*lv)->dlg2mol2(atom_typing,def_file.c_str(),true,verb,true,false,"X",max_ring_members,
                                    no_planar_free_rot,0,0,prot_acids,prot_guanidin,prot_amidin,
                                    prot_amin,prot_phosphate,prot_sulfate);
                    (*lv)->rename_atoms(false);
                    string namebuf;
                    namebuf = *fit;
                    remove_ext(namebuf);
                    namebuf = namebuf + (*lv)->name;
                    (*lv)->name = namebuf;
                }
            } else {
                if (sourcefile == "X") {
                    cerr << c_message<cERROR>("you have to specify a reference file by --s=your_ref.mol2") << endl;
                    exit(1);
                }
                STRUCTURE *rstructure;
                PARSER *rparser;
                rstructure = new STRUCTURE();
                rparser = new PARSER(rstructure,verb,true,true,true,parse_comments);
                if (!(rparser->read_mol2(sourcefile.c_str()))) {
                    cerr << c_message<cERROR>("could not read '") << sourcefile << "' as mol2 file" << endl;
                    continue;
                }
                if (rstructure->ligands.size() == 0) {
                    cerr << c_message<cERROR>("found no ligand in '") << sourcefile << endl;
                    delete rstructure;
                    delete rparser;
                    continue;
                }
                string namebuf = rstructure->ligands[0]->name;
                rstructure->ligands[0]->get_bonded_only();
                for (ligands_vec lv=mstructure->ligands.begin(); lv!=mstructure->ligands.end(); ++lv) {
                    (*lv)->dlg2mol2(rstructure->ligands[0]);
                    string nb = namebuf + (*lv)->name;
                    (*lv)->name = nb;
                }
                delete rstructure;
                delete rparser;
            }
            if (dlg2singlemol2) {
                if (high > mstructure->ligands.size()) high = mstructure->ligands.size();
                cout << "writing solutions " << low << "-" << high << " from " << *fit << endl;
                for (ligands_vec lig2=mstructure->ligands.begin(); lig2!=mstructure->ligands.end(); ++lig2) {
                    if (sourcefiles.size() > 1) {
                        tfile = *fit;
                        if (targetfile != "X") {
                            ostringstream os;
                            os << "_" << targetfile << "_" << (*lig2)->name << ".mol2";
                            replace_ext(tfile,".dlg",os.str());
                        } else {
                            replace_ext(tfile,".dlg","_");
                            tfile += (*lig2)->name;
                            tfile += ".mol2";
                        }
                    } else {
                        if (targetfile != "X") {
                            tfile = targetfile + "_" + (*lig2)->name + ".mol2";
                        } else {
                            tfile = *fit;
                            replace_ext(tfile,".dlg","_");
                            tfile += (*lig2)->name;
                            tfile += ".mol2";
                        }
                    }
                    string_fu::remove_char(tfile);
                    inparser->write_mol2_ligand(*lig2,tfile.c_str());
                }
            } else {
                if (sourcefiles.size() > 1) {
                    tfile = *fit;
                    if (targetfile != "X") {
                        ostringstream os;
                        os << "_" << targetfile << ".mol2";
                        replace_ext(tfile,".dlg",os.str());
                    } else replace_ext(tfile,"dlg","mol2");
                } else {
                    if (targetfile != "X") tfile = targetfile;
                    else {
                        tfile = *fit;
                        replace_ext(tfile,"dlg","mol2");
                    }
                }
                string_fu::remove_char(tfile);
                if (high > mstructure->ligands.size()) high = low + mstructure->ligands.size() - 1;
                cout << "writing solutions " << low << "-" << high << " from " << *fit << " to " << tfile << endl;
                mparser->write_mol2(tfile.c_str());
            }
            delete mstructure;
            delete mparser;

        } else if (mode == "cif") {
            unsigned int total_mol_read = 0;
            int block_size = max_mol2_block;
            bool go_on = true;
            must_append = false;
            while (go_on) {
                int n_mol_read = 0;
                n_mol_read = inparser->read_next_cif((*fit).c_str(),block_size);
                
                if (n_mol_read < 1) break;
                total_mol_read += n_mol_read;
            
                if (n_mol_read < block_size) go_on = false;
            
                if (convert) {
                    for (ligands_vec lig=instruct->ligands.begin(); lig!=instruct->ligands.end(); ++lig) {

                        // cerr << (*lig)->name << endl;

                        (*lig)->get_atom_typing(atom_typing,false,def_file.c_str(),true,
                                                verb,true,false,"X",max_ring_members,
                                                true,0,0,prot_acids,prot_guanidin,prot_amidin,prot_amin,prot_phosphate,
                                                prot_sulfate,kekulize_aromatics,kekulize_charged,allow_charged_aromatics);
                    }
                    if (sourcefiles.size() > 1) {
                        tfile = *fit;
                        if (targetfile != "X") {
                            ostringstream os;
                            os << "_" << targetfile << ".mol2";
                            replace_ext(tfile,".cif",os.str());
                        } else replace_ext(tfile,"cif","mol2");
                    } else {
                        if (targetfile != "X") tfile = targetfile;
                        else {tfile = *fit; replace_ext(tfile,"cif","mol2");}
                    }
                    string_fu::remove_char(tfile);
                    if (go_on || must_append) {
                        if (int(instruct->ligands.size()) >= max_mol2_block && tfile == *fit) {
                            cerr << c_message<cERROR>("if there are more than ") << max_mol2_block-1 
                            << " molecules in your mol2 (the current blocksize of fconv), you have to specify a target with '--t'" << endl;
                            exit(0);
                        }
                        inparser->write_next_mol2(instruct->ligands,tfile.c_str());
                    } else inparser->write_mol2(tfile.c_str());
                    if (go_on) must_append = true;
                    
                    instruct->clear();
                    continue;
                }
                
                if (build_crystal) {
                    int lig_num = 0 + total_mol_read - n_mol_read;
                    
                    for (ligands_vec lig=instruct->ligands.begin(); lig!=instruct->ligands.end(); ++lig) {
                        lig_num++;
                        CRYSTALIZER *mycrys;
                        mycrys = new CRYSTALIZER(*lig);
                        
                    //    if (unit_cell_multi > 0) mycrys->build_crystal2(unit_cell_multi);
                        if (unit_cell_multi > 0) mycrys->build_crystal2_splitted(unit_cell_multi);
                        else mycrys->build_crystal_splitted(unit_cell_exp_rad);
                    //    else mycrys->build_crystal(unit_cell_exp_rad);
                        
                        PARSER *pars2;
                        pars2 = new PARSER(mycrys->structure,1,true,true,true,parse_comments);
                    
                        if (instruct->ligands.size() > 1) {
                            tfile = *fit;
                            if (targetfile == "X") {
                                ostringstream os;
                                os << "_" << lig_num << "_crys.mol2";
                                replace_ext(tfile,".mol2",os.str());
                            } else {
                                tfile = targetfile;
                                ostringstream os;
                                os << "_" << lig_num << ".mol2";
                                replace_ext(tfile,".mol2",os.str());
                            }
                        } else {
                            tfile = *fit;
                            if (targetfile == "X") {
                                ostringstream os;
                                os << "_crys.mol2";
                                replace_ext(tfile,".mol2",os.str());
                            } else {
                                tfile = targetfile;
                            }
                        }
                        
                        for (ligands_vec lig=mycrys->structure->splitted_ligands.begin(); 
                                lig!=mycrys->structure->splitted_ligands.end(); ++lig) {
                            (*lig)->get_atom_typing(atom_typing,false,def_file.c_str(),true,0,true,false,"X",
                                                    max_ring_members,true,0,0,prot_acids,prot_guanidin,prot_amidin,
                                                    prot_amin,prot_phosphate,prot_sulfate,kekulize_charged,allow_charged_aromatics);
                        }
                        string_fu::remove_char(tfile);
                        pars2->write_mol2(mycrys->structure->splitted_ligands,tfile.c_str());
                    //    pars2->write_mol2(tfile.c_str());
                        delete pars2;
                        delete mycrys;
                    }
                    
                    instruct->clear();
                    continue;
                }
                
                instruct->clear_ligands();
            }
        }
        if (must_append && (!merge_mol2_in1mol)) inparser->write_next_mol2_end();
        if (!merge_mol2_in1mol) instruct->clear();
    }
    
    if (merge_mol2 || (add_mol2 && !(auto_dlg))) {
        inparser->write_next_mol2_end();
    //    inparser->write_mol2(tfile.c_str());
    } else if (merge_mol2_in1mol) {
        unsigned int offs = 0;
        unsigned int last_id = 0;
        unsigned int bnd_offs = 0;
        unsigned int last_bnd = 0;
        unsigned int res_offs = 0;
        unsigned int last_res = 0;
        for (ligands_vec lig=instruct->ligands.begin(); lig!=instruct->ligands.end(); ++lig) {
            for (atoms_vec at=(*lig)->atoms.begin(); at!=(*lig)->atoms.end(); ++at) {
                (*at)->intern_id += offs;
                (*at)->id += offs;
                (*at)->res_number += res_offs;
                last_id = (*at)->intern_id;
                last_res = (*at)->res_number;
            }
            for (bonds_vec bt=(*lig)->bonds.begin(); bt!=(*lig)->bonds.end(); ++bt) {
                (*bt)->id += bnd_offs;
                last_bnd = (*bt)->id;
            }
            offs = last_id;
            bnd_offs = last_bnd;
            res_offs = last_res;
            
            for (atoms_vec at=(*lig)->atoms.begin(); at!=(*lig)->atoms.end(); ++at) {
                stl_ptr<ATOM> atm(new ATOM(**at));
                m2_lig->atoms.push_back(atm);
            }
            for (bonds_vec bt=(*lig)->bonds.begin(); bt!=(*lig)->bonds.end(); ++bt) {
                stl_ptr<BOND> bnd(new BOND());
                bnd->id = (*bt)->id;
                bnd->free_rot = (*bt)->free_rot;
                bnd->type = (*bt)->type;
                for (atoms_vec at=m2_lig->atoms.begin(); at!=m2_lig->atoms.end(); ++at) {
                    if ((*at)->intern_id == (*bt)->from->intern_id) bnd->from = (*at);
                    else if ((*at)->intern_id == (*bt)->to->intern_id) bnd->to = (*at);
                }
                m2_lig->bonds.push_back(bnd);
            }
        }
        
        m2_lig->name = "merged_mol";
        vector<stl_ptr<LIGAND> > lw;
        lw.push_back(m2_lig);
        string_fu::remove_char(tfile);
        inparser->write_mol2(lw,tfile.c_str());
        m2_lig.kill();
    }
    
    delete instruct;
    delete inparser;
    
    return 0;
}

