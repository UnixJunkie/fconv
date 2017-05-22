
//============================================================================
// atom_properties_GN.cpp -*- C++ -*-; contains type specific atom properties
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


#include"atom_properties_GN.h"

//==================================================================================================
//Atomtypen:
//==================================================================================================

//!Anzahl der internen Atomtypen:
const int n_intern_types = 160;

//!Versionsstring fuer die atom types:
const string a_t_version = "250911_2";

const int n_metal_list = 19;

const string metal_list[n_metal_list] = { // Si, Se and As not considered as metals
    "LI","BE","NA","MG","AL"," K","CA","CR","MN","FE","CO","NI","CU","ZN","AG","SN","CE","PT","AU"
};

const string metal_list_ra[n_metal_list] = { // Si, Se and As not considered as metals
    " LI"," BE"," NA"," MG"," AL","  K"," CA"," CR"," MN"," FE"," CO"," NI"," CU"," ZN"," AG"," SN"," CE"," PT"," AU"
};

//!es folgen die internen Atomtypen (mode 0):
const string i_t[n_intern_types] = {
    "H.ac","H.onh","H.n","H.o","H.0",                                //5
    
    "C.ar6p","C.ar6x","C.ar6","C.arp","C.arx","C.ar","C.2r3o","C.2r3x","C.3r3x","C.2r3","C.3r3",
    "C.1n","C.1p","C.1s","C.co2h","C.co2","C.es","C.hal","C.am","C.o","C.s","C.gu","C.guh",
    "C.mi","C.mih","C.n","C.2p","C.2s","C.2t","C.et","C.ohp","C.ohs","C.oht","C.3n","C.3p",
    "C.3s","C.3t","C.3q",                                        //38
    
    "N.ar6p","N.ar6","N.arp","N.ar2","N.ar3","N.ar3h","N.r3","N.az","N.1","N.o2","N.ohac","N.oh","N.ims",
    "N.imt","N.amp","N.ams","N.amt","N.samp","N.sams","N.samt","N.gu1","N.gu2","N.guh","N.mi1",
    "N.mi2","N.mih","N.aap","N.aas2","N.aas3","N.aat2","N.aat3","N.2n","N.2p","N.2s","N.2t","N.3n","N.3p","N.3s",
    "N.3t","N.4q","N.4h",                                            //41
    
    "O.ar","O.r3","O.h2o","O.n","O.noh","O.2co2","O.2es","O.2hal","O.am","O.co2","O.2po","O.2so",
    "O.2p","O.2s","O.3po","O.3so","O.carb","O.o","O.3ac","O.ph",
    "O.3oh","O.3es","O.3eta","O.3et",                                //24
    
    "S.ar","S.r3","S.thi","S.o","S.o2h","S.o3h","S.o4h","S.o2","S.o3",
    "S.o4","S.2","S.sh","S.s","S.3",                                //14
    
    "P.r3","P.o","P.o2h","P.o3h","P.o4h","P.o2","P.o3","P.o4","P.3",                //9
    
    "F.0","Cl.0","Br.0","I.0","F.i","Cl.i","Br.i","I.i",                        //8
    
    "Li","Na","Mg","Al","Si","K","Ca","Cr.th","Cr.oh","Mn","Fe","Co","Cu","Zn","Se","Mo","Sn",
    "Ni","Hg","B","As"                                            //21
};

//!es folgen die entsprechenden Sybyltypen (mode 1):
const string mode_1[n_intern_types] = {
    "H","H","H","H","H",                                        //5
    
    "C.ar","C.ar","C.ar","C.2","C.2","C.ar","C.2","C.2","C.3","C.2","C.3",
    "C.1","C.1","C.1","C.2","C.2","C.2","C.2","C.2","C.2","C.2","C.2","C.cat",
    "C.2","C.2","C.2","C.2","C.2","C.2","C.3","C.3","C.3","C.3","C.3","C.3",
    "C.3","C.3","C.3",                                        //38
    
    "N.pl3","N.ar","N.pl3","N.2","N.pl3","N.pl3","N.3","N.1","N.1","N.pl3","N.am","N.3","N.am",
    "N.am","N.am","N.am","N.am","N.am","N.am","N.am","N.2","N.pl3","N.pl3","N.2",
    "N.pl3","N.pl3","N.pl3","N.pl3","N.3","N.pl3","N.3","N.2","N.2","N.2","N.pl3","N.3","N.3","N.3",
    "N.3","N.4","N.4",                                            //41
    
    "O.3","O.3","O.3","O.2","O.3","O.2","O.2","O.2","O.2","O.co2","O.2","O.2",
    "O.co2","O.co2","O.3","O.3","O.2","O.3","O.3","O.3",
    "O.3","O.3","O.3","O.3",                                    //24
    
    "S.3","S.3","S.2","S.o","S.o2","S.o2","S.o2","S.o2","S.o2","S.o2","S.2","S.3","S.3","S.3",    //14
    
    "P.3","P.3","P.3","P.3","P.3","P.3","P.3","P.3","P.3",                        //9
    
    "F","Cl","Br","I","F","Cl","Br","I",                                //8
    
    "Li","Na","Mg","Al","Si","K","Ca","Cr.th","Cr.oh","Mn","Fe","Co","Cu","Zn","Se","Mo","Sn",
    "Ni","Hg","B","As"                                            //21
};

//!es folgen die modifizierten Sybyltypen (mode 2):
const string mode_2[n_intern_types] = {
    "H","H","H","H","H",                                        //5
    
    "C.ar","C.ar","C.ar","C.ar","C.ar","C.ar","C.2","C.2","C.3","C.2","C.3",
    "C.1","C.1","C.1","C.2","C.2","C.2","C.2","C.2","C.2","C.2","C.2","C.cat",
    "C.2","C.2","C.2","C.2","C.2","C.2","C.3","C.3","C.3","C.3","C.3","C.3",
    "C.3","C.3","C.3",                                        //38
    
    "N.pl3","N.ar","N.pl3","N.ar","N.pl3","N.pl3","N.3","N.1","N.1","N.pl3","N.am","N.3","N.am",
    "N.am","N.am","N.am","N.am","N.am","N.am","N.am","N.2","N.pl3","N.pl3","N.2",
    "N.pl3","N.pl3","N.pl3","N.pl3","N.3","N.pl3","N.3","N.2","N.2","N.2","N.pl3","N.3","N.3","N.3",
    "N.3","N.4","N.4",                                            //41
    
    "O.2","O.3","O.3","O.2","O.3","O.2","O.2","O.2","O.2","O.co2","O.2","O.2",
    "O.co2","O.co2","O.3","O.3","O.2","O.3","O.3","O.3",
    "O.3","O.3","O.3","O.3",                                    //24
    
    "S.2","S.3","S.2","S.o","S.o2","S.o2","S.o2","S.o2","S.o2","S.o2","S.2","S.3","S.3","S.3",    //14
    
    "P.3","P.3","P.3","P.3","P.3","P.3","P.3","P.3","P.3",                        //9
    
    "F","Cl","Br","I","F","Cl","Br","I",                                //8
    
    "Li","Na","Mg","Al","Si","K","Ca","Cr.th","Cr.oh","Mn","Fe","Co","Cu","Zn","Se","Mo","Sn",
    "Ni","Hg","B","As"                                            //21
};

//!Jetzt die Anzahl der Valenzen fuer jeden Atomtyp:
const int atom_valence[n_intern_types] = {
    1,1,1,1,1,                                            //5
    
    3,3,3,3,3,3,3,3,4,3,4,
    2,2,2,3,3,3,3,3,3,3,3,3,
    3,3,3,3,3,3,4,4,4,4,4,4,
    4,4,4,                                                //38
    
    3,2,3,2,3,3,3,2,1,3,3,3,3,
    3,3,3,3,3,3,3,2,3,3,2,
    3,3,3,3,3,3,3,2,2,2,3,3,3,3,
    3,4,4,                                                //41
    
    2,2,2,1,2,1,1,1,1,1,1,1,
    1,1,2,2,1,2,2,2,
    2,2,2,2,                                            //24
    
    2,2,1,3,4,4,4,4,4,
    4,1,2,2,2,                                            //14
    
    4,4,4,4,4,4,4,4,4,                                        //9
    
    1,1,1,1,0,0,0,0,                                        //8
    
    0,0,0,0,4,0,0,0,0,0,0,0,0,0,4,0,0,
    0,0,3,4                                                //21
};

//!Bindungsgeometrien:
const int atom_geometrie[n_intern_types] = { //!0:nichts / 1:linear / 2:trigonal planar / 3:tetraedrisch / 4:unknown
    0,0,0,0,0,                                            //5
    
    2,2,2,2,2,2,2,2,3,2,3,
    1,1,1,2,2,2,2,2,2,2,2,2,
    2,2,2,2,2,2,3,3,3,3,3,3,
    3,3,3,                                                //38
    
    2,2,2,2,2,2,3,1,1,2,2,3,2,
    2,2,2,2,2,2,2,2,3,2,2,
    3,2,3,2,3,2,3,2,2,2,2,3,3,3,
    3,3,3,                                                //41
    
    2,3,3,2,3,2,2,2,2,2,2,2,
    2,2,3,3,2,3,3,3,
    3,3,3,3,                                            //24
    
    2,3,3,2,3,3,3,3,3,
    3,2,3,3,3,                                            //14
    
    3,3,3,3,3,3,3,3,3,                                        //9
    
    0,0,0,0,0,0,0,0,                                        //8
    
    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
    4,4,4,4                                                //21
};

//!Hybridisierungen:
const int atom_hybridizations[n_intern_types] = {
    0,0,0,0,0,                                            //5
    
    2,2,2,2,2,2,2,2,3,2,3,
    1,1,1,2,2,2,2,2,2,2,2,2,
    2,2,2,2,2,2,3,3,3,3,3,3,
    3,3,3,                                                //38
    
    2,2,2,2,2,2,3,1,1,2,2,3,2,
    2,2,2,2,2,2,2,2,3,2,2,
    3,2,3,2,3,2,3,2,2,2,2,3,3,3,
    3,3,3,                                                //41
    
    2,3,3,2,3,2,2,2,2,2,2,2,
    2,2,3,3,2,3,3,3,
    3,3,3,3,                                            //24
    
    2,3,2,2,2,2,2,2,2,
    2,2,3,3,3,                                            //14
    
    3,2,2,2,2,2,2,2,3,                                        //9
    
    0,0,0,0,0,0,0,0,                                        //8
    
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0                                                //21
};

//! X--H  Bindungslaengen:
const float atom_hb_length[n_intern_types] = { //!Typen die keine Protonen haben koennen haben dummy-werte!!!
    1.08061,1.08061,1.08061,1.08061,1.08061,                            //5
    
    0.956379,0.955857,0.960605,0.956379,0.952796,0.960605,0.968379,0.968379,0.979254,0.968379,0.976759,
    0.952907,0.952907,0.965793,0.991346,1.01374,0.955684,0.955684,0.982986,0.980125,0.978913,0.950231,0.950231,
    0.963721,0.950231,0.966764,0.972083,0.962652,0.962652,0.982888,0.989513,0.991182,0.991182,0.98106,0.975382,
    0.984575,0.985716,0.985716,                                    //38
    
    0.913131,0.913131,0.902918,0.899805,0.899805,0.899805,0.901294,0.897547,0.897547,0.909496,0.881489,0.909496,0.8869,
    0.8869,0.899901,0.886293,0.886293,0.865833,0.864477,0.864477,0.907462,0.887531,0.882072,0.884302,
    0.897229,0.889454,0.894858,0.889982,0.889982,0.889982,0.889982,0.899218,0.892104,0.902016,0.902016,0.902304,0.894173,0.889304,
    0.889304,0.922775,0.922775,                                        //41
    
    0.853137,1.10295,0.876543,0.914642,0.914642,0.92502,0.853137,0.853137,0.853137,0.92502,0.930701,0.949969,
    0.853137,0.752269,1.36144,1.36144,0.853137,0.917998,0.92502,0.892335,
    0.872473,0.92502,1.39311,1.10295,                                //24
    
    1.07966,1.23896,1.23896,1.23896,1.23896,1.23896,1.23896,1.23896,1.23896,
    1.23896,1.23896,1.2022,1.23896,1.23896,                                //14
    
    1.23896,1.23896,1.23896,1.23896,1.23896,1.23896,1.23896,1.23896,1.23896,            //9
    
    0.973445,1.24552,1.4,1.6,0.973445,1.24552,1.4,1.6,                        //8
    
    1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,
    1.,1.,1.,1.                                            //21
};

const bool prot_acids_default = false;
const bool prot_guanidin_default = true;
const bool prot_amidin_default = true;
const bool prot_amin_default = false;
const bool prot_phosphate_default = false;
const bool prot_sulfate_default = false;
const bool get_bonds_default = true;
const bool kekulize_aromatics_default = false;
const bool kekulize_charged_default = false;
const bool allow_charged_aromatics_default = true;
const int max_ring_members_default = 10;

namespace def_map {
    //!und hier die Zuweisungsmap fuer die Atomtypen:
    tr1::unordered_map<string,string> a_t_map;
    
    //!eine zweite map fuer anwendungen in denen immer zwischen 2 maps gewechselt wird:
    tr1::unordered_map<string,string> a_t_map2;
    
    //!ein flag, mit welchem def_file die map bereits befuellt ist:
    string curr_a_t_map = "X"; //X / def_file / m0 / m1 / m2
    
    string curr_a_t_map2 = "X";

    bool curr_prot_acids_default = false;
    bool curr_prot_guanidin_default = true;
    bool curr_prot_amidin_default = true;
    bool curr_prot_amin_default = false;
    bool curr_prot_phosphate_default = false;
    bool curr_prot_sulfate_default = false;
    bool curr_get_bonds_default = true;
    bool curr_kekulize_aromatics_default = false;
    bool curr_kekulize_charged_default = false;
    bool curr_allow_charged_aromatics_default = true;
    int curr_max_ring_members_default = 10;

    bool curr_prot_acids_default2 = false;
    bool curr_prot_guanidin_default2 = true;
    bool curr_prot_amidin_default2 = true;
    bool curr_prot_amin_default2 = false;
    bool curr_prot_phosphate_default2 = false;
    bool curr_prot_sulfate_default2 = false;
    bool curr_get_bonds_default2 = true;
    bool curr_kekulize_aromatics_default2 = false;
    bool curr_kekulize_charged_default2 = false;
    bool curr_allow_charged_aromatics_default2 = true;
    int curr_max_ring_members_default2 = 10;
    
    int last_a_t_map = 0;
}


namespace atom_properties {
    const string normal_aacids[31] = {"ALA","ARG","ASN","ASP","ASH","ASX","CYS","CYX","CYM","GLN","GLU",
                                      "GLH","GLX","GLY","ILE","LEU","LYS","LYN","MET","PHE","PRO","SER",
                                      "THR","TRP","TYR","TYM","VAL","HIS","HID","HIE","HIP"};
    const string modified_aacids[11] = {"ABA","CGU","CME","CSD","MEN","MLY","MSE","PCA","PTR","SEP","TPO"};
    
    const string nucleic_acids[11] = {"  A","  T","  G","  C","  U"," +U"," DA"," DT"," DG"," DC"," DU"};
    
    const string known_elements[29] = {
        "H","C","N","O","S","P","F","Cl","Br","I","Li","Na","Mg","Al","Si","K","Ca","Cr",
        "Mn","Fe","Co","Cu","Zn","Se","Mo","Sn","Ni","Hg","B"
    };
    const int number_of_elements = 25;
    const string vdW_elements[number_of_elements]= {
        "H", "C", "N", "O", "S", "P", "F", "Cl", "Br", "I",
        "Li", "Na", "Mg", "Al", "Si", "K", "Ca", "B", "X", "Zn",
        "Fe", "Cu", "Ni",
        "Co", "Mn"
    };
    const float vdW_radii[number_of_elements]= { //united atom modell (wasserstoffe mit drin)
        1.20, //! H
        1.77, //! C
        1.35, //! N (um 0.2 erniedrigt um moegliche H-Bond zu beruecksichtigen)
        1.32, //! O (um 0.2 erniedrigt um moegliche H-Bond zu beruecksichtigen)
        1.65, //! S (um 0.2 erniedrigt um moegliche H-Bond zu beruecksichtigen)
        1.80, //! P
        1.47, //! F
        1.75, //! Cl
        1.85, //! Br
        1.98, //! I
        0.46, //! Li+   bei Koordinationszahl 6 (vdW - 0.3)
        0.72, //! Na+   bei Koordinationszahl 6 (vdW - 0.3)
        0.42, //! Mg++  bei Koordinationszahl 6 (vdW - 0.3)
        0.24, //! Al+++ bei Koordinationszahl 6 (vdW - 0.3)
        2.10, //! Si
        1.08, //! K+    bei Koordinationszahl 6 (vdW - 0.3)
        0.70, //! Ca++  bei Koordinationszahl 6 (vdW - 0.3)
        0.80, //! B
        1.60, //! X
        0.44, //! Zn++  bei Koordinationszahl 6 (vdW - 0.3)
        0.31, //! Fe++  bei Koordinationszahl 6 (vdW - 0.3)
        0.43, //! Cu++  bei Koordinationszahl 6 (vdW - 0.3)
        0.39, //! Ni++  bei Koordinationszahl 6 (vdW - 0.3)
        0.35, //! Co++  bei Koordinationszahl 6 (vdW - 0.3)
        0.53  //! Mn++  bei Koordinationszahl 6 (vdW - 0.3)
    };
    /*
    const float vdW_radii2[number_of_elements]= {
        //! http://de.wikipedia.org/wiki/Van-der-Waals-Radius
        1.10, //! H
        1.70, //! C
        1.35, //! N (um 0.2 erniedrigt um moegliche H-Bond zu beruecksichtigen)
        1.32, //! O (um 0.2 erniedrigt um moegliche H-Bond zu beruecksichtigen)
        1.60, //! S (um 0.2 erniedrigt um moegliche H-Bond zu beruecksichtigen)
        1.80, //! P
        1.47, //! F
        1.75, //! Cl
        1.85, //! Br
        1.98, //! I
        0.46, //! Li+   bei Koordinationszahl 6 (vdW - 0.3)
        0.72, //! Na+   bei Koordinationszahl 6 (vdW - 0.3)
        0.42, //! Mg++  bei Koordinationszahl 6 (vdW - 0.3)
        0.24, //! Al+++ bei Koordinationszahl 6 (vdW - 0.3)
        2.10, //! Si
        1.08, //! K+    bei Koordinationszahl 6 (vdW - 0.3)
        0.70, //! Ca++  bei Koordinationszahl 6 (vdW - 0.3)
        1.92, //! B
        1.60, //! X
        0.44, //! Zn++  bei Koordinationszahl 6 (vdW - 0.3)
        0.31, //! Fe++  bei Koordinationszahl 6 (vdW - 0.3)
        0.43, //! Cu++  bei Koordinationszahl 6 (vdW - 0.3)
        0.39, //! Ni++  bei Koordinationszahl 6 (vdW - 0.3)
        0.35, //! Co++  bei Koordinationszahl 6 (vdW - 0.3)
        0.53  //! Mn++  bei Koordinationszahl 6 (vdW - 0.3)
    };
    */
    const float clash_radii[number_of_elements]= { //Atomradien, bzw. Ionenradien
        //! http://www.uniterra.de/rutherford/
        0.373, //! H
        0.772, //! C
        0.710, //! N
        0.604, //! O
        1.040, //! S
        0.930, //! P
        0.709, //! F
        0.994, //! Cl
        1.145, //! Br
        1.331, //! I
        0.780, //! Li+
        0.980, //! Na+
        0.780, //! Mg++
        0.570, //! Al+++
        1.170, //! Si
        1.330, //! K+
        1.060, //! Ca++
        0.830, //! B
        0.600, //! X
        0.830, //! Zn++
        0.670, //! Fe+++
        0.720, //! Cu++
        0.620, //! Ni+++
        0.640, //! Co++
        0.520  //! Mn+++
    };
    const float mweights[number_of_elements] = {
        1.008,  //! H
        12.011, //! C
        14.007, //! N
        15.999, //! O
        32.070, //! S
        30.974, //! P
        18.998, //! F
        35.453, //! Cl
        79.900, //! Br
        126.90, //! I
        6.941,  //! Li
        22.990, //! Na
        24.305, //! Mg
        26.982, //! Al
        28.086, //! Si
        39.100, //! K
        40.080, //! Ca
        10.810, //! B
        0.0,    //! X
        65.38,  //! Zn
        55.85,  //! Fe
        63.55,  //! Cu
        58.70,  //! Ni
        58.93,  //! Co
        54.93   //! Mn
    };
    
    tr1::unordered_map<string,float> vdW_map; //!vdW-Radien
    tr1::unordered_map<string,float> square_vdW_map;
    tr1::unordered_map<string,float> clash_map; //!Atom(/Ionen)-Radien
    tr1::unordered_map<string,float> square_clash_map;
    tr1::unordered_map<string,float> mol_weight;
    
    tr1::unordered_set<string> known_aacids;
    tr1::unordered_set<string> known_mod_aacids;
    tr1::unordered_set<string> known_nucleic_acids;
    
    void initialize() {
        if (vdW_map.size() > 0) return;
        for (int i=0; i<number_of_elements; ++i) {
            vdW_map[vdW_elements[i]] = vdW_radii[i];
            square_vdW_map[vdW_elements[i]] = vdW_radii[i] * vdW_radii[i];
            clash_map[vdW_elements[i]] = clash_radii[i];
            square_clash_map[vdW_elements[i]] = clash_radii[i] * clash_radii[i];
            mol_weight[vdW_elements[i]] = mweights[i];
        }
        for (int i=0; i<31; ++i) known_aacids.insert(normal_aacids[i]);
        for (int i=0; i<11; ++i) known_mod_aacids.insert(modified_aacids[i]);
        for (int i=0; i<11; ++i) known_nucleic_acids.insert(nucleic_acids[i]);
    }
}

