
//============================================================================
// crystallizer_GN.cpp -*- C++ -*-; building crystal packings
//
// Copyright (C) 2007, 2008, 2009, 2010 Gerd Neudert
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
// This library provides routines for building crystal packings from
// information contained in CIF files or in CRYSIN entries of MOL2 files.
//============================================================================

#include"crystallizer_GN.h"


//==============================================================================================
//Definitionen fr CRYSTALIZER:
//==============================================================================================

CRYSTALIZER::CRYSTALIZER(LIGAND &ligand,bool vis,bool verb) {
    structure = new STRUCTURE(); //bekommt ein neues STRUCTURE-Object
    stl_ptr<LIGAND> lig(new LIGAND(ligand)); //mit einer echten Kopie des "zu kristallisierenden" Liganden
    lig->main_structure = structure;
    structure->ligands.push_back(lig);
    l = lig;
    visualize = vis;
    verbosity = verb;
    sym_ele = 0;
    cart_x = vec3d<float>(1.,0.,0.);
    cart_y = vec3d<float>(0.,1.,0.);
    cart_z = vec3d<float>(0.,0.,1.);
    if (!(ligand.crysin.zero())) {
        cart_x *= l->crysin->cryst2cart;
        cart_y *= l->crysin->cryst2cart;
        cart_z *= l->crysin->cryst2cart;
    }
}

CRYSTALIZER::CRYSTALIZER(stl_ptr<LIGAND> &ligand,bool vis,bool verb) {
    structure = new STRUCTURE(); //bekommt ein neues STRUCTURE-Object
    stl_ptr<LIGAND> lig(new LIGAND(*ligand)); //mit einer echten Kopie des "zu kristallisierenden" Liganden
    lig->main_structure = structure;
    structure->ligands.push_back(lig);
    l = lig;
    visualize = vis;
    verbosity = verb;
    sym_ele = 0;
    cart_x = vec3d<float>(1.,0.,0.);
    cart_y = vec3d<float>(0.,1.,0.);
    cart_z = vec3d<float>(0.,0.,1.);
    if (!(ligand->crysin.zero())) {
        cart_x *= l->crysin->cryst2cart;
        cart_y *= l->crysin->cryst2cart;
        cart_z *= l->crysin->cryst2cart;
    }
}

CRYSTALIZER::~CRYSTALIZER() {
    delete structure;
}

//!PROBLEM:
//! 1.) Ein Molekuel, welches auf einer Flaeche der Einheitszelle liegt und durch eine Operation auf der
//!     gegenueberliegenden Flaeche zu liegen kommt, wird unter Umstaenden nicht als identisch erkannt!!!
//!     Im Prinzip muesste hier bei der Vervielfaeltigung nochmal die Ueberpruefung auf Duplizitaet kommen!!!
//!     (z.B. bei QIHCIT)
//! 2.) Diese Ueberpruefung hier macht das Kristallbauen extrem langsam fuer Molekuele mit grosser Anzahl
//!     an Atomen. => Was besseres ausdenken!!!
bool CRYSTALIZER::is_copy(stl_ptr<LIGAND> &lig) {
    //!Bei symmetrischen Molekuelen gehoert meist nur ein Teil des Molekuels zur asymmetrischen Einheit
    //!=> Wenn z.B. nur das halbe Molekuel zur asym. Einheit gehoert ergeben 8 Symmetrieoperationen nur
    //!   4 Molekuele
    //!Da die Operationen auf ganze Molekuele angewendet werden, muessen die doppelten hier rausgefiltert werden
    //vec3d<float> buf;
    
    vec3d<float> ref[27];
    
    ref[0] = lig->atoms[0]->coord;
    
    ref[1] = lig->atoms[0]->coord + cart_x; ref[2] = lig->atoms[0]->coord - cart_x;
    ref[3] = lig->atoms[0]->coord + cart_y; ref[4] = lig->atoms[0]->coord - cart_y;
    ref[5] = lig->atoms[0]->coord + cart_z; ref[6] = lig->atoms[0]->coord - cart_z;
    
    ref[7] = ref[1] + cart_y; ref[8] = ref[1] - cart_y;
    ref[9] = ref[1] + cart_z; ref[10] = ref[1] - cart_z;
    ref[11] = ref[2] + cart_y; ref[12] = ref[2] - cart_y;
    ref[13] = ref[2] + cart_z; ref[14] = ref[2] - cart_z;
    
    ref[15] = ref[3] + cart_z; ref[16] = ref[3] - cart_z;
    ref[17] = ref[4] + cart_z; ref[18] = ref[4] - cart_z;
    
    ref[19] = ref[7] + cart_z; ref[20] = ref[11] + cart_z;
    ref[21] = ref[7] - cart_z; ref[22] = ref[11] - cart_z;
    
    ref[23] = ref[8] + cart_z; ref[24] = ref[12] + cart_z;
    ref[24] = ref[8] - cart_z; ref[26] = ref[12] - cart_z;
    
    //cerr << "KOPIE:";
    
    for (ligands_vec lt=structure->ligands.begin(); lt!=structure->ligands.end(); ++lt) {
        //wenn eines der Atome von lig mit einem Atom eines anderen Liganden uebereinstimmt wird davon
        //ausgegangen, dass es das gleiche Molekuel (spiegelsymmetrisch) ist
        
        //for (atoms_vec at=lig->atoms.begin(); at!=lig->atoms.end(); ++at) {
            for (atoms_vec sat=(*lt)->atoms.begin(); sat!=(*lt)->atoms.end(); ++sat) {
                //if (get_square_distance((*sat)->coord,(*at)->coord) < 0.025) return true;
                //if (get_square_distance((*sat)->coord,lig->atoms[0]->coord) < 0.025) return true;
                for (int i=0; i<27; ++i) {
                    if (get_square_distance((*sat)->coord,ref[i]) < 0.025) {/*cerr << " JA" << endl;*/ return true;}
                }
            }
        //}
    }
    
    return false;
}

void CRYSTALIZER::correct_setting(vec3d<float> &v) {
    if (l->crysin->group < 3 || l->crysin->setting > 142) {
        v *= l->crysin->cryst2cart;
        return;
    }
    if (l->crysin->group < 16) { //Monokline Systeme
        float buf;
        switch (l->crysin->setting) {
            case 1: v *= l->crysin->cryst2cart; return;
            case 2: buf = v[0]; v[0] = v[2]; v[2] = buf; v *= l->crysin->cryst2cart; break;
            case 3: buf = v[1]; v[1] = v[2]; v[2] = -1.*buf; v *= l->crysin->cryst2cart; break;
            case 4: buf = v[0]; v[0] = -1.*v[1]; v[1] = v[2]; v[2] = buf; v *= l->crysin->cryst2cart; break;
        }
    } else if (l->crysin->group < 75) { //Orthorhombische Systeme
        float buf;
        switch (l->crysin->setting) {
            case 1: v *= l->crysin->cryst2cart; return;
            case 2: buf = v[1]; v[1] = v[0]; v[0] = v[2]; v[2] = buf; v *= l->crysin->cryst2cart; break;
            case 3: buf = v[0]; v[0] = v[1]; v[1] = v[2]; v[2] = buf; v *= l->crysin->cryst2cart; break;
            case 4: buf = v[1]; v[1] = -1.*v[2]; v[2] = buf; v *= l->crysin->cryst2cart; break;
            case 5: buf = v[0]; v[0] = v[1]; v[1] = buf; v[2] = -1.*v[2]; v *= l->crysin->cryst2cart; break;
            case 6: buf = v[0]; v[0] = -1.*v[2]; v[2] = buf; v *= l->crysin->cryst2cart; break;
        }
    } else { //Tetragonale Systeme
        float buf;
        switch (l->crysin->setting) {
            case 1: v *= l->crysin->cryst2cart; return;
            case 2: buf = v[0]; v[0] += v[1]; v[1] -= buf; v *= l->crysin->cryst2cart; break;
        }
    }
}

/*
void CRYSTALIZER::change_origin_choice(vec3d<float> &v,vec3d<float> &sym1,vec3d<float> &sym2) {
    //!Bei einigen Raumgruppen kann ein unterschiedlicher Ursprung fuer die Einheitszelle
    //!gewaehlt werden. Bei Origin Choice 2 ist der Ursprung um (0.25,0.25,0.25) verschoben.
    //!=> Die Symmetrieelemente muessen entsprechend der neuen relativen Lage verschoben werden.
    //!Die Vektoren sym1 und sym2 sollen den Symmetrieelementen entsprechen (sym 2 fuer glides,
    //!sonst mit sym1 identisch). Komponenten, die in sym1 oder sym2 von Null verschieden sind
    //!werden NICHT geaendert. !ACHTUNG: Das ist falsch, siehe z.B. RG 224 !!!
    
    //!Das vostehende ist so nicht richtig!!!
    //!Ich muss noch die allgemeine Regel finden, nach der die Lage der Symmetrieelemente
    //!geaendert wird.
    
    if (l->crysin->group == 48 || l->crysin->group == 50 || l->crysin->group == 59 ||
        l->crysin->group == 68 || l->crysin->group == 70 || l->crysin->group == 85 ||
        l->crysin->group == 86 || l->crysin->group == 88 || l->crysin->group == 125 ||
        l->crysin->group == 126 || l->crysin->group == 129 || l->crysin->group == 130 ||
        l->crysin->group == 133 || l->crysin->group == 134 || l->crysin->group == 137 ||
        l->crysin->group == 138 || l->crysin->group == 141 || l->crysin->group == 142 ||
        l->crysin->group == 201 || l->crysin->group == 203 || l->crysin->group == 222 ||
        l->crysin->group == 224 || l->crysin->group == 227 || l->crysin->group == 228) {
        
        //!Hier muss ein check hin, ob es sich um Origin Choice 2 handelt
        //!Wenn nicht, dann return!!!
        
        v *= l->crysin->cart2cryst;
        if (sym1[0] == 0. && sym2[0] == 0.) v[0] -= 0.25;
        if (sym1[1] == 0. && sym2[1] == 0.) v[1] -= 0.25;
        if (sym1[2] == 0. && sym2[2] == 0.) v[2] -= 0.25;
        v *= l->crysin->cryst2cart;
    }
    
}
*/

bool CRYSTALIZER::origin_one() {
    //! Soll True liefern, wenn Origin Choice 1 vorliegt und False, wenn
    //! OC2 vorliegt.
    n_mol = structure->ligands.size();
    for (int i=0; i<n_mol; ++i) {
        for (int j=i+1; j<n_mol; ++j) {
            for (atoms_vec at=structure->ligands[i]->atoms.begin(); at!=structure->ligands[i]->atoms.end(); ++at) {
                for (atoms_vec at2=structure->ligands[j]->atoms.begin(); at2!=structure->ligands[j]->atoms.end(); ++at2) {
                    if ((*at)->sybyl_type[0] != 'C' && (*at)->sybyl_type[0] != 'N' &&
                        (*at)->sybyl_type[0] != 'O' && (*at)->sybyl_type[0] != 'S' && (*at)->sybyl_type[0] != 'H') continue;
                    if ((*at2)->sybyl_type[0] != 'C' && (*at2)->sybyl_type[0] != 'N' &&
                        (*at2)->sybyl_type[0] != 'O' && (*at2)->sybyl_type[0] != 'S' && (*at2)->sybyl_type[0] != 'H') continue;
                    if (get_square_distance((*at)->coord,(*at2)->coord) < 1.5) return false;
                }
            }
        }
    }
    return true;
}

void CRYSTALIZER::kill_oc1() {
    n_mol = structure->ligands.size();
    for (int i=1; i<n_mol; ++i) {
        structure->ligands[i].kill();
        --n_ligs;
    }
    structure->ligands.clear();
    structure->ligands.push_back(l);
    sym_ele = 0;
    sym_name.clear();
    
    if (visualize) {
        f_out.close();
        f_out.clear();
        f_out.open("sym_vis.py");
        f_out << "# This visualization file was automatically created by 'structure_GN.cpp'" << "\n";
        f_out << "# Start pymol and use 'run this_files_name.py' to start the visualization" << "\n" << "\n";
        
        vec3d<float> ucx(1.,0.,0.); vec3d<float> ucy(0.,1.,0.); vec3d<float> ucz(0.,0.,1.);
        vec3d<float> ucxy(1.,1.,0.); vec3d<float> ucxz(1.,0.,1.); vec3d<float> ucxyz(1.,1.,1.); vec3d<float> ucyz(0.,1.,1.);
        
        ucx *= l->crysin->cryst2cart;
        ucy *= l->crysin->cryst2cart;
        ucz *= l->crysin->cryst2cart;
        ucxy *= l->crysin->cryst2cart;
        ucxz *= l->crysin->cryst2cart;
        ucxyz *= l->crysin->cryst2cart;
        ucyz *= l->crysin->cryst2cart;
        
        for (int ih=0; ih<3; ++ih) {
            if (fabs(ucx[ih]) < 0.0001) ucx[ih] = 0.; //!weil pymol die nomenklatur '1.234e-11' nicht kann!!!
            if (fabs(ucy[ih]) < 0.0001) ucy[ih] = 0.;
            if (fabs(ucz[ih]) < 0.0001) ucz[ih] = 0.;
            if (fabs(ucxy[ih]) < 0.0001) ucxy[ih] = 0.;
            if (fabs(ucxz[ih]) < 0.0001) ucxz[ih] = 0.;
            if (fabs(ucxyz[ih]) < 0.0001) ucxyz[ih] = 0.;
            if (fabs(ucyz[ih]) < 0.0001) ucyz[ih] = 0.;
        }
        
        //!namen und anhaengen !!!!!
        ostringstream name;
        name << "unit_cell";
        sym_name.push_back(name.str());
        
        f_out << name.str() << "=[9.0," << "0.,0.,0.,"
              << ucx[0] << "," << ucx[1] << "," << ucx[2] << ",0.04,1.,0.,0.,1.,0.,0.,9.0,"
              << "0.,0.,0.,"
              << ucy[0] << "," << ucy[1] << "," << ucy[2] << ",0.04,0.,1.,0.,0.,1.,0.,9.0,"
              << "0.,0.,0.,"
              << ucz[0] << "," << ucz[1] << "," << ucz[2] << ",0.04,0.,0.,1.,0.,0.,1.,9.0,"
              
              << ucy[0] << "," << ucy[1] << "," << ucy[2] << ","
              << ucxy[0] << "," << ucxy[1] << "," << ucxy[2] << ",0.02,0.6,0.6,0.6,0.6,0.6,0.6,9.0,"
              << ucx[0] << "," << ucx[1] << "," << ucx[2] << ","
              << ucxy[0] << "," << ucxy[1] << "," << ucxy[2] << ",0.02,0.6,0.6,0.6,0.6,0.6,0.6,9.0,"
              << ucx[0] << "," << ucx[1] << "," << ucx[2] << ","
              << ucxz[0] << "," << ucxz[1] << "," << ucxz[2] << ",0.02,0.6,0.6,0.6,0.6,0.6,0.6,9.0,"
              << ucz[0] << "," << ucz[1] << "," << ucz[2] << ","
              << ucxz[0] << "," << ucxz[1] << "," << ucxz[2] << ",0.02,0.6,0.6,0.6,0.6,0.6,0.6,9.0,"
              << ucy[0] << "," << ucy[1] << "," << ucy[2] << ","
              << ucyz[0] << "," << ucyz[1] << "," << ucyz[2] << ",0.02,0.6,0.6,0.6,0.6,0.6,0.6,9.0,"
              << ucyz[0] << "," << ucyz[1] << "," << ucyz[2] << ","
              << ucxyz[0] << "," << ucxyz[1] << "," << ucxyz[2] << ",0.02,0.6,0.6,0.6,0.6,0.6,0.6,9.0,"
              << ucxz[0] << "," << ucxz[1] << "," << ucxz[2] << ","
              << ucxyz[0] << "," << ucxyz[1] << "," << ucxyz[2] << ",0.02,0.6,0.6,0.6,0.6,0.6,0.6,9.0,"
              << ucz[0] << "," << ucz[1] << "," << ucz[2] << ","
              << ucyz[0] << "," << ucyz[1] << "," << ucyz[2] << ",0.02,0.6,0.6,0.6,0.6,0.6,0.6,9.0,"
              << ucxy[0] << "," << ucxy[1] << "," << ucxy[2] << ","
              << ucxyz[0] << "," << ucxyz[1] << "," << ucxyz[2] << ",0.02,0.6,0.6,0.6,0.6,0.6,0.6]" << endl;
    }
}

void CRYSTALIZER::correct_pos(stl_ptr<LIGAND> &lig) {
    vec3d<float> buf;
    buf = lig->atoms[0]->coord;
    buf *= l->crysin->cart2cryst;
    for (int i=0; i<3; ++i) {
        if (buf[i] > 1.) {
            for (atoms_vec at=lig->atoms.begin(); at!=lig->atoms.end(); ++at) {
                (*at)->coord *= l->crysin->cart2cryst;
                (*at)->coord[i] -= 1.;
                (*at)->coord *= l->crysin->cryst2cart;
            }
        } else if (buf[i] < 0.) {
            for (atoms_vec at=lig->atoms.begin(); at!=lig->atoms.end(); ++at) {
                (*at)->coord *= l->crysin->cart2cryst;
                (*at)->coord[i] += 1.;
                (*at)->coord *= l->crysin->cryst2cart;
            }
        }
    }
}

//!Die Vektoren koennen nicht per reference uebergeben werden, weil sie eventuell
//!durch correct_setting() veraendert werden!!!
void CRYSTALIZER::gl(vec3d<float> v1,vec3d<float> v2,vec3d<float> posi,vec3d<float> direction) { //glide
    
    correct_setting(v1); correct_setting(v2); correct_setting(posi); correct_setting(direction);
    
    //!TEST:
//    change_origin_choice(posi,v1,v2);
    
    //Spiegelmatrix:  A = E - 2 * n*n^T   (Achtung: n*n^T => n*n-Matrix   /   die Normalenvectoren muessen normiert sein)
    vec3d<float> norm = vectorproduct(v1,v2);
    norm.norm(); //normierter Normalenvector der Spiegelebene
    matrix<float> m(1.,0.,0.,
                    0.,1.,0.,
                    0.,0.,1.);
    m[0][0] -= 2. * norm[0] * norm[0]; m[0][1] -= 2. * norm[0] * norm[1]; m[0][2] -= 2. * norm[0] * norm[2];
    m[1][0] -= 2. * norm[1] * norm[0]; m[1][1] -= 2. * norm[1] * norm[1]; m[1][2] -= 2. * norm[1] * norm[2];
    m[2][0] -= 2. * norm[2] * norm[0]; m[2][1] -= 2. * norm[2] * norm[1]; m[2][2] -= 2. * norm[2] * norm[2];
    stl_ptr<LIGAND> lig(new LIGAND(*l)); lig->main_structure = structure;
    ostringstream hos; hos << n_ligs; lig->name += hos.str();
    
    //!neu: wenn die Verschiebung nicht orthogonal zur Spiegelebene ist, muss der Verschiebungsvektor neu
    //!berechnet werden:
    vec3d<float> rv(norm);
    float pbuf = posi.value();
    if (pbuf > 0.) rv *= pbuf * cos(angle(norm,posi));
    else rv = posi;
    rv *= 2.;
    //vec3d<float> addv(posi); addv *= 2.;
    
    for (atoms_vec at=lig->atoms.begin(); at!=lig->atoms.end(); ++at) {
    //    (*at)->coord *= l->crysin->cart2cryst;
        (*at)->coord *= m; //Spiegelung an Ursprungsebene
    
    //    (*at)->coord += addv; //Verschiebung fr nicht im Ursprung liegende Spiegelebene
        (*at)->coord += rv; //Verschiebung fr nicht im Ursprung liegende Spiegelebene
        
        (*at)->coord += direction; //Translationsteil
    //    (*at)->coord *= l->crysin->cryst2cart;
    }
    
    if (visualize) {
        sym_ele++;
        ostringstream name;
        if (direction.value() < 0.0001) name << "mirror_" << sym_ele;
        else name << "glide_" << sym_ele;
        sym_name.push_back(name.str());
        vec3d<float> dv1(l->atoms[0]->coord);
        vec3d<float> dv3(lig->atoms[0]->coord);
    //    dv1 *= l->crysin->cart2cryst;
        vec3d<float> dv2(dv1);
        dv2 *= m;
        
    //    dv2 += addv;
        dv2 += rv;
        
        vec3d<float> ev1(v1); vec3d<float> ev2(v2);
        vec3d<float> ev3(ev1); vec3d<float> ev4(ev2);
        ev3 *= -1.; ev4 *= -1.;
        ev1 += posi; ev2 += posi; ev3 += posi; ev4 += posi;
    /*    
        dv1 *= l->crysin->cryst2cart;
        dv2 *= l->crysin->cryst2cart;
        ev1 *= l->crysin->cryst2cart;
        ev2 *= l->crysin->cryst2cart;
        ev3 *= l->crysin->cryst2cart;
        ev4 *= l->crysin->cryst2cart;
    */    
        f_out << name.str() << "=[9.0," << dv1[0] << "," << dv1[1] << "," << dv1[2] << ","
              << dv2[0] << "," << dv2[1] << "," << dv2[2] << ",0.03,0.6,0.6,0.6,0.6,0.6,0.6,9.0," 
              << ev1[0] << "," << ev1[1] << "," << ev1[2] << ","
              << ev2[0] << "," << ev2[1] << "," << ev2[2] << ",0.03,0.8,0.8,0.8,0.8,0.8,0.8,9.0,"
              << ev1[0] << "," << ev1[1] << "," << ev1[2] << ","
              << ev3[0] << "," << ev3[1] << "," << ev3[2] << ",0.03,0.8,0.8,0.8,0.8,0.8,0.8,9.0,"
              << ev1[0] << "," << ev1[1] << "," << ev1[2] << ","
              << ev4[0] << "," << ev4[1] << "," << ev4[2] << ",0.03,0.8,0.8,0.8,0.8,0.8,0.8,9.0,"
              << ev2[0] << "," << ev2[1] << "," << ev2[2] << ","
              << ev3[0] << "," << ev3[1] << "," << ev3[2] << ",0.03,0.8,0.8,0.8,0.8,0.8,0.8,9.0,"
              << ev2[0] << "," << ev2[1] << "," << ev2[2] << ","
              << ev4[0] << "," << ev4[1] << "," << ev4[2] << ",0.03,0.8,0.8,0.8,0.8,0.8,0.8,9.0,"
              << ev3[0] << "," << ev3[1] << "," << ev3[2] << ","
              << ev4[0] << "," << ev4[1] << "," << ev4[2] << ",0.03,0.8,0.8,0.8,0.8,0.8,0.8,9.0,"
              << dv2[0] << "," << dv2[1] << "," << dv2[2] << ","
              << dv3[0] << "," << dv3[1] << "," << dv3[2] << ",0.03,0.6,0.6,0.6,0.6,0.6,0.6]" << endl;
    }
    correct_pos(lig);
    if (!is_copy(lig)) {structure->ligands.push_back(lig); ++n_ligs;}
    else lig.kill();
}

void CRYSTALIZER::sc(vec3d<float> axis,float angle,vec3d<float> posi,vec3d<float> direction) { //!achsen in fractional space
    
    correct_setting(axis); correct_setting(posi); correct_setting(direction);
    
    //!TEST:
//    change_origin_choice(posi,axis,axis);
    
    vec3d<float> bax(axis);
    bax.norm();
    matrix<float> m = rotmatrix(bax,angle);
    stl_ptr<LIGAND> lig(new LIGAND(*l)); lig->main_structure = structure;
    ostringstream hos; hos << n_ligs; lig->name += hos.str();
    
    
    //!PROBLEM: Die Berechnung der Rotationsmatrix beruht auf der Annahme eines rechtwinkligen Koordinatensystems
    //!=> nicht die Atomkoordinaten in Kristallkoordinaten umrechnen, sondern stattdessen die Symmetrieelemente 
    //!   in den kartesischen Raum ueberfuehren!
    
    for (atoms_vec at=lig->atoms.begin(); at!=lig->atoms.end(); ++at) {
    //    (*at)->coord *= l->crysin->cart2cryst;
        (*at)->coord -= posi;
        (*at)->coord *= m;
        (*at)->coord += posi;
        (*at)->coord += direction;
    //    (*at)->coord *= l->crysin->cryst2cart;
    }
    
    if (visualize) {
        sym_ele++;
        int hdw = int((angle * 360. / 6.2831853) + 0.5);
        if (hdw > 180) hdw -= 360;
        ostringstream name;
        if (direction.value() < 0.0001) name << "rot";
        else name << "screw";
        if (hdw < 0) {name << "I"; hdw *= -1;}
        name << hdw << "_" << sym_ele;
        sym_name.push_back(name.str());
        vec3d<float> v1(0.,0.,0.), v2(axis);
        vec3d<float> a(axis);
        if (a[0] < 0.) a[0] *= -1.;
        if (a[1] < 0.) a[1] *= -1.;
        if (a[2] < 0.) a[2] *= -1.;
        v1 += posi; v2 += posi;
        
        vec3d<float> hv1(l->atoms[0]->coord);
    //    hv1 *= l->crysin->cart2cryst;
        vec3d<float> dv1 = hv1 - v1;
        vec3d<float> dv2(dv1);
        dv2 *= m;
        dv2 += direction;
        dv1 += posi;
        dv2 += posi;
        vec3d<float> dv02 = v1 + direction;
    //    dv1 *= l->crysin->cryst2cart;
    //    dv2 *= l->crysin->cryst2cart;
    //    dv02 *= l->crysin->cryst2cart;
        
    //    v1 *= l->crysin->cryst2cart;
    //    v2 *= l->crysin->cryst2cart;
        f_out << name.str() << "=[9.0," << v1[0] << "," << v1[1] << "," << v1[2] << ","
              << v2[0] << "," << v2[1] << "," << v2[2] << ",0.07," << a[0] << "," << a[1] << "," << a[2] << "," 
              << a[0] << "," << a[1] << "," << a[2] 
              << ",9.0," << v1[0] << "," << v1[1] << "," << v1[2] << ","
              << dv1[0] << "," << dv1[1] << "," << dv1[2] << ",0.03," << a[0]*0.8 << "," << a[1]*0.8 << "," << a[2]*0.8 << "," 
              << a[0]*0.8 << "," << a[1]*0.8 << "," << a[2]*0.8
              << ",9.0," << dv02[0] << "," << dv02[1] << "," << dv02[2] << ","
              << dv2[0] << "," << dv2[1] << "," << dv2[2] << ",0.03," << a[0]*0.8 << "," << a[1]*0.8 << "," << a[2]*0.8 << "," 
              << a[0]*0.8 << "," << a[1]*0.8 << "," << a[2]*0.8
              << ",9.0," << dv1[0] << "," << dv1[1] << "," << dv1[2] << ","
              << dv2[0] << "," << dv2[1] << "," << dv2[2] << ",0.03," << a[0]*0.8 << "," << a[1]*0.8 << "," << a[2]*0.8 << "," 
              << a[0]*0.8 << "," << a[1]*0.8 << "," << a[2]*0.8 << "]" << endl;
    }
    
    correct_pos(lig);
    if (!is_copy(lig)) {structure->ligands.push_back(lig); ++n_ligs;}
    else lig.kill();
}

void CRYSTALIZER::sci(vec3d<float> axis,float angle,vec3d<float> posi,vec3d<float> direction,vec3d<float> point) { //!Rotoinversion
    
    correct_setting(axis); correct_setting(posi); correct_setting(direction); correct_setting(point);
    
    vec3d<float> bax(axis);
    bax.norm();
    matrix<float> m = rotmatrix(bax,angle);
    stl_ptr<LIGAND> lig(new LIGAND(*l)); lig->main_structure = structure;
    ostringstream hos; hos << n_ligs; lig->name += hos.str();
    for (atoms_vec at=lig->atoms.begin(); at!=lig->atoms.end(); ++at) {
    //    (*at)->coord *= l->crysin->cart2cryst;
        (*at)->coord -= posi;
        (*at)->coord *= m;
        (*at)->coord += posi;
        (*at)->coord += direction; //!Achtung: direction darf nur die Verschiebung durch Screw beinhalten (nicht Inversion)
        (*at)->coord -= point;
        (*at)->coord *= -1.;
        (*at)->coord += point;
    //    (*at)->coord *= l->crysin->cryst2cart;
    }
    
    if (visualize) {
        sym_ele++;
        int hdw = int((angle * 360. / 6.2831853) + 0.5);
        if (hdw > 180) hdw -= 360;
        ostringstream name;
        name << "rotoinversion";
        if (hdw < 0) {name << "I"; hdw *= -1;}
        name << hdw << "_" << sym_ele;
        sym_name.push_back(name.str());
        vec3d<float> v1(0.,0.,0.), v2(axis);
        vec3d<float> a(axis);
        if (a[0] < 0.) a[0] *= -1.;
        if (a[1] < 0.) a[1] *= -1.;
        if (a[2] < 0.) a[2] *= -1.;
        v1 += posi; v2 += posi;
        
        vec3d<float> hv1(l->atoms[0]->coord);
    //    hv1 *= l->crysin->cart2cryst;
        vec3d<float> dv1 = hv1 - v1;
        vec3d<float> dv2(dv1);
        dv2 *= m;
        dv2 += direction;
        dv1 += posi;
        dv2 += posi;
        vec3d<float> dv02 = v1 + direction;
        
        vec3d<float> dv3(point);
        dv3 -= dv2;
        dv3 += point;
    /*    
        dv1 *= l->crysin->cryst2cart;
        dv2 *= l->crysin->cryst2cart;
        dv02 *= l->crysin->cryst2cart;
        dv3 *= l->crysin->cryst2cart;
    */    
        vec3d<float> hpoint(point);
    //    hpoint *= l->crysin->cryst2cart;
        
    //    v1 *= l->crysin->cryst2cart;
    //    v2 *= l->crysin->cryst2cart;
        
        f_out << name.str() << "=[9.0," << v1[0] << "," << v1[1] << "," << v1[2] << ","
              << v2[0] << "," << v2[1] << "," << v2[2] << ",0.07," << a[0] << "," << a[1] << "," << a[2] << "," 
              << a[0] << "," << a[1] << "," << a[2] << ",7.0," << hpoint[0] << "," 
              << hpoint[1] << "," << hpoint[2] << ",0.4"
              << ",9.0," << v1[0] << "," << v1[1] << "," << v1[2] << ","
              << dv1[0] << "," << dv1[1] << "," << dv1[2] << ",0.03," << a[0]*0.8 << "," << a[1]*0.8 << "," << a[2]*0.8 << "," 
              << a[0]*0.8 << "," << a[1]*0.8 << "," << a[2]*0.8
              << ",9.0," << dv02[0] << "," << dv02[1] << "," << dv02[2] << ","
              << dv2[0] << "," << dv2[1] << "," << dv2[2] << ",0.03," << a[0]*0.8 << "," << a[1]*0.8 << "," << a[2]*0.8 << "," 
              << a[0]*0.8 << "," << a[1]*0.8 << "," << a[2]*0.8
              << ",9.0," << dv1[0] << "," << dv1[1] << "," << dv1[2] << ","
              << dv2[0] << "," << dv2[1] << "," << dv2[2] << ",0.03," << a[0]*0.8 << "," << a[1]*0.8 << "," << a[2]*0.8 << "," 
              << a[0]*0.8 << "," << a[1]*0.8 << "," << a[2]*0.8
              << ",9.0," << dv2[0] << "," << dv2[1] << "," << dv2[2] << ","
              << dv3[0] << "," << dv3[1] << "," << dv3[2] << ",0.02," << a[0]*0.8 << "," << a[1]*0.8 << "," << a[2]*0.8 << "," 
              << a[0]*0.8 << "," << a[1]*0.8 << "," << a[2]*0.8
              <<"]" << endl;
    }
    correct_pos(lig);
    if (!is_copy(lig)) {structure->ligands.push_back(lig); ++n_ligs;}
    else lig.kill();
}

void CRYSTALIZER::in(vec3d<float> point) { //!point ist coordinate im fractional space
    
    correct_setting(point);
    
    //!TEST:
//    vec3d<float> v0(0.,0.,0.);
//    change_origin_choice(point,v0,v0);
    
    stl_ptr<LIGAND> lig(new LIGAND(*l)); lig->main_structure = structure;
    ostringstream hos; hos << n_ligs; lig->name += hos.str();
    for (atoms_vec at=lig->atoms.begin(); at!=lig->atoms.end(); ++at) {
    //    (*at)->coord *= l->crysin->cart2cryst;
        (*at)->coord -= point;
        (*at)->coord *= -1.;
        (*at)->coord += point;
    //    (*at)->coord *= l->crysin->cryst2cart;
    }
    
    if (visualize) {
        sym_ele++;
        ostringstream name;
        name << "inversion_" << sym_ele;
        sym_name.push_back(name.str());
        vec3d<float> dv1(l->atoms[0]->coord);
    //    dv1 *= l->crysin->cart2cryst;
        vec3d<float> dv3(point);
        dv3 -= dv1;
        dv3 += point;
    //    dv1 *= l->crysin->cryst2cart;
    //    dv3 *= l->crysin->cryst2cart;
        vec3d<float> hpoint(point);
    //    hpoint *= l->crysin->cryst2cart;
        //!--------------------------
        f_out << name.str() << "=[9.0," << dv1[0] << "," << dv1[1] << "," << dv1[2] << ","
              << dv3[0] << "," << dv3[1] << "," << dv3[2] << ",0.03,0.8,0.8,0.8,0.8,0.8,0.8,7.0," << hpoint[0] << "," 
              << hpoint[1] << "," << hpoint[2] << ",0.4]" << endl;
    }
    correct_pos(lig);
    if (!is_copy(lig)) {structure->ligands.push_back(lig); ++n_ligs;}
    else lig.kill();
}

void CRYSTALIZER::tr(vec3d<float> trans) {
    
    correct_setting(trans);
    
    stl_ptr<LIGAND> lig(new LIGAND(*l)); lig->main_structure = structure;
    ostringstream hos; hos << n_ligs; lig->name += hos.str();
    for (atoms_vec at=lig->atoms.begin(); at!=lig->atoms.end(); ++at) {
    //    (*at)->coord *= l->crysin->cart2cryst;
        (*at)->coord += trans;
    //    (*at)->coord *= l->crysin->cryst2cart;
    }
    
    if (visualize) {
        sym_ele++;
        ostringstream name;
        name << "translation_" << sym_ele;
        sym_name.push_back(name.str());
        vec3d<float> dv1(l->atoms[0]->coord);
        vec3d<float> dv2(lig->atoms[0]->coord);
        f_out << name.str() << "=[9.0," << dv1[0] << "," << dv1[1] << "," << dv1[2] << ","
              << dv2[0] << "," << dv2[1] << "," << dv2[2] << ",0.03,0.8,0.8,0.8,0.8,0.8,0.8]" << endl;
    }
    correct_pos(lig);
    if (!is_copy(lig)) {structure->ligands.push_back(lig); ++n_ligs;}
    else lig.kill();
}

/*
void CRYSTALIZER::vis_asym(float xb, float yb, float zb) {
    ofstream f_out;
    f_out.open("cell_axes.py");

    f_out << "# This visualization file was automatically created by 'structure_GN.cpp'" << "\n";
    f_out << "# Start pymol and use 'run this_files_name.py' to start the visualization" << "\n" << "\n";
    
    f_out << "cell_cgo = [";
    
    vec3d<float> orig(0.,0.,0.);
    vec3d<float> x(xb,0.,0.), y(0.,yb,0.), z(0.,0.,zb);
    x *= l->crysin->cryst2cart;
    y *= l->crysin->cryst2cart;
    z *= l->crysin->cryst2cart;
    
    float rad = 0.05;
    
    f_out << "9.0," << orig[0] << "," << orig[1] << "," << orig[2] << ","
                    << x[0] << "," << x[1] << "," << x[2] << "," << rad << ",1.0,0.0,0.0,1.0,0.0,0.0,"
          << "9.0," << orig[0] << "," << orig[1] << "," << orig[2] << ","
                    << y[0] << "," << y[1] << "," << y[2] << "," << rad << ",0.0,1.0,0.0,0.0,1.0,0.0,"
          << "9.0," << orig[0] << "," << orig[1] << "," << orig[2] << ","
                    << z[0] << "," << z[1] << "," << z[2] << "," << rad << ",0.0,0.0,1.0,0.0,0.0,1.0]" << "\n";
    
    f_out << "cmd.load_cgo(cell_cgo, 'asym', 1)" << endl;
    f_out.close();
}
*/

void CRYSTALIZER::t(stl_ptr<CRYSIN_POSITION> &p) {
    
    //cerr << "calling with: " << p->const_add << " " << p->x_add << " " << p->y_add << " " << p->z_add << endl;
    
    stl_ptr<LIGAND> lig(new LIGAND(*l)); lig->main_structure = structure;
    ostringstream hos; hos << n_ligs; lig->name += hos.str();
    
    for (atoms_vec at=lig->atoms.begin(); at!=lig->atoms.end(); ++at) {
        (*at)->coord *= l->crysin->cart2cryst;
        float x = (*at)->coord[0]; float y = (*at)->coord[1]; float z = (*at)->coord[2];
        (*at)->coord[0] = p->const_add[0] + p->x_add[0] * x + p->y_add[0] * y + p->z_add[0] * z;
        (*at)->coord[1] = p->const_add[1] + p->x_add[1] * x + p->y_add[1] * y + p->z_add[1] * z;
        (*at)->coord[2] = p->const_add[2] + p->x_add[2] * x + p->y_add[2] * y + p->z_add[2] * z;
        (*at)->coord *= l->crysin->cryst2cart;
    }
    
    correct_pos(lig);
    if (!is_copy(lig)) {structure->ligands.push_back(lig); ++n_ligs;}
    else lig.kill();
}

bool CRYSTALIZER::build_unit_cell_from_positions() {
    for (vector<stl_ptr<CRYSIN_POSITION> >::iterator it=l->crysin->positions.begin(); it!=l->crysin->positions.end(); ++it) {
        t(*it);
    }
    return true;
}

bool CRYSTALIZER::build_unit_cell() {
    
    if (l->crysin->positions.size() > 0) return build_unit_cell_from_positions();
    
    if (visualize) {
        f_out.open("sym_vis.py");
        f_out << "# This visualization file was automatically created by 'structure_GN.cpp'" << "\n";
        f_out << "# Start pymol and use 'run this_files_name.py' to start the visualization" << "\n" << "\n";
        
        vec3d<float> ucx(1.,0.,0.); vec3d<float> ucy(0.,1.,0.); vec3d<float> ucz(0.,0.,1.);
        vec3d<float> ucxy(1.,1.,0.); vec3d<float> ucxz(1.,0.,1.); vec3d<float> ucxyz(1.,1.,1.); vec3d<float> ucyz(0.,1.,1.);
        
        ucx *= l->crysin->cryst2cart;
        ucy *= l->crysin->cryst2cart;
        ucz *= l->crysin->cryst2cart;
        ucxy *= l->crysin->cryst2cart;
        ucxz *= l->crysin->cryst2cart;
        ucxyz *= l->crysin->cryst2cart;
        ucyz *= l->crysin->cryst2cart;
        
        for (int ih=0; ih<3; ++ih) {
            if (fabs(ucx[ih]) < 0.0001) ucx[ih] = 0.; //!weil pymol die nomenklatur '1.234e-11' nicht kann!!!
            if (fabs(ucy[ih]) < 0.0001) ucy[ih] = 0.;
            if (fabs(ucz[ih]) < 0.0001) ucz[ih] = 0.;
            if (fabs(ucxy[ih]) < 0.0001) ucxy[ih] = 0.;
            if (fabs(ucxz[ih]) < 0.0001) ucxz[ih] = 0.;
            if (fabs(ucxyz[ih]) < 0.0001) ucxyz[ih] = 0.;
            if (fabs(ucyz[ih]) < 0.0001) ucyz[ih] = 0.;
        }
        
        //!namen und anh�ngen !!!!!
        ostringstream name;
        name << "unit_cell";
        sym_name.push_back(name.str());
        
        f_out << name.str() << "=[9.0," << "0.,0.,0.,"
              << ucx[0] << "," << ucx[1] << "," << ucx[2] << ",0.04,1.,0.,0.,1.,0.,0.,9.0,"
              << "0.,0.,0.,"
              << ucy[0] << "," << ucy[1] << "," << ucy[2] << ",0.04,0.,1.,0.,0.,1.,0.,9.0,"
              << "0.,0.,0.,"
              << ucz[0] << "," << ucz[1] << "," << ucz[2] << ",0.04,0.,0.,1.,0.,0.,1.,9.0,"
              
              << ucy[0] << "," << ucy[1] << "," << ucy[2] << ","
              << ucxy[0] << "," << ucxy[1] << "," << ucxy[2] << ",0.02,0.6,0.6,0.6,0.6,0.6,0.6,9.0,"
              << ucx[0] << "," << ucx[1] << "," << ucx[2] << ","
              << ucxy[0] << "," << ucxy[1] << "," << ucxy[2] << ",0.02,0.6,0.6,0.6,0.6,0.6,0.6,9.0,"
              << ucx[0] << "," << ucx[1] << "," << ucx[2] << ","
              << ucxz[0] << "," << ucxz[1] << "," << ucxz[2] << ",0.02,0.6,0.6,0.6,0.6,0.6,0.6,9.0,"
              << ucz[0] << "," << ucz[1] << "," << ucz[2] << ","
              << ucxz[0] << "," << ucxz[1] << "," << ucxz[2] << ",0.02,0.6,0.6,0.6,0.6,0.6,0.6,9.0,"
              << ucy[0] << "," << ucy[1] << "," << ucy[2] << ","
              << ucyz[0] << "," << ucyz[1] << "," << ucyz[2] << ",0.02,0.6,0.6,0.6,0.6,0.6,0.6,9.0,"
              << ucyz[0] << "," << ucyz[1] << "," << ucyz[2] << ","
              << ucxyz[0] << "," << ucxyz[1] << "," << ucxyz[2] << ",0.02,0.6,0.6,0.6,0.6,0.6,0.6,9.0,"
              << ucxz[0] << "," << ucxz[1] << "," << ucxz[2] << ","
              << ucxyz[0] << "," << ucxyz[1] << "," << ucxyz[2] << ",0.02,0.6,0.6,0.6,0.6,0.6,0.6,9.0,"
              << ucz[0] << "," << ucz[1] << "," << ucz[2] << ","
              << ucyz[0] << "," << ucyz[1] << "," << ucyz[2] << ",0.02,0.6,0.6,0.6,0.6,0.6,0.6,9.0,"
              << ucxy[0] << "," << ucxy[1] << "," << ucxy[2] << ","
              << ucxyz[0] << "," << ucxyz[1] << "," << ucxyz[2] << ",0.02,0.6,0.6,0.6,0.6,0.6,0.6]" << endl;
    }
    
    
    //!Momentan alles nur fuer setting == 1 (origin choice 1)!!!!!
    float pi53 = 5.2359878; //300
    float pi32 = 4.7123890; //270
    float pi43 = 4.1887902; //240
    float pi = 3.1415927;  //180
    float pi23 = 2.0943951;//120
    float pi2 = 1.5707963; //90
    float pi3 = 1.0471976; //60
    //float pi4 = 0.7853981; //45 
    
    vec3d<float> v0(0.,0.,0.); vec3d<float> v1(1.,1.,1.);

    vec3d<float> x(1.,0.,0.); vec3d<float> y(0.,1.,0.); vec3d<float> z(0.,0.,1.);
    vec3d<float> x34(0.75,0.,0.); vec3d<float> y34(0.,0.75,0.); vec3d<float> z34(0.,0.,0.75);
    vec3d<float> x23(0.66666666,0.,0.); vec3d<float> y23(0.,0.66666666,0.); vec3d<float> z23(0.,0.,0.66666666);
    vec3d<float> x2(0.5,0.,0.); vec3d<float> y2(0.,0.5,0.); vec3d<float> z2(0.,0.,0.5);
    vec3d<float> x38(0.375,0.,0.); vec3d<float> y38(0.,0.375,0.); vec3d<float> z38(0.,0.,0.375);
    vec3d<float> x3(0.33333333,0.,0.); vec3d<float> y3(0.,0.333333333,0.); vec3d<float> z3(0.,0.,0.33333333);
    vec3d<float> x4(0.25,0.,0.); vec3d<float> y4(0.,0.25,0.); vec3d<float> z4(0.,0.,0.25);
    vec3d<float> x8(0.125,0.,0.); vec3d<float> y8(0.,0.125,0.); vec3d<float> z8(0.,0.,0.125);
    
    vec3d<float> xy(1.,1.,0.); vec3d<float> yz(0.,1.,1.); vec3d<float> xz(1.,0.,1.);
    vec3d<float> xiy(1.,-1.,0.); vec3d<float> yiz(0.,1.,-1.); vec3d<float> xiz(1.,0.,-1.);
    vec3d<float> ixy(-1.,1.,0.); vec3d<float> iyz(0.,-1.,1.); vec3d<float> ixz(-1.,0.,1.);
    vec3d<float> x34y34(0.75,0.75,0.); vec3d<float> y34z34(0.,0.75,0.75); vec3d<float> x34z34(0.75,0.,0.75);
    vec3d<float> x2y2(0.5,0.5,0.); vec3d<float> y2z2(0.,0.5,0.5); vec3d<float> x2z2(0.5,0.,0.5);
    vec3d<float> x4y4(0.25,0.25,0.); vec3d<float> y4z4(0.,0.25,0.25); vec3d<float> x4z4(0.25,0.,0.25);
    vec3d<float> x8y8(0.125,0.125,0.); vec3d<float> y8z8(0.,0.125,0.125); vec3d<float> x8z8(0.125,0.,0.125);
    
    vec3d<float> x2y2z2(0.5,0.5,0.5); vec3d<float> x4y4z4(0.25,0.25,0.25); vec3d<float> x8y8z8(0.125,0.125,0.125); 
    
    vec3d<float> x2yz(0.5,1.,1.); vec3d<float> xy2z(1.,0.5,1.); vec3d<float> xyz2(1.,1.,0.5);
    
    vec3d<float> x2y4(0.5,0.25,0.); vec3d<float> y2z4(0.,0.5,0.25); vec3d<float> x2z4(0.5,0.,0.25);
    vec3d<float> x4y2(0.25,0.5,0.); vec3d<float> y4z2(0.,0.25,0.5); vec3d<float> x4z2(0.25,0.,0.5);
    vec3d<float> x4y8(0.25,0.125,0.); vec3d<float> y4z8(0.,0.25,0.125); vec3d<float> x4z8(0.25,0.,0.125);
    vec3d<float> x8y4(0.125,0.25,0.); vec3d<float> y8z4(0.,0.125,0.25); vec3d<float> x8z4(0.125,0.,0.25);
    vec3d<float> x4y38(0.25,0.375,0.); vec3d<float> y4z38(0.,0.25,0.375); vec3d<float> x4z38(0.25,0.,0.375);
    vec3d<float> x38y4(0.375,0.25,0.); vec3d<float> y38z4(0.,0.375,0.25); vec3d<float> x38z4(0.375,0.,0.25);
    vec3d<float> x4y34(0.25,0.75,0.); vec3d<float> x4z34(0.25,0.,0.75);
    vec3d<float> x34y4(0.75,0.25,0.); vec3d<float> y4z34(0.,0.25,0.75);
    vec3d<float> x34z4(0.75,0.,0.25); vec3d<float> y34z4(0.,0.75,0.25);
    
    vec3d<float> ix4y4(-0.25,0.25,0.); vec3d<float> x4iy4(0.25,-0.25,0.);
    vec3d<float> ix4z4(-0.25,0.,0.25); vec3d<float> x4iz4(0.25,0.,-0.25);
    vec3d<float> iy4z4(0.,-0.25,0.25); vec3d<float> y4iz4(0.,0.25,-0.25);
    
    vec3d<float> ix4y4z4(-0.25,0.25,0.25); vec3d<float> x4y4z34(0.25,0.25,0.75);
    vec3d<float> x4iy4z34(0.25,-0.25,0.75); vec3d<float> ix4(-0.25,0.,0.);
    
    vec3d<float> x8y38z38(0.125,0.375,0.375); vec3d<float> x38y8z38(0.375,0.125,0.375); vec3d<float> x38y38z8(0.375,0.375,0.125);
    
    vec3d<float> x3y3(0.33333333,0.33333333,0.); vec3d<float> y3z3(0.,0.33333333,0.33333333); vec3d<float> x3z3(0.33333333,0.,0.33333333);
    vec3d<float> x23y3(0.66666666,0.33333333,0.); vec3d<float> y23z3(0.,0.66666666,0.33333333);
    vec3d<float> x3y23(0.33333333,0.66666666,0.); vec3d<float> y3z23(0.,0.33333333,0.66666666);
    vec3d<float> x23z3(0.66666666,0.,0.33333333); vec3d<float> x3z23(0.33333333,0.,0.66666666);
    vec3d<float> x3y3z3(0.33333333,0.33333333,0.33333333); vec3d<float> x23y3z3(0.66666666,0.33333333,0.33333333);
    vec3d<float> x3y23z3(0.33333333,0.66666666,0.33333333); vec3d<float> x3y3z23(0.33333333,0.33333333,0.66666666);
    
    vec3d<float> x23y23z3(0.66666666,0.66666666,0.33333333); vec3d<float> x3y23z23(0.33333333,0.66666666,0.66666666);
    vec3d<float> x23y3z23(0.66666666,0.33333333,0.66666666); vec3d<float> x23y23z23(0.66666666,0.66666666,0.66666666);
    
    vec3d<float> ix3y3(-0.33333333,0.33333333,0.); vec3d<float> iy3z3(0.,-0.33333333,0.33333333); vec3d<float> ix3z3(-0.33333333,0.,0.33333333);
    vec3d<float> x3iy3(0.33333333,-0.33333333,0.); vec3d<float> y3iz3(0.,0.33333333,-0.33333333); vec3d<float> x3iz3(0.33333333,0.,-0.33333333);
    
    vec3d<float> dxy(2.,1.,0.); vec3d<float> dxz(2.,0.,1.);
    vec3d<float> xdy(1.,2.,0.); vec3d<float> xdz(1.,0.,2.);
    vec3d<float> dyz(0.,2.,1.); vec3d<float> ydz(0.,1.,2.);
    
    vec3d<float> ix4y4z34(-0.25,0.25,0.75); vec3d<float> x4iy4z4(0.25,-0.25,0.25); vec3d<float> x3y6z6(0.33333333,0.16666666,0.16666666);
    vec3d<float> x3iy3z6(0.33333333,-0.33333333,0.16666666); vec3d<float> x3y23z6(0.33333333,0.66666666,0.16666666);
    vec3d<float> x6y3z3(0.16666666,0.33333333,0.33333333); vec3d<float> ix3y3z3(-0.33333333,0.33333333,0.33333333);
    vec3d<float> z6(0.,0.,0.16666666); vec3d<float> iy6z6(0.,-0.16666666,0.16666666); vec3d<float> y6z6(0.,0.16666666,0.16666666);
    vec3d<float> x3z6(0.33333333,0.,0.16666666); vec3d<float> y6z3(0.,0.16666666,0.33333333); vec3d<float> x6z3(0.16666666,0.,0.33333333);
    vec3d<float> x6iy6z3(0.16666666,-0.16666666,0.33333333); vec3d<float> iy2(0.,-0.5,0.); vec3d<float> ix2(-0.5,0.,0.);
    vec3d<float> ix6y6z23(-0.16666666,0.16666666,0.66666666); vec3d<float> x3y6z23(0.33333333,0.16666666,0.66666666);
    vec3d<float> x6iy6z56(0.16666666,-0.16666666,0.83333333); vec3d<float> x6y3z56(0.16666666,0.33333333,0.83333333);
    vec3d<float> x23y3z56(0.66666666,0.33333333,0.83333333); vec3d<float> ix6y6z6(-0.16666666,0.16666666,0.16666666);
    vec3d<float> iy6z512(0.,-0.16666666,0.41666666); vec3d<float> y6z512(0.,0.16666666,0.41666666);
    vec3d<float> x3z512(0.33333333,0.,0.41666666); vec3d<float> y6z12(0.,0.16666666,0.08333333);
    vec3d<float> y3z12(0.,0.33333333,0.08333333); vec3d<float> x6z12(0.16666666,0.,0.08333333); vec3d<float> z56(0.,0.,0.83333333);
    vec3d<float> z512(0.,0.,0.41666666); vec3d<float> z12(0.,0.,0.08333333);
    
    vec3d<float> ix4iy4iz4(-0.25,-0.25,-0.25); vec3d<float> ix4iy4(-0.25,-0.25,0.);
    vec3d<float> x38y38(0.375,0.375,0.); vec3d<float> x38z38(0.375,0.,0.375); vec3d<float> y38z38(0.,0.375,0.375);
    vec3d<float> x38y8(0.375,0.125,0.); vec3d<float> x38z8(0.375,0.,0.125); vec3d<float> x8y38(0.125,0.375,0.);
    vec3d<float> y38z8(0.,0.375,0.125); vec3d<float> x8z38(0.125,0.,0.375); vec3d<float> y8z38(0.,0.125,0.375);
    
    vec3d<float> x2y4z38(0.5,0.25,0.375); vec3d<float> ix4y2(-0.25,0.5,0.); vec3d<float> x2iy4(0.5,-0.25,0.);
    vec3d<float> x2iy4z8(0.5,-0.25,0.125); vec3d<float> y34z38(0.,0.75,0.375); vec3d<float> iy4z8(0.,-0.25,0.125);
    vec3d<float> x2iy4z38(0.5,-0.25,0.375); vec3d<float> y34z8(0.,0.75,0.125); vec3d<float> x2y4z8(0.5,0.25,0.125);
    vec3d<float> x34y34z4(0.75,0.75,0.25); vec3d<float> iy4z38(0.,-0.25,0.375); vec3d<float> x34y34z34(0.75,0.75,0.75);
    
    //!Ab hier Vectoren fuer die Vertieferraumgruppen:
    vec3d<float> xyz(1.,1.,1.); vec3d<float> ixyiz(-1.,1.,-1.); vec3d<float> xiyiz(1.,-1.,-1.); vec3d<float> ixiy(-1.,-1.,0.);
    vec3d<float> ixiyz(-1.,-1.,1.); vec3d<float> ix3iy6(-0.33333333,-0.16666666,0.); vec3d<float> x3iy6(0.33333333,-0.16666666,0.);
    vec3d<float> ix6y6(-0.16666666,0.16666666,0.); vec3d<float> x6y6(0.16666666,0.16666666,0.); vec3d<float> ix2y2(-0.5,0.5,0.);
    vec3d<float> x6iy6(0.16666666,-0.16666666,0.); vec3d<float> x3iy3z3(0.33333333,-0.33333333,0.33333333);
    vec3d<float> x23iy3(0.66666666,-0.33333333,0.); vec3d<float> x2iy2(0.5,-0.5,0.); vec3d<float> ix6iy3(-0.16666666,-0.33333333,0.);
    vec3d<float> ix6y3(-0.16666666,0.33333333,0.); vec3d<float> x6y3(0.16666666,0.33333333,0.);
    vec3d<float> x3y3iz3(0.33333333,0.33333333,-0.33333333); vec3d<float> x3y6(0.33333333,0.16666666,0.);
    vec3d<float> x6iy6iz6(0.16666666,-0.16666666,-0.16666666); vec3d<float> x6iy6z6(0.16666666,-0.16666666,0.16666666);
    vec3d<float> x6y6iz6(0.16666666,0.16666666,-0.16666666); vec3d<float> ix6iy6z6(-0.16666666,-0.16666666,0.16666666);
    vec3d<float> ix6y6iz6(-0.16666666,0.16666666,-0.16666666); vec3d<float> ix3y23(-0.33333333,0.66666666,0.);
    vec3d<float> x4y34iz4(0.25,0.75,-0.25); vec3d<float> ix4y34z4(-0.25,0.75,0.25); vec3d<float> x34iy4z4(0.75,-0.25,0.25);
    vec3d<float> x34y4iz4(0.75,0.25,-0.25); vec3d<float> xy2(1.,0.5,0.); vec3d<float> ix2iy2(-0.5,-0.5,0.);
    vec3d<float> iy2z2(0.,-0.5,0.5); vec3d<float> x2y(0.5,1.,0.); vec3d<float> ix2y(-0.5,1.,0.); vec3d<float> ixy2(-1.,0.5,0.);
    vec3d<float> y2iz2(0.,0.5,-0.5); vec3d<float> xiy2(1.,-0.5,0.); vec3d<float> x2iz2(0.5,0.,-0.5);
    vec3d<float> ix2z2(-0.5,0.,0.5); vec3d<float> x2iy(0.5,-1.,0.);
    vec3d<float> ix2y4(-0.5,0.25,0.); vec3d<float> x4y(0.25,1.,0.); vec3d<float> x2y34iz4(0.5,0.75,-0.25);
    vec3d<float> x34iy4(0.75,-0.25,0.); vec3d<float> x54iy(1.25,-1.,0.); vec3d<float> x2iy4z34(0.5,-0.25,0.75);
    vec3d<float> ixy4(-1.,0.25,0.); vec3d<float> x34iz4(0.75,0.,-0.25); vec3d<float> x6y512(0.16666666,0.41666666,0.);
    vec3d<float> x34iy2(0.75,-0.5,0.); vec3d<float> ix6y712(-0.16666666,0.58333333,0.); vec3d<float> ixy54(-1.,1.25,0.);
    vec3d<float> ix4y2z34(-0.25,0.5,0.75); vec3d<float> ix4y(-0.25,1.,0.); vec3d<float> ix8y8z38(-0.125,0.125,0.375);
    vec3d<float> x8y38iz8(0.125,0.375,-0.125); vec3d<float> x38iy8z8(0.375,-0.125,0.125); vec3d<float> x8iy8z38(0.125,-0.125,0.375);
    vec3d<float> ix8y38z8(-0.125,0.375,0.125); vec3d<float> x38y8iz8(0.375,0.125,-0.125); vec3d<float> x8y58z8(0.125,0.625,0.125);
    vec3d<float> ix32y(-1.5,1.,0.); vec3d<float> ix58y8z78(0.625,0.125,0.875); vec3d<float> x8y78iz8(0.125,0.875,-0.125);
    vec3d<float> x32y2(1.5,0.5,0.); vec3d<float> x78iy8z58(0.875,-0.125,0.625); vec3d<float> x8y8z58(0.125,0.125,0.625);
    vec3d<float> x8iy8z78(0.125,-0.125,0.875); vec3d<float> ix58y78z8(-0.625,0.875,0.125); vec3d<float> x78y58iz8(0.875,0.625,-0.125);
    vec3d<float> x58y8z8(0.625,0.125,0.125); vec3d<float> x78iy58z8(0.875,-0.625,0.125); vec3d<float> xiy32(1.,-1.5,0.);
    vec3d<float> xy4(1.,0.25,0.); vec3d<float> ix4y34(-0.25,0.75,0.); vec3d<float> x34y2iz4(0.75,0.5,-0.25);
    vec3d<float> ixy34(-1.,0.75,0.); vec3d<float> ix4iy4z4(-0.25,-0.25,0.25); vec3d<float> ix4y4iz4(-0.25,0.25,-0.25);
    vec3d<float> x4iy4iz4(0.25,-0.25,-0.25); vec3d<float> x8y34(0.125,0.75,0.); vec3d<float> iy4z2(0.,-0.25,0.5);
    vec3d<float> x38y34(0.375,0.75,0.); vec3d<float> ix4y8(-0.25,0.125,0.); vec3d<float> x34y38(0.75,0.375,0.);
    vec3d<float> x34y8(0.75,0.125,0.); vec3d<float> ix4z2(-0.25,0.,0.5); vec3d<float> x38y2(0.375,0.5,0.);
    vec3d<float> x38iy4(0.375,-0.25,0.); vec3d<float> x4y54(0.25,1.25,0.); vec3d<float> y2z38(0.,0.5,0.375);
    vec3d<float> iy4z34(0.,-0.25,0.75); vec3d<float> x54y4(1.25,0.25,0.); vec3d<float> x712y512(0.58333333,0.41666666,0.);
    vec3d<float> x512y712(0.41666666,0.58333333,0.); vec3d<float> x34iz34(0.75,0.,-0.75); vec3d<float> iy34z34(0.,-0.75,0.75);
    vec3d<float> x2y38(0.5,0.375,0.); vec3d<float> x2iz4(0.5,0.,-0.25); vec3d<float> y2iz4(0.,0.5,-0.25);
    vec3d<float> x4y4iz4(0.25,0.25,-0.25); vec3d<float> x6y6z6(0.16666666,0.16666666,0.16666666); 
    vec3d<float> ix4y34z2(-0.25,0.75,0.5); vec3d<float> x34iy(0.75,-1.,0.); vec3d<float> x34iy4z2(0.75,-0.25,0.5);
    vec3d<float> ix4z34(-0.25,0.,0.75); vec3d<float> x4iy2(0.25,0.,-0.5); vec3d<float> y34iz4(0.,0.75,-0.25);
    vec3d<float> ix2y34(-0.5,0.75,0.); vec3d<float> x512y6(0.41666666,0.16666666,0.); vec3d<float> x712iy6(0.58333333,-0.16666666,0.);
    vec3d<float> x32iy34(1.5,-0.75,0.); vec3d<float> x34iy32(0.75,-1.25,0.); vec3d<float> x34iy34(0.75,-0.75,0.);
    vec3d<float> x32iy2(1.5,-0.5,0.); vec3d<float> x78y8iz58(0.875,0.125,-0.625); vec3d<float> ixy32(-1.,1.5,0.);
    vec3d<float> ix34y34(-0.75,0.75,0.); vec3d<float> y34iz34(0.,0.75,-0.75); vec3d<float> ix34y32(-0.75,1.5,0.);
    vec3d<float> ix34z34(-0.75,0.,0.75); vec3d<float> ix32y34(-1.5,0.75,0.); vec3d<float> ix8y78z8(-0.125,0.875,0.125);
    vec3d<float> x58iy8z78(0.625,-0.125,0.875); vec3d<float> x8y78iz58(0.125,0.875,-0.625); vec3d<float> x32iy(1.5,-1.,0.);
    vec3d<float> x78iy8z8(0.875,-0.125,0.125); vec3d<float> ix8y58z78(-0.125,0.625,0.875); vec3d<float> x2y32(0.5,1.5,0.);
    vec3d<float> ix8y8z78(-0.125,0.125,0.875); vec3d<float> x58y78iz8(0.625,0.875,-0.125); vec3d<float> x8iy58z78(0.125,-0.625,0.875);
    vec3d<float> ix8y78z58(-0.125,0.875,0.625); vec3d<float> x78y8iz8(0.875,0.125,-0.125); vec3d<float> ix2y32(-0.5,1.5,0.);
    
    bool ret_false = false;
    
    switch (l->crysin->group) {
        //!Rotation wird als Spezialfall von Screw behandelt und Mirrorplane als Spezialfall von Glide
        //!Folgende Synthaxregeln gelten fuer die Symmetrieoperationen:
        //! Inversion     :  in(inversionspunkt)
        //! Translation   :  tr(verschiebung)
        //! Screw         :  sc(drehachse,drehwinkel,verschiebung_der_drehachse,verschiebung)
        //! Rotoinversion :  sci(drehachse,drehwinkel,verschiebung_der_drehachse,verschiebung,inversionspunkt)
        //! Glide         :  gl(achse1,achse2,verschiebung_der_ebene,verschiebung)
        //!
        //!ACHTUNG!: Fuer a-, b- oder c-Glideebenen muss die entsprechende Verschiebung (x2, y2 oder z2) mitgegeben werden!!!
        
        //!--------------------------------------------------------------------------------------------------------
        //!Monoklin:
        case 1: break;
        
        case 2: in(v0); break; //P1_inv    Inversion
        
        case 3: sc(y,pi,v0,v0); break; //P2    180 Rotation um y 
        
        case 4: sc(y,pi,v0,y2); break; //P2_1    180 Screw um y mit Verschiebung 1/2
        
        case 5: sc(y,pi,v0,v0); tr(x2y2); sc(y,pi,x4,y2); break; //C2    jeweils Screw fr Origin-point und C-Flaechenpunkt
        
        case 6: gl(x,z,v0,v0); break; //Pm    Spiegelung an xz-Ebene
        
        case 7: gl(x,z,v0,z2); break; //Pc    Glide-Spiegelung an xz-Ebene mit Translation Richtung z
        
        case 8: gl(x,z,v0,v0); tr(x2y2); gl(x,z,y4,x2); break; //Cm
        
        case 9: gl(x,z,v0,z2); tr(x2y2); gl(x,z,y4,x2z2); break; //Cc
        
        case 10: sc(y,pi,v0,v0); in(v0); gl(x,z,v0,v0); break; //P2/m
        
        case 11: sc(y,pi,v0,y2); in(v0); gl(x,z,y4,v0); break; //P2_1/m
        
        case 12: sc(y,pi,v0,v0); in(v0); gl(x,z,v0,v0); tr(x2y2); sc(y,pi,x4,y2); in(x4y4); gl(x,z,y4,x2); break; //C2/m
        
        case 13: sc(y,pi,z4,v0); in(v0); gl(x,z,v0,z2); break; //P2/c
        
        case 14: sc(y,pi,z4,y2); in(v0); gl(x,z,y4,z2); break; //P2_1/c
        //case 14: t(-1.,0.5,1.,0.5,-1.,0.5); t(-1.,0.,-1.,0.,-1.,0.); t(1.,0.5,-1.,0.5,1.,0.5); break;
        //case 14: t(-1.,0.5,1.,0.5,-1.,0.5); t(-1.,0.,-1.,0.,-1.,0.); t(1.,-0.5,-1.,-0.5,1.,-0.5); break;
        
        case 15: sc(y,pi,z4,v0); in(v0); gl(x,z,v0,z2); tr(x2y2); sc(y,pi,x4z4,y2); in(x4y4); gl(x,z,y4,x2z2); break; //C2/c
        
        //!--------------------------------------------------------------------------------------------------------
        //!Orthorhombisch:
        case 16: sc(z,pi,v0,v0); sc(y,pi,v0,v0); sc(x,pi,v0,v0); break; //P222
        
        case 17: sc(z,pi,v0,z2); sc(y,pi,z4,v0); sc(x,pi,v0,v0); break; //P222_1
        
        case 18: sc(z,pi,v0,v0); sc(y,pi,x4,y2); sc(x,pi,y4,x2); break; //P2_1 2_1 2
        
        case 19: sc(z,pi,x4,z2); sc(y,pi,z4,y2); sc(x,pi,y4,x2); break; //P2_1 2_1 2_1
        
        case 20: sc(z,pi,v0,z2); sc(y,pi,z4,v0); sc(x,pi,v0,v0); tr(x2y2);
                 sc(z,pi,x4y4,z2); sc(y,pi,x4z4,y2); sc(x,pi,y4,x2); break; //C222_1
        
        case 21: sc(z,pi,v0,v0); sc(y,pi,v0,v0); sc(x,pi,v0,v0); tr(x2y2);
                 sc(z,pi,x4y4,v0); sc(y,pi,x4,y2); sc(x,pi,y4,x2); break; //C222
        
        case 22: sc(z,pi,v0,v0); sc(y,pi,v0,v0); sc(x,pi,v0,v0); tr(y2z2);
                 sc(z,pi,y4,z2); sc(y,pi,z4,y2); sc(x,pi,y4z4,v0); tr(x2z2);
                 sc(z,pi,x4,z2); sc(y,pi,x4z4,v0); sc(x,pi,z4,x2); tr(x2y2);
                 sc(z,pi,x4y4,v0); sc(y,pi,x4,y2); sc(x,pi,y4,x2); break; //F222
        
        case 23: sc(z,pi,v0,v0); sc(y,pi,v0,v0); sc(x,pi,v0,v0); tr(x2y2z2);
                 sc(z,pi,x4y4,z2); sc(y,pi,x4z4,y2); sc(x,pi,y4z4,x2); break; //I222
        
        case 24: sc(z,pi,x4,z2); sc(y,pi,z4,y2); sc(x,pi,y4,x2); tr(x2y2z2);
                 sc(z,pi,y4,v0); sc(y,pi,x4,v0); sc(x,pi,z4,v0); break;//I2_12_12_1
        
        case 25: sc(z,pi,v0,v0); gl(x,z,v0,v0); gl(y,z,v0,v0); break; //Pmm2
        
        case 26: sc(z,pi,v0,z2); gl(x,z,v0,z2); gl(y,z,v0,v0); break; //Pmc2_1
        
        case 27: sc(z,pi,v0,v0); gl(x,z,v0,z2); gl(y,z,v0,z2); break; //Pcc2
        
        case 28: sc(z,pi,v0,v0); gl(x,z,v0,x2); gl(y,z,x4,v0); break; //Pma2
        
        case 29: sc(z,pi,v0,z2); gl(x,z,v0,x2); gl(y,z,x4,z2); break; //Pca2_1
        
        case 30: sc(z,pi,v0,v0); gl(x,z,y4,z2); gl(y,z,v0,y2z2); break; //Pnc2
        
        case 31: sc(z,pi,x4,z2); gl(x,z,v0,x2z2); gl(y,z,v0,v0); break; //Pmn2_1
        
        case 32: sc(z,pi,v0,v0); gl(x,z,y4,x2); gl(y,z,x4,y2); break; //Pba2
        
        case 33: sc(z,pi,v0,z2); gl(x,z,y4,x2); gl(y,z,x4,y2z2); break; //Pna2_1
        
        case 34: sc(z,pi,v0,v0); gl(x,z,y4,x2z2); gl(y,z,x4,y2z2); break; //Pnn2
        
        case 35: sc(z,pi,v0,v0); gl(x,z,v0,v0); gl(y,z,v0,v0); tr(x2y2);
                 sc(z,pi,x4y4,v0); gl(x,z,y4,x2); gl(y,z,x4,y2); break; //Cmm2
        
        case 36: sc(z,pi,v0,z2); gl(x,z,v0,z2); gl(y,z,v0,v0); tr(x2y2);
                 sc(z,pi,x4y4,z2); gl(x,z,y4,x2z2); gl(y,z,x4,y2); break; //Cmc2_1
        
        case 37: sc(z,pi,v0,v0); gl(x,z,v0,z2); gl(y,z,v0,z2); tr(x2y2);
                 sc(z,pi,x4y4,v0); gl(x,z,y4,x2z2); gl(y,z,x4,y2z2); break; //Ccc2
        
        case 38: sc(z,pi,v0,v0); gl(x,z,v0,v0); gl(y,z,v0,v0); tr(y2z2);
                 sc(z,pi,y4,z2); gl(x,z,y4,z2); gl(y,z,v0,y2z2); break; //Amm2
        
        case 39: sc(z,pi,v0,v0); gl(x,z,y4,v0); gl(y,z,v0,y2); tr(y2z2);
                 sc(z,pi,y4,z2); gl(x,z,v0,z2); gl(y,z,v0,z2); break; //Abm2
        
        case 40: sc(z,pi,v0,v0); gl(x,z,v0,x2); gl(y,z,x4,v0); tr(y2z2);
                 sc(z,pi,y4,z2); gl(x,z,y4,x2z2); gl(y,z,x4,y2z2); break; //Ama2
        
        case 41: sc(z,pi,v0,v0); gl(x,z,y4,x2); gl(y,z,x4,y2); tr(y2z2);
                 sc(z,pi,y4,z2); gl(x,z,v0,x2z2); gl(y,z,x4,z2); break; //Aba2
        
        case 42: sc(z,pi,v0,v0); gl(x,z,v0,v0); gl(y,z,v0,v0); tr(y2z2);
                         sc(z,pi,y4,z2); gl(x,z,y4,z2); gl(y,z,v0,y2z2); tr(x2z2);
                         sc(z,pi,x4,z2); gl(x,z,v0,x2z2); gl(y,z,x4,z2); tr(x2y2);
                         sc(z,pi,x4y4,v0); gl(x,z,y4,x2); gl(y,z,x4,y2); break; //Fmm2
        
        case 43: sc(z,pi,v0,v0); gl(x,z,y8,x4z4); gl(y,z,x8,y4z4); tr(y2z2);
                         sc(z,pi,y4,z2); gl(x,z,y38,x4z34); gl(y,z,x8,y34z34); tr(x2z2);
                         sc(z,pi,x4,z2); gl(x,z,y8,x34z34); gl(y,z,x38,y4z34); tr(x2y2);
                         sc(z,pi,x4y4,v0); gl(x,z,y38,x34z4); gl(y,z,x38,y34z4); break; //Fdd2
        
        case 44: sc(z,pi,v0,v0); gl(x,z,v0,v0); gl(y,z,v0,v0); tr(x2y2z2);
                 sc(z,pi,x4y4,z2); gl(x,z,y4,x2z2); gl(y,z,x4,y2z2); break; //Imm2
        
        case 45: sc(z,pi,v0,v0); gl(x,z,y4,x2); gl(y,z,x4,y2); tr(x2y2z2);
                 sc(z,pi,x4y4,z2); gl(x,z,v0,z2); gl(y,z,v0,z2); break; //Iba2
        
        case 46: sc(z,pi,v0,v0); gl(x,z,v0,x2); gl(y,z,x4,v0); tr(x2y2z2);
                 sc(z,pi,x4y4,z2); gl(x,z,y4,z2); gl(y,z,v0,y2z2); break; //Ima2
        
        case 47: sc(z,pi,v0,v0); sc(y,pi,v0,v0); sc(x,pi,v0,v0); in(v0);
                 gl(x,y,v0,v0); gl(x,z,v0,v0); gl(y,z,v0,v0); break; //Pmmm
        
        //korrekte Origin Choice ermitteln:
        case 48: sc(z,pi,v0,v0); sc(y,pi,v0,v0); sc(x,pi,v0,v0); in(x4y4z4);
                 gl(x,y,z4,x2y2); gl(x,z,y4,x2z2); gl(y,z,x4,y2z2); //Pnnn
                 if (origin_one()) break;
                 else {
                kill_oc1();
                sc(z,pi,x4y4,v0); sc(y,pi,x4z4,v0); sc(x,pi,y4z4,v0); in(v0);
                gl(x,y,v0,x2y2); gl(x,z,v0,x2z2); gl(y,z,v0,y2z2); break; //Pnnn
                 }
        
        case 49: sc(z,pi,v0,v0); sc(y,pi,z4,v0); sc(x,pi,z4,v0); in(v0);
                 gl(x,y,v0,v0); gl(x,z,v0,z2); gl(y,z,v0,z2); break; //Pccm
        
        case 50: sc(z,pi,v0,v0); sc(y,pi,v0,v0); sc(x,pi,v0,v0); in(x4y4);
                 gl(x,y,v0,x2y2); gl(x,z,y4,x2); gl(y,z,x4,y2); //Pban
                 if (origin_one()) break;
                 else {
                kill_oc1();
                sc(z,pi,x4y4,v0); sc(y,pi,x4,v0); sc(x,pi,y4,v0); in(v0);
                gl(x,y,v0,x2y2); gl(x,z,v0,x2); gl(y,z,v0,y2); break; //Pban
                 }
        
        
        case 51: sc(z,pi,x4,v0); sc(y,pi,v0,v0); sc(x,pi,v0,x2); in(v0);
                 gl(x,y,v0,x2); gl(x,z,v0,v0); gl(y,z,x4,v0); break; //Pmma
        //case 51: sc(y,pi,z4,v0); sc(x,pi,v0,v0); sc(z,pi,v0,x2); in(v0);
        //         gl(x,z,v0,z2); gl(y,z,v0,v0); gl(x,y,z4,v0); break; //Pmma
        
        case 52: sc(z,pi,x4,v0); sc(y,pi,x4z4,y2); sc(x,pi,y4z4,v0); in(v0);
                 gl(x,y,v0,x2); gl(x,z,y4,x2z2); gl(y,z,v0,y2z2); break; //Pnna
        
        case 53: sc(z,pi,x4,z2); sc(y,pi,x4z4,v0); sc(x,pi,v0,v0); in(v0);
                 gl(x,y,z4,x2); gl(x,z,v0,x2z2); gl(y,z,v0,v0); break; //Pmna
        
        case 54: sc(z,pi,x4,v0); sc(y,pi,z4,v0); sc(x,pi,z4,x2); in(v0);
                 gl(x,y,v0,x2); gl(x,z,v0,z2); gl(y,z,x4,z2); break; //Pcca
        
        case 55: sc(z,pi,v0,v0); sc(y,pi,x4,y2); sc(x,pi,y4,x2); in(v0);
                 gl(x,y,v0,v0); gl(x,z,y4,x2); gl(y,z,x4,y2); break; //Pbam
        
        case 56: sc(z,pi,x4y4,v0); sc(y,pi,z4,y2); sc(x,pi,z4,x2); in(v0);
                 gl(x,y,v0,x2y2); gl(x,z,y4,z2); gl(y,z,x4,z2); break; //Pccn
        
        case 57: sc(z,pi,v0,z2); sc(y,pi,z4,y2); sc(x,pi,y4,v0); in(v0);
                 gl(x,y,z4,v0); gl(x,z,y4,z2); gl(y,z,v0,y2); break; //Pbcm
        
        case 58: sc(z,pi,v0,v0); sc(y,pi,x4z4,y2); sc(x,pi,y4z4,x2); in(v0);
                 gl(x,y,v0,v0); gl(x,z,y4,x2z2); gl(y,z,x4,y2z2); break; //Pnnm
        
        case 59: sc(z,pi,v0,v0); sc(y,pi,x4,y2); sc(x,pi,y4,x2); in(x4y4);
                 gl(x,y,v0,x2y2); gl(x,z,v0,v0); gl(y,z,v0,v0); //Pmmn
                 if (origin_one()) break;
                 else {
                kill_oc1();
                sc(z,pi,x4y4,v0); sc(y,pi,v0,y2); sc(x,pi,v0,x2); in(v0);
                gl(x,y,v0,x2y2); gl(x,z,y4,v0); gl(y,z,x4,v0); break; //Pmmn
            }
        
        case 60: sc(z,pi,x4y4,z2); sc(y,pi,z4,v0); sc(x,pi,y4,x2); in(v0);
                 gl(x,y,z4,x2y2); gl(x,z,v0,z2); gl(y,z,x4,y2); break; //Pbcn
        
        case 61: sc(z,pi,x4,z2); sc(y,pi,z4,y2); sc(x,pi,y4,x2); in(v0);
                 gl(x,y,z4,x2); gl(x,z,y4,z2); gl(y,z,x4,y2); break; //Pbca
        
        case 62: sc(z,pi,x4,z2); sc(y,pi,v0,y2); sc(x,pi,y4z4,x2); in(v0);
                 gl(x,y,z4,x2); gl(x,z,y4,v0); gl(y,z,x4,y2z2); break; //Pnma
        
        case 63: sc(z,pi,v0,z2); sc(y,pi,z4,v0); sc(x,pi,v0,v0); in(v0); gl(x,y,z4,v0);
                 gl(x,z,v0,z2); gl(y,z,v0,v0); tr(x2y2); sc(z,pi,x4y4,z2); sc(y,pi,x4z4,y2);
                 sc(x,pi,y4,x2); in(x4y4); gl(x,y,z4,x2y2); gl(x,z,y4,x2z2); gl(y,z,x4,y2); break; //Cmcm
        
        case 64: sc(z,pi,y4,z2); sc(y,pi,z4,y2); sc(x,pi,v0,v0); in(v0); gl(x,y,z4,y2);
                 gl(x,z,y4,z2); gl(y,z,v0,v0); tr(x2y2); sc(z,pi,x4,z2); sc(y,pi,x4z4,v0);
                 sc(x,pi,y4,x2); in(x4y4); gl(x,y,z4,x2); gl(x,z,v0,x2z2); gl(y,z,x4,y2); break; //Cmca
        
        case 65: sc(z,pi,v0,v0); sc(y,pi,v0,v0); sc(x,pi,v0,v0); in(v0); gl(x,y,v0,v0);
                 gl(x,z,v0,v0); gl(y,z,v0,v0); tr(x2y2); sc(z,pi,x4y4,v0); sc(y,pi,x4,y2);
                 sc(x,pi,y4,x2); in(x4y4); gl(x,y,v0,x2y2); gl(x,z,y4,x2); gl(y,z,x4,y2); break; //Cmmm
        
        case 66: sc(z,pi,v0,v0); sc(y,pi,z4,v0); sc(x,pi,z4,v0); in(v0); gl(x,y,v0,v0);
                 gl(x,z,v0,z2); gl(y,z,v0,z2); tr(x2y2); sc(z,pi,x4y4,v0); sc(y,pi,x4z4,y2); sc(x,pi,y4z4,x2);
                 in(x4y4); gl(x,y,v0,x2y2); gl(x,z,y4,x2z2); gl(y,z,x4,y2z2); break; //Cccm
        
        case 67: sc(z,pi,y4,v0); sc(y,pi,v0,y2); sc(x,pi,v0,v0); in(v0); gl(x,y,v0,y2); gl(x,z,y4,v0); gl(y,z,v0,v0); tr(x2y2);
                         sc(z,pi,x4,v0); sc(y,pi,x4,v0); sc(x,pi,y4,x2); in(x4y4); gl(x,y,v0,x2); gl(x,z,v0,x2); gl(y,z,x4,y2); break; //Cmma
        
        case 68: sc(z,pi,x4y4,v0); sc(y,pi,v0,v0); sc(x,pi,y4,x2); in(y4z4); gl(x,y,z4,x2); gl(x,z,y4,z2); gl(y,z,x4,z2); tr(x2y2);
                         sc(z,pi,v0,v0); sc(y,pi,x4,y2); sc(x,pi,v0,v0); in(x4z4); gl(x,y,z4,y2); gl(x,z,v0,x2z2); gl(y,z,v0,y2z2); //Ccca
                         if (origin_one()) break;
                 else {
                kill_oc1();
                sc(z,pi,x4,v0); sc(y,pi,z4,v0); sc(x,pi,z4,x2); in(v0); gl(x,y,v0,x2); gl(x,z,v0,z2); gl(y,z,x4,z2); tr(x2y2);
                sc(z,pi,y4,v0); sc(y,pi,x4z4,y2); sc(x,pi,y4z4,v0); in(x4y4); gl(x,y,v0,y2); gl(x,z,y4,x2z2); gl(y,z,v0,y2z2);
                break;
                 }
        
        case 69: sc(z,pi,v0,v0); sc(y,pi,v0,v0); sc(x,pi,v0,v0); in(v0); gl(x,y,v0,v0); gl(x,z,v0,v0); gl(y,z,v0,v0); tr(y2z2);
                         sc(z,pi,y4,z2); sc(y,pi,z4,y2); sc(x,pi,y4z4,v0); in(y4z4); gl(x,y,z4,y2); gl(x,z,y4,z2); gl(y,z,v0,y2z2); tr(x2z2);
                         sc(z,pi,x4,z2); sc(y,pi,x4z4,v0); sc(x,pi,z4,x2); in(x4z4); gl(x,y,z4,x2); gl(x,z,v0,x2z2); gl(y,z,x4,z2); tr(x2y2);
                         sc(z,pi,x4y4,v0); sc(y,pi,x4,y2); sc(x,pi,y4,x2); in(x4y4); gl(x,y,v0,x2y2); gl(x,z,y4,x2); gl(y,z,x4,y2); break; //Fmmm 
        
        case 70: sc(z,pi,v0,v0); sc(y,pi,v0,v0); sc(x,pi,v0,v0); in(x8y8z8); gl(x,y,z8,x4y4); gl(x,z,y8,x4z4);
                 gl(y,z,x8,y4z4); tr(y2z2); sc(z,pi,y4,z2); sc(y,pi,z4,y2); sc(x,pi,y4z4,v0); in(x8y38z38); gl(x,y,z38,x4y34);
                 gl(x,z,y38,x4z34); gl(y,z,x8,y34z34); tr(x2z2); sc(z,pi,x4,z2); sc(y,pi,x4z4,v0); sc(x,pi,z4,x2); in(x38y8z38);
                 gl(x,y,z38,x34y4); gl(x,z,y8,x34z34); gl(y,z,x38,y4z34); tr(x2y2); sc(z,pi,x4y4,v0); sc(y,pi,x4,y2); sc(x,pi,y4,x2);
                 in(x38y38z8); gl(x,y,z8,x34y34); gl(x,z,y38,x34z4); gl(y,z,x38,y34z4); //Fddd
                 if (origin_one()) break;
                 else {
                kill_oc1();
                sc(z,pi,x38y38,v0); sc(y,pi,x38z38,v0); sc(x,pi,y38z38,v0); in(v0); gl(x,y,v0,x4y4); gl(x,z,v0,x4z4);
                gl(y,z,v0,y4z4); tr(y2z2); sc(z,pi,x38y8,z2); sc(y,pi,x38z8,y2); sc(x,pi,y8z8,v0); in(y4z4); gl(x,y,z4,x4y34);
                gl(x,z,y4,x4z34); gl(y,z,v0,y34z34); tr(x2z2); sc(z,pi,x8y38,z2); sc(y,pi,x8z8,v0); sc(x,pi,y38z8,x2); in(x4z4);
                gl(x,y,z4,x34y4); gl(x,z,v0,x34z34); gl(y,z,x4,y4z34); tr(x2y2); sc(z,pi,x8y8,v0); sc(y,pi,x8z38,y2); 
                sc(x,pi,y8z38,x2); in(x4y4); gl(x,y,v0,x34y34); gl(x,z,y4,x34z4); gl(y,z,x4,y34z4); //Fddd
                break;
            }
        
        case 71: sc(z,pi,v0,v0); sc(y,pi,v0,v0); sc(x,pi,v0,v0); in(v0); gl(x,y,v0,v0); gl(x,z,v0,v0); gl(y,z,v0,v0);
                 tr(x2y2z2); sc(z,pi,x4y4,z2); sc(y,pi,x4z4,y2); sc(x,pi,y4z4,x2); in(x4y4z4); gl(x,y,z4,x2y2);
                 gl(x,z,y4,x2z2); gl(y,z,x4,y2z2); break; //Immm
        
        case 72: sc(z,pi,v0,v0); sc(y,pi,x4,y2); sc(x,pi,y4,x2); in(v0); gl(x,y,v0,v0); gl(x,z,y4,x2); gl(y,z,x4,y2); tr(x2y2z2);
                         sc(z,pi,x4y4,z2); sc(y,pi,z4,v0); sc(x,pi,z4,v0); in(x4y4z4); gl(x,y,z4,x2y2); gl(x,z,v0,z2); gl(y,z,v0,z2); break; //Ibam
        
        case 73: sc(z,pi,x4,z2); sc(y,pi,z4,y2); sc(x,pi,y4,x2); in(v0); gl(x,y,z4,x2); gl(x,z,y4,z2); gl(y,z,x4,y2); tr(x2y2z2);
                         sc(z,pi,y4,v0); sc(y,pi,x4,v0); sc(x,pi,z4,v0); in(x4y4z4); gl(x,y,v0,y2); gl(x,z,v0,x2); gl(y,z,v0,z2); break; //Ibca
        
        case 74: sc(z,pi,y4,v0); sc(y,pi,v0,y2); sc(x,pi,v0,v0); in(v0); gl(x,y,v0,y2); gl(x,z,y4,v0);
                 gl(y,z,v0,v0); tr(x2y2z2); sc(z,pi,x4,z2); sc(y,pi,x4z4,v0); sc(x,pi,y4z4,x2); in(x4y4z4);
                 gl(x,y,z4,x2); gl(x,z,v0,x2z2); gl(y,z,x4,y2z2); break; //Imma
        
        //!--------------------------------------------------------------------------------------------------------
        //!Tetragonal:
        case 75: sc(z,pi,v0,v0); sc(z,pi2,v0,v0); sc(z,pi32,v0,v0); break; //P4   3x um 90 drehen um z
        
        case 76: sc(z,pi,v0,z2); sc(z,pi2,v0,z4); sc(z,pi32,v0,z34); break; //P4_1
        
        case 77: sc(z,pi,v0,v0); sc(z,pi2,v0,z2); sc(z,pi32,v0,z2); break; //P4_2
        
        case 78: sc(z,pi,v0,z2); sc(z,pi2,v0,z34); sc(z,pi32,v0,z4); break; //P4_3
        
        case 79: sc(z,pi,v0,v0); sc(z,pi2,v0,v0); sc(z,pi32,v0,v0); tr(x2y2z2);
                 sc(z,pi,x4y4,z2); sc(z,pi2,y2,z2); sc(z,pi32,x2,z2); break; //I4
        
        case 80: sc(z,pi,x4y4,z2); sc(z,pi2,ix4y4,z4); sc(z,pi32,x4iy4,z34); tr(x2y2z2);
                 sc(z,pi,v0,v0); sc(z,pi2,x4y4,z34); sc(z,pi32,x4y4,z4); break; //I4_1 
        
        case 81: sc(z,pi,v0,v0); sci(z,pi2,v0,v0,v0); sci(z,pi32,v0,v0,v0); break; //P4_inv  -- Rotoinversion
        
        case 82: sc(z,pi,v0,v0); sci(z,pi2,v0,v0,v0); sci(z,pi32,v0,v0,v0); tr(x2y2z2);
                 sc(z,pi,x4y4,z2); sci(z,pi2,x2,v0,x2z4); sci(z,pi32,y2,v0,y2z4); break; //I4_inv
        
        case 83: sc(z,pi,v0,v0); sc(z,pi2,v0,v0); sc(z,pi32,v0,v0); in(v0); gl(x,y,v0,v0);
                 sci(z,pi2,v0,v0,v0); sci(z,pi32,v0,v0,v0); break; //P4/m
        
        case 84: sc(z,pi,v0,v0); sc(z,pi2,v0,z2); sc(z,pi32,v0,z2); in(v0); gl(x,y,v0,v0);
                 sci(z,pi2,v0,v0,z4); sci(z,pi32,v0,v0,z4); break; //P4_2/m
        
        case 85: sc(z,pi,v0,v0); sc(z,pi2,y2,v0); sc(z,pi32,x2,v0); in(x4y4); gl(x,y,v0,x2y2);
                 sci(z,pi2,v0,v0,v0); sci(z,pi32,v0,v0,v0); //P4/n
                 if (origin_one()) break;
                 else {
                kill_oc1();
                sc(z,pi,x4y4,v0); sc(z,pi2,x4y4,v0); sc(z,pi32,x4y4,v0); in(v0); gl(x,y,v0,x2y2);
                    sci(z,pi2,x4iy4,v0,x4iy4); sci(z,pi32,ix4y4,v0,ix4y4); //P4/n
                    break;
            }
        
        case 86: sc(z,pi,v0,v0); sc(z,pi2,y2,z2); sc(z,pi32,x2,z2); in(x4y4z4); gl(x,y,z4,x2y2);
                 sci(z,pi2,v0,v0,v0); sci(z,pi32,v0,v0,v0); //P4_2/n
                 if (origin_one()) break;
                 else {
                kill_oc1();
                sc(z,pi,x4y4,v0); sc(z,pi2,ix4y4,z2); sc(z,pi32,x4iy4,z2); in(v0); gl(x,y,v0,x2y2);
                sci(z,pi2,x4y4,v0,x4y4z4); sci(z,pi32,x4y4,v0,x4y4z4); //P4_2/n
                    break;
            }
        
        case 87: sc(z,pi,v0,v0); sc(z,pi2,v0,v0); sc(z,pi32,v0,v0); in(v0); gl(x,y,v0,v0); sci(z,pi2,v0,v0,v0);
                 sci(z,pi32,v0,v0,v0); tr(x2y2z2); sc(z,pi,x4y4,z2); sc(z,pi2,y2,z2); sc(z,pi32,x2,z2); in(x4y4z4);
                 gl(x,y,z4,x2y2); sci(z,pi2,x2,v0,x2z4); sci(z,pi32,y2,v0,y2z4); break; //I4/m
        
        case 88: sc(z,pi,x4y4,z2); sc(z,pi2,ix4y4,z4); sc(z,pi32,x4iy4,z34); in(y4z8); gl(x,y,z38,x2); sci(z,pi2,v0,v0,v0);
                 sci(z,pi32,y2,v0,y2z4); tr(x2y2z2); sc(z,pi,v0,v0); sc(z,pi2,x4y4,z34); sc(z,pi32,x4y4,z4);
                 in(x4z38); gl(x,y,z8,y2); sci(z,pi2,x2,v0,x2z4); sci(z,pi32,v0,v0,v0); //I4_1/a
                 if (origin_one()) break;
                 else {
                kill_oc1();
                sc(z,pi,x4,z2); sc(z,pi2,x4y2,z4); sc(z,pi32,x34,z34); in(v0); gl(x,y,z4,x2); sci(z,pi2,x2y4,v0,x2y4z38);
                sci(z,pi32,y4,v0,y4z8); tr(x2y2z2); sc(z,pi,y4,v0); sc(z,pi2,ix4y2,z34); sc(z,pi32,x4,z4);
                in(x4y4z4); gl(x,y,v0,y2); sci(z,pi2,x2iy4,v0,x2iy4z8); sci(z,pi32,y34,v0,y34z38); //I4_1/a
                break;
            }
        
        case 89: sc(z,pi,v0,v0); sc(z,pi2,v0,v0); sc(z,pi32,v0,v0); sc(y,pi,v0,v0); sc(x,pi,v0,v0);
                 sc(xy,pi,v0,v0); sc(xiy,pi,v0,v0);  break;//P422
        
        case 90: sc(z,pi,v0,v0); sc(z,pi2,y2,v0); sc(z,pi32,x2,v0); sc(y,pi,x4,y2); sc(x,pi,y4,x2);
                 sc(xy,pi,v0,v0); sc(xiy,pi,v0,v0); break; //P42_1
        
        case 91: sc(z,pi,v0,z2); sc(z,pi2,v0,z4); sc(z,pi32,v0,z34); sc(y,pi,v0,v0); sc(x,pi,z4,v0);
                 sc(xy,pi,z38,v0); sc(xiy,pi,z8,v0); break; //P4_122
        
        case 92: sc(z,pi,v0,z2); sc(z,pi2,y2,z4); sc(z,pi32,x2,z34); sc(y,pi,x4z8,y2); sc(x,pi,y4z38,x2);
                 sc(xy,pi,v0,v0); sc(xiy,pi,z4,v0); break; //P4_12_12
        
        case 93: sc(z,pi,v0,v0); sc(z,pi2,v0,z2); sc(z,pi32,v0,z2); sc(y,pi,v0,v0); sc(x,pi,v0,v0);
                 sc(xy,pi,z4,v0); sc(xiy,pi,z4,v0); break; //P4_222
        
        case 94: sc(z,pi,v0,v0); sc(z,pi2,y2,z2); sc(z,pi32,x2,z2); sc(y,pi,x4z4,y2); sc(x,pi,y4z4,x2);
                 sc(xy,pi,v0,v0); sc(xiy,pi,v0,v0); break; //P4_22_12
        
        case 95: sc(z,pi,v0,z2); sc(z,pi2,v0,z34); sc(z,pi32,v0,z4); sc(y,pi,v0,v0); sc(x,pi,z4,v0);
                 sc(xy,pi,z8,v0); sc(xiy,pi,z38,v0); break; //P4_322
        
        case 96: sc(z,pi,v0,z2); sc(z,pi2,y2,z34); sc(z,pi32,x2,z4); sc(y,pi,x4z38,y2); sc(x,pi,y4z8,x2);
                 sc(xy,pi,v0,v0); sc(xiy,pi,z4,v0); break; //P4_32_12
        
        case 97: sc(z,pi,v0,v0); sc(z,pi2,v0,v0); sc(z,pi32,v0,v0); sc(y,pi,v0,v0); sc(x,pi,v0,v0); sc(xy,pi,v0,v0); sc(xiy,pi,v0,v0); 
             tr(x2y2z2); sc(z,pi,x4y4,z2); sc(z,pi2,y2,z2); sc(z,pi32,x2,z2); sc(y,pi,x4z4,y2); sc(x,pi,y4z4,x2); sc(xy,pi,z4,x2y2);
             sc(xiy,pi,y2z4,v0); break; //I422
        
        case 98: sc(z,pi,x4y4,z2); sc(z,pi2,ix4y4,z4); sc(z,pi32,x4iy4,z34); sc(y,pi,x4z38,v0); sc(x,pi,y4z8,v0);
                 sc(xy,pi,z4,x2y2); sc(xiy,pi,v0,v0); tr(x2y2z2); sc(z,pi,v0,v0); sc(z,pi2,x4y4,z34); sc(z,pi32,x4y4,z4);
                 sc(y,pi,z8,y2); sc(x,pi,z38,x2); sc(xy,pi,v0,v0); sc(xiy,pi,y2z4,v0); break; //I4_122
        
        case 99: sc(z,pi,v0,v0); sc(z,pi2,v0,v0); sc(z,pi32,v0,v0); gl(x,z,v0,v0); gl(y,z,v0,v0);
                 gl(xiy,z,v0,v0); gl(xy,z,v0,v0); break; //P4mm
        
        case 100: sc(z,pi,v0,v0); sc(z,pi2,v0,v0); sc(z,pi32,v0,v0); gl(x,z,y4,x2); gl(y,z,x4,y2);
                  gl(xiy,z,x2,v0); gl(xy,z,v0,x2y2); break; //P4bm
        
        case 101: sc(z,pi,v0,v0); sc(z,pi2,v0,z2); sc(z,pi32,v0,z2); gl(x,z,v0,z2); gl(y,z,v0,z2);
                  gl(xiy,z,v0,v0); gl(xy,z,v0,v0); break; //P4_2cm
        
        case 102: sc(z,pi,v0,v0); sc(z,pi2,y2,z2); sc(z,pi32,x2,z2); gl(x,z,y4,x2z2); gl(y,z,x4,y2z2);
                  gl(xiy,z,v0,v0); gl(xy,z,v0,v0); break; //P4_2nm
        
        case 103: sc(z,pi,v0,v0); sc(z,pi2,v0,v0); sc(z,pi32,v0,v0); gl(x,z,v0,z2); gl(y,z,v0,z2);
                  gl(xiy,z,v0,z2); gl(xy,z,v0,z2); break; //P4cc
        
        case 104: sc(z,pi,v0,v0); sc(z,pi2,v0,v0); sc(z,pi32,v0,v0); gl(x,z,y4,x2z2); gl(y,z,x4,y2z2);
                  gl(xiy,z,x2,z2); gl(xy,z,v0,x2y2z2); break; //P4nc
        
        case 105: sc(z,pi,v0,v0); sc(z,pi2,v0,z2); sc(z,pi32,v0,z2); gl(x,z,v0,v0); gl(y,z,v0,v0);
                  gl(xiy,z,v0,z2); gl(xy,z,v0,z2); break; //P4_2mc
        
        case 106: sc(z,pi,v0,v0); sc(z,pi2,v0,z2); sc(z,pi32,v0,z2); gl(x,z,y4,x2); gl(y,z,x4,y2);
                  gl(xiy,z,x2,z2); gl(xy,z,v0,x2y2z2); break; //P4_2bc
        
        case 107: sc(z,pi,v0,v0); sc(z,pi2,v0,v0); sc(z,pi32,v0,v0); gl(x,z,v0,v0); gl(y,z,v0,v0); gl(xiy,z,v0,v0); gl(xy,z,v0,v0); 
                  tr(x2y2z2); sc(z,pi,x4y4,z2); sc(z,pi2,y2,z2); sc(z,pi32,x2,z2); gl(x,z,y4,x2z2); gl(y,z,x4,y2z2); gl(xiy,z,x2,z2); 
                  gl(xy,z,v0,x2y2z2); break; //I4mm
        
        case 108: sc(z,pi,v0,v0); sc(z,pi2,v0,v0); sc(z,pi32,v0,v0); gl(x,z,v0,z2); gl(y,z,v0,z2); gl(xiy,z,v0,z2); gl(xy,z,v0,z2); 
                  tr(x2y2z2); sc(z,pi,x4y4,z2); sc(z,pi2,y2,z2); sc(z,pi32,x2,z2); gl(x,z,y4,x2); gl(y,z,x4,y2); gl(xiy,z,x2,v0); 
                  gl(xy,z,v0,x2y2); break; //I4cm
        
        case 109: sc(z,pi,x4y4,z2); sc(z,pi2,ix4y4,z4); sc(z,pi32,x4iy4,z34); gl(x,z,v0,v0); gl(y,z,x4,y2z2); gl(xiy,z,x4,ix4y4z4);
                  gl(xy,z,x4,x4y4z34); tr(x2y2z2); sc(z,pi,v0,v0); sc(z,pi2,x4y4,z34); sc(z,pi32,x4y4,z4); gl(x,z,y4,x2z2); gl(y,z,v0,v0);
                  gl(xiy,z,x4,x4iy4z34); gl(xy,z,ix4,x4y4z4); break; //I4_1md
        
        case 110: sc(z,pi,x4y4,z2); sc(z,pi2,ix4y4,z4); sc(z,pi32,x4iy4,z34); gl(x,z,v0,z2); gl(y,z,x4,y2); gl(xiy,z,x4,ix4y4z34);
                  gl(xy,z,x4,x4y4z4); tr(x2y2z2); sc(z,pi,v0,v0); sc(z,pi2,x4y4,z34); sc(z,pi32,x4y4,z4); gl(x,z,y4,x2);
                  gl(y,z,v0,z2); gl(xiy,z,x4,x4iy4z4); gl(xy,z,ix4,x4y4z34); break; //I4_1cd
        
        case 111: sc(z,pi,v0,v0); sci(z,pi2,v0,v0,v0); sci(z,pi32,v0,v0,v0); sc(y,pi,v0,v0); sc(x,pi,v0,v0); gl(xiy,z,v0,v0);
                  gl(xy,z,v0,v0); break; //P4_inv2m
        
        case 112: sc(z,pi,v0,v0); sci(z,pi2,v0,v0,v0); sci(z,pi32,v0,v0,v0); sc(y,pi,z4,v0); sc(x,pi,z4,v0); gl(xiy,z,v0,z2);
                  gl(xy,z,v0,z2); break; //P4_inv2c
        
        case 113: sc(z,pi,v0,v0); sci(z,pi2,v0,v0,v0); sci(z,pi32,v0,v0,v0); sc(y,pi,x4,y2); sc(x,pi,y4,x2); gl(xiy,z,x2,v0);
                  gl(xy,z,v0,x2y2); break; //P4_inv2_1m
        
        case 114: sc(z,pi,v0,v0); sci(z,pi2,v0,v0,v0); sci(z,pi32,v0,v0,v0); sc(y,pi,x4z4,y2); sc(x,pi,y4z4,x2); gl(xiy,z,x2,z2);
                  gl(xy,z,v0,x2y2z2); break; //P4_inv2_1c
        
        case 115: sc(z,pi,v0,v0); sci(z,pi2,v0,v0,v0); sci(z,pi32,v0,v0,v0); gl(x,z,v0,v0); gl(y,z,v0,v0); sc(xy,pi,v0,v0);
                  sc(xiy,pi,v0,v0); break; //P4_invm2
        
        case 116: sc(z,pi,v0,v0); sci(z,pi2,v0,v0,v0); sci(z,pi32,v0,v0,v0); gl(x,z,v0,z2); gl(y,z,v0,z2); sc(xy,pi,z4,v0);
                  sc(xiy,pi,z4,v0); break; //P4_invc2
        
        case 117: sc(z,pi,v0,v0); sci(z,pi2,v0,v0,v0); sci(z,pi32,v0,v0,v0); gl(x,z,y4,x2); gl(y,z,x4,y2); sc(xy,pi,v0,x2y2);
                  sc(xiy,pi,y2,v0); break; //P4_invb2
        
        case 118: sc(z,pi,v0,v0); sci(z,pi2,v0,v0,v0); sci(z,pi32,v0,v0,v0); gl(x,z,y4,x2z2); gl(y,z,x4,y2z2); sc(xy,pi,z4,x2y2);
                  sc(xiy,pi,y2z4,v0); break; //P4_invn2
        
        case 119: sc(z,pi,v0,v0); sci(z,pi2,v0,v0,v0); sci(z,pi32,v0,v0,v0); gl(x,z,v0,v0); gl(y,z,v0,v0); sc(xy,pi,v0,v0);
                  sc(xiy,pi,v0,v0); tr(x2y2z2); sc(z,pi,x4y4,z2); sci(z,pi2,x2,v0,x2z4); sci(z,pi32,y2,v0,y2z4); gl(x,z,y4,x2z2);
                  gl(y,z,x4,y2z2); sc(xy,pi,z4,x2y2); sc(xiy,pi,y2z4,v0);break; //I4_invm2
        
        case 120: sc(z,pi,v0,v0); sci(z,pi2,v0,v0,v0); sci(z,pi32,v0,v0,v0); gl(x,z,v0,z2); gl(y,z,v0,z2); sc(xy,pi,z4,v0);
                  sc(xiy,pi,z4,v0); tr(x2y2z2); sc(z,pi,x4y4,z2); sci(z,pi2,x2,v0,x2z4); sci(z,pi32,y2,v0,y2z4); gl(x,z,y4,x2);
                  gl(y,z,x4,y2); sc(xy,pi,v0,x2y2); sc(xiy,pi,y2,v0);break; //I4_invc2
        
        case 121: sc(z,pi,v0,v0); sci(z,pi2,v0,v0,v0); sci(z,pi32,v0,v0,v0); sc(y,pi,v0,v0); sc(x,pi,v0,v0); gl(xiy,z,v0,v0);
                  gl(xy,z,v0,v0); tr(x2y2z2); sc(z,pi,x4y4,z2); sci(z,pi2,x2,v0,x2z4); sci(z,pi32,y2,v0,y2z4); sc(y,pi,x4z4,y2);
                  sc(x,pi,y4z4,x2); gl(xiy,z,x2,z2); gl(xy,z,v0,x2y2z2); break; //I4_inv2m
        
        case 122: sc(z,pi,v0,v0); sci(z,pi2,v0,v0,v0); sci(z,pi32,v0,v0,v0); sc(y,pi,x4z38,v0); sc(x,pi,z38,x2); gl(xiy,z,x4,x4iy4z34);
                  gl(xy,z,x4,x4y4z34); tr(x2y2z2); sc(z,pi,x4y4,z2); sci(z,pi2,x2,v0,x2z4); sci(z,pi32,y2,v0,y2z4); sc(y,pi,z8,y2);
                  sc(x,pi,y4z8,v0); gl(xiy,z,x4,ix4y4z4); gl(xy,z,ix4,x4y4z4); break; //I4_inv2d
        
        case 123: sc(z,pi,v0,v0); sc(z,pi2,v0,v0); sc(z,pi32,v0,v0); sc(y,pi,v0,v0); sc(x,pi,v0,v0); sc(xy,pi,v0,v0); sc(xiy,pi,v0,v0);
                  in(v0); gl(x,y,v0,v0); sci(z,pi2,v0,v0,v0); sci(z,pi32,v0,v0,v0); gl(x,z,v0,v0); gl(y,z,v0,v0); gl(xiy,z,v0,v0);
                  gl(xy,z,v0,v0); break; //P4/mmm
        
        case 124: sc(z,pi,v0,v0); sc(z,pi2,v0,v0); sc(z,pi32,v0,v0); sc(y,pi,z4,v0); sc(x,pi,z4,v0); sc(xy,pi,z4,v0); sc(xiy,pi,z4,v0);
                  in(v0); gl(x,y,v0,v0); sci(z,pi2,v0,v0,v0); sci(z,pi32,v0,v0,v0); gl(x,z,v0,z2); gl(y,z,v0,z2); gl(xiy,z,v0,z2);
                  gl(xy,z,v0,z2); break; //P4/mcc
        
        case 125: sc(z,pi,v0,v0); sc(z,pi2,v0,v0); sc(z,pi32,v0,v0); sc(y,pi,v0,v0); sc(x,pi,v0,v0); sc(xy,pi,v0,v0); sc(xiy,pi,v0,v0);
                  in(x4y4); gl(x,y,v0,x2y2); sci(z,pi2,x2,v0,x2); sci(z,pi32,y2,v0,y2); gl(x,z,y4,x2); gl(y,z,x4,y2); gl(xiy,z,x2,v0);
                  gl(xy,z,v0,x2y2); //P4/nbm
                  if (origin_one()) break;
                  else {
                kill_oc1();
                sc(z,pi,x4y4,v0); sc(z,pi2,x4y4,v0); sc(z,pi32,x4y4,v0); sc(y,pi,x4,v0); sc(x,pi,y4,v0); sc(xy,pi,v0,v0); 
                sc(xiy,pi,y2,v0); in(v0); gl(x,y,v0,x2y2); sci(z,pi2,x4iy4,v0,x4iy4); sci(z,pi32,ix4y4,v0,ix4y4);
                gl(x,z,v0,x2); gl(y,z,v0,y2); gl(xiy,z,v0,v0); gl(xy,z,v0,x2y2); //P4/nbm
                break;
            }
        
        case 126: sc(z,pi,v0,v0); sc(z,pi2,v0,v0); sc(z,pi32,v0,v0); sc(y,pi,v0,v0); sc(x,pi,v0,v0); sc(xy,pi,v0,v0); sc(xiy,pi,v0,v0);
                  in(x4y4z4); gl(x,y,z4,x2y2); sci(z,pi2,x2,v0,x2z4); sci(z,pi32,y2,v0,y2z4); gl(x,z,y4,x2z2); gl(y,z,x4,y2z2);
                  gl(xiy,z,x2,z2); gl(xy,z,v0,x2y2z2); //P4/nnc
                  if (origin_one()) break;
                  else {
                kill_oc1();
                sc(z,pi,x4y4,v0); sc(z,pi2,x4y4,v0); sc(z,pi32,x4y4,v0); sc(y,pi,x4z4,v0); sc(x,pi,y4z4,v0); sc(xy,pi,z4,v0); 
                sc(xiy,pi,y2z4,v0); in(v0); gl(x,y,v0,x2y2); sci(z,pi2,x4iy4,v0,x4iy4); sci(z,pi32,ix4y4,v0,ix4y4); 
                gl(x,z,v0,x2z2); gl(y,z,v0,y2z2); gl(xiy,z,v0,z2); gl(xy,z,v0,x2y2z2);
                break;
            }
        
        case 127: sc(z,pi,v0,v0); sc(z,pi2,v0,v0); sc(z,pi32,v0,v0); sc(y,pi,x4,y2); sc(x,pi,y4,x2); sc(xy,pi,v0,x2y2); sc(xiy,pi,y2,v0);
                  in(v0); gl(x,y,v0,v0); sci(z,pi2,v0,v0,v0); sci(z,pi32,v0,v0,v0); gl(x,z,y4,x2); gl(y,z,x4,y2); gl(xiy,z,x2,v0);
                  gl(xy,z,v0,x2y2); break; //P4/mbm
        
        case 128: sc(z,pi,v0,v0); sc(z,pi2,v0,v0); sc(z,pi32,v0,v0); sc(y,pi,x4z4,y2); sc(x,pi,y4z4,x2); sc(xy,pi,z4,x2y2);
                  sc(xiy,pi,y2z4,v0); in(v0); gl(x,y,v0,v0); sci(z,pi2,v0,v0,v0); sci(z,pi32,v0,v0,v0); gl(x,z,y4,x2z2);
                  gl(y,z,x4,y2z2); gl(xiy,z,x2,z2); gl(xy,z,v0,x2y2z2); break; //P4/mnc
        
        case 129: sc(z,pi,v0,v0); sc(z,pi2,y2,v0); sc(z,pi32,x2,v0); sc(y,pi,x4,y2); sc(x,pi,y4,x2); sc(xy,pi,v0,v0); sc(xiy,pi,v0,v0);
                  in(x4y4); gl(x,y,v0,x2y2); sci(z,pi2,v0,v0,v0); sci(z,pi32,v0,v0,v0); gl(x,z,v0,v0); gl(y,z,v0,v0); gl(xiy,z,x2,v0);
                  gl(xy,z,v0,x2y2); //P4/nmm
                  if (origin_one()) break;
                  else {
                kill_oc1();
                sc(z,pi,x4y4,v0); sc(z,pi2,x4y4,v0); sc(z,pi32,x4y4,v0); sc(y,pi,v0,y2); sc(x,pi,v0,x2); sc(xy,pi,v0,x2y2); 
                sc(xiy,pi,v0,v0); in(v0); gl(x,y,v0,x2y2); sci(z,pi2,x4iy4,v0,x4iy4); sci(z,pi32,ix4y4,v0,ix4y4);
                gl(x,z,y4,v0); gl(y,z,x4,v0); gl(xiy,z,x2,v0); gl(xy,z,v0,v0); //P4/nmm
                break;
            }
        
        case 130: sc(z,pi,v0,v0); sc(z,pi2,y2,v0); sc(z,pi32,x2,v0); sc(y,pi,x4z4,y2); sc(x,pi,y4z4,x2); sc(xy,pi,z4,v0); sc(xiy,pi,z4,v0);
                  in(x4y4); gl(x,y,v0,x2y2); sci(z,pi2,v0,v0,v0); sci(z,pi32,v0,v0,v0); gl(x,z,v0,z2); gl(y,z,v0,z2); gl(xiy,z,x2,z2);
                  gl(xy,z,v0,x2y2z2); //P4/ncc
                  if (origin_one()) break;
                  else {
                kill_oc1();
                sc(z,pi,x4y4,v0); sc(z,pi2,x4y4,v0); sc(z,pi32,x4y4,v0); sc(y,pi,z4,y2); sc(x,pi,z4,x2); sc(xy,pi,z4,x2y2); 
                sc(xiy,pi,z4,v0); in(v0); gl(x,y,v0,x2y2); sci(z,pi2,x4iy4,v0,x4iy4); sci(z,pi32,ix4y4,v0,ix4y4); gl(x,z,y4,z2); 
                gl(y,z,x4,z2); gl(xiy,z,x2,z2); gl(xy,z,v0,z2); //P4/ncc
                break;
            }
        
        case 131: sc(z,pi,v0,v0); sc(z,pi2,v0,z2); sc(z,pi32,v0,z2); sc(y,pi,v0,v0); sc(x,pi,v0,v0); sc(xy,pi,z4,v0); sc(xiy,pi,z4,v0);
                  in(v0); gl(x,y,v0,v0); sci(z,pi2,v0,v0,z4); sci(z,pi32,v0,v0,z4); gl(x,z,v0,v0); gl(y,z,v0,v0); gl(xiy,z,v0,z2);
                  gl(xy,z,v0,z2); break; //P4_2/mmc
        
        case 132: sc(z,pi,v0,v0); sc(z,pi2,v0,z2); sc(z,pi32,v0,z2); sc(y,pi,z4,v0); sc(x,pi,z4,v0); sc(xy,pi,v0,v0); sc(xiy,pi,v0,v0);
                  in(v0); gl(x,y,v0,v0); sci(z,pi2,v0,v0,z4); sci(z,pi32,v0,v0,z4); gl(x,z,v0,z2); gl(y,z,v0,z2); gl(xiy,z,v0,v0);
                  gl(xy,z,v0,v0); break; //P4_2/mcm
        
        case 133: sc(z,pi,v0,v0); sc(z,pi2,y2,z2); sc(z,pi32,x2,z2); sc(y,pi,z4,v0); sc(x,pi,z4,v0); sc(xy,pi,v0,x2y2); sc(xiy,pi,y2,v0);
                  in(x4y4z4); gl(x,y,z4,x2y2); sci(z,pi2,v0,v0,v0); sci(z,pi32,v0,v0,v0); gl(x,z,y4,x2); gl(y,z,x4,y2); gl(xiy,z,v0,z2);
                  gl(xy,z,v0,z2); //P4_2/nbc
                  if (origin_one()) break;
                  else {
                kill_oc1();
                sc(z,pi,x4y4,v0); sc(z,pi2,x4y4,z2); sc(z,pi32,x4y4,z2); sc(y,pi,x4,v0); sc(x,pi,y4,v0); sc(xy,pi,z4,v0); 
                sc(xiy,pi,y2z4,v0); in(v0); gl(x,y,v0,x2y2); sci(z,pi2,x4iy4,v0,x4iy4z4); sci(z,pi32,ix4y4,v0,ix4y4z4); 
                gl(x,z,v0,x2); gl(y,z,v0,y2); gl(xiy,z,v0,z2); gl(xy,z,v0,x2y2z2); //P4_2/nbc
                break;
            }
        
        case 134: sc(z,pi,v0,v0); sc(z,pi2,y2,z2); sc(z,pi32,x2,z2); sc(y,pi,v0,v0); sc(x,pi,v0,v0); sc(xy,pi,z4,x2y2); sc(xiy,pi,y2z4,v0);
                  in(x4y4z4); gl(x,y,z4,x2y2); sci(z,pi2,v0,v0,v0); sci(z,pi32,v0,v0,v0); gl(x,z,y4,x2z2); gl(y,z,x4,y2z2); gl(xiy,z,v0,v0);
                  gl(xy,z,v0,v0); //P4_2/nnm
                  if (origin_one()) break;
                  else {
                kill_oc1();
                sc(z,pi,x4y4,v0); sc(z,pi2,x4y4,z2); sc(z,pi32,x4y4,z2); sc(y,pi,x4z4,v0); sc(x,pi,y4z4,v0); sc(xy,pi,v0,v0); 
                sc(xiy,pi,y2,v0); in(v0); gl(x,y,v0,x2y2); sci(z,pi2,x4iy4,v0,x4iy4z4); sci(z,pi32,ix4y4,v0,ix4y4z4); 
                gl(x,z,v0,x2z2); gl(y,z,v0,y2z2); gl(xiy,z,v0,v0); gl(xy,z,v0,x2y2); //P4_2/nnm
                break;
            }
        
        case 135: sc(z,pi,v0,v0); sc(z,pi2,v0,z2); sc(z,pi32,v0,z2); sc(y,pi,x4,y2); sc(x,pi,y4,x2); sc(xy,pi,z4,x2y2); sc(xiy,pi,y2z4,v0);
                  in(v0); gl(x,y,v0,v0); sci(z,pi2,v0,v0,z4); sci(z,pi32,v0,v0,z4); gl(x,z,y4,x2); gl(y,z,x4,y2); gl(xiy,z,x2,z2);
                  gl(xy,z,v0,x2y2z2); break; //P4_2/mbc
        
        case 136: sc(z,pi,v0,v0); sc(z,pi2,y2,z2); sc(z,pi32,x2,z2); sc(y,pi,x4z4,y2); sc(x,pi,y4z4,x2); sc(xy,pi,v0,v0); sc(xiy,pi,v0,v0);
                  in(v0); gl(x,y,v0,v0); sci(z,pi2,x2,v0,x2z4); sci(z,pi32,y2,v0,y2z4); gl(x,z,y4,x2z2); gl(y,z,x4,y2z2); gl(xiy,z,v0,v0);
                  gl(xy,z,v0,v0); break; //P4_2/mnm
        
        case 137: sc(z,pi,v0,v0); sc(z,pi2,y2,z2); sc(z,pi32,x2,z2); sc(y,pi,x4z4,y2); sc(x,pi,y4z4,x2); sc(xy,pi,v0,v0); sc(xiy,pi,v0,v0);
                  in(x4y4z4); gl(x,y,z4,x2y2); sci(z,pi2,v0,v0,v0); sci(z,pi32,v0,v0,v0); gl(x,z,v0,v0); gl(y,z,v0,v0); gl(xiy,z,x2,z2);
                  gl(xy,z,v0,x2y2z2); //P4_2/nmc
                  if (origin_one()) break;
                  else {
                kill_oc1();
                sc(z,pi,x4y4,v0); sc(z,pi2,x4y4,z2); sc(z,pi32,x4y4,z2); sc(y,pi,v0,y2); sc(x,pi,v0,x2); sc(xy,pi,z4,x2y2); 
                sc(xiy,pi,z4,v0); in(v0); gl(x,y,v0,x2y2); sci(z,pi2,x4iy4,v0,x4iy4z4); sci(z,pi32,ix4y4,v0,ix4y4z4); 
                gl(x,z,y4,v0); gl(y,z,x4,v0); gl(xiy,z,x2,z2); gl(xy,z,v0,z2); //P4_2/nmc
                break;
            }
        
        case 138: sc(z,pi,v0,v0); sc(z,pi2,y2,z2); sc(z,pi32,x2,z2); sc(y,pi,x4,y2); sc(x,pi,y4,x2); sc(xy,pi,z4,v0); sc(xiy,pi,z4,v0);
                  in(x4y4z4); gl(x,y,z4,x2y2); sci(z,pi2,v0,v0,v0); sci(z,pi32,v0,v0,v0); gl(x,z,v0,z2); gl(y,z,v0,z2); gl(xiy,z,x2,v0);
                  gl(xy,z,v0,x2y2); //P4_2/ncm
                  if (origin_one()) break;
                  else {
                kill_oc1();
                sc(z,pi,x4y4,v0); sc(z,pi2,x4y4,z2); sc(z,pi32,x4y4,z2); sc(y,pi,z4,y2); sc(x,pi,z4,x2); sc(xy,pi,v0,x2y2); 
                sc(xiy,pi,v0,v0); in(v0); gl(x,y,v0,x2y2); sci(z,pi2,x4iy4,v0,x4iy4z4); sci(z,pi32,ix4y4,v0,ix4y4z4); 
                gl(x,z,y4,z2); gl(y,z,x4,z2); gl(xiy,z,x2,v0); gl(xy,z,v0,v0); //P4_2/ncm
                break;
            }
        
        case 139: sc(z,pi,v0,v0); sc(z,pi2,v0,v0); sc(z,pi32,v0,v0); sc(y,pi,v0,v0); sc(x,pi,v0,v0); sc(xy,pi,v0,v0); sc(xiy,pi,v0,v0);
                  in(v0); gl(x,y,v0,v0); sci(z,pi2,v0,v0,v0); sci(z,pi32,v0,v0,v0); gl(x,z,v0,v0); gl(y,z,v0,v0); gl(xiy,z,v0,v0);
                  gl(xy,z,v0,v0); tr(x2y2z2); sc(z,pi,x4y4,z2); sc(z,pi2,y2,z2); sc(z,pi32,x2,z2); sc(y,pi,x4z4,y2); sc(x,pi,y4z4,x2);
                  sc(xy,pi,z4,x2y2); sc(xiy,pi,y2z4,v0); in(x4y4z4); gl(x,y,z4,x2y2); sci(z,pi2,x2,v0,x2z4); sci(z,pi32,y2,v0,y2z4);
                  gl(x,z,y4,x2z2); gl(y,z,x4,y2z2); gl(xiy,z,x2,z2); gl(xy,z,v0,x2y2z2); break; //I4/mmm
        
        case 140: sc(z,pi,v0,v0); sc(z,pi2,v0,v0); sc(z,pi32,v0,v0); sc(y,pi,z4,v0); sc(x,pi,z4,v0); sc(xy,pi,z4,v0); sc(xiy,pi,z4,v0);
                  in(v0); gl(x,y,v0,v0); sci(z,pi2,v0,v0,v0); sci(z,pi32,v0,v0,v0); gl(x,z,v0,z2); gl(y,z,v0,z2); gl(xiy,z,v0,z2);
                  gl(xy,z,v0,z2); tr(x2y2z2); sc(z,pi,x4y4,z2); sc(z,pi2,y2,z2); sc(z,pi32,x2,z2); sc(y,pi,x4,y2); sc(x,pi,y4,x2);
                  sc(xy,pi,v0,x2y2); sc(xiy,pi,y2,v0); in(x4y4z4); gl(x,y,z4,x2y2); sci(z,pi2,x2,v0,x2z4); sci(z,pi32,y2,v0,y2z4);
                  gl(x,z,y4,x2); gl(y,z,x4,y2); gl(xiy,z,x2,v0); gl(xy,z,v0,x2y2); break; //I4/mcm
        
        case 141: sc(z,pi,x4y4,z2); sc(z,pi2,ix4y4,z4); sc(z,pi32,x4iy4,z34); sc(y,pi,x4z38,v0); sc(x,pi,y4z8,v0); sc(xy,pi,z4,x2y2);
                  sc(xiy,pi,v0,v0); in(y4z8); gl(x,y,z38,x2); sci(z,pi2,v0,v0,v0); sci(z,pi32,y2,v0,y2z4); gl(x,z,y4,x2z2);
                  gl(y,z,v0,v0); gl(xiy,z,x4,x4iy4z34); gl(xy,z,ix4,x4y4z4); tr(x2y2z2); sc(z,pi,v0,v0); sc(z,pi2,x4y4,z34);
                  sc(z,pi32,x4y4,z4); sc(y,pi,z8,y2); sc(x,pi,z38,x2); sc(xy,pi,v0,v0); sc(xiy,pi,y2z4,v0); in(x4z38);
                  gl(x,y,z8,y2); sci(z,pi2,x2,v0,x2z4); sci(z,pi32,v0,v0,v0); gl(x,z,v0,v0); gl(y,z,x4,y2z2);
                  gl(xiy,z,x4,ix4y4z4); gl(xy,z,x4,x4y4z34); //I4_1/amd
                  if (origin_one()) break;
                  else {
                kill_oc1();
                sc(z,pi,x4,z2); sc(z,pi2,ix4y2,z4); sc(z,pi32,x4,z34); sc(y,pi,x4z4,v0); sc(x,pi,v0,v0); sc(xy,pi,y4z8,x2y2);
                sc(xiy,pi,y4z38,v0); in(v0); gl(x,y,z4,x2); sci(z,pi2,x2iy4,v0,x2iy4z38); sci(z,pi32,y34,v0,y34z8); gl(x,z,v0,x2z2);
                gl(y,z,v0,v0); gl(xiy,z,x2,x4iy4z34); gl(xy,z,v0,x34y34z4); tr(x2y2z2); sc(z,pi,y4,v0); sc(z,pi2,x4y2,z34);
                sc(z,pi32,x34,z4); sc(y,pi,v0,y2); sc(x,pi,y4z4,x2); sc(xy,pi,iy4z38,x2y2); sc(xiy,pi,y34z8,v0); in(x4y4z4);
                gl(x,y,v0,y2); sci(z,pi2,x2y4,v0,x2y4z8); sci(z,pi32,y4,v0,y4z38); gl(x,z,y4,v0); gl(y,z,x4,y2z2);
                gl(xiy,z,x2,ix4y4z4); gl(xy,z,v0,x4y4z34); //I4_1/amd
                break;
            }
        
        case 142: sc(z,pi,x4y4,z2); sc(z,pi2,ix4y4,z4); sc(z,pi32,x4iy4,z34); sc(y,pi,x4z8,v0); sc(x,pi,y4z38,v0); sc(xy,pi,v0,x2y2);
                  sc(xiy,pi,z4,v0); in(y4z8); gl(x,y,z38,x2); sci(z,pi2,v0,v0,v0); sci(z,pi32,y2,v0,y2z4); gl(x,z,y4,x2);
                  gl(y,z,v0,z2); gl(xiy,z,x4,x4iy4z4); gl(xy,z,ix4,x4y4z34); tr(x2y2z2); sc(z,pi,v0,v0); sc(z,pi2,x4y4,z34);
                  sc(z,pi32,x4y4,z4); sc(y,pi,z38,y2); sc(x,pi,z8,x2); sc(xy,pi,z4,v0); sc(xiy,pi,y2,v0); in(x4z38);
                  gl(x,y,z8,y2); sci(z,pi2,x2,v0,x2z4); sci(z,pi32,v0,v0,v0); gl(x,z,v0,z2); gl(y,z,x4,y2);
                  gl(xiy,z,x4,ix4y4z34); gl(xy,z,x4,x4y4z4); //I4_1/acd
                  if (origin_one()) break;
                  else {
                kill_oc1();
                sc(z,pi,x4,z2); sc(z,pi2,ix4y2,z4); sc(z,pi32,x4,z34); sc(y,pi,x4,v0); sc(x,pi,z4,v0); sc(xy,pi,y4z38,x2y2);
                sc(xiy,pi,y4z8,v0); in(v0); gl(x,y,z4,x2); sci(z,pi2,x2iy4,v0,x2iy4z38); sci(z,pi32,y34,v0,y34z8); gl(x,z,v0,x2);
                gl(y,z,v0,z2); gl(xiy,z,x2,x4iy4z4); gl(xy,z,v0,x34y34z34); tr(x2y2z2); sc(z,pi,y4,v0); sc(z,pi2,x4y2,z34);
                sc(z,pi32,x34,z4); sc(y,pi,z4,y2); sc(x,pi,y4,x2); sc(xy,pi,iy4z8,x2y2); sc(xiy,pi,y34z38,v0); in(x4y4z4);
                gl(x,y,v0,y2); sci(z,pi2,x2y4,v0,x2y4z8); sci(z,pi32,y4,v0,y4z38); gl(x,z,y4,z2); gl(y,z,x4,y2);
                gl(xiy,z,x2,ix4y4z34); gl(xy,z,v0,x4y4z4); //I4_1/acd
                break;
            }
        
        //!--------------------------------------------------------------------------------------------------------
        //!Trigonal:
        case 143: sc(z,pi23,v0,v0); sc(z,pi43,v0,v0); break; //P3
        
        case 144: sc(z,pi23,v0,z3); sc(z,pi43,v0,z23); break; //P3_1
        
        case 145: sc(z,pi23,v0,z23); sc(z,pi43,v0,z3); break; //P3_2
        
        //Hexagonale Achsen:
        case 146: if (l->crysin->rhombo) { //Rhombohedrale Achsen
                sc(xyz,pi23,v0,v0); sc(xyz,pi43,v0,v0); break; //R3
            } else { //Hexagonale Achsen
                sc(z,pi23,v0,v0); sc(z,pi43,v0,v0); tr(x23y3z3); sc(z,pi23,x3y3,z3); sc(z,pi43,x3,z3);
                tr(x3y23z23); sc(z,pi23,y3,z23); sc(z,pi43,x3y3,z23); break; //R3
            }
        
        case 147: sc(z,pi23,v0,v0); sc(z,pi43,v0,v0); in(v0); sci(z,pi23,v0,v0,v0); sci(z,pi43,v0,v0,v0); break; //P3_inv
        
        case 148: if (l->crysin->rhombo) { //Rhombohedrale Achsen
                sc(xyz,pi23,v0,v0); sc(xyz,pi43,v0,v0); in(v0); sci(xyz,pi23,v0,v0,v0);
                sci(xyz,pi43,v0,v0,v0); break; //R3_inv
            } else {
                sc(z,pi23,v0,v0); sc(z,pi43,v0,v0); in(v0); sci(z,pi23,v0,v0,v0); sci(z,pi43,v0,v0,v0);
                tr(x23y3z3); sc(z,pi23,x3y3,z3); sc(z,pi43,x3,z3); in(x3y6z6); sci(z,pi23,x3iy3,v0,x3iy3z6);
                sci(z,pi43,x3y23,v0,x3y23z6); tr(x3y23z23); sc(z,pi23,y3,z23); sc(z,pi43,x3y3,z23);
                in(x6y3z3); sci(z,pi23,x23y3,v0,x23y3z3); sci(z,pi43,ix3y3,v0,ix3y3z3); break; //R3_inv
            }
        
        case 149: sc(z,pi23,v0,v0); sc(z,pi43,v0,v0); sc(xiy,pi,v0,v0); sc(xdy,pi,v0,v0); sc(dxy,pi,v0,v0); break; //P312
        
        case 150: sc(z,pi23,v0,v0); sc(z,pi43,v0,v0); sc(xy,pi,v0,v0); sc(x,pi,v0,v0); sc(y,pi,v0,v0); break; //P321
        
        case 151: sc(z,pi23,v0,z3); sc(z,pi43,v0,z23); sc(xiy,pi,z3,v0); sc(xdy,pi,z6,v0); sc(dxy,pi,v0,v0); break; //P3_112
        
        case 152: sc(z,pi23,v0,z3); sc(z,pi43,v0,z23); sc(xy,pi,v0,v0); sc(x,pi,z3,v0); sc(y,pi,z6,v0); break; //P3_121
        
        case 153: sc(z,pi23,v0,z23); sc(z,pi43,v0,z3); sc(xiy,pi,z6,v0); sc(xdy,pi,z3,v0); sc(dxy,pi,v0,v0); break; //P3_212
        
        case 154: sc(z,pi23,v0,z23); sc(z,pi43,v0,z3); sc(xy,pi,v0,v0); sc(x,pi,z6,v0); sc(y,pi,z3,v0); break; //P3_221
        
        case 155: if (l->crysin->rhombo) { //Rhombohedrale Achsen
                sc(xyz,pi23,v0,v0); sc(xyz,pi43,v0,v0); sc(xiy,pi,v0,v0); sc(yiz,pi,v0,v0); sc(ixz,pi,v0,v0); break;
            } else {
                sc(z,pi23,v0,v0); sc(z,pi43,v0,v0); sc(xy,pi,v0,v0); sc(x,pi,v0,v0); sc(y,pi,v0,v0);
                tr(x23y3z3); sc(z,pi23,x3y3,z3); sc(z,pi43,x3,z3); sc(xy,pi,iy6z6,x2y2); sc(x,pi,y6z6,x2);
                sc(y,pi,x3z6,v0); tr(x3y23z23); sc(z,pi23,y3,z23); sc(z,pi43,x3y3,z23);
                sc(xy,pi,y6z3,x2y2); sc(x,pi,y3z3,v0); sc(y,pi,x6z3,y2); break; //R32
            }
        
        case 156: sc(z,pi23,v0,v0); sc(z,pi43,v0,v0); gl(xy,z,v0,v0); gl(xdy,z,v0,v0); gl(dxy,z,v0,v0); break; //P3m1
        
        case 157: sc(z,pi23,v0,v0); sc(z,pi43,v0,v0); gl(xy,z,v0,v0); gl(x,z,v0,v0); gl(y,z,v0,v0); break; //P31m
        
        case 158: sc(z,pi23,v0,v0); sc(z,pi43,v0,v0); gl(xy,z,v0,z2); gl(xdy,z,v0,z2); gl(dxy,z,v0,z2); break; //P3c1
        
        case 159: sc(z,pi23,v0,v0); sc(z,pi43,v0,v0); gl(xy,z,v0,z2); gl(x,z,v0,z2); gl(y,z,v0,z2); break; //P31c
        
        case 160: if (l->crysin->rhombo) { //Rhombohedrale Achsen
                sc(xyz,pi23,v0,v0); sc(xyz,pi43,v0,v0); gl(xiy,z,v0,v0); gl(x,yz,v0,v0);
                gl(xz,y,v0,v0); break;
            } else {
                sc(z,pi23,v0,v0); sc(z,pi43,v0,v0); gl(xiy,z,v0,v0); gl(xdy,z,v0,v0); gl(dxy,z,v0,v0);
                tr(x23y3z3); sc(z,pi23,x3y3,z3); sc(z,pi43,x3,z3); gl(xiy,z,x2,x6iy6z3); gl(xdy,z,iy2,x6y3z3);
                gl(dxy,z,v0,x23y3z3); tr(x3y23z23); sc(z,pi23,y3,z23); sc(z,pi43,x3y3,z23); gl(xiy,z,x2,ix6y6z23);
                gl(xdy,z,v0,x3y23z23); gl(dxy,z,ix2,x3y6z23); break; //R3m
            }
        
        case 161: if (l->crysin->rhombo) { //Rhombohedrale Achsen
                sc(xyz,pi23,v0,v0); sc(xyz,pi43,v0,v0); gl(xy,z,v0,x2y2z2); gl(x,yz,v0,x2y2z2);
                gl(xz,y,v0,x2y2z2); break;
            } else {
                sc(z,pi23,v0,v0); sc(z,pi43,v0,v0); gl(xiy,z,v0,z2); gl(xdy,z,v0,z2); gl(dxy,z,v0,z2);
                tr(x23y3z3); sc(z,pi23,x3y3,z3); sc(z,pi43,x3,z3); gl(xiy,z,x2,x6iy6z56); gl(xdy,z,iy2,x6y3z56);
                gl(dxy,z,v0,x23y3z56); tr(x3y23z23); sc(z,pi23,y3,z23); sc(z,pi43,x3y3,z23); gl(xiy,z,x2,ix6y6z6);
                gl(xdy,z,v0,x3y23z6); gl(dxy,z,ix2,x3y6z6); break; //R3c
            }
        
        case 162: sc(z,pi23,v0,v0); sc(z,pi43,v0,v0); sc(xiy,pi,v0,v0); sc(xdy,pi,v0,v0); sc(dxy,pi,v0,v0);
                  in(v0); sci(z,pi23,v0,v0,v0); sci(z,pi43,v0,v0,v0); gl(xy,z,v0,v0); gl(x,z,v0,v0);
                  gl(y,z,v0,v0); break; //P3_inv1m
        
        case 163: sc(z,pi23,v0,v0); sc(z,pi43,v0,v0); sc(xiy,pi,z4,v0); sc(xdy,pi,z4,v0); sc(dxy,pi,z4,v0);
                  in(v0); sci(z,pi23,v0,v0,v0); sci(z,pi43,v0,v0,v0); gl(xy,z,v0,z2); gl(x,z,v0,z2);
                  gl(y,z,v0,z2); break; //P3_inv1c
        
        case 164: sc(z,pi23,v0,v0); sc(z,pi43,v0,v0); sc(xy,pi,v0,v0); sc(x,pi,v0,v0); sc(y,pi,v0,v0);
                  in(v0); sci(z,pi23,v0,v0,v0); sci(z,pi43,v0,v0,v0); gl(xiy,z,v0,v0); gl(xdy,z,v0,v0);
                  gl(dxy,z,v0,v0); break; //P3_invm1
        
        case 165: sc(z,pi23,v0,v0); sc(z,pi43,v0,v0); sc(xy,pi,z4,v0); sc(x,pi,z4,v0); sc(y,pi,z4,v0);
                  in(v0); sci(z,pi23,v0,v0,v0); sci(z,pi43,v0,v0,v0); gl(xiy,z,v0,z2); gl(xdy,z,v0,z2);
                  gl(dxy,z,v0,z2); break; //P3_invc1
        
        case 166: if (l->crysin->rhombo) { //Rhombohedrale Achsen
                sc(xyz,pi23,v0,v0); sc(xyz,pi43,v0,v0); sc(xiy,pi,v0,v0); sc(yiz,pi,v0,v0);
                sc(ixz,pi,v0,v0); in(v0); sci(xyz,pi23,v0,v0,v0); sci(xyz,pi43,v0,v0,v0);
                gl(xy,z,v0,v0); gl(x,yz,v0,v0); gl(xz,y,v0,v0); break;
            } else {
                sc(z,pi23,v0,v0); sc(z,pi43,v0,v0); sc(xy,pi,v0,v0); sc(x,pi,v0,v0); sc(y,pi,v0,v0);
                in(v0); sci(z,pi23,v0,v0,v0); sci(z,pi43,v0,v0,v0); gl(xiy,z,v0,v0); gl(xdy,z,v0,v0);
                gl(dxy,z,v0,v0); tr(x23y3z3); sc(z,pi23,x3y3,z3); sc(z,pi43,x3,z3); sc(xy,pi,iy6z6,x2y2);
                sc(x,pi,y6z6,x2); sc(y,pi,x3z6,v0); in(x3y6z6); sci(z,pi23,x3iy3,v0,x3iy3z6); sci(z,pi43,x3y23,v0,x3y23z6);
                gl(xiy,z,x2,x6iy6z3); gl(xdy,z,iy2,x6y3z3); gl(dxy,z,v0,x23y3z3); tr(x3y23z23); sc(z,pi23,y3,z23);
                sc(z,pi43,x3y3,z23); sc(xy,pi,y6z3,x2y2); sc(x,pi,y3z3,v0); sc(y,pi,x6z3,y2); in(x6y3z3);
                sci(z,pi23,x23y3,v0,x23y3z3); sci(z,pi43,ix3y3,v0,ix3y3z3); gl(xiy,z,x2,ix6y6z23);
                gl(xdy,z,v0,x3y23z23); gl(dxy,z,ix2,x3y6z23); break; //R3_invm
            }
        
        case 167: if (l->crysin->rhombo) { //Rhombohedrale Achsen
                sc(xyz,pi23,v0,v0); sc(xyz,pi43,v0,v0); sc(xiy,pi,y2z4,v0); sc(yiz,pi,x4z2,v0);
                sc(ixz,pi,x2y4,v0); in(v0); sci(xyz,pi23,v0,v0,v0); sci(xyz,pi43,v0,v0,v0);
                gl(xy,z,v0,x2y2z2); gl(x,yz,v0,x2y2z2); gl(xz,y,v0,x2y2z2); break;
            } else {
                sc(z,pi23,v0,v0); sc(z,pi43,v0,v0); sc(xy,pi,z4,v0); sc(x,pi,z4,v0); sc(y,pi,z4,v0);
                in(v0); sci(z,pi23,v0,v0,v0); sci(z,pi43,v0,v0,v0); gl(xiy,z,v0,z2); gl(xdy,z,v0,z2);
                gl(dxy,z,v0,z2); tr(x23y3z3); sc(z,pi23,x3y3,z3); sc(z,pi43,x3,z3); sc(xy,pi,iy6z512,x2y2);
                sc(x,pi,y6z512,x2); sc(y,pi,x3z512,v0); in(x3y6z6); sci(z,pi23,x3iy3,v0,x3iy3z6); sci(z,pi43,x3y23,v0,x3y23z6);
                gl(xiy,z,x2,x6iy6z56); gl(xdy,z,iy2,x6y3z56); gl(dxy,z,v0,x23y3z56); tr(x3y23z23); sc(z,pi23,y3,z23);
                sc(z,pi43,x3y3,z23); sc(xy,pi,y6z12,x2y2); sc(x,pi,y3z12,v0); sc(y,pi,x6z12,y2); in(x6y3z3);
                sci(z,pi23,x23y3,v0,x23y3z3); sci(z,pi43,ix3y3,v0,ix3y3z3); gl(xiy,z,x2,ix6y6z6);
                gl(xdy,z,v0,x3y23z6); gl(dxy,z,ix2,x3y6z6); break; //R3_invc
            }
        
        //!--------------------------------------------------------------------------------------------------------
        //!Hexagonal:
        case 168: sc(z,pi23,v0,v0); sc(z,pi43,v0,v0); sc(z,pi,v0,v0); sc(z,pi53,v0,v0); sc(z,pi3,v0,v0); break; //P6
        
        case 169: sc(z,pi23,v0,z3); sc(z,pi43,v0,z23); sc(z,pi,v0,z2); sc(z,pi53,v0,z56); sc(z,pi3,v0,z6); break; //P6_1
        
        case 170: sc(z,pi23,v0,z23); sc(z,pi43,v0,z3); sc(z,pi,v0,z2); sc(z,pi53,v0,z6); sc(z,pi3,v0,z56); break; //P6_5
        
        case 171: sc(z,pi23,v0,z23); sc(z,pi43,v0,z3); sc(z,pi,v0,v0); sc(z,pi53,v0,z23); sc(z,pi3,v0,z3); break; //P6_2
        
        case 172: sc(z,pi23,v0,z3); sc(z,pi43,v0,z23); sc(z,pi,v0,v0); sc(z,pi53,v0,z3); sc(z,pi3,v0,z23); break; //P6_4
        
        case 173: sc(z,pi23,v0,v0); sc(z,pi43,v0,v0); sc(z,pi,v0,z2); sc(z,pi53,v0,z2); sc(z,pi3,v0,z2); break; //P6_3
        
        case 174: sc(z,pi23,v0,v0); sc(z,pi43,v0,v0); gl(x,y,v0,v0); sci(z,pi53,v0,v0,v0); sci(z,pi3,v0,v0,v0); break; //P6_inv
        
        case 175: sc(z,pi23,v0,v0); sc(z,pi43,v0,v0); sc(z,pi,v0,v0); sc(z,pi53,v0,v0); sc(z,pi3,v0,v0); in(v0);
                  sci(z,pi23,v0,v0,v0); sci(z,pi43,v0,v0,v0); gl(x,y,v0,v0); sci(z,pi53,v0,v0,v0); sci(z,pi3,v0,v0,v0); break; //P6/m
        
        case 176: sc(z,pi23,v0,v0); sc(z,pi43,v0,v0); sc(z,pi,v0,z2); sc(z,pi53,v0,z2); sc(z,pi3,v0,z2); in(v0);
                  sci(z,pi23,v0,v0,v0); sci(z,pi43,v0,v0,v0); gl(x,y,z4,v0); sci(z,pi53,v0,v0,z4); sci(z,pi3,v0,v0,z4); break; //P6_3/m
        
        case 177: sc(z,pi23,v0,v0); sc(z,pi43,v0,v0); sc(z,pi,v0,v0); sc(z,pi53,v0,v0); sc(z,pi3,v0,v0); sc(xy,pi,v0,v0);
                  sc(x,pi,v0,v0); sc(y,pi,v0,v0); sc(xiy,pi,v0,v0); sc(xdy,pi,v0,v0); sc(dxy,pi,v0,v0); break; //P622
        
        case 178: sc(z,pi23,v0,z3); sc(z,pi43,v0,z23); sc(z,pi,v0,z2); sc(z,pi53,v0,z56); sc(z,pi3,v0,z6); sc(xy,pi,z6,v0);
                  sc(x,pi,v0,v0); sc(y,pi,z3,v0); sc(xiy,pi,z512,v0); sc(xdy,pi,z4,v0); sc(dxy,pi,z12,v0); break; //P6_122
        
        case 179: sc(z,pi23,v0,z23); sc(z,pi43,v0,z3); sc(z,pi,v0,z2); sc(z,pi53,v0,z6); sc(z,pi3,v0,z56); sc(xy,pi,z3,v0);
                  sc(x,pi,v0,v0); sc(y,pi,z6,v0); sc(xiy,pi,z12,v0); sc(xdy,pi,z4,v0); sc(dxy,pi,z512,v0); break; //P6_522
        
        case 180: sc(z,pi23,v0,z23); sc(z,pi43,v0,z3); sc(z,pi,v0,v0); sc(z,pi53,v0,z23); sc(z,pi3,v0,z3); sc(xy,pi,z3,v0);
                  sc(x,pi,v0,v0); sc(y,pi,z6,v0); sc(xiy,pi,z3,v0); sc(xdy,pi,v0,v0); sc(dxy,pi,z6,v0); break; //P6_222
        
        case 181: sc(z,pi23,v0,z3); sc(z,pi43,v0,z23); sc(z,pi,v0,v0); sc(z,pi53,v0,z3); sc(z,pi3,v0,z23); sc(xy,pi,z6,v0);
                  sc(x,pi,v0,v0); sc(y,pi,z3,v0); sc(xiy,pi,z6,v0); sc(xdy,pi,v0,v0); sc(dxy,pi,z3,v0); break; //P6_422
        
        case 182: sc(z,pi23,v0,v0); sc(z,pi43,v0,v0); sc(z,pi,v0,z2); sc(z,pi53,v0,z2); sc(z,pi3,v0,z2); sc(xy,pi,v0,v0);
                  sc(x,pi,v0,v0); sc(y,pi,v0,v0); sc(xiy,pi,z4,v0); sc(xdy,pi,z4,v0); sc(dxy,pi,z4,v0); break; //P6_322
        
        case 183: sc(z,pi23,v0,v0); sc(z,pi43,v0,v0); sc(z,pi,v0,v0); sc(z,pi53,v0,v0); sc(z,pi3,v0,v0); 
                  gl(xiy,z,v0,v0); gl(xdy,z,v0,v0); gl(dxy,z,v0,v0); gl(xy,z,v0,v0); gl(x,z,v0,v0); gl(y,z,v0,v0); break;
        
        case 184: sc(z,pi23,v0,v0); sc(z,pi43,v0,v0); sc(z,pi,v0,v0); sc(z,pi53,v0,v0); sc(z,pi3,v0,v0);
                  gl(xiy,z,v0,z2); gl(xdy,z,v0,z2); gl(dxy,z,v0,z2); gl(xy,z,v0,z2); gl(x,z,v0,z2); gl(y,z,v0,z2); break;
        
        case 185: sc(z,pi23,v0,v0); sc(z,pi43,v0,v0); sc(z,pi,v0,z2); sc(z,pi53,v0,z2); sc(z,pi3,v0,z2);
                  gl(xiy,z,v0,z2); gl(xdy,z,v0,z2); gl(dxy,z,v0,z2); gl(xy,z,v0,v0); gl(x,z,v0,v0); gl(y,z,v0,v0); break;
        
        case 186: sc(z,pi23,v0,v0); sc(z,pi43,v0,v0); sc(z,pi,v0,z2); sc(z,pi53,v0,z2); sc(z,pi3,v0,z2);
                  gl(xiy,z,v0,v0); gl(xdy,z,v0,v0); gl(dxy,z,v0,v0); gl(xy,z,v0,z2); gl(x,z,v0,z2); gl(y,z,v0,z2); break;
        
        case 187: sc(z,pi23,v0,v0); sc(z,pi43,v0,v0); gl(x,y,v0,v0); sci(z,pi53,v0,v0,v0); sci(z,pi3,v0,v0,v0);
                  gl(xiy,z,v0,v0); gl(xdy,z,v0,v0); gl(dxy,z,v0,v0); sc(xiy,pi,v0,v0); sc(xdy,pi,v0,v0); sc(dxy,pi,v0,v0); break;
        
        case 188: sc(z,pi23,v0,v0); sc(z,pi43,v0,v0); gl(x,y,z4,v0); sci(z,pi53,v0,v0,z4); sci(z,pi3,v0,v0,z4);
                  gl(xiy,z,v0,z2); gl(xdy,z,v0,z2); gl(dxy,z,v0,z2); sc(xiy,pi,v0,v0); sc(xdy,pi,v0,v0); sc(dxy,pi,v0,v0); break;
        
        case 189: sc(z,pi23,v0,v0); sc(z,pi43,v0,v0); gl(x,y,v0,v0); sci(z,pi53,v0,v0,v0); sci(z,pi3,v0,v0,v0);
                  sc(xy,pi,v0,v0); sc(x,pi,v0,v0); sc(y,pi,v0,v0); gl(xy,z,v0,v0); gl(x,z,v0,v0); gl(y,z,v0,v0); break;
        
        case 190: sc(z,pi23,v0,v0); sc(z,pi43,v0,v0); gl(x,y,z4,v0); sci(z,pi53,v0,v0,z4); sci(z,pi3,v0,v0,z4);
                  sc(xy,pi,v0,v0); sc(x,pi,v0,v0); sc(y,pi,v0,v0); gl(xy,z,v0,z2); gl(x,z,v0,z2); gl(y,z,v0,z2); break;
        
        case 191: sc(z,pi23,v0,v0); sc(z,pi43,v0,v0); sc(z,pi,v0,v0); sc(z,pi53,v0,v0); sc(z,pi3,v0,v0);
                  sc(xy,pi,v0,v0); sc(x,pi,v0,v0); sc(y,pi,v0,v0); sc(xiy,pi,v0,v0); sc(xdy,pi,v0,v0); sc(dxy,pi,v0,v0);
                  in(v0); sci(z,pi23,v0,v0,v0); sci(z,pi43,v0,v0,v0); gl(x,y,v0,v0); sci(z,pi53,v0,v0,v0); sci(z,pi3,v0,v0,v0);
                  gl(xiy,z,v0,v0); gl(xdy,z,v0,v0); gl(dxy,z,v0,v0); gl(xy,z,v0,v0); gl(x,z,v0,v0); gl(y,z,v0,v0); break;
        
        case 192: sc(z,pi23,v0,v0); sc(z,pi43,v0,v0); sc(z,pi,v0,v0); sc(z,pi53,v0,v0); sc(z,pi3,v0,v0);
                  sc(xy,pi,z4,v0); sc(x,pi,z4,v0); sc(y,pi,z4,v0); sc(xiy,pi,z4,v0); sc(xdy,pi,z4,v0); sc(dxy,pi,z4,v0);
                  in(v0); sci(z,pi23,v0,v0,v0); sci(z,pi43,v0,v0,v0); gl(x,y,v0,v0); sci(z,pi53,v0,v0,v0); sci(z,pi3,v0,v0,v0);
                  gl(xiy,z,v0,z2); gl(xdy,z,v0,z2); gl(dxy,z,v0,z2); gl(xy,z,v0,z2); gl(x,z,v0,z2); gl(y,z,v0,z2); break;
        
        case 193: sc(z,pi23,v0,v0); sc(z,pi43,v0,v0); sc(z,pi,v0,z2); sc(z,pi53,v0,z2); sc(z,pi3,v0,z2);
                  sc(xy,pi,z4,v0); sc(x,pi,z4,v0); sc(y,pi,z4,v0); sc(xiy,pi,v0,v0); sc(xdy,pi,v0,v0); sc(dxy,pi,v0,v0);
                  in(v0); sci(z,pi23,v0,v0,v0); sci(z,pi43,v0,v0,v0); gl(x,y,z4,v0); sci(z,pi53,v0,v0,z4); sci(z,pi3,v0,v0,z4);
                  gl(xiy,z,v0,z2); gl(xdy,z,v0,z2); gl(dxy,z,v0,z2); gl(xy,z,v0,v0); gl(x,z,v0,v0); gl(y,z,v0,v0); break;
        
        case 194: sc(z,pi23,v0,v0); sc(z,pi43,v0,v0); sc(z,pi,v0,z2); sc(z,pi53,v0,z2); sc(z,pi3,v0,z2);
                  sc(xy,pi,v0,v0); sc(x,pi,v0,v0); sc(y,pi,v0,v0); sc(xiy,pi,z4,v0); sc(xdy,pi,z4,v0); sc(dxy,pi,z4,v0);
                  in(v0); sci(z,pi23,v0,v0,v0); sci(z,pi43,v0,v0,v0); gl(x,y,z4,v0); sci(z,pi53,v0,v0,z4); sci(z,pi3,v0,v0,z4);
                  gl(xiy,z,v0,v0); gl(xdy,z,v0,v0); gl(dxy,z,v0,v0); gl(xy,z,v0,z2); gl(x,z,v0,z2); gl(y,z,v0,z2);break;
        
        //!--------------------------------------------------------------------------------------------------------
        //!Kubisch:
        case 195: sc(z,pi,v0,v0); sc(y,pi,v0,v0); sc(x,pi,v0,v0); sc(xyz,pi23,v0,v0); sc(ixyiz,pi23,v0,v0);
                  sc(xiyiz,pi23,v0,v0); sc(ixiyz,pi23,v0,v0); sc(xyz,pi43,v0,v0); sc(xiyiz,pi43,v0,v0); sc(ixiyz,pi43,v0,v0);
                  sc(ixyiz,pi43,v0,v0); break;
        
        
        case 196: sc(z,pi,v0,v0); sc(y,pi,v0,v0); sc(x,pi,v0,v0); sc(xyz,pi23,v0,v0); sc(ixyiz,pi23,v0,v0);
                  sc(xiyiz,pi23,v0,v0); sc(ixiyz,pi23,v0,v0); sc(xyz,pi43,v0,v0); sc(xiyiz,pi43,v0,v0); sc(ixiyz,pi43,v0,v0);
                  sc(ixyiz,pi43,v0,v0); tr(y2z2); sc(z,pi,y4,z2); sc(y,pi,z4,y2); sc(x,pi,y4z4,v0); sc(xyz,pi23,ix3iy6,x3y3z3);
                  sc(ixyiz,pi23,y2,v0); sc(xiyiz,pi23,x3iy6,ix3y3z3); sc(ixiyz,pi23,y2,v0); sc(xyz,pi43,ix6y6,x3y3z3);
                  sc(xiyiz,pi43,x6y6,ix3y3z3); sc(ixiyz,pi43,x2y2,v0); sc(ixyiz,pi43,ix2y2,v0); tr(x2z2); sc(z,pi,x4,z2);
                  sc(y,pi,x4z4,v0); sc(x,pi,z4,x2); sc(xyz,pi23,x6iy6,x3y3z3); sc(ixyiz,pi23,x6y6,x3iy3z3); sc(xiyiz,pi23,x2iy2,v0);
                  sc(ixiyz,pi23,x2y2,v0); sc(xyz,pi43,ix6iy3,x3y3z3); sc(xiyiz,pi43,x2,v0); sc(ixiyz,pi43,x2,v0);
                  sc(ixyiz,pi43,ix6y3,x3iy3z3); tr(x2y2); sc(z,pi,x4y4,v0); sc(y,pi,x4,y2); sc(x,pi,y4,x2); sc(xyz,pi23,x6y3,x3y3z3);
                  sc(ixyiz,pi23,x2,v0); sc(xiyiz,pi23,x2,v0); sc(ixiyz,pi23,x6y3,x3y3iz3); sc(xyz,pi43,x3y6,x3y3z3);
                  sc(xiyiz,pi43,y2,v0); sc(ixiyz,pi43,x3y6,x3y3iz3); sc(ixyiz,pi43,y2,v0); break;
        
        case 197: sc(z,pi,v0,v0); sc(y,pi,v0,v0); sc(x,pi,v0,v0); sc(xyz,pi23,v0,v0); sc(ixyiz,pi23,v0,v0);
                  sc(xiyiz,pi23,v0,v0); sc(ixiyz,pi23,v0,v0); sc(xyz,pi43,v0,v0); sc(xiyiz,pi43,v0,v0); sc(ixiyz,pi43,v0,v0);
                  sc(ixyiz,pi43,v0,v0); tr(x2y2z2); sc(z,pi,x4y4,z2); sc(y,pi,x4z4,y2); sc(x,pi,y4z4,x2); sc(xyz,pi23,v0,x2y2z2);
                  sc(ixyiz,pi23,x3y3,x6iy6z6); sc(xiyiz,pi23,x23iy3,ix6y6z6); sc(ixiyz,pi23,x3y23,x6y6iz6); sc(xyz,pi43,v0,x2y2z2);
                  sc(xiyiz,pi43,x3y3,ix6y6z6); sc(ixiyz,pi43,x23y3,x6y6iz6); sc(ixyiz,pi43,ix3y23,x6iy6z6); break;
        
        
        case 198: sc(z,pi,x4,z2); sc(y,pi,z4,y2); sc(x,pi,y4,x2); sc(xyz,pi23,v0,v0); sc(ixyiz,pi23,x2,v0); sc(xiyiz,pi23,x2iy2,v0);
                  sc(ixiyz,pi23,y2,v0); sc(xyz,pi43,v0,v0); sc(xiyiz,pi43,x6y6,ix3y3z3); sc(ixiyz,pi43,x3y6,x3y3iz3); 
                  sc(ixyiz,pi43,ix6y3,x3iy3z3); break;
        
        case 199: sc(z,pi,x4,z2); sc(y,pi,z4,y2); sc(x,pi,y4,x2); sc(xyz,pi23,v0,v0); sc(ixyiz,pi23,x2,v0); sc(xiyiz,pi23,x2iy2,v0);
                  sc(ixiyz,pi23,y2,v0); sc(xyz,pi43,v0,v0); sc(xiyiz,pi43,x6y6,ix3y3z3); sc(ixiyz,pi43,x3y6,x3y3iz3); 
                  sc(ixyiz,pi43,ix6y3,x3iy3z3); tr(x2y2z2); sc(z,pi,y4,v0); sc(y,pi,x4,v0); sc(x,pi,z4,v0); sc(xyz,pi23,v0,x2y2z2); 
                  sc(ixyiz,pi23,ix6y3,x6iy6z6); sc(xiyiz,pi23,x6y6,ix6y6z6); sc(ixiyz,pi23,x3y6,x6y6iz6); sc(xyz,pi43,v0,x2y2z2); 
                  sc(xiyiz,pi43,x6y6,x6iy6iz6); sc(ixiyz,pi43,x3y6,ix6iy6z6); sc(ixyiz,pi43,ix6y3,ix6y6iz6); break; 
        
        case 200: sc(z,pi,v0,v0); sc(y,pi,v0,v0); sc(x,pi,v0,v0); sc(xyz,pi23,v0,v0); sc(ixyiz,pi23,v0,v0);
                  sc(xiyiz,pi23,v0,v0); sc(ixiyz,pi23,v0,v0); sc(xyz,pi43,v0,v0); sc(xiyiz,pi43,v0,v0); sc(ixiyz,pi43,v0,v0);
                  sc(ixyiz,pi43,v0,v0); in(v0); gl(x,y,v0,v0); gl(x,z,v0,v0); gl(y,z,v0,v0); sci(xyz,pi23,v0,v0,v0); 
                  sci(ixyiz,pi23,v0,v0,v0); sci(xiyiz,pi23,v0,v0,v0); sci(ixiyz,pi23,v0,v0,v0); sci(xyz,pi43,v0,v0,v0);
                  sci(xiyiz,pi43,v0,v0,v0); sci(ixiyz,pi43,v0,v0,v0); sci(ixyiz,pi43,v0,v0,v0); break;
        
        case 201: sc(z,pi,v0,v0); sc(y,pi,v0,v0); sc(x,pi,v0,v0); sc(xyz,pi23,v0,v0); sc(ixyiz,pi23,v0,v0);
                  sc(xiyiz,pi23,v0,v0); sc(ixiyz,pi23,v0,v0); sc(xyz,pi43,v0,v0); sc(xiyiz,pi43,v0,v0); sc(ixiyz,pi43,v0,v0);
                  sc(ixyiz,pi43,v0,v0); in(x4y4z4); gl(x,y,z4,x2y2); gl(x,z,y4,x2z2); gl(y,z,x4,y2z2); sci(xyz,pi23,v0,v0,x4y4z4);
                  sci(ixyiz,pi23,ixy,v0,ix4y4z34);//xquer mit +1 siehe bis incl. 204
                  sci(xiyiz,pi23,y,v0,x4y34iz4); sci(ixiyz,pi23,x,v0,x34iy4z4); sci(xyz,pi43,v0,v0,x4y4z4); sci(xiyiz,pi43,xiy,v0,x4iy4z34);
                  sci(ixiyz,pi43,y,v0,ix4y34z4); sci(ixyiz,pi43,x,v0,x34y4iz4);
                      if (origin_one()) break;
                  else {
                kill_oc1();
                sc(z,pi,x4y4,v0); sc(y,pi,x4z4,v0); sc(x,pi,y4z4,v0); sc(xyz,pi23,v0,v0); sc(ixyiz,pi23,y2,v0);
                sc(xiyiz,pi23,x2,v0); sc(ixiyz,pi23,x2y2,v0); sc(xyz,pi43,v0,v0); sc(xiyiz,pi43,x2,v0); sc(ixiyz,pi43,x2y2,v0);
                sc(ixyiz,pi43,y2,v0); in(v0); gl(x,y,v0,x2y2); gl(x,z,v0,x2z2); gl(y,z,v0,y2z2); sci(xyz,pi23,v0,v0,v0);
                sci(ixyiz,pi23,ixy2,v0,ix2z2); sci(xiyiz,pi23,ix2y,v0,y2iz2); sci(ixiyz,pi23,x2iy2,v0,x2iy2);
                sci(xyz,pi43,v0,v0,v0); sci(xiyiz,pi43,x2iy,v0,iy2z2); sci(ixiyz,pi43,ix2y2,v0,ix2y2);
                sci(ixyiz,pi43,xiy2,v0,x2iz2);
                break;
            }
        
        case 202: sc(z,pi,v0,v0); sc(y,pi,v0,v0); sc(x,pi,v0,v0); sc(xyz,pi23,v0,v0); sc(ixyiz,pi23,v0,v0);
                  sc(xiyiz,pi23,v0,v0); sc(ixiyz,pi23,v0,v0); sc(xyz,pi43,v0,v0); sc(xiyiz,pi43,v0,v0); sc(ixiyz,pi43,v0,v0);
                  sc(ixyiz,pi43,v0,v0); in(v0); gl(x,y,v0,v0); gl(x,z,v0,v0); gl(y,z,v0,v0); sci(xyz,pi23,v0,v0,v0); 
                  sci(ixyiz,pi23,v0,v0,v0); sci(xiyiz,pi23,v0,v0,v0); sci(ixiyz,pi23,v0,v0,v0); sci(xyz,pi43,v0,v0,v0);
                  sci(xiyiz,pi43,v0,v0,v0); sci(ixiyz,pi43,v0,v0,v0); sci(ixyiz,pi43,v0,v0,v0); tr(y2z2); sc(z,pi,y4,z2);
                  sc(y,pi,z4,y2); sc(x,pi,y4z4,v0); sc(xyz,pi23,ix3iy6,x3y3z3); sc(ixyiz,pi23,y2,v0); sc(xiyiz,pi23,x3iy6,ix3y3z3);
                  sc(ixiyz,pi23,y2,v0); sc(xyz,pi43,ix6y6,x3y3z3); sc(xiyiz,pi43,x6y6,ix3y3z3); sc(ixiyz,pi43,x2y2,v0); 
                  sc(ixyiz,pi43,ix2y2,v0); in(y4z4); gl(x,y,z4,y2); gl(x,z,y4,z2); gl(y,z,v0,y2z2); sci(xyz,pi23,y2,v0,y2); 
                  sci(ixyiz,pi23,ixy2,v0,ix2z2); sci(xiyiz,pi23,y2,v0,y2); sci(ixiyz,pi23,xy2,v0,x2z2); sci(xyz,pi43,ix2iy2,v0,z2); 
                  sci(xiyiz,pi43,x2iy2,v0,z2); sci(ixiyz,pi43,ix2y2,v0,ix2y2); sci(ixyiz,pi43,x2y2,v0,x2y2); tr(x2z2); sc(z,pi,x4,z2); 
                  sc(y,pi,x4y4,v0); sc(x,pi,z4,x2); sc(xyz,pi23,x6iy6,x3y3z3); sc(ixyiz,pi23,x6y6,x3iy3z3); sc(xiyiz,pi23,x2iy2,v0); 
                  sc(ixiyz,pi23,x2y2,v0); sc(xyz,pi43,ix6iy3,x3y3z3); sc(xiyiz,pi43,x2,v0); sc(ixiyz,pi43,x2,v0); 
                  sc(ixyiz,pi43,ix6y3,x3iy3z3); in(x4z4); gl(x,y,z4,x2); gl(x,z,v0,x2z2); gl(y,z,x4,z2); sci(xyz,pi23,ix2iy2,v0,z2); 
                  sci(ixyiz,pi23,ix2y2,v0,z2); sci(xiyiz,pi23,x2y2,v0,x2y2); sci(ixiyz,pi23,x2iy2,v0,x2iy2); sci(xyz,pi43,x2,v0,x2); 
                  sci(xiyiz,pi43,x2iy,v0,iy2z2); sci(ixiyz,pi43,x2y,v0,y2z2); sci(ixyiz,pi43,x2,v0,x2); tr(x2y2); sc(z,pi,x4y4,v0); 
                  sc(y,pi,x4,y2); sc(x,pi,y4,x2); sc(xyz,pi23,x6y3,x3y3z3); sc(ixyiz,pi23,x2,v0); sc(xiyiz,pi23,x2,v0); 
                  sc(ixiyz,pi23,x6y3,x3y3iz3); sc(xyz,pi43,x3y6,x3y3z3); sc(xiyiz,pi43,y2,v0); sc(ixiyz,pi43,x3y6,x3y3iz3); 
                  sc(ixyiz,pi43,y2,v0); in(x4y4); gl(x,y,v0,x2y2); gl(x,z,y4,x2); gl(y,z,x4,y2); sci(xyz,pi23,x2,v0,x2); 
                  sci(ixyiz,pi23,ix2y,v0,y2z2); sci(xiyiz,pi23,ix2y,v0,y2iz2); sci(ixiyz,pi23,x2,v0,x2); sci(xyz,pi43,y2,v0,y2); 
                  sci(xiyiz,pi43,xiy2,v0,x2z2); sci(ixiyz,pi43,y2,v0,y2); sci(ixyiz,pi43,xiy2,v0,x2iz2); break;
        
        case 203: sc(z,pi,v0,v0); sc(y,pi,v0,v0); sc(x,pi,v0,v0); sc(xyz,pi23,v0,v0); sc(ixyiz,pi23,v0,v0);
                  sc(xiyiz,pi23,v0,v0); sc(ixiyz,pi23,v0,v0); sc(xyz,pi43,v0,v0); sc(xiyiz,pi43,v0,v0); sc(ixiyz,pi43,v0,v0);
                  sc(ixyiz,pi43,v0,v0); in(x8y8z8); gl(x,y,z8,x4y4); gl(x,z,y8,x4z4); gl(y,z,x8,y4z4); sci(xyz,pi23,v0,v0,x8y8z8);
                  sci(ixyiz,pi23,ix2y2,v0,ix8y8z38); sci(xiyiz,pi23,y2,v0,x8y38iz8); sci(ixiyz,pi23,x2,v0,x38iy8z8);
                  sci(xyz,pi43,v0,v0,x8y8z8); sci(xiyiz,pi43,x2iy2,v0,x8iy8z38); sci(ixiyz,pi43,y2,v0,ix8y38z8); 
                  sci(ixyiz,pi43,x2,v0,x38y8iz8); tr(y2z2); sc(z,pi,y4,z2); sc(y,pi,z4,y2); sc(x,pi,y4z4,v0); sc(xyz,pi23,ix3iy6,x3y3z3); 
                  sc(ixyiz,pi23,y2,v0); sc(xiyiz,pi23,x3iy6,ix3y3z3); sc(ixiyz,pi23,y2,v0); sc(xyz,pi43,ix6y6,x3y3z3); 
                  sc(xiyiz,pi43,x6y6,ix3y3z3); sc(ixiyz,pi43,x2y2,v0); sc(ixyiz,pi43,ix2y2,v0); in(x8y38z38); gl(x,y,z38,x4y34);
                  gl(x,z,y38,x4z34); gl(y,z,x8,y34z34); sci(xyz,pi23,y2,v0,x8y58z8); sci(ixyiz,pi23,ix32y,v0,ix58y8z78); 
                  sci(xiyiz,pi23,y,v0,x8y78iz8); sci(ixiyz,pi23,x32y2,v0,x78iy8z58); sci(xyz,pi43,ix2iy2,v0,x8y8z58); 
                  sci(xiyiz,pi43,xiy,v0,x8iy8z78); sci(ixiyz,pi43,ix2y,v0,ix58y78z8); sci(ixyiz,pi43,xy2,v0,x78y58iz8);
                  tr(x2z2); sc(z,pi,x4,z2); sc(y,pi,x4z4,v0); sc(x,pi,z4,x2); sc(xyz,pi23,x6iy6,x3y3z3); sc(ixyiz,pi23,x6y6,x3iy3z3); 
                  sc(xiyiz,pi23,x2iy2,v0); sc(ixiyz,pi23,x2y2,v0); sc(xyz,pi43,ix6iy3,x3y3z3); sc(xiyiz,pi43,x2,v0); sc(ixiyz,pi43,x2,v0); 
                  sc(ixyiz,pi43,ix6y3,x3iy3z3); in(x38y8z38); gl(x,y,z38,x34y4); gl(x,z,y8,x34z34); gl(y,z,x38,y4z34); 
                  sci(xyz,pi23,ix2iy2,v0,x8y8z58); sci(ixyiz,pi23,ixy,v0,ix8y8z78); sci(xiyiz,pi23,x2y,v0,x58y78iz8); 
                  sci(ixiyz,pi23,xiy2,v0,x78iy58z8); sci(xyz,pi43,x2,v0,x58y8z8); sci(xiyiz,pi43,xiy32,v0,x8iy58z78);
                  sci(ixiyz,pi43,x2y32,v0,ix8y78z58); sci(ixyiz,pi43,x,v0,x78y8iz8); tr(x2y2); sc(z,pi,x4y4,v0); sc(y,pi,x4,y2); 
                  sc(x,pi,y4,x2); sc(xyz,pi23,x6y3,x3y3z3); sc(ixyiz,pi23,x2,v0); sc(xiyiz,pi23,x2,v0); sc(ixiyz,pi23,x6y3,x3y3iz3); 
                  sc(xyz,pi43,x3y6,x3y3z3); sc(xiyiz,pi43,y2,v0); sc(ixiyz,pi43,x3y6,x3y3iz3); sc(ixyiz,pi43,y2,v0); in(x38y38z8); 
                  gl(x,y,z8,x34y34); gl(x,z,y38,x34z4); gl(y,z,x38,y34z4); sci(xyz,pi23,x2,v0,x58y8z8); sci(ixyiz,pi23,ixy32,v0,ix8y58z78); 
                  sci(xiyiz,pi23,ix2y32,v0,x8y78iz58); sci(ixiyz,pi23,x,v0,x78iy8z8); sci(xyz,pi43,y2,v0,x8y58z8); 
                  sci(xiyiz,pi43,x32iy,v0,x58iy8z78); sci(ixiyz,pi43,y,v0,ix8y78z8); sci(ixyiz,pi43,x32iy2,v0,x78y8iz58);
                  if (origin_one()) break;
                  else {
                kill_oc1();
                sc(z,pi,x8y8,v0); sc(y,pi,x8z8,v0); sc(x,pi,y8z8,v0); sc(xyz,pi23,v0,v0); sc(ixyiz,pi23,y4,v0);
                sc(xiyiz,pi23,x4,v0); sc(ixiyz,pi23,x4y4,v0); sc(xyz,pi43,v0,v0); sc(xiyiz,pi43,x4,v0); sc(ixiyz,pi43,x4y4,v0);
                sc(ixyiz,pi43,y4,v0); in(v0); gl(x,y,v0,x34y34); gl(x,z,v0,x34z34); gl(y,z,v0,y34z34); sci(xyz,pi23,v0,v0,v0);
                sci(ixyiz,pi23,ix32y34,v0,ix34z34); sci(xiyiz,pi23,ix34y32,v0,y34iz34); sci(ixiyz,pi23,x34iy34,v0,x34iy34);
                sci(xyz,pi43,v0,v0,v0); sci(xiyiz,pi43,x34iy32,v0,iy34z34); sci(ixiyz,pi43,ix34y34,v0,ix34y34);
                sci(ixyiz,pi43,x32iy34,v0,x34iz34);
                tr(y2z2); sc(z,pi,x8y38,z2); sc(y,pi,x8z38,y2); sc(x,pi,y38z38,v0); 
                sc(xyz,pi23,ix3iy6,x3y3z3); sc(ixyiz,pi23,y34,v0); sc(xiyiz,pi23,x712iy6,ix3y3z3); sc(ixiyz,pi23,x4y34,v0); 
                sc(xyz,pi43,ix6y6,x3y3z3); sc(xiyiz,pi43,x512y6,ix3y3z3); sc(ixiyz,pi43,x34y34,v0); sc(ixyiz,pi43,ix2y34,v0);
                in(y4z4); gl(x,y,z4,x34y4); gl(x,z,y4,x34z4); gl(y,z,v0,y4z4); sci(xyz,pi23,y2,v0,y2);
                sci(ixyiz,pi23,ix2y4,v0,ix4z4); sci(xiyiz,pi23,x4y,v0,x2y34iz4); sci(ixiyz,pi23,x34iy4,v0,x34iy4);
                sci(xyz,pi43,ix2iy2,v0,z2); sci(xiyiz,pi43,x54iy,v0,x2iy4z34); sci(ixiyz,pi43,ix4y4,v0,ix4y4);
                sci(ixyiz,pi43,ixy4,v0,x34iz4);
                tr(x2z2); sc(z,pi,x38y8,z2); sc(y,pi,x38z38,v0); sc(x,pi,y8z38,x2); 
                sc(xyz,pi23,x6iy6,x3y3z3); sc(ixyiz,pi23,x6y512,x3iy3z3); sc(xiyiz,pi23,x34iy2,v0); sc(ixiyz,pi23,x34y34,v0);
                sc(xyz,pi43,ix6iy3,x3y3z3); sc(xiyiz,pi43,x34,v0); sc(ixiyz,pi43,x34y4,v0); sc(ixyiz,pi43,ix6y712,x3iy3z3); 
                in(x4z4); gl(x,y,z4,x4y34); gl(x,z,v0,x4z4); gl(y,z,x4,y34z4); sci(xyz,pi23,ix2iy2,v0,z2);
                sci(ixyiz,pi23,ixy54,v0,ix4y2z34); sci(xiyiz,pi23,ix4y,v0,y34iz4); sci(ixiyz,pi23,x4iy4,v0,x4iy4); 
                sci(xyz,pi43,x2,v0,x2); sci(xiyiz,pi43,x4iy2,v0,iy4z4); sci(ixiyz,pi43,ix4y34,v0,ix4y34);
                sci(ixyiz,pi43,xy4,v0,x34y2iz4);
                tr(x2y2); sc(z,pi,x34y34,v0); sc(y,pi,x38z8,y2); sc(x,pi,y38z8,x2);
                sc(xyz,pi23,x6y3,x3y3z3); sc(ixyiz,pi23,x2y4,v0); sc(xiyiz,pi23,x34,v0); sc(ixiyz,pi23,x512y712,x3y3iz3);
                sc(xyz,pi43,x3y6,x3y3z3); sc(xiyiz,pi43,x4y2,v0); sc(ixiyz,pi43,x712y512,x3y3iz3); sc(ixyiz,pi43,y34,v0); 
                in(x4y4); gl(x,y,v0,x4y4); gl(x,z,y4,x4z34); gl(y,z,x4,y4z34); sci(xyz,pi23,x2,v0,x2);
                sci(ixyiz,pi23,ixy34,v0,ix4z34); sci(xiyiz,pi23,ix4y2,v0,y4iz4); sci(ixiyz,pi23,x54y4,v0,x34iy4z2);
                sci(xyz,pi43,y2,v0,y2); sci(xiyiz,pi43,x34iy,v0,iy4z34); sci(ixiyz,pi43,x4y54,v0,ix4y34z2);
                sci(ixyiz,pi43,x2iy4,v0,x4iz4); break;
            }
        
        case 204: sc(z,pi,v0,v0); sc(y,pi,v0,v0); sc(x,pi,v0,v0); sc(xyz,pi23,v0,v0); sc(ixyiz,pi23,v0,v0);
                  sc(xiyiz,pi23,v0,v0); sc(ixiyz,pi23,v0,v0); sc(xyz,pi43,v0,v0); sc(xiyiz,pi43,v0,v0); sc(ixiyz,pi43,v0,v0);
                  sc(ixyiz,pi43,v0,v0); in(v0); gl(x,y,v0,v0); gl(x,z,v0,v0); gl(y,z,v0,v0); sci(xyz,pi23,v0,v0,v0); 
                  sci(ixyiz,pi23,v0,v0,v0); sci(xiyiz,pi23,v0,v0,v0); sci(ixiyz,pi23,v0,v0,v0); sci(xyz,pi43,v0,v0,v0);
                  sci(xiyiz,pi43,v0,v0,v0); sci(ixiyz,pi43,v0,v0,v0); sci(ixyiz,pi43,v0,v0,v0); tr(x2y2z2); sc(z,pi,x4y4,z2);
                  sc(y,pi,x4z4,y2); sc(x,pi,y4z4,x2); sc(xyz,pi23,v0,x2y2z2); sc(ixyiz,pi23,x3y3,x6iy6z6); sc(xiyiz,pi23,x23iy3,ix6y6z6);
                  sc(ixiyz,pi23,x3y23,x6y6z6); sc(xyz,pi43,v0,x2y2z2); sc(xiyiz,pi43,x3y3,ix6y6z6); sc(ixiyz,pi43,x23y3,x6y6iz6);
                  sc(ixyiz,pi43,ix3y23,x6iy6z6); in(x4y4z4); gl(x,y,z4,x2y2); gl(x,z,y4,x2z2); gl(y,z,x4,y2z2); sci(xyz,pi23,v0,v0,x4y4z4);
                  sci(ixyiz,pi23,ixy,v0,ix4y4z34); sci(xiyiz,pi23,y,v0,x4y34iz4); sci(ixiyz,pi23,x,v0,x34iy4z4); sci(xyz,pi43,v0,v0,x4y4z4);
                  sci(xiyiz,pi43,xiy,v0,x4iy4z34); sci(ixiyz,pi43,y,v0,ix4y34z4); sci(ixyiz,pi43,x,v0,x34y4iz4);break;
        
        case 205: sc(z,pi,x4,z2); sc(y,pi,z4,y2); sc(x,pi,y4,x2); sc(xyz,pi23,v0,v0); sc(ixyiz,pi23,x2,v0); sc(xiyiz,pi23,x2iy2,v0);
                  sc(ixiyz,pi23,y2,v0); sc(xyz,pi43,v0,v0); sc(xiyiz,pi43,x6y6,ix3y3z3); sc(ixiyz,pi43,x3y6,x3y3iz3); 
                  sc(ixyiz,pi43,ix6y3,x3iy3z3); in(v0); gl(x,y,z4,x2); gl(x,z,y4,z2); gl(y,z,x4,y2); sci(xyz,pi23,v0,v0,v0); 
                  sci(ixyiz,pi23,ix2y,v0,y2z2); sci(xiyiz,pi23,x2y2,v0,x2y2); sci(ixiyz,pi23,xy2,v0,x2z2); sci(xyz,pi43,v0,v0,v0); 
                  sci(xiyiz,pi43,x2iy2,v0,z2); sci(ixiyz,pi43,y2,v0,y2); sci(ixyiz,pi43,x2,v0,x2); break;
        
        case 206: sc(z,pi,x4,z2); sc(y,pi,z4,y2); sc(x,pi,y4,x2); sc(xyz,pi23,v0,v0); sc(ixyiz,pi23,x2,v0);
                  sc(xiyiz,pi23,x2iy2,v0); sc(ixiyz,pi23,y2,v0); sc(xyz,pi43,v0,v0); sc(xiyiz,pi43,x6y6,ix3y3z3);
                  sc(ixiyz,pi43,x3y6,x3y3iz3); sc(ixyiz,pi43,ix6y3,x3iy3z3); in(v0); sc(z,pi,y4,v0); sc(y,pi,x4,v0);
                  sc(x,pi,z4,v0); sc(xyz,pi23,v0,x2y2z2); sc(ixyiz,pi23,ix6y3,x6iy6z6); sc(xiyiz,pi23,x6y6,ix6y6z6);
                  sc(ixiyz,pi23,x3y6,x6y6iz6); sc(xyz,pi43,v0,x2y2z2); sc(xiyiz,pi43,x6y6,x6iy6iz6); sc(ixiyz,pi43,x3y6,ix6iy6z6);
                  sc(ixyiz,pi43,ix6y3,ix6y6iz6); in(x4y4z4); gl(x,y,v0,y2); gl(x,z,v0,x2); gl(y,z,v0,z2); sci(xyz,pi23,v0,v0,x4y4z4);
                  sci(ixyiz,pi23,ix2,v0,ix4iy4z4); sci(xiyiz,pi23,ix2y2,v0,ix4y4iz4); sci(ixiyz,pi23,iy2,v0,x4iy4iz4);
                  sci(xyz,pi43,v0,v0,x4y4z4); sci(xiyiz,pi43,x2iy2,v0,x4iy4z4); sci(ixiyz,pi43,y2,v0,ix4y4z4); 
                  sci(ixyiz,pi43,x2,v0,x4y4iz4); break;
        
        case 207: sc(z,pi,v0,v0); sc(y,pi,v0,v0); sc(x,pi,v0,v0); sc(xyz,pi23,v0,v0); sc(ixyiz,pi23,v0,v0);
                  sc(xiyiz,pi23,v0,v0); sc(ixiyz,pi23,v0,v0); sc(xyz,pi43,v0,v0); sc(xiyiz,pi43,v0,v0); sc(ixiyz,pi43,v0,v0);
                  sc(ixyiz,pi43,v0,v0); sc(xy,pi,v0,v0); sc(xiy,pi,v0,v0); sc(z,pi32,v0,v0); sc(z,pi2,v0,v0); sc(x,pi32,v0,v0);
                  sc(yz,pi,v0,v0); sc(yiz,pi,v0,v0); sc(x,pi2,v0,v0); sc(y,pi32,v0,v0); sc(xz,pi,v0,v0); sc(y,pi32,v0,v0);
                  sc(ixz,pi,v0,v0); break;
        
        case 208: sc(z,pi,v0,v0); sc(y,pi,v0,v0); sc(x,pi,v0,v0); sc(xyz,pi23,v0,v0); sc(ixyiz,pi23,v0,v0);
                  sc(xiyiz,pi23,v0,v0); sc(ixiyz,pi23,v0,v0); sc(xyz,pi43,v0,v0); sc(xiyiz,pi43,v0,v0); sc(ixiyz,pi43,v0,v0);
                  sc(ixyiz,pi43,v0,v0); sc(xy,pi,z4,x2y2); sc(xiy,pi,y2z4,v0); sc(z,pi32,x2,z2); sc(z,pi2,y2,z2); sc(x,pi32,y2,x2);
                  sc(yz,pi,x4,y2z2); sc(yiz,pi,x4y2,v0); sc(x,pi2,z2,x2); sc(y,pi2,x2,y2); sc(xz,pi,y4,x2z2); sc(y,pi32,z2,y2);
                  sc(ixz,pi,x2y4,v0); break;
        
        case 209: sc(z,pi,v0,v0); sc(y,pi,v0,v0); sc(x,pi,v0,v0); sc(xyz,pi23,v0,v0); sc(ixyiz,pi23,v0,v0);
                  sc(xiyiz,pi23,v0,v0); sc(ixiyz,pi23,v0,v0); sc(xyz,pi43,v0,v0); sc(xiyiz,pi43,v0,v0); sc(ixiyz,pi43,v0,v0);
                  sc(ixyiz,pi43,v0,v0); sc(xy,pi,v0,v0); sc(xiy,pi,v0,v0); sc(z,pi32,v0,v0); sc(z,pi2,v0,v0); sc(x,pi32,v0,v0);
                  sc(yz,pi,v0,v0); sc(yiz,pi,v0,v0); sc(x,pi2,v0,v0); sc(y,pi32,v0,v0); sc(xz,pi,v0,v0); sc(y,pi32,v0,v0);
                  sc(ixz,pi,v0,v0); tr(y2z2); sc(z,pi,y4,z2); sc(y,pi,z4,y2); sc(x,pi,y4z4,v0); sc(xyz,pi23,ix3iy6,x3y3z3);
                  sc(ixyiz,pi23,y2,v0); sc(xiyiz,pi23,x3iy6,ix3y3z3); sc(ixiyz,pi23,y2,v0); sc(xyz,pi43,ix6y6,x3y3z3);
                  sc(xiyiz,pi43,x6y6,ix3y3z3); sc(ixiyz,pi43,x2y2,v0); sc(ixyiz,pi43,ix2y2,v0); sc(xy,pi,y4z4,x4y4);
                  sc(xiy,pi,y4z4,ix4y4); sc(z,pi32,x4y4,z2); sc(z,pi2,ix4y4,z2); sc(x,pi32,y2,v0); sc(yz,pi,v0,y2z2);
                  sc(yiz,pi,y2,v0); sc(x,pi2,z2,v0); sc(y,pi2,x4z4,y2); sc(xz,pi,ix4y4,x4z4); sc(y,pi32,ix4y4,y2);
                  sc(ixz,pi,x4y4,ix4z4); tr(x2z2); sc(z,pi,x4,z2); sc(y,pi,x4z4,v0); sc(x,pi,z4,x2); sc(xyz,pi23,x6iy6,x3y3z3);
                  sc(ixyiz,pi23,x6y6,x3iy3z3); sc(xiyiz,pi23,x2iy2,v0); sc(ixiyz,pi23,x2y2,v0); sc(xyz,pi43,ix6iy3,x3y3z3);
                  sc(xiyiz,pi43,x2,v0); sc(ixiyz,pi43,x2,v0); sc(ixyiz,pi43,ix6y3,x3iy3z3); sc(xy,pi,iy4z4,x4y4);
                  sc(xiy,pi,y4z4,x4iy4); sc(z,pi32,x4iy4,z2); sc(z,pi2,x4y4,z2); sc(x,pi32,y4z4,x2); sc(yz,pi,x4iy4,x4z4);
                  sc(yiz,pi,x4y4,iy4z4); sc(x,pi2,iy4z4,x2); sc(y,pi2,x2,v0); sc(xz,pi,v0,x2z2); sc(y,pi32,z2,v0);
                  sc(ixz,pi,x2,v0); tr(x2y2); sc(z,pi,x4y4,v0); sc(y,pi,x4,y2); sc(x,pi,y4,x2); sc(xyz,pi23,x6y3,x3y3z3);
                  sc(ixyiz,pi23,x2,v0); sc(xiyiz,pi23,x2,v0); sc(ixiyz,pi23,x6y3,x3y3iz3); sc(xyz,pi43,x3y6,x3y3z3);
                  sc(xiyiz,pi43,y2,v0); sc(ixiyz,pi43,x3y6,x3y3iz3); sc(ixyiz,pi43,y2,v0); sc(xy,pi,v0,x2y2);
                  sc(xiy,pi,y2,v0); sc(z,pi32,x2,v0); sc(z,pi2,y2,v0); sc(x,pi32,y4iz4,x2); sc(yz,pi,x4y4,y4z4);
                  sc(yiz,pi,x4y4,y4iz4); sc(x,pi2,y4z4,x2); sc(y,pi2,x4iy4,y2); sc(xz,pi,x4y4,x4z4); sc(y,pi32,x4z4,y2);
                  sc(ixz,pi,x4y4,x4iz4); break;
        
        case 210: tr(y2z2); sc(z,pi,v0,v0); sc(y,pi,x4z4,v0); sc(x,pi,y4,x2); sc(xyz,pi23,ix3iy6,x3y3z3);
                  sc(ixyiz,pi23,x2,v0); sc(xiyiz,pi23,v0,v0); sc(ixiyz,pi23,x2y2,v0); sc(xyz,pi43,ix6y6,x3y3z3);
                  sc(xiyiz,pi43,x2,v0); sc(ixiyz,pi43,x3y6,x3y3iz3); sc(ixyiz,pi43,v0,v0); sc(xy,pi,z8,x34y34);
                  sc(xiy,pi,y2z38,ix4y4); sc(z,pi32,x4,z4); sc(z,pi2,x4y2,z34); sc(x,pi32,y2iz4,x34); sc(yz,pi,x38iy4,y2z2);
                  sc(yiz,pi,x8y34,v0); sc(x,pi2,z4,x4); sc(y,pi2,x2iz4,y34); sc(xz,pi,y8,x4z4); sc(y,pi32,z34,y4);
                  sc(ixz,pi,x2y38,ix4z4); tr(x2z2); sc(z,pi,x4y4,v0); sc(y,pi,z4,y2); sc(x,pi,v0,v0); sc(xyz,pi23,x6iy6,x3y3z3);
                  sc(ixyiz,pi23,v0,v0); sc(xiyiz,pi23,x2,v0); sc(ixiyz,pi23,y2,v0); sc(xyz,pi43,ix6iy3,x3y3z3); sc(xiyiz,pi43,x6y6,ix3y3z3);
                  sc(ixiyz,pi43,v0,v0); sc(ixyiz,pi43,y2,v0); sc(xy,pi,z8,x4y4); sc(xiy,pi,y2z38,x4iy4); sc(z,pi32,x34,z4);
                  sc(z,pi2,ix4y2,z34); sc(x,pi32,y4,x4); sc(yz,pi,x8,y34z34); sc(yiz,pi,x38y2,iy4z4); sc(x,pi2,y4z2,x34);
                  sc(y,pi2,x4,y4); sc(xz,pi,x4y38,x2z2); sc(y,pi32,ix4z2,y34); sc(ixz,pi,x34y8,v0); tr(x2y2); sc(z,pi,x4,z2);
                  sc(y,pi,v0,v0); sc(x,pi,y4z4,v0); sc(xyz,pi23,x6y3,x3y3z3); sc(ixyiz,pi23,y2,v0); sc(xiyiz,pi23,x2iy2,v0);
                  sc(ixiyz,pi23,v0,v0); sc(xyz,pi43,x3y6,x3y3z3); sc(xiyiz,pi43,v0,v0); sc(ixiyz,pi43,x2y2,v0); 
                  sc(ixyiz,pi43,ix6y3,x3iy3z3); sc(xy,pi,y4z38,x2y2); sc(xiy,pi,y34z8,v0); sc(z,pi32,x2iy4,z34); sc(z,pi2,y4,z4); 
                  sc(x,pi32,y34,x4); sc(yz,pi,x8,y4z4); sc(yiz,pi2,x38y2,y2iz4); sc(x,pi2,iy4z2,x34); sc(y,pi2,x2z4,y34); 
                  sc(xz,pi,y8,x34z34); sc(y,pi32,z4,y4); sc(ixz,pi,x2y38,x4iz4); break;
        
        case 211: sc(z,pi,v0,v0); sc(y,pi,v0,v0); sc(x,pi,v0,v0); sc(xyz,pi23,v0,v0); sc(ixyiz,pi23,v0,v0);
                  sc(xiyiz,pi23,v0,v0); sc(ixiyz,pi23,v0,v0); sc(xyz,pi43,v0,v0); sc(xiyiz,pi43,v0,v0); sc(ixiyz,pi43,v0,v0);
                  sc(ixyiz,pi43,v0,v0); sc(xy,pi,v0,v0); sc(xiy,pi,v0,v0); sc(z,pi32,v0,v0); sc(z,pi2,v0,v0); sc(x,pi32,v0,v0);
                  sc(yz,pi,v0,v0); sc(yiz,pi,v0,v0); sc(x,pi2,v0,v0); sc(y,pi32,v0,v0); sc(xz,pi,v0,v0); sc(y,pi32,v0,v0);
                  sc(ixz,pi,v0,v0);tr(x2y2z2); sc(z,pi,x4y4,z2); sc(y,pi,x4z4,y2); sc(x,pi,y4z4,x2); sc(xyz,pi23,v0,x2y2z2);
                  sc(ixyiz,pi23,x3y3,x6iy6z6); sc(xiyiz,pi23,x23iy3,ix6y6z6); sc(ixiyz,pi23,x3y23,x6y6iz6); sc(xyz,pi43,v0,x2y2z2);
                  sc(xiyiz,pi43,x3y3,ix6y6z6); sc(ixiyz,pi43,x23y3,x6y6iz6); sc(ixyiz,pi43,ix3y23,x6iy6z6); sc(xy,pi,z4,x2y2);
                  sc(yz,pi,x4,y2z2); sc(yiz,pi,x4y2,v0); sc(x,pi2,z2,x2); sc(y,pi2,x2,y2); sc(xz,pi,y4,x2z2); sc(y,pi32,z2,y2);
                  sc(ixz,pi,x2y4,v0); break;
        
        case 212: sc(z,pi,x4,z2); sc(y,pi,z4,y2); sc(x,pi,y4,x2); sc(xyz,pi23,v0,v0); sc(ixyiz,pi23,x2,v0); sc(xiyiz,pi23,x2iy2,v0);
                  sc(ixiyz,pi23,y2,v0); sc(xyz,pi43,v0,v0); sc(xiyiz,pi43,x6y6,ix3y3z3); sc(ixiyz,pi43,x3y6,x3y3iz3); 
                  sc(ixyiz,pi43,ix6y3,x3iy3z3); sc(xy,pi,y4z38,x2y2); sc(xiy,pi,y4z8,v0); sc(z,pi32,x34,z4); sc(z,pi2,x4y2,z34); 
                  sc(x,pi32,y34,x4); sc(yz,pi,x38iy4,y2z2); sc(yiz,pi,x8y4,v0); sc(x,pi2,y4z2,x34); sc(y,pi2,x2z4,y34); 
                  sc(xz,pi,x4y38,x2z2); sc(y,pi32,z34,y4); sc(ixz,pi,x4y8,v0); break;
        
        case 213: sc(z,pi,x4,z2); sc(y,pi,z4,y2); sc(x,pi,y4,x2); sc(xyz,pi23,v0,v0); sc(ixyiz,pi23,x2,v0); sc(xiyiz,pi23,x2iy2,v0);
                  sc(ixiyz,pi23,y2,v0); sc(xyz,pi43,v0,v0); sc(xiyiz,pi43,x6y6,ix3y3z3); sc(ixiyz,pi43,x3y6,x3y3iz3); 
                  sc(ixyiz,pi43,ix6y3,x3iy3z3); sc(xy,pi,iy4z8,x2y2); sc(xiy,pi,y34z38,v0); sc(z,pi32,x4,z34); sc(z,pi2,ix4y2,z2); 
                  sc(x,pi32,y4,x34); sc(yz,pi,x8y4,y2z2); sc(yiz,pi,x38y34,v0); sc(x,pi2,iy4z2,x4); sc(y,pi2,x2iz4,y4); 
                  sc(xz,pi,ix4y8,x2z2); sc(y,pi32,z4,y34); sc(ixz,pi,x34y38,v0);break;
        
        default: if (l->crysin->group < 1 || l->crysin->group > 230) {
                if (verbosity) cerr << c_message<cWARNING>("CRYSTALIZER::build_unit_cell():  --> Spacegroup ") << l->crysin->group
                     << " does not exist" << endl; ret_false = true; break;
            } else {
                if (verbosity) cerr << c_message<cWARNING>("CRYSTALIZER::build_unit_cell():  --> Spacegroup ") << l->crysin->group
                     << " is not yet implemented!" << endl; ret_false = true; break;
            }
    }
    
    //cerr << "unit is ready" << endl;
    
    if (visualize) {
        for (int i=0; i<=sym_ele; ++i) {
            f_out << "cmd.load_cgo(" << sym_name[i] << ",'" << sym_name[i] << "')" << endl;
        }
        f_out.close();
    }
    
    if (ret_false) return false;
    else return true;
    
}

void CRYSTALIZER::duplicate_unit_cell_splitted() {
    vec3d<float> dv;
    float tmp_dst;
    bool take_it;
    for (int x=-3; x<4; ++x) {
        for (int y=-3; y<4; ++y) {
            for (int z=-3; z<4; ++z) {
                if (x==0 && y==0 && z==0) continue;
                dv[0] = x; dv[1] = y; dv[2] = z;
                for (int i=0; i<n_mol; ++i) {
                    l = structure->splitted_ligands[i];
                    stl_ptr<LIGAND> lig(new LIGAND(*l)); lig->main_structure = structure;
                    ostringstream hos; hos << n_ligs; lig->name += hos.str();
                    take_it = false;
                    for (atoms_vec at=lig->atoms.begin(); at!=lig->atoms.end(); ++at) {
                        (*at)->coord *= l->crysin->cart2cryst;
                        (*at)->coord += dv;
                        (*at)->coord *= l->crysin->cryst2cart;
                        if (!take_it) for (int nnn=0; nnn<n_splitted; ++nnn) {
                            for (atoms_vec bt=structure->splitted_ligands[nnn]->atoms.begin();
                                       bt!=structure->splitted_ligands[nnn]->atoms.end(); ++bt) {
                                tmp_dst = get_square_distance((*bt)->coord,(*at)->coord);
                                if (tmp_dst <= smax_d) {take_it = true; break;}
                            }
                            if (take_it) break;
                        }
                    }
                
                    //!=== neu: 22.05.2009 ====================================
                    //! Pruefen auf doppelt besetzte Lage:
                    //! Leider kann es passieren, dass ein Ligand tatsaechlich aus mehreren Molekuelen
                    //! besteht, von denen eventuell nur eines doppelt ist (weil es auf eine Grenzflaeche
                    //! der Einheitszelle liegt). => Deshalb muessen leider alle Atome gegen alle verglichen werden:
                    //! Umsetzung siehe oben!!!
                    if (take_it) {
                        for (atoms_vec at=lig->atoms.begin(); at!=lig->atoms.end(); ++at) {
                            for (ligands_vec vlig=structure->splitted_ligands.begin(); 
                                             vlig!=structure->splitted_ligands.end(); ++vlig) {
                                for (atoms_vec bt=(*vlig)->atoms.begin(); bt!=(*vlig)->atoms.end(); ++bt) {
                                    if (get_square_distance((*bt)->coord,(*at)->coord) < 0.25) {
                                        take_it = false;
                                        break;
                                    }
                                }
                                if (!take_it) break;
                            }
                            if (!take_it) break;
                        }
                    }
                    //!========================================================
        
                    if (take_it) {
                        structure->splitted_ligands.push_back(lig); ++n_ligs;
                    } else {
                        lig.kill();
                    }
                }
            }
        }
    }
}

void CRYSTALIZER::duplicate_unit_cell2_splitted(int number) {
    vec3d<float> dv;
    bool take_it;
    for (int x=-number; x<=number; ++x) {
        for (int y=-number; y<=number; ++y) {
            for (int z=-number; z<=number; ++z) {
                if (x==0 && y==0 && z==0) continue;
                dv[0] = x; dv[1] = y; dv[2] = z;
                for (int i=0; i<n_mol; ++i) {
                    l = structure->splitted_ligands[i];
                    stl_ptr<LIGAND> lig(new LIGAND(*l)); lig->main_structure = structure;
                    ostringstream hos; hos << n_ligs; lig->name += hos.str();
                    for (atoms_vec at=lig->atoms.begin(); at!=lig->atoms.end(); ++at) {
                        (*at)->coord *= l->crysin->cart2cryst;
                        (*at)->coord += dv;
                        (*at)->coord *= l->crysin->cryst2cart;
                    }
                    take_it = true;
                    for (atoms_vec at=lig->atoms.begin(); at!=lig->atoms.end(); ++at) {
                        for (ligands_vec vlig=structure->splitted_ligands.begin(); 
                                    vlig!=structure->splitted_ligands.end(); ++vlig) {
                            for (atoms_vec bt=(*vlig)->atoms.begin(); bt!=(*vlig)->atoms.end(); ++bt) {
                                if (get_square_distance((*bt)->coord,(*at)->coord) < 0.25) {
                                    take_it = false;
                                    break;
                                }
                            }
                            if (!take_it) break;
                        }
                        if (!take_it) break;
                    }
                    
                    if (take_it) {
                        structure->splitted_ligands.push_back(lig); ++n_ligs;
                    } else {
                        lig.kill();
                    }
                }
            }
        }
    }
}

int CRYSTALIZER::build_crystal_splitted(float max_dst) {
    if (l->crysin.zero()) {if (verbosity) cerr << c_message<cERROR>("CRYSTALIZER::build_crystal --> no CRYSIN-entry for ") 
                                << l->name << endl; return false;}
    if (l->crysin->group == 0) {if (verbosity) cerr << c_message<cERROR>("CRYSTALIZER::build_crystal --> no CRYSIN-entry for ") 
                                     << l->name << endl; return false;}
    max_d = max_dst;
    smax_d = max_dst * max_dst;
    
    //! Achtung: split() braucht ein atom typing!!!
    l->get_atom_typing(0,true,"X",false,0,false,false,"X");
    n_splitted = l->split();
    
    //!ligands leer machen:
    for (ligands_vec lig=structure->ligands.begin(); lig!=structure->ligands.end(); ++lig) {
        lig->kill();
    }
    structure->ligands.clear();
    
    //! Kristallbau muss fuer alle gesplitteten durchgefuehrt werden
    n_ligs = 1;
    for (int ilig=0; ilig<n_splitted; ++ilig) {
        l = structure->splitted_ligands[ilig];
        structure->ligands.push_back(structure->splitted_ligands[ilig]);
        if (build_unit_cell()) {
            //! von ligands nach splitted ligands kopieren:
            for (unsigned int i=1; i<structure->ligands.size(); ++i) {
                structure->splitted_ligands.push_back(structure->ligands[i]);
            }
            structure->ligands.clear();
        } else {
            structure->ligands.clear();
            return 0;
        }
    }
    
    n_mol = structure->splitted_ligands.size();
    duplicate_unit_cell_splitted();
    
    return n_splitted;
}

int CRYSTALIZER::build_crystal2_splitted(int number) {
    if (l->crysin.zero()) {if (verbosity) cerr << c_message<cERROR>("CRYSTALIZER::build_crystal --> no CRYSIN-entry for ") 
                                << l->name << endl; return false;}
    if (l->crysin->group == 0) {if (verbosity) cerr << c_message<cERROR>("CRYSTALIZER::build_crystal --> no CRYSIN-entry for ") 
                                     << l->name << endl; return false;}
    
    //! Achtung: split() braucht ein atom typing!!!
    l->get_atom_typing(0,true,"X",false,0,false,false,"X");
    n_splitted = l->split();
    
    //!ligands leer machen:
    for (ligands_vec lig=structure->ligands.begin(); lig!=structure->ligands.end(); ++lig) {
        lig->kill();
    }
    structure->ligands.clear();
    
    //! Kristallbau muss fuer alle gesplitteten durchgefuehrt werden
    n_ligs = 1;
    for (int ilig=0; ilig<n_splitted; ++ilig) {
        l = structure->splitted_ligands[ilig];
        structure->ligands.push_back(structure->splitted_ligands[ilig]);
        if (build_unit_cell()) {
            //! von ligands nach splitted ligands kopieren:
            for (unsigned int i=1; i<structure->ligands.size(); ++i) {
                structure->splitted_ligands.push_back(structure->ligands[i]);
            }
            structure->ligands.clear();
        } else {
            structure->ligands.clear();
            return 0;
        }
    }
    
    n_mol = structure->splitted_ligands.size();
    duplicate_unit_cell2_splitted(number);
    
    return n_splitted;
}
