
//============================================================================
// structure_additional_GN.h -*- C++ -*-; CRYSIN objects
//
// Copyright (C) 2007, 2008, 2009 Gerd Neudert
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
// This library implements objects for CRYSIN representation.
//============================================================================


#ifndef ADDITIONALGN_H
#define ADDITIONALGN_H

#include<stdlib.h>
#include<iostream>
#include<string>
#include<sstream>
#include<vector>
#include<map>
#include"linalg_GN.hpp"
#include"stl_ptr_GN.hpp"
#include"message_GN.hpp"


using namespace std;
using namespace TXT;

class CRYSIN_POSITION;
class CRYSIN;

//==================================================================================================
//Deklaration der Klasse CRYSIN:
//==================================================================================================

class CRYSIN_POSITION {
    public:
        vec3d<float> const_add;
        vec3d<float> x_add;
        vec3d<float> y_add;
        vec3d<float> z_add;
        CRYSIN_POSITION();
        CRYSIN_POSITION(CRYSIN_POSITION &ref);
        ~CRYSIN_POSITION();
};

class CRYSIN {
    public:
        float a; float b; float c;
        float wa; float wb; float wc;
        matrix<float> cryst2cart; //Transformationsmatrix um Kristallkoordinaten in kartesische Koordinaten umzuwandeln
        matrix<float> cart2cryst; //the other way around
        int group;
        int setting;
        bool rhombo;
        vector<stl_ptr<CRYSIN_POSITION> > positions; //!positions aus cif-files geparsed
        CRYSIN();
        CRYSIN(CRYSIN &crys);
        ~CRYSIN();
};

#endif
