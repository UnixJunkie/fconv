
//============================================================================
// structure_additional_GN.cpp -*- C++ -*-; CRYSIN objects
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


#include"structure_additional_GN.h"


//==============================================================================================
//Definitionen fuer CRYSIN:
//==============================================================================================
CRYSIN_POSITION::CRYSIN_POSITION() {}

CRYSIN_POSITION::CRYSIN_POSITION(CRYSIN_POSITION &ref) : const_add(ref.const_add), x_add(ref.x_add),
                                                         y_add(ref.y_add), z_add(ref.z_add) {}

CRYSIN_POSITION::~CRYSIN_POSITION() {}

CRYSIN::CRYSIN() : group(0), setting(1) {}

CRYSIN::CRYSIN(CRYSIN &crys) : a(crys.a),b(crys.b),c(crys.c),wa(crys.wa),wb(crys.wb),wc(crys.wc),cryst2cart(crys.cryst2cart),
                               cart2cryst(crys.cart2cryst), group(crys.group), setting(crys.setting), rhombo(crys.rhombo) {
    for (vector<stl_ptr<CRYSIN_POSITION> >::iterator it=crys.positions.begin(); it!= crys.positions.end(); ++it) {
        stl_ptr<CRYSIN_POSITION> ncp = new CRYSIN_POSITION(**it);
        positions.push_back(ncp);
    }
}

CRYSIN::~CRYSIN() {
    for (vector<stl_ptr<CRYSIN_POSITION> >::iterator it=positions.begin(); it!=positions.end(); ++it) it->kill();
    positions.clear();
}
