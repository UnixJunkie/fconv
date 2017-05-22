
//============================================================================
// timer_GN.hpp -*- C++ -*-; simple timer
//
// Copyright (C) 2009 Gerd Neudert
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
//============================================================================


#ifndef __TIMERGN
#define __TIMERGN

#include<time.h>

using namespace std;

template<class T>
class TIMER_GN {
    private:
        clock_t tstart;
        T tcurr;
    public:
        TIMER_GN():tcurr(0.) {}
        ~TIMER_GN() {}
        
        void start_t() {
            tcurr = 0.;
            tstart = clock();
        }
        
        void stop_t() {
            clock_t tend = clock();
            tcurr += tend - tstart;
        }
        
        void continue_t() {tstart = clock();}
        
        T get_time() {return (tcurr/CLOCKS_PER_SEC);}
};

#endif
