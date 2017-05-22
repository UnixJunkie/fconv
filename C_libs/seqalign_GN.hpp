
//============================================================================
// seqalign_GN.hpp -*- C++ -*-; sequence alignment
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
// This library implements sequence alignment by dynamic programming
// following Needleman Wunsch and Smith Waterman.
//============================================================================


#ifndef __SEQ_ALIGN
#define __SEQ_ALIGN

#include<iostream>
#include<vector>
#include<algorithm>
#include"stl_ptr_GN.hpp"

using namespace std;

template<typename T,typename T2>
class SA_CONTAINER {
    public:
        T at;
        T2 comparer;
        bool operator==(SA_CONTAINER &rechts) {if (rechts.comparer == comparer) return true; else return false;}
};

template<typename T>
class SA_MATRIX_ELEMENT {
    public:
        vector<stl_ptr<SA_MATRIX_ELEMENT<T> > > prev;
        vector<stl_ptr<SA_MATRIX_ELEMENT<T> > > next;
        stl_ptr<T> ival;
        stl_ptr<T> jval;
        int score;
        unsigned int ni;
        unsigned int nj;
        
        SA_MATRIX_ELEMENT(unsigned int i,unsigned int j);
        ~SA_MATRIX_ELEMENT();
};

template<typename T>
class SEQALIGN {
    private:
        vector<T> refvec;
        vector<T> compvec;
        vector<vector<stl_ptr<SA_MATRIX_ELEMENT<T> > > > smatrix;
        T gap_val;
        unsigned int ni; //Anzahl der Zeilen von smatrix
        unsigned int nj; //Anzahl der Spalten von smatrix
        unsigned int max_sm_val; //maximaler score in der smith-waterman matrix
        unsigned int max_val_i;
        unsigned int max_val_j;
        
        void generate_smatrix(bool gap_penalty,bool smith_water = false);
        void traceback();
        void sm_traceback();
    public:
        vector<T> ires;
        vector<T> jres;
        SEQALIGN(vector<T> &rv,vector<T> &cv,T &gv);
        ~SEQALIGN();
        void needleman_wunsch1();
        void needleman_wunsch2();
        void smith_waterman();
};

template<typename T>
SA_MATRIX_ELEMENT<T>::SA_MATRIX_ELEMENT(unsigned int i,unsigned int j): ni(i),nj(j) {}

template<typename T>
SA_MATRIX_ELEMENT<T>::~SA_MATRIX_ELEMENT() {}


template<typename T>
SEQALIGN<T>::SEQALIGN(vector<T> &rv,vector<T> &cv,T &gv): refvec(rv),compvec(cv),gap_val(gv) {
    ni = 0;
    nj = 0;
}

template<typename T>
SEQALIGN<T>::~SEQALIGN() {
    for (unsigned int i=0; i<ni; ++i) {
        for (unsigned int j=0; j<nj; ++j) smatrix[i][j].kill();
    }
}

template<typename T>
void SEQALIGN<T>::generate_smatrix(bool gap_penalty,bool smith_water) {
    ni = refvec.size() + 1;
    nj = compvec.size() + 1;
    
    max_sm_val = 0;
    max_val_i = 0;
    max_val_j = 0;
    
    for (unsigned int i=0; i<ni; ++i) {
        smatrix.push_back(vector<stl_ptr<SA_MATRIX_ELEMENT<T> > >());
        for (unsigned int j=0; j<nj; ++j) {
            stl_ptr<SA_MATRIX_ELEMENT<T> > curr_elem = new SA_MATRIX_ELEMENT<T>(i,j);
            
            if (i != 0) curr_elem->ival = &(refvec[i-1]);
            else curr_elem->ival = gap_val;
            if (j != 0) curr_elem->jval = &(compvec[j-1]);
            else curr_elem->jval = gap_val;
            
            if (i == 0 || j == 0) {
                if (gap_penalty && (!smith_water)) {
                    curr_elem->score = (i+j) * (-2);
                } else curr_elem->score = 0;
            } else {
                int s1 = smatrix[i-1][j-1]->score;
                int s2 = smatrix[i-1][j]->score;
                int s3 = smatrix[i][j-1]->score;
                
                if (s1 >= s2 && s1 >= s3) {    
                    curr_elem->prev.push_back(smatrix[i-1][j-1]);
                    smatrix[i-1][j-1]->next.push_back(curr_elem);
                } else if (s2 > s3) {
                    curr_elem->prev.push_back(smatrix[i-1][j]);
                    smatrix[i-1][j]->next.push_back(curr_elem);
                } else if (s3 > s2) {
                    curr_elem->prev.push_back(smatrix[i][j-1]);
                    smatrix[i][j-1]->next.push_back(curr_elem);
                } else { //2 Wege
                    curr_elem->prev.push_back(smatrix[i-1][j]);
                    smatrix[i-1][j]->next.push_back(curr_elem);
                //    curr_elem.prev.push_back(&(smatrix[i][j-1]));
                //    smatrix[i][j-1].next.push_back(&(curr_elem));
                }
                
                if (curr_elem->ival == curr_elem->jval) s1 += 1;
                else s1 -= 1;
                if (gap_penalty) {s2 -= 2; s3 -= 2;}
                
                if (s1 >= s2 && s1 >= s3) {
                    curr_elem->score = s1;
                } else if (s2 > s3) {
                    curr_elem->score = s2;
                } else if (s3 > s2) {
                    curr_elem->score = s3;
                } else {
                    curr_elem->score = s2;
                }
                if (smith_water) {
                    if (curr_elem->score < 0) curr_elem->score = 0;
                    if (curr_elem->score > int(max_sm_val)) {
                        max_sm_val = curr_elem->score;
                        max_val_i = i;
                        max_val_j = j;
                    }
                }
            }
            smatrix[i].push_back(curr_elem);
        }
    }
}

template<typename T>
void SEQALIGN<T>::traceback() {
    unsigned int ci = ni-1;
    unsigned int cj = nj-1;
    for (unsigned int maxiter=0; maxiter<(ni+nj); ++maxiter) {
        
        stl_ptr<SA_MATRIX_ELEMENT<T> > curr_elem = smatrix[ci][cj];
        
        if (curr_elem->prev.size() == 0) {
            ires.push_back(*(curr_elem->ival));
            jres.push_back(*(curr_elem->jval));
            break;
        }
        if (curr_elem->prev[0]->ni < ci && curr_elem->prev[0]->nj < cj) {
            ires.push_back(*(curr_elem->ival));
            jres.push_back(*(curr_elem->jval));
        } else if (curr_elem->prev[0]->ni < ci) {
            ires.push_back(*(curr_elem->ival));
            jres.push_back(gap_val);
        } else {
            ires.push_back(gap_val);
            jres.push_back(*(curr_elem->jval));
        }
        
        ci = curr_elem->prev[0]->ni;
        cj = curr_elem->prev[0]->nj;
    }
}

template<typename T>
void SEQALIGN<T>::sm_traceback() {
    unsigned int ci = max_val_i;
    unsigned int cj = max_val_j;
    
    for (unsigned int filli=ni-1; filli > ci; --filli) ires.push_back(gap_val);
    for (unsigned int fillj=nj-1; fillj > cj; --fillj) jres.push_back(gap_val);
    
    for (unsigned int maxiter=0; maxiter<(ni+nj); ++maxiter) {
        
        stl_ptr<SA_MATRIX_ELEMENT<T> > curr_elem = smatrix[ci][cj];
        
        if (curr_elem->score == 0) break;
        
        if (ci < max_val_i )
        if (curr_elem->prev.size() == 0) {
            ires.push_back(*(curr_elem->ival));
            jres.push_back(*(curr_elem->jval));
            break;
        }
        if (curr_elem->prev[0]->ni < ci && curr_elem->prev[0]->nj < cj) {
            ires.push_back(*(curr_elem->ival));
            jres.push_back(*(curr_elem->jval));
        } else if (curr_elem->prev[0]->ni < ci) {
            ires.push_back(*(curr_elem->ival));
            jres.push_back(gap_val);
        } else {
            ires.push_back(gap_val);
            jres.push_back(*(curr_elem->jval));
        }
        
        ci = curr_elem->prev[0]->ni;
        cj = curr_elem->prev[0]->nj;
    }
    
    for (unsigned int filli=ci; filli > 0; --filli) ires.push_back(gap_val);
    for (unsigned int fillj=cj; fillj > 0; --fillj) jres.push_back(gap_val);
}

template<typename T>
void SEQALIGN<T>::needleman_wunsch1() {
    generate_smatrix(false);
    traceback();
    reverse(ires.begin(),ires.end());
    reverse(jres.begin(),jres.end());
}

template<typename T>
void SEQALIGN<T>::needleman_wunsch2() {
    generate_smatrix(true);
    traceback();
    reverse(ires.begin(),ires.end());
    reverse(jres.begin(),jres.end());
}

template<typename T>
void SEQALIGN<T>::smith_waterman() {
    generate_smatrix(true,true);
    sm_traceback();
    reverse(ires.begin(),ires.end());
    reverse(jres.begin(),jres.end());
}

#endif

