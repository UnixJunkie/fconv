
//============================================================================
// bron_kerbosch_GN.hpp -*- C++ -*-; implements the Bron Kerbosch algorithm
//
// Copyright (C) 2008, 2009, 2010 Gerd Neudert
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
// This library is an generic implementation of the well known Bron Kerbosch
// algorithm. It uses an matrix representation for edges if less then 65000
// nodes and otherwise it uses a slower variant with sets.
//============================================================================


#ifndef __BRONKERBOSCHGN
#define __BRONKERBOSCHGN

#include<string.h>
#include<stdint.h>
#include<iostream>
#include<vector>
#include<map>
#include<tr1/unordered_set>
#include<time.h>

using namespace std;

//!===============================================================================
//!Forward Deklarationen:
template<class T> class BK_NODE; //!Diese Klasse repraesentiert einen Knoten im assoziativen Graphen
template<class T> class BK_SOLVER; //!Diese Klasse ermittelt die groessten Cliquen im ass. Graphen
//!===============================================================================


//!===============================================================================
//!Deklaration der Klasse BK_NODE:
template<class T>
class BK_NODE {
    public:
        T obj1; //! Zeiger auf Knoten aus G1
        T obj2; //! Zeiger auf Knoten aus G2
        tr1::unordered_set<BK_NODE<T>*> edges;
        BK_NODE(T const& o1,T const& o2);
        ~BK_NODE();
};
//!===============================================================================


//!===============================================================================
//!Deklaration der Klasse BK_SOLVER:
template<class T> 
class BK_SOLVER {
    private:
        unsigned int curr_n_cliques;
        bool all_cliques;
        unsigned int max_n_cliques; //! nach dieser Anzahl gefundener Cliquen gleicher Groesse wird min_clique_size erhoeht
        unsigned int min_clique_size;
        double max_time;
        double start;
        double ende;
        unsigned int curr_BK_C_size;
        uint32_t *edge_matrix;
        int *g2BK_C;
        
        void get_clique_matrix(int *BK_N,unsigned int BK_N_size,vector<int> &BK_P);
        
        void get_clique_sets(int *BK_N,unsigned int BK_N_size,vector<int> &BK_P);
        
        void solve_with_matrix();
        
        void solve_with_sets();
    public:
        unsigned int max_clique; //! Groesse der groessten gefundenen Clique
        vector<BK_NODE<T> * > nodes; //! Alle Knoten von G
        multimap<unsigned int,vector<BK_NODE<T> * > > cliques; //! Hier kommen alle gefundenen Cliquen rein
        
        BK_SOLVER(vector<BK_NODE<T> * > &bkn); //! Initialisieren mit fertigem vector mit BK_NODEs
                                                       //! => Speicherverwaltung liegt beim User
        ~BK_SOLVER();
        
        void solve(unsigned int min_cs = 1,                  //! Mindestgroesse fuer die Cliquen / sollen auch kleinere gesucht werden /
                   bool ac = false,unsigned int mnc = 1000,  //! Nach welcher Anzahl gleich grosser Cliquen in Folge abbrechen
                   double m_t = -1.);                        //! Wenn mehr Sekunden zwischen 2 new_clique vergeht wird abgebrochen
};
//!===============================================================================


//!===============================================================================
//!Definition der Klasse BK_NODE:
template<class T>
BK_NODE<T>::BK_NODE(T const& o1,T const& o2):obj1(o1),obj2(o2) {}

template<class T>
BK_NODE<T>::~BK_NODE() {}
//!===============================================================================


//!===============================================================================
//!Definition der Klasse BK_SOLVER:
template<class T>
BK_SOLVER<T>::BK_SOLVER(vector<BK_NODE<T> * > &bkn):max_clique(0),nodes(bkn) {}

template<class T>
BK_SOLVER<T>::~BK_SOLVER() {}

template<class T>
void BK_SOLVER<T>::get_clique_matrix(int *BK_N,unsigned int BK_N_size,vector<int> &BK_P) {
    if (curr_n_cliques > max_n_cliques) {
        min_clique_size = max_clique + 1;
        curr_n_cliques -= 2;
    }
    
    //! gucken, ob ein Knoten aus P mit allen Knoten aus N verbunden ist (Abbruchkriterium):
    for (vector<int>::iterator it=BK_P.begin(); it!=BK_P.end(); ++it) {
        bool abbruch = true;
        for (unsigned int jt=0; jt<BK_N_size; ++jt) {
            uint32_t offset = BK_N[jt] * nodes.size() + *it;
            uint32_t index = offset / 32;
            if (!(edge_matrix[index] & (0x01 << (offset % 32)))) {
                abbruch = false;
                break;
            }
        }
        if (abbruch) return; //! Es kann keine maximale Clique mehr gefunden werden
    }
    
    //! Jetzt Knoten fuer Knoten aus BK_N BK_C zufuegen:
    for (unsigned int i=0; i<BK_N_size; ++i) {
        int v = BK_N[i];
        g2BK_C[v] = 1;
        ++curr_BK_C_size;
        
        int e_size = nodes[v]->edges.size() - curr_BK_C_size + 1;
        if (e_size > 0) {
            int *BK_NN = new int[e_size];
            
            unsigned int BK_NN_size = 0;
            for (unsigned int jt=i+1; jt<BK_N_size; ++jt) {
                uint32_t offset = BK_N[jt] * nodes.size() + v;
                uint32_t index = offset / 32;
                if (!(edge_matrix[index] & (0x01 << (offset % 32)))) continue;
                BK_NN[BK_NN_size] = BK_N[jt]; ++BK_NN_size;
            }
        
            vector<int> BK_PN; //! neue Liste bereits zur Erweiterung genutzter Knoten
            for (vector<int>::iterator jt=BK_P.begin(); jt!=BK_P.end(); ++jt) {
                uint32_t offset = (*jt)*nodes.size() + v;
                uint32_t index = offset / 32;
                if (!(edge_matrix[index] & (0x01 << (offset % 32)))) continue;
                BK_PN.push_back(*jt);
            }
        
            if (BK_PN.size() == 0 && BK_NN_size == 0) {
                //! Clique gefunden:
                //! Die folgende Bedingung size() >= max_clique  stellt sicher, dass wenn bereits eine "gute Reihenfolge"
                //! in G1 und G2 vorgegeben ist und somit die groessten Cliquen als erstes gefunden werden, nicht extra
                //! eine neue new_clique Liste angelegt werden muss
                if (all_cliques && (curr_BK_C_size >= min_clique_size)) { //! Es sollen alle Cliquen gefunden werden
                    
                    ende = clock();
                    if ((ende - start) > max_time) min_clique_size = 999999;
                    start = clock();
                    
                    if (curr_BK_C_size > max_clique) {
                        curr_n_cliques = 0;
                        max_clique = curr_BK_C_size;
                    }
                    ++curr_n_cliques;
                    vector<BK_NODE<T> * > new_clique;
                    for (unsigned int nc=0; nc<nodes.size(); ++nc) {
                        if (g2BK_C[nc]) new_clique.push_back(nodes[nc]);
                    }
                    cliques.insert(pair<unsigned int,vector<BK_NODE<T> * > >(new_clique.size(),new_clique));
                } else if ((curr_BK_C_size >= max_clique) && (curr_BK_C_size >= min_clique_size)) { //! keine kleineren als max_clique mehr anlegen
                    
                    ende = clock();
                    if ((ende - start) > max_time) min_clique_size = 999999;
                    start = clock();
                    
                    if (curr_BK_C_size > max_clique) {
                        curr_n_cliques = 0;
                        max_clique = curr_BK_C_size;
                    }
                    ++curr_n_cliques;
                    vector<BK_NODE<T> * > new_clique;
                    for (unsigned int nc=0; nc<nodes.size(); ++nc) {
                        if (g2BK_C[nc]) new_clique.push_back(nodes[nc]);
                    }
                    
                    cliques.insert(pair<unsigned int,vector<BK_NODE<T> * > >(new_clique.size(),new_clique));
                }
            } else {
                if ((curr_BK_C_size + BK_NN_size) >= min_clique_size) {
                    if (all_cliques || (curr_BK_C_size + BK_NN_size >= max_clique)) get_clique_matrix(BK_NN,BK_NN_size,BK_PN);
                }
            }
            delete[] BK_NN;
        } else {
            //! Clique gefunden:
            //! Die folgende Bedingung size() >= max_clique  stellt sicher, dass wenn bereits eine "gute Reihenfolge"
            //! in G1 und G2 vorgegeben ist und somit die groessten Cliquen als erstes gefunden werden, nicht extra
            //! eine neue new_clique Liste angelegt werden muss
            if (all_cliques && (curr_BK_C_size >= min_clique_size)) { //! Es sollen alle Cliquen gefunden werden
                
                ende = clock();
                if ((ende - start) > max_time) min_clique_size = 999999;
                start = clock();
                    
                if (curr_BK_C_size > max_clique) {
                    curr_n_cliques = 0;
                    max_clique = curr_BK_C_size;
                }
                ++curr_n_cliques;
                vector<BK_NODE<T> * > new_clique;
                for (unsigned int nc=0; nc<nodes.size(); ++nc) {
                    if (g2BK_C[nc]) new_clique.push_back(nodes[nc]);
                }
                cliques.insert(pair<unsigned int,vector<BK_NODE<T> * > >(new_clique.size(),new_clique));
            } else if ((curr_BK_C_size >= max_clique) && (curr_BK_C_size >= min_clique_size)) { //! keine kleineren als max_clique mehr anlegen
                
                ende = clock();
                if ((ende - start) > max_time) min_clique_size = 999999;
                start = clock();
                    
                if (curr_BK_C_size > max_clique) {
                    curr_n_cliques = 0;
                    max_clique = curr_BK_C_size;
                }
                ++curr_n_cliques;
                vector<BK_NODE<T> * > new_clique;
                for (unsigned int nc=0; nc<nodes.size(); ++nc) {
                    if (g2BK_C[nc]) new_clique.push_back(nodes[nc]);
                }
                
                cliques.insert(pair<unsigned int,vector<BK_NODE<T> * > >(new_clique.size(),new_clique));
            }
        }
        
        //! von BK_C nach BK_P verschieben:
        g2BK_C[v] = 0;
        --curr_BK_C_size;
        BK_P.push_back(v);
    }
}

template<class T>
void BK_SOLVER<T>::get_clique_sets(int *BK_N,unsigned int BK_N_size,vector<int> &BK_P) {
    if (curr_n_cliques > max_n_cliques) {
        min_clique_size = max_clique + 1;
        curr_n_cliques -= 2;
    }
    
    //! gucken, ob ein Knoten aus P mit allen Knoten aus N verbunden ist (Abbruchkriterium):
    for (vector<int>::iterator it=BK_P.begin(); it!=BK_P.end(); ++it) {
        bool abbruch = true;
        for (unsigned int jt=0; jt<BK_N_size; ++jt) {
            if (nodes[*it]->edges.find(nodes[BK_N[jt]]) == nodes[*it]->edges.end()) {
                abbruch = false;
                break;
            }
        }
        if (abbruch) return; //! Es kann keine maximale Clique mehr gefunden werden
    }
    
    //! Jetzt Knoten fuer Knoten aus BK_N BK_C zufuegen:
    for (unsigned int i=0; i<BK_N_size; ++i) {
        int v = BK_N[i];
        g2BK_C[v] = 1;
        ++curr_BK_C_size;
        
        int e_size = nodes[v]->edges.size() - curr_BK_C_size + 1;
        if (e_size > 0) {
            int *BK_NN = new int[e_size];
            unsigned int BK_NN_size = 0;
            for (unsigned int jt=i+1; jt<BK_N_size; ++jt) {
                if (nodes[v]->edges.find(nodes[BK_N[jt]]) == nodes[v]->edges.end()) continue;
                BK_NN[BK_NN_size] = BK_N[jt]; ++BK_NN_size;
            }
        
            vector<int> BK_PN; //! neue Liste bereits zur Erweiterung genutzter Knoten
            for (vector<int>::iterator jt=BK_P.begin(); jt!=BK_P.end(); ++jt) {
                if (nodes[v]->edges.find(nodes[*jt]) == nodes[v]->edges.end()) continue;
                BK_PN.push_back(*jt);
            }
        
            if (BK_PN.size() == 0 && BK_NN_size == 0) {
                //! Clique gefunden:
                //! Die folgende Bedingung size() >= max_clique  stellt sicher, dass wenn bereits eine "gute Reihenfolge"
                //! in G1 und G2 vorgegeben ist und somit die groessten Cliquen als erstes gefunden werden, nicht extra
                //! eine neue new_clique Liste angelegt werden muss
                if (all_cliques && (curr_BK_C_size >= min_clique_size)) { //! Es sollen alle Cliquen gefunden werden
                    
                    ende = clock();
                    if ((ende - start) > max_time) min_clique_size = 999999;
                    start = clock();
                    
                    if (curr_BK_C_size > max_clique) {
                        curr_n_cliques = 0;
                        max_clique = curr_BK_C_size;
                    }
                    ++curr_n_cliques;
                    vector<BK_NODE<T> * > new_clique;
                    for (unsigned int nc=0; nc<nodes.size(); ++nc) {
                        if (g2BK_C[nc]) new_clique.push_back(nodes[nc]);
                    }
                    cliques.insert(pair<unsigned int,vector<BK_NODE<T> * > >(new_clique.size(),new_clique));
                } else if ((curr_BK_C_size >= max_clique) && (curr_BK_C_size >= min_clique_size)) { //! keine kleineren als max_clique mehr anlegen
                    
                    ende = clock();
                    if ((ende - start) > max_time) min_clique_size = 999999;
                    start = clock();
                    
                    if (curr_BK_C_size > max_clique) {
                        curr_n_cliques = 0;
                        max_clique = curr_BK_C_size;
                    }
                    ++curr_n_cliques;
                    vector<BK_NODE<T> * > new_clique;
                    for (unsigned int nc=0; nc<nodes.size(); ++nc) {
                        if (g2BK_C[nc]) new_clique.push_back(nodes[nc]);
                    }
                    cliques.insert(pair<unsigned int,vector<BK_NODE<T> * > >(new_clique.size(),new_clique));
                }
            } else {
                if ((curr_BK_C_size + BK_NN_size) >= min_clique_size) {
                    if (all_cliques || (curr_BK_C_size + BK_NN_size >= max_clique)) get_clique_sets(BK_NN,BK_NN_size,BK_PN);
                }
            }
            delete[] BK_NN;
        } else {
            //! Clique gefunden:
            //! Die folgende Bedingung size() >= max_clique  stellt sicher, dass wenn bereits eine "gute Reihenfolge"
            //! in G1 und G2 vorgegeben ist und somit die groessten Cliquen als erstes gefunden werden, nicht extra
            //! eine neue new_clique Liste angelegt werden muss
            if (all_cliques && (curr_BK_C_size >= min_clique_size)) { //! Es sollen alle Cliquen gefunden werden
                
                ende = clock();
                if ((ende - start) > max_time) min_clique_size = 999999;
                start = clock();
                    
                if (curr_BK_C_size > max_clique) {
                    curr_n_cliques = 0;
                    max_clique = curr_BK_C_size;
                }
                ++curr_n_cliques;
                vector<BK_NODE<T> * > new_clique;
                for (unsigned int nc=0; nc<nodes.size(); ++nc) {
                    if (g2BK_C[nc]) new_clique.push_back(nodes[nc]);
                }
                cliques.insert(pair<unsigned int,vector<BK_NODE<T> * > >(new_clique.size(),new_clique));
            } else if ((curr_BK_C_size >= max_clique) && (curr_BK_C_size >= min_clique_size)) { //! keine kleineren als max_clique mehr anlegen
                
                ende = clock();
                if ((ende - start) > max_time) min_clique_size = 999999;
                start = clock();
                
                if (curr_BK_C_size > max_clique) {
                    curr_n_cliques = 0;
                    max_clique = curr_BK_C_size;
                }
                ++curr_n_cliques;
                vector<BK_NODE<T> * > new_clique;
                for (unsigned int nc=0; nc<nodes.size(); ++nc) {
                    if (g2BK_C[nc]) new_clique.push_back(nodes[nc]);
                }
                cliques.insert(pair<unsigned int,vector<BK_NODE<T> * > >(new_clique.size(),new_clique));
            }
        }
        
        //! von BK_C nach BK_P verschieben:
        g2BK_C[v] = 0;
        --curr_BK_C_size;
        BK_P.push_back(v);
    }
}

template<class T>
void BK_SOLVER<T>::solve_with_matrix() {
    g2BK_C = new int[nodes.size()];
    curr_BK_C_size = 0;
    int *BK_N = new int[nodes.size()];
    for (int i=0; i<int(nodes.size()); ++i) BK_N[i] = i;
    vector<int> BK_P; //! bereits in frueheren Rekursionen zur Erweiterung benutzt
    
    //!Bitmatrix fuer die Kanten erstellen:
    uint32_t ne = nodes.size() * nodes.size();
    ne /= 32; ne += 1;
    edge_matrix = new uint32_t[ne];
    memset(edge_matrix,0,ne*4);
    for (int i=0; i<int(nodes.size()); ++i) {
        g2BK_C[i] = 0;
        for (int j=i+1; j<int(nodes.size()); ++j) {
            if (nodes[i]->edges.find(nodes[j]) != nodes[i]->edges.end()) {
                uint32_t offset = j*nodes.size() + i;
                uint32_t index = offset / 32;
                uint32_t shift = offset % 32;
                uint32_t setB = 0x01; setB <<= shift;
                edge_matrix[index] = (edge_matrix[index] | setB);
                offset = i*nodes.size() + j;
                index = offset / 32;
                shift = offset % 32;
                setB = 0x01; setB <<= shift;
                edge_matrix[index] = (edge_matrix[index] | setB);
            }
        }
    }
    
    start = clock();
    
    get_clique_matrix(BK_N,nodes.size(),BK_P);
    delete[] edge_matrix;
    delete[] BK_N;
    delete[] g2BK_C;
}

template<class T>
void BK_SOLVER<T>::solve_with_sets() {
    g2BK_C = new int[nodes.size()];
    curr_BK_C_size = 0;
    int *BK_N = new int[nodes.size()];
    for (int i=0; i<int(nodes.size()); ++i) {
        BK_N[i] = i;
        g2BK_C[i] = 0;
    }
    vector<int> BK_P; //! bereits in frueheren Rekursionen zur Erweiterung benutzt
    
    start = clock();
    
    get_clique_sets(BK_N,nodes.size(),BK_P);
    delete[] BK_N;
    delete[] g2BK_C;
}

template<class T>
void BK_SOLVER<T>::solve(unsigned int min_cs,bool ac,unsigned int mnc,double m_t) {
    all_cliques = ac;
    max_n_cliques = mnc;
    min_clique_size = min_cs;
    if (m_t > 0.) max_time = m_t * CLOCKS_PER_SEC;
    else max_time = 3600. * CLOCKS_PER_SEC;
    curr_n_cliques = 0;
    
    //! Mit Matrix gehts schneller, aber bei zuvielen Knoten wird der Speicher vollgeballert
    //! => Bei ueber 65000 Knoten die langsamere Variante
    //! (Speicherbedarf der edge_matrix:  nodes*nodes/8 Bytes)
    if (nodes.size() < 65000) solve_with_matrix();
    else {
        if (m_t > 0.) max_time /= 10.;
        solve_with_sets();
    }
}
//!===============================================================================


#endif
