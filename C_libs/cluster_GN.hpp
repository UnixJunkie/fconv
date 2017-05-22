
//============================================================================
// cluster_GN.hpp -*- C++ -*-; generic implementation of hierarchic clustering
//
// Copyright (C) 2009, 2010 Gerd Neudert
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
// This library is an generic implementation for hierarchic clustering and
// visualization.
//============================================================================

#ifndef __CLUSTERGN
#define __CLUSTERGN

#include<string.h>
#include<list>
#include<vector>
#include<fstream>
#include<string>

using namespace std;

template<typename T>
class D_ELEMENT {
    public:
        T element;
        D_ELEMENT* next;
        
        D_ELEMENT():element(-1),next(0) {}
        D_ELEMENT(T const& ele,D_ELEMENT* n):element(ele),next(n) {}
        ~D_ELEMENT() {};
};

template<typename T> //! T sollte ein integer-typ sein
class DENDRO_NODE {
    public:
        bool old;
        float distance; //! Distanz zwischen left und right
        float x_coord;
        DENDRO_NODE* prev; //! Wurzel => NULL, wenn root-Node
        DENDRO_NODE* left;
        DENDRO_NODE* right;
        D_ELEMENT<T>* elements; //! Nur bei root-Nodes gefuellt!!!
        D_ELEMENT<T>* last;
        int n_elements;
        
        DENDRO_NODE():old(false),distance(0.),x_coord(-1.),prev(0),left(0),right(0),elements(0),last(0),n_elements(0) {}
        DENDRO_NODE(T const& ele):old(false),distance(0.),x_coord(-1.),prev(0),left(0),right(0) {
            elements = new D_ELEMENT<T>(ele,0);
            last = elements;
            n_elements = 1;
        }
        DENDRO_NODE(float const& dst,DENDRO_NODE* l,DENDRO_NODE* r):old(false),distance(dst),x_coord(-1.),prev(0),left(l),right(r) {
            n_elements = l->n_elements + r->n_elements;
            elements = l->elements;
            l->last->next = r->elements;
            last = r->last;
            l->prev = this;
            r->prev = this;
        }
        
        ~DENDRO_NODE() {if (n_elements == 1) delete elements;}
};

template<typename T>
class HIERA_CLUST {
    private:
        float curr_min;
        float* werte;
        bool kill_werte;
        bool has_werte;
        float get_CC_dist_complete(DENDRO_NODE<int>* c1,DENDRO_NODE<int>* c2);
        float get_CC_dist_single(DENDRO_NODE<int>* c1,DENDRO_NODE<int>* c2);
        float get_CC_dist_average(DENDRO_NODE<int>* c1,DENDRO_NODE<int>* c2);
        int dimension;
        vector<DENDRO_NODE<int>*> base;
        
        void calc_werte();
        void get_base(DENDRO_NODE<int>* root);
    public:
        typedef list<DENDRO_NODE<int>*>::iterator cluster_iter;
        typedef D_ELEMENT<int>* ele_ptr;
        vector<T> ini_vec;
        vector<DENDRO_NODE<int>*> dnodes;
        list<DENDRO_NODE<int>*> clusters;
        float (*get_elem_distance)(T&,T&);
        
        vector<string> names_vec; //! wenn ein matrix_values file mit namen eingelesen wird kommen die hier rein
        
        HIERA_CLUST(vector<T>& dn,float (*gd)(T&,T&)); //! zu Clusternde Objekte und Zeiger auf Distanzfunktion
        HIERA_CLUST(char const* fname); //! Name eines Files, welches mit write_matrix_values geschrieben wurde
        HIERA_CLUST(float* w, int dim); //! Wertematrix und Dimension => Nummerierung anhand der Reihenfolge in der Matrix
        ~HIERA_CLUST() {
            for (vector<DENDRO_NODE<int>*>::iterator it=dnodes.begin(); it!=dnodes.end(); ++it) delete *it;
            if (werte && kill_werte) delete[] werte;
        }
        
        void cluster(float const& d_max,float (HIERA_CLUST::*dist_func)(DENDRO_NODE<int>* c1,DENDRO_NODE<int>* c2));
        void cluster_single_link(float const& d_max = -1.);
        void cluster_complete_link(float const& d_max = -1.);
        void cluster_average_link(float const& d_max = -1.);
        float get_wert(int const& a,int const& b) {
            int index = a*dimension + b;
            return (werte[index]);
        }
        
        void write_dendrogram(char const* fname,string (*get_elem_label)(T const&) = 0,float const& height_scale = 15.);
        void write_matrix(char const* fname,string (*get_elem_label)(T const&) = 0);
        
        void write_clusters(char const* fname,string (*get_elem_label)(T const&) = 0);
        void write_matrix_values(char const* fname,string (*get_elem_label)(T const&) = 0);
};

template<typename T>
HIERA_CLUST<T>::HIERA_CLUST(vector<T>& dn,float (*gd)(T&,T&)) {
    kill_werte = true;
    ini_vec = dn;
    get_elem_distance = gd;
    DENDRO_NODE<int>* curr;
    for (int i=0; i<int(dn.size()); ++i) {
        curr = new DENDRO_NODE<int>(i);
        dnodes.push_back(curr);
        clusters.push_back(curr);
    }
    dimension = dn.size();
    has_werte = false;
    werte = 0;
}

template<typename T>
HIERA_CLUST<T>::HIERA_CLUST(char const* fname) {
    kill_werte = true;
    ifstream fin;
    fin.open(fname);
    if (!fin) {
        cerr << "error: HIERA_CLUST::HIERA_CLUST -> could not open " << fname << endl;
        return;
    }
    string row,dummy;
    dimension = -1;
    bool has_names = false;
    while (!fin.eof()) {
        row = "@@@";
        getline(fin,row);
        if (row == "@@@") continue;
        istringstream is;
        is.str(row);
        
        is >> dummy;
        if (dummy == "n_elements") is >> dummy >> dimension;
        else if (dummy == "has_names") {
            is >> dummy >> dummy;
            if (dummy == "yes") has_names = true;
            break;
        }
    }
    if (dimension < 0) {
        cerr << "error: HIERA_CLUST::HIERA_CLUST -> no number of elements given in " << fname << endl;
        return;
    }
    if (dimension > 10000) {
        cerr << "warning: HIERA_CLUST::HIERA_CLUST  -> matrix of " << dimension << " objects will need more than "
             << int((dimension * dimension + dimension ) * 0.000004 + 1) << " MB of memory!" << endl;
    }
    
    for (int i=0; i<dimension; ++i) {
        DENDRO_NODE<int>* curr = new DENDRO_NODE<int>(i);
        dnodes.push_back(curr);
        clusters.push_back(curr);
    }
    
    werte = new float[dimension*dimension];
    for (int i=0; i<dimension; ++i) {
        row = "@@@";
        getline(fin,row);
        if (row == "@@@") {
            cerr << "error: HIERA_CLUST::HIERA_CLUST -> missing lines in " << fname << endl;
            return;
        }
        istringstream is;
        is.str(row);
        int checker = -1;
        is >> checker;
        if (checker != i) {
            cerr << "error: HIERA_CLUST::HIERA_CLUST -> expected element " << i << " but found " << checker << " in " << fname << endl;
            return;
        }
        is >> dummy >> dummy;
        names_vec.push_back(dummy);
        for (int j=i+1; j<dimension; ++j) {
            int index = i*dimension + j;
            is >> werte[index];
            if (is.fail()) {
                cerr << "error: HIERA_CLUST::HIERA_CLUST -> missing columns in " << fname << endl;
                return;
            }
            int index2 = j*dimension + i;
            werte[index2] = werte[index];
        }
    }
    has_werte = true;
    
    fin.close();
}

template<typename T>
HIERA_CLUST<T>::HIERA_CLUST(float* w, int dim) {
    dimension = dim;
    werte = w;
    kill_werte = false;
    has_werte = true;
    
    for (int i=0; i<dimension; ++i) {
        DENDRO_NODE<int>* curr = new DENDRO_NODE<int>(i);
        dnodes.push_back(curr);
        clusters.push_back(curr);
    }
}

template<typename T>
float HIERA_CLUST<T>::get_CC_dist_complete(DENDRO_NODE<int>* c1,DENDRO_NODE<int>* c2) {
    float max_CC = 0.;
    ele_ptr e1 = c1->elements;
    while (e1) {
        ele_ptr e2 = c2->elements;
        while (e2) {
            int index = e1->element*dimension + e2->element;
            if (werte[index] > max_CC) {
                max_CC = werte[index];
                if (max_CC > curr_min) return 999999999.;
            }
            e2 = e2->next;
        }
        e1 = e1->next;
    }
    return max_CC;
}

template<typename T>
float HIERA_CLUST<T>::get_CC_dist_single(DENDRO_NODE<int>* c1,DENDRO_NODE<int>* c2) {
    float min_CC = 999999999.;
    ele_ptr e1 = c1->elements;
    while (e1) {
        ele_ptr e2 = c2->elements;
        while (e2) {
            int index = e1->element*dimension + e2->element;
            if (werte[index] < min_CC) {
                min_CC = werte[index];
            }
            e2 = e2->next;
        }
        e1 = e1->next;
    }
    return min_CC;
}

template<typename T>
float HIERA_CLUST<T>::get_CC_dist_average(DENDRO_NODE<int>* c1,DENDRO_NODE<int>* c2) {
    float av_CC = 0.;
    unsigned int n_CC = 0;
    ele_ptr e1 = c1->elements;
    while (e1) {
        ele_ptr e2 = c2->elements;
        while (e2) {
            int index = e1->element*dimension + e2->element;
            av_CC += werte[index];
            ++n_CC;
            e2 = e2->next;
        }
        e1 = e1->next;
    }
    av_CC /= n_CC;
    return (av_CC);
}

template<typename T>
void HIERA_CLUST<T>::calc_werte() {
    //! Distanzmatrix berechnen:
    if (dimension > 10000) {
        cerr << "warning: HIERA_CLUST::cluster  -> clustering of " << dimension << " objects will need more than "
             << int((dimension * dimension + dimension ) * 0.000004 + 1) << " MB of memory!" << endl;
    }
    werte = new float[dimension*dimension];
    for (int i=0; i<int(dnodes.size()); ++i) {
        for (int j=i+1; j<int(dnodes.size()); ++j) {
            int index = i*dimension + j;
            float dist = get_elem_distance(ini_vec[i],ini_vec[j]);
            werte[index] = dist;
            index = j*dimension + i;
            werte[index] = dist;    
        }
    }
    has_werte = true;
}

template<typename T>
void HIERA_CLUST<T>::cluster(float const& d_max,float (HIERA_CLUST::*dist_func)(DENDRO_NODE<int>* c1,DENDRO_NODE<int>* c2)) {
    if (!has_werte) calc_werte();
    
    float val;
    while (clusters.size() > 1) {
        curr_min = 999999999.;
        DENDRO_NODE<int>* cj1 = 0;
        DENDRO_NODE<int>* cj2 = 0;
        for (list<DENDRO_NODE<int>*>::iterator it=clusters.begin(); it!=clusters.end(); ++it) {
            list<DENDRO_NODE<int>*>::iterator jt=it; ++jt;
            for (; jt!=clusters.end(); ++jt) {
                val = (this->*dist_func)(*it,*jt);
                if (val < curr_min) {
                    curr_min = val;
                    cj1 = *it;
                    cj2 = *jt;
                }
            }
        }
        
        if (curr_min > d_max && d_max > -0.5) break;
        
        DENDRO_NODE<int>* curr = new DENDRO_NODE<int>(curr_min,cj1,cj2);
        clusters.remove(cj1);
        clusters.remove(cj2);
        clusters.push_back(curr);
        dnodes.push_back(curr);
    }
}

template<typename T>
void HIERA_CLUST<T>::cluster_complete_link(float const& d_max) {cluster(d_max,&HIERA_CLUST::get_CC_dist_complete);}
    
template<typename T>
void HIERA_CLUST<T>::cluster_single_link(float const& d_max) {cluster(d_max,&HIERA_CLUST::get_CC_dist_single);}
    
template<typename T>
void HIERA_CLUST<T>::cluster_average_link(float const& d_max) {cluster(d_max,&HIERA_CLUST::get_CC_dist_average);}
    
template<typename T>
void HIERA_CLUST<T>::get_base(DENDRO_NODE<int>* root) {
    //! Diese rekursive Funktion bringt die Blaetter in die richtige Reihenfolge
    //! fuer die Darstellung als Dendrogramm
    DENDRO_NODE<int>* curr = root;
    while (curr->left) curr = curr->left;
    base.push_back(curr);
    if (!curr->prev) return;
    while (curr->prev->right == curr) {
        curr = curr->prev;
        if (!curr->prev) return;
    }
    curr = curr->prev;
    if (curr->right) get_base(curr->right);
}

template<typename T>
void HIERA_CLUST<T>::write_dendrogram(char const* fname,string (*get_elem_label)(T const&),float const& height_scale) {
    float width, w_scale;
    float height, h_scale;
    
    width = dimension * 20.;
    
    w_scale = width / dimension;
    width += 90.;
    
    float max_d_dst = 0.;
    for (list<DENDRO_NODE<int>*>::iterator it=clusters.begin(); it!=clusters.end(); ++it) {
        if ((*it)->distance > max_d_dst) max_d_dst = (*it)->distance;
    }
    
    for (list<DENDRO_NODE<int>*>::iterator it=clusters.begin(); it!=clusters.end(); ++it) {
        get_base(*it);
    }
    
    int n_steps = 1;
    for (vector<DENDRO_NODE<int>*>::iterator it=base.begin(); it!=base.end(); ++it) {
        int cns = 0;
        DENDRO_NODE<int>* cn = *it;
        while (cn->prev) {
            ++cns;
            cn = cn->prev;
        }
        if (cns > n_steps) n_steps = cns;
    }
    
    
    height = n_steps * height_scale;
    h_scale = height / max_d_dst;
    height += 115.;
    
    ofstream fout;
    fout.open(fname);
    
    fout << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    fout << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\"\n";
    fout << "\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n";
    fout << "<svg width=\"" << width << "px\" height=\"" << height << "px\"\n";
    fout << "xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\">\n";
    
    fout << "<!-- This SVG is generated by cluster_GN.hpp -->\n";
    fout << "<!-- " << clusters.size() << " clusters and " << dimension << " elements -->\n";
    
    fout << "<!-- 1.) generate distance scale: -->\n";
    float s0 = height - 105.;
    float s1 = s0 - h_scale * max_d_dst;
    fout << "<g fill=\"none\" stroke=\"red\" stroke-width=\"2\">\n";
    fout << "<line x1=\"10\" y1=\"" << s0 << "\" x2=\"10\" y2=\"" << s1 << "\"/>\n";
    fout << "</g>\n";
    fout << "<g fill=\"none\" stroke=\"red\" stroke-width=\"1\">\n";
    fout << "<line x1=\"5\" y1=\"" << s0 << "\" x2=\"15\" y2=\"" << s0 << "\"/>\n";
    fout << "<line x1=\"5\" y1=\"" << s1 << "\" x2=\"15\" y2=\"" << s1 << "\"/>\n";
    fout << "</g>\n";
    fout << "<text x=\"18\" y=\"" << s0+4 << "\" font-size=\"10\" font-family=\"Verdana\" fill=\"red\">\n";
    fout << "0.0\n";
    fout << "</text>\n";
    fout << "<text x=\"18\" y=\"" << s1+4 << "\" font-size=\"10\" font-family=\"Verdana\" fill=\"red\">\n";
    fout << max_d_dst << "\n";
    fout << "</text>\n";
    
    
    fout << "<!-- 2.) generate labels: -->\n";
    float curr_x = 60.;
    float curr_y = height - 100.;
    for (vector<DENDRO_NODE<int>*>::iterator it=base.begin(); it!=base.end(); ++it) {
        (*it)->x_coord = curr_x + 3.;
        fout << "<g transform=\"translate(" << curr_x << "," << curr_y << ")\">\n";
        fout << "<g transform=\"rotate(90)\">\n";
        fout << "<text x=\"0\" y=\"0\" font-size=\"10\" font-family=\"Verdana\" fill=\"blue\">\n";
        fout << (*it)->elements->element;
        if (get_elem_label) fout << " : " << get_elem_label(ini_vec[(*it)->elements->element]) << "\n";
        else fout << "\n";
        fout << "</text>\n";
        fout << "</g>\n</g>\n";
        curr_x += 20.;
    }
    
    
    fout << "<!-- 3.) generate tree: -->\n";
    float gcurr_y = height - 105.;
    float wy,wy2;
    fout << "<g fill=\"none\" stroke=\"black\" stroke-width=\"1\">\n";
    vector<DENDRO_NODE<int>*> new_base;
    while (base.size() > clusters.size()) {
        for (vector<DENDRO_NODE<int>*>::iterator it=base.begin(); it!=base.end(); ++it) {
            if (!(*it)->prev) continue;
            wy = gcurr_y - h_scale * (*it)->prev->distance;
            wy2 = gcurr_y - h_scale * (*it)->distance;
            
            if (!(*it)->old) {
                fout << "<line x1=\"" << (*it)->x_coord << "\" y1=\"" << wy << "\" x2=\"" 
                     << (*it)->x_coord << "\" y2=\"" << wy2 << "\"/>\n";
                (*it)->old = true;
            }
        }
        for (vector<DENDRO_NODE<int>*>::iterator it=base.begin(); it!=base.end(); ++it) {
            if (!(*it)->prev) continue;
            wy = gcurr_y - h_scale * (*it)->prev->distance;
            wy2 = gcurr_y - h_scale * (*it)->distance;
        
            if ((*it)->prev->left == *it) {
                if ((*it)->prev->right->x_coord > -0.5) {
                    fout << "<line x1=\"" << (*it)->prev->right->x_coord << "\" y1=\"" << wy
                         << "\" x2=\"" << (*it)->x_coord << "\" y2=\"" << wy << "\"/>\n";
                    (*it)->prev->x_coord = ((*it)->prev->right->x_coord - (*it)->x_coord) / 2. + (*it)->x_coord;
                    new_base.push_back((*it)->prev);
                    (*it)->prev->right->x_coord = -1.;
                    (*it)->x_coord = -1.;
                    break;
                }
            }
        }
        
        for (vector<DENDRO_NODE<int>*>::iterator it=base.begin(); it!=base.end(); ++it) {
            if ((*it)->x_coord > -0.5) new_base.push_back(*it);
        }
        base.clear();
        for (vector<DENDRO_NODE<int>*>::iterator it=new_base.begin(); it!=new_base.end(); ++it) {
            base.push_back(*it);
        }
        new_base.clear();
    }
    fout << "</g>\n";
    
    fout << "</svg>\n";
    fout.close();
}


template<typename T>
void HIERA_CLUST<T>::write_matrix(char const* fname,string (*get_elem_label)(T const&)) {
    base.clear();
    
    float width;
    float height;
    
    width = dimension * 20.;
    width += 145.;
    height = width;
    
    float min_wert = 999999999.;
    float max_wert = 0.;
    for (int i=0; i<dimension; ++i) {
        for (int j=i+1; j<dimension; ++j) {
            int index = i*dimension + j;
            if (werte[index] < min_wert) min_wert = werte[index];
            else if (werte[index] > max_wert) max_wert = werte[index];
        }
    }
    float col_scale = max_wert - min_wert;
    if (col_scale < 0.00001) col_scale = 0.00001;
//    col_scale = 253. / col_scale;
    col_scale = 500. / col_scale;
    
    for (list<DENDRO_NODE<int>*>::iterator it=clusters.begin(); it!=clusters.end(); ++it) {
        get_base(*it);
    }
    
    ofstream fout;
    fout.open(fname);
    
    fout << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    fout << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\"\n";
    fout << "\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n";
    fout << "<svg width=\"" << (width + 100.) << "px\" height=\"" << height << "px\"\n";
    fout << "xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\">\n";
    
    fout << "<!-- This SVG is generated by cluster_GN.hpp -->\n";
    
    fout << "<!-- 1.) generate matrix: -->\n";
    float curr_x = 10.;
    for (vector<DENDRO_NODE<int>*>::iterator it=base.begin(); it!=base.end(); ++it) {
        float curr_y = 10.;
        for (vector<DENDRO_NODE<int>*>::iterator jt=base.begin(); jt!=base.end(); ++jt) {
            int rot = 255;
            int gruen = 255;
            int blau = 255;
            if ((*it)->elements->element != (*jt)->elements->element) {
                int index = (*it)->elements->element * dimension + (*jt)->elements->element;
                int cval = int((werte[index]-min_wert) * col_scale + 0.5);
                if (cval < 256) {
                    gruen -= cval;
                    blau -= cval;
                } else {
                    gruen = 0;
                    blau = 0;
                    rot = 510 - cval;
                }
            }
            
            fout << "<rect x=\"" << curr_x << "\" y=\"" << curr_y << "\" width=\"20\" height=\"20\" fill=\"rgb("
                 << rot << "," << gruen << "," << blau << ")\" stroke=\"black\" stroke-width=\"1\" />\n";
            curr_y += 20.;
        }
        curr_x += 20.;
    }
    
    fout << "<!-- 2.) generate labels: -->\n";
    curr_x = 15.;
    float curr_y = 23.;
    float fixed = width - 130.;
    for (vector<DENDRO_NODE<int>*>::iterator it=base.begin(); it!=base.end(); ++it) {
        fout << "<g transform=\"translate(" << curr_x << "," << fixed << ")\">\n";
        fout << "<g transform=\"rotate(90)\">\n";
        fout << "<text x=\"0\" y=\"0\" font-size=\"10\" font-family=\"Verdana\" fill=\"blue\">\n";
        fout << (*it)->elements->element;
        if (get_elem_label) fout << " : " << get_elem_label(ini_vec[(*it)->elements->element]) << "\n";
        else fout << "\n";
        fout << "</text>\n";
        fout << "</g>\n</g>\n";
        curr_x += 20.;
        
        fout << "<text x=\"" << fixed << "\" y=\"" << curr_y << "\" font-size=\"10\" font-family=\"Verdana\" fill=\"blue\">\n";
        fout << (*it)->elements->element;
        if (get_elem_label) fout << " : " << get_elem_label(ini_vec[(*it)->elements->element]) << "\n";
        else fout << "\n";
        fout << "</text>\n";
        curr_y += 20.;
    }
    
    fout << "<!-- 3.) generate scale: -->\n";
    curr_x = width;
    fout << "<rect x=\"" << curr_x << "\" y=\"10\" width=\"30\" height=\"30\" fill=\"rgb("
                 << /*"255,0,0"*/"10,0,0" << ")\" stroke=\"black\" stroke-width=\"2\" />\n";
    fout << "<text x=\"" << curr_x+36. << "\" y=\"28\" font-size=\"10\" font-family=\"Verdana\" fill=\"red\">\n";
    fout << max_wert << "\n";
    fout << "</text>\n";
    
    fout << "<rect x=\"" << curr_x << "\" y=\"40\" width=\"30\" height=\"30\" fill=\"rgb("
                 << /*"255,51,51"*/"110,0,0" << ")\" stroke=\"black\" stroke-width=\"2\" />\n";
    fout << "<text x=\"" << curr_x+36. << "\" y=\"58\" font-size=\"10\" font-family=\"Verdana\" fill=\"red\">\n";
    fout << min_wert + 4. * (max_wert - min_wert) / 5. << "\n";
    fout << "</text>\n";
    
    fout << "<rect x=\"" << curr_x << "\" y=\"70\" width=\"30\" height=\"30\" fill=\"rgb("
                 << /*"255,102,102"*/"210,0,0" << ")\" stroke=\"black\" stroke-width=\"2\" />\n";
    fout << "<text x=\"" << curr_x+36. << "\" y=\"88\" font-size=\"10\" font-family=\"Verdana\" fill=\"red\">\n";
    fout << min_wert + 3. * (max_wert - min_wert) / 5. << "\n";
    fout << "</text>\n";
    
    fout << "<rect x=\"" << curr_x << "\" y=\"100\" width=\"30\" height=\"30\" fill=\"rgb("
                 << /*"255,153,153"*/"255,55,55" << ")\" stroke=\"black\" stroke-width=\"2\" />\n";
    fout << "<text x=\"" << curr_x+36. << "\" y=\"118\" font-size=\"10\" font-family=\"Verdana\" fill=\"red\">\n";
    fout << min_wert + 2. * (max_wert - min_wert) / 5. << "\n";
    fout << "</text>\n";
    
    fout << "<rect x=\"" << curr_x << "\" y=\"130\" width=\"30\" height=\"30\" fill=\"rgb("
                 << /*"255,204,204"*/"255,155,155" << ")\" stroke=\"black\" stroke-width=\"2\" />\n";
    fout << "<text x=\"" << curr_x+36. << "\" y=\"148\" font-size=\"10\" font-family=\"Verdana\" fill=\"red\">\n";
    fout << min_wert + (max_wert - min_wert) / 5. << "\n";
    fout << "</text>\n";
    
    fout << "<rect x=\"" << curr_x << "\" y=\"160\" width=\"30\" height=\"30\" fill=\"rgb("
                 << "255,255,255" << ")\" stroke=\"black\" stroke-width=\"2\" />\n";
    fout << "<text x=\"" << curr_x+36. << "\" y=\"178\" font-size=\"10\" font-family=\"Verdana\" fill=\"red\">\n";
    fout << "0.0 to " << min_wert << "\n";
    fout << "</text>\n";
    
    fout << "</svg>\n";
    fout.close();
}

template<typename T>
void HIERA_CLUST<T>::write_clusters(char const* fname,string (*get_elem_label)(T const&)) {
    base.clear();
    
    ofstream fout;
    fout.open(fname);
    
    fout << "@GENERAL\n";
    fout << "n_elements " << dimension << "\n";
    fout << "n_clusters " << dnodes.size() << "\n";
    fout << "n_final_clusters " << clusters.size() << "\n";
    
    float max_dst = 0.;
    float min_dst = 999999999.;
    for (list<DENDRO_NODE<int>*>::iterator it=clusters.begin(); it!=clusters.end(); ++it) {
        if ((*it)->distance > max_dst) max_dst = (*it)->distance;
    }
    for (list<DENDRO_NODE<int>*>::iterator it=clusters.begin(); it!=clusters.end(); ++it) {
        get_base(*it);
    }
    for (vector<DENDRO_NODE<int>*>::iterator it=base.begin(); it!=base.end(); ++it) {
        if ((*it)->prev) {
            if ((*it)->prev->distance < min_dst) min_dst = (*it)->prev->distance;
        }
    }
    fout << "min_dst " << min_dst << "\n";
    fout << "max_dst " << max_dst << "\n";
    
    for (vector<DENDRO_NODE<int>*>::iterator it=base.begin(); it!=base.end(); ++it) {
        fout << "@NEW_BASE_CLUSTER\n";
        fout << "cluster_id " << *it << "\n";
        fout << "element_id " << (*it)->elements->element << "\n";
        fout << "element_name ";
        if (get_elem_label) fout << get_elem_label(ini_vec[(*it)->elements->element]) << "\n";
        else fout << "\n";
    }
    
    
    //!Jetzt von unten nach oben jeweils die naechste Zusammenfassung angeben:
    list<DENDRO_NODE<int>*> new_base;
    for (vector<DENDRO_NODE<int>*>::iterator it=base.begin(); it!=base.end(); ++it) {
        new_base.push_back(*it);
    }
    
    while (new_base.size() > clusters.size()) {
        list<DENDRO_NODE<int>*>::iterator min_base;// = new_base.begin();
        for (list<DENDRO_NODE<int>*>::iterator it=new_base.begin(); it!=new_base.end(); ++it) {
            if ((*it)->prev) {
                min_base = it;
                break;
            }
        }
        
        for (list<DENDRO_NODE<int>*>::iterator it=new_base.begin(); it!=new_base.end(); ++it) {
            if (!(*it)->prev) continue;
            if ((*it)->prev->distance < (*min_base)->prev->distance) min_base = it;
        }
        
        if ((*min_base)->prev) new_base.push_back((*min_base)->prev);
        new_base.erase(min_base);
        
        if ((*min_base)->prev->left == *min_base) {
            new_base.remove((*min_base)->prev->right);
        } else {
            new_base.remove((*min_base)->prev->left);
        }
        
        fout << "@NEW_CLUSTER\n";
        fout << "cluster_id " << (*min_base)->prev << "\n";
        fout << "n_elements_in_cluster " << (*min_base)->prev->n_elements << "\n";
        fout << "distance " << (*min_base)->prev->distance << "\n";
        fout << "child1_id " << (*min_base)->prev->left << "\n";
        fout << "child2_id " << (*min_base)->prev->right << "\n";
        fout << "n_total_clusters " << new_base.size() << "\n";
        /*
        D_ELEMENT<int>* ele = (*min_base)->prev->elements;
        for (int i=0; i<(*min_base)->prev->n_elements; ++i) {
            fout << ele->element;
            if (get_elem_label) fout << " : " << get_elem_label(ini_vec[ele->element]) << "\n";
            else fout << "\n";
            ele = ele->next;
        }
        */
        //! Jetzt noch einen Repraesentanten fuer den Cluster bestimmen:
        D_ELEMENT<int>* ele = (*min_base)->prev->elements;
        D_ELEMENT<int>* min_ele = ele;
        float min_ele_val = 999999999.;
        for (int i=0; i<(*min_base)->prev->n_elements; ++i) {
            float c_val = 0.;
            D_ELEMENT<int>* o_ele = (*min_base)->prev->elements;
            for (int j=0; j<(*min_base)->prev->n_elements; ++j) {
                if (i == j) continue;
                int index = i*dimension + j;
                c_val += werte[index];
                o_ele = o_ele->next;
            }
            if (c_val < min_ele_val) {
                min_ele_val = c_val;
                min_ele = ele;
            }
            ele = ele->next;
        }
        fout << "average_element " << min_ele->element << "\n";
    }
    fout.close();
}

template<typename T>
void HIERA_CLUST<T>::write_matrix_values(char const* fname,string (*get_elem_label)(T const&)) {
    if (!has_werte) {
        cerr << "error: HIERA_CLUST::write_matrix_values -> no similarity values available" << endl;
        return;
    }
    
    ofstream fout;
    fout.open(fname);
    
    fout << "n_elements = " << dimension << "\n";
    fout << "has_names = ";
    if (get_elem_label) fout << "yes\n";
    else fout << "no\n";
    for (int i=0; i<dimension; ++i) {
        fout << i << " : ";
        if (get_elem_label) fout << get_elem_label(ini_vec[dnodes[i]->elements->element]);
        for (int j=i+1; j<dimension; ++j) {
            int index = i*dimension + j;
            fout << " " << werte[index];
        }
        fout << "\n";
    }
    
    fout.close();
}

#endif
