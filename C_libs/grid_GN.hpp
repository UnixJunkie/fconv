
//============================================================================
// grid_GN.hpp -*- C++ -*-; n-dimensional container
//
// Copyright (C) 2006, 2007, 2008 Gerd Neudert
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
// This library implements and n-dimensional container, together with
// some iterators.
//============================================================================


#ifndef __GRIDGN
#define __GRIDGN

#include<vector>
#include<cmath>
#include<fstream>
#include"message_GN.hpp"

using namespace std;
using namespace TXT;

//==================================================================================================
//!Forward-Deklarationen:
//==================================================================================================

template<int n_dim,class T> class GRID; //Feld beliebiger Dimension
template<int n_dim,class T> class GRIDERATOR; //Iterator fr die GRID-Klasse
template<int n_dim,class T> class FIXED_GRIDERATOR; //Iteraor, der eine Dimension festhaelt
template<int n_dim,class T> class SHELL_GRIDERATOR; //Iterator fuer eine n-dimensionale Schicht
template<int n_dim,class T> class SUB_GRIDERATOR; //Iterator fuer ein Sub-Grid
template<int n_dim,int cur_dim,class T> class GRIDINDEXER; //um den []-Operator zu realisieren


//!*************************************************************************************************
//!                                   KLASSEN-DEKLARATIONEN                                        *
//!*************************************************************************************************

//==================================================================================================
//!Deklaration der Klasse GRID:
//==================================================================================================

template<int n_dim,class T> //'n_dim'-dimensionales Gitter mit den Ausdehnungen 'sizes' und den Gitterelementen 'T'
class GRID {
    private:
        int sizes[n_dim]; //Ausdehnung des Gitters
        vector<T> val; //Vektor mit Gitterpunkten
        int n_elements; //Anzahl der Gitterpunkte
        int read_index; //Hilfsvariable
        int mulsizes[n_dim]; //mulsizes[i] = Produkt ber alle sizes[n] fr n=0 bis i
        GRIDERATOR<n_dim,T> *begin_iter; //schon berechneter Iterator (verhindert Neukonstruktion bei wiederholter Abfrage)
        GRIDERATOR<n_dim,T> *rbegin_iter;//=>Beschleunigung 
        GRIDERATOR<n_dim,T> *end_iter;
        GRIDERATOR<n_dim,T> *rend_iter;
        FIXED_GRIDERATOR<n_dim,T> *fix_begin_iter[n_dim]; //fixed_iterator (verhindert Neukonstruktion bei wiederholter Abfrage)
        FIXED_GRIDERATOR<n_dim,T> *fix_rbegin_iter[n_dim];//(wie z.B. check auf end() in Schleifen)
        FIXED_GRIDERATOR<n_dim,T> *fix_end_iter[n_dim];
        FIXED_GRIDERATOR<n_dim,T> *fix_rend_iter[n_dim];
        SHELL_GRIDERATOR<n_dim,T> *shell_begin_iter; //fuer shell_iterator vorerst nur begin und end und NICHT fuer
        SHELL_GRIDERATOR<n_dim,T> *shell_end_iter;   //alle shell's einzeln!
        
        inline void calc_mulsizes();
        inline void clear_iter(); //Iteratoren nach einem resize zurcksetzen
        inline void get_max_shell(); //max_shell berechnen
    public:
        typedef T value_type;
        typedef GRIDERATOR<n_dim,T> iterator; //Iterator ber alle Gitterelemente (ungeprft!)
        typedef FIXED_GRIDERATOR<n_dim,T> fixed_iterator; //Iterator ber alle Gitterelemente mit einer festgelegten Dimension
        typedef SUB_GRIDERATOR<n_dim,T> sub_iterator; //Iterator ber ein SubGrid (z.B. x1,y1,z1 bis x2,y2,z2)
        typedef SHELL_GRIDERATOR<n_dim,T> shell_iterator; //Iterator fuer z.B. eine Wuerfelschicht
        int max_shell; //fuer n_dim = 3 waere dies z.B. die Anzahl der Wuerfelschalen
        
        GRID();  //Standardkonstruktor fr Verwendung in Containern
        explicit GRID(int *s); //mit array fr die Ausdehnungen der verschiedenen Dimensionen anlegen
        explicit GRID(int *s,T standard); //mit Standardwert fllen
        explicit GRID(int s); //gleiche Groesse in jede Richtung
        explicit GRID(int s,T standard); //...
        explicit GRID(GRID<n_dim,T> &ref); //! tiefe Kopie
        
        inline GRID<n_dim,T>& operator=(GRID<n_dim,T> &rechts);

        ~GRID();
        
        //es folgen die Mehtoden zum Iteratoren erzeugen
        inline iterator& begin(); //Zeiger (Iterator) auf das erste Gitterelement
        inline iterator& end();   //Iterator hinter das letzte Element
        inline iterator& rbegin(); //Iterator auf das letzte Gitterelement
        inline iterator& rend(); //Iterator vor das erste Gitterelement
        
        inline fixed_iterator& begin(int fixed_dim);
        inline fixed_iterator& end(int fixed_dim);
        inline fixed_iterator& rbegin(int fixed_dim);
        inline fixed_iterator& rend(int fixed_dim);
        
        inline shell_iterator& sbegin(int shell);
        inline shell_iterator& send(int shell);
        
        //es folgen die Methoden zur Indexumwandlung
        inline int get_direct_index(int *index); //direkten index fuer *index berechnen
        inline int* get_coords(int *coord,int index); //andersrum (int coord[n_dim] muss uebergeben werden)
                                                      //der zurckgegebene Zeiger zeigt also auf coord
                                  //!=> lieber mal ein void draus machen
        inline int get_new_index(int *oldmulsizes,int *newmulsizes,int oldindex); //internen Index fr Gitter mit 
                                                                                  //verschiedenen sizes[n_dim] umrechnen 
                                                                                          //(zum Umkopieren bei Gr�en�derung)
        inline int* get_mulsizes(int *muls,int *s); //mulsizes fr die sizes s berechnen
        
        //es folgen die Methoden zur Groessenaenderung
        inline void resize(int *s); //!Achtung: altes array wird geloescht
        inline void resize(int s); //...
        inline void resize(int *s,T standard); //Nach Groessenaenderung mit Standardelementen fllen
        inline void resize(int s,T standard); //...
        inline void smart_resize(int *newsize); //!Elemente bleiben erhalten und behalten auch ihren index
        inline void smart_resize(int newsize); //!Iteratoren werden ungltig (bzw. zeigen auf falsche Elemente
        inline void smart_resize(int dim,int newsize); //Groesse der Dimension dim auf newsize aendern
        inline void smart_resize(int *newsize,T standard); //Die neuen Zellen mit Standard fuellen
        inline void smart_resize(int dim,int newsize,T standard); //...

        //es folgen Zugriffmethoden
        inline T& value(int *index); //Gitterwert lesen und Referenz liefern (also auch zum schreiben!!!)
        inline T& dvalue(int index); //direkte Angabe des entsprechenden Index von val (wird von den Iteratoren genutzt)
        inline T* get_pos(int *index); //Zeiger auf die entsprechende Gitterposition
        inline T* dget_pos(int index); //direkte Angabe des entsprechenden Index von val (wird von den Iteratoren genutzt)
        inline T& front(); //Referenz auf das erste Gitterelement
        inline T& back(); //Referenz auf das letzte Gitterelement
        
        inline void push_back(T& value); //einfach an den intern verwendeten vector anhaengen => mulsizes usw.
                         //neu berechnen / Es wird der letzten Dimension hinzugefuegt!!!
        
        inline GRIDINDEXER<n_dim,n_dim-1,T> operator[](int index); //liefert ein Objekt, auf das  wieder der
                                                                   // []-Operator anwendbar ist
        
        //sonstige
        inline int size(); //Anzahl der Gitterelemente liefern
        inline int get_n_elements(); //Speicherbedarf angeben
        inline int* get_sizes(); //Zeiger auf sizes
        inline void clear(); //alle elemente loeschen, aber erhaelt die Kapazitaet! (Unterschied zu stl-Containern)
        inline void full_clear(); //wie clear fuer stl-Container
        inline void set_all(T &standard); //alle Elemente auf den Wert standard setzen
        
        //Methoden zum Lesen und Speichern eines Gitters
        inline void save(const char *filename,bool binary = true); //Gitter abspeichern
        inline void load(const char *filename,bool binary = true); //Gitter laden (sizes wird automatisch angepasst
                                                              //stimmt die Dimension nicht gibt es einen error
        
        friend class GRIDERATOR<n_dim,T>; //braucht Zugriff auf sizes[n_dim]
        friend class FIXED_GRIDERATOR<n_dim,T>; //braucht Zugriff auf sizes[n_dim]
        friend class SUB_GRIDERATOR<n_dim,T>; //...
        friend class SHELL_GRIDERATOR<n_dim,T>;
        template<int n,int cur_dim,class T2> friend class GRIDINDEXER; //greift auf mulsizes[] und val zu
};


//==================================================================================================
//!Deklaration der Klasse GRIDERATOR:
//==================================================================================================

template<int n_dim,class T>
class GRIDERATOR {
    protected:
        GRID<n_dim,T> *grid; //Zeiger auf das Gitter zu dem der Iterator gehoert
        int index; //aktueller index fr GRID<n_dim,T>::val
    public:
        //Es folgen die Konstruktoren, Zuweisungsoperatoren und Destruktor
        explicit GRIDERATOR(GRID<n_dim,T> *grd,int i); //Iterator fr das ganze GRID mit direktem Startindex
        explicit GRIDERATOR(GRID<n_dim,T> *grd,int *i); //mit Startindex int index[n_dim]
        
        inline GRIDERATOR<n_dim,T>& operator=(GRIDERATOR<n_dim,T> &rechts); //Zuweisungsoperator
        inline GRIDERATOR<n_dim,T>& operator=(int *i); //Iterator auf entsprechenden Index setzen
        
        virtual ~GRIDERATOR(); //virtual, weil Basisklasse fr weitere Iteratoren

        //es folgen die Zugriffsoperatoren
        inline T& operator*(); //Dereferenzierungsoperator liefert eine Referenz auf das aktuelle Gitterelement
        inline T* operator->(); //Dereferenzierungsoperator liefert Zeiger auf aktuelles Gitterelement
        inline int get_index();

        //es folgen die Vergleichsoperatoren
        inline bool operator==(GRIDERATOR<n_dim,T> rechts); //!kann hier nicht als Referenz bergeben werden, damit die
        inline bool operator!=(GRIDERATOR<n_dim,T> rechts); //!Anwendung von begin() usw. funktioniert (das Objekt wird
                                                            //!in diesem Fall trotzdem nur einmal an der entsprechenden
                                                                    //!Stelle erzeugt
        
        //es folgen die Operatoren fr Zeigerarithmetik
        inline GRIDERATOR<n_dim,T>& operator++(); //Pr�ix
        inline GRIDERATOR<n_dim,T> operator++(int); //Postfix
        inline GRIDERATOR<n_dim,T>& operator--();
        inline GRIDERATOR<n_dim,T> operator--(int);
        inline GRIDERATOR<n_dim,T>& operator+=(int inc);
        inline GRIDERATOR<n_dim,T>& operator-=(int dec);
};


//==================================================================================================
//!Deklaration der Klasse FIXED_GRIDERATOR:
//==================================================================================================

template<int n_dim,class T>
class FIXED_GRIDERATOR : public GRIDERATOR<n_dim,T> {
    protected:
        typedef GRIDERATOR<n_dim,T> base;
        
        int dim; //!von 0 bis n_dim-1
        int n_step; //wie oft kann der Index um 1 erh�t werden
        int n_jump; //wie weit mu�gesprungen werden, wenn cur_step == n_step
        int cur_step; //wie oft wurde bereits um 1 erh�t
    public:
        //es folgen die Konstruktoren, Zuweisungsoperatoren und Destruktor
        explicit FIXED_GRIDERATOR(GRID<n_dim,T> *grd,int i,int fixed_dim);
        explicit FIXED_GRIDERATOR(GRID<n_dim,T> *grd,int *index,int fixed_dim);
        
        inline FIXED_GRIDERATOR<n_dim,T>& operator=(FIXED_GRIDERATOR<n_dim,T> &rechts); //Zuweisungsoperator
        inline FIXED_GRIDERATOR<n_dim,T>& operator=(int *i); //Iterator auf entsprechenden Index setzen
        
        ~FIXED_GRIDERATOR();
        
        //es folgen die Vergleichsoperatoren
        inline bool operator==(FIXED_GRIDERATOR<n_dim,T> rechts); //!siehe oben
        inline bool operator!=(FIXED_GRIDERATOR<n_dim,T> rechts);
        
        //es folgen die Operatoren fr Zeigerarithmetik
        inline FIXED_GRIDERATOR<n_dim,T>& operator++();
        inline FIXED_GRIDERATOR<n_dim,T> operator++(int);
        inline FIXED_GRIDERATOR<n_dim,T>& operator--();
        inline FIXED_GRIDERATOR<n_dim,T> operator--(int);
        inline FIXED_GRIDERATOR<n_dim,T>& operator+=(int inc);
        inline FIXED_GRIDERATOR<n_dim,T>& operator-=(int dec);
};


//==================================================================================================
//!Deklaration der Klasse SHELL_GRIDERATOR:
//==================================================================================================

template<int n_dim,class T>
class SHELL_GRIDERATOR : public GRIDERATOR<n_dim,T> {
    protected:
        typedef GRIDERATOR<n_dim,T> base;
        int op_ind[n_dim];
    public:
        int shell;
        //es folgen die Konstruktoren, Zuweisungsoperatoren und Destruktor
        explicit SHELL_GRIDERATOR(GRID<n_dim,T> *grd,int i,int sh);
        explicit SHELL_GRIDERATOR(GRID<n_dim,T> *grd,int *index,int sh);
        
        inline SHELL_GRIDERATOR<n_dim,T>& operator=(SHELL_GRIDERATOR<n_dim,T> &rechts); //Zuweisungsoperator
        inline SHELL_GRIDERATOR<n_dim,T>& operator=(int *i); //Iterator auf entsprechenden Index setzen
        
        ~SHELL_GRIDERATOR();
        
        //es folgen die Vergleichsoperatoren
        inline bool operator==(SHELL_GRIDERATOR<n_dim,T> rechts);
        inline bool operator!=(SHELL_GRIDERATOR<n_dim,T> rechts);
        
        //es folgen die Operatoren fr Zeigerarithmetik
        inline SHELL_GRIDERATOR<n_dim,T>& operator++();
        inline SHELL_GRIDERATOR<n_dim,T> operator++(int);
        inline SHELL_GRIDERATOR<n_dim,T>& operator--();
        inline SHELL_GRIDERATOR<n_dim,T> operator--(int);
        inline SHELL_GRIDERATOR<n_dim,T>& operator+=(int inc);
        inline SHELL_GRIDERATOR<n_dim,T>& operator-=(int dec);
};


//==================================================================================================
//!Deklaration der Klasse SUB_GRIDERATOR:
//==================================================================================================

template<int n_dim,class T>
class SUB_GRIDERATOR : public GRIDERATOR<n_dim,T> {
    protected:
        typedef GRIDERATOR<n_dim,T> base;
        
    //    int n_step; //wie oft kann der Index um 1 erh�t werden
    //    int n_jump; //wie weit mu�gesprungen werden, wenn cur_step == n_step
    //    int cur_step; //wie oft wurde bereits um 1 erh�t
    public:
        //es folgen die Konstruktoren, Zuweisungsoperatoren und Destruktor
        explicit SUB_GRIDERATOR(GRID<n_dim,T> *grd,int i_start,int i_end);
    //    explicit SUB_GRIDERATOR(GRID<n_dim,T> *grd,int *i_start,int *i_end);
        
    //    inline SUB_GRIDERATOR<n_dim,T>& operator=(SUB_GRIDERATOR<n_dim,T> &rechts); //Zuweisungsoperator
    //    inline SUB_GRIDERATOR<n_dim,T>& operator=(int *i); //Iterator auf entsprechenden Index setzen
        
        ~SUB_GRIDERATOR();
        
        //es folgen die Vergleichsoperatoren
    //    inline bool operator==(FIXED_GRIDERATOR<n_dim,T> rechts); //!siehe oben
    //    inline bool operator!=(FIXED_GRIDERATOR<n_dim,T> rechts);
        
        //es folgen die Operatoren fr Zeigerarithmetik
    //    inline FIXED_GRIDERATOR<n_dim,T>& operator++();
    //    inline FIXED_GRIDERATOR<n_dim,T> operator++(int);
    //    inline FIXED_GRIDERATOR<n_dim,T>& operator--();
    //    inline FIXED_GRIDERATOR<n_dim,T> operator--(int);
    //    inline FIXED_GRIDERATOR<n_dim,T>& operator+=(int inc);
    //    inline FIXED_GRIDERATOR<n_dim,T>& operator-=(int dec);
};


//==================================================================================================
//!Deklaration der Klasse GRIDINDEXER:
//==================================================================================================

template<int n_dim,int cur_dim,class T>
class GRIDINDEXER { //Hilfsklasse fr GRID - liefert ein Objekt, auf das wieder der []-Operator anwendbar ist
    private:
        GRID<n_dim,T> &grid;
    public:
        GRIDINDEXER (GRID<n_dim,T> &grd);
        ~GRIDINDEXER();

        inline GRIDINDEXER<n_dim,cur_dim-1,T> operator[](int index);
};

template<int n_dim,class T>
class GRIDINDEXER<n_dim,1,T> { //Spezialisierung fr cur_dim == 1   => liefert die eigentliche Referenz auf T
    private:
        GRID<n_dim,T> &grid;
    public:
        GRIDINDEXER (GRID<n_dim,T> &grd);
        ~GRIDINDEXER();

        inline T& operator[](int index);
};



//!*************************************************************************************************
//!                                   KLASSEN-DEFINITIONEN                                         *
//!*************************************************************************************************

//==================================================================================================
//!Definition der Klasse GRID:
//==================================================================================================

template<int n_dim,class T> 
GRID<n_dim,T>::GRID():n_elements(1),begin_iter(NULL),rbegin_iter(NULL),end_iter(NULL),rend_iter(NULL),
                      shell_begin_iter(NULL), shell_end_iter(NULL) {
    for (int i = 0; i<n_dim; ++i) {
        sizes[i] = 1;
        fix_begin_iter[i] = NULL;
        fix_rbegin_iter[i] = NULL;
        fix_end_iter[i] = NULL;
        fix_rend_iter[i] = NULL;
    }
    calc_mulsizes();
    val.resize(n_elements);
    get_max_shell();
}

template<int n_dim,class T> 
GRID<n_dim,T>::GRID(int *s):n_elements(1),begin_iter(NULL),rbegin_iter(NULL),end_iter(NULL),rend_iter(NULL),
                            shell_begin_iter(NULL), shell_end_iter(NULL) {
    for (int i = 0; i<n_dim; ++i) {
        n_elements *= *(s+i);
        sizes[i] = *(s+i);
        fix_begin_iter[i] = NULL;
        fix_rbegin_iter[i] = NULL;
        fix_end_iter[i] = NULL;
        fix_rend_iter[i] = NULL;
    }
    calc_mulsizes();
    val.resize(n_elements);
    get_max_shell();
}

template<int n_dim,class T> 
GRID<n_dim,T>::GRID(int *s,T standard):n_elements(1),begin_iter(NULL),rbegin_iter(NULL),end_iter(NULL),rend_iter(NULL),
                                       shell_begin_iter(NULL), shell_end_iter(NULL) {
    for (int i = 0; i<n_dim; ++i) {
        n_elements *= *(s+i);
        sizes[i] = *(s+i);
        fix_begin_iter[i] = NULL;
        fix_rbegin_iter[i] = NULL;
        fix_end_iter[i] = NULL;
        fix_rend_iter[i] = NULL;
    }
    calc_mulsizes();
    val.resize(n_elements,standard);
    get_max_shell();
}

template<int n_dim,class T> 
GRID<n_dim,T>::GRID(int s):begin_iter(NULL),rbegin_iter(NULL),end_iter(NULL),rend_iter(NULL),
                           shell_begin_iter(NULL), shell_end_iter(NULL) {
    for (int i = 0; i<n_dim; ++i) {
        sizes[i] = s;
        fix_begin_iter[i] = NULL;
        fix_rbegin_iter[i] = NULL;
        fix_end_iter[i] = NULL;
        fix_rend_iter[i] = NULL;
    }
    calc_mulsizes();
    n_elements = int(pow(float(s),n_dim) + 0.5);
    val.resize(n_elements);
    get_max_shell();
}

template<int n_dim,class T> 
GRID<n_dim,T>::GRID(int s,T standard):begin_iter(NULL),rbegin_iter(NULL),end_iter(NULL),rend_iter(NULL),
                                      shell_begin_iter(NULL), shell_end_iter(NULL) {
    for (int i = 0; i<n_dim; ++i) {
        sizes[i] = s;
        fix_begin_iter[i] = NULL;
        fix_rbegin_iter[i] = NULL;
        fix_end_iter[i] = NULL;
        fix_rend_iter[i] = NULL;
    }
    calc_mulsizes();
    n_elements = int(pow(float(s),n_dim) + 0.5);
    val.resize(n_elements,standard);
    get_max_shell();
}

template<int n_dim,class T>
GRID<n_dim,T>::GRID(GRID<n_dim,T> &ref):begin_iter(NULL),rbegin_iter(NULL),end_iter(NULL),rend_iter(NULL),
                                        shell_begin_iter(NULL), shell_end_iter(NULL) {
    for (int i = 0; i<n_dim; ++i) {
        sizes[i] = ref.sizes[i];
        fix_begin_iter[i] = NULL;
        fix_rbegin_iter[i] = NULL;
        fix_end_iter[i] = NULL;
        fix_rend_iter[i] = NULL;
    }
    calc_mulsizes();
    n_elements = ref.n_elements;
    val = ref.val; //tiefe Kopie
    get_max_shell();
}

template<int n_dim,class T> 
GRID<n_dim,T>::~GRID() {
    if (begin_iter) delete begin_iter;
    if (rbegin_iter) delete rbegin_iter;
    if (end_iter) delete end_iter;
    if (rend_iter) delete rend_iter;
    if (shell_begin_iter) delete shell_begin_iter;
    if (shell_end_iter) delete shell_end_iter;
    for (int i = 0; i<n_dim; ++i) {
        if (fix_begin_iter[i]) delete fix_begin_iter[i];
        if (fix_rbegin_iter[i]) delete fix_rbegin_iter[i];
        if (fix_end_iter[i]) delete fix_end_iter[i];
        if (fix_rend_iter[i]) delete fix_rend_iter[i];
    }
}


//--------------------------------------------------------------------------------------------------
//Operatoren fr GRID
//--------------------------------------------------------------------------------------------------

template<int n_dim,class T> 
GRID<n_dim,T>& GRID<n_dim,T>::operator=(GRID<n_dim,T> &rechts) {
    if (this == &rechts) return *this; //Zeit sparen bei Selbstzuweisung!
    for (int i = 0; i<n_dim; ++i) {sizes[i] = rechts.sizes[i];}
    calc_mulsizes();
    n_elements = rechts.n_elements;
    val = rechts.val;
    max_shell = rechts.max_shell;
    return *this;
}


//--------------------------------------------------------------------------------------------------
//Methoden fr GRID
//--------------------------------------------------------------------------------------------------

template<int n_dim,class T> 
GRIDERATOR<n_dim,T>& GRID<n_dim,T>::begin() {
    if (!begin_iter) begin_iter = new GRIDERATOR<n_dim,T>(this,0);
    return *begin_iter;
}

template<int n_dim,class T> 
GRIDERATOR<n_dim,T>& GRID<n_dim,T>::end() {
    if (!end_iter) end_iter = new GRIDERATOR<n_dim,T>(this,val.size());
    return *end_iter;
}

template<int n_dim,class T> 
GRIDERATOR<n_dim,T>& GRID<n_dim,T>::rbegin() {
    if (!rbegin_iter) rbegin_iter = new GRIDERATOR<n_dim,T>(this,val.size()-1);
    return *rbegin_iter;
}

template<int n_dim,class T> 
GRIDERATOR<n_dim,T>& GRID<n_dim,T>::rend() {
    if (!rend_iter) rend_iter = new GRIDERATOR<n_dim,T>(this,-1);
    return *rend_iter;
}

template<int n_dim,class T> 
FIXED_GRIDERATOR<n_dim,T>& GRID<n_dim,T>::begin(int fixed_dim) {
    if (!fix_begin_iter[fixed_dim]) fix_begin_iter[fixed_dim] = new FIXED_GRIDERATOR<n_dim,T>(this,0,fixed_dim);
    return *fix_begin_iter[fixed_dim];
}

template<int n_dim,class T> 
FIXED_GRIDERATOR<n_dim,T>& GRID<n_dim,T>::end(int fixed_dim) {
    if (!fix_end_iter[fixed_dim]) fix_end_iter[fixed_dim] = new FIXED_GRIDERATOR<n_dim,T>(this,val.size(),fixed_dim);
    return *fix_end_iter[fixed_dim];
}

template<int n_dim,class T> 
FIXED_GRIDERATOR<n_dim,T>& GRID<n_dim,T>::rbegin(int fixed_dim) {
    if (!fix_rbegin_iter[fixed_dim]) fix_rbegin_iter[fixed_dim] = new FIXED_GRIDERATOR<n_dim,T>(this,val.size()-1,fixed_dim);
    return *fix_rbegin_iter[fixed_dim];
}

template<int n_dim,class T> 
FIXED_GRIDERATOR<n_dim,T>& GRID<n_dim,T>::rend(int fixed_dim) {
    if (!fix_rend_iter[fixed_dim]) fix_rend_iter[fixed_dim] = new FIXED_GRIDERATOR<n_dim,T>(this,-1,fixed_dim);
    return *fix_rend_iter[fixed_dim];
}

template<int n_dim,class T> 
SHELL_GRIDERATOR<n_dim,T>& GRID<n_dim,T>::sbegin(int shell) {
    if (shell_begin_iter) {
        if (shell_begin_iter->shell != shell) {
            delete shell_begin_iter;
            shell_begin_iter = new SHELL_GRIDERATOR<n_dim,T>(this,0,shell);
        }
    } else shell_begin_iter = new SHELL_GRIDERATOR<n_dim,T>(this,0,shell);
    return *shell_begin_iter;
}

template<int n_dim,class T> 
SHELL_GRIDERATOR<n_dim,T>& GRID<n_dim,T>::send(int shell) {
    if (shell_end_iter) {
        if (shell_end_iter->shell != shell) {
            delete shell_end_iter;
            shell_end_iter = new SHELL_GRIDERATOR<n_dim,T>(this,val.size(),shell);
        }
    } else shell_end_iter = new SHELL_GRIDERATOR<n_dim,T>(this,val.size(),shell);
    return *shell_end_iter;
}

template<int n_dim,class T>
void GRID<n_dim,T>::calc_mulsizes() {
    mulsizes[0] = sizes[0];
    for (int i = 1; i < n_dim; ++i) {
        mulsizes[i] = mulsizes[i-1]*sizes[i];
    }
}

template<int n_dim,class T>
void GRID<n_dim,T>::get_max_shell() {
    max_shell = sizes[0];
    for (int i = 1; i < n_dim; ++i) {
        if (sizes[i] < max_shell) max_shell = sizes[i];
    }
    max_shell = int((max_shell / 2.) + 0.6) - 1;
}

template<int n_dim,class T>
int* GRID<n_dim,T>::get_mulsizes(int *muls,int *s) {
    muls[0] = s[0];
    for (int i = 1; i < n_dim; ++i) {
        muls[i] = muls[i-1]*s[i];
    }
    return muls;
}

template<int n_dim,class T>
int GRID<n_dim,T>::get_direct_index(int *index) {
    read_index = *index;
    for (int i = 1; i<n_dim; ++i) { //Index berechnen
        read_index += (*(index+i))*mulsizes[i-1];
    }
    return read_index;
}

template<int n_dim,class T>
int* GRID<n_dim,T>::get_coords(int *coord,int index) {
    int help = index;
    int old = index;
    for (int i = n_dim-1; i > 0; --i) {
        help %= mulsizes[i-1];
        coord[i] = (old - help) / mulsizes[i-1];
        old = help;
    }
    coord[0] = help;
    return coord;
}

template<int n_dim,class T>
int GRID<n_dim,T>::get_new_index(int *oldmulsizes,int *newmulsizes,int oldindex) {
    int coord[n_dim];
    int help = oldindex;
    int old = oldindex;
    for (int i = n_dim-1; i > 0; --i) { //alten Index in Koordinaten umrechnen
        help %= oldmulsizes[i-1];
        coord[i] = (old - help) / oldmulsizes[i-1];
        old = help;
    }
    read_index = help;
    for (int i = 1; i<n_dim; ++i) { //neuen Index aus Koordinaten berechnen
        read_index += coord[i]*newmulsizes[i-1];
    }
    return read_index;
}


template<int n_dim,class T>
void GRID<n_dim,T>::clear_iter() {
    //!Da alle Iteratoren nach einem resize ungltig werden mssen die bereits gespeicherten
    //!wieder gel�cht werden
    if (begin_iter) {delete begin_iter; begin_iter = NULL;}
    if (end_iter) {delete end_iter; end_iter = NULL;}
    if (shell_begin_iter) {delete shell_begin_iter; shell_begin_iter = NULL;}
    if (shell_end_iter) {delete shell_end_iter; shell_end_iter = NULL;}
    if (rbegin_iter) {delete rbegin_iter; rbegin_iter = NULL;}
    if (rend_iter) {delete rend_iter; rend_iter = NULL;}
    for (int i = 0; i<n_dim; ++i) {
        if (fix_begin_iter[i]) {delete fix_begin_iter[i]; fix_begin_iter[i] = NULL;}
        if (fix_end_iter[i]) {delete fix_end_iter[i]; fix_end_iter[i] = NULL;}
        if (fix_rbegin_iter[i]) {delete fix_rbegin_iter[i]; fix_rbegin_iter[i] = NULL;}
        if (fix_rend_iter[i]) {delete fix_rend_iter[i]; fix_rend_iter[i] = NULL;}
    }
}

template<int n_dim,class T>
void GRID<n_dim,T>::resize(int *s) {
    n_elements = 1;
    for (int i = 0; i<n_dim; ++i) {
        n_elements *= *(s+i);
        sizes[i] = *(s+i);
    }
    calc_mulsizes();
    val.clear();
    val.resize(n_elements);
    clear_iter();
    get_max_shell();
}

template<int n_dim,class T>
void GRID<n_dim,T>::resize(int s) {
    for (int i = 0; i<n_dim; ++i) {sizes[i] = s;}
    calc_mulsizes();
    n_elements = int(pow(float(s),n_dim) + 0.5);
    val.clear();
    val.resize(n_elements);
    clear_iter();
    get_max_shell();
}

template<int n_dim,class T>
void GRID<n_dim,T>::resize(int *s,T standard) {
    n_elements = 1;
    for (int i = 0; i<n_dim; ++i) {
        n_elements *= *(s+i);
        sizes[i] = *(s+i);
    }
    calc_mulsizes();
    val.clear();
    val.resize(n_elements,standard);
    clear_iter();
    get_max_shell();
}

template<int n_dim,class T>
void GRID<n_dim,T>::resize(int s,T standard) {
    for (int i = 0; i<n_dim; ++i) {sizes[i] = s;}
    calc_mulsizes();
    n_elements = int(pow(float(s),n_dim) + 0.5);
    val.clear();
    val.resize(n_elements,standard);
    clear_iter();
    get_max_shell();
}

template<int n_dim,class T>
void GRID<n_dim,T>::smart_resize(int *newsizes) {
    int oldmulsizes[n_dim];
    int old_n_elements = n_elements;
    n_elements = 1;
    int help = -1;
    for (int i = 0; i < n_dim; ++i) {
        n_elements *= *(newsizes+i);
        oldmulsizes[i] = mulsizes[i];
        if ((sizes[i] != *(newsizes+i))&&(help == -1)) help = i; 
        sizes[i] = *(newsizes+i);
    }
    vector<T> buf(val);
    val.clear();
    val.resize(n_elements);
    calc_mulsizes();
    //!jetzt den alten (internen) index in den neuen index umrechnen und rberkopieren
    if (help < 2) help = 0; else help = oldmulsizes[help-1];
    for (int i = 0; i < old_n_elements; ++i) {
        if (i < sizes[0]) val[i] = buf[i]; //bis hier �dert sich nichts
        else if (i < help) val[i] = buf[i]; //auch bis hier �dert sich nichts (spart viel Zeit, wenn nur eine hohe dim. ge�dert wurde!
        else val[get_new_index(oldmulsizes,mulsizes,i)] = buf[i];
    }
    clear_iter();
    get_max_shell();
}

template<int n_dim,class T>
void GRID<n_dim,T>::smart_resize(int newsize) {
    int oldmulsizes[n_dim];
    int old_n_elements = n_elements;
    int help = -1;
    n_elements = int(pow(float(newsize),n_dim) + 0.5);
    for (int i = 0; i < n_dim; ++i) {
        oldmulsizes[i] = mulsizes[i];
        if ((sizes[i] != newsize)&&(help == -1)) help = i; 
        sizes[i] = newsize;
    }
    vector<T> buf(val);
    val.clear();
    val.resize(n_elements);
    calc_mulsizes();
    if (help < 2) help = 0; else help = oldmulsizes[help-1];
    for (int i = 0; i < old_n_elements; ++i) {
        if (i < sizes[0]) val[i] = buf[i];
        else if (i < help) val[i] = buf[i];
        else val[get_new_index(oldmulsizes,mulsizes,i)] = buf[i];
    }
    clear_iter();
    get_max_shell();
}

template<int n_dim,class T>
void GRID<n_dim,T>::smart_resize(int dim,int newsize) {
    int oldmulsizes[n_dim];
    int help = -1;
    int old_n_elements = n_elements;
    n_elements = newsize*n_elements/sizes[dim];
    sizes[dim] = newsize;
    for (int i = 0; i < n_dim; ++i) {
        oldmulsizes[i] = mulsizes[i]; 
    }
    vector<T> buf(val);
    val.clear();
    val.resize(n_elements);
    calc_mulsizes();
    if (dim < 2) help = 0; else help = oldmulsizes[help-1];
    for (int i = 0; i < old_n_elements; ++i) {
        if (i < sizes[0]) val[i] = buf[i]; //bis hier �dert sich nichts
        else if (i < help) val[i] = buf[i]; //auch bis hier �dert sich nichts (spart viel Zeit, wenn nur eine hohe dim. ge�dert wurde!
        else val[get_new_index(oldmulsizes,mulsizes,i)] = buf[i];
    }
    clear_iter();
    get_max_shell();
}

template<int n_dim,class T>
void GRID<n_dim,T>::smart_resize(int *newsizes,T standard) {
    int oldmulsizes[n_dim];
    int old_n_elements = n_elements;
    n_elements = 1;
    int help = -1;
    for (int i = 0; i < n_dim; ++i) {
        n_elements *= *(newsizes+i);
        oldmulsizes[i] = mulsizes[i];
        if ((sizes[i] != *(newsizes+i))&&(help == -1)) help = i; 
        sizes[i] = *(newsizes+i);
    }
    vector<T> buf(val);
    val.clear();
    val.resize(n_elements);
    calc_mulsizes();
    if (help < 2) help = 0; else help = oldmulsizes[help-1];
    for (int i = 0; i < old_n_elements; ++i) {
        if (i < sizes[0]) val[i] = buf[i]; //bis hier �dert sich nichts
        else if (i < help) val[i] = buf[i]; //auch bis hier �dert sich nichts (spart viel Zeit, wenn nur eine hohe dim. ge�dert wurde!
        else val[get_new_index(oldmulsizes,mulsizes,i)] = buf[i];
    }
    clear_iter();
    get_max_shell();
}

template<int n_dim,class T>
void GRID<n_dim,T>::smart_resize(int dim,int newsize,T standard) {
    int oldmulsizes[n_dim];
    int help = -1;
    int old_n_elements = n_elements;
    n_elements = newsize*n_elements/sizes[dim];
    sizes[dim] = newsize;
    for (int i = 0; i < n_dim; ++i) {
        oldmulsizes[i] = mulsizes[i]; 
    }
    vector<T> buf(val);
    val.clear();
    val.resize(n_elements);
    calc_mulsizes();
    if (dim < 2) help = 0; else help = oldmulsizes[help-1];
    for (int i = 0; i < old_n_elements; ++i) {
        if (i < sizes[0]) val[i] = buf[i]; //bis hier �dert sich nichts
        else if (i < help) val[i] = buf[i]; //auch bis hier �dert sich nichts (spart viel Zeit, wenn nur eine hohe dim. ge�dert wurde!
        else val[get_new_index(oldmulsizes,mulsizes,i)] = buf[i];
    }
    clear_iter();
    get_max_shell();
}

template<int n_dim,class T>
T& GRID<n_dim,T>::value(int *index) { //zum Lesen und Schreiben!
    return val[get_direct_index(index)];
}

template<int n_dim,class T>
T& GRID<n_dim,T>::dvalue(int index) {
    return val[index];
}

template<int n_dim,class T>
T* GRID<n_dim,T>::get_pos(int *index) { //Zeiger auf entspr. Position
    return &(val[get_direct_index(index)]);
}

template<int n_dim,class T>
T* GRID<n_dim,T>::dget_pos(int index) {
    return &(val[index]);
}

template<int n_dim,class T>
T& GRID<n_dim,T>::front() {
    return *(val.begin());
}

template<int n_dim,class T>
T& GRID<n_dim,T>::back() {
    return *(val.rbegin());
}

template<int n_dim,class T>
void GRID<n_dim,T>::push_back(T& value) {
    val.push_back(value);
    sizes[n_dim-1] += 1; //!es wird die letzte Dimension geaendert!!!
    calc_mulsizes();
    n_elements = val.size();
    clear_iter();
    get_max_shell();
}

template<int n_dim,class T>
GRIDINDEXER<n_dim,n_dim-1,T> GRID<n_dim,T>::operator[](int index) {
    read_index = index;
    return GRIDINDEXER<n_dim,n_dim-1,T>(*this);
}

template<int n_dim,class T>
int GRID<n_dim,T>::size() {
    return n_elements;
}

template<int n_dim,class T>
int GRID<n_dim,T>::get_n_elements() {
    return n_elements;
}

template<int n_dim,class T>
int* GRID<n_dim,T>::get_sizes() {
    return sizes;
}

template<int n_dim,class T>
void GRID<n_dim,T>::clear() {
    val.clear();
    val.resize(n_elements);
}

template<int n_dim,class T>
void GRID<n_dim,T>::full_clear() {
    val.resize(0);
}

template<int n_dim,class T>
void GRID<n_dim,T>::set_all(T &standard) {
    val.clear();
    val.resize(n_elements,standard);
}

template<int n_dim,class T>
void GRID<n_dim,T>::save(const char *filename,bool binary) { //!nur fr primitive Typen T
    ofstream f_out;
    if (!binary) { //Textdatei schreiben
        f_out.open(filename);
        if (!f_out) {cerr << c_message<cERROR>("GRID::save --> could not open ") << filename << " for writing!" << endl; return;}
        //zun�hst den Kopf schreiben
        f_out << n_dim << " " << n_elements << " ";
        for (int i=0; i<n_dim; ++i) {
            f_out << sizes[i] << " ";
        }
        for (typename vector<T>::iterator it=val.begin(); it!=val.end(); ++it) {
            f_out << *it << " ";
        }
    } else { //Bin�datei schreiben
        f_out.open(filename,ios::out|ios::binary);
        if (!f_out) {cerr << c_message<cERROR>("GRID::save --> could not open ") << filename << " for writing!" << endl; return;}
        int buf = n_dim;
        f_out.write(reinterpret_cast<char*>(&buf),sizeof(buf));
        f_out.write(reinterpret_cast<char*>(&n_elements),sizeof(n_elements));
        f_out.write(reinterpret_cast<char*>(&(sizes[0])), sizeof(sizes[0])*buf);
        f_out.write(reinterpret_cast<char*>(&(val[0])), sizeof(val[0])*n_elements);
    }
    f_out.close();
}

template<int n_dim,class T>
void GRID<n_dim,T>::load(const char *filename,bool binary) {
    ifstream f_in;
    if (!binary) { //Textdatei lesen
        f_in.open(filename);
        if (!f_in) {cerr << c_message<cERROR>("GRID::load --> could not open ") << filename << " for reading in text-mode!" << endl; return;}
        int buf;
        f_in >> buf;
        if (buf!=n_dim) {
            cerr << c_message<cERROR>("GRID::load --> can not load a grid with dimension ") << buf << " in a grid with dimension " << n_dim << endl;
        }
        f_in >> n_elements;
        for (int i=0; i<n_dim; ++i) {
            f_in >> sizes[i];
        }
        resize(sizes); //!GN:04.08.06
        for (int i=0; i<n_elements; ++i) {
            f_in >> val[i];
        }
    } else { //Bin�datei lesen
        f_in.open(filename,ios::in | ios::binary);
        if (!f_in) {cerr << c_message<cERROR>("GRID::load --> could not open ") << filename << " for reading in binary-mode!" << endl; return;}
        int buf;
        f_in.read(reinterpret_cast<char*>(&buf),sizeof(buf));
        if (buf!=n_dim) {
            cerr << c_message<cERROR>("GRID::load --> can not load a grid with dimension ") << buf << " in a grid with dimension " << n_dim << endl;
        } else
        f_in.read(reinterpret_cast<char*>(&n_elements),sizeof(n_elements));
        f_in.read(reinterpret_cast<char*>(&(sizes[0])), sizeof(sizes[0])*buf);
        resize(sizes); //!GN:04.08.06
        f_in.read(reinterpret_cast<char*>(&(val[0])), sizeof(val[0])*n_elements);
    }
    f_in.close();
}


//==================================================================================================
//!Definitionen der Klasse GRIDERATOR:
//==================================================================================================

template<int n_dim,class T> 
GRIDERATOR<n_dim,T>::GRIDERATOR(GRID<n_dim,T> *grd,int i):grid(grd),index(i) {}

template<int n_dim,class T> 
GRIDERATOR<n_dim,T>::GRIDERATOR(GRID<n_dim,T> *grd,int *i):grid(grd) {
    index = grid->get_direct_index(i);
}

template<int n_dim,class T> 
GRIDERATOR<n_dim,T>::~GRIDERATOR() {}


//--------------------------------------------------------------------------------------------------
//Operatoren fr GRIDERATOR
//--------------------------------------------------------------------------------------------------

template<int n_dim,class T> 
GRIDERATOR<n_dim,T>& GRIDERATOR<n_dim,T>::operator=(GRIDERATOR<n_dim,T> &rechts) {
    if (this == &rechts) return *this; //Zeit sparen bei Selbstzuweisung!
    grid = rechts.grid;
    index = rechts.index;
    return *this;
}

template<int n_dim,class T> 
GRIDERATOR<n_dim,T>& GRIDERATOR<n_dim,T>::operator=(int *i) {
    index = grid->get_direct_index(i);
    return *this;
}

template<int n_dim,class T> 
T& GRIDERATOR<n_dim,T>::operator*() {
    return grid->dvalue(index);
}

template<int n_dim,class T> 
T* GRIDERATOR<n_dim,T>::operator->() {
    return grid->dget_pos(index);
}

template<int n_dim,class T> 
int GRIDERATOR<n_dim,T>::get_index() {
    return index;
}

template<int n_dim,class T>
bool GRIDERATOR<n_dim,T>::operator==(GRIDERATOR<n_dim,T> rechts) {
    if ((index == rechts.index)&&(grid == rechts.grid)) return true;
    return false;
}

template<int n_dim,class T>
bool GRIDERATOR<n_dim,T>::operator!=(GRIDERATOR<n_dim,T> rechts) {
    if ((index == rechts.index)&&(grid == rechts.grid)) return false;
    return true;
}

template<int n_dim,class T>
GRIDERATOR<n_dim,T>& GRIDERATOR<n_dim,T>::operator++() {
    ++index;
    return *this;
}

template<int n_dim,class T>
GRIDERATOR<n_dim,T> GRIDERATOR<n_dim,T>::operator++(int) {
    ++index;
    return *this;
}

template<int n_dim,class T>
GRIDERATOR<n_dim,T>& GRIDERATOR<n_dim,T>::operator--() {
    --index;
    return *this;
}

template<int n_dim,class T>
GRIDERATOR<n_dim,T> GRIDERATOR<n_dim,T>::operator--(int) {
    --index;
    return *this;
}

template<int n_dim,class T>
GRIDERATOR<n_dim,T>& GRIDERATOR<n_dim,T>::operator+=(int inc) {
    index += inc;
    return *this;
}

template<int n_dim,class T>
GRIDERATOR<n_dim,T>& GRIDERATOR<n_dim,T>::operator-=(int dec) {
    index -= dec;
    return *this;
}


//==================================================================================================
//!Definitionen der Klasse FIXED_GRIDERATOR:
//==================================================================================================

template<int n_dim,class T> 
FIXED_GRIDERATOR<n_dim,T>::FIXED_GRIDERATOR(GRID<n_dim,T> *grd,int i,int fixed_dim) : base(grd,i),
                                        dim(fixed_dim),n_step(1) {
    for (int j=0; j<dim; ++j) {
        n_step *= grd->sizes[j];
    }
    if (dim == 0) cur_step = base::index;
    else cur_step = base::index % n_step;
    n_jump = n_step * grd->sizes[dim];
    n_step -= 1; //gibt an wie oft um eins hochgez�lt wird bevor ein Sprung erfolgt
    n_jump -= n_step; //gibt die Sprungweite bis zur n�hsten Sequenz
    if (cur_step > n_step) cur_step = 0;
}

template<int n_dim,class T> 
FIXED_GRIDERATOR<n_dim,T>::FIXED_GRIDERATOR(GRID<n_dim,T> *grd,int *index,int fixed_dim) : base(grd,index),
                                            dim(fixed_dim),n_step(1) {
//    base::grid = grd;
//    base::index = grd->get_direct_index(index);
    for (int j=0; j<dim; ++j) {
        n_step *= grd->sizes[j];
    }
    if (dim == 0) cur_step = base::index;
    else cur_step = base::index % n_step;
    n_jump = n_step * grd->sizes[dim];
    n_step -= 1; //gibt an wie oft um eins hochgez�lt wird bevor ein Sprung erfolgt
    n_jump -= n_step; //gibt die Sprungweite bis zur n�hsten Sequenz
    if (cur_step > n_step) cur_step = 0;
}

template<int n_dim,class T> 
FIXED_GRIDERATOR<n_dim,T>::~FIXED_GRIDERATOR() {}


//--------------------------------------------------------------------------------------------------
//Operatoren fr FIXED_GRIDERATOR
//--------------------------------------------------------------------------------------------------

template<int n_dim,class T> 
FIXED_GRIDERATOR<n_dim,T>& FIXED_GRIDERATOR<n_dim,T>::operator=(FIXED_GRIDERATOR<n_dim,T> &rechts) {
    base::grid = rechts.grid;
    base::index = rechts.index;
    dim = rechts.dim;
    n_step = rechts.n_step;
    n_jump = rechts.n_jump;
    cur_step = rechts.cur_step;
    return *this;
}

template<int n_dim,class T> 
FIXED_GRIDERATOR<n_dim,T>& FIXED_GRIDERATOR<n_dim,T>::operator=(int *i) {
    base::index = base::grid->get_direct_index(i);
    cur_step = base::index % (n_step+1);
    return *this;
}

template<int n_dim,class T>
bool FIXED_GRIDERATOR<n_dim,T>::operator==(FIXED_GRIDERATOR<n_dim,T> rechts) {
    if (base::grid == rechts.grid) {
        if (base::index == rechts.index) return true;
        if ((rechts.index == int(base::grid->val.size()))&&(base::index > rechts.index)) {
            return true; //index liegt bereits hinter dem Ende (fr ==end(dim) -Abfrage)
        }
        if ((rechts.index == -1)&&(base::index < rechts.index)) {
            return true; //index liegt bereits vor dem rEnde (fr ==rend(dim) -Abfrage)
        }
    }
    return false;
}

template<int n_dim,class T>
bool FIXED_GRIDERATOR<n_dim,T>::operator!=(FIXED_GRIDERATOR<n_dim,T> rechts) {
    if ((base::index == rechts.index)&&(base::grid == rechts.grid)) return false;
    if ((rechts.index == int(base::grid->val.size()))&&(base::index > rechts.index)) {
        return false; //index liegt bereits hinter dem Ende (fr !=end(dim) -Abfrage)
    }
    if ((rechts.index == -1)&&(base::index < rechts.index)) {
        return false; //index liegt bereits vor dem rEnde (fr !=rend(dim) -Abfrage)
    }
    return true;
}


template<int n_dim,class T>
FIXED_GRIDERATOR<n_dim,T>& FIXED_GRIDERATOR<n_dim,T>::operator++() {
    if (cur_step == n_step) {
        base::index += n_jump;
        cur_step = 0;
    } else {
        base::index += 1;
        ++cur_step;
    }
    return *this;
}

template<int n_dim,class T>
FIXED_GRIDERATOR<n_dim,T> FIXED_GRIDERATOR<n_dim,T>::operator++(int) {
    if (cur_step == n_step) {
        base::index += n_jump;
        cur_step = 0;
    } else {
        base::index += 1;
        ++cur_step;
    }
    return *this;
}

template<int n_dim,class T>
FIXED_GRIDERATOR<n_dim,T>& FIXED_GRIDERATOR<n_dim,T>::operator--() {
    if (cur_step == 0) {
        base::index -= n_jump;
        cur_step = n_step;
    } else {
        base::index -= 1;
        --cur_step;
    }
    return *this;
}

template<int n_dim,class T>
FIXED_GRIDERATOR<n_dim,T> FIXED_GRIDERATOR<n_dim,T>::operator--(int) {
    if (cur_step == 0) {
        base::index -= n_jump;
        cur_step = n_step;
    } else {
        base::index -= 1;
        --cur_step;
    }
    return *this;
}

template<int n_dim,class T>
FIXED_GRIDERATOR<n_dim,T>& FIXED_GRIDERATOR<n_dim,T>::operator+=(int inc) {
    for (int i=0; i<inc; ++i) ++*this;
    return *this;
}

template<int n_dim,class T>
FIXED_GRIDERATOR<n_dim,T>& FIXED_GRIDERATOR<n_dim,T>::operator-=(int dec) {
    for (int i=0; i<dec; ++i) --*this;
    return *this;
}



//==================================================================================================
//!Definitionen der Klasse SHELL_GRIDERATOR:
//==================================================================================================

template<int n_dim,class T> 
SHELL_GRIDERATOR<n_dim,T>::SHELL_GRIDERATOR(GRID<n_dim,T> *grd,int i,int sh) : base(grd,i),
                                        shell(sh) {
    for (int j=0; j<n_dim; ++j) {
        op_ind[j] = grd->sizes[j] - shell - 1;
    }
    base::index--;
    ++*this; //um den naechsten gueltigen Index zu bekommen (falls i nicht zur shell gehoert)
}

template<int n_dim,class T> 
SHELL_GRIDERATOR<n_dim,T>::SHELL_GRIDERATOR(GRID<n_dim,T> *grd,int *index,int sh) : base(grd,index),
                                            shell(sh) {
    for (int j=0; j<n_dim; ++j) {
        op_ind[j] = grd->sizes[j] - shell - 1;
    }
    base::index--;
    ++*this;
}

template<int n_dim,class T> 
SHELL_GRIDERATOR<n_dim,T>::~SHELL_GRIDERATOR() {}


//--------------------------------------------------------------------------------------------------
//Operatoren fr SHELL_GRIDERATOR
//--------------------------------------------------------------------------------------------------

template<int n_dim,class T> 
SHELL_GRIDERATOR<n_dim,T>& SHELL_GRIDERATOR<n_dim,T>::operator=(SHELL_GRIDERATOR<n_dim,T> &rechts) {
    base::grid = rechts.grid;
    base::index = rechts.index;
    shell = rechts.shell;
    for (int j=0; j<n_dim; ++j) {
        op_ind[j] = rechts.op_ind[j];
    }
    return *this;
}

template<int n_dim,class T> 
SHELL_GRIDERATOR<n_dim,T>& SHELL_GRIDERATOR<n_dim,T>::operator=(int *i) {
    base::index = base::grid->get_direct_index(i);
    return *this;
}

template<int n_dim,class T>
bool SHELL_GRIDERATOR<n_dim,T>::operator==(SHELL_GRIDERATOR<n_dim,T> rechts) {
    if (base::grid == rechts.grid) {
        if (base::index == rechts.index) return true;
        if ((rechts.index == int(base::grid->val.size()))&&(base::index > rechts.index)) {
            return true; //index liegt bereits hinter dem Ende (fr ==end(dim) -Abfrage)
        }
        if ((rechts.index == -1)&&(base::index < rechts.index)) {
            return true; //index liegt bereits vor dem rEnde (fr ==rend(dim) -Abfrage)
        }
    }
    return false;
}

template<int n_dim,class T>
bool SHELL_GRIDERATOR<n_dim,T>::operator!=(SHELL_GRIDERATOR<n_dim,T> rechts) {
    if ((base::index == rechts.index)&&(base::grid == rechts.grid)) return false;
    if ((rechts.index == int(base::grid->val.size()))&&(base::index > rechts.index)) {
        return false; //index liegt bereits hinter dem Ende (fr !=end(dim) -Abfrage)
    }
    if ((rechts.index == -1)&&(base::index < rechts.index)) {
        return false; //index liegt bereits vor dem rEnde (fr !=rend(dim) -Abfrage)
    }
    return true;
}

template<int n_dim,class T>
SHELL_GRIDERATOR<n_dim,T>& SHELL_GRIDERATOR<n_dim,T>::operator++() {
    bool breaker;
    bool conti;
    int tester;    
    do {
        base::index++;
        int help = base::index;
        int old = base::index;
        breaker = false;
        conti = false;
        for (int i = n_dim-1; i > 0; --i) {
            help %= base::grid->mulsizes[i-1];
            tester = (old - help) / base::grid->mulsizes[i-1];
            if ((tester < shell) || (tester > op_ind[i])) {
                conti = true;
            }
            if ((tester == shell) || (tester == op_ind[i])) {
                breaker = true;
            }
            old = help;
        }
        if ((help < shell) || (help > op_ind[0])) conti = true;
        if ((help == shell) || (help == op_ind[0])) breaker = true;
        if (conti) continue;
        if (breaker) break;
    } while(base::index < int(base::grid->val.size()));
    return *this;
}

template<int n_dim,class T>
SHELL_GRIDERATOR<n_dim,T> SHELL_GRIDERATOR<n_dim,T>::operator++(int) {
    bool breaker;
    bool conti;
    int tester;    
    do {
        base::index++;
        int help = base::index;
        int old = base::index;
        breaker = false;
        conti = false;
        for (int i = n_dim-1; i > 0; --i) {
            help %= base::grid->mulsizes[i-1];
            tester = (old - help) / base::grid->mulsizes[i-1];
            if ((tester < shell) || (tester > op_ind[i])) {
                conti = true;
            }
            if ((tester == shell) || (tester == op_ind[i])) {
                breaker = true;
            }
            old = help;
        }
        if ((help < shell) || (help > op_ind[0])) conti = true;
        if ((help == shell) || (help == op_ind[0])) breaker = true;
        if (conti) continue;
        if (breaker) break;
    } while(base::index < int(base::grid->val.size()));
    return *this;
}

template<int n_dim,class T>
SHELL_GRIDERATOR<n_dim,T>& SHELL_GRIDERATOR<n_dim,T>::operator--() {
    bool breaker;
    bool conti;
    int tester;    
    do {
        base::index--;
        int help = base::index;
        int old = base::index;
        breaker = false;
        conti = false;
        for (int i = n_dim-1; i > 0; --i) {
            help %= base::grid->mulsizes[i-1];
            tester = (old - help) / base::grid->mulsizes[i-1];
            if ((tester < shell) || (tester > op_ind[i])) {
                conti = true;
            }
            if ((tester == shell) || (tester == op_ind[i])) {
                breaker = true;
            }
            old = help;
        }
        if ((help < shell) || (help > op_ind[0])) conti = true;
        if ((help == shell) || (help == op_ind[0])) breaker = true;
        if (conti) continue;
        if (breaker) break;
    } while(base::index < int(base::grid->val.size()));
    return *this;
}

template<int n_dim,class T>
SHELL_GRIDERATOR<n_dim,T> SHELL_GRIDERATOR<n_dim,T>::operator--(int) {
    bool breaker;
    bool conti;
    int tester;    
    do {
        base::index--;
        int help = base::index;
        int old = base::index;
        breaker = false;
        conti = false;
        for (int i = n_dim-1; i > 0; --i) {
            help %= base::grid->mulsizes[i-1];
            tester = (old - help) / base::grid->mulsizes[i-1];
            if ((tester < shell) || (tester > op_ind[i])) {
                conti = true;
            }
            if ((tester == shell) || (tester == op_ind[i])) {
                breaker = true;
            }
            old = help;
        }
        if ((help < shell) || (help > op_ind[0])) conti = true;
        if ((help == shell) || (help == op_ind[0])) breaker = true;
        if (conti) continue;
        if (breaker) break;
    } while(base::index < int(base::grid->val.size()));
    return *this;
}

template<int n_dim,class T>
SHELL_GRIDERATOR<n_dim,T>& SHELL_GRIDERATOR<n_dim,T>::operator+=(int inc) {
    for (int i=0; i<inc; ++i) ++*this;
    return *this;
}

template<int n_dim,class T>
SHELL_GRIDERATOR<n_dim,T>& SHELL_GRIDERATOR<n_dim,T>::operator-=(int dec) {
    for (int i=0; i<dec; ++i) --*this;
    return *this;
}



//!bisher nur FIXED_GRIDERATOR kopiert und das FIXED durch SUB ersetzt
/*
//==================================================================================================
//!Definitionen der Klasse SUB_GRIDERATOR:
//==================================================================================================

template<int n_dim,class T> 
SUB_GRIDERATOR<n_dim,T>::SUB_GRIDERATOR(GRID<n_dim,T> *grd,int i_start,int i_end) : base(grd,i_start) {
    n_step = 1;
    for (int j=0; j<dim; ++j) {
        n_step *= grd->sizes[j];
    }
    if (dim == 0) cur_step = base::index;
    else cur_step = base::index % n_step;
    n_jump = n_step * grd->sizes[dim];
    n_step -= 1; //gibt an wie oft um eins hochgez�lt wird bevor ein Sprung erfolgt
    n_jump -= n_step; //gibt die Sprungweite bis zur n�hsten Sequenz
}

template<int n_dim,class T> 
SUB_GRIDERATOR<n_dim,T>::SUB_GRIDERATOR(GRID<n_dim,T> *grd,int *index,int fixed_dim) {
    base::grid = grd;
    base::index = grd->get_direct_index(index);
    dim = fixed_dim;
    n_step = 1;
    for (int j=0; j<dim; ++j) {
        n_step *= grd->sizes[j];
    }
    if (dim == 0) cur_step = base::index;
    else cur_step = base::index % n_step;
    n_jump = n_step * grd->sizes[dim];
    n_step -= 1; //gibt an wie oft um eins hochgez�lt wird bevor ein Sprung erfolgt
    n_jump -= n_step; //gibt die Sprungweite bis zur n�hsten Sequenz
}

template<int n_dim,class T> 
SUB_GRIDERATOR<n_dim,T>::~SUB_GRIDERATOR() {}


//--------------------------------------------------------------------------------------------------
//Operatoren fr SUB_GRIDERATOR
//--------------------------------------------------------------------------------------------------

template<int n_dim,class T> 
SUB_GRIDERATOR<n_dim,T>& SUB_GRIDERATOR<n_dim,T>::operator=(SUB_GRIDERATOR<n_dim,T> &rechts) {
    base::grid = rechts.grid;
    base::index = rechts.index;
    dim = rechts.dim;
    n_step = rechts.n_step;
    n_jump = rechts.n_jump;
    cur_step = rechts.cur_step;
    return *this;
}

template<int n_dim,class T> 
SUB_GRIDERATOR<n_dim,T>& SUB_GRIDERATOR<n_dim,T>::operator=(int *i) {
    base::index = base::grid->get_direct_index(i);
    cur_step = base::index % (n_step+1);
    return *this;
}

template<int n_dim,class T>
bool SUB_GRIDERATOR<n_dim,T>::operator==(SUB_GRIDERATOR<n_dim,T> rechts) {
    if (base::grid == rechts.grid) {
        if (base::index == rechts.index) return true;
        if ((rechts.index == base::grid->val.size())&&(base::index > rechts.index)) {
            return true; //index liegt bereits hinter dem Ende (fr ==end(dim) -Abfrage)
        }
        if ((rechts.index == -1)&&(base::index < rechts.index)) {
            return true; //index liegt bereits vor dem rEnde (fr ==rend(dim) -Abfrage)
        }
    }
    return false;
}

template<int n_dim,class T>
bool SUB_GRIDERATOR<n_dim,T>::operator!=(SUB_GRIDERATOR<n_dim,T> rechts) {
    if ((base::index == rechts.index)&&(base::grid == rechts.grid)) return false;
    if ((rechts.index == base::grid->val.size())&&(base::index > rechts.index)) {
        return false; //index liegt bereits hinter dem Ende (fr !=end(dim) -Abfrage)
    }
    if ((rechts.index == -1)&&(base::index < rechts.index)) {
        return false; //index liegt bereits vor dem rEnde (fr !=rend(dim) -Abfrage)
    }
    return true;
}


template<int n_dim,class T>
SUB_GRIDERATOR<n_dim,T>& SUB_GRIDERATOR<n_dim,T>::operator++() {
    if (cur_step == n_step) {
        base::index += n_jump;
        cur_step = 0;
    } else {
        base::index += 1;
        ++cur_step;
    }
    return *this;
}

template<int n_dim,class T>
SUB_GRIDERATOR<n_dim,T> SUB_GRIDERATOR<n_dim,T>::operator++(int) {
    if (cur_step == n_step) {
        base::index += n_jump;
        cur_step = 0;
    } else {
        base::index += 1;
        ++cur_step;
    }
    return *this;
}

template<int n_dim,class T>
SUB_GRIDERATOR<n_dim,T>& SUB_GRIDERATOR<n_dim,T>::operator--() {
    if (cur_step == 0) {
        base::index -= n_jump;
        cur_step = n_step;
    } else {
        base::index -= 1;
        --cur_step;
    }
    return *this;
}

template<int n_dim,class T>
SUB_GRIDERATOR<n_dim,T> SUB_GRIDERATOR<n_dim,T>::operator--(int) {
    if (cur_step == 0) {
        base::index -= n_jump;
        cur_step = n_step;
    } else {
        base::index -= 1;
        --cur_step;
    }
    return *this;
}

template<int n_dim,class T>
SUB_GRIDERATOR<n_dim,T>& SUB_GRIDERATOR<n_dim,T>::operator+=(int inc) {
    for (int i=0; i<inc; ++i) ++*this;
    return *this;
}

template<int n_dim,class T>
SUB_GRIDERATOR<n_dim,T>& SUB_GRIDERATOR<n_dim,T>::operator-=(int dec) {
    for (int i=0; i<dec; ++i) --*this;
    return *this;
}
*/



//==================================================================================================
//!Definitionen der Klasse GRIDINDEXER:
//==================================================================================================

template<int n_dim,int cur_dim,class T> 
GRIDINDEXER<n_dim,cur_dim,T>::GRIDINDEXER(GRID<n_dim,T> &grd) : grid(grd) {}

template<int n_dim,int cur_dim,class T> 
GRIDINDEXER<n_dim,cur_dim,T>::~GRIDINDEXER() {}


//--------------------------------------------------------------------------------------------------
//Operatoren fr GRIDINDEXER
//--------------------------------------------------------------------------------------------------

template<int n_dim,int cur_dim,class T>
GRIDINDEXER<n_dim,cur_dim-1,T> GRIDINDEXER<n_dim,cur_dim,T>::operator[](int index) {
    grid.read_index += index*grid.mulsizes[n_dim-cur_dim-1];
    return GRIDINDEXER<n_dim,cur_dim-1,T>(grid);
}


//--------------------------------------------------------------------------------------------------
//Spezialisierung von GRIDINDEXER fr cur_dim == 1
//--------------------------------------------------------------------------------------------------

template<int n_dim,class T> 
GRIDINDEXER<n_dim,1,T>::GRIDINDEXER(GRID<n_dim,T> &grd) : grid(grd) {}

template<int n_dim,class T> 
GRIDINDEXER<n_dim,1,T>::~GRIDINDEXER() {}

template<int n_dim,class T>
T& GRIDINDEXER<n_dim,1,T>::operator[](int index) {
    grid.read_index += index*grid.mulsizes[n_dim-2];
    return grid.val[grid.read_index];
}


#endif //__GRIDGN

