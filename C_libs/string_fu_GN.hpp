
//============================================================================
// string_fu_GN.hpp -*- C++ -*-; some string manipulations
//
// Copyright (C) 2007 Gerd Neudert
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
// This library implements some string manipulations, but is nearly obsolete
// as many of those functions are not safe and even not effective. I had just
// no time to reimplement this.
//============================================================================


#ifndef __STRINGFUGN
#define __STRINGFUGN
#include<string>
#include<sstream>
#include<vector>
#include<iostream>

using namespace std;

namespace string_fu {

//!============================================================================================
//!s2v(string,zahl):
//!Wandelt den string in eine Zahl um und
//!liefert bei erfolgreicher Umwandlung 'true' zurueck (sonst 'false'):
template<typename T>
inline bool s2v(string const &s, T &v) {
    istringstream is;
    is.str(s); is >> v;
    if (is.fail()) return false;
    else return true;
}
//!============================================================================================


//!============================================================================================
//!add2s(string,T):
//!Konkateniert einen beliebigen Typ an den string
//!Absichtlich keine Fehlerabfrage!!!
template<typename T>
inline void add2s(string &s,T const &v) {
    ostringstream os;
    os << s << v;
    s = os.str();
}
//!============================================================================================


//!============================================================================================
//!mysplit(string,string_vector,trennzeichen):
//!Teilt den uebergebenen String in durch das
//!Trennzeichen getrennte Teilstrings auf und packt diese in den uebergebenen Vector.
//!Zurueckgegeben wird die Anzahl der Teilstrings. Ist das Trennzeichen nicht in dem String
//!enthalten, so enthaelt der Vector den gesamten String als einziges Element. Besteht der
//!String NUR aus dem Trennzeichen, so wird Null zurueckgeliefert!
inline int mysplit(string const &s,vector<string> &vec,char const split_char = ' ') {
    //vec.clear(); //besser dem user ueberlassen!!!
    string tmp;
    string::size_type i1 = s.find(split_char);
    if (i1 == string::npos) {
        vec.push_back(s);
        return 1;
    }
    tmp = s.substr(0,i1);
    if (tmp.size() > 0) vec.push_back(tmp);
    i1 += 1;
    string::size_type i2;
    while ((i2 = s.find(split_char,i1)) != string::npos) {
        tmp = s.substr(i1,(i2-i1));
        if (tmp.size() > 0) vec.push_back(tmp);
        i1 = i2 + 1;
    }
    if (s.size() >= i1) {
        tmp = s.substr(i1);
        if (tmp.size() > 0) vec.push_back(tmp);
    }
    return vec.size();
}
//!============================================================================================


//!============================================================================================
//!remove_char(string,char):
//!Entfernt ein bestimmtes Zeichen aus einem String
//!Absichtlich keine Fehlerabfrage!!!
inline void remove_char(string &s,char const rc = ' ') {
    vector<string> tv;
    mysplit(s,tv,rc);
    ostringstream os;
    for (vector<string>::iterator it=tv.begin(); it!=tv.end(); ++it) {
        os << *it;
    }
    s = os.str();
}
//!============================================================================================


//!============================================================================================
//!replace_char(string,old_char,new_char):
//!Ersetzt ein bestimmtes Zeichen in einem String
//!durch ein anderes. ACHTUNG: Zeichen am Anfang und Ende des Strings werden NICHT ersetzt!!!
//!Absichtlich keine Fehlerabfrage!!!
inline void replace_char(string &s,char const c_old,char const c_new) {
    vector<string> tv;
    mysplit(s,tv,c_old);
    ostringstream os;
    bool nf = false;
    for (vector<string>::iterator it=tv.begin(); it!=tv.end(); ++it) {
        if (nf) os << c_new;
        else nf =true;
        os << *it;
    }
    s = os.str();
}
//!============================================================================================


//!============================================================================================
//!get_ext(string[,char]): 
//!Liefert die Dateierweiterung des uebergebenen Strings, also den letzten
//!Punkt im String und das was folgt (der Punkt wird also mit zurueckgeliefert!!!)
//!Alternativ kann auch ein anderes Zeichen als der Punkt mitgegeben werden, um die
//!Funktion allgemein zu halten.
inline string get_ext(string const &name,char lc = '.') {
    string::size_type posi = name.rfind(lc);
    if (posi == string::npos) return "";
    else return name.substr(posi);
}
//!============================================================================================


//!============================================================================================
//!get_text_after(string[,char]): 
//!Liefert alles was nach dem letzten Vorkommen des gegebenen Zeichens kommt.
//!Liefert den kompletten String, wenn das Zeichen nicht enthalten ist.
inline string get_text_after(string const &name,char lc = '.') {
    string::size_type posi = name.rfind(lc);
    if (posi == string::npos) return name;
    posi += 1;
    if (posi == string::npos) return name;
    else return name.substr(posi);
}
//!============================================================================================


//!============================================================================================
//!replace_ext(target_string,current_ext,new_ext):
//!Ersetzt das Ende des Strings (als current_ext uebergeben) durch eine neue Endung
//!(als new_ext uebergeben). Wie der Name schon sagt, ist die Funktion urspruenglich
//!dazu gedacht Dateierweiterungen zu aendern.
//!ACHTUNG!!!: Es wird nicht geprueft, ob der String wirklich die angegebene
//!Endung hat! Diese wird lediglich benutzt, um zu ermitteln, wieviel vom
//!uebergebenen String abgeschnitten werden muss!
//!current_ext kann natuerlich auch als get_ext(target_string) uebergeben werden, wobei
//!dann zu beachten ist, dass new_ext auch den Punkt (z.b. ".mol2") enthalten muss.
inline void replace_ext(string &tfile,string const &ext,string const &by_ext) {
    string v_file(tfile,0,tfile.size()-ext.size());
    tfile = v_file + by_ext;
}
//!============================================================================================


//!============================================================================================
//!remove_ext(target_string[,char]):
//!Schneidet den String vor dem letzten vorkommen von char ab.
inline void remove_ext(string &name,char lc = '.') {
    replace_ext(name,get_ext(name,lc),"");
}
//!============================================================================================


//!============================================================================================
//!replace_substr(target_string,substr,replace_string):
inline void replace_substr(string &name,string const &substr,string const &newsubstr) {
    string::size_type posi = name.find(substr);
    name.replace(posi,posi+substr.size()-1,newsubstr);
}
//!============================================================================================


//!============================================================================================
//!regex_compare(regex,reference):
//!Unterstuetzt bisher lediglich '*' als Wildcard !!!
inline bool regex_compare(string &regex,string const &ref) {
    string::size_type ref_von = 0;
    string::size_type von = 0;
    string::size_type bis = 0;
    string::size_type tmp;
    
    von = regex.find("*");
    if (von == string::npos) {
        if (regex == ref) return true;
        else return false;
    }
    if (von > 0) { //! wenn keine Wildcard am Anfang steht
        string sub = regex.substr(0,von);
        tmp = ref.find(sub,ref_von);
        if (tmp != ref_von) return false;
        ref_von = tmp + sub.size();
    }
    
    while (ref_von < ref.size()) {
        von = regex.find("*",von);
        if (von == string::npos) if (ref_von != ref.size()) return false; //!fehlt ein abschliessendes '*'
        bis = regex.find("*",von+1);
        
        string sub;
        if (bis == string::npos) {
            if (von == regex.size()-1) return true;
            sub = regex.substr(von+1);
            tmp = ref.find(sub,ref_von);
            if (tmp == string::npos) return false;
        } else {
            sub = regex.substr(von+1,bis-von-1);
            tmp = ref.find(sub,ref_von);
            if (tmp == string::npos) return false;
        }
        
        ref_von = tmp + sub.size();
        von = bis;
    }
    return true;
}
//!============================================================================================

} //Ende Namespace

#endif //__STRINGFUGN
