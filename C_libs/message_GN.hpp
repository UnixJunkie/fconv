
//============================================================================
// message_GN.hpp -*- C++ -*-; vector and matrix calculations
//
// Copyright (C) 2006, 2007, 2008, 2009 Gerd Neudert
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
// This library provides some standard output functions. If you compile on
// Windows or MacOS, you should compile without '-D_LINUX_OS'
//============================================================================


#include<sstream>
#include<string>

#ifndef __MESSAGEGN
#define __MESSAGEGN

//!nuetzliche Escapesequenzen, die bisher nicht benutzt werden:
//! \033[1;31m  => rot in fett (funktioniert aber nicht richtig - Farbe wird stattdessen heller
//! \033[4;31m  => rot und Unterstrichen
//! \033[5;31m  => rot und blinkend
//! ( \033[0;31m   => standard rot)

using namespace std;

namespace TXT {
namespace { //!Muss nochmal in einen anonymen Namensbereich, weil der Compiler sonst wegen multiplen
            //!Deklarationen der Konstanten meckert!!!
    #if defined (_LINUX_OS)
    //Textfarben
    extern char const cNORM[] = "\033[0m";   //!weil intern gelinkte Objekte nicht als Template-Parameter zulaessig sind
    extern char const cBLACK[] = "\033[30m"; //!soll sich im zukuenftigen C-Standard aendern => dann wieder const char* name = ""
    extern char const cRED[] = "\033[31m";
    extern char const cGREEN[] = "\033[32m";
    extern char const cYELLOW[] = "\033[33m";
    extern char const cBLUE[] = "\033[34m";
    extern char const cMAGENTA[] = "\033[35m";
    extern char const cCYAN[] = "\033[36m";
    extern char const cGREY[] = "\033[37m";
    
    //Backgroundfarben
    extern char const bgcNORM[] = "\033[0m";
    extern char const bgcBLACK[] = "\033[40m";
    extern char const bgcRED[] = "\033[41m";
    extern char const bgcGREEN[] = "\033[42m";
    extern char const bgcYELLOW[] = "\033[43m";
    extern char const bgcBLUE[] = "\033[44m";
    extern char const bgcMAGENTA[] = "\033[45m";
    extern char const bgcCYAN[] = "\033[46m";
    extern char const bgcGREY[] = "\033[47m";
    
    //einige Standardtextfragmente:
    extern char const cERROR[] = "\033[31merror: ";
    extern char const cWARNING[] = "\033[33mwarning: ";
    
    //oder als Funktionen:
    template<const char* T>
    inline const char* c_message(const char* msg) {
        ostringstream os;
        static string res;
        os << T << msg << "\033[0m";
        res = os.str();
        return res.c_str();
    }
    //Spezialisierung fuer errors
    template<>
    inline const char* c_message<cERROR>(const char* msg) {
        ostringstream os;
        static string res;
        os << cERROR << "\033[0m" << msg;
        res = os.str();
        return res.c_str();
    }
    //Spezialisierung fuer warnings
    template<>
    inline const char* c_message<cWARNING>(const char* msg) {
        ostringstream os;
        static string res;
        os << cWARNING << "\033[0m" << msg;
        res = os.str();
        return res.c_str();
    }
    //Fuer strings:
    template<const char* T>
    inline const char* c_message(string &msg) {
        ostringstream os;
        static string res;
        os << T << msg << "\033[0m";
        res = os.str();
        return res.c_str();
    }
    //Spezialisierung fuer errors
    template<>
    inline const char* c_message<cERROR>(string &msg) {
        ostringstream os;
        static string res;
        os << cERROR << "\033[0m" << msg;
        res = os.str();
        return res.c_str();
    }
    //Spezialisierung fuer warnings
    template<>
    inline const char* c_message<cWARNING>(string &msg) {
        ostringstream os;
        static string res;
        os << cWARNING << "\033[0m" << msg;
        res = os.str();
        return res.c_str();
    }
    
    
    //Cursor bewegen:
    inline const char* c_left(int n = 1) {
        ostringstream os;
        static string res;
        os << "\033[" << n << "D";
        res = os.str();
        return res.c_str();
    }
    inline const char* c_right(int n = 1) {
        ostringstream os;
        static string res;
        os << "\033[" << n << "C";
        res = os.str();
        return res.c_str();
    }
    inline const char* c_up(int n = 1) {
        ostringstream os;
        static string res;
        os << "\033[" << n << "A";
        res = os.str();
        return res.c_str();
    }
    inline const char* c_down(int n = 1) {
        ostringstream os;
        static string res;
        os << "\033[" << n << "B";
        res = os.str();
        return res.c_str();
    }
    
    inline const char* c_work() {
        static int wind = 0;
        static string res;
        ostringstream os;
        if (wind > 3) wind = 0;
        switch(wind) {
            case 0: os << "\033[1D" << "|"; break;
            case 1: os << "\033[1D" << "/"; break;
            case 2: os << "\033[1D" << "-"; break;
            case 3: os << "\033[1D" << "\\"; break;
        }
        wind++;
        res = os.str();
        return res.c_str();
    }
    
    #else
    //Textfarben
    extern char const cNORM[] = "";  //!weil intern gelinkte Objekte nicht als Template-Parameter zulaessig sind
    extern char const cBLACK[] = ""; //!soll sich im zukuenftigen C-Standard aendern => dann wieder const char* name = ""
    extern char const cRED[] = "";
    extern char const cGREEN[] = "";
    extern char const cYELLOW[] = "";
    extern char const cBLUE[] = "";
    extern char const cMAGENTA[] = "";
    extern char const cCYAN[] = "";
    extern char const cGREY[] = "";
    
    //Backgroundfarben
    extern char const bgcNORM[] = "";
    extern char const bgcBLACK[] = "";
    extern char const bgcRED[] = "";
    extern char const bgcGREEN[] = "";
    extern char const bgcYELLOW[] = "";
    extern char const bgcBLUE[] = "";
    extern char const bgcMAGENTA[] = "";
    extern char const bgcCYAN[] = "";
    extern char const bgcGREY[] = "";
    
    //einige Standardtextfragmente:
    extern char const cERROR[] = "error: ";
    extern char const cWARNING[] = "warning: ";
    
    //oder als Funktionen:
    template<const char* T>
    inline const char* c_message(const char* msg) {
        ostringstream os;
        static string res;
        os << T << msg;
        res = os.str();
        return res.c_str();
    }
    //Spezialisierung fuer errors
    template<>
    inline const char* c_message<cERROR>(const char* msg) {
        ostringstream os;
        static string res;
        os << cERROR << msg;
        res = os.str();
        return res.c_str();
    }
    //Spezialisierung fuer warnings
    template<>
    inline const char* c_message<cWARNING>(const char* msg) {
        ostringstream os;
        static string res;
        os << cWARNING << msg;
        res = os.str();
        return res.c_str();
    }
    //Fuer strings:
    template<const char* T>
    inline const char* c_message(string &msg) {
        ostringstream os;
        static string res;
        os << T << msg;
        res = os.str();
        return res.c_str();
    }
    //Spezialisierung fuer errors
    template<>
    inline const char* c_message<cERROR>(string &msg) {
        ostringstream os;
        static string res;
        os << cERROR << msg;
        res = os.str();
        return res.c_str();
    }
    //Spezialisierung fuer warnings
    template<>
    inline const char* c_message<cWARNING>(string &msg) {
        ostringstream os;
        static string res;
        os << cWARNING << msg;
        res = os.str();
        return res.c_str();
    }
    #endif
    
    inline const char* ascii_eule() {
        static string res;
        res = "\n,___,\n{0,o}\n/)__)\n-\"-\"-";
        return res.c_str();
    }
    
    inline const char* ascii_bunny() {
        static string res;
        res = "\n(\\ /)\n( . .)\nc('')('')";
        return res.c_str();
    }
    
    inline const char* ascii_flower() {
        static string res;
        res = "\nvVVVv\n(___)\n\\~Y~/\n\\\\|//";
        return res.c_str();
    }
    
    inline const char* ascii_bunny_flower() {
        static string res;
        res = "\n          vVVVv\n(\\ /)     (___)\n( . .)    \\~Y~/\nc('')('') \\\\|//\n^^^^^^^^^^^^^^^^";
        return res.c_str();
    }
    
    inline const char* ascii_bunny_flower_bee() {
        static string res;
        res = "\n           vVVVv      __         .' '.\n(\\ /)      (___)    _/__)        .   .       .\n( . .)     \\~Y~/   (8|)_}}- .      .        .\nc('')('')  \\\\|//    `\\__)    '. . ' ' .  . '\n^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^";
        return res.c_str();
    }
    
    inline const char* ascii_bee() {
        static string res;
        res = "\n    __         .' '.\n  _/__)        .   .       .\n (8|)_}}- .      .        .\n  `\\__)    '. . ' ' .  . '\n";
        return res.c_str();
    }
    
    inline const char* ascii_garfield() {
        static string res;
        res = "\n         __ __\n        ,;::\\::\\\n      ,'/' `/'`/\n  _\\,: '.,-'.-':.\n -./\"'  :    :  :\\/,\n  ::.  ,:____;__; :-\n  :\"  ( .`-*'o*',);\n   \\.. ` `---'`' /\n    `:._..-   _.'\n    ,;  .     `.\n   /\"'| |       \\\n  ::. ) :        :\n  |\" (   \\       |\n  :.(_,  :       ;\n   \\'`-'_/      /\n    `...   , _,'\n     |,|  : |\n     |`|  | |\n     |,|  | |\n ,--.;`|  | '..--.\n/;' \"' ;  '..--. ))\n\\:.___(___   ) ))'\n             `-'-";
        return res.c_str();
    }
    
    inline const char* ascii_bender() {
        static string res;
        res = "\n      _\n     ( )\n      H\n      H\n     _H_\n  .-'-.-'-.\n /         \\\n|           |\n|   .-------'._\n|  / /  '.' '. \\\n|  \\ \\ @   @ / /\n|   '---------'\n|    _______|\n|  .'-+-+-+|\n|  '.-+-+-+|\n|    \"\"\"\"\"\" |\n'-.__   __.-'\n     \"\"\"";
        return res.c_str();
    }
}
}

#endif

