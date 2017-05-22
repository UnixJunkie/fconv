
//============================================================================
// optalign_GN.hpp -*- C++ -*-; generic coordinate superposition
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
// This library is a generic implementation for optimal superposition of a
// set of coordinates, following:
// Berthold K.P. Horn, "Closed-form solution of absolute orientation using
//                      unit quaternions", Journal of the Optical Society
//                      of America A, Vol. 4, page 629, April 1987
//============================================================================


#ifndef __OPT_ALIGN
#define __OPT_ALIGN

#include<iostream>
#include<cmath>
#include<complex>
#include<vector>
#include<algorithm>
#include"linalg_GN.hpp"
#include"stl_ptr_GN.hpp"
#include"message_GN.hpp"

using namespace std;
using namespace TXT;

const float max_imag_part4eigenvalue = 0.0005;
const float max_imag_part4eigenvalue2 = 0.1;

template<typename T>
inline bool coplanar_points(vector<vec3d<T> > &L,float threshold) {
    if (L.size() < 4) return true;
    for (unsigned int i=0; i<L.size()-3; ++i) {
        if (!(is_coplanar(L[i],L[i+1],L[i+2],L[i+3],threshold))) {
            return false;
        }
    }
    return true;
}

inline double get_max_eigenvalue(double c3,double c2,double c1,double c0) {
    //! Liefert die maximale Loesung zur quartischen Gleichung:
    //! x^4 + c3x^3 + c2x^2 + c1x + c0 = 0
    //! Loesung nach Ferrari:  http://en.wikipedia.org/wiki/Quartic_equation
    
    double c32 = c3 * c3;
    
    double a = (-3. * c32 / 8.) + c2;
    double b = (c32 * c3 / 8.) - (c2 * c3 / 2.) + c1;
    double g = (-3. * c32 * c32 / 256.) + (c2 * c32 / 16.) - (c3 * c1 / 4.) + c0;
    
    double a2 = a * a;
    
    if (b == 0.) {
        complex<double> sq1 = sqrt(complex<double>(a2 - 4. * g));
        complex<double> x1 = -c3 / 4.;
        complex<double> x2 = x1; complex<double> x3 = x1; complex<double> x4 = x1;
        x1 += sqrt((-a + sq1) / 2.);
        x2 += sqrt((-a - sq1) / 2.);
        x3 -= sqrt((-a + sq1) / 2.);
        x4 -= sqrt((-a - sq1) / 2.);
        
        vector<double> real_sol;
        if (abs(x1.imag()) < max_imag_part4eigenvalue) real_sol.push_back(x1.real());
        if (abs(x2.imag()) < max_imag_part4eigenvalue) real_sol.push_back(x2.real());
        if (abs(x3.imag()) < max_imag_part4eigenvalue) real_sol.push_back(x3.real());
        if (abs(x4.imag()) < max_imag_part4eigenvalue) real_sol.push_back(x4.real());
        
        if (real_sol.size() == 0) {
            if (abs(x1.imag()) < max_imag_part4eigenvalue2) real_sol.push_back(x1.real());
            if (abs(x2.imag()) < max_imag_part4eigenvalue2) real_sol.push_back(x2.real());
            if (abs(x3.imag()) < max_imag_part4eigenvalue2) real_sol.push_back(x3.real());
            if (abs(x4.imag()) < max_imag_part4eigenvalue2) real_sol.push_back(x4.real());
            if (real_sol.size() == 0) return -1.;
            else return *max_element(real_sol.begin(),real_sol.end());
        } else return *max_element(real_sol.begin(),real_sol.end());
    } else {
        double P = (-a2 / 12.) - g;
        double Q = (-a2 * a / 108.) + (a * g / 3.) - (b * b / 8.);
        
        complex<double> sq1 = sqrt(complex<double>((Q * Q / 4.) + (P * P * P / 27.)));
        
        complex<double> R1 = Q / 2.;
        R1 += sq1;
        
        complex<double> U1 = pow(R1,1. / 3.);
        
        complex<double> y = (-5. * a / 6.) - U1;
        if (U1.real() != 0.) y += P / (3. * U1);
        
        complex<double> W = sqrt(complex<double>(a + (2. * y)));
        
        complex<double> sq21 = sqrt(complex<double>(-3. * a - 2. * y - 2. * b / W));
        complex<double> sq22 = sqrt(complex<double>(-3. * a - 2. * y + 2. * b / W));
        
        complex<double> x1 = W; x1 += sq21; x1 /= 2.; x1 -= c3 / 4.;
        complex<double> x2 = W; x2 -= sq21; x2 /= 2.; x2 -= c3 / 4.;
        complex<double> x3 = -W; x3 += sq22; x3 /= 2.; x3 -= c3 / 4.;
        complex<double> x4 = -W; x4 -= sq22; x4 /= 2.; x4 -= c3 / 4.;
        
        vector<double> real_sol;
        
        if (abs(x1.imag()) < max_imag_part4eigenvalue) real_sol.push_back(x1.real());
        if (abs(x2.imag()) < max_imag_part4eigenvalue) real_sol.push_back(x2.real());
        if (abs(x3.imag()) < max_imag_part4eigenvalue) real_sol.push_back(x3.real());
        if (abs(x4.imag()) < max_imag_part4eigenvalue) real_sol.push_back(x4.real());
        
        if (real_sol.size() == 0) {
            if (abs(x1.imag()) < max_imag_part4eigenvalue2) real_sol.push_back(x1.real());
            if (abs(x2.imag()) < max_imag_part4eigenvalue2) real_sol.push_back(x2.real());
            if (abs(x3.imag()) < max_imag_part4eigenvalue2) real_sol.push_back(x3.real());
            if (abs(x4.imag()) < max_imag_part4eigenvalue2) real_sol.push_back(x4.real());
            if (real_sol.size() == 0) return -1.;
            else return *max_element(real_sol.begin(),real_sol.end());
        } else return *max_element(real_sol.begin(),real_sol.end());
    }
}

inline double get_max_eigenvalue_for_coplanar(double c2,double c0) {
    //! Fuer coplanare Punkte vereinfacht sich die Bestimmung des maximalen Eigenwertes
    //! Davon abgesehen liefert die Loesung nach Ferrari fuer coplanare Systeme auch meist
    //! ein fehlerhaftes Ergebnis! (frueher habe ich fuer diesen Fall extra Punkte
    //! eingefuegt, die nicht in der Ebene lagen)
    complex<double> wc2 = c2 * c2; wc2 /= 4.; wc2 -= c0; wc2 = ((-0.5) * c2) + sqrt(wc2);
    wc2 = sqrt(wc2);
    return wc2.real();
}

template<typename T>
inline void quaternion_mul(T *result,T *q1,T *q2) {
    //!Multiplikation fuer Quaternionen:
    result[0] = q1[0] * q2[0] - q1[1] * q2[1] - q1[2] * q2[2] - q1[3] * q2[3];
    result[1] = q1[1] * q2[0] + q1[0] * q2[1] + q1[2] * q2[3] - q1[3] * q2[2];
    result[2] = q1[2] * q2[0] + q1[0] * q2[2] + q1[3] * q2[1] - q1[1] * q2[3];
    result[3] = q1[3] * q2[0] + q1[0] * q2[3] + q1[1] * q2[2] - q1[2] * q2[1];
}

template<typename T>
inline void get_opt_align(vector<stl_ptr<vec3d<T> > > &L, vector<stl_ptr<vec3d<T> > > &R) {
    //! L und R muessen paarweise zugeordnete Koordinatensets sein
    //! Die Funktion liefert das optimale Alignment von R auf L
    //! Hierfuer werden zunaechst die Centroide der beiden Systeme ueberlagert.
    //! Dann wird das unit Quaternion bestimmt, welches   q^T * N * q   maximiert.
    //! Dieses Quaternion ist gleichzeitig der Eigenvektor zum maximalen Eigenwert ev von N
    //! det(N-evI) = 0   kann erweitert werden zu  ev^4+c3*ev^3+c2*ev^2+c1*ev+c0
    //! Diese quartische Gleichung liefert alle Eigenwerte
    //! Der Eigenvektor q (das gesuchte Quaternion) kann dann ueber   [N - ev_max*I]*q = 0   gefunden werden
    //! Alle Zeilenvektoren der zu  [N-ev_max*I] gehoerenden Kofaktormatrix sind parallel zu q
    //! Der Betragsmaessig groesste wird ausgewaehlt (um Rechengenauigkeit zu erhoehen) und normiert.
    //! nach:  
    //!  Berthold K.P. Horn, Closed-form solution of absolute orientation using unit quaternions, Journal of the
    //!  Optical Society of America A, Vol. 4, page 629, April 1987
    //! 
    //! Im Prinzip genau das Gleiche kann man auch nachlesen in:
    //!  Evangelos A. Coutsias et al., Using quaternions to Calculate RMSD, J Comput Chem 25: 1849-1857, 2004

    //! 1.) Centroide von L und R berechnen:
    vec3d<T> cL(0.,0.,0.); vec3d<T> cR(0.,0.,0.);
    typedef typename vector<stl_ptr<vec3d<T> > >::iterator piter;
    for (piter it=L.begin(); it!=L.end(); ++it) cL += **it;
    for (piter it=R.begin(); it!=R.end(); ++it) cR += **it;
    cL /= L.size();
    cR /= R.size();
    
    //! 2.) Centroide abziehen:
    for (piter it=L.begin(); it!=L.end(); ++it) **it -= cL;
    for (piter it=R.begin(); it!=R.end(); ++it) **it -= cR;
    
    //! 3.) Hilfsgroessen zur Berechnung der Matrix N berechnen:
    double S[3][3] = {{0.,0.,0.},
                      {0.,0.,0.},
                      {0.,0.,0.}};
    for (unsigned int i=0; i<L.size(); ++i) {
        for (unsigned int ix=0; ix<3; ++ix) {
            for (unsigned int iy=0; iy<3; ++iy) {
                S[ix][iy] += (*(L[i]))[ix] * (*(R[i]))[iy];
            }
        }
    }
    // => S = Summe_ueber_alle_coords_i(L_i * R_i_transponiert)
    
    //! 4.) N berechnen (Matrix N siehe 6.)
    double a = S[0][0] + S[1][1] + S[2][2];
    double b = S[0][0] - S[1][1] - S[2][2];
    double c = -S[0][0] + S[1][1] - S[2][2];
    double d = -S[0][0] - S[1][1] + S[2][2];
    double e = S[1][2] - S[2][1];
    double f = S[0][1] + S[1][0];
    double g = S[1][2] + S[2][1];
    double h = S[2][0] - S[0][2];
    double i = S[2][0] + S[0][2];
    double j = S[0][1] - S[1][0];
    
    //! 5.) maximalen Eigenwert berechnen:
    double c3 = a + b + c + d;
    
    double c2 = -2. * (S[0][0] * S[0][0] + S[0][1] * S[0][1] + S[0][2] * S[0][2] + 
                       S[1][0] * S[1][0] + S[1][1] * S[1][1] + S[1][2] * S[1][2] + 
                       S[2][0] * S[2][0] + S[2][1] * S[2][1] + S[2][2] * S[2][2]);
    double c1 = 8. * (S[0][0] * S[1][2] * S[2][1] + S[1][1] * S[2][0] * S[0][2] + S[2][2] * S[0][1] * S[1][0])
              - 8. * (S[0][0] * S[1][1] * S[2][2] + S[1][2] * S[2][0] * S[0][1] + S[2][1] * S[1][0] * S[0][2]);
    double c0 =  (a*b - e*e) * (c*d - g*g) + (e*h - a*f) * (f*d - g*i)
               + (a*i - e*j) * (f*g - c*i) + (e*f - b*h) * (h*d - g*j)
               + (b*j - e*i) * (h*g - c*j) + (h*i - f*j) * (h*i - f*j);
    
    double ev_max;
    if (coplanar_points(L,0.3) || coplanar_points(R,0.3)) ev_max = get_max_eigenvalue_for_coplanar(c2,c0);
    else ev_max = get_max_eigenvalue(c3,c2,c1,c0);
    
    if (ev_max < 0.) {
        if (coplanar_points(L,0.9) || coplanar_points(R,0.9)) {
            cerr << c_message<cWARNING>("optalign_GN::get_opt_align() -> no positive eigenvalue -> trying solution for planar systems")
                 << endl;
            ev_max = get_max_eigenvalue_for_coplanar(c2,c0);
        } else cerr << c_message<cERROR>("optalign_GN::get_opt_align() -> could not determine a positive eigenvalue") << endl;
    }
    
    //! 6.) korrespondierenden Eigenvektor berechnen (dieser entspricht dem Quaternion fuer die Rotation):
    matrix44<double> N(a-ev_max,e,h,j,
                       e,b-ev_max,f,i,
                       h,f,c-ev_max,g,
                       j,i,g,d-ev_max);
    
    double rq[4] = {0.,0.,0.,0.};
    double magni;
    double hmag = 0.;
    for (int ir=0; ir<4; ++ir) { //aus der Kofaktormatrix die Zeile mit der groessten Magnitude nehmen
        rq[0] += pow(-1.,1+ir) * N.sub_det(ir,0);
        rq[1] += pow(-1.,2+ir) * N.sub_det(ir,1);
        rq[2] += pow(-1.,3+ir) * N.sub_det(ir,2);
        rq[3] += pow(-1.,4+ir) * N.sub_det(ir,3);
    }
    hmag = rq[0] * rq[0] + rq[1] * rq[1] + rq[2] * rq[2] + rq[3] * rq[3];
    
    magni = sqrt(hmag);
    
    rq[0] /= magni; rq[1] /= -magni; rq[2] /= -magni; rq[3] /= -magni; //!warum hier Vorzeichenwechsel weiss ich nicht
                                                                       //!sonst aber falsche Drehrichtung
    double rq_con[4] = {rq[0],
                        -rq[1],
                        -rq[2],
                        -rq[3]};
    
    //! 7.) Koordinaten transformieren:
    for (piter it=R.begin(); it!=R.end(); ++it) {
        double result[4];
        double r[4] = {0.,
                       (*it)[0],
                       (*it)[1],
                       (*it)[2]};
        quaternion_mul(result,rq,r);
        quaternion_mul(r,result,rq_con);
        (*it)[0] = r[1];
        (*it)[1] = r[2];
        (*it)[2] = r[3];
    }
    
    //! 8.) L wieder zurueckschieben und R entsprechend verschieben:
    for (piter it=L.begin(); it!=L.end(); ++it) **it += cL;
    for (piter it=R.begin(); it!=R.end(); ++it) **it += cL;
}

template<typename T>
inline void get_align_matrix(vector<vec3d<T> > &L, vector<vec3d<T> > &R,matrix<T> &rotm,vec3d<T> &cL,vec3d<T> &cR) {
    //!Wie get_opt_align, aber es werden die uebergebenen Koordinaten nicht veraendert und stattdessen die Rotationsmatrix
    //!sowie die beiden Translationsvektoren (die Centroide) in die uebergebenen Argumente geschrieben:
    cL = vec3d<T>(0.,0.,0.); cR = vec3d<T>(0.,0.,0.);
    typedef typename vector<vec3d<T> >::iterator piter;
    for (piter it=L.begin(); it!=L.end(); ++it) cL += *it;
    for (piter it=R.begin(); it!=R.end(); ++it) cR += *it;
    cL /= L.size();
    cR /= R.size();
    
    for (piter it=L.begin(); it!=L.end(); ++it) *it -= cL;
    for (piter it=R.begin(); it!=R.end(); ++it) *it -= cR;
    
    double S[3][3] = {{0.,0.,0.},
                      {0.,0.,0.},
                      {0.,0.,0.}};
    for (unsigned int i=0; i<L.size(); ++i) {
        for (unsigned int ix=0; ix<3; ++ix) {
            for (unsigned int iy=0; iy<3; ++iy) {
                S[ix][iy] += L[i][ix] * R[i][iy];
            }
        }
    }
    
    double a = S[0][0] + S[1][1] + S[2][2];
    double b = S[0][0] - S[1][1] - S[2][2];
    double c = -S[0][0] + S[1][1] - S[2][2];
    double d = -S[0][0] - S[1][1] + S[2][2];
    double e = S[1][2] - S[2][1];
    double f = S[0][1] + S[1][0];
    double g = S[1][2] + S[2][1];
    double h = S[2][0] - S[0][2];
    double i = S[2][0] + S[0][2];
    double j = S[0][1] - S[1][0];
    
    //double c3 = a + b + c + d;
    double c3 = 0.; //!sollte immer Null sein
    
    double c2 = -2. * (S[0][0] * S[0][0] + S[0][1] * S[0][1] + S[0][2] * S[0][2] + 
                       S[1][0] * S[1][0] + S[1][1] * S[1][1] + S[1][2] * S[1][2] + 
                       S[2][0] * S[2][0] + S[2][1] * S[2][1] + S[2][2] * S[2][2]);
    double c1 = 8. * (S[0][0] * S[1][2] * S[2][1] + S[1][1] * S[2][0] * S[0][2] + S[2][2] * S[0][1] * S[1][0])
              - 8. * (S[0][0] * S[1][1] * S[2][2] + S[1][2] * S[2][0] * S[0][1] + S[2][1] * S[1][0] * S[0][2]);
    double c0 =  (a*b - e*e) * (c*d - g*g) + (e*h - a*f) * (f*d - g*i)
               + (a*i - e*j) * (f*g - c*i) + (e*f - b*h) * (h*d - g*j)
               + (b*j - e*i) * (h*g - c*j) + (h*i - f*j) * (h*i - f*j);
    
    double ev_max;
    
    if (coplanar_points(L,0.3) || coplanar_points(R,0.3)) ev_max = get_max_eigenvalue_for_coplanar(c2,c0);
    else ev_max = get_max_eigenvalue(c3,c2,c1,c0);
    
    if (ev_max < 0.) {
        if (coplanar_points(L,0.9) || coplanar_points(R,0.9)) {
            cerr << c_message<cWARNING>("optalign_GN::get_align_matrix() -> no positive eigenvalue -> trying solution for planar systems")
                 << endl;
            ev_max = get_max_eigenvalue_for_coplanar(c2,c0);
        } else cerr << c_message<cERROR>("optalign_GN::get_align_matrix() -> could not determine a positive eigenvalue") << endl;
    }
    
    matrix44<double> N(a-ev_max,e,h,j,
                       e,b-ev_max,f,i,
                       h,f,c-ev_max,g,
                       j,i,g,d-ev_max);
    
    double rq[4] = {0.,0.,0.,0.};
    double magni;
    double hmag = 0.;
    for (int ir=0; ir<4; ++ir) { //aus der Kofaktormatrix die Zeile mit der groessten Magnitude nehmen
        rq[0] += pow(-1.,1+ir) * N.sub_det(ir,0);
        rq[1] += pow(-1.,2+ir) * N.sub_det(ir,1);
        rq[2] += pow(-1.,3+ir) * N.sub_det(ir,2);
        rq[3] += pow(-1.,4+ir) * N.sub_det(ir,3);
    }
    hmag = rq[0] * rq[0] + rq[1] * rq[1] + rq[2] * rq[2] + rq[3] * rq[3];
    
    magni = sqrt(hmag);
    rq[0] /= magni; rq[1] /= -magni; rq[2] /= -magni; rq[3] /= -magni;
    
    //! Rotationsmatrix berechnen:
    T q02 = rq[0] * rq[0];
    T qx2 = rq[1] * rq[1];
    T qy2 = rq[2] * rq[2];
    T qz2 = rq[3] * rq[3];
    rotm[0][0] = q02 + qx2 - qy2 - qz2;
    rotm[0][1] = 2. * (rq[1] * rq[2] - rq[0] * rq[3]);
    rotm[0][2] = 2. * (rq[1] * rq[3] + rq[0] * rq[2]);
    
    rotm[1][0] = 2. * (rq[2] * rq[1] + rq[0] * rq[3]);
    rotm[1][1] = q02 - qx2 + qy2 - qz2;
    rotm[1][2] = 2. * (rq[2] * rq[3] - rq[0] * rq[1]);
    
    rotm[2][0] = 2. * (rq[3] * rq[1] - rq[0] * rq[2]);
    //rotm[2][1] = 2. * (rq[3] * rq[2] + rq[0] * rq[3]); //Fehler im Horn-Paper S.641
    rotm[2][1] = 2. * (rq[3] * rq[2] + rq[0] * rq[1]);
    rotm[2][2] = q02 - qx2 - qy2 + qz2;
    
    for (piter it=L.begin(); it!=L.end(); ++it) *it += cL;
    for (piter it=R.begin(); it!=R.end(); ++it) *it += cR;
}

template<typename T>
inline void get_weighted_align_matrix(vector<vec3d<T> > &L,vector<vec3d<T> > &R,vector<T> &weights,matrix<T> &rotm,vec3d<T> &cL,vec3d<T> &cR) {
    //! Hier die Variante mit gewichteten Punkten:
    
    cL = vec3d<T>(0.,0.,0.); cR = vec3d<T>(0.,0.,0.);
    
    T wsum = 0.; vec3d<T> tvec;
    for (unsigned int i=0; i<L.size(); ++i) {
        tvec = L[i]; tvec *= weights[i]; cL += tvec;
        tvec = R[i]; tvec *= weights[i]; cR += tvec;
        wsum += weights[i];
    }
    
    cL /= wsum;
    cR /= wsum;
    
    for (unsigned int i=0; i<L.size(); ++i) {
        L[i] -= cL;
        R[i] -= cR;
    }
    
    double S[3][3] = {{0.,0.,0.},
                      {0.,0.,0.},
                      {0.,0.,0.}};
    for (unsigned int i=0; i<L.size(); ++i) {
        for (unsigned int ix=0; ix<3; ++ix) {
            for (unsigned int iy=0; iy<3; ++iy) {
                S[ix][iy] += weights[i] * L[i][ix] * R[i][iy];
            }
        }
    }
    
    double a = S[0][0] + S[1][1] + S[2][2];
    double b = S[0][0] - S[1][1] - S[2][2];
    double c = -S[0][0] + S[1][1] - S[2][2];
    double d = -S[0][0] - S[1][1] + S[2][2];
    double e = S[1][2] - S[2][1];
    double f = S[0][1] + S[1][0];
    double g = S[1][2] + S[2][1];
    double h = S[2][0] - S[0][2];
    double i = S[2][0] + S[0][2];
    double j = S[0][1] - S[1][0];
    
    double c3 = 0.; //!sollte immer Null sein
    
    double c2 = -2. * (S[0][0] * S[0][0] + S[0][1] * S[0][1] + S[0][2] * S[0][2] + 
                       S[1][0] * S[1][0] + S[1][1] * S[1][1] + S[1][2] * S[1][2] + 
                       S[2][0] * S[2][0] + S[2][1] * S[2][1] + S[2][2] * S[2][2]);
    double c1 = 8. * (S[0][0] * S[1][2] * S[2][1] + S[1][1] * S[2][0] * S[0][2] + S[2][2] * S[0][1] * S[1][0])
              - 8. * (S[0][0] * S[1][1] * S[2][2] + S[1][2] * S[2][0] * S[0][1] + S[2][1] * S[1][0] * S[0][2]);
    double c0 =  (a*b - e*e) * (c*d - g*g) + (e*h - a*f) * (f*d - g*i)
               + (a*i - e*j) * (f*g - c*i) + (e*f - b*h) * (h*d - g*j)
               + (b*j - e*i) * (h*g - c*j) + (h*i - f*j) * (h*i - f*j);
    
    double ev_max;
    if (coplanar_points(L,0.3) || coplanar_points(R,0.3)) ev_max = get_max_eigenvalue_for_coplanar(c2,c0);
    else ev_max = get_max_eigenvalue(c3,c2,c1,c0);
    
    if (ev_max < 0.) {
        if (coplanar_points(L,0.9) || coplanar_points(R,0.9)) {
            cerr << c_message<cWARNING>("optalign_GN::get_align_matrix() -> no positive eigenvalue -> trying solution for planar systems")
                 << endl;
            ev_max = get_max_eigenvalue_for_coplanar(c2,c0);
        } else cerr << c_message<cERROR>("optalign_GN::get_align_matrix() -> could not determine a positive eigenvalue") << endl;
    }
    
    matrix44<double> N(a-ev_max,e,h,j,
                       e,b-ev_max,f,i,
                       h,f,c-ev_max,g,
                       j,i,g,d-ev_max);
    
    double rq[4] = {0.,0.,0.,0.};
    double magni;
    double hmag = 0.;
    for (int ir=0; ir<4; ++ir) { //aus der Kofaktormatrix die Zeile mit der groessten Magnitude nehmen
        rq[0] += pow(-1.,1+ir) * N.sub_det(ir,0);
        rq[1] += pow(-1.,2+ir) * N.sub_det(ir,1);
        rq[2] += pow(-1.,3+ir) * N.sub_det(ir,2);
        rq[3] += pow(-1.,4+ir) * N.sub_det(ir,3);
    }
    hmag = rq[0] * rq[0] + rq[1] * rq[1] + rq[2] * rq[2] + rq[3] * rq[3];
    
    magni = sqrt(hmag);
    rq[0] /= magni; rq[1] /= -magni; rq[2] /= -magni; rq[3] /= -magni;
    
    //! Rotationsmatrix berechnen:
    T q02 = rq[0] * rq[0];
    T qx2 = rq[1] * rq[1];
    T qy2 = rq[2] * rq[2];
    T qz2 = rq[3] * rq[3];
    rotm[0][0] = q02 + qx2 - qy2 - qz2;
    rotm[0][1] = 2. * (rq[1] * rq[2] - rq[0] * rq[3]);
    rotm[0][2] = 2. * (rq[1] * rq[3] + rq[0] * rq[2]);
    
    rotm[1][0] = 2. * (rq[2] * rq[1] + rq[0] * rq[3]);
    rotm[1][1] = q02 - qx2 + qy2 - qz2;
    rotm[1][2] = 2. * (rq[2] * rq[3] - rq[0] * rq[1]);
    
    rotm[2][0] = 2. * (rq[3] * rq[1] - rq[0] * rq[2]);
    //rotm[2][1] = 2. * (rq[3] * rq[2] + rq[0] * rq[3]); //Fehler im Horn-Paper S.641
    rotm[2][1] = 2. * (rq[3] * rq[2] + rq[0] * rq[1]);
    rotm[2][2] = q02 - qx2 - qy2 + qz2;
    
    
    for (unsigned int i=0; i<L.size(); ++i) {
        L[i] += cL;
        R[i] += cR;
    }
}

template<typename T>
inline void get_align_quaternion(vector<vec3d<T> > &L, vector<vec3d<T> > &R,T *quaternion,vec3d<T> &cL,vec3d<T> &cR) {
    //!Wie get_opt_align, aber es werden die uebergebenen Koordinaten nicht veraendert und stattdessen die Rotationsmatrix
    //!Sowie die beiden Translationsvektoren (die Centroide) in die uebergebenen Argumente geschrieben:
    
    cL = vec3d<T>(0.,0.,0.); cR = vec3d<T>(0.,0.,0.);
    typedef typename vector<vec3d<T> >::iterator piter;
    for (piter it=L.begin(); it!=L.end(); ++it) cL += *it;
    for (piter it=R.begin(); it!=R.end(); ++it) cR += *it;
    cL /= L.size();
    cR /= R.size();
    
    for (piter it=L.begin(); it!=L.end(); ++it) *it -= cL;
    for (piter it=R.begin(); it!=R.end(); ++it) *it -= cR;
    
    double S[3][3] = {{0.,0.,0.},
                      {0.,0.,0.},
                      {0.,0.,0.}};
    for (unsigned int i=0; i<L.size(); ++i) {
        for (unsigned int ix=0; ix<3; ++ix) {
            for (unsigned int iy=0; iy<3; ++iy) {
                S[ix][iy] += L[i][ix] * R[i][iy];
            }
        }
    }
    
    double a = S[0][0] + S[1][1] + S[2][2];
    double b = S[0][0] - S[1][1] - S[2][2];
    double c = -S[0][0] + S[1][1] - S[2][2];
    double d = -S[0][0] - S[1][1] + S[2][2];
    double e = S[1][2] - S[2][1];
    double f = S[0][1] + S[1][0];
    double g = S[1][2] + S[2][1];
    double h = S[2][0] - S[0][2];
    double i = S[2][0] + S[0][2];
    double j = S[0][1] - S[1][0];
    
    double c3 = a + b + c + d;
    
    double c2 = -2. * (S[0][0] * S[0][0] + S[0][1] * S[0][1] + S[0][2] * S[0][2] + 
                       S[1][0] * S[1][0] + S[1][1] * S[1][1] + S[1][2] * S[1][2] + 
                       S[2][0] * S[2][0] + S[2][1] * S[2][1] + S[2][2] * S[2][2]);
    double c1 = 8. * (S[0][0] * S[1][2] * S[2][1] + S[1][1] * S[2][0] * S[0][2] + S[2][2] * S[0][1] * S[1][0])
              - 8. * (S[0][0] * S[1][1] * S[2][2] + S[1][2] * S[2][0] * S[0][1] + S[2][1] * S[1][0] * S[0][2]);
    double c0 =  (a*b - e*e) * (c*d - g*g) + (e*h - a*f) * (f*d - g*i)
               + (a*i - e*j) * (f*g - c*i) + (e*f - b*h) * (h*d - g*j)
               + (b*j - e*i) * (h*g - c*j) + (h*i - f*j) * (h*i - f*j);
    
    double ev_max;
    if (coplanar_points(L,0.3) || coplanar_points(R,0.3)) ev_max = get_max_eigenvalue_for_coplanar(c2,c0);
    else ev_max = get_max_eigenvalue(c3,c2,c1,c0);
    
    if (ev_max < 0.) {
        if (coplanar_points(L,0.9) || coplanar_points(R,0.9)) {
            cerr << c_message<cWARNING>("optalign_GN::get_align_quaternion() -> no positive eigenvalue -> trying solution for planar systems")
                 << endl;
            ev_max = get_max_eigenvalue_for_coplanar(c2,c0);
        } else cerr << c_message<cERROR>("optalign_GN::get_align_quaternion() -> could not determine a positive eigenvalue") << endl;
    }
    
    matrix44<double> N(a-ev_max,e,h,j,
                       e,b-ev_max,f,i,
                       h,f,c-ev_max,g,
                       j,i,g,d-ev_max);
    
    double rq[4] = {0.,0.,0.,0.};
    double magni;
    double hmag = 0.;
    for (int ir=0; ir<4; ++ir) { //aus der Kofaktormatrix die Zeile mit der groessten Magnitude nehmen
        rq[0] += pow(-1.,1+ir) * N.sub_det(ir,0);
        rq[1] += pow(-1.,2+ir) * N.sub_det(ir,1);
        rq[2] += pow(-1.,3+ir) * N.sub_det(ir,2);
        rq[3] += pow(-1.,4+ir) * N.sub_det(ir,3);
    }
    hmag = rq[0] * rq[0] + rq[1] * rq[1] + rq[2] * rq[2] + rq[3] * rq[3];
    
    magni = sqrt(hmag);
    
    quaternion[0] = rq[0] / magni; quaternion[1] = rq[1] / -magni;
    quaternion[2] = rq[2] / -magni; quaternion[3] = rq[3] / -magni;
    
    for (piter it=L.begin(); it!=L.end(); ++it) *it += cL;
    for (piter it=R.begin(); it!=R.end(); ++it) *it += cR;
}

#endif
