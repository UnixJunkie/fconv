
//============================================================================
// linalg_GN.hpp -*- C++ -*-; vector and matrix calculations
//
// Copyright (C) 2006, 2007, 2008, 2009, 2010 Gerd Neudert
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
// This library provides objects for vector and matrix representations and
// implements the most important calculations.
// I invested much time for many of those functions and for most of them
// you should think twice before changing them (of course your always welcome
// to discuss things with me)
//============================================================================


#ifndef __LINALGGN
#define __LINALGGN
#include<iostream>
#include<cmath>
#include<vector>

using namespace std;

//==================================================================================================
//!Forward-Deklarationen:
//==================================================================================================

template<class T> class vec3d;
template<class T> class matrix;
template<class T> class matrix44;
template<class T> class QUATERNION;
template<class T> class INDEXER;
template<class T> class INDEXER44;
template<class T> const vec3d<T> operator+(vec3d<T> const& links,vec3d<T> const& rechts); //addiert 2 Vektoren (basiert auf += =>mu�nicht friend sein)
template<class T> const vec3d<T> operator-(vec3d<T> const& links,vec3d<T> const& rechts); //subtrahiert 2 Vektoren
template<class T> const vec3d<T> operator*(vec3d<T> const& links,vec3d<T> const& rechts); //Vektorprodukt
template<class T> const vec3d<T> operator*(vec3d<T> const& links,T const& rechts); //Multiplikation mit Skalar
template<class T> const vec3d<T> operator*(T const& links,vec3d<T> const& rechts); //        -"-
template<class T> const vec3d<T> operator*(matrix<T> const& links,vec3d<T> const& rechts); //Vektor mit Matrix transformieren
template<class T> const vec3d<T> operator*(vec3d<T> const& links,matrix<T> const& rechts); //Vektor mit Matrix transformieren (falschrum, aber aus bequemlichkeit auch moeglich)
template<class T> const vec3d<T> operator/(vec3d<T> const& links,T const& rechts); //Division mit Skalar
template<class T> inline const T skalarproduct(vec3d<T> const& links,vec3d<T> const& rechts); //bildet das Skalarprodukt aus 2 Vektoren
template<class T> inline const vec3d<T> vectorproduct(vec3d<T> const& links,vec3d<T> const& rechts); //bildet das Vektorprodukt aus 2 Vektoren
template<class T> inline const T angle(vec3d<T> const& links,vec3d<T> const& rechts); //liefert den Winkel zwischen 2 Vektoren
template<class T> inline const T angle_for_normed(vec3d<T> const& links,vec3d<T> const& rechts);
template<class T> inline const T angle(vec3d<T> const& p1,vec3d<T> const& p2,vec3d<T> const& p3); //Winkel ueber 3 Koordinaten
template<class T> inline const T get_distance(vec3d<T> const& v1,vec3d<T> const& v2); //Distanz zwischen 2 Punkten
template<class T> inline const T get_square_distance(vec3d<T> const& v1,vec3d<T> const& v2); //Quadrat der Distanz zwischen 2 Punkten
template<class T> const matrix<T> rotmatrix(vec3d<T> const& vec,double const& angle); //liefert die matrix fr Rotation um vec
template<class T> inline const T determinant(matrix<T> const& m); //Determinante einer Matrix bestimmen
template<class T> inline bool is_almost_equal(matrix<T> const& m1,matrix<T> const& m2);
template<class T> ostream &operator<<(ostream &os,vec3d<T> const& vec); //Ausgabeoperator fr vec3d-Obj.
template<class T> ostream &operator<<(ostream &os,matrix<T> const& m); //Ausgabeoperator fr matrix-Obj.
template<class T> ostream &operator<<(ostream &os,matrix44<T> const& m);

template<class T>
class vec3d {
    private:
        T vector[3];
    public:
        vec3d(T const& x=0,T const& y=0,T const& z=0);
        vec3d(const T vec[]);
        vec3d(vec3d<T> const& vec);
        template<class T2> vec3d(vec3d<T2> const& vec);
        ~vec3d();
        
        inline void add(vec3d<T> const& vec); //Vektor addieren
        inline void sub(vec3d<T> const& vec); //Vektor subtrahieren
        inline void vector_product(vec3d<T> const& vec); //Vektorprodukt bilden
        inline void quat_mul(QUATERNION<T> const& q);
        inline void quat_mul(T const *const q);
        inline const T skalar_product(vec3d<T> const& vec); //Skalarprodukt bilden
        template<class T2> inline void mul_skalar(T2 const& skalar); //Skalar multiplizieren (Template spart casts)
        template<class T2> inline void div_skalar(T2 const& skalar); //Division mit Skalar
        inline const T value() const; //Betrag des Vektors
        inline void norm(); //Vektor normieren
        inline void transform(matrix<T> const& mat); //Vector transformieren
        inline void transform(matrix44<T> const& mat); //Vector transformieren
    
        inline T* const get_vec() const; //liefert Zeiger auf das vector-array   //EXPERIMENTELL
        inline T* get_vec();
        
        inline const vec3d<T>& operator=(vec3d<T> const& rechts); //Zuweisungsoperator
        inline const vec3d<T>& operator=(T rechts[]); //Zuweisung von array (keine Prfung!)
        inline const vec3d<T>& operator+=(vec3d<T> const& rechts); // += Operator (der + Operator basiert darauf)
        inline const vec3d<T>& operator-=(vec3d<T> const& rechts); // -= Operator (- beruht darauf) fr Vektorprodukt
        template<class T2> const vec3d<T>& operator*=(T2 const& rechts); // *= Operator (* beruht darauf) fr Multipl. mit Skalar
        template<class T2> const vec3d<T>& operator/=(T2 const& rechts); // /= Operator fr Division durch Skalar
        inline const vec3d<T>& operator*=(vec3d<T> const& rechts); // *= fr Vektorprodukt
        inline const vec3d<T>& operator*=(QUATERNION<T> const& rechts);
        inline const vec3d<T>& operator*=(matrix<T> const& rechts); //fr Multiplikation mit Matrix
        inline const vec3d<T>& operator*=(matrix44<T> const& rechts);
        inline T& operator[](int const& index); // [] Operator fr Lesen und Schreiben der Vektorkomponenten
        inline const T& operator[](int const& index) const;
        inline operator T*(); //fuer Cast zum entsprechenden Array
        inline bool operator<(vec3d<T> const& rechts); //!hier wird nicht der Betrag verglichen, sondern ob alle Komponenten kleiner
        inline bool operator>(vec3d<T> const& rechts); //!bzw. groesser sind
        inline bool operator<=(vec3d<T> const& rechts);
        inline bool operator>=(vec3d<T> const& rechts);

        friend const T skalarproduct<>(vec3d<T> const& links,vec3d<T> const& rechts);
        friend const vec3d<T> vectorproduct<>(vec3d<T> const& links,vec3d<T> const& rechts); // * ist berladen
        friend const T angle<>(vec3d<T> const& links,vec3d<T> const& rechts); //Winkel zwischen 2 Vektoren
        friend const T angle_for_normed<>(vec3d<T> const& links,vec3d<T> const& rechts);
        friend const T angle<>(vec3d<T> const& p1,vec3d<T> const& p2,vec3d<T> const& p3);
        friend const T get_distance<>(vec3d<T> const& v1,vec3d<T> const& v2);
        friend const T get_square_distance<>(vec3d<T> const& v1,vec3d<T> const& v2);
        friend const matrix<T> rotmatrix<>(vec3d<T> const& vec,double const& angle); //Rotationsmatrix zurckgeben
        friend ostream &operator<< <>(ostream &os,vec3d<T> const& vec); //Ausgabeoperator
        template<class T2> friend class vec3d; //!hier darf kein <T2> angegeben werden!!!
};


template<class T>
class matrix {
    private:
        T mat[3][3];
    public:
        matrix() {mat[0][0] = 0.; mat[0][1] = 0.; mat[0][2] = 0.;
                  mat[1][0] = 0.; mat[1][1] = 0.; mat[1][2] = 0.;
                  mat[2][0] = 0.; mat[2][1] = 0.; mat[2][2] = 0.;}
        matrix(T const& a00,T const& a01,T const& a02,
               T const& a10,T const& a11,T const& a12,
               T const& a20,T const& a21,T const& a22);
        matrix(matrix<T> const& m);
        ~matrix();
        
        inline const T det();
        inline bool invert();
        inline void transpone();
        inline void mul_matrix(matrix<T> const& rechts);
        
        inline bool diagonalise();
        
        inline bool gaussj_invert();
        inline bool sqrt_gj();
        
        inline void get_rot_v(T a,vec3d<T> *rv);
    //    inline bool check_a_v(T a,vec3d<T> &b);
    //    inline bool full_check_a_v(T a,vec3d<T> &b);
        inline T get_rot_angle();
        inline vec3d<T> get_rot_vector();
        
        inline void get_cardan_angles(float &a,float &b,float &c);
        inline void get_fixed_angles(float &x,float &y,float &z);
        
        inline const matrix<T>& operator=(matrix<T> const& rechts);
        inline INDEXER<T> operator[](int const& index);
        inline const matrix<T>& operator*=(matrix<T> const& rechts);
        
        friend void vec3d<T>::transform(matrix<T> const& mat);
        friend const matrix<T> rotmatrix<>(vec3d<T> const& vec,double const& angle);
        friend const T determinant<>(matrix<T> const& m);
        friend bool is_almost_equal<>(matrix<T> const& m1,matrix<T> const& m2);
        friend ostream &operator<< <>(ostream &os,matrix<T> const& m);
        friend T& INDEXER<T>::operator[](int const& xindex);
};


template<class T>
class matrix44 {
    private:
        T mat[4][4];
    public:
        matrix44() {mat[0][0] = 0.; mat[0][1] = 0.; mat[0][2] = 0.; mat[0][3] = 0.;
                    mat[1][0] = 0.; mat[1][1] = 0.; mat[1][2] = 0.; mat[1][3] = 0.;
                    mat[2][0] = 0.; mat[2][1] = 0.; mat[2][2] = 0.; mat[2][3] = 0.;
                    mat[3][0] = 0.; mat[3][1] = 0.; mat[3][2] = 0.; mat[3][3] = 0.;}
        matrix44(T const& a00,T const& a01,T const& a02,T const& a03,
                 T const& a10,T const& a11,T const& a12,T const& a13,
                 T const& a20,T const& a21,T const& a22,T const& a23,
                 T const& a30,T const& a31,T const& a32,T const& a33);
        matrix44(matrix44<T> const& m);
        ~matrix44();
        
        inline const T sub_det(int row,int column);
        inline const T det();
        inline bool invert();
        
        inline void mul_matrix(matrix44<T> const& rechts);
        
        inline const matrix44<T>& operator*=(matrix44<T> const& rechts);
        inline const matrix44<T>& operator=(matrix44<T> const& rechts);
        inline INDEXER44<T> operator[](int const& index);
        
        friend void vec3d<T>::transform(matrix44<T> const& mat);
        friend ostream &operator<< <>(ostream &os,matrix44<T> const& m);
        friend T& INDEXER44<T>::operator[](int const& xindex);
};


template<class T>
class matrix55 {
    private:
        T mat[5][5];
    public:
        matrix55() {mat[0][0] = 0.; mat[0][1] = 0.; mat[0][2] = 0.; mat[0][3] = 0.; mat[0][4] = 0.;
                    mat[1][0] = 0.; mat[1][1] = 0.; mat[1][2] = 0.; mat[1][3] = 0.; mat[1][4] = 0.;
                    mat[2][0] = 0.; mat[2][1] = 0.; mat[2][2] = 0.; mat[2][3] = 0.; mat[2][4] = 0.;
                    mat[3][0] = 0.; mat[3][1] = 0.; mat[3][2] = 0.; mat[3][3] = 0.; mat[3][4] = 0.;
                    mat[4][0] = 0.; mat[4][1] = 0.; mat[4][2] = 0.; mat[4][3] = 0.; mat[4][4] = 0.;}
        matrix55(T const& a00 = 0.,T const& a01 = 0.,T const& a02 = 0.,T const& a03 = 0.,T const& a04 = 0.,
                 T const& a10 = 0.,T const& a11 = 0.,T const& a12 = 0.,T const& a13 = 0.,T const& a14 = 0.,
                 T const& a20 = 0.,T const& a21 = 0.,T const& a22 = 0.,T const& a23 = 0.,T const& a24 = 0.,
                 T const& a30 = 0.,T const& a31 = 0.,T const& a32 = 0.,T const& a33 = 0.,T const& a34 = 0.,
                 T const& a40 = 0.,T const& a41 = 0.,T const& a42 = 0.,T const& a43 = 0.,T const& a44 = 0.);
        matrix55(matrix55<T> const& m);
        ~matrix55();
        
        inline const T sub_det(int row,int column);
        inline const T det();
        
        inline INDEXER<T> operator[](int const& index);
        
        friend T& INDEXER<T>::operator[](int const& xindex);
};


template<class T>
class QUATERNION {
    public:
        T quat[4];
        
        QUATERNION(T const& w=1.,T const& x=0.,T const& y=0.,T const& z=0.) {
            quat[0] = w; quat[1] = x; quat[2] = y; quat[3] = z;
        }
        QUATERNION(QUATERNION<T> const& q) {memcpy(quat,q.quat,4*sizeof(T));}
        QUATERNION(vec3d<T> const& vec) {quat[0] = 0; quat[1] = vec[0]; quat[2] = vec[1]; quat[3] = vec[2];}
        QUATERNION(vec3d<T> const& v,T const& angle) {
            //! Initialisierung mit Drehachse und Winkel
            //! ACHTUNG: die Drehachse MUSS bereits normiert sein!!!
            quat[0] = cos(angle/2.);
            T s = sin(angle/2.);
            quat[1] = v[0] * s;
            quat[2] = v[1] * s;
            quat[3] = v[2] * s;
            norm();
        }
        ~QUATERNION() {}
        void norm() {
            T normer = sqrt(quat[0]*quat[0]+quat[1]*quat[1]+quat[2]*quat[2]+quat[3]*quat[3]);
            quat[0] /= normer;
            quat[1] /= normer;
            quat[2] /= normer;
            quat[3] /= normer;
        }
        void invert() {
            quat[1] = -quat[1];
            quat[2] = -quat[2];
            quat[3] = -quat[3];
        }
        void quat_mul(QUATERNION<T> const& b) {
            T q[4] = {quat[0]*b.quat[0] - quat[1]*b.quat[1]  -  quat[2]*b.quat[2]  -  quat[3]*b.quat[3],
                      quat[0]*b.quat[1] + quat[1]*b.quat[0]  +  quat[2]*b.quat[3]  -  quat[3]*b.quat[2],
                      quat[0]*b.quat[2] + quat[2]*b.quat[0]  +  quat[3]*b.quat[1]  -  quat[0]*b.quat[3],
                      quat[0]*b.quat[3] + quat[3]*b.quat[0]  +  quat[1]*b.quat[2]  -  quat[2]*b.quat[1]};
            quat[0] = q[0];
            quat[1] = q[1];
            quat[2] = q[2];
            quat[3] = q[3];
        }
        const QUATERNION<T>& operator*=(QUATERNION<T> const& rechts) {
            quat_mul(rechts);
            return *this;
        }
        void lerp(QUATERNION<T> const& a,QUATERNION<T> const& b,T const& w2) {
            //! von a nach b interpolieren
            //! w2 = [0.,1.]
            float w1 = 1. - w2;
            quat[0] = a.quat[0]*w1 + b.quat[0]*w2;
            quat[1] = a.quat[1]*w1 + b.quat[1]*w2;
            quat[2] = a.quat[2]*w1 + b.quat[2]*w2;
            quat[3] = a.quat[3]*w1 + b.quat[3]*w2;
            norm();
        }
        void slerp(QUATERNION<T> const& a,QUATERNION<T> const& b,T const& t) {
            //! ebenfalls von a nach b interpolieren, aber mit Fuehrung des
            //! Quaternions auf einer Kreisbahn => immer gleiche Geschwindigkeit
            //! normierung nicht noetig <- doch?
            T w1, w2;
            T cosa(quat[1]*b.quat[1] + quat[2]*b.quat[2] + quat[3]*b.quat[3] + quat[0]*b.quat[0]);
            T wa(acos(cosa));
            T sina(sin(wa));
            if (sina > 0.001) {
                w1 = sin(wa*(1.-t)) / sina;
                w2 = sin(wa*t) / sina;
            } else {
                w1 = 1. - t;
                w2 = t;
            }

            quat[0] = w1*a.quat[0] + w2*b.quat[0];
            quat[1] = w1*a.quat[1] + w2*b.quat[1];
            quat[2] = w1*a.quat[2] + w2*b.quat[2];
            quat[3] = w1*a.quat[3] + w2*b.quat[3];
        }
};


template<class T>
inline void quat_lerp(T *result,QUATERNION<T> const& a,QUATERNION<T> const& b,const T& w2) {
    QUATERNION<T> q; q.lerp(a,b,w2);
    result[0] = q.quat[0];
    result[1] = q.quat[1];
    result[2] = q.quat[2];
    result[3] = q.quat[3];
}

template<class T>
inline void quat_slerp(T *result,QUATERNION<T> const& a,QUATERNION<T> const& b,const T& w2) {
    QUATERNION<T> q; q.slerp(a,b,w2);
    result[0] = q.quat[0];
    result[1] = q.quat[1];
    result[2] = q.quat[2];
    result[3] = q.quat[3];
}

template<class T>
void vec3d<T>::quat_mul(QUATERNION<T> const& q) {
    T Q1(q.quat[1]*vector[0] + q.quat[2]*vector[1] + q.quat[3]*vector[2]);
    T Q2(q.quat[3]*vector[0] + q.quat[0]*vector[1] - q.quat[1]*vector[2]);
    T Q3(q.quat[0]*vector[0] - q.quat[3]*vector[1] + q.quat[2]*vector[2]);
    T Q4((-q.quat[2])*vector[0] + q.quat[1]*vector[1] + q.quat[0]*vector[2]);
    vector[0] = Q3*q.quat[0] + Q1*q.quat[1] + Q4*q.quat[2] - Q2*q.quat[3];
    vector[1] = Q2*q.quat[0] - Q4*q.quat[1] + Q1*q.quat[2] + Q3*q.quat[3];
    vector[2] = Q4*q.quat[0] + Q2*q.quat[1] - Q3*q.quat[2] + Q1*q.quat[3];
}

template<class T>
void vec3d<T>::quat_mul(T const *const q) {
    T Q1(q[1]*vector[0] + q[2]*vector[1] + q[3]*vector[2]);
    T Q2(q[3]*vector[0] + q[0]*vector[1] - q[1]*vector[2]);
    T Q3(q[0]*vector[0] - q[3]*vector[1] + q[2]*vector[2]);
    T Q4((-q[2])*vector[0] + q[1]*vector[1] + q[0]*vector[2]);
    vector[0] = Q3*q[0] + Q1*q[1] + Q4*q[2] - Q2*q[3];
    vector[1] = Q2*q[0] - Q4*q[1] + Q1*q[2] + Q3*q[3];
    vector[2] = Q4*q[0] + Q2*q[1] - Q3*q[2] + Q1*q[3];
}

template<class T>
inline float get_quat_dist(T *result,QUATERNION<T> const& a,QUATERNION<T> const& b) {
    // euklidische Distanz b-a:
    // ACHTUNG: Es gibt immer 2 Quaternionen, die zur gleichen Drehung fuehren (q und -q)
    //         => Es ist natuerlich der geringere Abstand gewuenscht
    //         => Die Quaternionen sollten deshalb immer auf positive w von 0 bis 1
    //            beschraenkt werden (entsprechend Drehwinkel 0 bis PI)
    //         -> Ist dies nicht der Fall kann get_quat_dist_with_check genutzt werden
    result[0] = b.quat[0] - a.quat[0];
    result[1] = b.quat[1] - a.quat[1];
    result[2] = b.quat[2] - a.quat[2];
    result[3] = b.quat[3] - a.quat[3];
    return(sqrt(result[0]*result[0]+result[1]*result[1]+result[2]*result[2]+result[3]*result[3]));
}

template<class T>
inline float get_quat_dist(T *result,T const * const a,T const * const b) {
    // euklidische Distanz b-a:
    // ACHTUNG: Es gibt immer 2 Quaternionen, die zur gleichen Drehung fuehren (q und -q)
    //         => Es ist natuerlich der geringere Abstand gewuenscht
    //         => Die Quaternionen sollten deshalb immer auf positive w von 0 bis 1
    //            beschraenkt werden (entsprechend Drehwinkel 0 bis PI)
    //         -> Ist dies nicht der Fall kann get_quat_dist_with_check genutzt werden
    result[0] = b[0] - a[0];
    result[1] = b[1] - a[1];
    result[2] = b[2] - a[2];
    result[3] = b[3] - a[3];
    return(sqrt(result[0]*result[0]+result[1]*result[1]+result[2]*result[2]+result[3]*result[3]));
}

template<class T>
inline float get_quat_dist_with_check(T *result,QUATERNION<T> const& a,QUATERNION<T> const& b) {
    // euklidische Distanz:
    // Wenn der Winkel zwischen den Vektoren > PI/2, dann eines der QUATERNIONen negieren (NICHT invertieren):
    // (Das negierte Quat liefert die gleiche Rotation, aber der Abstand ist geringer)
    if (fabs(b.quat[0]-a.quat[0]) > 1.) {
        result[0] = b.quat[0] + a.quat[0];
        result[1] = b.quat[1] + a.quat[1];
        result[2] = b.quat[2] + a.quat[2];
        result[3] = b.quat[3] + a.quat[3];
    } else {
        result[0] = b.quat[0] - a.quat[0];
        result[1] = b.quat[1] - a.quat[1];
        result[2] = b.quat[2] - a.quat[2];
        result[3] = b.quat[3] - a.quat[3];
    }
    return(sqrt(result[0]*result[0]+result[1]*result[1]+result[2]*result[2]+result[3]*result[3]));
}


//==================================================================================================
//!Deklaration der Klasse INDEXER:
//==================================================================================================

template<class T>
class INDEXER {
    private:
        int y;
        matrix<T> *m;
    public:
        INDEXER(int yindex, matrix<T> *matr);
        ~INDEXER();
        
        inline T& operator[](int const& xindex);
};

template<class T>
class INDEXER44 {
    private:
        int y;
        matrix44<T> *m;
    public:
        INDEXER44(int yindex, matrix44<T> *matr);
        ~INDEXER44();
        
        inline T& operator[](int const& xindex);
};


//==============================================================================================
//!Definitionen fr vec3d:
//==============================================================================================

template<class T> 
vec3d<T>::vec3d(T const& x,T const& y,T const& z) {
    vector[0] = x; vector[1] = y; vector[2] = z;
}

template<class T> 
vec3d<T>::vec3d(const T vec[]) {
    vector[0] = vec[0]; vector[1] = vec[1]; vector[2] = vec[2];
}

template<class T> 
vec3d<T>::vec3d(vec3d<T> const& vec) {
    vector[0] = vec.vector[0]; vector[1] = vec.vector[1]; vector[2] = vec.vector[2];
}

template<class T> template<class T2>
vec3d<T>::vec3d(vec3d<T2> const& vec) {
    vector[0] = vec.vector[0]; vector[1] = vec.vector[1]; vector[2] = vec.vector[2];
}

template<class T> 
vec3d<T>::~vec3d() {}


template<class T> 
void vec3d<T>::add(vec3d<T> const& vec) { //Addition
    vector[0] += vec.vector[0]; vector[1] += vec.vector[1]; vector[2] += vec.vector[2];
}

template<class T> 
void vec3d<T>::sub(vec3d<T> const& vec) { //Subtraktion
    vector[0] -= vec.vector[0]; vector[1] -= vec.vector[1]; vector[2] -= vec.vector[2];
}

template<class T> 
const T vec3d<T>::value() const { //Betrag des Vektors
    return sqrt(vector[0]*vector[0]+vector[1]*vector[1]+vector[2]*vector[2]);
}

template<class T> 
void vec3d<T>::norm() { //Vektor normieren
    T buf = value();
    vector[0] /= buf; vector[1] /= buf; vector[2] /= buf;
}

template<class T> template<class T2> 
void vec3d<T>::mul_skalar(T2 const& skalar) { //Vektor mit Skalar multiplizieren
    vector[0] *= skalar; vector[1] *= skalar; vector[2] *= skalar;
}

template<class T> template<class T2> 
void vec3d<T>::div_skalar(T2 const& skalar) { //Vektor durch Skalar teilen
    vector[0] /= skalar; vector[1] /= skalar; vector[2] /= skalar;
}

template<class T> 
void vec3d<T>::vector_product(vec3d<T> const& vec) { //Vektorprodukt bilden
    T x = vector[1]*vec.vector[2]-vector[2]*vec.vector[1];
    T y = vector[2]*vec.vector[0]-vector[0]*vec.vector[2];
    vector[2] = vector[0]*vec.vector[1]-vector[1]*vec.vector[0];
    vector[0] = x; vector[1] = y;
}

template<class T> 
const T vec3d<T>::skalar_product(vec3d<T> const& vec) { //Skalarprodukt zurckgeben
    return (vector[0]*vec.vector[0] + vector[1]*vec.vector[1] + vector[2]*vec.vector[2]);
}

template<class T> 
void vec3d<T>::transform(matrix<T> const& mat) { //Vektor mit matrix multiplizieren
    T r[3] = {0.,0.,0.};
    r[0]+=mat.mat[0][0]*vector[0]; r[0]+=mat.mat[0][1]*vector[1]; r[0]+=mat.mat[0][2]*vector[2];
    r[1]+=mat.mat[1][0]*vector[0]; r[1]+=mat.mat[1][1]*vector[1]; r[1]+=mat.mat[1][2]*vector[2];
    r[2]+=mat.mat[2][0]*vector[0]; r[2]+=mat.mat[2][1]*vector[1]; r[2]+=mat.mat[2][2]*vector[2];
    vector[0] = r[0]; vector[1] = r[1]; vector[2] = r[2];
}

template<class T> 
void vec3d<T>::transform(matrix44<T> const& mat) { //Vektor mit matrix multiplizieren
    T r[3] = {0.,0.,0.};
    r[0]+=mat.mat[0][0]*vector[0]; r[0]+=mat.mat[0][1]*vector[1]; r[0]+=mat.mat[0][2]*vector[2]; r[0]+=mat.mat[0][3];
    r[1]+=mat.mat[1][0]*vector[0]; r[1]+=mat.mat[1][1]*vector[1]; r[1]+=mat.mat[1][2]*vector[2]; r[1]+=mat.mat[1][3];
    r[2]+=mat.mat[2][0]*vector[0]; r[2]+=mat.mat[2][1]*vector[1]; r[2]+=mat.mat[2][2]*vector[2]; r[2]+=mat.mat[2][3];
    vector[0] = r[0]; vector[1] = r[1]; vector[2] = r[2];
}

template<class T>
T* const vec3d<T>::get_vec() const {
    return vector;
}

template<class T>
T* vec3d<T>::get_vec() {
    return vector;
}


template<class T> 
const vec3d<T>& vec3d<T>::operator=(vec3d<T> const& rechts) { //Zuweisungsoperator berladen
    if (this == &rechts) return *this; //Zeit sparen bei Selbstzuweisung!
    vector[0] = rechts.vector[0]; vector[1] = rechts.vector[1]; vector[2] = rechts.vector[2];
    return *this;
}

template<class T> 
const vec3d<T>& vec3d<T>::operator=(T rechts[]) {
    vector[0] = rechts[0]; vector[1] = rechts[1]; vector[2] = rechts[2];
    return *this;
}

template<class T> 
const vec3d<T>& vec3d<T>::operator+=(vec3d<T> const& rechts) {
    add(rechts);
    return *this;
}

template<class T> 
const vec3d<T>& vec3d<T>::operator-=(vec3d<T> const& rechts) {
    sub(rechts);
    return *this;
}

template<class T> template<class T2> 
const vec3d<T>& vec3d<T>::operator*=(T2 const& rechts) {
    mul_skalar(rechts);
    return *this;
}

template<class T> template<class T2> 
const vec3d<T>& vec3d<T>::operator/=(T2 const& rechts) {
    div_skalar(rechts);
    return *this;
}

template<class T>
const vec3d<T>& vec3d<T>::operator*=(vec3d<T> const& rechts) {
    vector_product(rechts);
    return *this;
}

template<class T>
const vec3d<T>& vec3d<T>::operator*=(QUATERNION<T> const& rechts) {
    quat_mul(rechts);
    return *this;
}

template<class T> 
const vec3d<T>& vec3d<T>::operator*=(matrix<T> const& rechts) {
    transform(rechts);
    return *this;
}

template<class T> 
const vec3d<T>& vec3d<T>::operator*=(matrix44<T> const& rechts) {
    transform(rechts);
    return *this;
}

template<class T> 
T& vec3d<T>::operator[](int const& index) { // [] Operator fuer Lesen und Schreiben
    return vector[index];
}

template<class T> 
const T& vec3d<T>::operator[](int const& index) const { // [] Operator fuer Lesen und Schreiben
    return vector[index];
}

template<class T> 
vec3d<T>::operator T*() {
    return vector;
}

template<class T> 
bool vec3d<T>::operator<(vec3d<T> const& rechts) {
    if (vector[0] < rechts.vector[0] && vector[1] < rechts.vector[1] && vector[2] < rechts.vector[2]) return true;
    return false;
}

template<class T> 
bool vec3d<T>::operator>(vec3d<T> const& rechts) {
    if (vector[0] > rechts.vector[0] && vector[1] > rechts.vector[1] && vector[2] > rechts.vector[2]) return true;
    return false;
}

template<class T> 
bool vec3d<T>::operator<=(vec3d<T> const& rechts) {
    if (vector[0] <= rechts.vector[0] && vector[1] <= rechts.vector[1] && vector[2] <= rechts.vector[2]) return true;
    return false;
}

template<class T> 
bool vec3d<T>::operator>=(vec3d<T> const& rechts) {
    if (vector[0] >= rechts.vector[0] && vector[1] >= rechts.vector[1] && vector[2] >= rechts.vector[2]) return true;
    return false;
}


//==============================================================================================
//!Definitionen fuer matrix:
//==============================================================================================

template<class T> 
matrix<T>::matrix(T const& a00,T const& a01,T const& a02,
                  T const& a10,T const& a11,T const& a12,
                  T const& a20,T const& a21,T const& a22) {
    mat[0][0] = a00; mat[0][1] = a01; mat[0][2] = a02;
    mat[1][0] = a10; mat[1][1] = a11; mat[1][2] = a12;
    mat[2][0] = a20; mat[2][1] = a21; mat[2][2] = a22;
}

template<class T> 
matrix<T>::matrix(matrix<T> const& m) {
    mat[0][0] = m.mat[0][0]; mat[0][1] = m.mat[0][1]; mat[0][2] = m.mat[0][2];
    mat[1][0] = m.mat[1][0]; mat[1][1] = m.mat[1][1]; mat[1][2] = m.mat[1][2];
    mat[2][0] = m.mat[2][0]; mat[2][1] = m.mat[2][1]; mat[2][2] = m.mat[2][2];
    
}

template<class T> 
matrix<T>::~matrix() {}


template<class T> 
const T matrix<T>::det() { //Determinante berechnen
    return ((mat[0][0]*mat[1][1]*mat[2][2])+
            (mat[0][1]*mat[1][2]*mat[2][0])+
            (mat[0][2]*mat[1][0]*mat[2][1])-
            (mat[2][0]*mat[1][1]*mat[0][2])-
            (mat[2][2]*mat[1][0]*mat[0][1])-
            (mat[0][0]*mat[2][1]*mat[1][2]));
}


template<class T>
bool matrix<T>::invert() { //Matrix invertieren
    T determinante = det();
    if (fabs(determinante) < 0.000001) return false;
    
    matrix<T> dummy(*this);
    mat[0][0] = (dummy.mat[1][1] * dummy.mat[2][2] - dummy.mat[1][2] * dummy.mat[2][1]) / determinante;
    mat[0][1] = (dummy.mat[0][2] * dummy.mat[2][1] - dummy.mat[0][1] * dummy.mat[2][2]) / determinante;
    mat[0][2] = (dummy.mat[0][1] * dummy.mat[1][2] - dummy.mat[0][2] * dummy.mat[1][1]) / determinante;
    mat[1][0] = (dummy.mat[1][2] * dummy.mat[2][0] - dummy.mat[1][0] * dummy.mat[2][2]) / determinante;
    mat[1][1] = (dummy.mat[0][0] * dummy.mat[2][2] - dummy.mat[0][2] * dummy.mat[2][0]) / determinante;
    mat[1][2] = (dummy.mat[0][2] * dummy.mat[1][0] - dummy.mat[0][0] * dummy.mat[1][2]) / determinante;
    mat[2][0] = (dummy.mat[1][0] * dummy.mat[2][1] - dummy.mat[1][1] * dummy.mat[2][0]) / determinante;
    mat[2][1] = (dummy.mat[0][1] * dummy.mat[2][0] - dummy.mat[0][0] * dummy.mat[2][1]) / determinante;
    mat[2][2] = (dummy.mat[0][0] * dummy.mat[1][1] - dummy.mat[0][1] * dummy.mat[1][0]) / determinante;
    
    return true;
}

template<class T>
bool matrix<T>::gaussj_invert() { //Matrix invertieren
    //! invert() wird sehr ungenau, wenn die Determinante nahe Null geht
    //! Daher alternativ diese deutlich umstaendlichere Variante
    //! (Gauss-Jordan Eliminierung)
    
    //! Die Schleifen bei Gelegenheit alle ersetzen, da der Algo nur fuer den
    //! 3x3 Spezialfall gebraucht wird !!!
    
    double a[3][3] = {{mat[0][0],mat[0][1],mat[0][2]},
                      {mat[1][0],mat[1][1],mat[1][2]},
                      {mat[2][0],mat[2][1],mat[2][2]}};
    /*
    double b[3][3] = {{1.,0.,0.},
                      {0.,1.,0.},
                      {0.,0.,1.}};
    */
    int indxc[3],indxr[3],ipiv[3];
    int i = 0;
    int icol = 0;
    int irow = 0;
    int j = 0;
    int k = 0;
    int l = 0;
    int ll = 0;
    double big,dum,pivinv,temp;
    
    for (j=0; j<3; ++j) ipiv[j] = 0;
    for (i=0; i<3; ++i) {
        big = 0.;
        for (j=0; j<3; ++j) {
            if (ipiv[j] != 1) for (k=0; k<3; ++k) {
                if (ipiv[k] == 0) {
                    if (fabs(a[j][k]) >= big) {
                        big = fabs(a[j][k]);
                        irow = j;
                        icol = k;
                    }
                } else if (ipiv[k] > 1) return false;
            }
        }
        ++(ipiv[icol]);
        
        if (irow != icol) {
            for (l=0; l<3; ++l) {temp = a[irow][l]; a[irow][l] = a[icol][l]; a[icol][l] = temp;}
        //    for (l=0; l<3; ++l) {temp = b[irow][l]; b[irow][l] = b[icol][l]; b[icol][l] = temp;}
        }
        
        indxr[i] = irow;
        indxc[i] = icol;
        if (a[icol][icol] == 0.) return false;
        pivinv = 1. / a[icol][icol];
        a[icol][icol] = 1.;
        for (l=0; l<3; ++l) a[icol][l] *= pivinv;
    //    for (l=0; l<3; ++l) b[icol][l] *= pivinv;
        
        for (ll=0; ll<3; ++ll) {
            if (ll != icol) {
                dum = a[ll][icol];
                a[ll][icol] = 0.;
                for (l=0; l<3; ++l) a[ll][l] -= a[icol][l] * dum;
            //    for (l=0; l<3; ++l) b[ll][l] -= b[icol][l] * dum;
            }
        }
    }
    
    for (l=2; l>=0; --l) {
        if (indxr[l] != indxc[l]) for (k=0; k<3; ++k) {
            temp = a[k][indxr[l]]; a[k][indxr[l]] = a[k][indxc[l]]; a[k][indxc[l]] = temp;
        }
    }
    
    mat[0][0] = a[0][0]; mat[0][1] = a[0][1]; mat[0][2] = a[0][2];
    mat[1][0] = a[1][0]; mat[1][1] = a[1][1]; mat[1][2] = a[1][2];
    mat[2][0] = a[2][0]; mat[2][1] = a[2][1]; mat[2][2] = a[2][2];
    
    
    return true;
}

template<class T>
bool matrix<T>::sqrt_gj() { //Wurzel einer Matrix nach Denman-Beavers
    int iter = 0;
    matrix<T> Y(*this);
    matrix<T> Z(1.,0.,0.,
                0.,1.,0.,
                0.,0.,1.);
    matrix<T> invY,invZ,Y1,Z1;
    matrix<T> tmat;
    do {
        invY = Y;
        invZ = Z;
        if (!invY.gaussj_invert()) break;//return false;
        if (!invZ.gaussj_invert()) break;//return false;
        
        //!hier gehts weiter:
        Y1.mat[0][0] = 0.5 * (Y.mat[0][0] + invZ[0][0]);
        Z1.mat[0][0] = 0.5 * (Z.mat[0][0] + invY[0][0]);
        
        Y1.mat[0][1] = 0.5 * (Y.mat[0][1] + invZ[0][1]);
        Z1.mat[0][1] = 0.5 * (Z.mat[0][1] + invY[0][1]);
        
        Y1.mat[0][2] = 0.5 * (Y.mat[0][2] + invZ[0][2]);
        Z1.mat[0][2] = 0.5 * (Z.mat[0][2] + invY[0][2]);
        
        Y1.mat[1][0] = 0.5 * (Y.mat[1][0] + invZ[1][0]);
        Z1.mat[1][0] = 0.5 * (Z.mat[1][0] + invY[1][0]);
        
        Y1.mat[1][1] = 0.5 * (Y.mat[1][1] + invZ[1][1]);
        Z1.mat[1][1] = 0.5 * (Z.mat[1][1] + invY[1][1]);
        
        Y1.mat[1][2] = 0.5 * (Y.mat[1][2] + invZ[1][2]);
        Z1.mat[1][2] = 0.5 * (Z.mat[1][2] + invY[1][2]);
        
        Y1.mat[2][0] = 0.5 * (Y.mat[2][0] + invZ[2][0]);
        Z1.mat[2][0] = 0.5 * (Z.mat[2][0] + invY[2][0]);
        
        Y1.mat[2][1] = 0.5 * (Y.mat[2][1] + invZ[2][1]);
        Z1.mat[2][1] = 0.5 * (Z.mat[2][1] + invY[2][1]);
        
        Y1.mat[2][2] = 0.5 * (Y.mat[2][2] + invZ[2][2]);
        Z1.mat[2][2] = 0.5 * (Z.mat[2][2] + invY[2][2]);
        
        Y = Y1;
        Z = Z1;
        
        tmat = Y;
        tmat.mul_matrix(Y);
        
        ++iter;
    } while (!is_almost_equal(*this,tmat) && iter < 20);
    
    *this = Y;
    return true;
}

template<class T>
void matrix<T>::transpone() { //Matrix transponieren
    matrix<T> dummy(*this);
    mat[0][0] = dummy.mat[0][0]; mat[0][1] = dummy.mat[1][0]; mat[0][2] = dummy.mat[2][0];
    mat[1][0] = dummy.mat[0][1]; mat[1][1] = dummy.mat[1][1]; mat[1][2] = dummy.mat[2][1];
    mat[2][0] = dummy.mat[0][2]; mat[2][1] = dummy.mat[1][2]; mat[2][2] = dummy.mat[2][2];
}

template<class T> 
void matrix<T>::mul_matrix(matrix<T> const& rechts) { //Matrixmultiplikation
    matrix<T> l(*this);
    mat[0][0] = l.mat[0][0] * rechts.mat[0][0] + l.mat[0][1] * rechts.mat[1][0] + l.mat[0][2] * rechts.mat[2][0];
    mat[0][1] = l.mat[0][0] * rechts.mat[0][1] + l.mat[0][1] * rechts.mat[1][1] + l.mat[0][2] * rechts.mat[2][1];
    mat[0][2] = l.mat[0][0] * rechts.mat[0][2] + l.mat[0][1] * rechts.mat[1][2] + l.mat[0][2] * rechts.mat[2][2];
    mat[1][0] = l.mat[1][0] * rechts.mat[0][0] + l.mat[1][1] * rechts.mat[1][0] + l.mat[1][2] * rechts.mat[2][0];
    mat[1][1] = l.mat[1][0] * rechts.mat[0][1] + l.mat[1][1] * rechts.mat[1][1] + l.mat[1][2] * rechts.mat[2][1];
    mat[1][2] = l.mat[1][0] * rechts.mat[0][2] + l.mat[1][1] * rechts.mat[1][2] + l.mat[1][2] * rechts.mat[2][2];
    mat[2][0] = l.mat[2][0] * rechts.mat[0][0] + l.mat[2][1] * rechts.mat[1][0] + l.mat[2][2] * rechts.mat[2][0];
    mat[2][1] = l.mat[2][0] * rechts.mat[0][1] + l.mat[2][1] * rechts.mat[1][1] + l.mat[2][2] * rechts.mat[2][1];
    mat[2][2] = l.mat[2][0] * rechts.mat[0][2] + l.mat[2][1] * rechts.mat[1][2] + l.mat[2][2] * rechts.mat[2][2];
}

template<class T> 
void matrix<T>::get_rot_v(T a,vec3d<T> *rv) {
    T helper = 1. - cos(a);
    rv[0][0] = sqrt((mat[0][0] - cos(a)) / helper);
    rv[0][1] = sqrt((mat[1][1] - cos(a)) / helper);
    rv[0][2] = sqrt((mat[2][2] - cos(a)) / helper);
    
    for (int i=1; i<8; ++i) rv[i] = rv[0];
    
    rv[1][0] *= -1.; //-++
    rv[2][1] *= -1.; //+-+
    rv[3][2] *= -1.; //++-
    rv[4][1] *= -1.; rv[4][2] *= -1.; //+--
    rv[5][0] *= -1.; rv[5][2] *= -1.; //-+-
    rv[6][0] *= -1.; rv[6][1] *= -1.; //--+
    rv[7][0] *= -1.; rv[7][1] *= -1.; rv[7][2] *= -1.; //---
}

template<class T>
T matrix<T>::get_rot_angle() {
    return acos((mat[0][0] + mat[1][1] + mat[2][2] - 1.) / 2.);
}

template<class T>
vec3d<T> matrix<T>::get_rot_vector() {
    T a1 = acos((mat[0][0] + mat[1][1] + mat[2][2] - 1.) / 2.);
    //Jetzt rausfinden welcher von beiden Winkeln der richtige ist:
    vec3d<T> cv[8];
    get_rot_v(a1,cv);
    int besta1 = 0;
    T a1diff = 999999999.;
    for (int i=0; i<8; ++i) {
        T curr_diff = 0.;
        matrix<T> mt = rotmatrix(cv[i],a1);
        for (int x=0; x<3; ++x) {
            for (int y=0; y<3; ++y) {
                curr_diff += abs(mt[x][y] - mat[x][y]);
            }
        }
        if (curr_diff < a1diff) {a1diff = curr_diff; besta1 = i;}
    }
    return cv[besta1];
}

template<class T>
void matrix<T>::get_cardan_angles(float &a,float &b,float &c) {
    //! Berechnet zur gegebenen Rotationsmatrix 3 Kardanwinkel
    //! Drehung gegen den Uhrzeigersinn nacheinander um z,y,x (c,b,a) sollte der gleichen
    //! Rotation wie durch die Rotationsmatrix entsprechen.
    //! Nach: http://mca.bv.tu-berlin.de/Sites/German_Sites/seminararbeiten/helmert.pdf
    //! Funktioniert irgendwie nicht!
    //! Am besten generell die Finger lassen von Cardan und Fixed angles
    //!   -> immer das Problem mit Singularitaeten und der Mehrdeutigkeit
    //!   -> lieber nur mit Rotationsmatrizen oder Quaternionen arbeiten
    b = asin(mat[2][0]);
    if (cos(b) < 0.0001) {
        c = 0.;
        a = asin(mat[0][1]);
    } else {
        a = atan(-mat[2][0]/mat[2][2]);
        c = atan(-mat[1][0]/mat[0][0]);
    }
}

template<class T>
void matrix<T>::get_fixed_angles(float &x,float &y,float &z) {
    //! Berechnet zur gegebenen Rotationsmatrix 3 fixed angle fuer x,y,z Drehung
    //! Nach: Handbook of Robotics Part A 1.2
    /*
    cerr << mat[0][0]+mat[1][0] << endl;
    y = atan2(-mat[2][0],sqrt(mat[0][0]+mat[1][0]));
    float cosy = cos(y);
    if (cosy < 0.001) { // Singularitaet
        y = 0.001; cosy = 1.5697963;
    }
    x = atan2(mat[1][0]/cosy,mat[0][0]/cosy);
    z = atan2(mat[2][1]/cosy,mat[2][2]/cosy);
    */
    //! Auch hier kommt eher Schwachsinn raus => lieber nur mit
    //! Rotationsmatrizen oder Quaternionen arbeiten
    x = atan2(mat[2][1],mat[2][2]);
    z = atan2(mat[1][0],mat[0][0]);
    y = atan2(-mat[2][0],cos(z)*mat[0][0]+sin(z)*mat[1][0]);
}


template<class T> 
const matrix<T>& matrix<T>::operator=(matrix<T> const& rechts) {
    if (this == &rechts) return *this; //Zeit sparen bei Selbstzuweisung!
    mat[0][0] = rechts.mat[0][0]; mat[0][1] = rechts.mat[0][1]; mat[0][2] = rechts.mat[0][2];
    mat[1][0] = rechts.mat[1][0]; mat[1][1] = rechts.mat[1][1]; mat[1][2] = rechts.mat[1][2];
    mat[2][0] = rechts.mat[2][0]; mat[2][1] = rechts.mat[2][1]; mat[2][2] = rechts.mat[2][2];
    return *this;
}

template<class T> 
INDEXER<T> matrix<T>::operator[](int const& index) { // [] Operator fuer Lesen und Schreiben
    return INDEXER<T>(index,this);
}

template<class T> 
const matrix<T>& matrix<T>::operator*=(matrix<T> const& rechts) {
    mul_matrix(rechts);
    return *this;
}

//==============================================================================================
//!Definitionen fuer matrix44:
//==============================================================================================
//--------------------------------------------------------------------------
//Konstruktoren und Destruktor fr matrix44:

template<class T> 
matrix44<T>::matrix44(const T& a00,const T& a01,const T& a02,const T& a03,
                      const T& a10,const T& a11,const T& a12,const T& a13,
                      const T& a20,const T& a21,const T& a22,const T& a23,
                      const T& a30,const T& a31,const T& a32,const T& a33) {
    mat[0][0] = a00; mat[0][1] = a01; mat[0][2] = a02; mat[0][3] = a03;
    mat[1][0] = a10; mat[1][1] = a11; mat[1][2] = a12; mat[1][3] = a13;
    mat[2][0] = a20; mat[2][1] = a21; mat[2][2] = a22; mat[2][3] = a23;
    mat[3][0] = a30; mat[3][1] = a31; mat[3][2] = a32; mat[3][3] = a33;
}

template<class T> 
matrix44<T>::matrix44(matrix44<T> const& m) {
    mat[0][0] = m.mat[0][0]; mat[0][1] = m.mat[0][1]; mat[0][2] = m.mat[0][2]; mat[0][3] = m.mat[0][3];
    mat[1][0] = m.mat[1][0]; mat[1][1] = m.mat[1][1]; mat[1][2] = m.mat[1][2]; mat[1][3] = m.mat[1][3];
    mat[2][0] = m.mat[2][0]; mat[2][1] = m.mat[2][1]; mat[2][2] = m.mat[2][2]; mat[2][3] = m.mat[2][3];
    mat[3][0] = m.mat[3][0]; mat[3][1] = m.mat[3][1]; mat[3][2] = m.mat[3][2]; mat[3][3] = m.mat[3][3];
}

template<class T> 
matrix44<T>::~matrix44() {}


template<class T> 
const T matrix44<T>::sub_det(int row,int column) { //Unterdeterminante berechnen
    int r_i[3]; int c_i[3];
    int ir = 0;
    for (int i=0; i<4; ++i) {
        if (i == row) continue;
        r_i[ir] = i;
        ++ir;
    }
    ir = 0;
    for (int i=0; i<4; ++i) {
        if (i == column) continue;
        c_i[ir] = i;
        ++ir;
    }
    return ((mat[r_i[0]][c_i[0]]*mat[r_i[1]][c_i[1]]*mat[r_i[2]][c_i[2]])+
            (mat[r_i[0]][c_i[1]]*mat[r_i[1]][c_i[2]]*mat[r_i[2]][c_i[0]])+
            (mat[r_i[0]][c_i[2]]*mat[r_i[1]][c_i[0]]*mat[r_i[2]][c_i[1]])-
            (mat[r_i[2]][c_i[0]]*mat[r_i[1]][c_i[1]]*mat[r_i[0]][c_i[2]])-
            (mat[r_i[2]][c_i[2]]*mat[r_i[1]][c_i[0]]*mat[r_i[0]][c_i[1]])-
            (mat[r_i[0]][c_i[0]]*mat[r_i[2]][c_i[1]]*mat[r_i[1]][c_i[2]]));
}

template<class T> 
const T matrix44<T>::det() { //Determinante berechnen
    return ((mat[0][0] * sub_det(0,0)) -
            (mat[0][1] * sub_det(0,1)) +
            (mat[0][2] * sub_det(0,2)) -
            (mat[0][3] * sub_det(0,3)));
}

template<class T>
bool matrix44<T>::invert() {
    T determinante = det();
    if (fabs(determinante) < 0.00001) return false;
    
    matrix44<T> dummy(*this);
    
    mat[0][0] = dummy.sub_det(0,0) / determinante;
    mat[0][1] = dummy.sub_det(1,0) / determinante;
    mat[0][2] = dummy.sub_det(2,0) / determinante;
    mat[0][3] = dummy.sub_det(3,0) / determinante;
    
    mat[1][0] = dummy.sub_det(0,1) / determinante;
    mat[1][1] = dummy.sub_det(1,1) / determinante;
    mat[1][2] = dummy.sub_det(2,1) / determinante;
    mat[1][3] = dummy.sub_det(3,1) / determinante;
    
    mat[2][0] = dummy.sub_det(0,2) / determinante;
    mat[2][1] = dummy.sub_det(1,2) / determinante;
    mat[2][2] = dummy.sub_det(2,2) / determinante;
    mat[2][3] = dummy.sub_det(3,2) / determinante;
    
    mat[3][0] = dummy.sub_det(0,3) / determinante;
    mat[3][1] = dummy.sub_det(1,3) / determinante;
    mat[3][2] = dummy.sub_det(2,3) / determinante;
    mat[3][3] = dummy.sub_det(3,3) / determinante;
    
    return true;
}

template<class T> 
void matrix44<T>::mul_matrix(matrix44<T> const& rechts) { //Matrixmultiplikation
    matrix44<T> l(*this);
    mat[0][0] = l.mat[0][0]*rechts.mat[0][0]+l.mat[0][1]*rechts.mat[1][0]+l.mat[0][2]*rechts.mat[2][0]+l.mat[0][3]*rechts.mat[3][0];
    mat[0][1] = l.mat[0][0]*rechts.mat[0][1]+l.mat[0][1]*rechts.mat[1][1]+l.mat[0][2]*rechts.mat[2][1]+l.mat[0][3]*rechts.mat[3][1];
    mat[0][2] = l.mat[0][0]*rechts.mat[0][2]+l.mat[0][1]*rechts.mat[1][2]+l.mat[0][2]*rechts.mat[2][2]+l.mat[0][3]*rechts.mat[3][2];
    mat[0][3] = l.mat[0][0]*rechts.mat[0][3]+l.mat[0][1]*rechts.mat[1][3]+l.mat[0][2]*rechts.mat[2][3]+l.mat[0][3]*rechts.mat[3][3];
    
    mat[1][0] = l.mat[1][0]*rechts.mat[0][0]+l.mat[1][1]*rechts.mat[1][0]+l.mat[1][2]*rechts.mat[2][0]+l.mat[1][3]*rechts.mat[3][0];
    mat[1][1] = l.mat[1][0]*rechts.mat[0][1]+l.mat[1][1]*rechts.mat[1][1]+l.mat[1][2]*rechts.mat[2][1]+l.mat[1][3]*rechts.mat[3][1];
    mat[1][2] = l.mat[1][0]*rechts.mat[0][2]+l.mat[1][1]*rechts.mat[1][2]+l.mat[1][2]*rechts.mat[2][2]+l.mat[1][3]*rechts.mat[3][2];
    mat[1][3] = l.mat[1][0]*rechts.mat[0][3]+l.mat[1][1]*rechts.mat[1][3]+l.mat[1][2]*rechts.mat[2][3]+l.mat[1][3]*rechts.mat[3][3];
    
    mat[2][0] = l.mat[2][0]*rechts.mat[0][0]+l.mat[2][1]*rechts.mat[1][0]+l.mat[2][2]*rechts.mat[2][0]+l.mat[2][3]*rechts.mat[3][0];
    mat[2][1] = l.mat[2][0]*rechts.mat[0][1]+l.mat[2][1]*rechts.mat[1][1]+l.mat[2][2]*rechts.mat[2][1]+l.mat[2][3]*rechts.mat[3][1];
    mat[2][2] = l.mat[2][0]*rechts.mat[0][2]+l.mat[2][1]*rechts.mat[1][2]+l.mat[2][2]*rechts.mat[2][2]+l.mat[2][3]*rechts.mat[3][2];
    mat[2][3] = l.mat[2][0]*rechts.mat[0][3]+l.mat[2][1]*rechts.mat[1][3]+l.mat[2][2]*rechts.mat[2][3]+l.mat[2][3]*rechts.mat[3][3];
    
    mat[3][0] = l.mat[3][0]*rechts.mat[0][0]+l.mat[3][1]*rechts.mat[1][0]+l.mat[3][2]*rechts.mat[2][0]+l.mat[3][3]*rechts.mat[3][0];
    mat[3][1] = l.mat[3][0]*rechts.mat[0][1]+l.mat[3][1]*rechts.mat[1][1]+l.mat[3][2]*rechts.mat[2][1]+l.mat[3][3]*rechts.mat[3][1];
    mat[3][2] = l.mat[3][0]*rechts.mat[0][2]+l.mat[3][1]*rechts.mat[1][2]+l.mat[3][2]*rechts.mat[2][2]+l.mat[3][3]*rechts.mat[3][2];
    mat[3][3] = l.mat[3][0]*rechts.mat[0][3]+l.mat[3][1]*rechts.mat[1][3]+l.mat[3][2]*rechts.mat[2][3]+l.mat[3][3]*rechts.mat[3][3];
}

//--------------------------------------------------------------------------
//Klassen-Operatoren fr matrix44:

template<class T> 
const matrix44<T>& matrix44<T>::operator=(matrix44<T> const& rechts) {
    if (this == &rechts) return *this; //Zeit sparen bei Selbstzuweisung!
    mat[0][0] = rechts.mat[0][0]; mat[0][1] = rechts.mat[0][1]; mat[0][2] = rechts.mat[0][2]; mat[0][3] = rechts.mat[0][3];
    mat[1][0] = rechts.mat[1][0]; mat[1][1] = rechts.mat[1][1]; mat[1][2] = rechts.mat[1][2]; mat[1][3] = rechts.mat[1][3];
    mat[2][0] = rechts.mat[2][0]; mat[2][1] = rechts.mat[2][1]; mat[2][2] = rechts.mat[2][2]; mat[2][3] = rechts.mat[2][3];
    mat[3][0] = rechts.mat[3][0]; mat[3][1] = rechts.mat[3][1]; mat[3][2] = rechts.mat[3][2]; mat[3][3] = rechts.mat[3][3];
    return *this;
}

template<class T> 
INDEXER44<T> matrix44<T>::operator[](int const& index) { // [] Operator fuer Lesen und Schreiben
    return INDEXER44<T>(index,this);
}


template<class T> 
const matrix44<T>& matrix44<T>::operator*=(matrix44<T> const& rechts) {
    mul_matrix(rechts);
    return *this;
}

//==============================================================================================
//!Definitionen fuer matrix55:
//==============================================================================================

template<class T> 
matrix55<T>::matrix55(const T& a00,const T& a01,const T& a02,const T& a03,const T& a04, 
                      const T& a10,const T& a11,const T& a12,const T& a13,const T& a14,
                      const T& a20,const T& a21,const T& a22,const T& a23,const T& a24,
                      const T& a30,const T& a31,const T& a32,const T& a33,const T& a34,
                      const T& a40,const T& a41,const T& a42,const T& a43,const T& a44) {
    mat[0][0] = a00; mat[0][1] = a01; mat[0][2] = a02; mat[0][3] = a03; mat[0][4] = a04;
    mat[1][0] = a10; mat[1][1] = a11; mat[1][2] = a12; mat[1][3] = a13; mat[1][4] = a14;
    mat[2][0] = a20; mat[2][1] = a21; mat[2][2] = a22; mat[2][3] = a23; mat[2][4] = a24;
    mat[3][0] = a30; mat[3][1] = a31; mat[3][2] = a32; mat[3][3] = a33; mat[3][4] = a34;
    mat[4][0] = a40; mat[4][1] = a41; mat[4][2] = a42; mat[4][3] = a43; mat[4][4] = a44;
}

template<class T> 
matrix55<T>::matrix55(matrix55<T> const& m) {
    mat[0][0] = m.mat[0][0]; mat[0][1] = m.mat[0][1]; mat[0][2] = m.mat[0][2]; mat[0][3] = m.mat[0][3]; mat[0][4] = m.mat[0][4];
    mat[1][0] = m.mat[1][0]; mat[1][1] = m.mat[1][1]; mat[1][2] = m.mat[1][2]; mat[1][3] = m.mat[1][3]; mat[1][4] = m.mat[1][4];
    mat[2][0] = m.mat[2][0]; mat[2][1] = m.mat[2][1]; mat[2][2] = m.mat[2][2]; mat[2][3] = m.mat[2][3]; mat[2][4] = m.mat[2][4];
    mat[3][0] = m.mat[3][0]; mat[3][1] = m.mat[3][1]; mat[3][2] = m.mat[3][2]; mat[3][3] = m.mat[3][3]; mat[3][4] = m.mat[3][4];
    mat[4][0] = m.mat[4][0]; mat[4][1] = m.mat[4][1]; mat[4][2] = m.mat[4][2]; mat[4][3] = m.mat[4][3]; mat[4][4] = m.mat[4][4];
}

template<class T> 
matrix55<T>::~matrix55() {}


template<class T> 
const T matrix55<T>::sub_det(int row,int column) { //Unterdeterminante berechnen
    int r_i[4]; int c_i[4];
    int ir = 0;
    for (int i=0; i<5; ++i) {
        if (i == row) continue;
        r_i[ir] = i;
        ++ir;
    }
    ir = 0;
    for (int i=0; i<5; ++i) {
        if (i == column) continue;
        c_i[ir] = i;
        ++ir;
    }
    return (matrix44<T>(mat[r_i[0]][c_i[0]],mat[r_i[0]][c_i[1]],mat[r_i[0]][c_i[2]],mat[r_i[0]][c_i[3]],
                        mat[r_i[1]][c_i[0]],mat[r_i[1]][c_i[1]],mat[r_i[1]][c_i[2]],mat[r_i[1]][c_i[3]],
                        mat[r_i[2]][c_i[0]],mat[r_i[2]][c_i[1]],mat[r_i[2]][c_i[2]],mat[r_i[2]][c_i[3]],
                        mat[r_i[3]][c_i[0]],mat[r_i[3]][c_i[1]],mat[r_i[3]][c_i[2]],mat[r_i[3]][c_i[3]]).det());
}

template<class T> 
const T matrix55<T>::det() { //Determinante berechnen
    return ((mat[0][0] * sub_det(0,0)) -
            (mat[0][1] * sub_det(0,1)) +
            (mat[0][2] * sub_det(0,2)) -
            (mat[0][3] * sub_det(0,3)) +
            (mat[0][4] * sub_det(0,4)));
}


template<class T> 
INDEXER<T> matrix55<T>::operator[](int const& index) { // [] Operator fuer Lesen und Schreiben
    return INDEXER<T>(index,this);
}


//==============================================================================================
//!Definitionen fuer INDEXER:
//==============================================================================================

template<class T> 
INDEXER<T>::INDEXER(int yindex, matrix<T> *matr) : y(yindex), m(matr) {}

template<class T> 
INDEXER<T>::~INDEXER() {}

template<class T>
T& INDEXER<T>::operator[](int const& xindex) {
    return m->mat[y][xindex];
}

template<class T> 
INDEXER44<T>::INDEXER44(int yindex, matrix44<T> *matr) : y(yindex), m(matr) {}

template<class T> 
INDEXER44<T>::~INDEXER44() {}

template<class T>
T& INDEXER44<T>::operator[](int const& xindex) {
    return m->mat[y][xindex];
}

//==============================================================================================
//!globale Definitionen:
//==============================================================================================

template<class T> 
inline const T skalarproduct(vec3d<T> const& links,vec3d<T> const& rechts) { //bildet das Skalarprodukt aus 2 Vektoren
    return (links.vector[0]*rechts.vector[0] +
            links.vector[1]*rechts.vector[1] + 
            links.vector[2]*rechts.vector[2]);
}

template<class T> 
inline const vec3d<T> vectorproduct(vec3d<T> const& links,vec3d<T> const& rechts) { //bildet das Vektorprodukt aus 2 Vektoren
    return (vec3d<T>(links) *= rechts);
}

template<class T> 
inline const T angle(vec3d<T> const& links,vec3d<T> const& rechts) { //liefert den Winkel zwischen 2 Vektoren
    return acos(skalarproduct(links,rechts)/(links.value()*rechts.value()));
}

template<class T> 
inline const T angle_for_normed(vec3d<T> const& links,vec3d<T> const& rechts) { //liefert den Winkel zwischen 2 Vektoren
    //! Fuer normierte Vektoren muss nicht nochmal extra value() aufgerufen werden!
    return acos(skalarproduct(links,rechts));
}

template<class T> 
inline const T angle(vec3d<T> const& p1,vec3d<T> const& p2,vec3d<T> const& p3) { //liefert den Winkel zwischen 2 Vektoren p2-p1 und p2-p3
    vec3d<T> v1(p1); v1 -= p2;
    vec3d<T> v2(p3); v2 -= p2;
    return acos(skalarproduct(v1,v2)/(v1.value()*v2.value()));
}

template<class T> 
inline const T get_distance(vec3d<T> const& v1,vec3d<T> const& v2) {
    T x = v1.vector[0] - v2.vector[0];
    T y = v1.vector[1] - v2.vector[1];
    T z = v1.vector[2] - v2.vector[2];
    return sqrt(x*x+y*y+z*z);
}

template<class T> 
inline const T get_square_distance(vec3d<T> const& v1,vec3d<T> const& v2) {
    T x = v1.vector[0] - v2.vector[0];
    T y = v1.vector[1] - v2.vector[1];
    T z = v1.vector[2] - v2.vector[2];
    return (x*x+y*y+z*z);
}

template<class T,class T2>
inline const T get_point2plane_distance(vec3d<T2> const& point,vec3d<T> const& p0,vec3d<T> const& p1,vec3d<T> const& p2) {
    vec3d<T> h1(p1); h1 -= p0;
    vec3d<T> h2(p2); h2 -= p0;
    h1 *= h2;
    h2 = point; h2 -= p0; T val = h2.value();
    val *= cos(angle(h1,h2));
    return abs(val);
}

template<class T> 
inline const matrix<T> rotmatrix(vec3d<T> const& vec,double const& angle) { //liefert die matrix fuer Rotation um vec
    vec3d<T> tvec = vec; tvec.norm(); //RotVektor muss normiert sein!!!
    T r[3][3]; //zum initialisieren der matrix
    //Werte vorberechnen:
    T v11=tvec.vector[0]*tvec.vector[0]; T v22=tvec.vector[1]*tvec.vector[1]; T v33=tvec.vector[2]*tvec.vector[2];
    T v12=tvec.vector[0]*tvec.vector[1]; T v13=tvec.vector[0]*tvec.vector[2]; T v23=tvec.vector[1]*tvec.vector[2];
    T c=cos(angle); T s=sin(angle); T x=1.0-c;
    T v1s=tvec.vector[0]*s; T v2s=tvec.vector[1]*s; T v3s=tvec.vector[2]*s;
    //Matrix berechnen
    r[0][0]=c+v11*x;
    r[0][1]=v12*x-v3s;
    r[0][2]=v13*x+v2s;
    r[1][0]=v12*x+v3s;
    r[1][1]=c+v22*x;
    r[1][2]=v23*x-v1s;
    r[2][0]=v13*x-v2s;
    r[2][1]=v23*x+v1s;
    r[2][2]=c+v33*x;
    return matrix<T>(r[0][0],r[0][1],r[0][2],
                     r[1][0],r[1][1],r[1][2],
                     r[2][0],r[2][1],r[2][2]);
}

template<class T> 
inline const T determinant(matrix<T> const& m) {
    return ((m.mat[1][1]*m.mat[2][2]*m.mat[3][3])+
            (m.mat[1][2]*m.mat[2][3]*m.mat[3][1])+
            (m.mat[1][3]*m.mat[2][1]*m.mat[3][2])-
            (m.mat[3][1]*m.mat[2][2]*m.mat[1][3])-
            (m.mat[3][3]*m.mat[2][1]*m.mat[1][2])-
            (m.mat[1][1]*m.mat[3][2]*m.mat[2][3]));
}

template<class T> 
inline const T dihedral(vec3d<T> const& r1,vec3d<T> const& r2,vec3d<T> const& r3,vec3d<T> const& r4) {
    //! liefert einen Winkel zwischen 0 und pi
    vec3d<T> v1(r2); v1 -= r1;
    vec3d<T> v2(r3); v2 -= r2;
    vec3d<T> v3(r4); v3 -= r3;
    v1 *= v2; v2 *= v3;
    return acos(skalarproduct(v1,v2)/(v1.value()*v2.value()));
}

template<class T>
inline const T p2p_angle(vec3d<T> const& e1p1,vec3d<T> const& e1p2,vec3d<T> const& e1p3,
                         vec3d<T> const& e2p1,vec3d<T> const& e2p2,vec3d<T> const& e2p3) {
    //! Winkel zwischen 2 Ebenen:
    //! ACHTUNG: Die Reihenfolge der Punkte entscheidet hier darueber welcher der 2 Winkel bestimmt wird!!!
    vec3d<T> v1 = e1p1; v1 -= e1p2;
    vec3d<T> v2 = e1p3; v2 -= e1p2;
    v1 *= v2;
    vec3d<T> v3 = e2p3; v3 -= e2p2;
    v2 = e2p1; v2 -= e2p2;
    v2 *= v3;
    return acos(skalarproduct(v1,v2)/(v1.value()*v2.value()));
}

template<class T>
inline const T circumcircle_radius(vec3d<T> const& p0,vec3d<T> const& p1,vec3d<T> const& p2) {
    //! Umkreisradius des Dreiecks: p0,p1,p2
    vec3d<T> a(p1); a -= p2;
    return (a.value() / (2. * sin(angle(p1,p0,p2))));
}

template<class T, class T2>
inline bool points_to_sphere(vec3d<T> const& p0,vec3d<T> const& p1,vec3d<T> const& p2,
                             vec3d<T> const& p3,vec3d<T> &center,T2 &sradius) {
    //!Der Funktion werden 4 Punkte, sowie die Referenz auf einen weiteren Punkt und eine Zahl uebergeben.
    //!Es wird die Umkugel fuer die 4 Punkte berechnet => Der Mittelpunkt wird in center geschrieben und 
    //!das Quadrat des Radius der Umkugel in radius.
    //!Wenn die Punkte koplanar, bzw. 3 Punkte kollinear sind, so wird false zurueckgeliefert!!!
    
    //!erstmal die Matrix zur Beschreibung der Kugel vorberechnen:
    matrix55<double> sphere_mat(1.,1.,1.,1.,1.,
                                (p0[0]*p0[0])+(p0[1]*p0[1])+(p0[2]*p0[2]), p0[0], p0[1], p0[2], 1.,
                                (p1[0]*p1[0])+(p1[1]*p1[1])+(p1[2]*p1[2]), p1[0], p1[1], p1[2], 1.,
                                (p2[0]*p2[0])+(p2[1]*p2[1])+(p2[2]*p2[2]), p2[0], p2[1], p2[2], 1.,
                                (p3[0]*p3[0])+(p3[1]*p3[1])+(p3[2]*p3[2]), p3[0], p3[1], p3[2], 1.);
    
    double m11 = sphere_mat.sub_det(0,0);
    if (m11 < 0.00000000001 && m11 > -0.00000000001) return false; //!erste Unterdeterminante == 0 
                                                                   //!koplanare, bzw. kollineare Punkte
    
    center[0] = 0.5 * sphere_mat.sub_det(0,1) / m11;
    center[1] = -0.5 * sphere_mat.sub_det(0,2) / m11;
    center[2] = 0.5 * sphere_mat.sub_det(0,3) / m11;
    sradius = (center[0] * center[0]) +
              (center[1] * center[1]) +
              (center[2] * center[2]) -
              (sphere_mat.sub_det(0,4) / m11);
    return true;
}

template<class T, class T2,class T3>
inline bool weighted_points_to_sphere(vec3d<T> const& p0,T2 const& w0,vec3d<T> const& p1,T2 const& w1,vec3d<T> const& p2,T2 const& w2,
                                      vec3d<T> const& p3,T2 const& w3,vec3d<T> &center,T3 &sradius) {
    //!Wie points_to_sphere, aber mit gewichteten Punkten
    //!Handelt es sich nicht um Punkte, sondern Kugeln, so kann das Quadrat des Radius dieser Kugeln
    //!als Gewicht genommen werden. Die resultierende Umkugel ist dann orthogonal zu den 4 Kugeln
    //!erstmal die Matrix zur Beschreibung der Kugel vorberechnen:
    //!Achtung: fuer gewichtete Punkte kann der Quadratradius auch negativ werden!!!
    matrix55<double> sphere_mat(1.,1.,1.,1.,1.,
                                (p0[0]*p0[0])+(p0[1]*p0[1])+(p0[2]*p0[2])-w0, p0[0], p0[1], p0[2], 1.,
                                (p1[0]*p1[0])+(p1[1]*p1[1])+(p1[2]*p1[2])-w1, p1[0], p1[1], p1[2], 1.,
                                (p2[0]*p2[0])+(p2[1]*p2[1])+(p2[2]*p2[2])-w2, p2[0], p2[1], p2[2], 1.,
                                (p3[0]*p3[0])+(p3[1]*p3[1])+(p3[2]*p3[2])-w3, p3[0], p3[1], p3[2], 1.);
    
    double m11 = sphere_mat.sub_det(0,0);
    if (m11 < 0.00000000001 && m11 > -0.00000000001) return false; //!erste Unterdeterminante == 0 
                                                                   //!koplanare, bzw. kollineare Punkte
    
    center[0] = 0.5 * sphere_mat.sub_det(0,1) / m11;
    center[1] = -0.5 * sphere_mat.sub_det(0,2) / m11;
    center[2] = 0.5 * sphere_mat.sub_det(0,3) / m11;
    sradius = (center[0] * center[0]) +
              (center[1] * center[1]) +
              (center[2] * center[2]) -
              (sphere_mat.sub_det(0,4) / m11);
    //! Hier kann auch ein negativer Wert fuer den sradius rauskommen, was aber kein Problem darstellt!
    //! Im Gegenteil -- so wird sichergestellt, dass auch bei einem Wert alpha=0 ein Tetraeder dessen Kugeln
    //! sich alle ueberlappen immer zum alpha-Komplex gehoert.
    return true;
}

template<class T> 
inline double is_P_in_sphere(vec3d<T> const& p0,vec3d<T> const& p1,vec3d<T> const& p2,vec3d<T> const& p3,vec3d<T> const& P) {
    //!Diese Funktion liefert die Determinante deren Vorzeichen entscheidet, ob P innerhalb
    //!der Umkugel um p0 bis p3 liegt, oder ausserhalb.
    //!Ist die Determinante negativ liegt P innerhalb -- ist sie positiv liegt P ausserhalb
    //!-- ist sie Null liegt P auf der Sphaere.
    //!Problem: p0 bis p3 muessen in allgemeiner Lage sein (nicht koplanar)!!!
    matrix55<double> sphere_mat(p0[0],p0[1],p0[2],(p0[0]*p0[0])+(p0[1]*p0[1])+(p0[2]*p0[2]),1.,
                                p1[0],p1[1],p1[2],(p1[0]*p1[0])+(p1[1]*p1[1])+(p1[2]*p1[2]),1.,
                                p2[0],p2[1],p2[2],(p2[0]*p2[0])+(p2[1]*p2[1])+(p2[2]*p2[2]),1.,
                                p3[0],p3[1],p3[2],(p3[0]*p3[0])+(p3[1]*p3[1])+(p3[2]*p3[2]),1.,
                                P[0],P[1],P[2],(P[0]*P[0])+(P[1]*P[1])+(P[2]*P[2]),1.);
    return sphere_mat.det();
    //!Um zu entscheiden werden die Punkte im 4D-Raum auf einen Paraboloid projeziert
    //!=> Die Koordinaten aller Punkte werden um (x^2 + y^2 + z^2) erweitert
    //!=> Liegt P auf dem Paraboloid tiefer als die restlichen Punkte, so bfindet er sich innerhalb.
    //!(Ein anschauliches Beispiel ist die Projektion einer Linie auf eine Parabel)
}

template<class T,class T2> 
inline double is_P_in_powersphere(vec3d<T> const& p0,T2 const& w0,vec3d<T> const& p1,T2 const& w1,
                                  vec3d<T> const& p2,T2 const& w2,vec3d<T> const& p3,T2 const& w3,vec3d<T> const& P,T2 const& wP) {
    //!Wie is_P_in_sphere ,  aber fuer gewichtete Punkte
    //!Problem: p0 bis p3 muessen in allgemeiner Lage sein (nicht koplanar)!!!
    matrix55<double> sphere_mat(p0[0],p0[1],p0[2],(p0[0]*p0[0])+(p0[1]*p0[1])+(p0[2]*p0[2])-w0,1.,
                                p1[0],p1[1],p1[2],(p1[0]*p1[0])+(p1[1]*p1[1])+(p1[2]*p1[2])-w1,1.,
                                p2[0],p2[1],p2[2],(p2[0]*p2[0])+(p2[1]*p2[1])+(p2[2]*p2[2])-w2,1.,
                                p3[0],p3[1],p3[2],(p3[0]*p3[0])+(p3[1]*p3[1])+(p3[2]*p3[2])-w3,1.,
                                P[0],P[1],P[2],(P[0]*P[0])+(P[1]*P[1])+(P[2]*P[2])-wP,1.);
    return sphere_mat.det();
}

template<class T>
inline bool trilaterate(vec3d<T> &sol1,vec3d<T> &sol2,vec3d<T> p1,vec3d<T> p2,vec3d<T> p3,T const& radius) {
    //! frei nach: http://en.wikipedia.org/wiki/Trilateration
    //! angepasst fuer den Spezialfall: r = r1 = r2 = r3
    //! p1,p2,p3 werden NICHT als Referenz uebergeben um erstens die Ruecktransformation
    //! zu sparen und zweitens, um die Koordinaten (durch Rechenungenauigkeit bei vielen
    //! Transformationen) nicht zu veraendern.
    
    //! 1.) Koordinatentransformation:
    vec3d<T> p1_old(p1); p1[0] = 0.; p1[1] = 0.; p1[2] = 0.; p2 -= p1_old; p3 -= p1_old;
    vec3d<T> n(p2); n *= p3; vec3d<T> z(0.,0.,1.);
    vec3d<T> axis(n); axis *= z; axis.norm();
    T alpha = angle(n,z);
    matrix<T> m1 = rotmatrix(axis,alpha); p2 *= m1; p3 *= m1;
    matrix<T> mr1 = rotmatrix(axis,-alpha); vec3d<T> x(1.,0.,0.);
    alpha = angle(p2,x);
    axis = p2; axis *= x; axis.norm();
    matrix<T> m2 = rotmatrix(axis,alpha); p2 *= m2; p3 *= m2;
    matrix<T> mr2 = rotmatrix(axis,-alpha);
    
    //! 2.) x-Wert:
    sol1[0] = p2[0] / 2.; sol2[0] = sol1[0];
    
    //! 3.) y-Wert:
    sol1[1] = (p3[0]*p3[0]+p3[1]*p3[1]) / (2.*p3[1]);
    sol1[1] -= p3[0] * sol1[0] / p3[1];
    sol2[1] = sol1[1];
    
    //! 4.) z-Werte:
    T deci = radius * radius - sol1[0] * sol1[0] - sol1[1] * sol1[1];
    if (deci < 0.) return false;
    sol1[2] = sqrt(deci);
    sol2[2] = -sol1[2];
    sol1 *= mr2; sol1 *= mr1; sol1 += p1_old;
    sol2 *= mr2; sol2 *= mr1; sol2 += p1_old;
    return true;
}


template<class T> 
inline T get_4point_det(vec3d<T> const& p0,vec3d<T> const& p1,vec3d<T> const& p2,vec3d<T> const& p3) {
    matrix44<double> check_det(p0[0],p0[1],p0[2],1,
                               p1[0],p1[1],p1[2],1,
                               p2[0],p2[1],p2[2],1,
                               p3[0],p3[1],p3[2],1);
    return check_det.det();
}

template<class T> 
inline bool is_coplanar(vec3d<T> const& p0,vec3d<T> const& p1,vec3d<T> const& p2,vec3d<T> const& p3,double const& epsilon = 0.1) {
    //!Die Funktion liefert true, wenn die Punkte koplanar sind
    if (fabs(get_4point_det(p0,p1,p2,p3)) < epsilon) return true;
    return false;
}

template<class T>
inline double triple_product(vec3d<T> const& p0,vec3d<T> const& p1,vec3d<T> const& p2,vec3d<T> const& p3) {
    //!Null wenn koplanare Punkte
    //!Positiv, wenn |p0,p1| |p0,p2| |p0,p3| ein Rechtssystem bilden
    //!Negativ fuer Linkssystem
    matrix<double> mat(p1[0]-p0[0],p1[1]-p0[1],p1[2]-p0[2],
                       p2[0]-p0[0],p2[1]-p0[1],p2[2]-p0[2],
                       p3[0]-p0[0],p3[1]-p0[1],p3[2]-p0[2]);
    return mat.det();
}

inline double triple_product(vec3d<double> const& p0,vec3d<double> const& p1,vec3d<double> const& p2,vec3d<float> const& p3) {
    //!Null wenn koplanare Punkte
    //!Positiv, wenn |p0,p1| |p0,p2| |p0,p3| ein Rechtssystem bilden
    //!Negativ fuer Linkssystem
    vec3d<double> pn(p3[0],p3[1],p3[2]);
    matrix<double> mat(p1[0]-p0[0],p1[1]-p0[1],p1[2]-p0[2],
                       p2[0]-p0[0],p2[1]-p0[1],p2[2]-p0[2],
                       pn[0]-p0[0],pn[1]-p0[1],pn[2]-p0[2]);
    return mat.det();
}

inline double triple_product(vec3d<float> const& p0,vec3d<float> const& p1,vec3d<float> const& p2,vec3d<double> const& p3) {
    //!Null wenn koplanare Punkte
    //!Positiv, wenn |p0,p1| |p0,p2| |p0,p3| ein Rechtssystem bilden
    //!Negativ fuer Linkssystem
    vec3d<double> pn(p3[0],p3[1],p3[2]);
    matrix<double> mat(p1[0]-p0[0],p1[1]-p0[1],p1[2]-p0[2],
                       p2[0]-p0[0],p2[1]-p0[1],p2[2]-p0[2],
                       pn[0]-p0[0],pn[1]-p0[1],pn[2]-p0[2]);
    return mat.det();
}

template<class T>
inline bool is_almost_equal(matrix<T> const& m1,matrix<T> const& m2) {
    //! Liefert true, wenn die zwei Matrizen nahezu identisch sind (epsilon = 0.001)
    if (fabs(m1.mat[0][0] - m2.mat[0][0]) > 0.001) return false;
    else if (fabs(m1.mat[0][1] - m2.mat[0][1]) > 0.001) return false;
    else if (fabs(m1.mat[0][2] - m2.mat[0][2]) > 0.001) return false;
    else if (fabs(m1.mat[1][0] - m2.mat[1][0]) > 0.001) return false;
    else if (fabs(m1.mat[1][1] - m2.mat[1][1]) > 0.001) return false;
    else if (fabs(m1.mat[1][2] - m2.mat[1][2]) > 0.001) return false;
    else if (fabs(m1.mat[2][0] - m2.mat[2][0]) > 0.001) return false;
    else if (fabs(m1.mat[2][1] - m2.mat[2][1]) > 0.001) return false;
    else if (fabs(m1.mat[2][2] - m2.mat[2][2]) > 0.001) return false;
    else return true;
}

inline float get_sphere_cap_volume(float const& r,float const& h) {return ((1.0471976 * h * h) * (3. * r - h));}

inline double get_sphere_cap_volume(double const& r,double const& h) {return ((1.0471976 * h * h) * (3. * r - h));}

template<class T>
inline float get_sphere_overlap(float const& r1,float const& r2,vec3d<T> const& m1,vec3d<T> const& m2) {
    float overlap = 0.;
    //! Berechnet das Ueberlappungsvolumen zweier Kugeln mit den Radien r1 und r2,
    //! Sowie den Mittelpunkten m1 und m2
    float rr = r1 + r2;
    float sd = get_square_distance(m1,m2);
    if ((rr*rr) <= sd) return 0.; // keine Ueberlappung
    float d = sqrt(sd);
    
    //! Nun erstmal das Kalotten-Volumen von Kugel 1 berechnen:
    float cos_alpha1 = (r1*r1+sd-r2*r2)/(2.*r1*d);
    float tmp = r1 * cos_alpha1;
    float a1 = r1 - tmp;
    // Fallunterscheidung: Wenn  a1 <= r1  dann ist der Volumenanteil die Kugelkappe
    //                     Wenn  a1 > r1   dann ist der Volumenanteil das Gesamtvolumen minus dem Volumen der Kappe
    if (a1 <= r1) overlap += get_sphere_cap_volume(r1,a1);
    else {
        if (d+r1 <= r2) return (4.1887902 * r1 * r1 * r1); // Kugel 1 liegt komplett innerhalb von Kugel 2
        a1 = r1 + tmp; // tmp ist in diesem Falle negativ (Schnittflaeche liegt hinter dem Mittelpunkt)
        overlap += (4.1887902 * r1 * r1 * r1) - get_sphere_cap_volume(r1,a1);
    }
    
    //! Und jetzt das Volumen der Kalotte von Kugel2 dazu:
    tmp = d - tmp;
    float a2 = r2 - tmp;
    if (a2 <= r2) overlap += get_sphere_cap_volume(r2,a2);
    else {
        if (d+r2 <= r1) return (4.1887902 * r2 * r2 * r2); // Kugel 2 liegt komplett innerhalb von Kugel 1
        a2 = r2 + tmp;
        overlap += (4.1887902 * r2 * r2 * r2) - get_sphere_cap_volume(r2,a2);
    }
    
    return overlap;
}

template<class T>
inline double get_sphere_overlap(double const& r1,double const& r2,vec3d<T> const& m1,vec3d<T> const& m2) {
    double overlap = 0.;
    //! Berechnet das Ueberlappungsvolumen zweier Kugeln mit den Radien r1 und r2,
    //! Sowie den Mittelpunkten m1 und m2
    double rr = r1 + r2;
    double sd = get_square_distance(m1,m2);
    if ((rr*rr) <= sd) return 0.; // keine Ueberlappung
    double d = sqrt(sd);
    
    //! Nun erstmal das Kalotten-Volumen von Kugel 1 berechnen:
    double cos_alpha1 = (r1*r1+sd-r2*r2)/(2.*r1*d);
    double tmp = r1 * cos_alpha1;
    double a1 = r1 - tmp;
    // Fallunterscheidung: Wenn  a1 <= r1  dann ist der Volumenanteil die Kugelkappe
    //                     Wenn  a1 > r1   dann ist der Volumenanteil das Gesamtvolumen minus dem Volumen der Kappe
    if (a1 <= r1) overlap += get_sphere_cap_volume(r1,a1);
    else {
        if (d+r1 <= r2) return (4.1887902 * r1 * r1 * r1); // Kugel 1 liegt komplett innerhalb von Kugel 2
        a1 = r1 + tmp; // tmp ist in diesem Falle negativ (Schnittflaeche liegt hinter dem Mittelpunkt)
        overlap += (4.1887902 * r1 * r1 * r1) - get_sphere_cap_volume(r1,a1);
    }
    
    //! Und jetzt das Volumen der Kalotte von Kugel2 dazu:
    tmp = d - tmp;
    double a2 = r2 - tmp;
    if (a2 <= r2) overlap += get_sphere_cap_volume(r2,a2);
    else {
        if (d+r2 <= r1) return (4.1887902 * r2 * r2 * r2); // Kugel 2 liegt komplett innerhalb von Kugel 1
        a2 = r2 + tmp;
        overlap += (4.1887902 * r2 * r2 * r2) - get_sphere_cap_volume(r2,a2);
    }
    
    return overlap;
}

template<class T>
inline T get_sphere3_overlap(T const& r1,T const& r2,T const& r3,
                             vec3d<T> const& m1,vec3d<T> const& m2,vec3d<T> const& m3) {
    T overlap = 0.;
    //! Berechnet das Ueberlappungsvolumen dreier Kugeln mit den Radien r1,r2 und r3,
    //! sowie den Mittelpunkten m1,m2 und m3:
    //! nach: "Volume of the Intersection of Three Spheres of Unequal Size"
    //!       von Gibson und Scheraga
    T a2 = get_square_distance(m2,m3); T a = sqrt(a2);
    T b2 = get_square_distance(m1,m3); T b = sqrt(b2);
    T c2 = get_square_distance(m1,m2); T c = sqrt(c2);
    
    T r12 = r1*r1; T r22 = r2*r2; T r32 = r3*r3;
    
    T e1 = (r22-r32) / a2;
    T e2 = (r32-r12) / b2;
    T e3 = (r12-r22) / c2;
    
    T w2 = (r12*a2+r22*b2+r32*c2) * (a2+b2+c2) - 2.*(r12*a2*a2+r22*b2*b2+r32*c2*c2) + a2*b2*c2*(e1*e2+e2*e3+e3*e1-1.);
    
    //! Wenn w2 == 0, dann beruehren sich die Kugeln in einem Punkt => Volumen gleich Null
    //! Wenn w2 < 0, dann gibt es keinen gemeinsamen Punkt, jedoch kann es sehr wohl ein Schnittvolumen geben:
    //!             1.) Wenn eine Kugel C den Schnitt AB der 2 anderen komplett beinhaltet => V_ABC = V_AB
    //!             2.) Wenn C komplett oder teilweise in AB enthalten ist => V_ABS = V_AB + V_AC - V_A
    //! Vorerst auch fuer diese Sonderfaelle Null zurueckgeben. Spaeter unbedingt anpassen!!!
    if (w2 <= 0.00001) return 0.;
    
    T w = sqrt(w2);
    
    T q1 = a * (b2+c2-a2+r22+r32-2.*r12+e1*(b2-c2));
    T q2 = b * (c2+a2-b2+r32+r12-2.*r22+e2*(c2-a2));
    T q3 = c * (a2+b2-c2+r12+r22-2.*r32+e3*(a2-b2));
    
    T c_PI = acos(-1.);
    T w_2 = w*2.;
    
    T t_1 = atan(w_2 / q1);
    if (t_1 < 0.) t_1 += c_PI;
    else if (t_1 > c_PI) t_1 -= c_PI;
    
    T t_2 = atan(w_2 / q2);
    if (t_2 < 0.) t_2 += c_PI;
    else if (t_2 > c_PI) t_2 -= c_PI;
    
    T t_3 = atan(w_2 / q3);
    if (t_3 < 0.) t_3 += c_PI;
    else if (t_3 > c_PI) t_3 -= c_PI;
    
    T t_4 = atan(b*w*(1.-e2)/(r1*q2));
    if (t_4 < 0.) t_4 += c_PI;
    else if (t_4 > c_PI) t_4 -= c_PI;
    
    T t_5 = atan(c*w*(1.+e3)/(r1*q3));
    if (t_5 < 0.) t_5 += c_PI;
    else if (t_5 > c_PI) t_5 -= c_PI;
    
    T t_6 = atan(c*w*(1.-e3)/(r2*q3));
    if (t_6 < 0.) t_6 += c_PI;
    else if (t_6 > c_PI) t_6 -= c_PI;
    
    T t_7 = atan(a*w*(1.+e1)/(r2*q1));
    if (t_7 < 0.) t_7 += c_PI;
    else if (t_7 > c_PI) t_7 -= c_PI;
    
    T t_8 = atan(a*w*(1.-e1)/(r3*q1));
    if (t_8 < 0.) t_8 += c_PI;
    else if (t_8 > c_PI) t_8 -= c_PI;
    
    T t_9 = atan(b*w*(1.+e2)/(r3*q2));
    if (t_9 < 0.) t_9 += c_PI;
    else if (t_9 > c_PI) t_9 -= c_PI;
    
    overlap += w/6.;
    overlap -= 0.5 * a * (r22+r32-a2*((1./6.)-(0.5*e1*e1))) * t_1;
    overlap -= 0.5 * b * (r32+r12-b2*((1./6.)-(0.5*e2*e2))) * t_2;
    overlap -= 0.5 * c * (r12+r22-c2*((1./6.)-(0.5*e3*e3))) * t_3;
    overlap += (2.*r1*r12/3.) * (t_4+t_5);
    overlap += (2.*r2*r22/3.) * (t_6+t_7);
    overlap += (2.*r3*r32/3.) * (t_8+t_9);
    
    return overlap;
}

template<class T>
inline T get_sphere4_overlap(T const& r1,T const& r2,T const& r3,T const& r4,
                             vec3d<T> const& m1,vec3d<T> const& m2,vec3d<T> const& m3,vec3d<T> const& m4) {
    T overlap = 0.;
    //! Berechnet das Ueberlappungsvolumen von 4 Kugeln mit den Radien r1,r2,r3 und r4,
    //! sowie den Mittelpunkten m1,m2,m3 und m4:
    //! nach: "Exact calculation of the volume and surface area of fused hard-sphere molecules with unequal radii"
    //!       von Gibson und Scheraga
    
    //! Wenn ein 3er Overlap Null ist, dann gibt es keinen 4er Overlap:
    T s31 = get_sphere3_overlap(r1,r2,r3,m1,m2,m3);
    if (s31 < 0.00001) return 0.;
    T s32 = get_sphere3_overlap(r1,r2,r4,m1,m2,m4);
    if (s32 < 0.00001) return 0.;
    T s33 = get_sphere3_overlap(r1,r3,r4,m1,m3,m4);
    if (s33 < 0.00001) return 0.;
    T s34 = get_sphere3_overlap(r2,r3,r4,m2,m3,m4);
    if (s34 < 0.00001) return 0.;
    
    
    T a2 = get_square_distance(m2,m3); T a = sqrt(a2);
    T b2 = get_square_distance(m1,m3); T b = sqrt(b2);
    T c2 = get_square_distance(m1,m2); T c = sqrt(c2);
    T g2 = get_square_distance(m2,m4); T g = sqrt(g2);
    T h2 = get_square_distance(m3,m4); T h = sqrt(h2);
    T f2 = get_square_distance(m1,m4); T f = sqrt(f2);
    
    T r12 = r1*r1; T r22 = r2*r2; T r32 = r3*r3; T r42 = r4*r4;
    
    matrix55<T> mw(0.,c2,b2,f2,1.,
                   c2,0.,a2,g2,1.,
                   b2,a2,0.,h2,1.,
                   f2,g2,h2,0.,1.,
                   1.,1.,1.,1.,0.);
    T w2 = 0.5 * mw.det();
    T w = sqrt(w2); //!check
    
    
    T s1 = a * (b2+c2-a2+g2+h2-2.*f2+(g2-h2)*(b2-c2)/a2);
    T s2 = b * (c2+a2-b2+h2+f2-2.*g2+(h2-f2)*(c2-a2)/b2);
    T s3 = c * (a2+b2-c2+f2+g2-2.*h2+(f2-g2)*(a2-b2)/c2);
    T s4 = h * (f2+b2-h2+a2+g2-2.*c2+(f2-b2)*(a2-g2)/h2);
    T s5 = g * (h2+a2-g2+c2+f2-2.*b2+(c2-f2)*(h2-a2)/g2);
    T s6 = f * (g2+c2-f2+b2+h2-2.*a2+(b2-h2)*(g2-c2)/f2);
    
    
    T c_PI = acos(-1.);
    T w_2 = w*2.;
    
    T t_1 = atan(w_2/s1);
    if (t_1 < 0.) t_1 += c_PI;
    else if (t_1 > c_PI) t_1 -= c_PI;
    t_1 *= ( (r22+r32)*a+(r22-r32)*(r22-r32)/(2.*a)-a2*a/6. );
    
    T t_2 = atan(w_2/s2);
    if (t_2 < 0.) t_2 += c_PI;
    else if (t_2 > c_PI) t_2 -= c_PI;
    t_2 *= ( (r12+r32)*b+(r12-r32)*(r12-r32)/(2.*b)-b2*b/6. );
    
    T t_3 = atan(w_2/s3);
    if (t_3 < 0.) t_3 += c_PI;
    else if (t_3 > c_PI) t_3 -= c_PI;
    t_3 *= ( (r12+r22)*c+(r12-r22)*(r12-r22)/(2.*c)-c2*c/6. );
    
    T t_4 = atan(w_2/s4);
    if (t_4 < 0.) t_4 += c_PI;
    else if (t_4 > c_PI) t_4 -= c_PI;
    t_4 *= ( (r32+r42)*h+(r32-r42)*(r32-r42)/(2.*h)-h2*h/6. );
    
    T t_5 = atan(w_2/s5);
    if (t_5 < 0.) t_5 += c_PI;
    else if (t_5 > c_PI) t_5 -= c_PI;
    t_5 *= ( (r22+r42)*g+(r22-r42)*(r22-r42)/(2.*g)-g2*g/6. );
    
    T t_6 = atan(w_2/s6);
    if (t_6 < 0.) t_6 += c_PI;
    else if (t_6 > c_PI) t_6 -= c_PI;
    t_6 *= ( (r12+r42)*f+(r12-r42)*(r12-r42)/(2.*f)-f2*f/6. );
    
    
    overlap -= w / 12.;
    overlap -= c_PI * (r12*r1+r22*r2+r32*r3+r42*r4) / 3.;
    overlap += 0.25 * t_1;
    overlap += 0.25 * t_2;
    overlap += 0.25 * t_3;
    overlap += 0.25 * t_4;
    overlap += 0.25 * t_5;
    overlap += 0.25 * t_6;
    overlap += 0.5 * s31;
    overlap += 0.5 * s32;
    overlap += 0.5 * s33;
    overlap += 0.5 * s34;
    
    if (overlap < 0.) return 0.;
    else return overlap;
}

template<class T>
inline T get_triangle_area(vec3d<T> const& p1,vec3d<T> const& p2,vec3d<T> const& p3) {
    //! Dreiecksflaeche nach Heron:
    vec3d<T> a(p2); a -= p1;
    vec3d<T> b(p3); b -= p2;
    vec3d<T> c(p1); b -= p3;
    T av = a.value(); T bv = b.value(); T cv = c.value();
    T s = av; s += bv; s+= cv; s /= 2.;
    return (sqrt(s * (s-av) * (s-bv) * (s-cv)));
}

template<class T>
inline T get_spherical_triangle_area(vec3d<T> const& center,T const& radius,vec3d<T> const& p1,
                                     vec3d<T> const& p2,vec3d<T> const& p3) {
    //! Diese Funktion berechnet die Flaeche eines sphaerischen Dreiecks, anhand von Kugelmittelpunkt,
    //! Radius und 3er Punkte. Die 3 Punkte muessen NICHT auf der Kugeloberflaeche liegen!!! Die
    //! Schnittpunkte der Geraden "Pn--Kugelmittelpunkt" mit der Kugeloberflaeche definieren das Dreieck.
    T a = angle(p1,center,p2);
    T b = angle(p1,center,p3);
    T c = angle(p2,center,p3);
    T alpha = a; if (fabs(sin(b)*sin(c)) > 0.00001) alpha = acos((cos(a)-cos(b)*cos(c))/(sin(b)*sin(c)));
    T beta = b; if (fabs(sin(a)*sin(c)) > 0.00001) beta = acos((cos(b)-cos(a)*cos(c))/(sin(a)*sin(c)));
    T gamma = c; if (fabs(sin(a)*sin(b)) > 0.00001) gamma = acos((cos(c)-cos(a)*cos(b))/(sin(a)*sin(b)));
    return((alpha + beta + gamma - acos(-1.)) * radius * radius);
}

template<class T>
inline T get_sphere_tetrahedron_intersection(vec3d<T> const& center,T const& radius,vec3d<T> const& p1,
                                             vec3d<T> const& p2,vec3d<T> const& p3) {
    //! Diese Funktion berechnet das Schnittvolumen einer Kugel mit einem Tetraeder, wobei ein Punkt
    //! des Tetraeders im Zentrum der Kugel liegt und die restlichen Punkte des Tetraeders ausserhalb
    //! der Kugel liegen muessen.
    /*
    T pi4 = 4. * acos(-1.) * radius * radius;
    T mod = get_spherical_triangle_area(center,radius,p1,p2,p3) / pi4;
    return (mod * radius * pi4 / 3.);
    */
    //! Leider weiss ich nicht mehr wirklich, wo ich obige Formel her habe, aber
    //! 1.) kommt exakt nur die Haelfte des korrekten Ergebnisses raus (ok, nun koennt man einfach
    //!      noch mit 2 multiplizieren, ist aber doof wenn ich die Grundlage nicht mehr kenne)
    //! 2.) geht es viel einfacher ;)
    //! einfach den Raumwinkel Omega bestimmen, den das Tetraeder aufspannt und dann das Kugelvolumen
    //! mit Omega/2PI multiplizieren:
    //! (Ich rufe mal nicht immer extra dihedral auf, um die Rechnung etwas effizienter zu machen)
    T vol = 4. * acos(-1.) * radius * radius * radius / 3.;
    vec3d<T> v1(center); v1 -= p1;
    vec3d<T> v2(center); v2 -= p2;
    vec3d<T> v3(center); v3 -= p3;
    vec3d<T> n1(v2); n1 *= v1; T val1 = n1.value(); // nach aussen
    v2 *= v3; T val2 = v2.value(); // nach innen
    vec3d<T> n2(v2); n2 *= -1.;
    v3 *= v1; T val3 = v3.value(); // nach innen
    return (((acos(skalarproduct(n1,v2)/(val1*val2))
              +acos(skalarproduct(n2,v3)/(val2*val3))
              +acos(skalarproduct(v3,n1)/(val3*val1))
              -acos(-1.))/(2.*acos(-1.))) * vol);
}

template<class T>
inline T get_sphere_tetrahedron_intersection(vec3d<T> const& center,vec3d<T> const& p1,
                                             vec3d<T> const& p2,vec3d<T> const& p3,T const& vol) {
    //! Variante zur vorhergehenden Funktion, bei der Volumen statt Radius mitgegeben wird:
    vec3d<T> v1(center); v1 -= p1;
    vec3d<T> v2(center); v2 -= p2;
    vec3d<T> v3(center); v3 -= p3;
    vec3d<T> n1(v2); n1 *= v1; T val1 = n1.value(); // nach aussen
    v2 *= v3; T val2 = v2.value(); // nach innen
    vec3d<T> n2(v2); n2 *= -1.;
    v3 *= v1; T val3 = v3.value(); // nach innen
    return (((acos(skalarproduct(n1,v2)/(val1*val2))
              +acos(skalarproduct(n2,v3)/(val2*val3))
              +acos(skalarproduct(v3,n1)/(val3*val1))
              -acos(-1.))/(2.*acos(-1.))) * vol);
}

//!============================================================================================================
//! Und nochmal einige spezielle Funktionen fuer Volumen und Flaechenberechnungen auf
//! Basis der Triangulierung. Hier gibt es einiges an Redundanz zu bereits vorhandenen
//! Funktionen, aber die Abtrennung schafft etwas Uebersicht.
//! Die Berechnungen basieren auf:
//! "Measuring Space Filling Diagrams and Voids" von Edelsbrunner und Fu
//! ACHTUNG:  Alle Winkel werden hier auf 0 bis 1 normiert. (Revolution angles)

const double DL_PI = acos(-1.);

template<class T> 
inline const T angle_dihedral(const vec3d<T> &s,const vec3d<T> &t,const vec3d<T> &u,const vec3d<T> &v) {
    //! liefert einen Winkel zwischen 0 und 1
    vec3d<T> us(u); us -= s;
    vec3d<T> ut(u); ut -= t;
    us *= ut; us /= sqrt(skalarproduct(us,us));
    vec3d<T> vs(v); vs -= s;
    vec3d<T> vt(v); vt -= t;
    vs *= vt; vs /= sqrt(skalarproduct(vs,vs));
    return (acos(skalarproduct(us,vs))/DL_PI);
}

template<class T> //! wird nur von sector_vol genutzt und ist somit eigentlich obsolet!!!
inline const T angle_solid(const vec3d<T> &i,const vec3d<T> &j,const vec3d<T> &k,const vec3d<T> &l) {
    T pj = angle_dihedral(i,j,l,k);
    T pk = angle_dihedral(i,k,j,l);
    T pl = angle_dihedral(i,l,k,j);
    
//    return (((pj+pk+pl)/2.)-0.25); ///FALSCH bei Edelsbrunner?
    //! Hmmm, nein -- eigentlich sollte ein voller solid angle 2Pi sein
    //! Beim Tetraeder gilt: solid_angle = da1 + da2 + da3 - PI
    //! Wenn da in Revolution ist => solid_angle = da1*PI + da2*PI + da3*PI - PI = PI * (da1 + da2 + da3 - 1)
    //! => Revolution: solid_angle /= 2PI  =>  solid_angle = 0.5 * (da1 + da2 + da3 - 1)
    return (0.5 * (pj+pk+pl-1.));
}

template<class T> //! Lahm:  lieber get_sphere_tetrahedron_intersection benutzen!!!
inline const T sector_vol(const vec3d<T> &i,const vec3d<T> &j,const vec3d<T> &k,const vec3d<T> &l,const T &ir) {
    return (angle_solid(i,j,k,l) * ir*ir*ir * 4. * DL_PI / 3.); //! ir = Radius der Kugel i
}

//!----------------------------------
template<class T> 
inline const T disk_length(const vec3d<T> &i,const vec3d<T> &j,const T &ir,const T &jr) {
    return (2. * DL_PI * disk_radius(i,j,ir,jr));
}

template<class T> 
inline const T disk_radius(const vec3d<T> &i,const vec3d<T> &j,const T &ir,const T &jr) {
    return sqrt(cap_height(i,j,ir,jr) * (2.*ir - cap_height(i,j,ir,jr)));
}

template<class T> 
inline const T disk_area(const vec3d<T> &i,const vec3d<T> &j,const T &ir,const T &jr) {
    return (disk_radius(i,j,ir,jr) * disk_length(i,j,ir,jr) / 2.);
}

template<class T> 
inline const vec3d<T> center2(const vec3d<T> &i,const vec3d<T> &j,const T &ir,const T &jr) {
    //! orthogonales Zentrum von 2 Kugeln bestimmen:
    vec3d<T> st(i); st -= j;
    T aux = 2. * skalarproduct(st,st);
    T lamda = 0.5 - (((ir * ir)-(jr * jr)) / aux);
    st = i; st *= lamda;
    vec3d<T> vj(j); vj *= (1. - lamda);
    st += vj;
    return st;
}

template<class T> //! Die Funktion sollte jetzt richtige Ergebnisse liefern:
inline const T cap_height(const vec3d<T> &i,const vec3d<T> &j,const T &ir,const T &jr) {
    vec3d<T> c2(center2(i,j,ir,jr));
    if (get_distance(i,j) < get_distance(j,c2)) return (ir + get_distance(i,c2));
    else return (ir - get_distance(i,c2));
}

template<class T> //! obsolet, weil nur von cap_vol verwendet
inline const T cap_area(const vec3d<T> &i,const vec3d<T> &j,const T &ir,const T &jr) {
    return (2. * DL_PI * ir * cap_height(i,j,ir,jr));
}

template<class T> //! obsolet, weil nur von ball2_vol verwendet
inline const T cap_vol(const vec3d<T> &i,const vec3d<T> &j,const T &ir,const T &jr) {
    T s = ir * cap_area(i,j,ir,jr) / 3.; //! ir = Radius von i
    T c = (ir - cap_height(i,j,ir,jr)) * disk_area(i,j,ir,jr) / 3.;
    return (s-c);
}

template<class T> //! obsolet, weil nur von wedge_vol verwendet -> lieber get_sphere_overlap() nutzen (siehe wedge_vol)
inline const T ball2_vol(const vec3d<T> &i,const vec3d<T> &j,const T &ir,const T &jr) {
    return (cap_vol(i,j,ir,jr) + cap_vol(j,i,jr,ir));
}

template<class T> //! lieber get_sphere_overlap() benutzen!!!
inline const T wedge_vol(const vec3d<T> &i,const vec3d<T> &j,const vec3d<T> &k,const vec3d<T> &l,const T &ir,const T &jr) {
    if (get_distance(i,j) >= ir+jr) return 0.;
    else return (angle_dihedral(i,j,k,l) * ball2_vol(i,j,ir,jr)); //! ir,jr = vdW-Radien der Kugeln
}
//!----------------------------------

template<class T> 
inline const vec3d<T> center3(const vec3d<T> &i,const vec3d<T> &j,const vec3d<T> &k,const T &ir,const T &jr,const T &kr) {
    matrix<T> ma1(i[1],i[2],1.,
                  j[1],j[2],1.,
                  k[1],k[2],1.);
    T a1 = ma1.det();
    matrix<T> ma2(i[2],i[0],1.,
                  j[2],j[0],1.,
                  k[2],k[0],1.);
    T a2 = ma2.det();
    matrix<T> ma3(i[0],i[1],1.,
                  j[0],j[1],1.,
                  k[0],k[1],1.);
    T a3 = ma3.det();
    matrix<T> ma4(i[0],i[1],i[2],
                  j[0],j[1],j[2],
                  k[0],k[1],k[2]);
    T a4 = ma4.det();
    
    T i0 = 0.5 * (ir*ir - i[0]*i[0] - i[1]*i[1] - i[2]*i[2]);
    T j0 = 0.5 * (jr*jr - j[0]*j[0] - j[1]*j[1] - j[2]*j[2]);
    T k0 = 0.5 * (kr*kr - k[0]*k[0] - k[1]*k[1] - k[2]*k[2]);
    
    matrix44<T> md0(i[0],i[1],i[2],1.,
                    j[0],j[1],j[2],1.,
                    k[0],k[1],k[2],1.,
                    a1,a2,a3,0.);
    T d0 = md0.det();
    matrix44<T> mdx(-i0,i[1],i[2],1.,
                    -j0,j[1],j[2],1.,
                    -k0,k[1],k[2],1.,
                    a4,a2,a3,0.);
    T dx = mdx.det();
    matrix44<T> mdy(i[0],-i0,i[2],1.,
                    j[0],-j0,j[2],1.,
                    k[0],-k0,k[2],1.,
                    a1,a4,a3,0.);
    T dy = mdy.det();
    matrix44<T> mdz(i[0],i[1],-i0,1.,
                    j[0],j[1],-j0,1.,
                    k[0],k[1],-k0,1.,
                    a1,a2,a4,0.);
    T dz = mdz.det();
    return vec3d<T>(dx/d0,dy/d0,dz/d0);
}

template<class T> 
inline const T segment_height(const vec3d<T> &i,const vec3d<T> &j,const vec3d<T> &k,const T &ir,const T &jr,const T &kr) {
    vec3d<T> y3 = center3(i,j,k,ir,jr,kr);
    vec3d<T> y2 = center2(i,j,ir,jr);
    if (get_distance(k,y2) < kr) return (disk_radius(i,j,ir,jr) + get_distance(y2,y3));
    else return (disk_radius(i,j,ir,jr) - get_distance(y2,y3));
}

template<class T> 
inline const T segment_length(const vec3d<T> &i,const vec3d<T> &j,const vec3d<T> &k,const T &ir,const T &jr,const T &kr) {
    return (segment_angle(i,j,k,ir,jr,kr) * disk_length(i,j,ir,jr));
}

template<class T> //! Muesste hier nicht das Ergebnis mit 0.5 multipliziert werden??? (weil dihedral nur auf PI normiert)
inline const T segment_angle(const vec3d<T> &i,const vec3d<T> &j,const vec3d<T> &k,const T &ir,const T &jr,const T &kr) {
    vec3d<T> jk = triangle_dual(i,j,k,ir,jr,kr);
    vec3d<T> kj = triangle_dual(i,k,j,ir,kr,jr);
    return (angle_dihedral(i,j,k,jk) + angle_dihedral(i,j,k,kj));
///    return (0.5*(angle_dihedral(i,j,k,jk) + angle_dihedral(i,j,k,kj)));
}

template<class T> 
inline const vec3d<T> triangle_dual(const vec3d<T> &i,const vec3d<T> &j,const vec3d<T> &k,const T &ir,const T &jr,const T &kr) {
    vec3d<T> y = center3(i,j,k,ir,jr,kr); vec3d<T> yy(y); y -= i;
    vec3d<T> n(j); n -= i; vec3d<T> n2(k); n2 -= i;
    n *= n2;
    T s1 = skalarproduct(y,n);
    T s2 = skalarproduct(n,n);
    T s3 = skalarproduct(y,y);
    T e = (sqrt(s1*s1-s3*s2+ir*ir*s2)-s1) / s2;
    n *= e; yy += n;
    return yy;
}

template<class T> 
inline const T segment_area(const vec3d<T> &i,const vec3d<T> &j,const vec3d<T> &k,const T &ir,const T &jr,const T &kr) {
    T s = 0.5 * disk_radius(i,j,ir,jr) * segment_length(i,j,k,ir,jr,kr);
    vec3d<T> jk = triangle_dual(i,j,k,ir,jr,kr);
    vec3d<T> kj = triangle_dual(i,k,j,ir,kr,jr);
    T h = disk_radius(i,j,ir,jr) - segment_height(i,j,k,ir,jr,kr);
    T t = 0.5 * h * get_distance(jk,kj);
    return (s-t);
}

template<class T> 
inline const T cap2_area(const vec3d<T> &i,const vec3d<T> &j,const vec3d<T> &k,const T &ir,const T &jr,const T &kr) {
    vec3d<T> pjk = triangle_dual(i,j,k,ir,jr,kr);
    vec3d<T> pkj = triangle_dual(i,k,j,ir,kr,jr);
    T lj = segment_angle(i,j,k,ir,jr,kr);
    T lk = segment_angle(i,k,j,ir,kr,jr);
    T jk = 0.5 - angle_dihedral(i,pjk,j,k);
    T kj = 0.5 - angle_dihedral(i,pkj,k,j);
    T a1 = 2.* DL_PI * ir*ir * (jk+kj);
    T a2 = 2. * DL_PI * ir * lj * (ir-cap_height(i,j,ir,jr));
    T a3 = 2. * DL_PI * ir * lk * (ir-cap_height(i,k,ir,kr));
    return (a1 - a2 - a3);
}

template<class T> 
inline const T cap2_vol(const vec3d<T> &i,const vec3d<T> &j,const vec3d<T> &k,const T &ir,const T &jr,const T &kr) {
    T s = ir * cap2_area(i,j,k,ir,jr,kr) / 3.;
    T cj = (ir - cap_height(i,j,ir,jr)) * segment_area(i,j,k,ir,jr,kr) / 3.;
    T ck = (ir - cap_height(i,k,ir,kr)) * segment_area(i,k,j,ir,kr,jr) / 3.;
    return (s-cj-ck);
}

template<class T> //! kommt Unfug raus :(  (lieber get_sphere3_overlap() nutzen)
inline const T ball3_vol(const vec3d<T> &i,const vec3d<T> &j,const vec3d<T> &k,const T &ir,const T &jr,const T &kr) {
    return (cap2_vol(i,j,k,ir,jr,kr) + cap2_vol(j,k,i,jr,kr,ir) + cap2_vol(k,i,j,kr,ir,jr));
}

template<class T> //! Leider kommt hier Schwachsinn raus
                  //! zumindest center3 liefert korrekte Ergebnisse
                  //! lieber 0.5 * get_sphere3_overlap()  nutzen!!!
inline const T pawn_vol(const vec3d<T> &i,const vec3d<T> &j,const vec3d<T> &k,const T &ir,const T &jr,const T &kr) {
    vec3d<T> y = center3(i,j,k,ir,jr,kr);
    if (get_distance(i,y) > ir ||
        get_distance(j,y) > jr ||
        get_distance(k,y) > kr) return 0.;
    else return (ball3_vol(i,j,k,ir,jr,kr) / 2.);
}

//!============================================================================================================

//--------------------------------------------------------------------------
//globale Operatoren:
template<class T> 
const vec3d<T> operator+(vec3d<T> const& links,vec3d<T> const& rechts) {
    return (vec3d<T>(links) += rechts);
}

template<class T> 
const vec3d<T> operator-(vec3d<T> const& links,vec3d<T> const& rechts) {
    return (vec3d<T>(links) -= rechts);
}

template<class T> 
const vec3d<T> operator*(vec3d<T> const& links,vec3d<T> const& rechts) { // * fr Vektorprodukt berladen
    return (vec3d<T>(links) *= rechts);
}

template<class T> 
const vec3d<T> operator*(vec3d<T> const& links,double const& rechts) { // * fr multipl. mit Skalar berladen
    return (vec3d<T>(links) *= rechts);
}

template<class T> 
const vec3d<T> operator*(double const& links,vec3d<T> const& rechts) {
    return (vec3d<T>(rechts) *= links);
}

template<class T> 
const vec3d<T> operator*(matrix<T> const& links,vec3d<T> const& rechts) { // * fr transform berladen
    return (vec3d<T>(rechts)*=links);
}

template<class T> 
const vec3d<T> operator*(vec3d<T> const& links,matrix<T> const& rechts) { // * fr transform berladen
    return (vec3d<T>(links)*=rechts);
}

template<class T> 
const vec3d<T> operator/(vec3d<T> const& links,double const& rechts) { // * fr Division mit Skalar berladen
    return (vec3d<T>(links) /= rechts);
}

template<class T> 
ostream &operator<<(ostream &os,vec3d<T> const& vec) { //<< ueberladen zur Ausgabe von vec3d-Objekten
    os << "[ " << vec.vector[0] << " , " << vec.vector[1] << " , " << vec.vector[2] << " ]";
    return os;
}

template<class T> 
ostream &operator<<(ostream &os,matrix<T> const& m) {
    os << endl << "[  ";
    os.width(6); os.setf(ios::fixed); os.precision(6);
    os << right << m.mat[0][0] << "  ";
    os.width(6); os << right << m.mat[0][1] << "  ";
    os.width(6); os << right << m.mat[0][2] << "  ]" << "\n" << "[  ";
    os.width(6); os << right << m.mat[1][0] << "  ";
    os.width(6); os << right << m.mat[1][1] << "  ";
    os.width(6); os << right << m.mat[1][2] << "  ]" << "\n" << "[  ";
    os.width(6); os << right << m.mat[2][0] << "  ";
    os.width(6); os << right << m.mat[2][1] << "  ";
    os.width(6); os << right << m.mat[2][2] << "  ]";
    os.unsetf(ios::fixed);
    return os;
}

template<class T> 
ostream &operator<<(ostream &os,matrix44<T> const& m) {
    os << endl << "[  ";
    os.width(6); os.setf(ios::fixed); os.precision(6);
    os << right << m.mat[0][0] << "  ";
    os.width(6); os << right << m.mat[0][1] << "  ";
    os.width(6); os << right << m.mat[0][2] << "  ";
    os.width(6); os << right << m.mat[0][3] << "  ]" << "\n" << "[  ";
    os.width(6); os << right << m.mat[1][0] << "  ";
    os.width(6); os << right << m.mat[1][1] << "  ";
    os.width(6); os << right << m.mat[1][2] << "  ";
    os.width(6); os << right << m.mat[1][3] << "  ]" << "\n" << "[  ";
    os.width(6); os << right << m.mat[2][0] << "  ";
    os.width(6); os << right << m.mat[2][1] << "  ";
    os.width(6); os << right << m.mat[2][2] << "  ";
    os.width(6); os << right << m.mat[2][3] << "  ]" << "\n" << "[  ";
    os.width(6); os << right << m.mat[3][0] << "  ";
    os.width(6); os << right << m.mat[3][1] << "  ";
    os.width(6); os << right << m.mat[3][2] << "  ";
    os.width(6); os << right << m.mat[3][3] << "  ]";
    os.unsetf(ios::fixed);
    return os;
}

#endif
