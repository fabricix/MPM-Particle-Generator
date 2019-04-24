/*
 * tensor.h
 *
 *  Created on: Oct 6, 2018
 *      Author: fabricio
 */


#ifndef INC_TENSOR_H
#define INC_TENSOR_H

#include <iostream>
#include <iomanip>
#include <cassert>
#include <limits>
#include <cmath>
#include <math.h>

const double pi=3.1415926535897932384626433832795028841971694;
const double MACHINE_TOL=std::numeric_limits<double>::epsilon();

/* ---------------------------------------------------------------------------*/

//template<class charT,class traits>
//   basic_ostream<charT,traits>&tab(basic_ostream<charT,traits>&os){os<<'\t'<<setw(10);return os;}
//template<class charT,class traits>
//   basic_ostream<charT,traits>&nwl(basic_ostream<charT,traits>&os){os<<'\n'<<setw(10);return os;}

/* ---------------------------------------------------------------------------*/

inline double Pceil(const double&d){return ceil(d-MACHINE_TOL);}

template<typename T>inline T abs(T a){return(a<0.?-a:a);}
template<typename T>inline T sgn(T a){return(a<0.?-1.:1.);}
template<bool>struct Assert;
template<>struct Assert<true>{}; // Assert<(1==2)>();

//template<typename T>T max(const vector<T>&v){double m=v[0];for(unsigned i=0;i<v.size();i+=1)m=(v[i]>m?v[i]:m);return m;}
//template<typename T>T min(const vector<T>&v){double m=v[0];for(unsigned i=0;i<v.size();i+=1)m=(v[i]<m?v[i]:m);return m;}

struct Matrix2;
struct Matrix3;

//-------------------------------------------------------------------------//
// Vector 2D
//-------------------------------------------------------------------------//

struct Vector2{

   double x,y;

   Vector2(){x=0.;y=0.;}
   Vector2(double a,double b){x=a;y=b;}
   Vector2(const Vector2&v){x=v.x;y=v.y;}

   Vector2&operator=(const double&a){x=a;y=a;return*this;}
   Vector2&operator=(const Vector2&v){x=v.x;y=v.y;return*this;}
   Vector2&operator+=(const Vector2&v){x+=v.x;y+=v.y;return*this;}
   Vector2&operator-=(const Vector2&v){x-=v.x;y-=v.y;return*this;}
   Vector2&operator*=(const Vector2&v){x*=v.x;y*=v.y;return*this;}
   Vector2&operator/=(const Vector2&v){x/=v.x;y/=v.y;return*this;}
   Vector2&operator+=(const double&v){x+=v;y+=v;return*this;}
   Vector2&operator-=(const double&v){x-=v;y-=v;return*this;}
   Vector2&operator*=(const double&v){x*=v;y*=v;return*this;}
   Vector2&operator/=(const double&v){x/=v;y/=v;return*this;}

   const Vector2 operator-(               )const{return Vector2(-x,-y);}
   const Vector2 operator+(const Vector2&v)const{return Vector2(x+v.x,y+v.y);}
   const Vector2 operator-(const Vector2&v)const{return Vector2(x-v.x,y-v.y);}
   const Vector2 operator*(const Vector2&v)const{return Vector2(x*v.x,y*v.y);}
   const Vector2 operator/(const Vector2&v)const{return Vector2(x/v.x,y/v.y);}
   const Vector2 operator+(const double&a)const{return Vector2(x+a,y+a);}
   const Vector2 operator-(const double&a)const{return Vector2(x-a,y-a);}
   const Vector2 operator*(const double&a)const{return Vector2(x*a,y*a);}
   const Vector2 operator/(const double&a)const{return Vector2(x/a,y/a);}

   bool operator==(const Vector2&v){return x==v.x && y==v.y;}
   bool operator==(const double &v){return x==v  && y==v  ;}
   bool operator!=(const Vector2&v){return x!=v.x || y!=v.y;}
   bool operator!=(const double &v){return x!=v  || y!=v  ;}

   inline double  inner(const Vector2&v)const;
   inline Vector2 inner(const Matrix2&m)const;
   inline Matrix2 outer(const Vector2&v)const;
};

inline const Vector2 operator+(const double&a,const Vector2&v){return Vector2(a+v.x,a+v.y);}
inline const Vector2 operator-(const double&a,const Vector2&v){return Vector2(a-v.x,a-v.y);}
inline const Vector2 operator*(const double&a,const Vector2&v){return Vector2(a*v.x,a*v.y);}
inline const Vector2 operator/(const double&a,const Vector2&v){return Vector2(a/v.x,a/v.y);}

//-------------------------------------------------------------------------//
// Vector 3D
//-------------------------------------------------------------------------//

struct Vector3{

   double x,y,z;

   Vector3(){x=0.0;y=0.0;z=0.0;}
   Vector3(double f){x=f;y=f;z=f;}
   Vector3(double vx,double vy, double vz){x=vx;y=vy;z=vz;}
   Vector3(const Vector3&v){x=v.x;y=v.y; z=v.z;}

   Vector3&operator=(const double&a)   {x=a;y=a;z=a;return*this;}
   Vector3&operator=(const Vector3&v)  {x=v.x;y=v.y;z=v.z;return*this;}
   Vector3&operator+=(const Vector3&v) {x+=v.x;y+=v.y;z+=v.z;return*this;}
   Vector3&operator-=(const Vector3&v) {x-=v.x;y-=v.y;z-=v.z;return*this;}
   Vector3&operator*=(const Vector3&v) {x*=v.x;y*=v.y;z*=v.z;return*this;}
   Vector3&operator/=(const Vector3&v) {x/=v.x;y/=v.y;z/=v.z;return*this;}
   Vector3&operator+=(const double&v)  {x+=v;y+=v;z+=v;return*this;}
   Vector3&operator-=(const double&v)  {x-=v;y-=v;z-=v;return*this;}
   Vector3&operator*=(const double&v)  {x*=v;y*=v;z*=v;return*this;}
   Vector3&operator/=(const double&v)  {x/=v;y/=v;z/=v;return*this;}

   const Vector3 operator-(               )const {return Vector3(-x,-y,-z);}
   const Vector3 operator+(const Vector3&v)const {return Vector3(x+v.x,y+v.y,z+v.z);}
   const Vector3 operator-(const Vector3&v)const {return Vector3(x-v.x,y-v.y,z-v.z);}
   const Vector3 operator*(const Vector3&v)const {return Vector3(x*v.x,y*v.y,z*v.z);}
   const Vector3 operator/(const Vector3&v)const {return Vector3(x/v.x,y/v.y,z/v.z);}
   const Vector3 operator+(const double&a)const  {return Vector3(x+a,y+a,z+a);}
   const Vector3 operator-(const double&a)const  {return Vector3(x-a,y-a,z-a);}
   const Vector3 operator*(const double&a)const  {return Vector3(x*a,y*a,z*a);}
   const Vector3 operator/(const double&a)const  {return Vector3(x/a,y/a,z/a);}

   bool operator==(const Vector3&v) {return x==v.x && y==v.y && z==v.z;}
   bool operator==(const double &v) {return x==v  && y==v && z==v;}
   bool operator!=(const Vector3&v) {return x!=v.x || y!=v.y || z!=v.z;}
   bool operator!=(const double &v) {return x!=v  || y!=v || z!=v  ;}

   inline double  inner(const Vector3&v)const;
   inline Vector3 inner(const Matrix3&m)const;
   inline Matrix3 outer(const Vector3&v)const;
};

inline const Vector3 operator+(const double&a,const Vector3&v){return Vector3(a+v.x,a+v.y,a+v.z);}
inline const Vector3 operator-(const double&a,const Vector3&v){return Vector3(a-v.x,a-v.y,a-v.z);}
inline const Vector3 operator*(const double&a,const Vector3&v){return Vector3(a*v.x,a*v.y,a*v.z);}
inline const Vector3 operator/(const double&a,const Vector3&v){return Vector3(a/v.x,a/v.y,a/v.z);}

//-------------------------------------------------------------------------//
// Matrix 2D
//-------------------------------------------------------------------------//

struct Matrix2{

   double xx,xy,yx,yy;

   Matrix2(){xx=0.;xy=0.;yx=0.;yy=0.;}
   Matrix2(double a){xx=a;xy=a;yx=a;yy=a;}
   Matrix2(double a,double b,double c,double d){xx=a;xy=b;yx=c;yy=d;}
   Matrix2(const Matrix2&m){xx=m.xx;xy=m.xy;yx=m.yx;yy=m.yy;}

   Matrix2&operator=(const double&a){xx=a;xy=a;yx=a;yy=a;return*this;}
   Matrix2&operator=(const Matrix2&m){xx=m.xx;xy=m.xy;yx=m.yx;yy=m.yy;return*this;}
   Matrix2&operator+=(const Matrix2&m){xx+=m.xx;xy+=m.xy;yx+=m.yx;yy+=m.yy;return*this;}
   Matrix2&operator-=(const Matrix2&m){xx-=m.xx;xy-=m.xy;yx-=m.yx;yy-=m.yy;return*this;}
   Matrix2&operator*=(const Matrix2&m){xx*=m.xx;xy*=m.xy;yx*=m.yx;yy*=m.yy;return*this;}
   Matrix2&operator/=(const Matrix2&m){xx/=m.xx;xy/=m.xy;yx/=m.yx;yy/=m.yy;return*this;}
   Matrix2&operator+=(const double&m){xx+=m;xy+=m;yx+=m;yy+=m;return*this;}
   Matrix2&operator-=(const double&m){xx-=m;xy-=m;yx-=m;yy-=m;return*this;}
   Matrix2&operator*=(const double&m){xx*=m;xy*=m;yx*=m;yy*=m;return*this;}
   Matrix2&operator/=(const double&m){xx/=m;xy/=m;yx/=m;yy/=m;return*this;}

   const Matrix2 operator+(const Matrix2&m)const{return Matrix2(xx+m.xx,xy+m.xy,yx+m.yx,yy+m.yy);}
   const Matrix2 operator-(const Matrix2&m)const{return Matrix2(xx-m.xx,xy-m.xy,yx-m.yx,yy-m.yy);}
   const Matrix2 operator*(const Matrix2&m)const{return Matrix2(xx*m.xx,xy*m.xy,yx*m.yx,yy*m.yy);}
   const Matrix2 operator/(const Matrix2&m)const{return Matrix2(xx/m.xx,xy/m.xy,yx/m.yx,yy/m.yy);}
   const Matrix2 operator+(const double&a)const{return Matrix2(xx+a,xy+a,yx+a,yy+a);}
   const Matrix2 operator-(const double&a)const{return Matrix2(xx-a,xy-a,yx-a,yy-a);}
   const Matrix2 operator*(const double&a)const{return Matrix2(xx*a,xy*a,yx*a,yy*a);}
   const Matrix2 operator/(const double&a)const{return Matrix2(xx/a,xy/a,yx/a,yy/a);}

   bool operator==(const Matrix2&m){return xx==m.xx&&xy==m.xy&&yx==m.yx&&yy==m.yy;}
   bool operator!=(const Matrix2&m){return xx!=m.xx||xy!=m.xy||yx!=m.yx||yy!=m.yy;}

   inline Vector2 inner(const Vector2&v)const;
   inline Matrix2 inner(const Matrix2&n)const;
};

inline const Matrix2 operator+(const double&a,const Matrix2&m){return Matrix2(a+m.xx,a+m.xy,a+m.yx,a+m.yy);}
inline const Matrix2 operator-(const double&a,const Matrix2&m){return Matrix2(a-m.xx,a-m.xy,a-m.yx,a-m.yy);}
inline const Matrix2 operator*(const double&a,const Matrix2&m){return Matrix2(a*m.xx,a*m.xy,a*m.yx,a*m.yy);}
inline const Matrix2 operator/(const double&a,const Matrix2&m){return Matrix2(a/m.xx,a/m.xy,a/m.yx,a/m.yy);}

//-------------------------------------------------------------------------//
// Matrix 3D
//-------------------------------------------------------------------------//

struct Matrix3 {

   // matrix components -->
   double xx,xy,xz,yx,yy,yz,zx,zy,zz;
   //<--

   // constructors -->
   Matrix3(){xx=0.0;xy=0.0;xz=0.0;yx=0.0;yy=0.0;yz=0.0;zx=0.0;zy=0.0;zz=0.0;}
   Matrix3(double a){xx=a;xy=a;xz=a;yx=a;yy=a;yz=a;zx=a;zy=a;zz=a;}
   Matrix3(double mxx,double mxy,double mxz,double myx,double myy,double myz,double mzx,double mzy,double mzz){
      xx=mxx; xy=mxy; xz=mxz; yx=myx; yy=myy; yz=myz; zx=mzx; zy=mzy; zz=mzz;}
   Matrix3(const Matrix3&m){xx=m.xx; xy=m.xy; xz=m.xz; yx=m.yx; yy=m.yy; yz=m.yz; zx=m.zx; zy=m.zy; zz=m.zz;}
   //<--

   // matrix operators -->
   Matrix3&operator =(const Matrix3&m)
   { xx=m.xx; xy=m.xy; xz=m.xz; yx=m.yx; yy=m.yy; yz=m.yz; zx=m.zx; zy=m.zy; zz=m.zz; return*this;}
   Matrix3&operator+=(const Matrix3&m)
   { xx+=m.xx; xy+=m.xy; xz+=m.xz; yx+=m.yx; yy+=m.yy; yz+=m.yz; zx+=m.zx; zy+=m.zy; zz+=m.zz; return*this;}
   Matrix3&operator-=(const Matrix3&m)
   { xx-=m.xx; xy-=m.xy; xz-=m.xz; yx-=m.yx; yy-=m.yy; yz-=m.yz; zx-=m.zx; zy-=m.zy; zz-=m.zz; return*this;}
   Matrix3&operator*=(const Matrix3&m)
   { xx*=m.xx; xy*=m.xy; xz*=m.xz; yx*=m.yx; yy*=m.yy; yz*=m.yz; zx*=m.zx; zy*=m.zy; zz*=m.zz; return*this;}
   Matrix3&operator/=(const Matrix3&m)
   { xx/=m.xx; xy/=m.xy; xz/=m.xz; yx/=m.yx; yy/=m.yy; yz/=m.yz; zx/=m.zx; zy/=m.zy; zz/=m.zz; return*this;}

   const Matrix3 operator+(const Matrix3&m)const
   {return Matrix3(xx+m.xx,xy+m.xy,xz+m.xz,yx+m.yx,yy+m.yy,yz+m.yz,zx+m.zx,zy+m.zy,zz+m.zz);}
   const Matrix3 operator-(const Matrix3&m)const
   {return Matrix3(xx-m.xx,xy-m.xy,xz-m.xz,yx-m.yx,yy-m.yy,yz-m.yz,zx-m.zx,zy-m.zy,zz-m.zz);}
   const Matrix3 operator*(const Matrix3&m)const
   {return Matrix3(xx*m.xx,xy*m.xy,xz*m.xz,yx*m.yx,yy*m.yy,yz*m.yz,zx*m.zx,zy*m.zy,zz*m.zz);}
   const Matrix3 operator/(const Matrix3&m)const
   {return Matrix3(xx/m.xx,xy/m.xy,xz/m.xz,yx/m.yx,yy/m.yy,yz/m.yz,zx/m.zx,zy/m.zy,zz/m.zz);}
   //<--

   // double operators
   Matrix3&operator+=(const double&m){xx+=m;xy+=m;xz+=m;yx+=m;yy+=m;yz+=m;zx+=m;zy+=m;zz+=m;return*this;}
   Matrix3&operator-=(const double&m){xx-=m;xy-=m;xz-=m;yx-=m;yy-=m;yz-=m;zx-=m;zy-=m;zz-=m;return*this;}
   Matrix3&operator*=(const double&m){xx*=m;xy*=m;xz*=m;yx*=m;yy*=m;yz*=m;zx*=m;zy*=m;zz*=m;return*this;}
   Matrix3&operator/=(const double&m){xx/=m;xy/=m;xz/=m;yx/=m;yy/=m;yz/=m;zx/=m;zy/=m;zz/=m;return*this;}

   const Matrix3 operator+(const double&a)const{return Matrix3(xx+a,xy+a,xz+a,yx+a,yy+a,yz+a,zx+a,zy+a,zz+a);}
   const Matrix3 operator-(const double&a)const{return Matrix3(xx-a,xy-a,xz-a,yx-a,yy-a,yz-a,zx-a,zy-a,zz-a);}
   const Matrix3 operator*(const double&a)const{return Matrix3(xx*a,xy*a,xz*a,yx*a,yy*a,yz*a,zx*a,zy*a,zz*a);}
   const Matrix3 operator/(const double&a)const{return Matrix3(xx/a,xy/a,xz/a,yx/a,yy/a,yz/a,zx/a,zy/a,zz/a);}
   //<--

   // boolean operations -->
   bool operator==(const Matrix3&m)
   {return xx==m.xx && xy==m.xy && xz==m.xz && yx==m.yx && yy==m.yy && yz==m.yz && zx==m.zx && zy==m.zy && zz==m.zz;}
   bool operator!=(const Matrix3&m)
   {return xx!=m.xx || xy!=m.xy || xz!=m.xz || yx!=m.yx || yy!=m.yy || yz!=m.yz || zx!=m.zx || zy!=m.zy || zz!=m.zz;}
   //<--

   inline Vector3 inner(const Vector3&v)const;
   inline Matrix3 inner(const Matrix3&m)const;
};

inline const Matrix3 operator+(const double&a,const Matrix3&m)
{return Matrix3(m.xx+a,m.xy+a,m.xz+a,m.yx+a,m.yy+a,m.yz+a,m.zx+a,m.zy+a,m.zz+a);}

inline const Matrix3 operator-(const double&a,const Matrix3&m)
{return Matrix3(m.xx-a,m.xy-a,m.xz-a,m.yx-a,m.yy-a,m.yz-a,m.zx-a,m.zy-a,m.zz-a);}

inline const Matrix3 operator*(const double&a,const Matrix3&m)
{return Matrix3(m.xx*a,m.xy*a,m.xz*a,m.yx*a,m.yy*a,m.yz*a,m.zx*a,m.zy*a,m.zz*a);}

inline const Matrix3 operator/(const double&a,const Matrix3&m)
{return Matrix3(m.xx*1/a,m.xy*1/a,m.xz*1/a,m.yx*1/a,m.yy*1/a,m.yz*1/a,m.zx*1/a,m.zy*1/a,m.zz*1/a);}

//-------------------------------------------------------------------------//
// unary operators
//-------------------------------------------------------------------------//

// square lenght
inline double len2(const Vector2&u){return u.x*u.x+u.y*u.y;}
inline double len2(const Vector3&u){return u.x*u.x+u.y*u.y+u.z*u.z;}

// lenght
inline double len(const Vector2&u){return sqrt(u.x*u.x+u.y*u.y);}
inline double len(const Vector3&u){return sqrt(u.x*u.x+u.y*u.y+u.z*u.z);}

// trace
inline double trace(const Matrix2&m){return m.xx+m.yy;}
inline double trace(const Matrix3&m){return m.xx+m.yy+m.zz;}

// mean stress
inline double mean(const Matrix2&m){return trace(m)/2.0;}
inline double mean(const Matrix3&m){return trace(m)/3.0;}

// determinans -->
inline double det(const Matrix2&m){return m.xx*m.yy-m.xy*m.yx;}
inline double det(const Matrix3&m){
   double m11=m.xx;double m12=m.xy;double m13=m.xz;
   double m21=m.yx;double m22=m.yy;double m23=m.yz;
   double m31=m.zx;double m32=m.zy;double m33=m.zz;
   return (m11*m22*m33 - m11*m23*m32 - m12*m21*m33 + m12*m23*m31 + m13*m21*m32 - m13*m22*m31);
}
inline double determinant(const Matrix3&m){ return det(m); }
//<--

// query
inline double absMax(const Matrix2&s){return std::max(std::max(abs(s.xx),abs(s.xy)),std::max(abs(s.yx),abs(s.yy)));}

// indentity -->
inline Matrix2 identity2D(const double f=1.){return Matrix2(f,0.,0.,f);}
inline Matrix3 identity3D(const double f=1.){return Matrix3(f,0.,0.,0.,f,0.,0.,0.,f);}
inline Matrix2 I2(const double f=1.){return identity2D(f);}
inline Matrix3 I3(const double f=1.){return identity3D(f);}
inline Matrix3 I (const double f=1.){return identity3D(f);}
//<--

// transpose
inline Matrix2 transpose(const Matrix2&m){return Matrix2(m.xx,m.yx,m.xy,m.yy);}
inline Matrix3 transpose(const Matrix3&m){
   double m11=m.xx;double m12=m.xy;double m13=m.xz;
   double m21=m.yx;double m22=m.yy;double m23=m.yz;
   double m31=m.zx;double m32=m.zy;double m33=m.zz;
   return Matrix3(m11, m21, m31,m12, m22, m32,m13, m23, m33);
}

// sum
inline Matrix2 sum(const Matrix2&m1, const Matrix2&m2)
{
   return Matrix2(m1.xx+m2.xx,m1.xy+m2.xy,m1.yx+m2.yx,m1.yy+m2.yy);
}

inline Matrix3 sum(const Matrix3&a, const Matrix3&b){
   double a11=a.xx;double a12=a.xy;double a13=a.xz;
   double a21=a.yx;double a22=a.yy;double a23=a.yz;
   double a31=a.zx;double a32=a.zy;double a33=a.zz;
   double b11=b.xx;double b12=b.xy;double b13=b.xz;
   double b21=b.yx;double b22=b.yy;double b23=b.yz;
   double b31=b.zx;double b32=b.zy;double b33=b.zz;
   return Matrix3(a11+b11,a12+b12,a13+b13,a21+b21,a22+b22,a23+b23,a31+b31,a32+b32,a33+b33);
}

inline double sum(const Vector3&v){ return (v.x+v.y+v.z);}
inline double sum(const Vector2&v){ return (v.x+v.y);}

// rest
inline Matrix2 sub(const Matrix2&m1, const Matrix2&m2)
{return Matrix2(m1.xx-m2.xx,m1.xy-m2.xy,m1.yx-m2.yx,m1.yy-m2.yy);}

inline Matrix3 sub(const Matrix3&a, const Matrix3&b){
   double a11=a.xx;double a12=a.xy;double a13=a.xz;
   double a21=a.yx;double a22=a.yy;double a23=a.yz;
   double a31=a.zx;double a32=a.zy;double a33=a.zz;
   double b11=b.xx;double b12=b.xy;double b13=b.xz;
   double b21=b.yx;double b22=b.yy;double b23=b.yz;
   double b31=b.zx;double b32=b.zy;double b33=b.zz;
   return Matrix3(a11-b11,a12-b12,a13-b13,a21-b21,a22-b22,a23-b23,a31-b31,a32-b32,a33-b33);
}

// inverse -->
inline Matrix2 inv(const Matrix2&m){const double d=det(m);return Matrix2(m.yy/d,-m.xy/d,-m.yx/d,m.xx/d);}

inline Matrix3 inv(const Matrix3&m){
   double m11=m.xx;double m12=m.xy;double m13=m.xz;
   double m21=m.yx;double m22=m.yy;double m23=m.yz;
   double m31=m.zx;double m32=m.zy;double m33=m.zz;
   double inv11 = (m22*m33-m23*m32)/(m11*m22*m33-m11*m23*m32-m12*m21*m33+m12*m23*m31+m13*m21*m32-m13*m22*m31);
   double inv12 =-(m12*m33-m13*m32)/(m11*m22*m33-m11*m23*m32-m12*m21*m33+m12*m23*m31+m13*m21*m32-m13*m22*m31);
   double inv13 = (m12*m23-m13*m22)/(m11*m22*m33-m11*m23*m32-m12*m21*m33+m12*m23*m31+m13*m21*m32-m13*m22*m31);
   double inv21 =-(m21*m33-m23*m31)/(m11*m22*m33-m11*m23*m32-m12*m21*m33+m12*m23*m31+m13*m21*m32-m13*m22*m31);
   double inv22 = (m11*m33-m13*m31)/(m11*m22*m33-m11*m23*m32-m12*m21*m33+m12*m23*m31+m13*m21*m32-m13*m22*m31);
   double inv23 =-(m11*m23-m13*m21)/(m11*m22*m33-m11*m23*m32-m12*m21*m33+m12*m23*m31+m13*m21*m32-m13*m22*m31);
   double inv31 = (m21*m32-m22*m31)/(m11*m22*m33-m11*m23*m32-m12*m21*m33+m12*m23*m31+m13*m21*m32-m13*m22*m31);
   double inv32 =-(m11*m32-m12*m31)/(m11*m22*m33-m11*m23*m32-m12*m21*m33+m12*m23*m31+m13*m21*m32-m13*m22*m31);
   double inv33 = (m11*m22-m12*m21)/(m11*m22*m33-m11*m23*m32-m12*m21*m33+m12*m23*m31+m13*m21*m32-m13*m22*m31);
   return Matrix3(inv11,inv12,inv13,inv21,inv22,inv23,inv31,inv32,inv33);
}
//<--

// multiplication -->
inline Matrix3 mmult(const Matrix3&a,const Matrix3&b)
{
   double a11=a.xx;double a12=a.xy;double a13=a.xz;
   double a21=a.yx;double a22=a.yy;double a23=a.yz;
   double a31=a.zx;double a32=a.zy;double a33=a.zz;
   double b11=b.xx;double b12=b.xy;double b13=b.xz;
   double b21=b.yx;double b22=b.yy;double b23=b.yz;
   double b31=b.zx;double b32=b.zy;double b33=b.zz;
   double ab11=a11*b11 + a12*b21 + a13*b31;
   double ab12=a11*b12 + a12*b22 + a13*b32;
   double ab13=a11*b13 + a12*b23 + a13*b33;
   double ab21=a21*b11 + a22*b21 + a23*b31;
   double ab22=a21*b12 + a22*b22 + a23*b32;
   double ab23=a21*b13 + a22*b23 + a23*b33;
   double ab31=a31*b11 + a32*b21 + a33*b31;
   double ab32=a31*b12 + a32*b22 + a33*b32;
   double ab33=a31*b13 + a32*b23 + a33*b33;
   return Matrix3(ab11,ab12,ab13,ab21,ab22,ab23,ab31,ab32,ab33);
}
//<--

// deviatoric -->
inline Matrix2 dev(const Matrix2&m){const double p=(m.xx+m.yy)/2.0;return Matrix2(m.xx-p,m.xy,m.yx,m.yy-p);}

inline Matrix3 deviatoric(const Matrix3&m){
   double p=trace(m)/3.0;
   return Matrix3(m.xx-p,m.xy,m.xz,m.yx,m.yy-p,m.yz,m.zx,m.zy,m.zz-p);
}
inline Matrix2 deviatoric(const Matrix2&m){
   return dev(m);
}
//<--

// Frobenius norm or Euclidean norm
inline double ddot(const Matrix2&m){return sqrt(m.xx*m.xx+m.xy*m.yx+m.yx*m.xy+m.yy*m.yy);}

// Frobenius norm or Euclidean norm -->
inline double frobeniusnorm(const Matrix3&m){

   double m11=m.xx;double m12=m.xy;double m13=m.xz;
   double m21=m.yx;double m22=m.yy;double m23=m.yz;
   double m31=m.zx;double m32=m.zy;double m33=m.zz;

   return std::sqrt(std::abs(m11)*std::abs(m11)+std::abs(m12)*std::abs(m12)+std::abs(m13)*std::abs(m13)+std::abs(m21)*std::abs(m21)+std::abs(m22)*std::abs(m22)+std::abs(m23)*std::abs(m23)+std::abs(m31)*std::abs(m31)+std::abs(m32)*std::abs(m32)+std::abs(m33)*std::abs(m33));
}
//<--

// double tensor contraction-->
inline double tensorddot(const Matrix2&a,const Matrix2&b ){
   return (a.xx*b.xx + a.xy*b.yx + a.yx*b.xy + a.yy*b.yy ) ;
}

/*(Holzapfel) A:B := Trace(A*transpose(B)) = A(ij)B(ij)*/
inline double tensorddot(const Matrix3&a,const Matrix3&b ){
   double a11=a.xx;double a12=a.xy;double a13=a.xz;
   double a21=a.yx;double a22=a.yy;double a23=a.yz;
   double a31=a.zx;double a32=a.zy;double a33=a.zz;
   double b11=b.xx;double b12=b.xy;double b13=b.xz;
   double b21=b.yx;double b22=b.yy;double b23=b.yz;
   double b31=b.zx;double b32=b.zy;double b33=b.zz;
   return (a11*b11 + a12*b12 + a13*b13 + a21*b21 + a22*b22 + a23*b23 + a31*b31 + a32*b32 + a33*b33);
}
//<--

// invariants -->
inline double invariant1(Matrix3 a){
	return trace(a);
}
inline double invariant2(Matrix3 a){
	return (a.xx*a.yy+a.yy*a.zz+a.zz*a.xx-(a.xy*a.yx+a.yz*a.zy+a.xz*a.zx));
}
inline double invariant3(Matrix3 a){
	return determinant(a);
}
//<--

// Printing operations -->
inline std::ostream& operator<<(std::ostream&os,const Vector2&v){os<<v.x<<", "<<v.y;return os;}
inline std::istream& operator>>(std::istream&is,Vector2&v){is>>v.x     >>v.y;return is;}
inline std::ostream& operator<<(std::ostream&os,const Matrix2&m){os<<m.xx<<", "<<m.xy<<", "<<m.yx<<", "<<m.yy;return os;}
inline std::istream& operator>>(std::istream&is,Matrix2&v){is>>v.xx    >>v.xy     >>v.yx     >>v.yy;return is;}
inline std::ostream& operator<<(std::ostream&os,const Vector3&v){os<<v.x<<", "<<v.y<<", "<<v.z;return os;}
inline std::istream& operator>>(std::istream&is,Vector3&v){is>>v.x>>v.y>>v.z;return is;}
inline std::ostream& operator<<(std::ostream&os,const Matrix3&m)
{os<<m.xx<<", "<<m.xy<<", "<<m.xz<<", "<<m.yx<<", "<<m.yy<<", "<<m.yz<<", "<<m.zx<<", "<<m.zy<<", "<<m.zz;return os;}
inline std::istream& operator>>(std::istream&is,Matrix3&v){is>>v.xx>>v.xy>>v.xz>>v.yx>>v.yy>>v.yz>>v.zx>>v.zy>>v.zz;return is;}
//<--

// magnitude of a vector -->
inline double magnitude(const Vector2&u){return std::sqrt(u.x*u.x+u.y*u.y);}
inline double magnitude(const Vector3&u){return std::sqrt(u.x*u.x+u.y*u.y+u.z*u.z);}
inline double mag(const Vector2&u){return magnitude(u);}
inline double mag(const Vector3&u){return magnitude(u);}
//<--

// distance -->
inline double EuclidDistance( Vector3 u, Vector3 v){
	return sqrt(pow((u.x-v.x),2.0) + pow((u.y-v.y),2.0) + pow((u.z-v.z),2.0));
}
inline double TwoPointsDistance( Vector3 u, Vector3 v){
   return EuclidDistance(u,v);
}
//<--

// utit vector -->
inline Vector2 unit(const Vector2&u){assert(abs(u.x)>MACHINE_TOL||abs(u.y)>MACHINE_TOL);return u/magnitude(u);}
inline Vector3 unit(const Vector3&u){assert(abs(u.x)>MACHINE_TOL||abs(u.y)>MACHINE_TOL||abs(u.z)>MACHINE_TOL);return u/magnitude(u);}
//<--

inline double radius(const Vector2&u,const Vector2&v){
   const double s=(u.x-v.x)*(u.x-v.x)+(u.y-v.y)*(u.y-v.y);
   return(s>MACHINE_TOL?std::sqrt(s):0.) ;
}

// return the effectove von mises stress -->
inline double vonMises(const Matrix3&s){
	Matrix3 sd  = deviatoric(s);
	return std::sqrt(3.0/2.0*tensorddot(sd,sd));
}
//<--

// principal values of a tensor-->
inline Vector3 principalValues(const Matrix3 s){

	using namespace std;

	Vector3 val(0.0);

	double I1 = invariant1(s);
	double I2 = invariant2(s);
	double I3 = invariant3(s);

	double Q  = (3.0*I2-I1*I1)/9.0;
	double R  = (2.0*(I1*I1*I1)-9.0*I1*I2+27.0*I3)/54.0;
	double th = acos(R/sqrt(-pow(Q,3.0)));

	val.x = 2.0*sqrt(-Q)*cos((th+0.0*pi)/3.0)+I1/3.0;
	val.y = 2.0*sqrt(-Q)*cos((th+2.0*pi)/3.0)+I1/3.0;
	val.z = 2.0*sqrt(-Q)*cos((th+4.0*pi)/3.0)+I1/3.0;

	val.x = isnan(val.x) ? 0.0:val.x;
	val.y = isnan(val.y) ? 0.0:val.y;
	val.z = isnan(val.z) ? 0.0:val.z;

	return val;
}
//<--

// principal values of a tensor
inline Vector3 pst(const Matrix3&s){
	return  principalValues(s);
}

// principal values of a tensor
inline Vector3 PincipalValuesGet(const Matrix3&s)
{
   return  principalValues(s);
}


// left-to-right operators where order matters -->
inline double  Vector2::inner(const Vector2&v)const{return x*v.x+y*v.y;}
inline double  Vector3::inner(const Vector3&v)const{return x*v.x+y*v.y+z*v.z;}
inline Vector2 Vector2::inner(const Matrix2&m)const{return Vector2(m.xx*x+m.yx*y,m.xy*x+m.yy*y);}
inline Vector3 Vector3::inner(const Matrix3&m)const{
   return  Vector3(m.xx*x + m.yx*y + m.zx*z, m.xy*x + m.yy*y + m.zy*z, m.xz*x + m.yz*y + m.zz*z);
}
inline Matrix2 Vector2::outer(const Vector2&v)const{return Matrix2(x*v.x,x*v.y,y*v.x,y*v.y);}
inline Vector2 Matrix2::inner(const Vector2&v)const{return Vector2(xx*v.x+xy*v.y,yx*v.x+yy*v.y);}
inline Vector3 Matrix3::inner(const Vector3&v)const{
	return Vector3( xx*v.x + xy*v.y + xz*v.z,yx*v.x + yy*v.y + yz*v.z,zx*v.x + zy*v.y + zz*v.z);}

inline Matrix2 Matrix2::inner(const Matrix2&n)const
{
	return Matrix2(xx*n.xx+xy*n.yx,xx*n.xy+xy*n.yy,yx*n.xx+yy*n.yx,yx*n.xy+yy*n.yy);
}

inline Matrix3 Matrix3::inner(const Matrix3&n)const{ Matrix3 a(xx,xy,xz,yx,yy,yz,zx,zy,zz); return mmult(a,n);};
//<--

void eigen(Matrix3 s, Vector3& eival, Matrix3& eivec);
void eigen(Matrix3 s, Vector3& eival,Vector3& n1,Vector3& n2,Vector3& n3);

#endif /* INC_TENSOR_H */
