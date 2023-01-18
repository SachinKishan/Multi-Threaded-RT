#pragma once

#include <cmath>
#include "VectorMath.h"

template<typename T>
class Matrix44 
{
public:

    T x[4][4] = { {1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1} };

    Matrix44()
    {
	    
    }

    Matrix44
       (T a, T b, T c, T d, 
        T e, T f, T g, T h,
        T i, T j, T k, T l, 
        T m, T n, T o, T p=1)
    {
        x[0][0] = a;
        x[0][1] = b;
        x[0][2] = c;
        x[0][3] = d;
        x[1][0] = e;
        x[1][1] = f;
        x[1][2] = g;
        x[1][3] = h;
        x[2][0] = i;
        x[2][1] = j;
        x[2][2] = k;
        x[2][3] = l;
        x[3][0] = m;
        x[3][1] = n;
        x[3][2] = o;
        x[3][3] = p;
    }

    const T* operator [] (uint8_t i) const { return x[i]; }
    T* operator [] (uint8_t i) { return x[i]; }

    static void multiply(const Matrix44<T>& a, const Matrix44& b, Matrix44& c)
    {
        for (uint8_t i = 0; i < 4; ++i) {
            for (uint8_t j = 0; j < 4; ++j) {
                c[i][j] = a[i][0] * b[0][j] + a[i][1] * b[1][j] +
                    a[i][2] * b[2][j] + a[i][3] * b[3][j];
            }
        }
    }

    friend Matrix44 operator *(Matrix44<double> a,Matrix44<double> b)
    {
        Matrix44 c;
        c.multiply(a, b, c);
        return c;
    }

    friend Matrix44 operator *(Matrix44<double> a,double f)
    {
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)a[i][j]=a[i][j] * f;
        return a;
    }
    void out()
    {
	    for(int i=0;i<4;i++)
	    {
		    for(int j=0;j<4;j++)
		    {
                std::cout << x[i][j]<<" ";
		    }
            std::cout << std::endl;
	    }
    }

    Matrix44 returnTransposed() const
    {
        return Matrix44(x[0][0],
            x[1][0],
            x[2][0],
            x[3][0],
            x[0][1],
            x[1][1],
            x[2][1],
            x[3][1],
            x[0][2],
            x[1][2],
            x[2][2],
            x[3][2],
            x[0][3],
            x[1][3],
            x[2][3],
            x[3][3]);
    }

    Matrix44& transpose()//operation
    {
        Matrix44 tmp(x[0][0],
            x[1][0],
            x[2][0],
            x[3][0],
            x[0][1],
            x[1][1],
            x[2][1],
            x[3][1],
            x[0][2],
            x[1][2],
            x[2][2],
            x[3][2],
            x[0][3],
            x[1][3],
            x[2][3],
            x[3][3]);
        *this = tmp;

        return *this;
    }
    
    void multiplyDirMatrix(const vec3& src, vec3& dst) const
    {
        double a, b, c;

        a = src[0] * x[0][0] + src[1] * x[1][0] + src[2] * x[2][0];
        b = src[0] * x[0][1] + src[1] * x[1][1] + src[2] * x[2][1];
        c = src[0] * x[0][2] + src[1] * x[1][2] + src[2] * x[2][2];

        dst[0] = a;
        dst[1] = b;
        dst[2] = c;
    }
   
    void multiplyVecMatrix(vec3& src, vec3& dst) const
    {
        double a, b, c, w;

        a = src[0] * x[0][0] + src[1] * x[1][0] + src[2] * x[2][0] + x[3][0];
        b = src[0] * x[0][1] + src[1] * x[1][1] + src[2] * x[2][1] + x[3][1];
        c = src[0] * x[0][2] + src[1] * x[1][2] + src[2] * x[2][2] + x[3][2];
        w = src[0] * x[0][3] + src[1] * x[1][3] + src[2] * x[2][3] + x[3][3];

        dst[0] = a / w;
        dst[1] = b / w;
        dst[2] = c / w;
    }

    vec3 multiplyVectorMatrix(vec3 src)//src=source vector
    {
        vec3 destination;
        double a, b, c, w;
        //vector is in row form. Vector * Matrix (V*M)
        a = src.x() * x[0][0] + src.y() * x[1][0] + src.z() * x[2][0] + x[3][0];
        b = src.x() * x[0][1] + src.y() * x[1][1] + src.z() * x[2][1] + x[3][1];
        c = src.x() * x[0][2] + src.y() * x[1][2] + src.z() * x[2][2] + x[3][2];
        w = src.x() * x[0][3] + src.y() * x[1][3] + src.z()* x[2][3] + x[3][3];

        destination[0] = a / w;
        destination[1] = b / w;
        destination[2] = c / w;
        
        
        return destination;

    }


    Matrix44<double> invert()
    {
        Matrix44<double> inverse;
        double determinant =
            x[0][0]
				*( AddAndMultiply(1,1,2,2,3,3,1,2,2,3,3,1,1,3,2,1,3,2,1,3,2,2,3,1,1,2,2,1,3,3,1,1,2,3,3,2));

        determinant -= x[1][0] *
				( AddAndMultiply(0,1,2,2,3,3,0,2,2,3,3,1,0,3,2,1,3,2,0,3,2,2,3,1,0,2,2,1,3,3,0,1,2,3,3,2));

        determinant+=x[2][0]*
				( AddAndMultiply(0,1,1,2,3,3,0,2,1,3,3,1,0,3,1,1,3,2,0,3,1,2,3,1,0,2,1,1,3,3,0,1,1,3,3,2));

        determinant-=x[3][0]*
				( AddAndMultiply(0,1,1,2,2,3,0,2,1,3,2,1,0,3,1,1,2,2,0,3,1,2,2,1,0,2,1,1,2,3,0,1,1,3,2,2));

    	inverse[0][0]=
				( x[1][1] * x[2][2] * x[3][3]
                + x[1][2] * x[2][3] * x[3][1]
                + x[1][3] * x[2][1] * x[3][2]
                - x[1][3] * x[2][2] * x[3][1]
                - x[1][2] * x[2][1] * x[3][3]
                - x[1][1] * x[2][3] * x[3][2]);
#pragma region
        inverse[0][0] = AddAndMultiply(1,1,2,2,3,3,1,2,2,3,3,1,1,3,2,1,3,2,1,3,2,2,3,1,1,2,2,1,3,3,1,1,2,3,3,2);
        inverse[0][1] -= AddAndMultiply(0, 1, 2, 2, 3, 3, 0, 2, 2, 3, 3, 1, 0, 3, 2, 1, 3, 2, 0, 3, 2, 2, 3, 1, 0, 2, 2, 1, 3, 3, 0, 1, 2, 3, 3, 2);
        inverse[0][2] = AddAndMultiply(0,1,1,2,3,3,0,2,1,3,3,1,0,3,1,1,3,2,0,3,1,2,3,1,0,2,1,1,3,3,0,1,1,3,3,2);
        inverse[0][3] -= AddAndMultiply(0,1,1,2,2,3,0,2,1,3,2,1,0,3,1,1,2,2,0,3,1,2,2,1,0,2,1,1,2,3,0,1,1,3,2,2);
        inverse[1][0] -= AddAndMultiply(0,1,2,2,3,3,1,2,2,3,3,0,1,3,2,0,3,2,1,3,2,2,3,0,1,2,2,0,3,3,1,0,2,3,3,2);
    	inverse[1][1] = AddAndMultiply(0, 0, 2, 2, 3, 3, 0, 2, 2, 3, 3, 0, 0, 3, 2, 0, 3, 2, 0, 3, 2, 2, 3, 0, 0, 2, 2, 0, 3, 3, 0, 0, 2, 3, 3, 2);
        inverse[1][2] -= AddAndMultiply(0,0,1,2,3,3,0,2,1,3,3,0,0,3,1,0,3,2,0,3,1,2,3,0,0,2,1,0,3,3,0,0,1,3,3,2);
        inverse[1][3] = AddAndMultiply(0,0,1,2,2,3,0,2,1,3,2,0,0,3,1,0,2,2,0,3,1,2,2,0,0,2,1,0,2,3,0,0,1,3,2,2);
        inverse[2][0] = AddAndMultiply(1,0,2,1,3,3,1,1,2,3,3,0,1,3,2,0,3,1,1,3,2,1,3,0,1,1,2,0,3,3,1,0,2,3,3,1);
        inverse[2][1] -= AddAndMultiply(0,0,2,1,3,3,0,1,2,3,3,0,0,3,2,0,3,1,0,3,2,1,3,0,0,1,2,0,3,3,0,0,2,3,3,1);
        inverse[2][2] = AddAndMultiply(0, 0, 1, 1, 3, 3, 0, 1, 1, 3, 3, 0, 0, 3, 1, 0, 3, 1, 0, 3, 1, 1, 3, 0, 0, 1, 1, 0, 3, 3, 0, 0, 1, 3, 3, 1);
    	inverse[2][3] -= AddAndMultiply(0, 0, 1, 1, 2, 3, 0, 1, 1, 3, 2, 0, 0, 3, 1, 0, 2, 1, 0, 3, 1, 1, 2, 0, 0, 1, 1, 0, 2, 3, 0, 0, 1, 3, 2, 1);
        inverse[3][0] -= AddAndMultiply(1,0,2,1,3,2,1,1,2,2,3,0,1,2,2,0,3,1,1,2,2,1,3,0,1,1,2,0,3,2,1,0,2,2,3,1);
        inverse[3][1] = AddAndMultiply(0, 0, 2, 1, 3, 2, 0, 1, 2, 2, 3, 0, 0, 2, 2, 0, 3, 1, 0, 2, 2, 1, 3, 0, 0, 1, 2, 0, 3, 2, 0, 0, 2, 2, 3, 1);
        inverse[3][2] -= AddAndMultiply(0, 0, 1, 1, 3, 2, 0, 1, 1, 2, 3, 0, 0, 2, 1, 0, 3, 1, 0, 2, 1, 1, 3, 0, 0, 1, 1, 0, 3, 2, 0, 0, 1, 2, 3, 1);
        inverse[3][3] = AddAndMultiply(0, 0, 1, 1, 2, 2, 0, 1, 1, 2, 2, 0, 0, 2, 1, 0, 2, 1, 0, 2, 1, 1, 2, 0, 0, 1, 1, 0, 2, 2, 0, 0, 1, 2, 2, 1);
        
#pragma endregion

        //std::cout << std::endl << determinant;
        inverse = inverse * determinant;
        return inverse;
    }
	double AddAndMultiply(
        int a1, int a2, int a3, int a4, int a5, int a6,
        int a7, int a8, int a9, int a10, int a11, int a12,
        int a13, int a14, int a15, int a16, int a17, int a18,
        int a19, int a20, int a21, int a22, int a23, int a24,
        int a25, int a26, int a27, int a28, int a29, int a30,
        int a31, int a32, int a33, int a34, int a35, int a36)
    {
        return x[a1][a2] * x[a3][a4] * x[a5][a6] +
            x[a7][a8] * x[a9][a10] * x[a11][a12] +
            x[a13][a14] * x[a15][a16] * x[a17][a18] -
            x[a19][a20] * x[a21][a22] * x[a23][a24] -
            x[a25][a26] * x[a27][a28] * x[a29][a30] -
            x[a31][a32] * x[a33][a34] * x[a35][a36];
    }

    void lookat(vec3 from, vec3 to, vec3 up)
    {
        vec3 forward = from - to;
        forward=unit_vector(forward);
        
        vec3 right =unit_vector(cross(up, forward));
        
        vec3 newup =unit_vector(cross(forward, right));

        x[0][0] = right.x(), x[0][1] = right.y(), x[0][2] = right.z();
        x[1][0] = newup.x(), x[1][1] = newup.y(), x[1][2] = newup.z();
        x[2][0] = forward.x(), x[2][1] = forward.y(), x[2][2] = forward.z();
        x[3][0] = -from.z(), x[3][1] = -from.y(), x[3][2] = -from.z();
    }
};



class TranslationMatrix :public Matrix44<double>
{
public:

    
    TranslationMatrix(double _x, double _y,double _z)
		:Matrix44<double>(1, 0, 0, 0,
						 0, 1, 0, 0,
						 0, 0, 1, 0,
						 _x, _y, _z)
    {
        
    }

    TranslationMatrix(vec3 translate)
        :Matrix44<double>(1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1, 0,
            translate.x(), translate.y(), translate.z())
    {

    }

};


class RotationMatrixX :public Matrix44<double>
{
public:
    RotationMatrixX() = default;
    RotationMatrixX(double degrees)
        :Matrix44<double>
			(1, 0, 0, 0,
            0, cos(degrees), -sin(degrees), 0,
            0, sin(degrees), cos(degrees), 0,
            0, 0, 0){}
};

class RotationMatrixY :public Matrix44<double>
{
public:

    RotationMatrixY() = default;
    RotationMatrixY(double degrees)
        :Matrix44<double>
        (cos(degrees), 0, sin(degrees), 0,
            0, 1,0, 0,
            -sin(degrees), 0, cos(degrees), 0,
            0, 0, 0) {}
};

class RotationMatrixZ :public Matrix44<double>
{
public:

    RotationMatrixZ() = default;
    RotationMatrixZ(double degrees)
        :Matrix44<double>
        (cos(degrees), -sin(degrees), 0, 0,
            sin(degrees), cos(degrees), 0, 0,
            0, 0, 1, 0,
            0, 0, 0) {}
};

class ScaleMatrix:public Matrix44<double>
{
public:
    ScaleMatrix() = default;
    ScaleMatrix(vec3 scale):
		Matrix44<double>(
        scale.x(), 0, 0, 0,
        0,scale.y(), 0, 0,
        0,0,scale.z(), 0,
        0,0,0)
    {
	    
    }
    ScaleMatrix(double scale) :
        Matrix44<double>(
            scale, 0, 0, 0,
            0, scale, 0, 0,
            0, 0, scale, 0,
            0, 0, 0)
    {

    }

    
};

