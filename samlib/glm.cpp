#include "glm.h"
#include <math.h>
#include <iostream>
using namespace std;

/**
 * @brief GLM::GLM
 */
GLM::GLM(){

}

/**
 * @brief GLM::~GLM
 */
GLM::~GLM(){

}

/**
 * @brief GLM::linkfunc
 * @param u
 * @return
 */
double GLM::linkfunc(double u){
    return log(u);
}

/**
 * @brief GLM::ddlinkfunc
 * @param u
 * @return
 */
double GLM::dlinkfunc(double u){
    if (u==0) return 0;
    return 1/u;
}

/**
 * @brief GLM::ddlinkfunc
 * @param u
 * @return
 */
double GLM::ddlinkfunc(double u){
    if (u==0) return 0;
    return -1/(u*u);
}

/**
 * @brief GLM::invlinkfunc
 * @param Xb
 * @return
 */
double GLM::invlinkfunc(double Xb){
    return exp(Xb);
}

/**
 * @brief GLM::irls
 * @param X
 * @param y
 * @param M
 * @param N
 * @param b
 * @param maxIter
 * @param epsilon
 * @return
 */
void GLM::irls(const double *X, const double *y, int N,
                 double *b, int maxIter, double epsilon){
    int iter=0;
    double  res=1e+100;
    double *eta=new double[N];
    double *z=new double[N];
    double *u=new double[N];
    double *Xb=new double[N];
    double *w=new double[N];

    for (int i=0; i<N; i++){
        eta[i]=linkfunc(y[i]);
    }

    for (int j=0; j<=1; j++){
        b[j]=0;
    }

    while (iter<maxIter && res>epsilon){
        // X*b
        for (int i=0; i<N; i++){
            Xb[i]=b[0];
            Xb[i]+=X[i]*b[1];
        }

        // z=inv_g(Xb)
        for (int i=0; i<N; i++){
            u[i]=invlinkfunc(Xb[i]);
        }

        // compute the pseudo-data zi
        for (int i=0; i<N; i++){
            z[i]=eta[i]+dlinkfunc(u[i])*(y[i]-u[i]);
        }
        // compute the weight diagonal
        for (int i=0; i<N; i++){
            w[i]=u[i];
        }

        // update b
        double a11=0, a12=0, a21=0, a22=0;
        double e=0, f=0;
        for (int i=0; i<N; i++){
            a11+=w[i];
            a12+=w[i]*X[i];
            a21+=X[i]*w[i];
            a22+=X[i]*w[i]*X[i];
        }
        for (int i=0; i<N; i++){
            e+=w[i]*z[i];
            f+=X[i]*w[i]*z[i];
        }
        double det=a11*a22-a12*a21;
        b[0]=a22*e-a12*f; b[0]/=det;
        b[1]=-a21*e+a11*f; b[1]/=det;

        // X*b
        for (int i=0; i<N; i++){
            Xb[i]=b[0];
            Xb[i]+=X[i]*b[1];
        }

        // z=inv_g(Xb)
        for (int i=0; i<N; i++){
            u[i]=invlinkfunc(Xb[i]);
        }

        // eta
        for (int i=0; i<N; i++){
            eta[i]=Xb[i];
        }

        // res
        res=0;
        for (int i=0; i<N; i++){
            res+=(y[i]-u[i])*(y[i]-u[i]);
        }

        iter++;
    }

    delete eta;
    delete z;
    delete u;
    delete Xb;
    delete w;
}

/**
 * @brief GLM::irls
 * @param X
 * @param y
 * @param N
 * @param b0
 * @param b
 * @param maxIter
 * @param epsilon
 */
void GLM::irls(const double *X, const double *y, int N, const double *b0, double *b, int maxIter, double epsilon){
    int iter=0;
    double  res=1e+100;
    double *eta=new double[N];
    double *z=new double[N];
    double *u=new double[N];
    double *Xb=new double[N];
    double *w=new double[N];

    for (int i=0; i<N; i++){
        eta[i]=linkfunc(y[i]);
    }

    for (int j=0; j<=1; j++){
        b[j]=b0[j];
    }

    while (iter<maxIter && res>epsilon){
        // X*b
        for (int i=0; i<N; i++){
            Xb[i]=b[0];
            Xb[i]+=X[i]*b[1];
        }

        // z=inv_g(Xb)
        for (int i=0; i<N; i++){
            u[i]=invlinkfunc(Xb[i]);
        }

        // compute the pseudo-data zi
        for (int i=0; i<N; i++){
            z[i]=eta[i]+dlinkfunc(u[i])*(y[i]-u[i]);
        }
        // compute the weight diagonal
        for (int i=0; i<N; i++){
            w[i]=u[i];
        }

        // update b
        double a11=0, a12=0, a21=0, a22=0;
        double e=0, f=0;
        for (int i=0; i<N; i++){
            a11+=w[i];
            a12+=w[i]*X[i];
            a21+=X[i]*w[i];
            a22+=X[i]*w[i]*X[i];
        }
        for (int i=0; i<N; i++){
            e+=w[i]*z[i];
            f+=X[i]*w[i]*z[i];
        }
        double det=a11*a22-a12*a21;
        b[0]=a22*e-a12*f; b[0]/=det;
        b[1]=-a21*e+a11*f; b[1]/=det;

        // X*b
        for (int i=0; i<N; i++){
            Xb[i]=b[0];
            Xb[i]+=X[i]*b[1];
        }

        // z=inv_g(Xb)
        for (int i=0; i<N; i++){
            u[i]=invlinkfunc(Xb[i]);
        }

        // eta
        for (int i=0; i<N; i++){
            eta[i]=Xb[i];
        }

        // res
        res=0;
        for (int i=0; i<N; i++){
            res+=(y[i]-u[i])*(y[i]-u[i]);
        }

        iter++;
    }

    delete eta;
    delete z;
    delete u;
    delete Xb;
    delete w;
}

/**
 * @brief GLM::sse
 * @param X
 * @param y
 * @param N
 * @param b
 * @return
 */
double GLM::sse(const double *X, const double *y, int N, const double *b){
    double res=0;
    double Xb[N];
    double u[N];

    // compute X*b
    for (int i=0; i<N; i++){
        Xb[i]=b[0]+b[1]*X[i];
    }

    // compute u
    for (int i=0; i<N; i++){
        u[i]=invlinkfunc(Xb[i]);
    }

    // compute score
    res=0;
    for (int i=0; i<N; i++){
        res+=(y[i]-u[i])*(y[i]-u[i]);
    }

    return res;
}
