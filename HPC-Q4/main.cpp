//
//  main.cpp
//  HPC-Q4
//
//  Created by hyo13 on 23/3/16.
//  Copyright (c) 2016 hyo13. All rights reserved.
//

#include <iostream>
#include <vector>
#include <math.h>
#include "TriMatrix.h"
#include <Accelerate/Accelerate.h>
using namespace std;

int main() {
    
    // INPUTS
    double L=1;
    int Nx=20;
    double T=5;
    double Nt=5000;
    double alpha=1;
    
    // INITIAL CALCULATIONS
    double dt=T/Nt;
    double dx=L/Nx;
    double v=alpha*dt/pow(dx,2);
    
    // X-POSIION VECTOR
    double *x;
    x=new double[Nx+1];
    for (int i=0;i<Nx+1;i++){
        x[i]=0+i*dx;
    }
    
    //INITIAL CONDITION
    double *u0;
    u0=new double[Nx+1];
    u0[1]=0;
    u0[Nx+1]=0;
    for (int i=1;i<Nx;i++){
        u0[i]=x[i]*(1-x[i]);
    }
    
    //IMPLICIT TIME INTEGRATION
    double arg=0.5;
    TriMatrix ML(-arg*v,Nx+1);
    TriMatrix MR((1-arg)*v,Nx+1);
    double *u0new;
    u0new=new double[Nx+1];
    for (int i=0;i<Nt;i++){
        u0new=ML/(MR*u0);
        u0=u0new;
    }
    
    //OUTPUT RESULTS
    for (int i=0;i<Nx+1;i++){
        cout<<u0new[i]<<endl;
    }
    
    return 0;
}
