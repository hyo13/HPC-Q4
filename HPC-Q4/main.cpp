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
    int Nx=10000;
    double T=0.01;
    double Nt=2000;
    double alpha=0.001;
    
    //calculate minimum input time step for Forward Euler to converge with v = or < 0.5
    double dx=L/Nx;
    int Ntmin=2*alpha*T/(pow(dx,2));
    
    //check that Nt is > or = Ntmax to make sure Forward Euler converges
    if (Nt<Ntmin){
        cout<<"ERROR: Input Nt is smaller than minimum time step allowed. Program Terminated."<<endl;
        terminate();
    }
    
    // INITIAL CALCULATIONS
    double dt=T/Nt;
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
    ML.LU();
    TriMatrix MR((1-arg)*v,Nx+1);
    MR.array();
    double *u0new;
    u0new=new double[Nx+1];
    for (int i=0;i<Nt;i++){
        u0new=ML/(MR*u0);
        u0=u0new;
        cout<<i<<endl;
    }
    
    //OUTPUT RESULTS
    for (int i=0;i<Nx+1;i++){
        cout<<u0new[i]<<endl;
    }
    
    return 0;
}
