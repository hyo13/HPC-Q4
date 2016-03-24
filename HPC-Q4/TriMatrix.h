//
//  TriMatrix.h
//  HPC-Q4
//
//  Created by hyo13 on 23/3/16.
//  Copyright (c) 2016 hyo13. All rights reserved.
//

#ifndef CLASS_TriMatrix
#define CLASS_TriMatrix

#include <stdio.h>
#include <iostream>
#include <vector>
#include <math.h>
#include <Accelerate/Accelerate.h>
using namespace std;

class TriMatrix{
    
private:
    double *diagm, *diagu, *diagl;
    int s;
    
public:
    TriMatrix(double v,int S){
        diagm = new double [S];
        diagu = new double [S-1];
        diagl = new double [S-1];
        
        //create diagm vector
        diagm[0]=1;
        diagm[S-1]=1;
        for (int i=1;i<S-1;i++){
            diagm[i]=1-2*v;
        }
        
        //create diagu vector
        diagu[0]=0;
        for (int i=1;i<S-1;i++){
            diagu[i]=v;
        }
        
        //create diagl vector
        diagl[S-2]=0;
        for (int i=0;i<S-2;i++){
            diagl[i]=v;
        }
        s=S;
    }
    
    //Operator Overload: Calculate Multiplication of TriMatrix to Vector X
    double *operator* (double *X){
        
        //diagm multiply X
        double Am[s];
        for (int i=0;i<s;i++){
            Am[i]=diagm[i]*X[i];
        }
        
        //diagu multiply X
        double Au[s];
        for (int i=0;i<s-1;i++){
            Au[i]=diagu[i]*X[i+1];
        }
        Au[s-1]=0;
        
        //diagl multiply X
        double Al[s];
        Al[0]=0;
        for (int i=1;i<s;i++){
            Al[i]=diagl[i-1]*X[i-1];
        }
        
        //Superposition of Results
        double *B;
        B=new double [s];
        for (int i=0;i<s;i++){
            B[i]=Am[i]+Au[i]+Al[i];
        }
        return B;
    }
    
    //Operator Overload: Matrix-Vector Solve Operation
    double *operator/ (double *B){
        double M[s];
        double U[s];
        double L[s];
        copy(diagm,diagm+s,M);
        copy(diagu,diagu+s,U);
        copy(diagl,diagl+s,L);
        
        
        //Forward Elimination
        for (int i=1;i<s;i++){
            double m=L[i-1]/M[i-1];
            M[i]=M[i]-m*U[i-1];
            B[i]=B[i]-m*B[i-1];
        }
        //Calculate Xn
        double *X;
        X=new double [s];
        X[s-1]=B[s-1]/M[s-1];
        
        //Backward Substitution
        for (int i=s-2;i>=0;i--){
            X[i]=(B[i]-U[i]*X[i+1])/M[i];
        }
        return X;
    }
    
};




#endif /* defined(__HPC_Q1__TriMatrix__) */
