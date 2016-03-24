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
    double *diagm,*diagu,*diagl,*A,*M,*U,*L,*U2;
    int s, *ipiv;
    
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

    //Function: Convert all column vectors of TriMatrix into array A for multiplitcation
    void array(){
        A=new double [s*s];
        //First column vector
        A[0]=diagm[0];
        A[1]=diagl[0];
        //Intermediate column vectors
        for (int i=1;i<s-1;i++){
            A[s*i+i-1]=diagu[i-1];
            A[s*i+i]=diagm[i];
            A[s*i+i+1]=diagl[i];
        }
        //Last column vector
        A[s*s-2]=diagu[s-2];
        A[s*s-1]=diagm[s-1];
    }
    
    //Operator Overload: Calculate Multiplication of TriMatrix to Vector X
    double *operator* (double *X){
        double *B;
        B = new double[s];
        double alpha=1;
        double beta=0;
        cblas_dgemv(CblasColMajor,CblasNoTrans,s,s,alpha,A,s,X,1,beta,B,1);
        return B;
    }
    
    //Function: LU Decomposition of TriMatrix for solve operation
    void LU(){
        M=new double [s];
        U=new double [s];
        L=new double [s];
        copy(diagm,diagm+s,M);
        copy(diagu,diagu+s,U);
        copy(diagl,diagl+s,L);
        ipiv=new int [s];
        U2=new double [s-2];
        int info;
        //LU Decomposition of TriMatrix
        dgttrf_(&s,L,M,U,U2,ipiv,&info);

    }
    
    //Operator Overload: Matrix-Vector Solve Operation
    double *operator/ (double *B){
        //Solve System
        int nrhs = 1;
        char trans = 'N';
        int info;
        dgttrs_(&trans,&s,&nrhs,L,M,U,U2,ipiv,B,&s,&info);
        return B;
    }
    
};




#endif /* defined(__HPC_Q1__TriMatrix__) */
