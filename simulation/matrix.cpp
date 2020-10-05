//
//  matrix.cpp
//  simulation
//
//  Created by Максим on 24.09.2020.
//  Copyright © 2020 Максим. All rights reserved.
//

#include "matrix.hpp"
#include "math.h"

//correct
matrix multMatrices(matrix A, matrix B){
    matrix mul;
    for(int i = 0; i<3; i++){
        for(int j = 0; j<3; j++){
            mul.pos[i][j] = 0;
            for(int k = 0; k<3; k++){
                mul.pos[i][j]+=A.pos[i][k]*B.pos[k][j];
            }
        }
    }
    return mul;
}

//correct
matrix transpose(matrix A){
    matrix transpose;
    for(int i = 0; i<3; i++){
        for(int j = 0; j<3; j++){
            transpose.pos[i][j] = A.pos[j][i];
        }
    }
    return transpose;
}

//correct
triple multMatVec(matrix A, triple v){
    triple res;
    res.x=A.pos[0][0]*v.x + A.pos[0][1]*v.y + A.pos[0][2]*v.z;
    res.y=A.pos[1][0]*v.x + A.pos[1][1]*v.y + A.pos[1][2]*v.z;
    res.z=A.pos[2][0]*v.x + A.pos[2][1]*v.y + A.pos[2][2]*v.z;
    return res;
}

//correct
matrix quaternion_to_matrix(const quaternion & q){
    matrix m;
    m.pos[0][0] = 1-2*q.j*q.j-2*q.k*q.k;
    m.pos[0][1] = 2*q.i*q.j-2*q.r*q.k;
    m.pos[0][2] = 2*q.i*q.k+2*q.r*q.j;
    m.pos[1][0] = 2*q.i*q.j+2*q.r*q.k;
    m.pos[1][1] = 1-2*q.i*q.i-2*q.k*q.k;
    m.pos[1][2] = 2*q.j*q.k-2*q.r*q.i;
    m.pos[2][0] = 2*q.i*q.k-2*q.r*q.j;
    m.pos[2][1] = 2*q.j*q.k+2*q.r*q.i;
    m.pos[2][2] = 1-2*q.i*q.i-2*q.j*q.j;
    return m;
}

//correct
quaternion normalize(const quaternion & q){
    quaternion res;
    double len = sqrt(q.r*q.r+q.i*q.i+q.j*q.j+q.k*q.k);
    res.r = q.r/len;
    res.i = q.i/len;
    res.j = q.j/len;
    res.k=q.k/len;
    return res;
}

//correct
quaternion matrix_to_quaternion(const matrix &m){
    quaternion q;
    double tr,s;
    tr = m.pos[0][0]+ m.pos[1][1] + m.pos[2][2];
    if(tr>=0){
        s = sqrt(tr + 1);
        q.r = 0.5*s;
        s = 0.5/s;
        q.i = (m.pos[2][1] - m.pos[1][2])*s;
        q.j = (m.pos[0][2] - m.pos[2][0])*s;
        q.k = (m.pos[1][0] - m.pos[0][1])*s;
    } else{
        int i = 0;
        if(m.pos[1][1] > m.pos[0][0]){
            i=1;
        }
        if(m.pos[2][2] > m.pos[i][i]) {
            i = 2;
        }
        switch(i){
            case 0 :
                s = sqrt((m.pos[0][0]-(m.pos[1][1]+m.pos[2][2]))+1);
                q.i = 0.5*s;
                s = 0.5/s;
                q.j = (m.pos[0][1]+ m.pos[1][0])*s;
                q.k = (m.pos[2][0]+ m.pos[0][2])*s;
                q.r = (m.pos[2][1] - m.pos[1][2])*s;
                break;
            case 1:
                s = sqrt((m.pos[1][1] -(m.pos[2][2]+m.pos[0][0]))+1);
                q.j = 0.5*s;
                s = 0.5/s;
                q.k = (m.pos[1][2]+ m.pos[2][1])*s;
                q.i = (m.pos[0][1] + m.pos[1][0])*s;
                q.r = (m.pos[0][2] - m.pos[2][0])*s;
                break;
            case 2:
                s = sqrt((m.pos[2][2] - (m.pos[0][0]+m.pos[1][1]))+1);
                q.k = 0.5*s;
                s = 0.5/s;
                q.i = (m.pos[2][0] + m.pos[0][2])*s;
                q.j = (m.pos[1][2] + m.pos[2][1])*s;
                q.k = (m.pos[1][0] - m.pos[0][1])*s;

        }

    }
    return q;
}

//correct
quaternion mult_quat_to_quat(const quaternion & q1, const quaternion & q2){
    quaternion res;
    res.r = q1.r*q2.r - q1.i*q2.i - q1.j*q2.j - q1.k*q2.k;
    
    res.i = q1.r*q2.i + q1.i*q2.r + q1.j*q2.k - q1.k*q2.j;
    res.j = q1.r*q2.j - q1.i*q2.k + q1.j*q2.r + q1.k*q2.i;
    res.k = q1.r*q2.k + q1.i*q2.j - q1.j*q2.i + q1.k*q2.r;
    return res;
}


