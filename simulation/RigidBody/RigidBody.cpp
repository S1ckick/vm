//
//  RigidBody.cpp
//  simulation
//
//  Created by Максим on 06.10.2020.
//  Copyright © 2020 Максим. All rights reserved.
//

#include "RigidBody.hpp"
#include <iostream>
namespace{
    const int STATE_SIZE = 13;
};

void State_to_Array(RigidBody *rb, double *y){

    //position vector at y[0],y[1],y[2]
    *y++ = rb->x.x;
    *y++ = rb->x.y;
    *y++ = rb->x.z;
    
    //quaternion at y[3],y[4],y[5],y[6]
    *y++ = rb->q.r;
    *y++ = rb->q.i;
    *y++ = rb->q.j;
    *y++ = rb->q.k;
    
    //linear momentum at y[7],y[8],y[9]
    *y++ = rb->P.x;
    *y++ = rb->P.y;
    *y++ = rb->P.z;

    //angular momentum at y[10],y[11],y[12]
    *y++ = rb->L.x;
    *y++ = rb->L.y;
    *y++ = rb->L.z;

}

void Array_to_State(RigidBody *rb, double *y){
    //take position vector
    rb->x.x = *y++;
    rb->x.y = *y++;
    rb->x.z = *y++;

    //take quaternion
    rb->q.r = *y++;
    rb->q.i = *y++;
    rb->q.j = *y++;
    rb->q.k = *y++;

    //take linear momentum
    rb->P.x = *y++;
    rb->P.y = *y++;
    rb->P.z = *y++;

    //take angular momentum
    rb->L.x = *y++;
    rb->L.y = *y++;
    rb->L.z = *y++;

    //compute auxiliary variables
    
    //v(t) = P(t)/M   triplet/double
    rb->v.x = rb->P.x/rb->mass;
    rb->v.y = rb->P.y/rb->mass;
    rb->v.z = rb->P.z/rb->mass;

    //produce R from quaternion
    rb->R = quaternion_to_matrix(normalize(rb->q));
    
    //I^-1(t)=R(t)*I^-1body*R(t)^T
    rb->Iinv = multMatrices(multMatrices(rb->R, rb->Ibodyinv),transpose(rb->R));

    //w(t) = I^-1(t)*L(t)
    rb->omega = multMatVec(rb->Iinv,rb->L);
    
}

void Compute_Force_and_Torque(RigidBody *rb){
   //force and torque are constants
}

void ddt_State_to_Array(RigidBody *rb, double *ydot){
    //copy d/dt*x(t) = v(t) into ydot
    *ydot++ = rb->v.x;
    *ydot++ = rb->v.y;
    *ydot++ = rb->v.z;

    //count resultive quaternion as 0.5 * [0, rb->omega]*q
    quaternion rotation_quat;
    rotation_quat.r = 0;
    rotation_quat.i = rb->omega.x;
    rotation_quat.j = rb->omega.y;
    rotation_quat.k = rb->omega.z;
    
    quaternion temp = mult_quat_to_quat(rotation_quat,rb->q);
    
    quaternion res;
    res.r=0.5*temp.r;
    res.i=0.5*temp.i;
    res.j=0.5*temp.j;
    res.k=0.5*temp.k;

    // copy resultive quaternion into ydot
    *ydot++ = res.r;
    *ydot++ = res.i;
    *ydot++ = res.j;
    *ydot++ = res.k;


    // d/dt*P(t) = F(t)
    *ydot++ = rb->force.x;
    *ydot++ = rb->force.y;
    *ydot++ = rb->force.z;

    // d/dt*L(t) = tau(t)
    *ydot++ = rb->torque.x;
    *ydot++ = rb->torque.y;
    *ydot++ = rb->torque.z;
}

void dydt(RigidBody *body,double *y, double *ydot){
    //put data in y[] into Body
    Array_to_State(body,y);
    Compute_Force_and_Torque(body);
    ddt_State_to_Array(body,ydot);
}



void MUL_Y_DOUBLE(double *y, double numb){
    for(int i = 0; i<STATE_SIZE; i++){
        y[i]*=numb;
    }
}

void SUM_Y_DOUBLE(double *y, double numb){
    for(int i = 0; i<STATE_SIZE; i++){
        y[i]+=numb;
    }
}

void SUM_Y_Y(double *y1, double *y2){
    for(int i = 0; i<STATE_SIZE; i++){
        y1[i]+=y2[i];
    }
}

void solveRungeKutta(RigidBody *body, double h){
    double yinit[STATE_SIZE];
    State_to_Array(body,yinit);
    
    //count k_1 = h*f(y)
    double k_1[STATE_SIZE];
    dydt(body,yinit,k_1);
    MUL_Y_DOUBLE(k_1, h);
    double temp[STATE_SIZE];
    for(int i = 0; i<STATE_SIZE; i++){
        temp[i] = k_1[i];
    }
    
    //count k_2 = h*f(y + k_1/2)
    MUL_Y_DOUBLE(temp, 1.0/2.0);
    SUM_Y_Y(temp, yinit);
    double k_2[STATE_SIZE];
    dydt(body,temp,k_2);
    MUL_Y_DOUBLE(k_2, h);
    
    //count k_3 = h*f(y + k_2/2)
    for(int i = 0; i<STATE_SIZE; i++){
        temp[i] = k_2[i];
    }
    MUL_Y_DOUBLE(temp, 1.0/2.0);
    SUM_Y_Y(temp, yinit);
    double k_3[STATE_SIZE];
    dydt(body,temp,k_3);
    MUL_Y_DOUBLE(k_3, h);
    
    //count k_4 = h*f(y + k_3)
    for(int i = 0; i<STATE_SIZE; i++){
        temp[i] = k_3[i];
    }
    SUM_Y_Y(temp, yinit);
    double k_4[STATE_SIZE];
    dydt(body,temp,k_4);
    MUL_Y_DOUBLE(k_4, h);
    
    //count (k_1)/6 , (k_2)/3, (k_3)/3, (k_4)/6
    MUL_Y_DOUBLE(k_1, 1.0/6.0);
    MUL_Y_DOUBLE(k_2, 1.0/3.0);
    MUL_Y_DOUBLE(k_3, 1.0/3.0);
    MUL_Y_DOUBLE(k_4, 1.0/6.0);
    
    double res[STATE_SIZE];
    for(int i = 0; i<STATE_SIZE; i++){
        res[i] = yinit[i];
    }
    //count y_n+1 = y_n + (k_1)/6 + (k_2)/3 + (k_3)/3 + (k_4)/6
    SUM_Y_Y(res, k_1);
    SUM_Y_Y(res, k_2);
    SUM_Y_Y(res, k_3);
    SUM_Y_Y(res, k_4);

    Array_to_State(body,res);
}

void InitStates(RigidBody *body){
    body->mass=2;
    
    //block
    double heightBlock = 2;
    double widthBlock = 2;
    double depthBlock = 4;
    
    matrix Ib_inv;
    //block
    Ib_inv.pos[0][0] = 12/((heightBlock*heightBlock + widthBlock*widthBlock)*body->mass);
    Ib_inv.pos[0][1] = 0;
    Ib_inv.pos[0][2] = 0;
    Ib_inv.pos[1][0] = 0;
    Ib_inv.pos[1][1] = 12/((heightBlock*heightBlock + depthBlock*depthBlock)*body->mass);
    Ib_inv.pos[1][2] = 0;
    Ib_inv.pos[2][0] = 0;
    Ib_inv.pos[2][1] = 0;
    Ib_inv.pos[2][2] = 12/((widthBlock*widthBlock + depthBlock*depthBlock)*body->mass);
    
    body->Ibodyinv = Ib_inv;
    
    body->x.x = 1;
    body->x.y = 1;
    body->x.z = 1;
    
    body->q.r = 0.5;
    body->q.i = 0.5;
    body->q.j = 0.5;
    body->q.k = 0.5;

    body->P.x = 0;
    body->P.y = 0;
    body->P.z = 0;

    body->L.x = 0.5;
    body->L.y = 4;
    body->L.z = 2;

    body->force.x=0;
    body->force.y=0;
    body->force.z=0;

    body->torque.x=0;
    body->torque.y=0;
    body->torque.z=0;
}
