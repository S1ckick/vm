//
//  RigidBody.cpp
//  simulation
//
//  Created by Максим on 06.10.2020.
//  Copyright © 2020 Максим. All rights reserved.
//

#include "RigidBody.hpp"

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

void Compute_Force_and_Torque(double t, RigidBody *rb){
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
