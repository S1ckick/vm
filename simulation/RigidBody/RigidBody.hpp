//
//  RigidBody.hpp
//  simulation
//
//  Created by Максим on 06.10.2020.
//  Copyright © 2020 Максим. All rights reserved.
//

#ifndef RigidBody_hpp
#define RigidBody_hpp
#include "../structsOperations/structsOperations.hpp"

struct RigidBody{
    //Constant quantities
    double mass; //mass M
    matrix Ibody; // I_body
    matrix Ibodyinv; //I^-1_body inverse of I_body

    //state variables
    triple x; //x(t)
    quaternion q; //q(t)
    triple P; //P(t)
    triple L; //L(t);

    //derived quantities
    matrix Iinv; //I^-1(t)
    matrix R;
    triple v; //v(t)
    triple omega; //w(t)

    //computed quantities
    triple force; //F(t)
    triple torque; //tau(t)
};

void State_to_Array(RigidBody *rb, double *y);
void Array_to_State(RigidBody *rb, double *y);
void Compute_Force_and_Torque(double t, RigidBody *rb);
void ddt_State_to_Array(RigidBody *rb, double *ydot);

#endif /* RigidBody_hpp */
