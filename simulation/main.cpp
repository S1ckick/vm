//
//  main.cpp
//  simulation
//
//  Created by Максим on 21.09.2020.
//  Copyright © 2020 Максим. All rights reserved.
//

#define GL_SILENCE_DEPRECATION
#include <iostream>
#include <GL/glew.h>
#include <GLUT/glsmap.h>
#include <GLUT/GLUT.h>
#include <GLFW/glfw3.h>

#include "matrix.hpp"

#include <cstdio>
#include <math.h>
using namespace std;

namespace{
    const int NBODIES = 1;
    const int STATE_SIZE = 13;
};

GLFWwindow* initWindow(const int resX, const int resY);
void controls(GLFWwindow* window, int key, int scancode, int action, int mods);
void drawCube();
matrix Star(triple a);
void InitStates();
void drawCone();
void drawCylinder();
void drawBlock();

void rotationMatrix(quaternion &q, GLfloat *mat){
    mat[0] = 1-2*q.j*q.j-2*q.k*q.k;
    mat[1] = 2*q.i*q.j-2*q.r*q.k;
    mat[2] = 2*q.i*q.k+2*q.r*q.j;
    
    mat[3] = 2*q.i*q.j+2*q.r*q.k;
    mat[4] = 1-2*q.i*q.i-2*q.k*q.k;
    mat[5] = 2*q.j*q.k-2*q.r*q.i;
    
    mat[6] = 2*q.i*q.k-2*q.r*q.j;
    mat[7] = 2*q.j*q.k+2*q.r*q.i;
    mat[8] = 1-2*q.i*q.i-2*q.j*q.j;
}

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

RigidBody Bodies[NBODIES];

void State_to_Array(RigidBody *rb, double *y){

    *y++ = rb->x.x;  //x component at position 0
    *y++ = rb->x.y;  //x component at position 1
    *y++ = rb->x.z; //x component at position 2

    *y++ = rb->q.r;
    *y++ = rb->q.i;
    *y++ = rb->q.j;
    *y++ = rb->q.k;

    *y++ = rb->P.x;
    *y++ = rb->P.y;
    *y++ = rb->P.z;

    *y++ = rb->L.x;
    *y++ = rb->L.y;
    *y++ = rb->L.z;

}

void Array_to_State(RigidBody *rb, double *y){
    rb->x.x = *y++;
    rb->x.y = *y++;
    rb->x.z = *y++;

    rb->q.r = *y++;
    rb->q.i = *y++;
    rb->q.j = *y++;
    rb->q.k = *y++;

    rb->P.x = *y++;
    rb->P.y = *y++;
    rb->P.z = *y++;

    rb->L.x = *y++;
    rb->L.y = *y++;
    rb->L.z = *y++;

    //compute auxiliary variables
    //v(t) = P(t)/M   triplet/double
    rb->v.x = rb->P.x/rb->mass;
    rb->v.y = rb->P.y/rb->mass;
    rb->v.z = rb->P.z/rb->mass;

    rb->R = quaternion_to_matrix(normalize(rb->q));
    
    //I^-1(t)=R(t)*I^-1body*R(t)^T
    rb->Iinv = multMatrices(multMatrices(rb->R, rb->Ibodyinv),transpose(rb->R));

    //w(t) = I^-1(t)*L(t)
    rb->omega = multMatVec(rb->Iinv,rb->L);
    
}

void Array_to_Bodies(double y[]){
    for(int i = 0; i< NBODIES; i++){
        Array_to_State(&Bodies[i],&y[i*STATE_SIZE]);
    }
}

void Bodies_to_Array(double y[]){
    for(int i = 0; i<NBODIES;i++){
        State_to_Array(&Bodies[i], &y[i*STATE_SIZE]);
    }
}

void Compute_Force_and_Torque(double t, RigidBody *rb){
   
    
    
}

void ddt_State_to_Array(RigidBody *rb, double *ydot){
    //copy d/dt*x(t) = v(t) into ydot
    *ydot++ = rb->v.x;
    *ydot++ = rb->v.y;
    *ydot++ = rb->v.z;

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

void dydt(double t, double y[], double ydot[]){
    //put data in y[] into Bodies[]
    Array_to_Bodies(y);
    for(int i = 0; i<NBODIES; i++){
        Compute_Force_and_Torque(t,&Bodies[i]);
        ddt_State_to_Array(&Bodies[i],&ydot[i*STATE_SIZE]);
    }
}

void display( GLFWwindow* window )
{
    /*
    GLfloat front_color[] = {0, 1, 0, 1};
    GLfloat back_color[] = {0, 0, 1, 1};
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, front_color);
    glMaterialfv(GL_BACK, GL_DIFFUSE, back_color);
    */
    
    double y0[STATE_SIZE*NBODIES];
    
    InitStates();
    double h=1./10.;
    for(double t = 0; t<60.0;t+=h){
        Bodies_to_Array(y0);
        double ydot1[STATE_SIZE*NBODIES];
        
        dydt(t/1000,y0,ydot1);
        
        for(int i = 0; i<13; i++){
            ydot1[i]*=h/2;
            ydot1[i]+=y0[i];
        }
        double ydot2[STATE_SIZE*NBODIES];
        dydt(t/1000,ydot1,ydot2);
        for(int i = 0; i<13; i++){
            ydot2[i]*=h;
            ydot2[i]+=y0[i];
        }
        
        quaternion q;
        q.r = ydot2[3];
        q.i = ydot2[4];
        q.j = ydot2[5];
        q.k = ydot2[6];

        GLfloat centerPosX = ydot2[0];
        GLfloat centerPosY = ydot2[1];
        GLfloat centerPosZ = ydot2[2];
        
        // Scale to window size
        GLint windowWidth, windowHeight;
        glfwGetFramebufferSize(window, &windowWidth, &windowHeight);
        glViewport(0, 0, windowWidth, windowHeight);
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        gluPerspective(30.0, ((float)windowWidth)/((float)windowHeight), 1.0, 1000.0);
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        gluLookAt(0, 50, 0, 0, 0, 0, 0, 0, 1);
        
        // Draw stuff
        glClearColor(0.0, 1.0, 0.0, 0.0);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        glColor3f (1.0, 1.0, 1.0);
        //glLoadIdentity();
        //glTranslatef(centerPosX-2,centerPosY,centerPosZ);
        
        glRotatef(2*acos(q.r)*180/M_PI, q.i, q.j, q.k);
        glPolygonMode(GL_FRONT, GL_FILL);
        
        glRotatef(90, 0, 1, 0);
        
        drawBlock();

        // Update Screen
        glfwSwapBuffers(window);

        // Check for any input, or window movement
        glfwPollEvents();
        
        
        
        //copy d/dt*Y(t+1/30) into state variables
        Array_to_Bodies(ydot2);
        
        double yShow[STATE_SIZE*NBODIES];
        Bodies_to_Array(yShow);
        cout << "quaternion" << yShow[3] << " " << yShow[4]  << " " << yShow[5]  << " " << yShow[6]  << "\n";
        cout << "pos" << yShow[0] << " " << yShow[1] << " " << yShow[2] << "\n";
        cout << "P: " << yShow[7] << " " << yShow[8] << " " << yShow[9] << "\n";
        cout << "L: " << yShow[10] << " " << yShow[11] << " " << yShow[12] << "\n";
        triple invariant= multMatVec(Bodies->Iinv, Bodies->omega);
        cout << "invariant" << invariant.x*Bodies->omega.x + invariant.y*Bodies->omega.y +invariant.z*Bodies->omega.z << " " << invariant.z*Bodies->omega.z<<  "\n" ;

    }

}



int main(int argc, char** argv)
{
    GLFWwindow* window = initWindow(1024, 620);
    if( NULL != window )
    {
        
        display( window );
    }
    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}


void InitStates(){
    Bodies->mass=2;
    //cube constant
    double sideLength = 2;
    //cone constants
    double coneRadius = 1;
    double coneHeight = 2;
    //cylinder constants
    double cylinderHeight = 2;
    double cylinderRadius = 1;
    //block
    double heightBlock = 2;
    double widthBlock = 2;
    double depthBlock = 4;
    
    matrix Ib;
    //cube
    
    Ib.pos[0][0] = Bodies->mass*sideLength*sideLength/6;
    Ib.pos[0][1] = 0;
    Ib.pos[0][2] = 0;
    Ib.pos[1][0] = 0;
    Ib.pos[1][1] = Bodies->mass*sideLength*sideLength/6;
    Ib.pos[1][2] = 0;
    Ib.pos[2][0] = 0;
    Ib.pos[2][1] = 0;
    Ib.pos[2][2] = Bodies->mass*sideLength*sideLength/6;
     
    //cone
    /*
    Ib.pos[0][0] = 3.0/5.0*Bodies->mass*(coneRadius*coneRadius/4+coneHeight*coneHeight);
    Ib.pos[0][1] = 0;
    Ib.pos[0][2] = 0;
    Ib.pos[1][0] = 0;
    Ib.pos[1][1] = 3.0/5.0*Bodies->mass*(coneRadius*coneRadius/4+coneHeight*coneHeight);
    Ib.pos[1][2] = 0;
    Ib.pos[2][0] = 0;
    Ib.pos[2][1] = 0;
    Ib.pos[2][2] = 3.0/10.0*Bodies->mass*coneRadius*coneRadius;
    */
    
    //cylinder
    /*
    Ib.pos[0][0] = 1.0/12.0*Bodies->mass*(3.0*cylinderRadius*cylinderRadius+cylinderHeight*cylinderHeight);
    Ib.pos[0][1] = 0;
    Ib.pos[0][2] = 0;
    Ib.pos[1][0] = 0;
    Ib.pos[1][1] = 1.0/12.0*Bodies->mass*(3.0*cylinderRadius*cylinderRadius+cylinderHeight*cylinderHeight);
    Ib.pos[1][2] = 0;
    Ib.pos[2][0] = 0;
    Ib.pos[2][1] = 0;
    Ib.pos[2][2] = Bodies->mass*cylinderRadius*cylinderRadius/2.0;
    */
    Bodies->Ibody=Ib;
    
    
    matrix Ib_inv;
    //cube
    /*
    Ib_inv.pos[0][0] = 6/(sideLength*sideLength*Bodies->mass);
    Ib_inv.pos[0][1] = 0;
    Ib_inv.pos[0][2] = 0;
    Ib_inv.pos[1][0] = 0;
    Ib_inv.pos[1][1] = 6/(sideLength*sideLength*Bodies->mass);
    Ib_inv.pos[1][2] = 0;
    Ib_inv.pos[2][0] = 0;
    Ib_inv.pos[2][1] = 0;
    Ib_inv.pos[2][2] = 6/(sideLength*sideLength*Bodies->mass);
    */
    //block
    
    Ib_inv.pos[0][0] = 12/((heightBlock*heightBlock + widthBlock*widthBlock)*Bodies->mass);
    Ib_inv.pos[0][1] = 0;
    Ib_inv.pos[0][2] = 0;
    Ib_inv.pos[1][0] = 0;
    Ib_inv.pos[1][1] = 6/((heightBlock*heightBlock + depthBlock*depthBlock)*sideLength*sideLength*Bodies->mass);
    Ib_inv.pos[1][2] = 0;
    Ib_inv.pos[2][0] = 0;
    Ib_inv.pos[2][1] = 0;
    Ib_inv.pos[2][2] = 6/((widthBlock*widthBlock + depthBlock*depthBlock)*sideLength*sideLength*Bodies->mass);
    
    
    //cone
    /*
    Ib_inv.pos[0][0] = 5.0/(3.0*Bodies->mass*(coneRadius*coneRadius/4+coneHeight*coneHeight));
    Ib_inv.pos[0][1] = 0;
    Ib_inv.pos[0][2] = 0;
    Ib_inv.pos[1][0] = 0;
    Ib_inv.pos[1][1] = 5.0/(3.0*Bodies->mass*(coneRadius*coneRadius/4+coneHeight*coneHeight));
    Ib_inv.pos[1][2] = 0;
    Ib_inv.pos[2][0] = 0;
    Ib_inv.pos[2][1] = 0;
    Ib_inv.pos[2][2] = 10.0/(3.0*Bodies->mass*coneRadius*coneRadius);
    */
    
    //cylinder
    /*
    Ib_inv.pos[0][0] = 12.0/(Bodies->mass*(3.0*cylinderRadius*cylinderRadius+cylinderHeight*cylinderHeight));
    Ib_inv.pos[0][1] = 0;
    Ib_inv.pos[0][2] = 0;
    Ib_inv.pos[1][0] = 0;
    Ib_inv.pos[1][1] = 12.0/(Bodies->mass*(3.0*cylinderRadius*cylinderRadius+cylinderHeight*cylinderHeight));
    Ib_inv.pos[1][2] = 0;
    Ib_inv.pos[2][0] = 0;
    Ib_inv.pos[2][1] = 0;
    Ib_inv.pos[2][2] = 2.0/(Bodies->mass*cylinderRadius*cylinderRadius);
    */
    
    Bodies->Ibodyinv = Ib_inv;
    Bodies->x.x = 0;
    Bodies->x.y = 0;
    Bodies->x.z = 0;
    
    Bodies->q.r = 0.5;
    Bodies->q.i = 0.5;
    Bodies->q.j = 0.5;
    Bodies->q.k = 0.5;

    Bodies->P.x = 0;
    Bodies->P.y = 0;
    Bodies->P.z = 0;

    Bodies->L.x = 1;
    Bodies->L.y = 1;
    Bodies->L.z = 1;

    Bodies->force.x=0;
    Bodies->force.y=0;
    Bodies->force.z=0;

    Bodies->torque.x=0;
    Bodies->torque.y=0;
    Bodies->torque.z=0;
}


GLFWwindow* initWindow(const int resX, const int resY)
{
    
    
    if(!glfwInit())
    {
        fprintf(stderr, "Failed to initialize GLFW\n");
        return NULL;
    }
    glfwWindowHint(GLFW_SAMPLES, 4); // 4x antialiasing

    // Open a window and create its OpenGL context
    GLFWwindow* window = glfwCreateWindow(resX, resY, "TEST", NULL, NULL);

    if(window == NULL)
    {
        fprintf(stderr, "Failed to open GLFW window.\n");
        glfwTerminate();
        return NULL;
    }

    glfwMakeContextCurrent(window);
    glfwSetKeyCallback(window, controls);

    // Get info of GPU and supported OpenGL version
    printf("Renderer: %s\n", glGetString(GL_RENDERER));
    printf("OpenGL version supported %s\n", glGetString(GL_VERSION));
    

    

    
    

    glEnable(GL_DEPTH_TEST); // Depth Testing
    glDepthFunc(GL_LEQUAL);
    glDisable(GL_CULL_FACE);
    glCullFace(GL_BACK);
    
    /*
    float pos[4] = {0, 3, 0.5, 1.5};
    float dir[3] = {-1, -1, -1};
    GLfloat mat_specular[] = {1, 1, 1, 1};
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glLightfv(GL_LIGHT0, GL_POSITION, pos);
    glLightfv(GL_LIGHT0, GL_SPOT_DIRECTION, dir);
    glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
    glMaterialf(GL_FRONT, GL_SHININESS, 128.0);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
    */
    return window;
}


void drawCube()
{
    glColor3f(1,0,0); //red
    glBegin (GL_QUADS);
        glVertex3f (1.0, 1.0, 1.0);
        glVertex3f (-1.0, 1.0, 1.0);
        glVertex3f (-1.0, -1.0, 1.0);
        glVertex3f (1.0, -1.0, 1.0);
    glEnd();
    glColor3f(0,1,0); //green
    glBegin (GL_QUADS);
        glVertex3f (1.0, 1.0, -1.0);
        glVertex3f (1.0, -1.0, -1.0);
        glVertex3f (-1.0, -1.0, -1.0);
        glVertex3f (-1.0, 1.0, -1.0);
    glEnd();
    glColor3f(0,0,1);
    glBegin (GL_QUADS);
        glVertex3f (-1.0, 1.0, 1.0);
        glVertex3f (-1.0, 1.0, -1.0);
        glVertex3f (-1.0, -1.0, -1.0);
        glVertex3f (-1.0, -1.0, 1.0);
    glEnd();
    glColor3f(1,1,0);
    glBegin (GL_QUADS);
        glVertex3f (1.0, 1.0, 1.0);
        glVertex3f (1.0, -1.0, 1.0);
        glVertex3f (1.0, -1.0, -1.0);
        glVertex3f (1.0, 1.0, -1.0);
    glEnd();
    glColor3f(1,0,1);
    glBegin (GL_QUADS);
        glVertex3f (-1.0, 1.0, -1.0);
        glVertex3f (-1.0, 1.0, 1.0);
        glVertex3f (1.0, 1.0, 1.0);
        glVertex3f (1.0, 1.0, -1.0);
    glEnd();
    glColor3f(1,1,1);
    glBegin(GL_QUADS);
        glVertex3f (-1.0, -1.0, -1.0);
        glVertex3f (1.0, -1.0, -1.0);
        glVertex3f (1.0, -1.0, 1.0);
        glVertex3f (-1.0, -1.0, 1.0);
    glEnd();
}

void drawBlock()
{
    glColor3f(1,0,0); //red
    glBegin (GL_QUADS);
        glVertex3f (2.0, 1.0, 1.0);
        glVertex3f (-2.0, 1.0, 1.0);
        glVertex3f (-2.0, -1.0, 1.0);
        glVertex3f (2.0, -1.0, 1.0);
    glEnd();
    glColor3f(0,1,0); //green
    glBegin (GL_QUADS);
        glVertex3f (2.0, 1.0, -1.0);
        glVertex3f (2.0, -1.0, -1.0);
        glVertex3f (-2.0, -1.0, -1.0);
        glVertex3f (-2.0, 1.0, -1.0);
    glEnd();
    glColor3f(0,0,1);
    glBegin (GL_QUADS);
        glVertex3f (-2.0, 1.0, 1.0);
        glVertex3f (-2.0, 1.0, -1.0);
        glVertex3f (-2.0, -1.0, -1.0);
        glVertex3f (-2.0, -1.0, 1.0);
    glEnd();
    glColor3f(1,1,0);
    glBegin (GL_QUADS);
        glVertex3f (2.0, 1.0, 1.0);
        glVertex3f (2.0, -1.0, 1.0);
        glVertex3f (2.0, -1.0, -1.0);
        glVertex3f (2.0, 1.0, -1.0);
    glEnd();
    glColor3f(1,0,1);
    glBegin (GL_QUADS);
        glVertex3f (-2.0, 1.0, -1.0);
        glVertex3f (-2.0, 1.0, 1.0);
        glVertex3f (2.0, 1.0, 1.0);
        glVertex3f (2.0, 1.0, -1.0);
    glEnd();
    glColor3f(1,1,1);
    glBegin(GL_QUADS);
        glVertex3f (-2.0, -1.0, -1.0);
        glVertex3f (2.0, -1.0, -1.0);
        glVertex3f (2.0, -1.0, 1.0);
        glVertex3f (-2.0, -1.0, 1.0);
    glEnd();
}


void drawCone(){
    glutSolidCone(1, 2, 100, 100);
    glTranslatef(0, 0, 1/2);
}

void drawCylinder(){
    GLUquadricObj *quadratic;
    quadratic = gluNewQuadric();
    gluCylinder(quadratic, 1, 1, 2, 100, 100);
}

matrix Star(triple a){
    matrix star;
    star.pos[0][0] = 0;
    star.pos[0][1] = -(a.z);
    star.pos[0][2] = a.y;
    star.pos[1][0] = a.z;
    star.pos[1][1] = 0;
    star.pos[1][2] = -(a.x);
    star.pos[2][0] = -(a.y);
    star.pos[2][1] = a.x;
    star.pos[2][2] =0;
    return star;
}

void controls(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    if(action == GLFW_PRESS)
        if(key == GLFW_KEY_ESCAPE)
            glfwSetWindowShouldClose(window, GL_TRUE);
}
