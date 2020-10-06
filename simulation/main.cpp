//
//  main.cpp
//  simulation
//
//  Created by Максим on 21.09.2020.
//  Copyright © 2020 Максим. All rights reserved.
//

#include <iostream>
#include <iomanip>
#include <GL/glew.h>
#include <GLUT/GLUT.h>
#include <GLFW/glfw3.h>

#include "structsOperations/structsOperations.hpp"
#include "RigidBody/RigidBody.hpp"

#include <math.h>
using namespace std;

namespace{
    const int NBODIES = 1;
    const int STATE_SIZE = 13;
};

GLFWwindow* initWindow(const int resX, const int resY);
void controls(GLFWwindow* window, int key, int scancode, int action, int mods);
void drawBlock();


RigidBody Bodies[NBODIES];

void InitStates(){
    
    Bodies->mass=2;
    
    //block
    double heightBlock = 2;
    double widthBlock = 2;
    double depthBlock = 4;
    
    matrix Ib;
    //block
    Ib.pos[0][0] = ((heightBlock*heightBlock + widthBlock*widthBlock)*Bodies->mass)/12;
    Ib.pos[0][1] = 0;
    Ib.pos[0][2] = 0;
    Ib.pos[1][0] = 0;
    Ib.pos[1][1] = ((heightBlock*heightBlock + depthBlock*depthBlock)*Bodies->mass)/12;
    Ib.pos[1][2] = 0;
    Ib.pos[2][0] = 0;
    Ib.pos[2][1] = 0;
    Ib.pos[2][2] = ((widthBlock*widthBlock + depthBlock*depthBlock)*Bodies->mass)/12;
    
    Bodies->Ibody=Ib;
    
    matrix Ib_inv;
    //block
    Ib_inv.pos[0][0] = 12/((heightBlock*heightBlock + widthBlock*widthBlock)*Bodies->mass);
    Ib_inv.pos[0][1] = 0;
    Ib_inv.pos[0][2] = 0;
    Ib_inv.pos[1][0] = 0;
    Ib_inv.pos[1][1] = 12/((heightBlock*heightBlock + depthBlock*depthBlock)*Bodies->mass);
    Ib_inv.pos[1][2] = 0;
    Ib_inv.pos[2][0] = 0;
    Ib_inv.pos[2][1] = 0;
    Ib_inv.pos[2][2] = 12/((widthBlock*widthBlock + depthBlock*depthBlock)*Bodies->mass);
    
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

    Bodies->L.x = 0.5;
    Bodies->L.y = 1;
    Bodies->L.z = 0.5;

    Bodies->force.x=0;
    Bodies->force.y=0;
    Bodies->force.z=0;

    Bodies->torque.x=0;
    Bodies->torque.y=0;
    Bodies->torque.z=0;
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
    double y0[STATE_SIZE*NBODIES];
    InitStates();
    
    //middle point method
    //x(t + h) = x(t) + h*f( x(t) + h/2*f( x(t) ) )
    double h=1./10.;
    for(double t = 0; t<60.0;t+=h){
        //push x(t) into y0
        Bodies_to_Array(y0);
        double ydot1[STATE_SIZE*NBODIES];
        
        //compute and push f( x(t) ) into ydot1
        dydt(t/1000,y0,ydot1);
        
        //multiply f(x(t)) by h/2 and add x(t)
        for(int i = 0; i<13; i++){
            ydot1[i]*=h/2;
            ydot1[i]+=y0[i];
        }
        
        double ydot2[STATE_SIZE*NBODIES];
        
        //compute f( x(t) + h/2*f(x(t))) and push it into ydot2
        dydt(t/1000,ydot1,ydot2);
        
        //multiply ydot2 by h and add x(t)
        for(int i = 0; i<13; i++){
            ydot2[i]*=h;
            ydot2[i]+=y0[i];
        }
        //ydot2 contains Y(t + h)
        
        //copy d/dt*Y(t+h) into state variables
        Array_to_Bodies(ydot2);
        
        //form quaternion to rotate an object
        quaternion q;
        q.r = ydot2[3];
        q.i = ydot2[4];
        q.j = ydot2[5];
        q.k = ydot2[6];

        //form position vector to translate an object
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
        glClearColor(0.0, 1.0, 0.0, 0.0);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        //translate an object
        glTranslatef(centerPosX,centerPosY,centerPosZ);
        
        //rotate an object
        glRotatef(2*acos(q.r)*180/M_PI, q.i, q.j, q.k);
        
        //draw an object
        glPolygonMode(GL_FRONT, GL_FILL);
        drawBlock();

        // Update Screen
        glfwSwapBuffers(window);

        // Check for any input, or window movement
        glfwPollEvents();
        
        // print logs
        double yShow[STATE_SIZE*NBODIES];
        Bodies_to_Array(yShow);
        cout << setprecision(20);
        cout << "quaternions: \n" << yShow[3] << " \n" << yShow[4]  << " \n" << yShow[5]  << " \n" << yShow[6]  << "\n";
        cout << "pos: \n" << yShow[0] << "\n" << yShow[1] << "\n" << yShow[2] << "\n";
        cout << "P: \n" << yShow[7] << "\n" << yShow[8] << "\n" << yShow[9] << "\n";
        cout << "L: \n" << yShow[10] << "\n" << yShow[11] << "\n" << yShow[12] << "\n";
        triple invariant= multMatVec(Bodies->Iinv, Bodies->omega);
        cout <<"invariant: " << invariant.x*Bodies->omega.x + invariant.y*Bodies->omega.y +invariant.z*Bodies->omega.z << "\n";
    }
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

    return window;
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
    glColor3f(0,0,0); //green
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

void controls(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    if(action == GLFW_PRESS)
        if(key == GLFW_KEY_ESCAPE)
            glfwSetWindowShouldClose(window, GL_TRUE);
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
