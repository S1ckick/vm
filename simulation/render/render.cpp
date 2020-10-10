//
//  render.cpp
//  simulation
//
//  Created by Максим on 10.10.2020.
//  Copyright © 2020 Максим. All rights reserved.
//

#include "render.hpp"
#include <iostream>
#include <iomanip>
#include "../structsOperations/structsOperations.hpp"
#include <math.h>
#include <ctime>
double DeltaTime;
RigidBody Body;
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


void Reshape(int W, int H) {
  glViewport(0, 0, W*2, H*2);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(30.0, ((float)W)/((float)H), 1.0, 1000.0);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  gluLookAt(0, 30, 0, 0, 0, 0, 0, 0, 1);
}


void Display(){
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    solveRungeKutta(&Body,DeltaTime);
    
    glPushMatrix();
    glTranslated(Body.x.x, Body.x.y, Body.x.z);
    glRotatef(2*acos(Body.q.r)*180/M_PI,  Body.q.i, Body.q.j, Body.q.k);
    glPolygonMode(GL_FRONT, GL_FILL);
    drawBlock();
    glPopMatrix();
    glFlush();
    glutSwapBuffers();
    
    //logs
    std::cout << std::setprecision(20);
    std::cout << "quaternion: " << Body.q.r << " " << Body.q.i << " " << Body.q.j << " "<<Body.q.k << "\n";
    triple invariant = multMatVec(Body.Iinv, Body.omega);
    std::cout << "invariant: " << invariant.x*Body.omega.x + invariant.y*Body.omega.y + invariant.z*Body.omega.z << "\n\n";
}

void Idle() {
  long Time;
  static long OldTime = -1;
  static long StartTime;

  if (OldTime == -1)
    StartTime = OldTime = clock();
  else {
    Time = clock();
    DeltaTime = (double)(Time - OldTime) / CLOCKS_PER_SEC;
    OldTime = clock();
  }
  glutPostRedisplay();
}

void deploy(){
    
    InitStates(&Body);
    
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
    glutInitWindowSize(1024, 620);
    glutInitWindowPosition(0, 0);
    glutCreateWindow("rotation");
        
    glutReshapeFunc(Reshape);
    glutDisplayFunc(Display);
    glutIdleFunc(Idle);
    
    glClearColor(0.45, 0.72, 0.65, 0);
    glEnable(GL_DEPTH_TEST);
    glShadeModel(GL_SMOOTH);
    
    glutMainLoop();
}
