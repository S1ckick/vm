//
//  structsOperations.hpp
//  simulation
//
//  Created by Максим on 24.09.2020.
//  Copyright © 2020 Максим. All rights reserved.
//

#ifndef structsOperations_hpp
#define structsOperations_hpp

struct quaternion{
    double r;
    double i;
    double j;
    double k;
    quaternion(){
        r = 0;
        i = 0;
        j = 0;
        k = 0;
    }
};

struct triple{
    double x;
    double y;
    double z;
    triple(){
        x = 0;
        y = 0;
        z = 0;
    }
};

struct matrix{
    double pos[3][3];

    matrix& operator=(matrix A){
        for(int i = 0; i<3; i++){
            for(int j = 0; j<3; j++){
                this->pos[i][j] = A.pos[i][j];
            }
        }
        return *this;
    }
};

matrix multMatrices(matrix A, matrix B);

matrix transpose(matrix A);
triple multMatVec(matrix A, triple v);
matrix quaternion_to_matrix(const quaternion & q);
quaternion normalize(const quaternion & q);
quaternion matrix_to_quaternion(const matrix &m);
quaternion mult_quat_to_quat(const quaternion & q1, const quaternion & q2);

#endif /* structsOperations_hpp */
