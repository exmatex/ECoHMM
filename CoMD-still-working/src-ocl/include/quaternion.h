#ifndef QUATERNION_H
#define QUATERNION_H

typedef struct quaternion { double x,y,z,w; } Quaternion;

extern void QuaternionMul(Quaternion* qr, Quaternion q1, Quaternion q2);

extern void QuaternionSetEulerAngles(Quaternion* qr, float pitch, float yaw, float roll);

extern void QuaternionGetRotMat(float* m, Quaternion q1);

#endif

