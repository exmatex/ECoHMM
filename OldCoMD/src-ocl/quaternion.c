
#include <math.h>

#include "quaternion.h"

void QuaternionMul(Quaternion* qr, Quaternion q1, Quaternion q2)
{
   float tx, ty, tz, tw;
   tx = q1.w*q2.x + q1.x*q2.w + q1.y*q2.z - q1.z*q2.y;
   ty = q1.w*q2.y + q1.y*q2.w + q1.z*q2.x - q1.x*q2.z;
   tz = q1.w*q2.z + q1.z*q2.w + q1.x*q2.y - q1.y*q2.x;
   tw = q1.w*q2.w - q1.x*q2.x - q1.y*q2.y - q1.z*q2.z;

   qr->x = tx; qr->y = ty; qr->z = tz; qr->w = tw;
}


void QuaternionSetEulerAngles(Quaternion* qr, float pitch, float yaw, float roll)
{
   qr->w = cos(pitch/2.0)*cos(yaw/2.0)*cos(roll/2.0) - sin(pitch/2.0)*sin(yaw/2.0)*sin(roll/2.0);
   qr->x = sin(pitch/2.0)*sin(yaw/2.0)*cos(roll/2.0) + cos(pitch/2.0)*cos(yaw/2.0)*sin(roll/2.0);
   qr->y = sin(pitch/2.0)*cos(yaw/2.0)*cos(roll/2.0) + cos(pitch/2.0)*sin(yaw/2.0)*sin(roll/2.0);
   qr->z = cos(pitch/2.0)*sin(yaw/2.0)*cos(roll/2.0) - sin(pitch/2.0)*cos(yaw/2.0)*sin(roll/2.0);

   float norm = sqrt(qr->x*qr->x + qr->y*qr->y + qr->z*qr->z + qr->w*qr->w);
   if (norm > 0.00001) 
   { 
      qr->x /= norm;  
      qr->y /= norm;  
      qr->z /= norm;  
      qr->w /= norm; 
   }
}


void QuaternionGetRotMat(float* m, Quaternion q1) 
{
   unsigned int i;
   float x = q1.x;  float y = q1.y;  float z = q1.z;  float w = q1.w;
   for (i=0; i<16; i++) m[i] = 0.0; m[15] = 1.0;
   m[0]  = 1 - 2*y*y - 2*z*z;  m[1]  = 2*x*y - 2*z*w;      m[2]  = 2*x*z + 2*y*w;
   m[4]  = 2*x*y + 2*z*w;      m[5]  = 1 - 2*x*x - 2*z*z;  m[6]  = 2*y*z - 2*x*w;
   m[8]  = 2*x*z - 2*y*w;      m[9]  = 2*y*z + 2*x*w;      m[10] = 1 - 2*x*x - 2*y*y;
}








