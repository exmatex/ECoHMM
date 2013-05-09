

#include <math.h>
#include <stdlib.h>

class Quaternion
{
public:
	Quaternion() { x = y = z = 0.0; w = 1.0; }
	Quaternion(double x, double y, double z, double w) : x(x), y(y), z(z), w(w) {};
	void set(double ax, double ay, double az, double aw) { x = ax; y = ay; z = az; w = aw; }
	void mul(Quaternion q);
	void setEulerAngles(float pitch, float yaw, float roll);
	void getRotMat(float* m) const;

	double x,y,z,w;
};
