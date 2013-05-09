#ifndef _ATOMVISUALIZE_H_
#define _ATOMVISUALIZE_H_

#include "pmdTypes.h"

#include <glut.h>

#ifdef USE_VTK
class vtkPoints; class vtkCellArray; class vtkPolyData; class vtkPolyDataMapper; class vtkUnsignedCharArray; class KeyPressInteractorStyle;
class vtkActor; class vtkRenderWindow; class vtkRenderer; class vtkRenderWindow; class vtkPointSource; class vtkRenderWindowInteractor; class vtkTimerCallback2;
#endif
class Quaternion;


extern SimFlat *blankSimulation(BasePotential *pot);
extern void allocDomains(SimFlat *s);

class AtomVisualize
{
public:
  AtomVisualize();
  ~AtomVisualize() {};
  static void initialize(SimFlat *atoms);
  static void updateData(SimFlat *atoms); 
  static void computeCentrosymmetry(SimFlat *atoms);
  static void colorMap(float min, float max, float value, unsigned char* color);
  static void renderData();
  static void interact();

#ifdef USE_VTK
  static vtkPoints *pts;
  static vtkCellArray *conn;
  static vtkPolyData *poly;
  static vtkUnsignedCharArray *colors;

  static vtkPointSource* pointSource;

  static vtkPolyDataMapper *inputMapper;
  static vtkActor *inputActor;
  static vtkRenderer *ren1;
  static vtkRenderWindow *renWin;
  static vtkRenderWindowInteractor *renWinInteractor;
  static KeyPressInteractorStyle *style;
  static vtkTimerCallback2 *cb;
#else
  static void keyboard(unsigned char key, int x, int y);
  static void idle();
  static void mouse(int button, int state, int x, int y);
  static void motion(int x, int y);

  static int mouse_old_x, mouse_old_y, mouse_buttons;
  static int displayList;
  static float zoomFactor;
#endif

  static SimFlat *sim;
  static bool newDataAvailable;
  static float minCentro;
  static float maxCentro;
  static float minX, maxX, minY, maxY, minZ, maxZ;
  static Quaternion* q;
};

#endif
