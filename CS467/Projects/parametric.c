#include "Tools/FPToolkit.c"
#include "Tools/M3d_matrix_tools.c"

// Global variables
int Tn, Ttypelist[100] ;
double Tvlist[100] ;
double x[20][1000], y[20][1000] ;
int Sn = 0 ;

/*
// Defines Vector2 to store both x and y
struct Vector2 {
  double x ;
  double y ;
};
typedef struct Vector2 Vector2;
*/

// Creates the screen
int initGraphics(int x, int y) {
  G_init_graphics(x, y) ;
  G_rgb(0, 0, 0) ;
  G_clear() ;
  G_rgb(1, 1, 1) ;
  
  return 0; 
}

/*
// Change of notation
Vector2 transformPoint(Vector2 p, double matrix[4][4]) {
  double vector[2] = {p.x, p.y};
  M3d_mat_mult_pt(vector, matrix, vector);
  p.x = vector[0];
  p.y = vector[1];
  return p;
}
*/

// Creates x and y of shapes
void xAy(int shape, double fx, double fy, double lo, double hi, double step) {
  double u ;
  int i = 0 ;
  
  for (u = lo; u <= hi; u += step) {
    x[shape][i] = fx ;
    y[shape][i] = fy ;
    i++;
  }
  
  Sn++;
}

// Plots the shape
int plot(int shapeNum, double v[4][4], double lo, double hi, double step) {
  double u ;
  int i = 0 ;

  for (u = lo ; u <= hi ; u += step) {
    M3d_mat_mult_pt(x[shapeNum], v, x[shapeNum]);
    M3d_mat_mult_pt(y[shapeNum], v, y[shapeNum]);
    G_pixel(x[shapeNum][i], y[shapeNum][i]) ;
    i++ ;
  }

  return 0 ;
}

// Creates the red partial circle
void circle() {
  double v[4][4] ;
  double vi[4][4] ;
  
  Tn = 0 ; 
  Ttypelist[Tn] = SX ; Tvlist[Tn] =   50.0 ; Tn++ ;
  Ttypelist[Tn] = SY ; Tvlist[Tn] =  100.0 ; Tn++ ;
  Ttypelist[Tn] = TX ; Tvlist[Tn] =  300.0 ; Tn++ ;
  Ttypelist[Tn] = TY ; Tvlist[Tn] =  500.0 ; Tn++ ;
  
  M3d_make_movement_sequence_matrix(v, vi, Tn, Ttypelist, Tvlist) ;
  G_rgb(1.0, 0.0, 0.0) ;
  xAy(Sn, cos, sin, 0.25*M_PI, 1.5*M_PI, 10) ;
  plot(Sn, v, 0.25*M_PI, 1.5*M_PI, 10) ;
  G_rgb(1.0, 1.0, 1.0) ;
}

//
void sum4() {
  double v[4][4] ;
  double vi[4][4] ;

  Tn = 0 ; 
  Ttypelist[Tn] = SX ; Tvlist[Tn] =   30.0 ; Tn++ ;
  Ttypelist[Tn] = SY ; Tvlist[Tn] =   30.0 ; Tn++ ;
  Ttypelist[Tn] = TX ; Tvlist[Tn] =  250.0 ; Tn++ ;
  Ttypelist[Tn] = TY ; Tvlist[Tn] =  170.0 ; Tn++ ;
  
  M3d_make_movement_sequence_matrix(v, vi, Tn, Ttypelist, Tvlist) ;
  G_rgb(0.0, 1.0, 0.0) ;
  //plot() ;
  G_rgb(1.0, 1.0, 1.0) ;
}

// Main function, calls the shapes
int main() {
  initGraphics(800, 800) ;
  
  circle() ;

  G_wait_key() ;
}
