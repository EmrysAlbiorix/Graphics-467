#include "Tools/FPToolkit.c"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int S_WIDTH, S_HEIGHT ;

int main() {
  // Initializes the screen
  S_WIDTH = 800 ; S_HEIGHT = 800 ;
  G_init_graphics(S_WIDTH, S_HEIGHT) ;
  G_rgb(0, 0, 0) ;
  G_clear() ;

  // Draws the circle
  G_rgb(0, 1, 0) ;
  G_circle(400, 400, 300) ;

  // Calculate the width and height of rectangle
  int p = 300/(sqrt(2)) ;
  int w = 2*p ;

  // Draws the rectangle
  G_rectangle(400-p, 400-p, w, w) ;

  // Draws another circle
  int n = 300 - p ;
  G_circle(400, 400, 300-n) ;

  // Draws another rectangle
  int m = (300-n)/(sqrt(2)) ;
  G_rectangle(400-m, 400-m, 2*m, 2*m) ;

  // Shows screen until key is pressed
  G_wait_key() ;
}

/** Notes:

    (R  (S  (T  DATA)))
- Instead you should do the Rotation, Scaling, and Translation first to save time complexity

    ((R  (S  T))  DATA)
- This method keeps all the movement as a 4x4 matrix before altering the DATA

- T(300, 40, 10) is a basic translation T
- Ti(-300, -40, -10) is the inverse matrix of T

- S(4, 2, 1/2) is a basic scale S
- Si(1/4, 1/2, 2) is the inverse of matrix S

- V is the matrix of all movement
- Vi is the inverse of matrix V
- V*Vi = I where I is the identity matrix

 **/
