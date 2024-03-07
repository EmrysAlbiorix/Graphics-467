#include "Tools/FPToolkit.c"
#include "Tools/M3d_matrix_tools.c"

double obmat[100][4][4] ;
double obinv[100][4][4] ;
double color[100][3] ;
int num_objects ;

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

void Draw_ellipsoid (int onum) {
  int n,i ;
  double t, xyz[3] ;
  double x,y ;

  G_rgb (color[onum][0],color[onum][1],color[onum][2]) ;
  
  n = 1000 ;
  for (i = 0 ; i < n ; i++) {
    t = i*2*M_PI/n ;
    xyz[0] = cos(t) ;
    xyz[1] = sin(t) ;
    xyz[2] = 0 ;
    M3d_mat_mult_pt(xyz, obmat[onum], xyz) ;
    x = xyz[0] ;
    y = xyz[1] ;
    G_point(x,y) ;
  }

}

// Draws the scene
void Draw_the_scene() {
  int onum ;
  for (onum = 0 ; onum < num_objects ; onum++) {
    Draw_ellipsoid(onum) ;
  }
}


/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

// Solves quadratic equation of the form ax^2 + bx + c = 0
double quadratic_solve(double a, double b, double c) {
  double x1, x2;

  // Calculate the solutions using the quadratic formula
  x1 = (-b + sqrt((b * b) - 4 * a * c)) / (2 * a);
  x2 = (-b - sqrt((b * b) - 4 * a * c)) / (2 * a);

  // Return the lower of the positive values of x1 and x2
  if (x1 > 0 && x1 <= x2) {
    return x1;
  } else if (x2 > 0 && x2 <= x1) {
    return x2;
  }
  // If both values are not positive or if none is lower, no value is returned
}

int test01() {
  double vm[4][4], vi[4][4] ;
  double Tvlist[100] ;
  int Tn, Ttypelist[100] ;
  double m[4][4], mi[4][4] ;
  double Rsource[3] ;
  double Rtip[3] ;
  double argb[3] ;

  //////////////////////////////////////////////////////////////////////
  M3d_make_identity(vm) ;    M3d_make_identity(vi) ; // OVERRIDE for 2d
  //////////////////////////////////////////////////////////////////////

  num_objects = 0 ;

  //////////////////////////////////////////////////////////////
  // Object 0: Green
  //////////////////////////////////////////////////////////////
  color[num_objects][0] = 0.0 ;
  color[num_objects][1] = 0.8 ; 
  color[num_objects][2] = 0.0 ;
	
  Tn = 0 ;
  Ttypelist[Tn] = SX ; Tvlist[Tn] =   60   ; Tn++ ;
  Ttypelist[Tn] = SY ; Tvlist[Tn] =  100   ; Tn++ ;
  Ttypelist[Tn] = RZ ; Tvlist[Tn] =   25   ; Tn++ ;
  Ttypelist[Tn] = TX ; Tvlist[Tn] =  300   ; Tn++ ;
  Ttypelist[Tn] = TY ; Tvlist[Tn] =  200   ; Tn++ ;
	
  M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
  M3d_mat_mult(obmat[num_objects], vm, m) ;
  M3d_mat_mult(obinv[num_objects], mi, vi) ;

  num_objects++ ; // don't forget to do this

  //////////////////////////////////////////////////////////////
  // Object 1: Orange
  //////////////////////////////////////////////////////////////
  color[num_objects][0] = 1.0 ;
  color[num_objects][1] = 0.3 ; 
  color[num_objects][2] = 0.0 ;
	
  Tn = 0 ;
  Ttypelist[Tn] = SX ; Tvlist[Tn] =  180   ; Tn++ ;
  Ttypelist[Tn] = SY ; Tvlist[Tn] =   40   ; Tn++ ;
  Ttypelist[Tn] = RZ ; Tvlist[Tn] =   60   ; Tn++ ;
  Ttypelist[Tn] = TX ; Tvlist[Tn] =  400   ; Tn++ ;
  Ttypelist[Tn] = TY ; Tvlist[Tn] =  550   ; Tn++ ;
	
  M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
  M3d_mat_mult(obmat[num_objects], vm, m) ;
  M3d_mat_mult(obinv[num_objects], mi, vi) ;

  num_objects++ ; // don't forget to do this
  
  //////////////////////////////////////////////////////////////
  // Object 2: Blue
  //////////////////////////////////////////////////////////////
  color[num_objects][0] = 0.3 ;
  color[num_objects][1] = 0.3 ; 
  color[num_objects][2] = 1.0 ;
	
  Tn = 0 ;
  Ttypelist[Tn] = SX ; Tvlist[Tn] =   75   ; Tn++ ;
  Ttypelist[Tn] = SY ; Tvlist[Tn] =   35   ; Tn++ ;
  Ttypelist[Tn] = RZ ; Tvlist[Tn] =  150   ; Tn++ ;
  Ttypelist[Tn] = TX ; Tvlist[Tn] =  360   ; Tn++ ;
  Ttypelist[Tn] = TY ; Tvlist[Tn] =  500   ; Tn++ ;
	
  M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
  M3d_mat_mult(obmat[num_objects], vm, m) ;
  M3d_mat_mult(obinv[num_objects], mi, vi) ;

  num_objects++ ; // don't forget to do this    
  
  //////////////////////////////////////////////////////////////
  // Object 3: Teal    
  //////////////////////////////////////////////////////////////
  color[num_objects][0] = 0.5 ;
  color[num_objects][1] = 1.0 ; 
  color[num_objects][2] = 1.0 ;
	
  Tn = 0 ;
  Ttypelist[Tn] = SX ; Tvlist[Tn] =  130   ; Tn++ ;
  Ttypelist[Tn] = SY ; Tvlist[Tn] =   30   ; Tn++ ;
  Ttypelist[Tn] = RZ ; Tvlist[Tn] =  -15   ; Tn++ ;
  Ttypelist[Tn] = TX ; Tvlist[Tn] =  100   ; Tn++ ;
  Ttypelist[Tn] = TY ; Tvlist[Tn] =  700   ; Tn++ ;
	
  M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
  M3d_mat_mult(obmat[num_objects], vm, m) ;
  M3d_mat_mult(obinv[num_objects], mi, vi) ;

  num_objects++ ; // don't forget to do this        
  //////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////
    
  G_rgb(0,0,0) ;
  G_clear() ;

  Draw_the_scene() ;
  
  // Draws the filter line
  Rsource[0] =  20 ;  Rsource[1] =  400 ;  Rsource[2] = 0 ;    
  G_rgb(1,0,1) ; G_fill_circle(Rsource[0], Rsource[1], 3) ;
  G_rgb(1,0,1) ; G_line(100,200,  100,600) ;
    
  G_wait_key() ;
    
  double ytip ;
  
  for (ytip = 200 ; ytip <= 600 ; ytip++) {
    // Set the tip of the ray
    Rtip[0] = 100; Rtip[1] = ytip; Rtip[2] = 0;   

    // Draw the ray from the source to the tip
    G_rgb(1, 1, 0); G_line(Rsource[0], Rsource[1], Rtip[0], Rtip[1]);
    //		ray (Rsource, Rtip, argb) ; 

    // Draw the scene
    Draw_the_scene();
    
    // Transform the ray endpoints into object-local coordinates for each object
    double RsourceT[3], RtipT[3];
    for (int objnum = 0; objnum < num_objects; objnum++) {
      // Transform the ray source and tip into object-local coordinates
      M3d_mat_mult_pt(RsourceT, obinv[objnum], Rsource);
      M3d_mat_mult_pt(RtipT, obinv[objnum], Rtip);

      // Calculate coefficients of the quadratic equation for the object
      double A, B, C;
      A = (RtipT[0] - RsourceT[0]) * (RtipT[0] - RsourceT[0]) +
          (RtipT[1] - RsourceT[1]) * (RtipT[1] - RsourceT[1]);
      B = 2 * RsourceT[0] * (RtipT[0] - RsourceT[0]) +
          2 * RsourceT[1] * (RtipT[1] - RsourceT[1]);
      C = (RsourceT[1] * RsourceT[1]) + (RsourceT[0] * RsourceT[0]) - 1;

      // Solve the quadratic equation for finding the intersections of the ray
      double t = quadratic_solve(A, B, C);

      // Calculate the intersection point in object-local coordinates
      double intersect[3];
      intersect[0] = RsourceT[0] + t * (RtipT[0] - RsourceT[0]);
      intersect[1] = RsourceT[1] + t * (RtipT[1] - RsourceT[1]);

      // Transform the intersection point back to world coordinates
      M3d_mat_mult_pt(intersect, obmat[objnum], intersect);

      // Draw the intersection point and the ray if it intersects with the object
      if (t > 0) {
          G_rgb(color[objnum][0], color[objnum][1], color[objnum][2]);
          G_line(Rsource[0], Rsource[1], intersect[0], intersect[1]);
          G_fill_circle(Rtip[0], Rtip[1], 1); // Mark the ray tip
      }
    }

    // Wait for user input before proceeding to the next ray
    G_wait_key();
  }

  G_rgb(1,1,1) ; G_draw_string("'q' to quit", 50,50) ;
  while (G_wait_key() != 'q') ;
  G_save_image_to_file("2d_Simple_Raytracer.xwd") ;
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

// Runs the program
int main() {
  G_init_graphics(800,800);
  test01() ;
}
