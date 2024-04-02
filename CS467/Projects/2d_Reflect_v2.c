#include "Tools/FPToolkit.c"
#include "Tools/M3d_matrix_tools.c"


double obmat[100][4][4] ;
double obinv[100][4][4] ;
double color[100][3] ;
int    num_objects ;
int    obtype[100] ; // 0 line segment, 1 circle, 2 hyperbola


/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////


// Draws the object
void Draw_object (int onum)
{
  int n,i ;
  double t, xyz[3] ;
  double x,y ;

  G_rgb (color[onum][0],color[onum][1],color[onum][2]) ;
  
  n = 1000 ;
  for (i = 0 ; i < n ; i++) {

    if (obtype[onum] == 0) {
      // horizontal line segment from (-1,0) to (1,0)
      t = -1 + 2.0*i/(n-1) ;
      xyz[0] = t ;
      xyz[1] = 0 ;
      xyz[2] = 0 ;
      M3d_mat_mult_pt(xyz, obmat[onum], xyz) ;
      x = xyz[0] ;
      y = xyz[1] ;
      G_point(x,y) ;
    } else if (obtype[onum] == 1) {
      // unit circle
      t = i*2*M_PI/n ;
      xyz[0] = cos(t) ;
      xyz[1] = sin(t) ;
      xyz[2] = 0 ;
      M3d_mat_mult_pt(xyz, obmat[onum], xyz) ;
      x = xyz[0] ;
      y = xyz[1] ;
      G_point(x,y) ;      
    } else if (obtype[onum] == 2) {
      // Trunkated Hyperbola
      t = -1 + 2.0*i/(n-1) ;
      
      xyz[0] = sqrt(1 + t*t) ;
      xyz[1] = t ;
      xyz[2] = 0 ;
      M3d_mat_mult_pt(xyz, obmat[onum], xyz) ;
      x = xyz[0] ;
      y = xyz[1] ;
      G_point(x,y) ;

      xyz[0] = -sqrt(1 + t*t) ;
      xyz[1] = t ;
      xyz[2] = 0 ;
      M3d_mat_mult_pt(xyz, obmat[onum], xyz) ;
      x = xyz[0] ;
      y = xyz[1] ;
      G_point(x,y) ;      
      
    }
    

  }

}


// Draws the Scene
void Draw_the_scene()
{
  int onum ;
  for (onum = 0 ; onum < num_objects ; onum++) {
    Draw_object(onum) ;
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

// Solves hyperbola's other side (returns the furthur x value)
double hyper_solve(double a, double b, double c) {
  double x1, x2;

  // Calculate the solutions using the quadratic formula
  x1 = (-b + sqrt((b * b) - 4 * a * c)) / (2 * a);
  x2 = (-b - sqrt((b * b) - 4 * a * c)) / (2 * a);

  // Return the lower of the positive values of x1 and x2
  if (x1 > 0 && x1 <= x2) {
    return x2;
  } else if (x2 > 0 && x2 <= x1) {
    return x1;
  }
  // If both values are not positive or if none is lower, no value is returned
}


/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

// Sets a standard length of scale
void normalizeVector(double v[], double res[], double scale){
  double magnitude = sqrt(v[0] * v[0] + v[1] * v[1]);
  
  res[0] = v[0]/magnitude * scale ;
  res[1] = v[1]/magnitude * scale ;
   
}

// Finds the cross product of v1 and v2
void crossProduct(double v1[], double v2[], double result[]) {
  result[0] = v1[1] * v2[2] - v1[2] * v2[1];
  result[1] = v1[2] * v2[0] - v1[0] * v2[2];
  result[2] = v1[0] * v2[1] - v1[1] * v2[0];
}

// Finds the dot product of v1 and v2 based on size
double dotProduct(double v1[], double v2[], int size) {
    double result = 0;
    for (int i = 0; i < size; i++) {
        result += v1[i] * v2[i];
    }
    return result;
}

//////////////////////////////////////////////////////////////
// Finds the normal vector
//////////////////////////////////////////////////////////////
void getNormal(int objnum, double normal[], double intersect[]) {
  // Moves back to obj-space
  double intersectObjSpace[3];
  M3d_mat_mult_pt(intersectObjSpace, obinv[objnum], intersect) ;
  
  // Calculate the normal vector in object-local coordinates
  //printf("%d, \t %lf, %lf\n", objnum, intersect[0], intersect[1]) ;
  switch(obtype[objnum]){
    case 0: // line
      normal[0] = intersectObjSpace[0] ;
      normal[1] = 1 ;
      normal[2] = 0 ;
    break;
    case 1: // circle
      normal[0] = 2 * intersectObjSpace[0] ;
      normal[1] = 2 * intersectObjSpace[1] ;
      normal[2] = 0 ;
    break;
    case 2: // hyperbola
      normal[0] = intersectObjSpace[0] ;
      normal[1] = 1 ;
      normal[2] = 0 ;
    break;
  }
  
  // Transoposes the shit
  double transpose[4][4] = {
  	obinv[objnum][0][0], obinv[objnum][1][0], obinv[objnum][2][0], 0,
  	obinv[objnum][0][1], obinv[objnum][1][1], obinv[objnum][2][1], 0,
  	obinv[objnum][0][2], obinv[objnum][1][2], obinv[objnum][2][2], 0,
  	0, 0, 0, 0} ;

  // Math
  M3d_mat_mult_pt(normal, transpose, normal) ;
  //printf("Normal: %lf, %lf\n", normal[0], normal[1]) ;
}

//////////////////////////////////////////////////////////////
// Finds Reflection Vector (Blue)
//////////////////////////////////////////////////////////////
void getReflect(double normal[], double incoming[], double reflection[]) {
    // Calculate reflection vector using incoming ray direction and surface normal
    // R = I - 2 * dot(N, I) * N
    double dotProd = dotProduct(normal, incoming, 3) ;
    
    normalizeVector(normal, normal, 25) ;
    
    reflection[0] = incoming[0] - 2 * dotProd * normal[0];
    reflection[1] = incoming[1] - 2 * dotProd * normal[1];
    reflection[2] = incoming[2] - 2 * dotProd * normal[2];
}

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
// Does the actual RayTracing shit
//////////////////////////////////////////////////////////////
int rayThing(double Rsource[], double Rtip[]) {
  // Stores the temporary distances and color
  double xBuff[3] ;
  double normal[3] ;
  xBuff[0] = -1 ;
  int RGBnum ;

  // Transform the ray endpoints into object-local coordinates for each object
  double RsourceT[3], RtipT[3];
  for (int objnum = 0; objnum < num_objects; objnum++) {
    // Transform the ray source and tip into object-local coordinates
    M3d_mat_mult_pt(RsourceT, obinv[objnum], Rsource);
    M3d_mat_mult_pt(RtipT, obinv[objnum], Rtip);

    // Calculate coefficients of the quadratic equation for the object
    double A, B, C;
    if (obtype[objnum] == 0 || obtype[objnum] == 1) {
      // Finds the intersections for circles and lines?
      A = (RtipT[0] - RsourceT[0]) * (RtipT[0] - RsourceT[0]) +
          (RtipT[1] - RsourceT[1]) * (RtipT[1] - RsourceT[1]);
      B = 2 * RsourceT[0] * (RtipT[0] - RsourceT[0]) +
          2 * RsourceT[1] * (RtipT[1] - RsourceT[1]);
      C = (RsourceT[1] * RsourceT[1]) + (RsourceT[0] * RsourceT[0]) - 1;
    } else {
      // Finds the intersections for the hyperbola
      A = (RtipT[0] - RsourceT[0]) * (RtipT[0] - RsourceT[0]) -
          (RtipT[1] - RsourceT[1]) * (RtipT[1] - RsourceT[1]);
      B = 2 * RsourceT[0] * (RtipT[0] - RsourceT[0]) -
          2 * RsourceT[1] * (RtipT[1] - RsourceT[1]);
      C = (RsourceT[0] * RsourceT[0]) - (RsourceT[1] * RsourceT[1]) - 1;
    }

    // Solve the quadratic equation for finding the intersections of the ray
    double t = 0, t1 = 0, t2 = 0 ;
    if (obtype[objnum] == 0) { // line
      double line[3] = {RtipT[0] - RsourceT[0], RtipT[1] - RsourceT[1], RtipT[2] - RsourceT[2]} ;
      t = -RsourceT[1]/line[1] ;
      double xt = RsourceT[0] + (line[0]*t) ;
      //printf("t: %lf,	xt: %lf\n", t, xt) ;
      if (t < 0 || fabs(xt) > 1) {
        continue ;
      }
    } else if (obtype[objnum] == 1) { // circle
      t = quadratic_solve(A, B, C) ;
    } else if (obtype[objnum] == 2) { // hyperbola
      double trunkY1 = 0, trunkY2 = 0 ;
      //hyper_solve(A,B,C,t1,t2);
      //t1 = (-B + sqrt((B * B) - 4 * A * C)) / (2 * A) ;
      t1 = quadratic_solve(A, B, C) ;
      //t2 = (-B - sqrt((B * B) - 4 * A * C)) / (2 * A) ;
      t2 = hyper_solve(A, B, C) ;
      trunkY1 = RsourceT[1] + (RtipT[1] - RsourceT[1]) * t1 ;
      trunkY2 = RsourceT[1] + (RtipT[1] - RsourceT[1]) * t2 ;
      //printf("%lf,%lf\n",trunkY1,trunkY2) ;
      if ((trunkY1 > 1 || trunkY1 < -1) && trunkY2 <= 1 && trunkY2 >= -1) {
        t = t2 ;
      } else if (trunkY1 <= 1 && trunkY1 >= -1 && (trunkY2 > 1 || trunkY2 < -1)) {
        t = t1 ;
      } else if (trunkY1 <= 1 && trunkY1 >= -1 && trunkY2 <= 1 && trunkY2 >= -1) {
        //t = (t1 < t2) ? t1 : t2 ; This is shorthand for the below if-else statement
        if (t1 < t2) {
          t = t1;
        } else {
          t = t2;
        }
      } else {
        continue ;
      }    
    } else {
      printf("%d is not an object type.\n", obtype[objnum]) ;
    }

    // Calculate the intersection point in object-local coordinates
    double intersect[3];
    intersect[0] = RsourceT[0] + t * (RtipT[0] - RsourceT[0]);
    intersect[1] = RsourceT[1] + t * (RtipT[1] - RsourceT[1]);

    // Transform the intersection point back to world coordinates
    M3d_mat_mult_pt(intersect, obmat[objnum], intersect);
    
    // Plots a point at the intersect to test
    G_rgb(1, 1, 1) ;
    G_fill_circle(intersect[0], intersect[1], 2) ;

    if (t > 0) {
      // Finds the distance between intersect and Rsource
      double xMath = (intersect[0] - Rsource[0]) * (intersect[0] - Rsource[0]) ;
      double yMath = (intersect[1] - Rsource[1]) * (intersect[1] - Rsource[1]) ;
      double distance = sqrt(xMath + yMath) ;
      
      // Saves the point that is closer to Rsource if intersects with objects
      if (intersect[0] < xBuff[0] || xBuff[0] < 0) {
        xBuff[0] = intersect[0] ;
        xBuff[1] = intersect[1] ;
        xBuff[2] = intersect[2] ;
        
        // Saves the correct color
        RGBnum = objnum ;
      }
    }
  }
  
  // Checks if < 0
  if (xBuff[0] < 0) {
    return 0 ;
  }
  
  // Find Normal
  getNormal(RGBnum, normal, xBuff) ;
  
  // Define the incoming ray direction
  double incoming[3] = {Rtip[0] - Rsource[0], Rtip[1] - Rsource[1], Rtip[2] - Rsource[2]};
  
  // Calculate reflection vector
  double reflection[3];
  getReflect(normal, incoming, reflection);
  
  // Normalize reflection vector
  normalizeVector(reflection, reflection, 50) ;

  
  // Does the drawing
  G_rgb(color[RGBnum][0], color[RGBnum][1], color[RGBnum][2]) ; // Shape color
  G_fill_circle(Rtip[0], Rtip[1], 1) ; // Mark the ray tip
  G_line(Rsource[0], Rsource[1], Rtip[0], Rtip[1]) ; // Inner Line
  G_rgb(.8, .8, .8) ;
  //G_line(Rtip[0], Rtip[1], xBuff[0], xBuff[1]) ; // Outer Line
  //G_line(xBuff[0], xBuff[1], xBuff[0] + normal[0], xBuff[1] + normal[1]) ; // Draws the Normal Vector
  G_rgb(1, 0, 0) ;
  G_line(xBuff[0] + 0.001*reflection[0], xBuff[1] + 0.001*reflection[1], xBuff[0] + reflection[0], xBuff[1] + reflection[1]); // Draws Reflection Vector
}
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


int test01()
{
  double vm[4][4], vi[4][4];
  double Tvlist[100];
  int Tn, Ttypelist[100];
  double m[4][4], mi[4][4];
  double Rsource[3];
  double Rtip[3];
  double argb[3] ;

    //////////////////////////////////////////////////////////////////////
    M3d_make_identity(vm) ;    M3d_make_identity(vi) ; // OVERRIDE for 2d
    //////////////////////////////////////////////////////////////////////

    num_objects = 0 ;

    //////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////
    obtype[num_objects] = 1 ; // circle
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
    obtype[num_objects] = 1 ; // circle    
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
    obtype[num_objects] = 1 ; // circle    
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
    obtype[num_objects] = 1 ; // circle    
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
    obtype[num_objects] = 0 ; // line segment
    color[num_objects][0] = 0.5 ;
    color[num_objects][1] = 0.5 ; 
    color[num_objects][2] = 1.0 ;
	
    Tn = 0 ;
    Ttypelist[Tn] = SX ; Tvlist[Tn] =   50   ; Tn++ ;
    Ttypelist[Tn] = RZ ; Tvlist[Tn] =  110   ; Tn++ ;
    Ttypelist[Tn] = TX ; Tvlist[Tn] =  300   ; Tn++ ;
    Ttypelist[Tn] = TY ; Tvlist[Tn] =  300   ; Tn++ ;
	
    M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
    M3d_mat_mult(obmat[num_objects], vm, m) ;
    M3d_mat_mult(obinv[num_objects], mi, vi) ;

    num_objects++ ; // don't forget to do this        
    //////////////////////////////////////////////////////////////
    obtype[num_objects] = 2 ; // hyperbola
    color[num_objects][0] = 0.4 ;
    color[num_objects][1] = 0.2 ; 
    color[num_objects][2] = 0.1 ;
	
    Tn = 0 ;
    Ttypelist[Tn] = SX ; Tvlist[Tn] =   15  ; Tn++ ;
    Ttypelist[Tn] = SY ; Tvlist[Tn] =   80   ; Tn++ ;    
    Ttypelist[Tn] = RZ ; Tvlist[Tn] =   -7  ; Tn++ ;
    Ttypelist[Tn] = TX ; Tvlist[Tn] =  200   ; Tn++ ;
    Ttypelist[Tn] = TY ; Tvlist[Tn] =  630   ; Tn++ ;
	
    M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
    M3d_mat_mult(obmat[num_objects], vm, m) ;
    M3d_mat_mult(obinv[num_objects], mi, vi) ;

    num_objects++ ; // don't forget to do this        
    //////////////////////////////////////////////////////////////    

    G_rgb(0,0,0) ;
    G_clear() ;

    Draw_the_scene() ;
    
    // Creates point at first click
    double P[2] ;
    G_wait_click(P) ;
    G_rgb(1,0,1) ;
    G_fill_circle(P[0], P[1], 2) ;
    
    // Creates point at second click
    double Q[2] ;
    G_wait_click(Q) ;
    G_rgb(1,0,1) ;
    G_fill_circle(Q[0], Q[1], 2) ;
    
    // Defines the points to send to rayThing
    Rsource[0] = P[0] ; Rsource[1] = P[1] ; Rsource[2] = 0 ;
    Rtip[0] = Q[0] ; Rtip[1] = Q[1] ; Rtip[2] = 0 ;
    
    // Draw the scene
    Draw_the_scene();
  
    // RayTrace thing, takes in two arrays (Still computes for further left on x axis)
    rayThing(Rsource, Rtip) ;

    // Wait for user input before proceeding to the next ray
    //G_wait_key();

    G_rgb(1,1,1) ; G_draw_string("'q' to quit", 50,50) ;
    while (G_wait_key() != 'q') ;
    G_save_image_to_file("2d_Reflect_v1.xwd") ;
}




//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////




int main()
{
  G_init_graphics(800,800);
  test01() ;
}
