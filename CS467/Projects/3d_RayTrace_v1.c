#include "Tools/FPToolkit.c"
#include "Tools/M3d_matrix_tools.c"

double obmat[100][4][4];
double obinv[100][4][4];
double color[100][3];
double refSource[3];
double refTip[3];
int num_objects;
int obtype[100]; // 0 sphere, 1 plane

/////////////////////////////////////////////////////////////////////////
// ... (functions from original code) ...
/////////////////////////////////////////////////////////////////////////
// Draws the object
void Draw_object(int onum) {
  int n, i;
  double t, xyz[3];
  double x, y, z;

  n = 1000;
  if (obtype[onum] == 0) { // sphere
    double center[3] = {0, 0, 0};
    M3d_mat_mult_pt(center, obmat[onum], center);
    double radius = 100; // Assuming radius of 100 for the sphere

    for (i = 0; i < n; i++) {
      double phi = 2 * M_PI * i / n;
      double theta = M_PI * i / n;

      xyz[0] = center[0] + radius * sin(theta) * cos(phi);
      xyz[1] = center[1] + radius * sin(theta) * sin(phi);
      xyz[2] = center[2] + radius * cos(theta);

      printf("Point: (%.2f, %.2f, %.2f)\n", xyz[0], xyz[1], xyz[2]);
    }
  } else if (obtype[onum] == 1) { // plane
    // Existing code for drawing the plane
    for (i = 0; i < n; i++) {
      t = -1 + 2.0 * i / (n - 1);
      xyz[0] = t;
      xyz[1] = t;
      xyz[2] = 0;
      M3d_mat_mult_pt(xyz, obmat[onum], xyz);
      x = xyz[0];
      y = xyz[1];
      z = xyz[2];
      printf("Point: (%.2f, %.2f, %.2f)\n", x, y, z);
    }
  }
}


// Draws the Scene
void Draw_the_scene() {
  int onum;
  for (onum = 0; onum < num_objects; onum++) {
    Draw_object(onum);
  }
}


/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

// Solves quadratic equation, returns number of solutions and saves solution in res array
int quadratic_solve(double a, double b, double c, double res[2]) {
  double d = b*b - 4*a*c ;

  if (d < 0) return 0 ;

  if (d == 0) { res[0] = -b/(2*a) ; return 1 ; }


  // Calculate the solutions using the quadratic formula
  res[0] = (-b + sqrt(d)) / (2 * a);
  res[1] = (-b - sqrt(d)) / (2 * a);

  return 2;

}

// Solves quadratic equation of the form ax^2 + bx + c = 0 (OLD)
double OLDquadratic_solve(double a, double b, double c) {
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
    normalizeVector(normal, normal, 1) ;
    normalizeVector(incoming, incoming, 1) ;
    
    double dotProd = dotProduct(normal, incoming, 3) ;
    
    reflection[0] = incoming[0] - 2 * dotProd * normal[0];
    reflection[1] = incoming[1] - 2 * dotProd * normal[1];
    reflection[2] = incoming[2] - 2 * dotProd * normal[2];
}

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
// Does the actual RayTracing
//////////////////////////////////////////////////////////////
int rayThing(double Rsource[], double Rtip[]) {
  double normal[3];
  int RGBnum;
  double RsourceT[3], RtipT[3];
  double tbest = 1e50;
  int obest = -1;
  int n, nn;
  double qres[2];
  double tq;

  for (int objnum = 0; objnum < num_objects; objnum++) {
    M3d_mat_mult_pt(RsourceT, obinv[objnum], Rsource);
    M3d_mat_mult_pt(RtipT, obinv[objnum], Rtip);

    double A, B, C;
    if (obtype[objnum] == 0) { // sphere
      A = (RtipT[0] - RsourceT[0]) * (RtipT[0] - RsourceT[0]) +
          (RtipT[1] - RsourceT[1]) * (RtipT[1] - RsourceT[1]) +
          (RtipT[2] - RsourceT[2]) * (RtipT[2] - RsourceT[2]);
      B = 2 * (RsourceT[0] * (RtipT[0] - RsourceT[0]) +
               RsourceT[1] * (RtipT[1] - RsourceT[1]) +
               RsourceT[2] * (RtipT[2] - RsourceT[2]));
      C = (RsourceT[0] * RsourceT[0]) + (RsourceT[1] * RsourceT[1]) +
          (RsourceT[2] * RsourceT[2]) - 1;
    } else { // plane
      A = (RtipT[0] - RsourceT[0]) * (RtipT[0] - RsourceT[0]) +
          (RtipT[1] - RsourceT[1]) * (RtipT[1] - RsourceT[1]) +
          (RtipT[2] - RsourceT[2]) * (RtipT[2] - RsourceT[2]);
      B = 2 * RsourceT[2] * (RtipT[2] - RsourceT[2]);
      C = RsourceT[2] * RsourceT[2];
    }

    double t = 1e50;
    if (obtype[objnum] == 0) { // sphere
      nn = quadratic_solve(A, B, C, qres);
      tq = 1e50;

      for (int ii = 0; ii < nn; ii++) {
        if (qres[ii] > 0 && qres[ii] < tq) {
          tq = qres[ii];
        }
      }

      if (tq < tbest) {
        tbest = tq;
        obest = objnum;
      }
    } else if (obtype[objnum] == 1) { // plane
      if (B == 0) {
        continue;
      }
      t = -C / B;
      if (t < 0) {
        continue;
      }

      if (t < tbest) {
        tbest = t;
        obest = objnum;
      }
    }
  }

  if (obest == -1) {
    printf("Error: obest == -1\n");
  } else {
    printf("obest = %d tbest = %lf\n", obest, tbest);
    RGBnum = obest;

    double intersect[3];
    intersect[0] = Rsource[0] + tbest * (Rtip[0] - Rsource[0]);
    intersect[1] = Rsource[1] + tbest * (Rtip[1] - Rsource[1]);
    intersect[2] = Rsource[2] + tbest * (Rtip[2] - Rsource[2]);

    M3d_mat_mult_pt(intersect, obmat[obest], intersect);

    // Check if the intersection point lies on the sphere
    if (obtype[obest] == 0) {
      double center[3] = {0, 0, 0};
      M3d_mat_mult_pt(center, obmat[obest], center);
      double radius = 100; // Assuming radius of 100 for the sphere
      double distFromCenter = sqrt((intersect[0] - center[0]) * (intersect[0] - center[0]) +
                                   (intersect[1] - center[1]) * (intersect[1] - center[1]) +
                                   (intersect[2] - center[2]) * (intersect[2] - center[2]));
      if (fabs(distFromCenter - radius) < 1e-6) {
        printf("Intersection point: (%.2f, %.2f, %.2f)\n", intersect[0], intersect[1], intersect[2]);
      }
    }

    getNormal(obest, normal, intersect);

    double len = sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);
    normal[0] /= len;
    normal[1] /= len;
    normal[2] /= len;

    printf("Normal vector: (%.2f, %.2f, %.2f)\n", normal[0], normal[1], normal[2]);

    double incoming[3];
    double reflec[3];

    incoming[0] = Rtip[0] - Rsource[0];
    incoming[1] = Rtip[1] - Rsource[1];
    incoming[2] = Rtip[2] - Rsource[2];

    getReflect(normal, incoming, reflec);

    refSource[0] = intersect[0] + 0.01 * reflec[0];
    refTip[0] = intersect[0] + 50 * reflec[0];
    refSource[1] = intersect[1] + 0.01 * reflec[1];
    refTip[1] = intersect[1] + 50 * reflec[1];
    refSource[2] = intersect[2] + 0.01 * reflec[2];
    refTip[2] = intersect[2] + 50 * reflec[2];

    printf("Reflection vector: (%.2f, %.2f, %.2f)\n", reflec[0], reflec[1], reflec[2]);

    printf("Color: (%.2f, %.2f, %.2f)\n", color[RGBnum][0], color[RGBnum][1], color[RGBnum][2]);
    printf("Ray: (%.2f, %.2f, %.2f) -> (%.2f, %.2f, %.2f)\n", Rsource[0], Rsource[1], Rsource[2], Rtip[0], Rtip[1], Rtip[2]);
  }
}

// ... (functions from original code) ...

// Sets the sphere color to blue and the plane color to green
void setColors() {
    color[0][0] = 0.0;
    color[0][1] = 0.0;
    color[0][2] = 1.0; // Blue sphere

    color[1][0] = 0.0;
    color[1][1] = 0.5;
    color[1][2] = 0.0; // Green plane
}

// Draws the scene with the sphere in the center and the reflective floor
/*void Draw_the_scene() {
  int onum;
  for (onum = 0; onum < num_objects; onum++) {
    Draw_object(onum);
  }
}*/

int main() {
    // Set up the scene with a sphere and a plane
    num_objects = 2;
    obtype[0] = 0; // Sphere
    obtype[1] = 1; // Plane

    // Set the colors
    setColors();

    // Set up the transformation matrices for the sphere and the plane
    // Here, we assume that the sphere is centered at (0, 0, 0) with a radius of 100
    // and the plane is at z = 0 with a size of 1000 x 1000
    M3d_make_translation(obmat[0], 0.0, 0.0, 0.0);
    M3d_make_scaling(obmat[1], 500.0, 500.0, 1.0);
    M3d_make_translation(obmat[1], 0.0, 0.0, -250.0);

    // Calculate the inverse matrices
    //M3d_mat_inverse(obinv[0], obmat[0]);
    //M3d_mat_inverse(obinv[1], obmat[1]);

    // Draw the scene
    Draw_the_scene();

    return 0;
}
