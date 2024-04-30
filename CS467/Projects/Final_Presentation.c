#include "Tools/FPToolkit.c"
#include "Tools/M3d_matrix_tools.c"
#include "../XWDhelp/XWD_TOOLS_03/xwd_tools_03.c"

// Constants for lighting model
#define AMBIENT 0.2
#define MAX_DIFFUSE 0.5
#define SPECPOW 50

#define hither 2
#define yon 100
#define REF_DEPTH 4
#define H 0.57735

// Light source position
double LightX = -8;
double LightY = 10;
double LightZ = -2;

// Arrays to store object information
double ARGB[3];
double obmat[100][4][4];
double obinv[100][4][4];
double color[100][3];
double reflectivity[100];
double refractive_index[100];
int obtype[100];
int num_objects;

// Light vectors in different spaces
double light_in_eye_space[3];
double light_in_world_space[3];

// Temporary variables used for calculations
double T;
int OB;

// Function to compute lighting model for a surface point
// Parameters:
//   irgb: Array containing the surface color (input)
//   s: Array containing the position of the observer (input)
//   p: Array containing the position of the surface point (input)
//   n: Array containing the surface normal at the point (input)
//   argb: Array to store the resulting color (output)
// Returns:
//   1 if successful, 0 otherwise
int Light_Model(double irgb[3], double s[3], double p[3], double n[3], double argb[3]) {
  // Calculate normalized normal vector N
  double len;
  double N[3];
  len = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
  if (len == 0) { return 0; }
  N[0] = n[0] / len;  N[1] = n[1] / len;  N[2] = n[2] / len;

  // Calculate normalized eye vector E
  double E[3];
  E[0] = s[0] - p[0];
  E[1] = s[1] - p[1];
  E[2] = s[2] - p[2];
  len = sqrt(E[0] * E[0] + E[1] * E[1] + E[2] * E[2]);
  
  if (len == 0) { 
  	return 0; 
  }
  
  E[0] /= len;  E[1] /= len;  E[2] /= len;
  double NdotE = N[0] * E[0] + N[1] * E[1] + N[2] * E[2];

  // Calculate normalized light vector L
  double L[3];
  L[0] = light_in_eye_space[0] - p[0];
  L[1] = light_in_eye_space[1] - p[1];
  L[2] = light_in_eye_space[2] - p[2];
  len = sqrt(L[0] * L[0] + L[1] * L[1] + L[2] * L[2]);
  
  if (len == 0) { 
  	return 0; 
  }
  
  L[0] /= len;  L[1] /= len;  L[2] /= len;
  double NdotL = N[0] * L[0] + N[1] * L[1] + N[2] * L[2];

  double max_ambient_and_diffuse = AMBIENT + MAX_DIFFUSE;

  // Initialize intensity
  double intensity;

  // Check if eye and light are on opposite sides of polygon
  if (NdotL * NdotE < 0) {
    intensity = AMBIENT;
  } else if ((NdotL < 0) && (NdotE < 0)) {
    // Check if eye and light are on same side but normal pointing "wrong" way
    N[0] *= (-1.0);    N[1] *= (-1.0);    N[2] *= (-1.0);
    NdotL *= (-1.0);
  }

  // Calculate reflection vector R
  double R[3];
  R[0] = 2 * NdotL * N[0] - L[0];
  R[1] = 2 * NdotL * N[1] - L[1];
  R[2] = 2 * NdotL * N[2] - L[2];
  double EdotR = E[0] * R[0] + E[1] * R[1] + E[2] * R[2];

  // Calculate diffuse component
  double diffuse = (NdotL <= 0.0) ? 0.0 : MAX_DIFFUSE * NdotL;

  // Calculate specular component
  double specular = (EdotR <= 0.0) ? 0.0 : (1.0 - max_ambient_and_diffuse) * pow(EdotR, SPECPOW);

  // Calculate total intensity
  intensity = AMBIENT + diffuse + specular;

  // Apply intensity to irgb color components
  double f, g;
  if (intensity <= max_ambient_and_diffuse) {
    f = intensity / max_ambient_and_diffuse;
    argb[0] = f * irgb[0];
    argb[1] = f * irgb[1];
    argb[2] = f * irgb[2];
      
  } else {
    f = (intensity - max_ambient_and_diffuse) / (1.0 - max_ambient_and_diffuse);
    g = 1.0 - f;
    argb[0] = g * irgb[0] + f;
    argb[1] = g * irgb[1] + f;
    argb[2] = g * irgb[2] + f;
  }

  return 1;
}

// Calculate lighting effects on an object using the Light_Model function
// irgb: inherent color of the object (input)
// p0, p1, p2: points defining the object
// argb: actual color of the object (output)
void light_model(double irgb[3], double p0[3], double p1[3], double p2[3], double argb[3]) {
  // Define the eye position (set to the origin for simplicity)
  double Eye[3] = {0, 0, 0};

  // Assign the first point of the object to P
  double P[3] = {p0[0], p0[1], p0[2]};

  // Calculate vectors a and b using the difference between points
  double a[3] = {p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]};
  double b[3] = {p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]};

  // Calculate the normal vector N of the object's surface
  double N[3];
  M3d_x_product(N, a, b); // Assuming M3d_x_product is a function to compute the cross product

  // Apply lighting effects using the Light_Model function
  Light_Model(irgb, Eye, P, N, argb);
}


// Calculate the magnitude of a 3D vector
// v: input vector
double magnitude(double v[3]) {
  // Calculate the magnitude using the Euclidean norm
  double mag = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);

  // Return the magnitude
  return mag;
}


// Find the reflection vector given the intersection point, light source, and surface normal
// intersection: coordinates of the intersection point
// source: coordinates of the light source
// N: surface normal vector
// R: reflection vector (output)
int find_reflection(double intersection[3], double source[3], double N[3], double R[3]) {
  // Calculate the direction vector of the incoming light ray (L)
  double L[3];
  L[0] = source[0] - intersection[0];
  L[1] = source[1] - intersection[1];
  L[2] = source[2] - intersection[2];

  // Normalize the direction vector of the incoming light ray (L)
  double mag = magnitude(L);
  L[0] /= mag;
  L[1] /= mag;
  L[2] /= mag;

  // Calculate the dot product of the surface normal (N) and the incoming light ray (L)
  double NdotL = L[0] * N[0] + L[1] * N[1] + L[2] * N[2];

  // Calculate the reflection vector (R)
  R[0] = 2 * NdotL * N[0] - L[0];
  R[1] = 2 * NdotL * N[1] - L[1];
  R[2] = 2 * NdotL * N[2] - L[2];

  // Return success
  return 1;
}


// Calculate the surface normal at the intersection point using the object matrix
// intersect: coordinates of the intersection point
// normal: surface normal vector (output)
// M: object transformation matrix
// obnum: object number
int find_normal(double intersect[3], double normal[3], double M[4][4], int obnum) {
  // Check the type of the object
  if (obtype[obnum] == 1) { // Sphere
    // Calculate the surface normal for a sphere
    normal[0] = 2 * intersect[0];
    normal[1] = 2 * intersect[1];
    normal[2] = 2 * intersect[2];

    // Apply the inverse object transformation matrix to the normal vector
    double nx = M[0][0] * normal[0] + M[1][0] * normal[1] + M[2][0] * normal[2];
    double ny = M[0][1] * normal[0] + M[1][1] * normal[1] + M[2][1] * normal[2];
    double nz = M[0][2] * normal[0] + M[1][2] * normal[1] + M[2][2] * normal[2];

    // Update the normal vector
    normal[0] = nx;
    normal[1] = ny;
    normal[2] = nz;

    // Normalize the normal vector
    double mag = magnitude(normal);
    normal[0] /= mag;
    normal[1] /= mag;
    normal[2] /= mag;

    return 1;
  } else if (obtype[obnum] == 0) { // Cylinder
    // Calculate the surface normal for a cylinder
    normal[0] = 2 * intersect[0];
    normal[1] = 0;
    normal[2] = 2 * intersect[2];

    // Apply the inverse object transformation matrix to the normal vector
    double nx = M[0][0] * normal[0] + M[1][0] * normal[1] + M[2][0] * normal[2];
    double ny = M[0][1] * normal[0] + M[1][1] * normal[1] + M[2][1] * normal[2];
    double nz = M[0][2] * normal[0] + M[1][2] * normal[1] + M[2][2] * normal[2];

    // Update the normal vector
    normal[0] = nx;
    normal[1] = ny;
    normal[2] = nz;

    // Normalize the normal vector
    double mag = magnitude(normal);
    normal[0] /= mag;
    normal[1] /= mag;
    normal[2] /= mag;

    return 1;
  } else if (obtype[obnum] == 2) { // Plane
    // Set the normal vector for a plane
    normal[0] = 0;
    normal[1] = 0;
    normal[2] = 1;

    // Apply the object transformation matrix to the normal vector
    double nx = M[0][2];
    double ny = M[1][2];
    double nz = M[2][2];

    // Update the normal vector
    normal[0] = nx;
    normal[1] = ny;
    normal[2] = nz;

    // Normalize the normal vector
    double mag = magnitude(normal);
    normal[0] /= mag;
    normal[1] /= mag;
    normal[2] /= mag;

    return 1;
  }
}


// Find intersection point between a ray and an object
// source: starting point of the ray
// tip: endpoint of the ray
// intersect: intersection point (output)
// obnum: object number
// Returns 1 if intersection is found, 0 otherwise
int find_intersection(double source[3], double tip[3], double intersect[3], int obnum) {
	// returns 1 if intersection is found; returns 0 if none

	double x0, y0, z0;
	double a, b, c;
	double t1, t2, t;

	double point_y1, point_y2, point_x;


	x0 = tip[0] - source[0];
	y0 = tip[1] - source[1];
	z0 = tip[2] - source[2];

	////////////////////////////////////////////////////////////////////////////////////

	if(obtype[obnum] == 1) {
		// quadratic equation for t
		a = x0*x0 + y0*y0 + z0*z0 ;
		b = 2*x0*source[0] + 2*y0*source[1] + 2*z0*source[2] ;
		c = source[0]*source[0] + source[1]*source[1] + source[2]*source[2] - 1 ;


		// if inside sqrt is negative return 0 (no intersections)
		if(b*b - 4*a*c > 0) {

			t1 = (-b + sqrt(b*b - 4*a*c))/(2*a) ;
			t2 = (-b - sqrt(b*b - 4*a*c))/(2*a) ;

			if(t1 <= t2) {
				t = t1 ;
			} else {
				t = t2 ;
			}

			if(t < T && t > 0) {

				T = t ;
				OB = obnum ;

				intersect[0] = source[0] + x0*T ;
				intersect[1] = source[1] + y0*T ;
				intersect[2] = source[2] + z0*T ;
			}

		} // end if/else

	} else if(obtype[obnum] == 0) {

		double point_y1, point_y2 ;
		// quadratic equation for t
		a = x0*x0 + z0*z0 ;
		b = 2*x0*source[0] + 2*z0*source[2] ;
		c = source[0]*source[0]+ source[2]*source[2] - 1 ;


		// if inside sqrt is negative return 0 (no intersections)
		if(b*b - 4*a*c > 0) {

			t1 = (-b + sqrt(b*b - 4*a*c))/(2*a) ;
			t2 = (-b - sqrt(b*b - 4*a*c))/(2*a) ;

			point_y1 = source[1] + t1*y0 ;
			point_y2 = source[1] + t2*y0 ;

			if(t1 <= t2) {
				t = t1 ;
			} else {
				t = t2 ;
			}

			if(t < T && t > 0) {

				T = t ;
				OB = obnum ;

				intersect[0] = source[0] + x0*T ;
				intersect[1] = source[1] + y0*T ;
				intersect[2] = source[2] + z0*T ;
			} 
		} // end if b

	} else if (obtype[obnum] == 2) { // Plane
    double t = -(source[2] / z0);

    if (t > 0 && t < T) {
      T = t;
      OB = obnum;

      intersect[0] = source[0] + x0 * T;
      intersect[1] = source[1] + y0 * T;
      intersect[2] = source[2] + z0 * T;
    }
	}

	if (T== 100000000) {
		return 0 ;
	} else {
		return 1 ;
	}
}


// Finds dot product
double dot_product(double N[3], double source[3]) {
	double result ;
	for (int i = 0; i < 3; i++) {
		result += N[i] * source[i] ;
	}
	
	return result ;
}

// Performs refraction of a ray using Snell's law
// I: Incident ray direction vector
// N: Surface normal vector
// eta: Refractive index ratio (n1 / n2), where n1 is the original medium and n2 is the new medium
// R: Refracted ray direction vector (output)
void refract(double I[3], double N[3], double eta, double R[3]) {
  // Calculate the refractive index ratio
  eta = 2.0 - eta;

  // Calculate the cosine of the angle of incidence
  double cosi = dot_product(N, I);

  // Calculate the direction of the refracted ray using Snell's law
  R[0] = I[0] * eta - N[0] * (eta * cosi - cosi);
  R[1] = I[1] * eta - N[1] * (eta * cosi - cosi);
  R[2] = I[2] * eta - N[2] * (eta * cosi - cosi);
}

//------------------------------------------------------------------------------------------------

// Test function to set up objects in scene
int test01() {
	int Tn, Ttypelist[100] ;
	double Tvlist[100] ;
	double Rsource[3] ;
	double Rtip[3] ;
	double transparency[100];
	double vm[4][4], vi[4][4] ;
	double m[4][4], mi[4][4] ;
	double irgb[100][3] ;
	double argb[100][3] ;
	double lrgb[100][3] ;

  // Set up light position
  light_in_eye_space[0] = LightX;
  light_in_eye_space[1] = LightY;
  light_in_eye_space[2] = LightZ;

  // Initialize matrices for view transformation
  M3d_make_identity(vm);
  M3d_make_identity(vi);

  // Initialize the number of objects in the scene
  num_objects = 0;
	
	//////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////
	// Goal: Make this into a background plane that is texured with forest

	obtype[num_objects] = 2; // Plane

	reflectivity[num_objects] = 0; // 0% reflective
	transparency[num_objects] = 0; // 0% transparency
	refractive_index[num_objects] = 0;

	irgb[num_objects][0] = 1;
	irgb[num_objects][1] = 1;
	irgb[num_objects][2] = 1; // White color (for testing)

	Tn = 0;
	Ttypelist[Tn] = SX ; Tvlist[Tn] = 	8.5		; Tn++ ;
	Ttypelist[Tn] = SY ; Tvlist[Tn] = 	8.5		; Tn++ ;
	Ttypelist[Tn] = SZ ; Tvlist[Tn] = 	1.0		; Tn++ ;
	Ttypelist[Tn] = TZ ; Tvlist[Tn] = 	10		; Tn++ ;

	M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
	M3d_mat_mult(obmat[num_objects], vm, m);
	M3d_mat_mult(obinv[num_objects], mi, vi);

	num_objects++; // don't forget to do this

	//////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////
	// Goal: Make this a semi-transparent refractive water droplet

	obtype[num_objects] = 1 ; // sphere

	reflectivity[num_objects] = 0 ; // 0% Reflective
	transparency[num_objects] = 1.0 ; // 100% transparency
	refractive_index[num_objects] = 0.5 ; // Water = 1.33
	
	irgb[num_objects][0] = 0.400;
	irgb[num_objects][1] = 0.730;
	irgb[num_objects][2] = 0.977; // lightish blue

	Tn = 0 ;
	Ttypelist[Tn] = SX ; Tvlist[Tn] =   0.5   ; Tn++ ;
	Ttypelist[Tn] = SY ; Tvlist[Tn] =   0.6   ; Tn++ ;
	Ttypelist[Tn] = SZ ; Tvlist[Tn] =   0.4   ; Tn++ ;
	Ttypelist[Tn] = TZ ; Tvlist[Tn] =   3.0   ; Tn++ ;
	Ttypelist[Tn] = TX ; Tvlist[Tn] =  -0.4   ; Tn++ ;
	Ttypelist[Tn] = TY ; Tvlist[Tn] =  -0.8   ; Tn++ ;

	M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
	M3d_mat_mult(obmat[num_objects], vm, m) ;
	M3d_mat_mult(obinv[num_objects], mi, vi) ;

	num_objects++ ; // don't forget to do this
	
	//////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////
	// Goal: Make this a semi-transparent refractive water streak

	obtype[num_objects] = 0 ; // cylinder

	reflectivity[num_objects] = 0 ; // 0% Reflective
	transparency[num_objects] = 1.0 ; // 80% transparency
	refractive_index[num_objects] = 0.7 ; // Water = 1.33
	
	irgb[num_objects][0] = 0.400;
	irgb[num_objects][1] = 0.730;
	irgb[num_objects][2] = 0.977; // lightish blue

	Tn = 0 ;
	Ttypelist[Tn] = SX ; Tvlist[Tn] =   0.2   ; Tn++ ;
	Ttypelist[Tn] = SY ; Tvlist[Tn] =   0.5   ; Tn++ ;
	Ttypelist[Tn] = SZ ; Tvlist[Tn] =   0.1   ; Tn++ ;
	Ttypelist[Tn] = TZ ; Tvlist[Tn] =   3.0   ; Tn++ ;
	Ttypelist[Tn] = TX ; Tvlist[Tn] =   1.4   ; Tn++ ;
	Ttypelist[Tn] = TY ; Tvlist[Tn] =   1.0   ; Tn++ ;

	M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
	M3d_mat_mult(obmat[num_objects], vm, m) ;
	M3d_mat_mult(obinv[num_objects], mi, vi) ;

	num_objects++ ; // don't forget to do this
	
	//////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////
	// Goal: Make this a semi-transparent refractive water droplet

	obtype[num_objects] = 1 ; // sphere

	reflectivity[num_objects] = 0 ; // 0% Reflective
	transparency[num_objects] = 1.0 ; // 100% transparency
	refractive_index[num_objects] = 0.8 ; // Water = 1.33
	
	irgb[num_objects][0] = 0.400;
	irgb[num_objects][1] = 0.730;
	irgb[num_objects][2] = 0.977; // lightish blue

	Tn = 0 ;
	Ttypelist[Tn] = SX ; Tvlist[Tn] =   0.2   ; Tn++ ;
	Ttypelist[Tn] = SY ; Tvlist[Tn] =   0.2   ; Tn++ ;
	Ttypelist[Tn] = SZ ; Tvlist[Tn] =   0.2   ; Tn++ ;
	Ttypelist[Tn] = TZ ; Tvlist[Tn] =   3.0   ; Tn++ ;
	Ttypelist[Tn] = TX ; Tvlist[Tn] =   0.7   ; Tn++ ;
	Ttypelist[Tn] = TY ; Tvlist[Tn] =   1.0   ; Tn++ ;

	M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
	M3d_mat_mult(obmat[num_objects], vm, m) ;
	M3d_mat_mult(obinv[num_objects], mi, vi) ;

	num_objects++ ; // don't forget to do this
	
	//////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////
	// Goal: Make this a semi-transparent refractive water droplet

	obtype[num_objects] = 1 ; // sphere

	reflectivity[num_objects] = 0 ; // 0% Reflective
	transparency[num_objects] = 1.0 ; // 100% transparency
	refractive_index[num_objects] = 0.8 ; // Water = 1.33
	
	irgb[num_objects][0] = 0.400;
	irgb[num_objects][1] = 0.730;
	irgb[num_objects][2] = 0.977; // lightish blue

	Tn = 0 ;
	Ttypelist[Tn] = SX ; Tvlist[Tn] =   0.2   ; Tn++ ;
	Ttypelist[Tn] = SY ; Tvlist[Tn] =   0.2   ; Tn++ ;
	Ttypelist[Tn] = SZ ; Tvlist[Tn] =   0.2   ; Tn++ ;
	Ttypelist[Tn] = TZ ; Tvlist[Tn] =   3.0   ; Tn++ ;
	Ttypelist[Tn] = TX ; Tvlist[Tn] =   0.1   ; Tn++ ;
	Ttypelist[Tn] = TY ; Tvlist[Tn] =  -0.1   ; Tn++ ;

	M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
	M3d_mat_mult(obmat[num_objects], vm, m) ;
	M3d_mat_mult(obinv[num_objects], mi, vi) ;

	num_objects++ ; // don't forget to do this
	
	//////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////
	// Goal: Make this a semi-transparent refractive water droplet

	obtype[num_objects] = 1 ; // sphere

	reflectivity[num_objects] = 0 ; // 0% Reflective
	transparency[num_objects] = 1.0 ; // 100% transparency
	refractive_index[num_objects] = 0.8 ; // Water = 1.33
	
	irgb[num_objects][0] = 0.400;
	irgb[num_objects][1] = 0.730;
	irgb[num_objects][2] = 0.977; // lightish blue

	Tn = 0 ;
	Ttypelist[Tn] = SX ; Tvlist[Tn] =   0.2   ; Tn++ ;
	Ttypelist[Tn] = SY ; Tvlist[Tn] =   0.2   ; Tn++ ;
	Ttypelist[Tn] = SZ ; Tvlist[Tn] =   0.2   ; Tn++ ;
	Ttypelist[Tn] = TZ ; Tvlist[Tn] =   3.0   ; Tn++ ;
	Ttypelist[Tn] = TX ; Tvlist[Tn] =   0.9   ; Tn++ ;
	Ttypelist[Tn] = TY ; Tvlist[Tn] =  -1.0   ; Tn++ ;

	M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
	M3d_mat_mult(obmat[num_objects], vm, m) ;
	M3d_mat_mult(obinv[num_objects], mi, vi) ;

	num_objects++ ; // don't forget to do this
	
	//////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////
	// Goal: Make this a semi-transparent refractive water droplet

	obtype[num_objects] = 1 ; // sphere

	reflectivity[num_objects] = 0 ; // 0% Reflective
	transparency[num_objects] = 1.0 ; // 100% transparency
	refractive_index[num_objects] = 0.9; // Water = 1.33
	
	irgb[num_objects][0] = 0.400;
	irgb[num_objects][1] = 0.730;
	irgb[num_objects][2] = 0.977; // lightish blue

	Tn = 0 ;
	Ttypelist[Tn] = SX ; Tvlist[Tn] =   0.3   ; Tn++ ;
	Ttypelist[Tn] = SY ; Tvlist[Tn] =   0.4   ; Tn++ ;
	Ttypelist[Tn] = SZ ; Tvlist[Tn] =   0.3   ; Tn++ ;
	Ttypelist[Tn] = TZ ; Tvlist[Tn] =   3.0   ; Tn++ ;
	Ttypelist[Tn] = TX ; Tvlist[Tn] =  -1.5   ; Tn++ ;
	Ttypelist[Tn] = TY ; Tvlist[Tn] =   1.6   ; Tn++ ;

	M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
	M3d_mat_mult(obmat[num_objects], vm, m) ;
	M3d_mat_mult(obinv[num_objects], mi, vi) ;

	num_objects++ ; // don't forget to do this
	
	//////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////
	// Goal: Make this a semi-transparent refractive water droplet

	obtype[num_objects] = 1 ; // sphere

	reflectivity[num_objects] = 0 ; // 0% Reflective
	transparency[num_objects] = 1.0 ; // 100% transparency
	refractive_index[num_objects] = 0.8 ; // Water = 1.33
	
	irgb[num_objects][0] = 0.400;
	irgb[num_objects][1] = 0.730;
	irgb[num_objects][2] = 0.977; // lightish blue

	Tn = 0 ;
	Ttypelist[Tn] = SX ; Tvlist[Tn] =   0.2   ; Tn++ ;
	Ttypelist[Tn] = SY ; Tvlist[Tn] =   0.2   ; Tn++ ;
	Ttypelist[Tn] = SZ ; Tvlist[Tn] =   0.2   ; Tn++ ;
	Ttypelist[Tn] = TZ ; Tvlist[Tn] =   3.0   ; Tn++ ;
	Ttypelist[Tn] = TX ; Tvlist[Tn] =  -1.3   ; Tn++ ;
	Ttypelist[Tn] = TY ; Tvlist[Tn] =  -0.3   ; Tn++ ;

	M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
	M3d_mat_mult(obmat[num_objects], vm, m) ;
	M3d_mat_mult(obinv[num_objects], mi, vi) ;

	num_objects++ ; // don't forget to do this
	
	//////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////
	// Goal: Make this a semi-transparent refractive water droplet

	obtype[num_objects] = 1 ; // sphere

	reflectivity[num_objects] = 0 ; // 0% Reflective
	transparency[num_objects] = 1.0 ; // 100% transparency
	refractive_index[num_objects] = 0.7 ; // Water = 1.33
	
	irgb[num_objects][0] = 0.400;
	irgb[num_objects][1] = 0.730;
	irgb[num_objects][2] = 0.977; // lightish blue

	Tn = 0 ;
	Ttypelist[Tn] = SX ; Tvlist[Tn] =   0.3   ; Tn++ ;
	Ttypelist[Tn] = SY ; Tvlist[Tn] =   0.3   ; Tn++ ;
	Ttypelist[Tn] = SZ ; Tvlist[Tn] =   0.3   ; Tn++ ;
	Ttypelist[Tn] = TZ ; Tvlist[Tn] =   3.0   ; Tn++ ;
	Ttypelist[Tn] = TX ; Tvlist[Tn] =  -0.4   ; Tn++ ;
	Ttypelist[Tn] = TY ; Tvlist[Tn] =   0.7   ; Tn++ ;

	M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
	M3d_mat_mult(obmat[num_objects], vm, m) ;
	M3d_mat_mult(obinv[num_objects], mi, vi) ;

	num_objects++ ; // don't forget to do this

	//////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////

	// Clear the graphics window
  G_rgb(0, 0, 0);
  G_clear();

  // Loop to render the scene and handle user input
  int key = 0 ;
  while (key != 'q') {
		// Initialize ray source and tip
    Rsource[0] = 0; Rsource[1] = 0; Rsource[2] = 0;
    Rtip[0] = -H; Rtip[1] = -H; Rtip[2] = 1;

		int obnum, insct, x, y, i ;
		double Rsource_new[3], Rtip_new[3], intersect[3] ;
		double ytip_save, xtip_save  ;

		ytip_save = Rtip[0] ;
		xtip_save = Rtip[0] ;

		double normal[3] ;
		double reflection[3] ;

		int objs_reflected[100] ;
		int num_intersections ;

		double proj_u, proj_v, u, v ;
		
		// Object texture mapping
		int id ;
		int dim[2] ;
		double width, height ;

		id = init_xwd_map_from_file("../XWDhelp/Change_xwd_styles/MistyForestJ.xwd") ;
		if(id == -1) {
			printf("failure: can't init map\n") ;
		}
		i = get_xwd_map_dimensions(id, dim) ;
		if(i == -1) {
			printf("failure: can't dimensions\n") ;
		}
		width = dim[0]; height = dim[1] ;

		// double for loop moves through each pt of "film"
		for(x = -400; x < 400; x++) {
			for(y = -400; y < 400; y++) {
				Rtip[0] = x*(H/400);
				Rtip[1] = y*(H/400);
				Rtip[2] = 1;

				Rsource[0] =  0 ;  Rsource[1] =  0 ;  Rsource[2] = 0 ;

				insct = 1;
				num_intersections = 0;

				while (insct == 1 && num_intersections < REF_DEPTH) {

					T = 100000000;
					OB = -1;

					for (obnum = 0; obnum < num_objects; obnum++) {
						M3d_mat_mult_pt(Rsource_new, obinv[obnum], Rsource);
						M3d_mat_mult_pt(Rtip_new, obinv[obnum], Rtip);

						insct = find_intersection(Rsource_new, Rtip_new, intersect, obnum);

					}

					if(OB > -1) {

						// find unit normal in obj space
						i = find_normal(intersect, normal, obinv[OB], OB);
						
						// texture mapping on plane...
						if (OB == 0) {
							double x = intersect[0];
							double y = intersect[1];
							double z = intersect[2];

							// Calculate texture coordinates based on intersection point
							u = (x + 1.0) / 2.0;
							v = (y + 1.0) / 2.0;

							// Scale texture coordinates to image dimensions
							proj_u = u * width;
							proj_v = v * height;

							// Get the color from the texture map
							i = get_xwd_map_color(id, proj_u, proj_v, irgb[OB]);
						}

						// send intersection to obj space and save it for later
						M3d_mat_mult_pt(intersect, obmat[OB], intersect);

						objs_reflected[num_intersections] = OB;
						num_intersections++;
						
						if (transparency[OB] > 0) {
		          // Refraction
							double eta = refractive_index[OB];
							double I[3];
							I[0] = intersect[0] - Rsource[0];
							I[1] = intersect[1] - Rsource[1];
							I[2] = intersect[2] - Rsource[2];

							// Normalize the incident ray direction
							double mag = magnitude(I);
							I[0] /= mag; I[1] /= mag; I[2] /= mag;

							double R[3];
							refract(I, normal, eta, R);

							// Update source and tip for refracted ray
							Rsource[0] = intersect[0] + 0.0001 * R[0];
							Rsource[1] = intersect[1] + 0.0001 * R[1];
							Rsource[2] = intersect[2] + 0.0001 * R[2];
							Rtip[0]    = intersect[0] + R[0];
							Rtip[1]    = intersect[1] + R[1];
							Rtip[2]    = intersect[2] + R[2];

							// Update the reflection vector with the refracted ray direction
							reflection[0] = R[0];
							reflection[1] = R[1];
							reflection[2] = R[2];
												
        		} else {
						// find unit reflection in obj space
						i = find_reflection(intersect, Rsource, normal, reflection);

						// light model
						i = Light_Model(irgb[OB], Rsource, intersect, normal, lrgb[OB]);

						// make intersection pt the new source; reflection pt the new tip
						Rsource[0] = intersect[0] + 0.0001 * reflection[0];
						Rsource[1] = intersect[1] + 0.0001 * reflection[1];
						Rsource[2] = intersect[2] + 0.0001 * reflection[2];
						Rtip[0]    = intersect[0] + reflection[0];
						Rtip[1]    = intersect[1] + reflection[1];
						Rtip[2]    = intersect[2] + reflection[2];

						}
					}
				} // end while

				if(num_intersections > 0) {
					obnum = objs_reflected[num_intersections - 1];
					ARGB[0] = lrgb[obnum][0];
					ARGB[1] = lrgb[obnum][1];
					ARGB[2] = lrgb[obnum][2];

					for(i = num_intersections - 1; i >= 0; i--){
						obnum = objs_reflected[i];
						
						if (transparency[obnum] > 0.0) {
		          // Refraction contribution
		          ARGB[0] = ARGB[0] * transparency[obnum] + (1 - transparency[obnum]) * lrgb[obnum][0];
		          ARGB[1] = ARGB[1] * transparency[obnum] + (1 - transparency[obnum]) * lrgb[obnum][1];
		          ARGB[2] = ARGB[2] * transparency[obnum] + (1 - transparency[obnum]) * lrgb[obnum][2];
        		} else {
        			// Reflection contribution
							ARGB[0] *= reflectivity[obnum];
							ARGB[0] += (1 - reflectivity[obnum])*lrgb[obnum][0];

							ARGB[1] *= reflectivity[obnum];
							ARGB[1] += (1 - reflectivity[obnum])*lrgb[obnum][1];

							ARGB[2] *= reflectivity[obnum];
							ARGB[2] += (1 - reflectivity[obnum])*lrgb[obnum][2];
						}
					}


				} else {
					ARGB[0] = 0; ARGB[1] = 0; ARGB[2] = 0;
				} 

				G_rgb(ARGB[0], ARGB[1], ARGB[2]);
				G_point(x + 400, y + 400) ;

			} // end for y
		} // end for x
		
		G_rgb(1,1,1) ; G_draw_string("'s' to save image, 'q' to quit", 50,50) ;
		
		// Saves as image file if 's' key pressed
		if (key == 's') {
			G_save_image_to_file("Final_Presentation.xwd") ;
			G_rgb(1,1,1) ; G_draw_string("Image saved as Final_Presentation.xwd", 560,750) ;
			key = G_wait_key() ;
			
			// Exit program if 'q' is pressed after saving
			if (key == 'q') {
				break ;
			}
		}
		
		key = G_wait_key();
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////

// Main function
int main() {
	G_init_graphics(800, 800);
	test01() ;
}
