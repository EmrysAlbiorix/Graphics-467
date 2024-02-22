#include "FPToolkit.c"
#include "M3d_matrix_tools.c"


//double x[10000],y[10000],z[10000];
double eye[3],coi[0],up[3];
int np;

//matricies for spheres of different colors
double RED[4][4],YELLOW[4][4],BLUE[4][4];

// To support the light model :
double light_in_eye_space[3] ;
double AMBIENT      = 0.2 ;
double MAX_DIFFUSE  = 0.5 ;
double SPECPOW      = 50 ;
double inherent_rgb[3];



int Light_Model (double irgb[3],
                 double s[3],
                 double p[3],
                 double n[3],
                 double argb[3])
// s,p,n in eyespace

// irgb == inherent color of object (input to this function)
// s = location of start of ray (probably the eye)
// p = point on object (input to this function)
// n = normal to the object at p (input to this function)
// argb == actual color of object (output of this function)
// globals : AMBIENT, MAX_DIFFUSE, SPECPOW, light_in_eye_space[3]

// return 1 if successful, 0 if error
{

  double len ;
  double N[3] ; 
  len = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]) ;
  if (len == 0) return 0 ;
  N[0] = n[0]/len ;  N[1] = n[1]/len ;  N[2] = n[2]/len ;

  double E[3] ;
  E[0] = s[0] - p[0] ; 
  E[1] = s[1] - p[1] ; 
  E[2] = s[2] - p[2] ; 
  len = sqrt(E[0]*E[0] + E[1]*E[1] + E[2]*E[2]) ;
  if (len == 0) return 0 ;
  E[0] /= len ;  E[1] /= len ;  E[2] /= len ;
  double NdotE = N[0]*E[0] + N[1]*E[1] + N[2]*E[2] ;

  double L[3] ;
  L[0] = light_in_eye_space[0] - p[0] ; 
  L[1] = light_in_eye_space[1] - p[1] ; 
  L[2] = light_in_eye_space[2] - p[2] ; 
  len = sqrt(L[0]*L[0] + L[1]*L[1] + L[2]*L[2]) ;
  if (len == 0) return 0 ;
  L[0] /= len ;  L[1] /= len ;  L[2] /= len ;
  double NdotL = N[0]*L[0] + N[1]*L[1] + N[2]*L[2] ;





  double max_ambient_and_diffuse = AMBIENT + MAX_DIFFUSE ;
     // this needs to occur BEFORE you possibly jump to LLL below




  double intensity ;
  if (NdotL*NdotE < 0) {
    // eye and light are on opposite sides of polygon
    intensity = AMBIENT ; 
    goto LLL ;
  } else if ((NdotL < 0) && (NdotE < 0)) {
    // eye and light on same side but normal pointing "wrong" way
    N[0] *= (-1.0) ;    N[1] *= (-1.0) ;    N[2] *= (-1.0) ; 
    NdotL *= (-1.0) ;
    NdotE *= (-1.0) ;   // don't use NdotE below, probably should eliminate this
  }


  // ignore Blinn's variant
  double R[3] ; // Reflection vector of incoming light
  R[0] = 2*NdotL*N[0] - L[0] ;
  R[1] = 2*NdotL*N[1] - L[1] ;
  R[2] = 2*NdotL*N[2] - L[2] ;

  double EdotR = E[0]*R[0] + E[1]*R[1] + E[2]*R[2] ;

  double diffuse ;
  if (NdotL <= 0.0) { diffuse = 0.0 ; }
  else { diffuse = MAX_DIFFUSE*NdotL ; }

  double specular ;
  if (EdotR <= 0.0) { specular = 0.0 ; }
  else { specular = (1.0 - max_ambient_and_diffuse)*pow(EdotR,SPECPOW) ;}

  // printf("%lf %lf\n",diffuse,specular) ;
  intensity = AMBIENT + diffuse + specular ;



 LLL : ;

  double f,g ;
  if (intensity <= max_ambient_and_diffuse) {
    f = intensity / max_ambient_and_diffuse ;
    argb[0] = f * irgb[0] ;
    argb[1] = f * irgb[1] ;
    argb[2] = f * irgb[2] ;
  } else {
    f = (intensity - max_ambient_and_diffuse) / 
                           (1.0 - max_ambient_and_diffuse) ;
    g = 1.0 - f ;
    argb[0] = g * irgb[0] + f ;
    argb[1] = g * irgb[1] + f ;
    argb[2] = g * irgb[2] + f ;
  }

  return 1 ;
}






void light_model (double irgb[3],
                  double *point, 
                  double *point2,
                  double *point3,
                  double argb[3])
// irgb == inherent color of object (input to this function)
// xx[],yy[],zz[] are points in the polygon
// argb == actual color of object (output of this function)
{
  double Eye[3] ;
  Eye[0] = 0 ; Eye[1] = 0 ; Eye[2] = 0 ; 

  double P[3]  ;
  P[0] = point[0] ;  P[1] = point[1] ;  P[2] = point[2] ;

  double a[3] ;
  //a[0] = xx[1] - xx[0] ;  a[1] = yy[1] - yy[0] ;  a[2] = zz[1] - zz[0] ;
  a[0] = point2[0] - point[0];  a[1] = point2[1] - point[1];  a[2] = point2[2] - point[2];

  double b[3] ;
  //b[0] = xx[2] - xx[0] ;  b[1] = yy[2] - yy[0] ;  b[2] = zz[2] - zz[0] ;
  b[0] = point3[0] - point[0] ;  b[1] = point3[1] - point[1] ;  b[2] = point3[2] - point[2] ;
 
  double N[3] ;
  M3d_x_product (N, a,b) ;

  Light_Model (irgb, Eye, P, N, argb) ;
}

double dist(double *point) {
	double x,y,z;
	
	x=point[0];
	y=point[1];
	z=point[2];
	
	//sqrt(xbar-x)^2+(ybar-y)^2+(zbar-z)^2
	
	return sqrt(x*x+y*y+z*z);


}


int plot (double (*f)(double u, double v), double (*g)(double u, double v), double (*l)(double u, double v), double ulo, double uhi,double vlo, double vhi) {
	
	double Tvlist[100];
	int Ttypelist[100],Tn;
	double V[4][4],Vi[4][4];
	double E[4][4],Ei[4][4];
	
	//movement sequence for RED sphere.
	Tn = 0; 
	//Ttypelist[Tn] = SY ; Tvlist[Tn] =   0.7 ; Tn++;
  	Ttypelist[Tn] = TZ ; Tvlist[Tn] =   5; Tn++;
	M3d_make_movement_sequence_matrix(V,Vi,Tn,Ttypelist,Tvlist);
	//M3d_view(E,Ei,eye,coi,up);
	
	
	double u,v;
	int i=0;
	double xbb,ybb,x,y,z;
	double halfangle;
	double H;
	double point[3],point2[3],point3[3],rgb[3],zbuffer[800][800];
	double pointdist,temp;
	int a,b;
	
	halfangle=40;
	H=tan(halfangle*(M_PI/180));
	pointdist=0;
	
	for (u = ulo; u <= uhi ; u+=0.005) {
    	for(v=vlo; v<= vhi; v+=0.005) {
   		point[0]=f(u,v);
   		point[1]=g(u,v);
   		point[2]=l(u,v);
   		
   		point2[0]=f(u+0.01,v);
   		point2[1]=g(u+0.01,v);
   		point2[2]=l(u+0.01,v);
   		
   		point3[0]=f(u,v+0.01);
   		point3[1]=g(u,v+0.01);
   		point3[2]=l(u,v+0.01);
   		
   		M3d_mat_mult_pt(point,V,point);
   		M3d_mat_mult_pt(point2,V,point2);
   		M3d_mat_mult_pt(point3,V,point3);
   		x=point[0];
   		y=point[1];
   		z=point[2];
   		xbb=((400/H)*(x/z)+400);
		ybb=((400/H)*(y/z)+400);
		
		
		temp=dist(point);
		if(pointdist>temp) {
			pointdist=pointdist;
		} else {
			pointdist=temp;
		}
		
		a=x/z;
		b=y/z;
		zbuffer[a][b]=pointdist;
		
		light_model (inherent_rgb,point,point2,point3,rgb);
		G_rgb(rgb[0], rgb[1], rgb[2]);
		G_point(xbb,ybb);
   		i++;
    	}		
  	}
  	return i;
}

double sphere_x(double u,double v) {
	double r = cos(v);
	return r*cos(u);
}

double sphere_y(double u,double v) {
	return sin(v);
}

double sphere_z(double u,double v) {
	double r = cos(v);
	return r*sin(u);
}

int main() {
	double u;
	double normal[2],mag;
	light_in_eye_space[0] =  10 ;
  	light_in_eye_space[1] =  10 ;
  	light_in_eye_space[2] =  0 ;
  
  	AMBIENT = 0.2 ;
  	MAX_DIFFUSE = 0.5 ;
  	SPECPOW = 30 ;
  	
	
	//initialize graphics pannel
	G_init_graphics(800,800);
	G_rgb(0.1,0.1,0.1);
	G_clear();
	
	inherent_rgb[0]=1;
	inherent_rgb[1]=0;
	inherent_rgb[2]=0;
	
	//initialize points for sphere
	np=plot(sphere_x,sphere_y,sphere_z,-1*M_PI,M_PI,-1*M_PI/2,M_PI/2);
	
	
	//wait for user to click key and end program.
	G_wait_key();
}

