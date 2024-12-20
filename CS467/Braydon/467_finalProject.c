#include "FPToolkit.c"
#include "M3d_matrix_tools.c"
#include "xwd_tools_03.c"


//double x[10000],y[10000],z[10000];
double zbuffer[800][800];
double eye[3],coi[3],up[3];
int np;

//matricies for spheres of different colors
double EARTH[4][4],SUN[4][4],STARS[4][4];
double EARTHi[4][4],SUNi[4][4],STARSi[4][4];
double EARTHobj[4][4],SUNobj[4][4],STARSobj[4][4];

// To support the light model :
double light_in_eye_space[3] ;
double AMBIENT      = 0.2 ;
double MAX_DIFFUSE  = 0.5 ;
double SPECPOW      = 50 ;
double inherent_rgb[3];

double specular;
double diffuse;
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

  diffuse =0 ;
  if (NdotL <= 0.0) { diffuse = 0.0 ; }
  else { diffuse = MAX_DIFFUSE*NdotL ; }

  specular=0;
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


int idEARTH,idSUN,idSTARS;
void getColor(double u,double v,double uhi,double vhi,double ulo,double vlo,double idA) {
	double e;
	int d[2],x,y;
	double width,height;
	double rgb[3];
	
  	e = get_xwd_map_dimensions(idA, d) ;
  	if (e == -1) { printf("failure\n") ;  exit(0) ; }
  	width = d[0] ; height = d[1];
  	
  	
  	x = width/(uhi-ulo)*(u-ulo);
  	y = height/(vhi-vlo)*(v-vlo);

	if (x < 0) printf("oops\n") ;
	
  	get_xwd_map_color(idA, x,y,rgb) ; // returns -1 on error, 1 if ok
  	
  	inherent_rgb[0]=rgb[0];
  	inherent_rgb[1]=rgb[1];
  	inherent_rgb[2]=rgb[2];
	
}


void plot (double (*f)(double u, double v), 
		  double (*g)(double u, double v), 
		  double (*l)(double u, double v), 
		  double ulo, double uhi,
		  double vlo, double vhi,
		  double MOVEMENT[4][4]) {
	
	double u,v;
	int j,k,i;
	double xbb,ybb,x,y,z;
	double halfangle;
	double H;
	double point[3],point2[3],point3[3],rgb[3];
	double pointdist,temp;
	int a,b;
	
	halfangle=40;
	H=tan(halfangle*(M_PI/180));
	pointdist=0;
	
	for (u = ulo; u <= uhi ; u+=0.0005) {
    	for(v=vlo; v<= vhi; v+=0.0005) {
   		point[0]=f(u,v);
   		point[1]=g(u,v);
   		point[2]=l(u,v);
		if(MOVEMENT==EARTHobj) {
			getColor(u,v,uhi,vhi,ulo,vlo,idEARTH);
		} else if(MOVEMENT==SUNobj) {
			getColor(u,v,uhi,vhi,ulo,vlo,idSUN);
		} else if(MOVEMENT==STARSobj) {
			getColor(u,v,uhi,vhi,ulo,vlo,idSTARS);
		}
		
   		point2[0]=f(u+0.01,v);
   		point2[1]=g(u+0.01,v);
   		point2[2]=l(u+0.01,v);
   		
   		point3[0]=f(u,v+0.01);
   		point3[1]=g(u,v+0.01);
   		point3[2]=l(u,v+0.01);
   		
   		M3d_mat_mult_pt(point,MOVEMENT,point);
   		M3d_mat_mult_pt(point2,MOVEMENT,point2);
   		M3d_mat_mult_pt(point3,MOVEMENT,point3);
   		x=point[0];
   		y=point[1];
   		z=point[2];

		double hither = 0.01 ;
                if (z < hither) { continue; }
                if (fabs(y/z) > H) { continue ; }
                if (fabs(x/z) > H) { continue ; }		
		
   		xbb=((400/H)*(x/z)+400);
			ybb=((400/H)*(y/z)+400);
			
		//zbuffer
		a=xbb;
		b=ybb;
		pointdist=z;
		if (a < 0 || a >= 800 || b < 0 || b >= 800) {
		  printf("a = %d  b = %d\n",a,b) ;
		  continue ;
		}
		
		if(pointdist<zbuffer[a][b]) {
		
		if(MOVEMENT==STARSobj) {SPECPOW=1000000000;} else {SPECPOW=30;}
			zbuffer[a][b]=pointdist;
			light_model (inherent_rgb,point,point2,point3,rgb);
			G_rgb(rgb[0], rgb[1], rgb[2]);
			G_point(xbb,ybb);
		} else if (pointdist>zbuffer[a][b]) {
			zbuffer[a][b]=zbuffer[a][b];
			}
    	}		
  	}
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

double square_x(double u) {
	double sgn,result;
	
	if(u>=0&&u<=M_PI/2) {
		sgn=-1;
		result = sgn*(cos(u))*cos(u); 
	} else if (u>=M_PI/2&&u<=M_PI){
		sgn=-1;
		result = sgn*(cos(u))*cos(u);  
	} else if (u>=M_PI&&u<=(3*M_PI)/2){
		sgn=1;
		result = sgn*(cos(u))*cos(u);  
	} else if (u>=(3*M_PI)/2&&u<=2*M_PI){
		sgn=1;
		result = sgn*(cos(u))*cos(u);  
	}

	return result;
	
}

double square_y(double u) {
	double sgn,result;
	
	if(u>=0&&u<=M_PI/2) {
		sgn=1;
		result = sgn*(sin(u))*sin(u); 
	} else if (u>=M_PI/2&&u<=M_PI){
		sgn=-1;
		result = sgn*(sin(u))*sin(u);  
	} else if (u>=M_PI&&u<=(3*M_PI)/2){
		sgn=-1;
		result = sgn*(sin(u))*sin(u);  
	} else if (u>=(3*M_PI)/2&&u<=2*M_PI){
		sgn=1;
		result = sgn*(sin(u))*sin(u);  
	}
	
	return result;
}


double torus_x(double u, double v) {
	double x = square_x(u)*cos(45)+square_y(u)*sin(45);
	return x;
}

double torus_y(double u, double v) {
	double y = square_y(u)*cos(45) - square_x(u)*sin(45);
	return cos(v)*(y+6);
}

double torus_z(double u, double v) {
double y = square_y(u)*cos(45) - square_x(u)*sin(45);
	return sin(v)*(y+6);
}

double disk_x(double u, double v) {
	double x = square_x(u)*cos(45)+square_y(u)*sin(45);
	return x;
}

double disk_y(double u, double v) {
	double y = square_y(u)*cos(45) - square_x(u)*sin(45);
	return cos(v)*(y+0.6);
}

double disk_z(double u, double v) {
double y = square_y(u)*cos(45) - square_x(u)*sin(45);
	return sin(v)*(y+0.6);
}

double cyl_x(double u, double v) {
	return cos(u);
}

double cyl_y(double u, double v) {
	return v;
}

double cyl_z(double u, double v) {
	return sin(u);
}


void init_zbuffer()
{
  int j,k ;
	//initialize zbuffer
	for(j=0;j<800;j++) {
		for(k=0;k<800;k++) {
			zbuffer[j][k]=100000;
		}
	}

}



int main() {
	int k,j;
	double u;
	double normal[2],mag;

	
	light_in_eye_space[0] = 1.4;
  light_in_eye_space[1] =  0.2;
  light_in_eye_space[2] = -0.7;
  
  AMBIENT = 0.2 ;
  MAX_DIFFUSE = 0.5 ;
  SPECPOW = 30;
	
	//initialize graphics pannel
	G_init_graphics(800,800);
	G_rgb(0,0,0);
	G_clear();
	
	init_zbuffer();
	
	double Tvlist[100];
	int Ttypelist[100],Tn;
	double V[4][4],Vi[4][4];
	double E[4][4],Ei[4][4];
	char idName[100];
	
	
	double t;
	int q,fnum;
	
	/*//VIEW MATRIX
	t=0.01*fnum;
  	eye[0] = 0; 
	eye[1] = 0; 
	eye[2] = 9-t; 
	
    coi[0] =  0;
    coi[1] =  0; 
    coi[2] =  0;

    up[0]  = eye[0] ; 
    up[1]  = eye[1] + 1 ;
    up[2]  = eye[2] ; 
	M3d_view(E,Ei,eye,coi,up);*/
	
	
	
	//movement sequence for MERCURY.
	Tn = 0; 
	Ttypelist[Tn] = RY ; Tvlist[Tn] =   -80; Tn++;
	Ttypelist[Tn] = SX ; Tvlist[Tn] =   .034; Tn++;
	Ttypelist[Tn] = SY ; Tvlist[Tn] =   .034; Tn++;
	Ttypelist[Tn] = SZ ; Tvlist[Tn] =   .034; Tn++;
  Ttypelist[Tn] = TY ; Tvlist[Tn] =   0; Tn++;
  Ttypelist[Tn] = TX ; Tvlist[Tn] =   0; Tn++;
  Ttypelist[Tn] = TZ ; Tvlist[Tn] =   2; Tn++;
	M3d_make_movement_sequence_matrix(EARTH,EARTHi,Tn,Ttypelist,Tvlist);
	
	//movement sequence for SUN.
	Tn = 0; 
	Ttypelist[Tn] = SX ; Tvlist[Tn] =   2; Tn++;
	Ttypelist[Tn] = SY ; Tvlist[Tn] =   2; Tn++;
	Ttypelist[Tn] = SZ ; Tvlist[Tn] =   2; Tn++;
	Ttypelist[Tn] = TX ; Tvlist[Tn] =   -2; Tn++;
	
	M3d_make_movement_sequence_matrix(SUN,SUNi,Tn,Ttypelist,Tvlist);
	
	//movement sequence for STARS.
	Tn = 0; 
	Ttypelist[Tn] = SX ; Tvlist[Tn] =   100; Tn++;
	Ttypelist[Tn] = SY ; Tvlist[Tn] =   100; Tn++;
	Ttypelist[Tn] = SZ ; Tvlist[Tn] =   100; Tn++;
	Ttypelist[Tn] = RY ; Tvlist[Tn] =   75; Tn++;
	
	M3d_make_movement_sequence_matrix(STARS,STARSi,Tn,Ttypelist,Tvlist);
	
	idEARTH = init_xwd_map_from_file ("mercuryB.xwd") ;// returns -1 on error, 1 if ok
  	if (idEARTH == -1) { printf("failure\n") ;  exit(0) ; }
	
	idSUN = init_xwd_map_from_file ("sunB.xwd") ;// returns -1 on error, 1 if ok
  	if (idSUN == -1) { printf("failure\n") ;  exit(0) ; }
  	
  	idSTARS = init_xwd_map_from_file ("starsB.xwd") ;// returns -1 on error, 1 if ok
  	if (idSUN == -1) { printf("failure\n") ;  exit(0) ; }
	
	//for movie, while less than desired number of frames
	fnum=0;
	q=45;
	double l=0;
	int count=0;
	char filename[100] = "FinalMovie/FinalMovie.xwd";
	//while(q != 'q') {
	while (q < 346) { //auto control: 346 for first movie - got to frame 45 start at this point next time
	  	init_zbuffer();

	  	t = 0.005*fnum;
	  	l = 0.0005*fnum;

			eye[0] = sin(l); 
			eye[1] = 0; 
			eye[2] = 3-t; 

		//printf("t = %lf   eye = %lf %lf %lf\n",t, eye[0],eye[1],eye[2]) ;

		coi[0] =  0;
		coi[1] =  0; 
		coi[2] =  0;

		up[0]  = eye[0]; 
		up[1]  = eye[1] + 1;
		up[2]  = eye[2]; 
		
		double rotate[4][4],rotatei[4][4];
		Tn=0;
		Ttypelist[Tn] = RY ; Tvlist[Tn] =   20; Tn++;
		M3d_make_movement_sequence_matrix(rotate,rotatei,Tn,Ttypelist,Tvlist);
		M3d_mat_mult(EARTHobj,  rotate,EARTH);

	  	// move ALL objects from WORLD SPACE into EYE SPACE:  
	  	M3d_view (E, Ei,eye,coi,up);
		
	  	M3d_mat_mult(EARTHobj,  E,EARTH);
	  	M3d_mat_mult(SUNobj,  E,SUN);
	  	M3d_mat_mult(STARSobj,  E,STARS);
	  	
	  	
	  
	  	G_rgb(0,0,0) ; 
	  	G_clear() ;
	  
		//draw MARS
		
		plot(sphere_x,sphere_y,sphere_z,-1*M_PI,M_PI,-1*M_PI/2,M_PI/2,EARTHobj);
		
		
		//draw SUN		
  		plot(sphere_x,sphere_y,sphere_z,-1*M_PI,M_PI,-1*M_PI/2,M_PI/2,SUNobj);
  	
  		//draw STARS
  		plot(sphere_x,sphere_y,sphere_z,-1*M_PI,M_PI,-1*M_PI/2,M_PI/2,STARSobj);
  	
		
	 	G_display_image();
		usleep(1);
	
		//manual control
	  	//q = G_wait_key(); 
	  	
	   sprintf(filename,"FinalMovie/FinalMovie%04d.xwd",q);
	   G_save_image_to_file(filename);
	  	
	   fnum++;
	   q++;
	   printf("%d\n",q);
  } //end while for camera movement
  
	G_wait_key();
	
}

