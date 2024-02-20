#include "Tools/FPToolkit.c"
#include "Tools/M3d_matrix_tools.c"


double x[1000],y[1000],z[1000];
double xv[1000],yv[1000];
int np;


int plot (double (*f)(double angle), double (*g)(double angle), double lo, double hi) {
	double angle;
	int i=0;
	
	for (angle = lo ; angle <= hi ; angle+=0.01) {
    	// draw a point
   		G_rgb (1.0, 0.0, 0.0) ; // red
   		x[i]=f(angle);
   		y[i]=g(angle);
   		xv[i]=f(angle+0.001);
   		yv[i]=g(angle+0.001);
   		i++;
  	}
  	return i;
}

void initGraphics() {
	G_init_graphics(800,800);
	G_rgb(0.1,0.1,0.1);
	G_clear();
}


double sum4_x(double u) {
	return u;
}

double sum4_y(double u) {
	double result;
	
	result = pow(1 - u*u*u*u, 0.25);
	
	return result;
}

double square1_x(double u) {
	
	double result,w;
	
	if(u>=0&&u<=1) {
		w=u;
		result=1-w;
	} else if (u>=1&&u<=2) {
		w=u-1;
		result = -w;
	} else if (u>=2&&u<=3) {
		w=u-2;
		result=w-1;
	} else if (u>=3&&u<=4) {
		w=u-3;
		result=w;
	}

	return result;
}

double square1_y(double u) {
	
	double result,w;
	
	if(u>=0&&u<=1) {
		w=u;
		result=w;
	} else if (u>=1&&u<=2) {
		w=u-1;
		result = 1-w;
	} else if (u>=2&&u<=3) {
		w=u-2;
		result=-w;
	} else if (u>=3&&u<=4) {
		w=u-3;
		result=w-1;
	}
	
	return result;
}


double square2_x(double u) {
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

double square2_y(double u) {
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

double ast_x(double u) {
	double sgn,result;
	
	if(u>=0&&u<=M_PI/2) {
		sgn=-1;
		result = sgn*(cos(u))*cos(u)*cos(u)*cos(u); 
	} else if (u>=M_PI/2&&u<=M_PI){
		sgn=-1;
		result = sgn*(cos(u))*cos(u)*cos(u)*cos(u); 
	} else if (u>=M_PI&&u<=(3*M_PI)/2){
		sgn=1;
		result = sgn*(cos(u))*cos(u)*cos(u)*cos(u);  
	} else if (u>=(3*M_PI)/2&&u<=2*M_PI){
		sgn=1;
		result = sgn*(cos(u))*cos(u)*cos(u)*cos(u);  
	}

	return result;
}

double ast_y(double u) {
	double sgn,result;
	
	if(u>=0&&u<=M_PI/2) {
		sgn=1;
		result = sgn*(sin(u))*sin(u)*sin(u)*sin(u); 
	} else if (u>=M_PI/2&&u<=M_PI){
		sgn=-1;
		result = sgn*(sin(u))*sin(u)*sin(u)*sin(u);
	} else if (u>=M_PI&&u<=(3*M_PI)/2){
		sgn=-1;
		result = sgn*(sin(u))*sin(u)*sin(u)*sin(u);
	} else if (u>=(3*M_PI)/2&&u<=2*M_PI){
		sgn=1;
		result = sgn*(sin(u))*sin(u)*sin(u)*sin(u);
	}
	
	return result;

}

double parab_y(double u) {
	return u*u;
}

double lemon_x(double u) {
	
	return cos(u)*cos(u)*cos(u);

}

double lemon_y(double u){
	return sin(u);
}

int main() {
	double u;
	double Tvlist[100];
	double v[4][4],vi[4][4];
	int Ttypelist[100],Tn;
	double normal[2];
	
	//initialize graphics pannel
	initGraphics();
	
	//-----------------circle-------------------------------------
	//initialize points for shape
	np=plot(cos,sin,0.25*M_PI, 1.5*M_PI);
	
	
	
	//circle
	Tn = 0 ; 
  	Ttypelist[Tn] = SX ; Tvlist[Tn] =   50.0 ; Tn++ ;
  	Ttypelist[Tn] = SY ; Tvlist[Tn] =  100.0 ; Tn++ ;
  	Ttypelist[Tn] = TX ; Tvlist[Tn] =  300.0 ; Tn++ ;
  	Ttypelist[Tn] = TY ; Tvlist[Tn] =  500.0 ; Tn++ ;
  	
  	M3d_make_movement_sequence_matrix(v,vi,Tn,Ttypelist,Tvlist);
	M3d_mat_mult_points(x,y,z,v,x,y,z,np);
	
	//M3d_x_product (normal, v1, v2);
	
	//draw points
	for(int i=0;i<np;i++) {
		G_rgb (1.0, 0.0, 0.0); // red
   		G_point (x[i], y[i]); // hard to see
   		normal[0]= (yv[i]-y[i]);
   		normal[1]=-1*(xv[i]-x[i]);
   		
   		G_line(x[i],y[i],normal[0],normal[1]);
	}
	
	
	
	
	
	//---------------------sum 4---------------------------------
	//initialize points for shape
	np=plot(sum4_x,sum4_y,-1,1);
	
	//sum4
	Tn = 0 ; 
  	Ttypelist[Tn] = SX ; Tvlist[Tn] =   30.0 ; Tn++ ;
  	Ttypelist[Tn] = SY ; Tvlist[Tn] =   30.0 ; Tn++ ;
  	Ttypelist[Tn] = TX ; Tvlist[Tn] =  250.0 ; Tn++ ;
  	Ttypelist[Tn] = TY ; Tvlist[Tn] =  170.0 ; Tn++ ;
  	
  	
  	M3d_make_movement_sequence_matrix(v,vi,Tn,Ttypelist,Tvlist);
	M3d_mat_mult_points(x,y,z,v,x,y,z,np);
	
	//draw points
	for(int i=0;i<np;i++) {
		G_rgb (0, 1, 0.0); // red
   		G_point (x[i], y[i]); // hard to see
	}
	//---------------------square 1---------------------------------
	//initialize points for shape
	np=plot(square1_x,square1_y,0,2*M_PI);
	
	Tn = 0 ; 
  	Ttypelist[Tn] = SX ; Tvlist[Tn] =  150.0 ; Tn++ ;
  	Ttypelist[Tn] = SY ; Tvlist[Tn] =   70.0 ; Tn++ ;
  	Ttypelist[Tn] = TX ; Tvlist[Tn] =  500.0 ; Tn++ ;
  	Ttypelist[Tn] = TY ; Tvlist[Tn] =  460.0 ; Tn++ ;
  	
  	M3d_make_movement_sequence_matrix(v,vi,Tn,Ttypelist,Tvlist);
	M3d_mat_mult_points(x,y,z,v,x,y,z,np);
	
	//draw points
	for(int i=0;i<np;i++) {
		G_rgb (0, 0.0, 1); // red
   		G_point (x[i], y[i]); // hard to see
	}
	
	
	//---------------------square 2-----------------------------------
	//initialize points for shape
	np=plot(square2_x,square2_y,0,2*M_PI);
	
	Tn = 0 ; 
  	Ttypelist[Tn] = SX ; Tvlist[Tn] =  150.0 ; Tn++ ;
  	Ttypelist[Tn] = SY ; Tvlist[Tn] =   70.0 ; Tn++ ;
  	Ttypelist[Tn] = TX ; Tvlist[Tn] =  500.0 ; Tn++ ;
  	Ttypelist[Tn] = TY ; Tvlist[Tn] =  670.0 ; Tn++ ;
  	
  	M3d_make_movement_sequence_matrix(v,vi,Tn,Ttypelist,Tvlist);
	M3d_mat_mult_points(x,y,z,v,x,y,z,np);
	
	//draw points
	for(int i=0;i<np;i++) {
		G_rgb (0, 0.5, 1); // red
   		G_point (x[i], y[i]); // hard to see
	}
	
	//-------------------astroid----------------------
	//initialize points for shape
	np=plot(ast_x,ast_y,0,2*M_PI);
	
	Tn = 0 ; 
  	Ttypelist[Tn] = SX ; Tvlist[Tn] =   80.0 ; Tn++ ;
  	Ttypelist[Tn] = SY ; Tvlist[Tn] =   40.0 ; Tn++ ;
  	Ttypelist[Tn] = RZ ; Tvlist[Tn] =   45.0 ; Tn++ ;
  	Ttypelist[Tn] = TX ; Tvlist[Tn] =  130.0 ; Tn++ ;
 	Ttypelist[Tn] = TY ; Tvlist[Tn] =  650.0 ; Tn++ ;
 	
 	M3d_make_movement_sequence_matrix(v,vi,Tn,Ttypelist,Tvlist);
	M3d_mat_mult_points(x,y,z,v,x,y,z,np);
  
  //draw points
	for(int i=0;i<np;i++) {
		G_rgb (1, 1, 0); // red
   		G_point (x[i], y[i]); // hard to see
	}
	
	//------------------------hyperbola-----------
	//initialize points for shape
	np=plot(cosh,sinh,-1,1.5);
	
	Tn = 0 ; 
  	Ttypelist[Tn] = SX ; Tvlist[Tn] =   70.0 ; Tn++ ;
  	Ttypelist[Tn] = SY ; Tvlist[Tn] =   70.0 ; Tn++ ;
  	Ttypelist[Tn] = TX ; Tvlist[Tn] =  250.0 ; Tn++ ;
  	Ttypelist[Tn] = TY ; Tvlist[Tn] =  150.0 ; Tn++ ;
  	
  	M3d_make_movement_sequence_matrix(v,vi,Tn,Ttypelist,Tvlist);
	M3d_mat_mult_points(x,y,z,v,x,y,z,np);
  
  //draw points
	for(int i=0;i<np;i++) {
		G_rgb (0, 1, 1); // red
   		G_point (x[i], y[i]); // hard to see
	}
	
	//------------------parabola------------------
	np=plot(sum4_x,parab_y,-1,2);
	
	Tn = 0 ; 
  	Ttypelist[Tn] = SX ; Tvlist[Tn] =  150.0 ; Tn++ ;
  	Ttypelist[Tn] = SY ; Tvlist[Tn] =   50.0 ; Tn++ ;
  	Ttypelist[Tn] = RZ ; Tvlist[Tn] =   60.0 ; Tn++ ;
  	Ttypelist[Tn] = TX ; Tvlist[Tn] =  140.0 ; Tn++ ;
  	Ttypelist[Tn] = TY ; Tvlist[Tn] =  200.0 ; Tn++ ;
	
  	M3d_make_movement_sequence_matrix(v,vi,Tn,Ttypelist,Tvlist);
	M3d_mat_mult_points(x,y,z,v,x,y,z,np);
  
  	//draw points
	for(int i=0;i<np;i++) {
		G_rgb (1, 0, 1); // red
   		G_point (x[i], y[i]); // hard to see
	}
	
	//--------------lemon-----------------------------
	np=plot(lemon_x,lemon_y,0,2*M_PI);
	
	Tn = 0 ; 
  	Ttypelist[Tn] = SX ; Tvlist[Tn] =  125.0 ; Tn++ ;
  	Ttypelist[Tn] = SY ; Tvlist[Tn] =  125.0 ; Tn++ ;
  	Ttypelist[Tn] = RZ ; Tvlist[Tn] =   60.0 ; Tn++ ;
  	Ttypelist[Tn] = TX ; Tvlist[Tn] =  620.0 ; Tn++ ;
  	Ttypelist[Tn] = TY ; Tvlist[Tn] =  210.0 ; Tn++ ;
  	
  	M3d_make_movement_sequence_matrix(v,vi,Tn,Ttypelist,Tvlist);
	M3d_mat_mult_points(x,y,z,v,x,y,z,np);
  	
	//draw points
	for(int i=0;i<np;i++) {
		G_rgb (1, 1, 1);
   		G_point (x[i], y[i]); // hard to see
	}
	
	
	G_wait_key();
}

