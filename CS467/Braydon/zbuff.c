#include "FPToolkit.c"
#include "M3d_matrix_tools.c"


double x[10000],y[10000],z[10000];
double xv[10000],yv[10000],zv[10000];
int np;

int plot (double (*f)(double u, double v), double (*g)(double u, double v), double (*l)(double u, double v), double ulo, double uhi,double vlo, double vhi) {
	double u,v;
	int i=0;
	
	for (u = ulo; u <= uhi ; u+=0.01) {
    	for(v=vlo; v<= vhi; v+=0.01) {
    		G_rgb (1.0, 0.0, 0.0) ; // red
   		x[i]=f(u,v);
   		y[i]=g(u,v);
   		z[i]=l(u,v);
   		i++;
   		v+=0.1;
    	}		
  	}
  	return i;
}

void initGraphics() {
	G_init_graphics(800,800);
	G_rgb(0.1,0.1,0.1);
	G_clear();
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

void draw_object() {
	double xbb[10000],ybb[10000];
	double halfangle;
	double H;
	int i;
	
	halfangle=40;
	H=tan(halfangle*(M_PI/180));
	
	for(i=0;i<np;i++) {
		xbb[i]=((400/H)*(x[i]/z[i])+400);
		ybb[i]=((400/H)*(y[i]/z[i])+400);
		G_rgb(0,1,0);
		G_point(xbb[i],ybb[i]);
	}
}

int main() {
	double u;
	double Tvlist[100];
	double v[4][4],vi[4][4];
	int Ttypelist[100],Tn;
	double normal[2],mag;
	
	//initialize graphics pannel
	G_init_graphics(800,800);
	G_rgb(0.1,0.1,0.1);
	G_clear();
	
	//initialize points for sphere
	np=plot(sphere_x,sphere_y,sphere_z,-1*M_PI,M_PI,-1*M_PI/2,M_PI/2);
	
	//translate and scale sphere
	Tn = 0; 
	Ttypelist[Tn] = SX ; Tvlist[Tn] =   50.0 ; Tn++;
	Ttypelist[Tn] = SY ; Tvlist[Tn] =   50.0 ; Tn++;
	Ttypelist[Tn] = SZ ; Tvlist[Tn] =   50.0 ; Tn++;
  	Ttypelist[Tn] = TZ ; Tvlist[Tn] =   10.0 ; Tn++;
  	
	M3d_make_movement_sequence_matrix(v,vi,Tn,Ttypelist,Tvlist);
	M3d_mat_mult_points(x,y,z,v,x,y,z,np);
	
	
	G_rgb(1,0,0);
	draw_object();
	
	//wait for user to click key and end program.
	G_wait_key();
}

