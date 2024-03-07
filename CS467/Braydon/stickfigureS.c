#include "FPToolkit.c"
#include "M3d_matrix_tools.c"

// stickfigure initially designed as centered for a 400x400 window :
double x[13] = {175,225,225,300,225,225,250,200,150,175,175,100,175} ; 
double y[13] = {300,300,250,225,225,200,100,175,100,200,225,225,250} ;
double z[13] = {0,0,0,0,0,0,0,0,0,0,0,0,0} ;
       // z[] values unimportant but should NOT be left uninitialized
       // as nan values WILL propagate through
int np = 13 ;



int main(int argc, char **argv) 
{

 if (argc != 3) {
    printf("usage : pgm_name   window_size  microseconds(30000)\n") ;
    exit(0) ;
 }
  
 	double winsize = atoi(argv[1]) ;
 	int u = atoi(argv[2]) ; 

	double v[4][4], vi[4][4], p[4][4], m[4][4];
  	int n;
  	int mtype[100] ;
  	double mparam[100] ;

 	G_init_graphics(winsize,winsize);
	G_rgb(0.2,0.2,0.2);
 	G_clear();

 // the original design was for a 400x400
 // window and the object is centered on 200,200
 // so we recenter it and make it larger
 // (you get to do this ... use the
 // M3d_make_movement_sequence_matrix  function :)
	
	
	//x width 200
	//y width 200
	
	//initial adjustments for middle
	n = 0 ;
  	mtype[n] = TX ;  mparam[n] =  -200 ; n++ ;
 	mtype[n] = TY ;  mparam[n] =  -200 ; n++ ;
 	mtype[n] = SX ;  mparam[n] =  (winsize/200)/2; n++ ;
 	mtype[n] = SY ;  mparam[n] =  (winsize/200)/2; n++ ;
 	mtype[n] = TX ;  mparam[n] =  winsize/2 ; n++ ;
 	mtype[n] = TY ;  mparam[n] =  winsize/2 ; n++ ;
	
	M3d_make_movement_sequence_matrix(v,vi,n,mtype,mparam);
	M3d_mat_mult_points(x,y,z,v,x,y,z,np);
	

	
	//make rotation and shrinking
	n = 0;
	mtype[n] = TX ;  mparam[n] =  -winsize/2 ; n++ ;
 	mtype[n] = TY ;  mparam[n] =  -winsize/2 ; n++ ;
  	mtype[n] = RZ ;  mparam[n] =  5; n++ ;
  	mtype[n] = SX ;  mparam[n] = .95; n++ ;
 	mtype[n] = SY ;  mparam[n] = .95; n++ ;
 	mtype[n] = TX ;  mparam[n] =  winsize/2 ; n++ ;
 	mtype[n] = TY ;  mparam[n] =  winsize/2 ; n++ ;
 
	M3d_make_movement_sequence_matrix(v,vi,n,mtype,mparam);
	
	//draw shape
	G_rgb(0,1,0);
	G_polygon(x,y,np);
	
	int key,i;
	
	//wait so user can initialize movement
	G_wait_key();
	//loop to spin and shrink
	i=0;
	while (i<70) {
		M3d_mat_mult_points(x,y,z,v,x,y,z,np);
 		G_rgb(0.2,0.2,0.2);
 		G_clear();
 		G_rgb(0,1,0);
		G_polygon(x,y,np);
		G_display_image();
		usleep(u);
		i++;
	}

	//make blood splatter
	G_rgb(1,0,0);
	G_fill_circle(winsize/2,winsize/2,100);
	//wait again so user can exit program
	G_wait_key();
	
}

