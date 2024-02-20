#include "Tools/FPToolkit.c"
#include "Tools/M3d_matrix_tools.c"

// stickfigure initially designed as centered for a 400x400 window :
double x[13] = {175,225,225,300,225,225,250,200,150,175,175,100,175} ;
double y[13] = {300,300,250,225,225,200,100,175,100,200,225,225,250} ;
double z[13] = {0,0,0,0,0,0,0,0,0,0,0,0,0} ;
       // z[] values unimportant but should NOT be left uninitialized
       // as nan values WILL propagate through
int n = 13 ;



int main(int argc, char **argv) {

 if (argc != 3) {
    printf("usage : pgm_name   window_size  microseconds(30000)\n") ;
    exit(0) ;
 }
 
 // User input parameters
 double winsize = atoi(argv[1]) ;
 int u = atoi(argv[2]) ; 

 // Creates the window
 G_init_graphics(winsize,winsize) ;
 G_rgb(0, 0, 0);
 G_clear();
 
 // Pauses for given time u
 usleep(u);

 // the original design was for a 400x400
 // window and the object is centered on 200,200
 // so we recenter it and make it larger
 // (you get to do this ... use the
 // M3d_make_movement_sequence_matrix  function :)

 // .....
 
 double halfsize = winsize/2 ;
 double scale = halfsize/200 ;
 double matrix[4][4], iMatrix[4][4] ;
 int tName[100] ; double tMove[100] ;
 
 int frame = 0 ;
 tName[frame] = TX ; tMove[frame] = -200 ;	frame++ ;
 tName[frame] = TY ; tMove[frame] = -200 ;	frame++ ;
 tName[frame] = SX ; tMove[frame] = scale ; 	frame++ ;
 tName[frame] = SY ; tMove[frame] = scale ;	frame++ ;
 tName[frame] = TX ; tMove[frame] = halfsize ;	frame++ ;
 tName[frame] = TY ; tMove[frame] = halfsize ;	frame++ ;
 
 M3d_make_movement_sequence_matrix(matrix, iMatrix, frame, tName, tMove) ;
 M3d_mat_mult_points(x, y, z, matrix, x, y, z, n) ;
 
 // Creates the figure
 G_rgb(0, 1, 0);
 G_fill_polygon(x, y, n);
 G_wait_key() ;

 
 // now make the movie the rotates and shrinks about the center :
 // .....
  
 scale = 0.95 ;
 double degree = 4 ;
 int num = 100 ;
 
 frame = 0 ;
 tName[frame] = TX ; tMove[frame] = -halfsize ;	frame++ ;
 tName[frame] = TY ; tMove[frame] = -halfsize ;	frame++ ;
 tName[frame] = SX ; tMove[frame] = scale ; 	frame++ ;
 tName[frame] = SY ; tMove[frame] = scale ;	frame++ ;
 tName[frame] = RZ ; tMove[frame] = degree ;	frame++ ;
 tName[frame] = TX ; tMove[frame] = halfsize ;	frame++ ;
 tName[frame] = TY ; tMove[frame] = halfsize ;	frame++ ;
 
 M3d_make_movement_sequence_matrix(matrix, iMatrix, frame, tName, tMove) ;

 // Loops num times
 for (int i = 0; i < num; i++) {
    M3d_mat_mult_points(x, y, z, matrix, x, y, z, n) ;
    G_rgb(0, 0, 0);
    G_clear();
    G_rgb(0, 1, 0);
    G_fill_polygon(x, y, n);
    G_display_image() ;
    usleep(u) ;
 }
 
 // Blood splatter
 G_rgb(0.8, 0, 0) ;
 G_fill_circle(halfsize, halfsize, winsize/32) ;
 G_wait_key() ;
}

