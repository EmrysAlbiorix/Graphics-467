#define MAX_FRAMES 100
#include <unistd.h>

#include "tools.c"

// stickfigure initially designed as centered for a 400x400 window :
double x[13] = {175, 225, 225, 300, 225, 225, 250,
                200, 150, 175, 175, 100, 175};
double y[13] = {300, 300, 250, 225, 225, 200, 100,
                175, 100, 200, 225, 225, 250};
double z[13] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
// z[] values unimportant but should NOT be left uninitialized
// as nan values WILL propagate through
int n = 13;

int main(int argc, char **argv) {
  if (argc != 3) {
    printf("usage : pgm_name   window_size  microseconds(30000)\n");
    exit(0);
  }

  // the original design was for a 400x400
  // window and the object is centered on 200,200
  // so we recenter it and make it larger
  // (you get to do this ... use the
  // M3d_make_movement_sequence_matrix  function :)
  // .....
  //   double scale = winsize / 400.0;

  /***
   * render stick figure
   * wait for u microseconds
   * loop:
   * scale by 95% rotate 5 degrees
   * sleep 5 microseconds
   */
  // render stick figure

  // get info from arguments
  double winsize = atoi(argv[1]);
  int u = atoi(argv[2]);

  // init drawing
  G_init_graphics(winsize, winsize);
  G_rgb(0, 0, 0);
  G_clear();

  // consts for loop
  double trans[4][4], itrans[4][4];
  double xAvg = 0, yAvg = 0;
  double scale = .9;
  int center = winsize / 2;
  for (int i = 0; i < n; i++) {
    xAvg += x[i];
    yAvg += y[i];
    printf("%lf\t%lf\n", x[i], y[i]);
  }
  xAvg /= n;
  yAvg /= n;
  printf("Average:\t%lf\t%lf\n", xAvg, yAvg);

  int initialMoveTypes[] = {TX, TY, SX, SY, TX, TY};
  double initialMoveParams[] = {-xAvg,           -yAvg,  winsize / 400.0,
                                winsize / 400.0, center, center};
  M3d_make_movement_sequence_matrix(trans, itrans, 6, initialMoveTypes,
                                    initialMoveParams);
  M3d_mat_mult_points(x, y, z, trans, x, y, z, n);

  // draw initial view
  G_rgb(0, 1, 0);
  int res = G_fill_polygon(x, y, n);
  printf("res: %d\n", res);

  G_wait_key();
  // wait for u microseconds
  printf("beforeSleep\n");
  usleep(u);
  printf("afterSleep\n");

  int moveTypes[] = {TX, TY, RZ, SX, SY, TX, TY};
  double moveParams[] = {-center, -center, 4, scale, scale, center, center};
  M3d_make_movement_sequence_matrix(trans, itrans, 7, moveTypes, moveParams);
  for (int frame = 0; frame < MAX_FRAMES; frame++) {
    // debug prints
    for (int i = 0; i < n; i++) {
      printf("%lf\t%lf\n", x[i], y[i]);
    }

    // clear canvas
    G_rgb(0, 0, 0);
    G_clear();

    // draw new
    G_rgb(0, 1, 0);
    G_fill_polygon(x, y, n);

    // do math
    // clang_format off
    M3d_mat_mult_points(x, y, z, trans, x, y, z, n);
    // clang_format on

    // sleep/wait key
    G_wait_key();
    usleep(30000);
    frame++;
  }

  G_rgb(0, 0, 0);
  G_clear();
  //   G_wait_mouse();
}
