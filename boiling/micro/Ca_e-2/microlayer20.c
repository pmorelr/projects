#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "view.h"
#include "contact.h"

// Dimensionless parameters :
#define RHOR 100 //[5.9; 1650]
#define MUR 10 //[3.2; 23.3]
#define Ca 0.1 //[0.001; 0.1]
#define R0 1
//#define tc 1e-2

char str[99];
int LEVEL = 7; //12

vector h[];

p[right]  = neumann(0);
p[left]  = neumann(0);
p[top]  = neumann(0);

h.t[left] = contact_angle (20*pi/180.);

int main() {
  size (2);
  //origin(-8, 0);
  init_grid (1 << LEVEL);

  f.height = h;

  rho1 = 1.;
  rho2 = 1./RHOR;
  mu1 = 1.; 
  mu2 = 1./(MUR);
  f.sigma = 1./Ca;

  run();
}


event init (t = 0) {
      sprintf (str, "Level = %i", LEVEL);
      refine (sq(x) + sq(y) - sq(R0*1.20) < 0 && level < LEVEL);
      fraction (f, sq(x) + sq(y) - sq(R0));
}

event boiling (t++) {
  foreach(){
    Sv[] = 3*t*(1 - f[])/(1e-6*(1 + t)); 
  }
}

event interface_position (i+= 10) {

  static FILE * fp2 = fopen ("xm", "w");

  int counter = 0;
  double xa [5];
  
  foreach() {
    if (f[] > 0. && f[] < 1. && y - Delta/2. <= 1e-6) {// - 1e-6
    face vector s; 
    coord segment[2];
    s.x.i = -1;
    coord n = facet_normal (point, f, s);
    double alpha = plane_alpha (f[], n);
    facets (n, alpha, segment);

    xa[counter++] = x + segment[0].x*Delta;
    }
  }
  if(xa[counter-1] < 1e6 && (xa[counter-1] - xa[0]) > 1e-6){
   fprintf (fp2, "%i %g %g %g\n", counter, t, xa[0], xa[counter-1]);
  }
}


event adapt (i++) {
  //double uemax = 5e-3;
  adapt_wavelet ({p,f}, (double[]){0.005, 0.01}, LEVEL, 5); //adapt_wavelet ({f,u}, (double[]){0.01,uemax,uemax}, LEVEL, 5);
}


event movies (t += 1.5/500.) {
  view (fov= 20, ty = -0.5, quat = {0, 0, -0.707, 0.707}, bg = {0,0,0});

  draw_vof ("f", lw = 4);
  squares ("p", linear = true);
  box();
  mirror ({0,1}) {
    draw_vof ("f", lw = 4);
    cells();
   // vectors("u", scale = 0.006);
  }
  draw_string (str, pos = 1, lw  = 2, size = 30);
  save ("movie.mp4");
}

event logo (i++) {
  fprintf (stdout, "%g \n", t);
  }


event end (t = 0.9) {
  output_facets (f, stdout);
  printf ("i = %d t = %g\n", i, t);
}