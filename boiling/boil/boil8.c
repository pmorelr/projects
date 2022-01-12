#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "view.h"

#define GRAVITY 0

// Dimensionless parameters :
#define RHOR 100 //[5.9; 1650]
#define MUR 10 //[3.2; 23.3]
#define Ca 0.1 //[0.001; 0.1]
#define R0 1

char str[99];
int LEVEL = 8; //12

p[right]  = neumann(0);
p[left]  = neumann(0);
p[top]  = neumann(0);
//u.t[right]  = neumann(0);
//u.t[left]  = neumann(0);
//u.t[top]  = neumann(0);

//uf.n[bottom] = 0;

/**
u.n[top] = neumann(0); 
u.t[top] = neumann(0);  

//u.n[bottom] = dirichlet(0); 
//u.t[bottom] = dirichlet(0);  

u.n[right] = neumann(0); 
u.n[left] = neumann(0); 
*/

int main() {
  size (4);
  origin(-2, 0);
  init_grid (1 << LEVEL);

  rho1 = 1.;
  rho2 = 1./RHOR;
  mu1 = 1.; 
  mu2 = 1./(MUR);
  f.sigma = 1./Ca;

  // Variation of the source term
  //Sv = 3*f/(R0 + t);

  run();
}

event init (t = 0) {
      sprintf (str, "Level = %i", LEVEL);
      refine (sq(x) + sq(y) - sq(R0*1.20) < 0 && level < LEVEL);
      fraction (f, sq(x) + sq(y) - sq(R0));
}


event boiling (t++) {
  foreach(){
    Sv[] = 10*(1 - f[]); //3*f[]/(R0 + t);
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

event u_errors (t+= 0.1) {

  static FILE * fp3 = fopen ("ue", "w");

  norm nx = normf(u.x);
  norm ny = normf(u.y);
  fprintf (fp3, "%g %g %g %g %g\n", t, nx.rms, nx.max, ny.rms, ny.max);
}


event adapt (i++) {
  //double uemax = 5e-3;
  adapt_wavelet ({p,f,Sv}, (double[]){0.01,0.01,0.01}, LEVEL, 5);
}

event movies (t += 3./500.) {
  view (fov= 18, quat = {0, 0, -0.707, 0.707}, bg = {0,0,0});

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


event end (t = 3) {
  printf ("i = %d t = %g\n", i, t);
}