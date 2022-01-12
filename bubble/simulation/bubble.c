
#include "axi.h"

#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "curvature.h"
#include "view.h"

//#include "navier-stokes/swirl.h"


#define mu(f)  (1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2))

#define RHOR 1000
#define MUR 100

#define Ga 30
#define Bo 100

#define radius 1


int LEVEL = 9; //12


u.n[top] = dirichlet(0); 
u.t[top] = dirichlet(0);  

u.n[right] = neumann(0); 
u.n[left] = neumann(0); 

//w[left] = dirichlet(y);
//w[right] = dirichlet(0);
//w[top] = dirichlet(0);

int main() {
  size (8);
  origin(-1.6, 0);
  init_grid (1 << 8);

  rho1 = 1.;
  rho2 = 1./RHOR;
  mu1 = 1./Ga; 
  mu2 = 1./(MUR*Ga);
  f.sigma = 1./Bo;

  //G.x = -1;

  run();
}

//bid shape;

event init (t = 0) {
      //mask (fabs(y) > 8 ? shape : none);
      refine (sq(x) + sq(y) - sq(radius*1.20) < 0 && level < 10);
      fraction (f, sq(x) + sq(y) - sq(radius));
}


event acceleration (i++) {
  face vector av = a;
  foreach_face(x)
    av.x[] = -1.;
}


event adapt (i++) {
  double uemax = 5e-3;
  adapt_wavelet ({f,u}, (double[]){0.01,uemax,uemax,uemax}, 9, 5);
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

event movies (t += 4./400.) {
  view (fov= 20, ty = -0.3, quat = {0, 0, -0.707, 0.707}, bg = {0,0,0});
  scalar vort[];
  vorticity (u, vort);

  draw_vof ("f", lw = 3);
  squares ("u.x", linear = true);
  mirror ({0,1}) {
    draw_vof ("f", lw = 3);
    squares ("u.y", linear = true);
   // vectors("u", scale = 0.006);
  }
  save ("movie.mp4");
}

event images (t = 4) {
  view (fov= 18, ty = -0.3, quat = {0, 0, -0.707, 0.707}, bg = {0,0,0});
  draw_vof ("f", lw = 3);
  squares ("u.x", linear = true);
  mirror ({0,1}) {
    draw_vof ("f", lw = 3);
    squares ("u.y", linear = true);
  }
  save ("photo1.png");
}


event logo (i++) {
  fprintf (stdout, "%g \n", t);
  }


event end (t = 5) {
  printf ("i = %d t = %g\n", i, t);
}