/**
# Axisymmetric pipe flow: embedded wall vs masked (bid) wall give different physics

This is a minimal demonstration of what appears to be a serious problem
with the **embed + axi** combination: the advective (inertia) term of the
momentum equation effectively acquires a spurious factor $y$ (the radius),
so that inertia vanishes toward the axis. Diffusion, pressure and mass
conservation are unaffected, so the fully-developed (Poiseuille) state is
exactly right -- but every inertia-controlled quantity is wrong. In
particular the laminar entrance length comes out ~4.4 times too short.

The same pipe-entrance problem is solved here twice, changing *only* the
way the wall at $y = 1$ is represented:

* default: a masked region with a `bid` boundary (staircase wall, exactly
  face-aligned since the domain size 8 is divisible by the radius 1);
* with `-DEMBED_WALL=1`: an embedded boundary, using the cm/fm metric
  updates required to combine *embed.h* with *axi.h* (the usual recipe,
  cf. several sandbox examples -- if there is a correct recipe we would
  love to know it).

Uniform inflow $u = 1$ enters a pipe of radius 1 with viscosity $1/Re$
($Re \equiv U a/\nu$, radius-based). The classical laminar entrance length
(Shah & London; Durst et al. 2005) is
$$ L_{99} = 0.056 \times Re_D \times D = 0.224 \times Re \times a. $$
With $Re = 25$ (default) this is $5.6\,a$, comfortably inside the
8-radii-long domain, so both runs cross the 99% line at a measurable
$z_{99}$: the bid run crosses where the classical correlation says it
should ($z_{99} \approx 5.6$), the embed run roughly $6\times$ earlier
($z_{99} \approx 0.9$).

The second, sharper diagnostic is a pointwise momentum balance evaluated
from the solution by finite differences,
$$ R = u u_z + v u_r + p_z - \frac{1}{Re}\Bigl[\frac1r (r u_r)_r + u_{zz}\Bigr], $$
printed at probe points in the entrance region. For the bid run $R \approx 0$
(discretisation noise). For the embed run $R \approx (1-y)\,(u u_z + v u_r)$:
the code behaves as if the advection term were multiplied by $y$.

Compile and run both:

~~~bash
qcc -disable-dimensions -O2 -Wall pipe-axi-embed-vs-bid.c -o pipe-bid -lm
qcc -disable-dimensions -O2 -Wall -DEMBED_WALL=1 pipe-axi-embed-vs-bid.c -o pipe-embed -lm
qcc -disable-dimensions -O2 -Wall -DEMBED_WALL=1 -DPATCH_ADVECTION=1 \
    pipe-axi-embed-vs-bid.c -o pipe-embed-fixed -lm
./pipe-bid 25 8 ; ./pipe-embed 25 8 ; ./pipe-embed-fixed 25 8   # args: Re LEVEL
~~~

(`-disable-dimensions` is required: the constants in this demo are not
dimension-annotated). Each run writes the final-time centerline profile
to `uc-<mode>-<Re>-<LEVEL>` and prints the residual table and a summary
to standard output. */

#if EMBED_WALL
# include "embed.h"
# if PATCH_ADVECTION
#  define MODESTR "embed-fixed"
# else
#  define MODESTR "embed"
# endif
#else
# define MODESTR "bid"
#endif
#include "axi.h"

/**
## The one-line fix (compile with -DEMBED_WALL=1 -DPATCH_ADVECTION=1)

The culprit is `update_tracer()` in *embed.h* (the "small cell" advection
update, to which *bcg.h* delegates the flux divergence when EMBED is
defined). Its full-cell branch divides the flux difference by $\Delta$
only, whereas the generic (non-embed) update in *bcg.h* divides by
$\Delta\,c_m$. A comment in *embed.h* states the limitation explicitly:
"the distinction should be made between $c_m$, the cell fraction metric,
and $c_s$, the embedded fraction. This is not done now so that embedded
boundaries cannot be combined with a metric yet."

With *axi.h*, $c_m = y\,c_s$: in full cells the update misses the factor
$1/y$ while the fluxes (built from $u_f$ which carries $f_m = y f_s$)
already carry a factor $y$ -- the advection term is effectively
multiplied by $y$, which is exactly the residual law measured below.

The replacement below is a verbatim copy of `update_tracer()` except
that the full-cell branch divides by `Delta*cm[]`. Since embed's own
metric event makes $c_m \equiv c_s$ when no other metric is present
(so $c_m = 1$ in full cells), the change is a no-op for pure embed and
repairs the embed+metric combination. */

#if EMBED_WALL && PATCH_ADVECTION
void update_tracer_metric (scalar f, face vector uf, face vector flux, double dt)
{
  scalar e[];
  foreach() {
    if (cs[] <= 0.)
      e[] = 0.;
    else if (cs[] >= 1.) {
      foreach_dimension()
	f[] += dt*(flux.x[] - flux.x[1])/(Delta*cm[]);   // was: /Delta
      e[] = 0.;
    }
    else {
      double umax = 0.;
      for (int i = 0; i <= 1; i++)
	foreach_dimension()
	  if (fabs(uf.x[i]) > umax)
	    umax = fabs(uf.x[i]);
      double dtmax = Delta*cm[]/(umax + SEPS);
      double F = 0.;
      foreach_dimension()
	F += flux.x[] - flux.x[1];
      F /= Delta*cm[];
      if (dt <= dtmax) {
	f[] += dt*F;
	e[] = 0.;
      }
      else {
	f[] += dtmax*F;
	double scs = 0.;
	foreach_neighbor(1)
	  scs += sq(cm[]);
	e[] = (dt - dtmax)*F*cm[]/scs;
      }
    }
  }
  foreach() {
    double se = 0.;
    foreach_neighbor(1)
      se += e[];
    f[] += cs[]*se;
  }
}
#define update_tracer update_tracer_metric
#endif // EMBED_WALL && PATCH_ADVECTION

#include "navier-stokes/centered.h"

double Re = 25.;
int LEVEL = 8;

const double LENGTH = 8.;    // divisible by the radius: face-aligned wall
face vector muv[];

#if !EMBED_WALL
bid pipewall;
#endif

int main (int argc, char * argv[])
{
  if (argc > 1) Re = atof (argv[1]);
  if (argc > 2) LEVEL = atoi (argv[2]);
  printf ("# mode %s Re %g LEVEL %d\n", MODESTR, Re, LEVEL);

  size (LENGTH);
  init_grid (1 << LEVEL);
  mu = muv;
  TOLERANCE = 1e-7;
  run();
}

event properties (i++) {
  foreach_face()
    muv.x[] = fm.x[]/Re;
}

/**
Uniform inflow on the left, outflow on the right, no-slip on the wall. */

u.n[left]  = dirichlet (y < 1. ? 1. : 0.);
p[left]    = neumann (0.);
pf[left]   = neumann (0.);

u.n[right] = neumann (0.);
p[right]   = dirichlet (0.);
pf[right]  = dirichlet (0.);

#if EMBED_WALL
u.x[embed] = dirichlet (0.);
u.y[embed] = dirichlet (0.);
#else
u.n[pipewall] = dirichlet (0.);
u.t[pipewall] = dirichlet (0.);
#endif

event init (t = 0)
{
#if EMBED_WALL
  double eps = L0/(1 << LEVEL)/1000.;
  for (scalar s in {u, p, pf})
    s.third = true;
  solid (cs, fs, 1. - y - eps);
  fractions_cleanup (cs, fs);

  /* the metric updates needed to combine embed with axi */
#if AXI
  cm_update (cm, cs, fs);
  fm_update (fm, cs, fs);
#if TREE
  cm.refine = cm.prolongation = refine_cm_axi;
  cs.refine = cs.prolongation = fraction_refine;
  fm.x.refine = refine_face_x_axi;
  fm.y.refine = refine_face_y_axi;
  metric_embed_factor = axi_factor;
#endif
  restriction ({cm, fm, cs, fs});
#endif
  foreach()
    u.x[] = cs[] ? (y < 1.) : 0.;
#else // bid wall
  mask (y > 1. ? pipewall : none);
  foreach()
    u.x[] = y < 1.;
#endif
  boundary (all);
}

/**
## Diagnostics at the final time $t = Re$

The flow is steady well before $t = Re$ (one radial diffusion time). */

event diagnostics (t = Re)
{
  double h = L0/(1 << LEVEL);

  /* centerline profile */
  char name[80];
  sprintf (name, "uc-%s-%g-%d", MODESTR, Re, LEVEL);
  FILE * fp = fopen (name, "w");
  double z99 = -1.;
  for (double z = h; z < L0; z += h) {
    double uc = interpolate (u.x, z, h/4.);
    fprintf (fp, "%g %g\n", z, uc);
    if (z99 < 0. && uc >= 0.99*2.)
      z99 = z;
  }
  fclose (fp);

  /* pointwise momentum balance in the entrance region */
  printf ("#     z     r      advec      resid   resid/((1-y)*advec)\n");
  double zp[6] = {1., 1., 2., 2., 3., 3.}, rp[6] = {0.3, 0.5, 0.3, 0.5, 0.5, 0.7};
  for (int k = 0; k < 6; k++) {
    double z0 = zp[k], r0 = rp[k];
    double uc = interpolate (u.x, z0, r0);
    double uzp = interpolate (u.x, z0 + h, r0), uzm = interpolate (u.x, z0 - h, r0);
    double urp = interpolate (u.x, z0, r0 + h), urm = interpolate (u.x, z0, r0 - h);
    double v   = interpolate (u.y, z0, r0);
    double dpz = (interpolate (p, z0 + h, r0) - interpolate (p, z0 - h, r0))/(2.*h);
    double uz = (uzp - uzm)/(2.*h), ur = (urp - urm)/(2.*h);
    double lap = (urp - 2.*uc + urm)/sq(h) + ur/r0 + (uzp - 2.*uc + uzm)/sq(h);
    double adv = uc*uz + v*ur;
    double res = adv + dpz - lap/Re;
    printf ("# %5g %5g %10.5f %10.5f %10.3f\n", z0, r0, adv, res,
	    res/((1. - r0)*adv));
  }
  printf ("# summary: mode %s Re %g LEVEL %d z99 %g u_c(exit) %g "
	  "(theory: L99 = 0.224*Re = %g)\n",
	  MODESTR, Re, LEVEL, z99, interpolate (u.x, L0 - h, h/4.), 0.224*Re);
}

/**
## Results

With the defaults ($Re = 25$, level 8) the two executables give, at $t = Re$:

* **bid**: $z_{99}$ close to the classical $L_{99} = 5.6$ (and to an
  independent Galerkin/parabolized solution of the same problem), with
  momentum residuals at the percent level of the advection term;
* **embed**: $z_{99} \approx 0.9$, roughly $6\times$ too short, and
  momentum residuals equal to $(1-y)$ times the advection term -- the
  last column reads $\approx 0.99$ at every probe -- i.e.\ the effective
  advection is $y\,(u\cdot\nabla)u$.

The same comparison at $Re = 50$ (run with `./pipe-bid 50 8` etc., where
$L_{99} = 11.2$ exceeds the domain) gives bid $u_c(8) = 1.947$ against
$1.945$ from the Galerkin solution, embed $z_{99} = 1.66$, and residual
ratios $0.984$--$0.991$.

With the one-line fix (`pipe-embed-fixed`), the embedded-wall run
recovers the bid/classical result: measured $z_{99} = 5.56$ (bid: 5.59,
theory: 5.6), the centerline profile agrees with the bid run to
$\lesssim 0.005$ everywhere, and the momentum residuals collapse from
$0.99\times(1-y)\,\mathrm{advec}$ to the percent level -- while the
patch remains bit-identical to the unpatched code for pure embed (no
metric) computations, since embed's metric event sets $c_m \equiv c_s$.

Because the spurious factor multiplies only the *inertia*, the developed
Poiseuille state (profile shape, $dp/dz = -8/Re_{\rm mean}$, flux) is
exactly right, which makes the problem easy to miss: entrance lengths,
and presumably all unsteady/inertial axisymmetric embed computations,
are affected while every steady diagnostic looks fine.

~~~gnuplot Centerline velocity: embed vs bid wall
set xlabel 'z / a'
set ylabel 'u(r=0)'
set key bottom right
set arrow from 5.6, 1.0 to 5.6, 1.98 nohead dt 2 lc rgb 'blue'
set label 'L_{99} = 0.224*Re' at 5.7, 1.3 tc rgb 'blue'
plot 'uc-bid-25-8' w l t 'bid wall', \
     'uc-embed-25-8' w l t 'embedded wall', \
     'uc-embed-fixed-25-8' w l t 'embedded wall (fixed)', \
     1.98 t '99% developed' lc rgb 'gray'
~~~ */
