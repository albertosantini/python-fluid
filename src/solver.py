""" Real-Time Fluid Dynamics for Games by Jos Stam (2003).
Parts of author's work are also protected
under U. S. patent #6,266,071 B1 [Patent].
"""
def set_bnd(N, b, x):
  """We assume that the fluid is contained in a
  box with solid walls: no flow should exit the walls. This simply means that the horizontal
  component of the velocity should be zero on the vertical walls, while the vertical component of
  the velocity should be zero on the horizontal walls. For the density and other fields considered
  in the code we simply assume continuity. The following code implements these conditions.
  """
  for i in range(1, N+1):
    if b == 1: x[0,i] = -x[1,i]
    else: x[0,i] = x[1,i]
    if b == 1: x[N+1,i] = -x[N,i]
    else: x[N+1,i] = x[N,i]
    if b == 2: x[i,0] = -x[i,1]
    else: x[i,0] = x[i,1]
    if b == 2: x[i,N+1] = -x[i,N]
    else: x[i,N+1] = x[i,N]

  x[0,0] = 0.5*(x[1,0]+x[0,1])
  x[0,N+1] = 0.5*(x[1,N+1]+x[0,N])
  x[N+1,0] = 0.5*(x[N,0]+x[N+1,1])
  x[N+1,N+1] = 0.5*(x[N,N+1]+x[N+1,N])

def lin_solve(N, b, x, x0, a, c):
  for k in range(0, 20):
    x[1:N+1,1:N+1] = (x0[1:N+1,1:N+1] + a*(x[0:N,1:N+1]+x[2:N+2,1:N+1]+x[1:N+1,0:N]+x[1:N+1,2:N+2]))/c
    set_bnd(N, b, x)

# Addition of forces: the density increases due to sources
def add_source(N, x, s, dt):
  size = (N+2)
  x[0:size,0:size] += dt*s[0:size,0:size]

# Diffusion: the density diffuses at a certain rate
def diffuse (N, b, x, x0, diff, dt):
  """ The basic idea behind our method is to find the densities which when diffused backward
  in time yield the densities we started with.
  The simplest iterative solver which works well in practice is Gauss-Seidel relaxation.
  """
  a = dt*diff*N*N
  lin_solve(N, b, x, x0, a, 1+4*a)

# Advection: the density follows the velocity field
def advect (N, b, d, d0, u, v, dt):
  """ The basic idea behind the advection step. Instead of moving the cell centers forward in
  time through the velocity field, we look for the particles which end up exactly at
  the cell centers by tracing backwards in time from the cell centers.
  """
  dt0 = dt*N;
  for i in range(1, N+1):
    for j in range(1, N+1):
      x = i-dt0*u[i,j]; y = j-dt0*v[i,j]
      if x<0.5: x=0.5
      if x>N+0.5: x=N+0.5
      i0 = int(x); i1 = i0+1
      if y<0.5: y=0.5
      if y>N+0.5: y=N+0.5
      j0 = int(y); j1 = j0+1
      s1 = x-i0; s0 = 1-s1; t1 = y-j0; t0 = 1-t1
      d[i,j] = s0*(t0*d0[i0,j0]+t1*d0[i0,j1])+s1*(t0*d0[i1,j0]+t1*d0[i1,j1])
  set_bnd (N, b, d)

def project(N, u, v, p, div):
  h = 1.0/N
  div[1:N+1,1:N+1] = -0.5*h*(u[2:N+2,1:N+1]-u[0:N,1:N+1]+v[1:N+1,2:N+2]-v[1:N+1,0:N])
  p[1:N+1,1:N+1] = 0
  set_bnd (N, 0, div)
  set_bnd (N, 0, p)
  lin_solve (N, 0, p, div, 1, 4)
  u[1:N+1,1:N+1] -= 0.5*(p[2:N+2,1:N+1]-p[0:N,1:N+1])/h # ??? not in the paper /h
  v[1:N+1,1:N+1] -= 0.5*(p[1:N+1,2:N+2]-p[1:N+1,0:N])/h # ??? not in the paper /h
  set_bnd (N, 1, u)
  set_bnd (N, 2, v)

# Evolving density: advection, diffusion, addition of sources
def dens_step (N, x, x0, u, v, diff, dt):
  add_source(N, x, x0, dt)
  x0, x = x, x0 # swap
  diffuse(N, 0, x, x0, diff, dt)
  x0, x = x, x0 # swap
  advect(N, 0, x, x0, u, v, dt)

# Evolving velocity: self-advection, viscous diffusion, addition of forces
def vel_step (N, u, v, u0, v0, visc, dt):
  add_source(N, u, u0, dt)
  add_source(N, v, v0, dt);
  u0, u = u, u0 # swap
  diffuse(N, 1, u, u0, visc, dt)
  v0, v = v, v0 # swap
  diffuse(N, 2, v, v0, visc, dt)
  project(N, u, v, u0, v0)
  u0, u = u, u0 # swap
  v0, v = v, v0 # swap
  advect(N, 1, u, u0, u0, v0, dt)
  advect(N, 2, v, v0, u0, v0, dt)
  project(N, u, v, u0, v0)

