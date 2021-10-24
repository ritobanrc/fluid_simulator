# Next Steps

* [X] Reset Button  
* JSON Serialization/Deserialization of settings
* Partio/VTK input/output
* Fix APIC Stability issue (invertibility?)
* Marching Cubes
* More Initial Conditions
    - Report suggested rest density for initial conditions 

* Refactor Constiutive Models into a separate file
* Better Camera Controller

## Longer Term
* Generic over 2D vs 3d?
* Other Degree Kernels for MPM
     https://www.desmos.com/calculator/gwhf6t3lqs


## PhysBAM Stuff
What classes in PhysBAM would you actually need to port to get most of it?

### Not necessary
* Vectors and Matricies (provided by Nalgebra)
* Arrays not necessary (Rust's `Vec` is probably good enough)
* Math_Tools not necessary (we may or may not want to re-implement RANGE)
* LOG -- just use rust's `log` library, maybe `slog` or `tracing` (including scope timers)
* Random Numbers -- just use `rand`

### Necessary
* GRID -- this is a really handy class
* Iterators -- RANGE_ITERATOR, NODE_ITERATOR, CELL_ITERATOR, FACE_ITERATOR (w/ FACE_INDEX), these can probably be written pretty cleanly with Rust's iterator (better than separate Next and Valid functions)
* Arrays Nd -- we can probably implement these ourselves better than ndarray
* Particles -- we can definitely come up with a better way to do this than PhysBAM's dynamic particle arrays
* Poisson Disk sampling
* Grid_PDE
    - Advection Trait w/ Semi lagrangian (maybe generic over a particle tracing API?), WENO, ENO, etc.
    - There has to be a more consistent way to handle boundaries
    - Craig's new Interpolation API is _very_ nice, porting that to rust would be cool (and not being tied down by backwards compatibility would be nice as well)
* Meshes -- this is necessary
* Output -- somehow use serde for this, if binary is necessary (which it probably is), serde can handle it -- maybe partio for particles
* Level Sets -- probably fine
    - what about analytic implicit objects? i'm not sure...
* Krylov Solvers
    - Don't know this area well enough, but I'm sure PhysBAM's API can be cleaned up -- LAPLACE, POISSON, PROJECTION, KRYLOV_SYSTEM_BASE
* Constitutive Models -- the rest of it probably isn't strictly necessary immediately
* and of course, a viewer

* Formats
    - Bincode for data serialization
    - Tracing/Tracing-subscriber for logging
