# ProfileMover
Implementation of Pixar's Profile Mover in C++/Python

**Build Instructions**
You will need to clone Polyscope into a folder called deps. I also haven't tested this setup on anything but Apple Silicon so the Cmakelists.txt might need some additional lines to support Windows.
Polyscope and libigl deps planned to be removed in the future.

```
git clone https://github.com/dcui-03/ProfileMover
cd ./ProfileMover
mkdir deps && cd ./deps
git clone --recurse-submodules https://github.com/nmwsharp/polyscope.git
git clone --recurse-submodules https://github.com/libigl/libigl.git
cd .. && mkdir build && cd ./build
cmake ..
make -j4
cd ..
./build/profile_mover ./data/small_sphere.obj
```

# Updates and Notes
**Update 7/9**

Finished the mesh splitting code. Yay!!!!! BUT haven't checked it yet so there could be logical bugs. Once the entire straightest geodesics pipeline is checked, I should probably do a bunch of reorganization and compartmentalization to start making modularity a priority (see below). After that, it's FINALLY time to finally move on to system assembly (i.e., the profilemover class).

More practically speaking, adding polyscope test support for visual tests is also loooong overdue. This would be a good thing to work on when I'm tired of other stuff, although maintaining polyscope-consistent spline indexing AND allowing for various tools to break, merge, and restructure splines will be quite the challenge. This will definitely take its own day.

Also, structs should be reorganized so that they're better compartmentalized... right now there's lots of loose variables that can be grouped together. This should also make it easier to switch modes (i.e., deformation vs. color (RGBA?) vs. scalar interpolation on cut-mesh).

**NOTES**

Important note about Eigen. For Eigen fixed-size containers that are a multiple of 16 in size (ex. Vector2d, Matrix4d), they cannot be directly placed into a std::vector<> if compiling with c++11 or 14, since they need to be placed at fixed sizes in memory. See utils.hpp or mesh.hpp for how to deal with this (i.e., allocator).

This does NOT affect Vector3d, Matrix3d, or dynamic sized (ex. MatrixXd) objects, though. In addition, c++17 handles this implicitly, so no need to handle if using c++17. However, it's good to put the allocators in for fixed-size 16 data types anyways for reliability if you happen to be below c++17.

**Big TODOs**:

- *Code Restructuring*: Some major restructures would be nice, but best saved for later. For one, templating the mesh class would make the cut-mesh class easier to interface with and avoid storing a bunch of unused data in the mesh class. Potentially add cutVertex and cutHE to the struct list to accommodate this. The same could be done for the dCurvenet class, where instead of only storing deformation info, it could be used to store a bunch of other info. Combining struct info would make the code much cleaner. Look also into where we can do parallelization. Ex. during runtime, intermediate steps can be pretty cleanly parallelized to assemble all matrices, spline discretization (both splines themself and sampling).

- *Polyscope Tests*: For running tests. Re-do polyscope front-end so visual debugging is enabled; this needs its own editable curve network class and converters from the new internal curvenet/dCN classes.

- *Profile Mover*: Integrate all components into the profilemover class and write function for actual runtime computation. This class should also have the construct the operators, and needs functions for the intermediate stages of computation (i.e., computing per-face deformed polygons, computing deformed vertex projections). Needs a few index maps to get in and out of the matrix indices.

- *Front End*: Figure out how to attach to Blender w/ all options that I want (Blender Bezier's may be a problem). Try Maya later. These both need separate API calls, and conversions from their internal types to a readable input type of profile mover (and vice versa).

- *Loose Ends*: A couple of extra things to try once its tested and stable: (1) use a better arclength curve sampling strategy than the naive one we currently have; (2) Maybe figure out how to do bounding boxes for faster projections? Not sure how this works... (3) Preliminary validity checks and unique error signatures (1. Test manifoldness + planarity of faces, 2. Test that no projected vertices collide, 3. Bound the sampling rate by the snapping criteria, 4. Test failure of the geodesics step, 5. Test validity of input splines at init AND at runtime i.e., no degenerate splines such as when start and endpoint have the same location and tangent vectors).

- *Projection Posing*: As mentioned in paper, compute cut-mesh from rest pose, then move all dcurvenet vertices and cut-mesh to the new configuration. Seems like curvenet is unusable like this??? Maybe only if you have pre-defined poses for the curvenet, you can transfer them? Maybe there's a scheme for editing the dcurvenet directly that I can look into?

**Small Steps**:
- Check entire geodesics pipeline, from projection to cutting
- Restructure code. Minor pointers: check that the sampling rate for curvenet curves is correct, and that the num_samples are actually applied to the curves. Also check to make sure that dCN anchors' halfedges are next and prev instead of dangling.
- Move proj vector from mesh verts to dCN verts and add a function to pre-compute the deformed proj during runtime.
- Add a function which checks if any edge chains are fully within a face, then removes them (and makes their associated dCN vertices inactive).
- Finish Polyscope front end mesh and curvenet classes.
- Add IO function to visualize edits in polyscope.
- See Appendix for DEC halfedge laplacian shortcut
