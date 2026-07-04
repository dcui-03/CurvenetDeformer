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
**Update 7/3**

Major update, mesh class and cut-mesh classes were once again split ugh. Cut-mesh class is now a child class of the mesh class and has its own set of functions to augment the mesh class's. This is to accommodate the geodesic computation and make expensive vertex projection parallelizable. It's now more accurate to the paper! But still needs a lotttt of reorganization to clear up issues. Also had classes friend each other to make inter-op easier.

Restarting straightest geodesics now, which operates on the original mesh and not on the cut-mesh. Shouldn't be horribly difficult given the new setup (famous last words). However, likely need to move this to the mesh class and then call it from the cut-mesh class when necessary.

**NOTES**

Important note about Eigen. For Eigen fixed-size containers that are a multiple of 16 in size (ex. Vector2d, Matrix4d), they cannot be directly placed into a std::vector<> if compiling with c++11 or 14, since they need to be placed at fixed sizes in memory. See utils.hpp or mesh.hpp for how to deal with this (i.e., allocator).

This does NOT affect Vector3d, Matrix3d, or dynamic sized (ex. MatrixXd) objects, though. In addition, c++17 handles this implicitly, so no need to handle if using c++17. However, it's good to put the allocators in for fixed-size 16 data types anyways for reliability if you happen to be below c++17.

**Big TODOs**:

- *Polyscope Tests*: For running tests. Re-do polyscope front-end so visual debugging is enabled; this needs its own editable curve network class and converters from the new internal curvenet/dCN classes.

- *Straightest Geodesic*: Yeah... Figure

- *Profile Mover*: Integrate all components into the profilemover class and write function for actual runtime computation. This class should also have the construct the operators, and needs functions for the intermediate stages of computation (i.e., computing per-face deformed polygons, computing deformed vertex projections). Needs a few index maps to get in and out of the matrix indices.

- *Front End*: Figure out how to attach to Blender w/ all options that I want (Blender Bezier's may be a problem). Try Maya later. These both need separate API calls, and conversions from their internal types to a readable input type of profile mover (and vice versa).

- *Speed*: After everything works, find shortcuts and cut down on unnecessary operations. Consider adding bounding boxes for faster projection (how to do this?)

- *Projection Posing*: As mentioned in paper, compute cut-mesh from rest pose, then move all dcurvenet vertices and cut-mesh to the new configuration. Seems like curvenet is unusable like this??? Maybe only if you have pre-defined poses for the curvenet, you can transfer them? Maybe there's a scheme for editing the dcurvenet directly that I can look into?

**Small Steps**:
- Clean up mesh class (needs a new initialization, iterators, etc.) and cut-mesh initialization
- Finish up mesh embedding. Currently missing a few parts, like t-vals for edge sorting, and halfedge sorting on the tangent plane.
- Work on mesh cutting. This should be mosly straightforward with some edge cases to be wary of (ex. boundaries).
- Finish Straightest Geodesics. Start with fast version and then move to slow version. Should be relatively well-structured at this point.
- Add a function which checks if any edge chains are fully within a face, then removes them (and makes their associated dCN vertices inactive).
- Finish Polyscope front end mesh and curvenet classes.
- Add IO function to visualize edits in polyscope.
- See Appendix for DEC halfedge laplacian shortcut
