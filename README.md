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
**Update 7/6**

Next walk direction code for straightest geodesics is almost done. Definitely needs t-values for edge splits (maybe make this an attribute of struct?), but also needs to handle walk direction recomputation where necessary. There are two ways to implement this: at the start of the function and at the end. If you compute at the start, you don't need to precompute to call straightest geodesics the first time. Maybe do everything directional at the start (i.e., compute the next walk face and direction at the start), then simply trace the geodesic on that next face, then pass along to next call. Also need some sort of treatment for degenerate next walk directions. Maybe fail this case and hope sampling density/mesh quality is good enough?

Also, would recommend moving projection info into a separate mesh struct to avoid these huge function inputs. ex. just a struct called meshElements which contains just the type and index, or even just treat is as a pair. This streamlines things and makes it clearer what the variable is doing.

**NOTES**

Important note about Eigen. For Eigen fixed-size containers that are a multiple of 16 in size (ex. Vector2d, Matrix4d), they cannot be directly placed into a std::vector<> if compiling with c++11 or 14, since they need to be placed at fixed sizes in memory. See utils.hpp or mesh.hpp for how to deal with this (i.e., allocator).

This does NOT affect Vector3d, Matrix3d, or dynamic sized (ex. MatrixXd) objects, though. In addition, c++17 handles this implicitly, so no need to handle if using c++17. However, it's good to put the allocators in for fixed-size 16 data types anyways for reliability if you happen to be below c++17.

**Big TODOs**:

- *Polyscope Tests*: For running tests. Re-do polyscope front-end so visual debugging is enabled; this needs its own editable curve network class and converters from the new internal curvenet/dCN classes.

- *Straightest Geodesic*: Yeah... Figure this out.

- *Profile Mover*: Integrate all components into the profilemover class and write function for actual runtime computation. This class should also have the construct the operators, and needs functions for the intermediate stages of computation (i.e., computing per-face deformed polygons, computing deformed vertex projections). Needs a few index maps to get in and out of the matrix indices.

- *Front End*: Figure out how to attach to Blender w/ all options that I want (Blender Bezier's may be a problem). Try Maya later. These both need separate API calls, and conversions from their internal types to a readable input type of profile mover (and vice versa).

- *Speed*: After everything works, find shortcuts and cut down on unnecessary operations. Consider adding bounding boxes for faster projection (how to do this?)

- *Projection Posing*: As mentioned in paper, compute cut-mesh from rest pose, then move all dcurvenet vertices and cut-mesh to the new configuration. Seems like curvenet is unusable like this??? Maybe only if you have pre-defined poses for the curvenet, you can transfer them? Maybe there's a scheme for editing the dcurvenet directly that I can look into?

**Small Steps**:
- t-values for edge insertion on cut-mesh
- Next walk direction recomputation mechanism
- Work on mesh cutting. This should be mosly straightforward with some edge cases to be wary of (ex. boundaries).
- Move proj vector from mesh verts to dCN verts and add a function to pre-compute the deformed proj during runtime.
- Add a function which checks if any edge chains are fully within a face, then removes them (and makes their associated dCN vertices inactive).
- Finish Polyscope front end mesh and curvenet classes.
- Add IO function to visualize edits in polyscope.
- See Appendix for DEC halfedge laplacian shortcut
