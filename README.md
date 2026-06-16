# ProfileMover
Implementation of Pixar's Profile Mover in C++/Python

You will need to clone Polyscope into a folder called deps. I also haven't tested this setup on anything but Apple Silicon so the Cmakelists.txt might need some additional lines to support Windows.
Polyscope and libigl deps planned to be removed in the future.

**Update 6/15**

Re-written curvenet class stripped to the bare minimum. The paradigm for the curvenet class is that it takes an existing curve network of bezier splines from some other editing software and statically builds it, meaning we don't need too many complicated functions for adding and removing, etc. These are the job of the front-end bezier curve data structure.

**NOTES**

Important note about Eigen. For Eigen fixed-size containers that are a multiple of 16 in size (ex. Vector2d, Matrix4d), they cannot be directly placed into a std::vector<> if compiling with c++11 or 14, since they need to be placed at fixed sizes in memory. See utils.hpp or mesh.hpp for how to deal with this (i.e., allocator).

This does NOT affect Vector3d, Matrix3d, or dynamic sized (ex. MatrixXd) objects, though. In addition, c++17 handles this implicitly, so no need to handle if using c++17. However, it's good to put the allocators in for fixed-size 16 data types anyways for reliability.

**Big TODO's**:

- *Discrete Curvenet*: Re-do discrete curvenet class. Treat as a halfedge data structure, which makes indices into the dCN on the cut-mesh easier.

- *Straightest Geodesic*: Re-do straightest geodesics. Include a fast (same-face check) version, and a slower (visibility-check) version. The first should be fairly straightforward to implement. The second will require extra Utils functions. Maybe look at other implementations for inspo? With face normals, should be fairly straightforward to just apply explicit rotations to a unit direction vector.

- *Polygon DEC*: Check polygon operators and try to find better shortcuts for constructing them (check Appendix?)

- *Profile Mover*: Integrate all components into the profilemover class and write function for actual runtime computation. 

- *Debugging and Cleaning*: Try to avoid constructing 2D basis + redundancy.

**Small Steps**:

- Debug and check the 

- Re-implement discrete curvenet as a "half-edge" type data structure where each "edge" has two directional halfedges and each halfedge can store its own scaled frame. Verts, half edges, edges, and splines can all be structs in this case, where each points to the others in a principled way. Provide functions for iterators (vertAdjHEs, vertAdjSplines), basic spline insertion, vertPairToHE map.

- Add IO function to visualize mesh edits in polyscope and debug mesh class.
