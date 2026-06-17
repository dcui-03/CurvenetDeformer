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
**Update 6/15**

Realized the original paper has "curves" that are distinct from splines (i.e., chains of splines such that the endpoints are intersections/anchors). Now knee-deep in dCurvenet initialization. Curvenet can be converted into a dcurvenet now, but still drafting functions for computing scaled local frames.

**NOTES**

Important note about Eigen. For Eigen fixed-size containers that are a multiple of 16 in size (ex. Vector2d, Matrix4d), they cannot be directly placed into a std::vector<> if compiling with c++11 or 14, since they need to be placed at fixed sizes in memory. See utils.hpp or mesh.hpp for how to deal with this (i.e., allocator).

This does NOT affect Vector3d, Matrix3d, or dynamic sized (ex. MatrixXd) objects, though. In addition, c++17 handles this implicitly, so no need to handle if using c++17. However, it's good to put the allocators in for fixed-size 16 data types anyways for reliability.

**Big TODOs**:
- *Discrete Curvenet*: Re-do discrete curvenet class

- *Straightest Geodesic*: Re-do straightest geodesics. Include a fast (same-face check) version, and a slower (visibility-check) version. The first should be fairly straightforward to implement. The second will require extra Utils functions. Maybe look at other implementations for inspo? With face normals, should be fairly straightforward to just apply explicit rotations to a unit direction vector.

- *Profile Mover*: Integrate all components into the profilemover class and write function for actual runtime computation.

- *Front End*: Figure out how to attach to Blender w/ all options that I want (Blender Bezier's may be a problem). Try Maya later.

- *Speed*: After everything works, find shortcuts and cut down on unnecessary operations. Consider adding bounding boxes for faster projection (how to do this?)

- *Projection Posing*: As mentioned in paper, compute cut-mesh from rest pose, then move all dcurvenet vertices and cut-mesh to the new configuration. Seems like curvenet is unusable like this??? Maybe only if you have pre-defined poses for the curvenet, you can transfer them? Maybe there's a scheme for editing the dcurvenet directly that I can look into?

**Small Steps**:

- Add a function which checks if any edge chains are fully within a face, then removes them (and makes their associated dCN vertices inactive).
- Implement the dCN scaled local frames.
- Check polygon DEC operators and try to find better shortcuts for constructing them (check Appendix?)
- Add IO function to visualize mesh edits in polyscope and debug mesh class.
- Save an initialized spline curvenetwork to some sparse file type (OBJ for splines???) so they can be loaded and reused.
