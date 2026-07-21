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
**Update 7/18**

I realized that if I switch to Release mode, it's very much realtime! Having OpenMP isn't even necessary, although it can be nice to have. Also tested with CG and it's remarkably MUCH slower, even with capped iterations/tolerance. Also, although I keep calling it "ARAP", the current system is more like a single block-coordinate ARAP solve, meaning we're not quite at a local minimum actually. I think that even on release mode, a true ARAP would be extremely expensive, to the point of becoming non-realtime.

I think that next I should refactor the front-end so that it can do closest-point queries on polygons. The easiest thing would just be to use the custom mesh class as a basis for this. I also need to start thinking about cleanup and shipping (i.e., which components to include like OpenMP, etc.)

Another thing: I need to make this spline-type agnostic. Maybe a good idea is to template the spline class/make a parent spline class that can handle the different variations in curve type (i.e., to support Catmull-Rom splines). It might be hard to allow mixing curve types though, since they each have different needs and therefore sizes. In this same 

It's also time to start thinking about the plug-in itself. It looks like Maya and Blender's Bezier curve types are not robust enough to handle multiple adjacent curves to a control.

**NOTE ON EIGEN**

Important note about Eigen. For Eigen fixed-size containers that are a multiple of 16 in size (ex. Vector2d, Matrix4d), they cannot be directly placed into a std::vector<> if compiling with c++11 or 14, since they need to be placed at fixed sizes in memory. This can be circumvented, but it's an extra headache and not worth it imo. This is not an issue for c++17, so for simplicity, I recommend we stick with c++17 and above. Note also that this does NOT affect Vector3d, Matrix3d, or dynamic sized (ex. MatrixXd) objects, though.

**Big TODOs**:

- *Code Cleanup*: Change checks into asserts and add error flags. Delete unused variables, uniformify naming.

- *Robustness*: A couple of extra things to try once its tested and stable: (1) use a better arclength curve sampling strategy than the naive one we currently have --> Do we need to? Naive samples with many more points is potentially just faster and simpler; (2) Preliminary validity checks and unique error signatures (1) Test manifoldness + planarity of faces, 2. Test that no projected vertices collide, 3. Bound the sampling rate by the snapping criteria, 4. Test failure modes of the geodesics step (do we catch most/all errors?) 4. Test that the input curvenets/deformed curvenets do not have degenerate tangents/normals, 5. Test validity of input splines at init AND at runtime i.e., no degenerate splines such as when start and endpoint have the same location and tangent vectors, 6. Test mesh orientability code. 7. Test validity of input meshes and curve networks. Prevent curvenets and meshes with loose verts/edges or clean them up before using.

- *Mid-level Speedups*: (1) Bounding boxes for projections: Create a hierarchical bounding box structure for subsets of faces (binary AABB is probably fine). Then test closest face for any box we are inside or is nearest to us (need to test both if so). (2) Assemble the halfedge Laplacian operator without building the other operators, if possible. (3) Identify places where we can parallelize and use OpenMP OR even CUDA... (ex. curvenet creation requires a bunch of projection. Could use this here)

- *Stretch Goals/Major Speedups*: (1) Parallelization with CUDA. There are a bunch of runtime operations (small/medium sized matrix multiplications) that could definitely use this. (3) Augment the cutmesh class to accept colors/scalars instead of deformation gradients. This lets us do discontinuous color interpolation, although figuring out how to interface this is hard... (4) Template the mesh class so that we can have various other attributes for free. This one seems tough though... (5) Projection posing mentioned in the paper.

- *Port to Blender + Maya*: This will require some python front-end and figuring out how to compile and interface.
