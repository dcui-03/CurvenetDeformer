# ProfileMover
Implementation of Pixar's Profile Mover in C++/Python

**Build Instructions**
You only need to clone Polyscope into a folder called deps; Eigen is fetched automatically by CMake, and libigl is no longer a dependency. This has been tested on Windows (MSVC), Linux, and macOS (Apple Silicon).

```
git clone https://github.com/dcui-03/ProfileMover
cd ./ProfileMover
mkdir deps && cd ./deps
git clone --recurse-submodules https://github.com/nmwsharp/polyscope.git
cd ..
```

*macOS / Linux*
```
cmake -S . -B build   -DCMAKE_CXX_COMPILER=$(brew --prefix llvm)/bin/clang++   -DOpenMP_ROOT=$(brew --prefix libomp)
cmake --build build -j4
./build/profile_mover ./data/small_sphere.obj
```
The `-DCMAKE_CXX_COMPILER`/`-DOpenMP_ROOT` flags are only needed on macOS, since AppleClang doesn't ship OpenMP (install both via `brew install llvm libomp`). On Linux, plain `cmake -S . -B build` is enough.

*Windows (MSVC)*
```
cmake -S . -B build
cmake --build build --config Release
.\build\Release\profile_mover.exe .\data\small_sphere.obj
```
MSVC finds OpenMP on its own, so no extra flags are needed.

# Updates and Notes
**Update 8/16**

Major reorganization of files. Cleaned up curvenet, dcurvenet, mesh, and cutmesh classes such that they are barebones structures I can build on. Meanwhile, all their deformation functionality has been moved to profilemover, which they all friend. This lets me cleanly separate out objects so I can reuse them flexibly!

Next few TODO's:
- Projection posing
- Face rigging using the profilerig class
- On the front end, my curve representation needs to be broad enough to handle catmull-rom and cubic bezier. Leave this for Blender/Maya API
- Code cleanup: It's getting messy, let's clean it up...

It's also time to start thinking about the plug-in itself. It looks like Maya and Blender's Bezier curve types are not robust enough to handle multiple adjacent curves to a control. May need to do some front-end magic to make this happen. Thoughts on Houdini as well?

**NOTE ON EIGEN**

Important note about Eigen. For Eigen fixed-size containers that are a multiple of 16 in size (ex. Vector2d, Matrix4d), they cannot be directly placed into a std::vector<> if compiling with c++11 or 14, since they need to be placed at fixed sizes in memory. This can be circumvented, but it's an extra headache and not worth it imo. This is not an issue for c++17, so for simplicity, I recommend we stick with c++17 and above. Note also that this does NOT affect Vector3d, Matrix3d, or dynamic sized (ex. MatrixXd) objects, though.

**Big TODOs**:

- *Code Cleanup*: Change checks into asserts and add error flags. Uniformify naming.

- *Curvenet-Adjacent Features*: (1) Add support for face rigging. This is less of a profilemover class problem, but perhaps a new class called facerig, which combines the mesh class with the curvenet and dcurvenet classes? Basically, add a bunch of robustness to the curvenet/dcurvenet classes and their surrounding features.

- *Robustness*: A couple of extra things to try once its tested and stable: (1) use a better arclength curve sampling strategy than the naive one we currently have --> Do we need to? Naive samples with many more points is potentially just faster and simpler; (2) Preliminary validity checks and unique error signatures (1) Test manifoldness + planarity of faces, 2. Test that no projected vertices collide, 3. Bound the sampling rate by the snapping criteria, 4. Test failure modes of the geodesics step (do we catch most/all errors?) 4. Test that the input curvenets/deformed curvenets do not have degenerate tangents/normals, 5. Test validity of input splines at init AND at runtime i.e., no degenerate splines such as when start and endpoint have the same location and tangent vectors, 6. Test mesh orientability code. 7. Test validity of input meshes and curve networks. Prevent curvenets and meshes with loose verts/edges or clean them up before using.

- *Mid-level Speedups*: (1) Bounding boxes for projections: Create a hierarchical bounding box structure for subsets of faces (binary AABB is probably fine). Then test closest face for any box we are inside or is nearest to us (need to test both if so).

- *Stretch Goals + Major Speedups*: (1) Parallelization with CUDA. There are a bunch of runtime operations (small/medium sized matrix multiplications) that could use this(?) (2) Augment the cutmesh class to accept colors/scalars instead of deformation gradients. This lets us do discontinuous color interpolation, although figuring out how to interface this is hard... (3) Template the mesh class so that we can have various other attributes for free. This one seems tough though...

- *Port to Blender + Maya*: This will require some python front-end and figuring out how to compile and interface.
