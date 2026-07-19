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

Ok, I matricized things and added IO functionality to the pscurvenet. The good thing is that matricizing speeds up face deformation by half! The bad news is that everything else takes pretty much the same amount of time, so overall it's still not "real-time" on larger meshes. So unless I incorporate something even faster, like CUDA, I'm pretty much capped here.

Turns out I did not in fact have OpenMP turned on. It's now fixed (changed CMakelists.txt). This does make it faster (cuts matrix compute time AND face def grads by half!), although I'm learning a few things about openmp that may be troublesome... Definitely check TODO's, esp in cutmesh and dcurvenet for more details. May still consider (1) CUDA for some computations and (2) an iterative solver like conjugate gradient instead of cholesky.

Next, I need to test meshes with boundaries, and also stress test polygonal mesh inputs and non-manifold inputs. ALSO I forgot to debug and clean up the utils.cpp file. There's lots to unused functions there that should just be discarded. After this, the big thing is cleanup and robustness. The code is actually quite messy right now and needs a lot more error messages + asserts.

Another question: What happens if we don't enforce the corner normal alignment? Is that better? I do notice there actually is that discontinuous "snapping" when we flip.

**NOTE ON EIGEN**

Important note about Eigen. For Eigen fixed-size containers that are a multiple of 16 in size (ex. Vector2d, Matrix4d), they cannot be directly placed into a std::vector<> if compiling with c++11 or 14, since they need to be placed at fixed sizes in memory. This can be circumvented, but it's an extra headache and not worth it imo. This is not an issue for c++17, so for simplicity, I recommend we stick with c++17 and above. Note also that this does NOT affect Vector3d, Matrix3d, or dynamic sized (ex. MatrixXd) objects, though.

**Big TODOs**:

- *Code Cleanup*: Change checks into asserts and add error flags. Delete unused variables, uniformify naming.

- *Robustness*: A couple of extra things to try once its tested and stable: (1) use a better arclength curve sampling strategy than the naive one we currently have --> Do we need to? Naive samples with many more points is potentially just faster and simpler; (2) Preliminary validity checks and unique error signatures (1) Test manifoldness + planarity of faces, 2. Test that no projected vertices collide, 3. Bound the sampling rate by the snapping criteria, 4. Test failure modes of the geodesics step (do we catch most/all errors?) 4. Test that the input curvenets/deformed curvenets do not have degenerate tangents/normals, 5. Test validity of input splines at init AND at runtime i.e., no degenerate splines such as when start and endpoint have the same location and tangent vectors, 6. Test mesh orientability code. 7. Test validity of input meshes and curve networks. Prevent curvenets and meshes with loose verts/edges or clean them up before using.

- *Mid-level Speedups*: (1) Bounding boxes for projections: Create a hierarchical bounding box structure for subsets of faces (binary AABB is probably fine). Then test closest face for any box we are inside or is nearest to us (need to test both if so). (2) Assemble the halfedge Laplacian operator without building the other operators, if possible. (3) Identify places where we can parallelize and use OpenMP OR even CUDA... (ex. curvenet creation requires a bunch of projection. Could use this here)

- *Stretch Goals/Major Speedups*: (1) Parallelization with CUDA. There are a bunch of runtime operations (small/medium sized matrix multiplications) that could definitely use this. (3) Augment the cutmesh class to accept colors/scalars instead of deformation gradients. This lets us do discontinuous color interpolation, although figuring out how to interface this is hard... (4) Template the mesh class so that we can have various other attributes for free. This one seems tough though... (5) Projection posing mentioned in the paper.

- *Port to Blender + Maya*: This will require some python front-end and figuring out how to compile and interface.
