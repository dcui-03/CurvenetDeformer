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

Fixed the edge start and vertex start bugs! Next, I need to test meshes with boundaries, and also stress test polygonal mesh inputs and non-manifold inputs. ALSO I forgot to debug and clean up the utils.cpp file. There's lots to unused functions there that should just be discarded.

A few things to think about: What happens if we don't enforce the corner normal alignment? Does the code just break? Also maybe I should time each component of the solve... I wonder what's taking up most of the compute?

It looks like the face deformation assembly is the bottleneck --> I wonder if I should just go ahead and make this a precomputed matrix operation. It seems like the problem is likely that each iteration, I need to reindex all adjacent vertices and halfedges, which is probably too expensive. I think, actually, that I should precompute operators for everything and just use those, but do not delete the mesh/cut-mesh objects. During the final release, I can have these "deleted" after creation.

Small TODOs: Change dCN verts to store only a single outgoing halfedge, and then use adjHE to get the rest. Add a safe normalization helper function to the utils so that we can safely normalize or catch errors.

Lastly, another thought on parallelization: We can potentially parallelize the projection function itself by parallelizing projection over faces. NOTE that this is mostly only good when we need to do single point projection queries, as we always want to run the batch of points to project in parallel and we do NOT want nested parallelism for safety reasons. We can potentially do this by enforcing nested parallelism depth of 1, thus making each projection single-threaded when doing a batch, and multi-threaded otherwise.

**NOTE ON EIGEN**

Important note about Eigen. For Eigen fixed-size containers that are a multiple of 16 in size (ex. Vector2d, Matrix4d), they cannot be directly placed into a std::vector<> if compiling with c++11 or 14, since they need to be placed at fixed sizes in memory. This can be circumvented, but it's an extra headache and not worth it imo. This is not an issue for c++17, so for simplicity, I recommend we stick with c++17 and above. Note also that this does NOT affect Vector3d, Matrix3d, or dynamic sized (ex. MatrixXd) objects, though.

**Big TODOs**:

- *Code Cleanup*: Change checks into asserts and add error flags. Delete unused variables, uniformify naming.

- *Robustness*: A couple of extra things to try once its tested and stable: (1) use a better arclength curve sampling strategy than the naive one we currently have; (2) Preliminary validity checks and unique error signatures (1) Test manifoldness + planarity of faces, 2. Test that no projected vertices collide, 3. Bound the sampling rate by the snapping criteria, 4. Test failure modes of the geodesics step (do we catch most/all errors?), 5. Test validity of input splines at init AND at runtime i.e., no degenerate splines such as when start and endpoint have the same location and tangent vectors, 6. Test mesh orientability code. 7. Test validity of input meshes and curve networks. Prevent curvenets and meshes with loose verts/edges or clean them up before using

- *Mid-level Speedups*: (1) Bounding boxes for projections: Create a hierarchical bounding box structure for subsets of faces (binary is probably fine). Then test closest face for any box we are inside or is nearest to us. (2) Assemble the halfedge Laplacian operator without building the other operators, if possible. (3) Identify places where we can parallelize and use OpenMP OR even CUDA... (ex. curvenet creation requires a bunch of projection. Could use this here)

- *Stretch Goals/Major Speedups*: (1) Parallelization with CUDA (2) After computing the cutmesh operators (and more), store those operators and values in profilemover, then discard both the mesh and cutmesh class. (Almost) everything then becomes matrix operations! Discarding the cutmesh class is a little bit scary for debugging though... maybe make a branch to test this out? It's probably better for the Blender/Maya ports anyways. I wrote the matrices out more explicitly in Goodnotes. (3) ARAP. This would definitely not be realtime, but it would be extremely cool! (4) Augment the cutmesh class to accept colors/scalars instead of deformation gradients. This lets us do discontinuous color interpolation, although figuring out how to interface this is hard... (5) Template the mesh class so that we can have various other attributes for free. This one seems tough though... (6) Projection posing mentioned in the paper.

- *Port to Blender + Maya*: This will require some python front-end and figuring out how to compile and interface.
