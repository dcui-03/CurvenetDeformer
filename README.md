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
**Update 7/14**

Finished checking mesh projection and geodesics code. Moving on to cutmesh.

An oddity to think about: Currently I forcibly align the corner normals in dCN frame computation with the normal at the intersection via dot product test (i.e., when the signed angle at the corner is larger than 180, it can cause the frame normal to "flip" w.r.t., the intersection normal, so we fix it via a sign change). This leads to pretty frames where pos. and neg. normals generally face the same way, but because this flip is binary, it can cause the frame to flip discontinuously as splines move. BUT if we don't do this, then the splines are no longer oriented w.r.t. the local normal, and therefore some frames don't intuitively describe local "stretch" on the surface anymore, which may also cause problems. Not sure what to do here, maybe change after doing testing.

TODO later: in dCN class, change adjHE to a single stored halfedge, and then add an iterator to get all adjacent halfedges; add a safe normalization function to check certain that inputs will not normalize to a NaN or 0.

One big thought: We technically don't need to store the mesh or the cutmesh after computation... once we have operators, it's entirely possible to just convert everything to operators and solve as a bunch of matrix multiplications, although assembling the matrices could be a challenge... The only step that seems non-matricizable is folding the flattened def grads into 3x3 def grads. There likely is a pair of operators that could do this though.....

Lastly, another thought on parallelization: We can potentially parallelize the projection function itself by parallelizing projection over faces. NOTE that this is mostly only good when we need to do single point projection queries, as we always want to run the batch of points to project in parallel and we do NOT want nested parallelism for safety reasons. We can potentially do this by enforcing nested parallelism depth of 1, thus making each projection single-threaded when doing a batch, and multi-threaded otherwise.

**NOTE ON EIGEN**

Important note about Eigen. For Eigen fixed-size containers that are a multiple of 16 in size (ex. Vector2d, Matrix4d), they cannot be directly placed into a std::vector<> if compiling with c++11 or 14, since they need to be placed at fixed sizes in memory. This can be circumvented, but it's an extra headache and not worth it imo. This is not an issue for c++17, so for simplicity, I recommend we stick with c++17 and above. Note also that this does NOT affect Vector3d, Matrix3d, or dynamic sized (ex. MatrixXd) objects, though.

**Big TODOs**:

- *Code Restructuring*: Templating mesh class for easy interfacing and reducing wasted storage? Look into parallelization options (are they even worth it?)

- *Code Cleanup*: Change checks into asserts and add error flags. Delete unused variables, uniformify naming.

- *Front End*: Figure out how to attach to Blender w/ all options that I want (Blender Bezier's may be a problem). Try Maya later. These both need separate API calls, and conversions from their internal types to a readable input type of profile mover (and vice versa).

- *Loose Ends*: A couple of extra things to try once its tested and stable: (1) use a better arclength curve sampling strategy than the naive one we currently have; (2) Maybe figure out how to do bounding boxes for faster projections? Not sure how this works... -> hierarchically find the nearest bounding box, but how do we deal with being inside one? (3) Preliminary validity checks and unique error signatures (1. Test manifoldness + planarity of faces, 2. Test that no projected vertices collide, 3. Bound the sampling rate by the snapping criteria, 4. Test failure modes of the geodesics step (do we catch most/all errors?), 5. Test validity of input splines at init AND at runtime i.e., no degenerate splines such as when start and endpoint have the same location and tangent vectors, 6. Maybe test if input mesh is made up of orientable components, and if so, orient the closed shapes s.t. normals are facing out --> but does the exact normal direc. really matter?). 7. Test validity of input meshes and curve networks. Prevent curvenets and meshes with loose verts/edges

- *Projection Posing*: As mentioned in paper, compute cut-mesh from rest pose, then move all dcurvenet vertices and cut-mesh to the new configuration. Seems like curvenet is unusable like this??? Maybe only if you have pre-defined poses for the curvenet, you can transfer them? Maybe there's a scheme for editing the dcurvenet directly that I can look into?

**Small Steps**:
- Check embedding and cutting
- Check that the sampling rate for curvenet curves is correct.
- Visualize cutmesh in polyscope
- See Appendix for DEC halfedge laplacian shortcut
