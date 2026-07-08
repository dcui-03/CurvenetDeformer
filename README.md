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
**Update 7/8**

Finished my first full draft of mesh cutting... although I think it'll need some extra debugging. However, it also doesn't add new boundary halfedges in, and simply leaves their twins dangling. I'm not sure if this is a problem... but it would definitely be better to have the halfedges than not, so as to prevent the mesh iterators from breaking (we may need them later). Here's my thought: After adding each split, record the start and end boundary halfedge of each (make a new one + new edge if needed), as well as tracking their previous and next. Then, let the start boundary he's prev be the end. Lastly, rewire the remaining boundaries. If a start boundary halfedge's twin's dest is the same as the previous split's end boundary halfedge, then connect them using prev/next. If not, then use the computed previous and next from the previous step. NOTE: if any pair of matching boundaries were caused by a split, then they will both require new halfedges. There will never be a case where we need to rewire with an old boundary and a new boundary. This is made easier by the sequential nature of this function, which fully let's us treat prior splits as boundaries.
NOTE: do NOT rewire the boundary halfedges until AFTER all the splits are done for this vertex. This makes it easier to iterate without fearing problems.

Also, structs should be reorganized so that they're better compartmentalized... right now there's lots of loose variables that can be grouped together. This should also make it easier to switch modes (i.e., deformation vs. color vs. scalar interpolation on cut-mesh)

I think there are some other organizational shortcuts that may be useful. It turns out that we don't really need a mapping from dCN verts to mesh verts, or even dCN halfedges to mesh halfedges. The paper circumvents this entirely by averaging everything onto dCN verts first. Problem is, this means we need to store 2 def grads per dCN verts, which is a bit annoying. Instead, maybe have each copied cut-mesh vertex reference a halfedge in  the dCN mesh. This halfedge + its next (and its destination vert) can be used to define a corner in the dCN, which let's us query corresponding def grads + positions much more easily. dCN vert to mesh vert mappings are therefore "implicit" via dCN halfedge destination. However, note that we absolutely do need cut-mesh halfedges to store their corresponding dCN halfedges. BUT this means that instead of storing a dCN halfedge for each cut-vert, we can actually store one of its adjacent cut-mesh halfedges instead, since these halfedges reference their dCN halfedge. i.e., cut-vert --> cut-halfedge --> dCN halfedge --> dCN corner

**NOTES**

Important note about Eigen. For Eigen fixed-size containers that are a multiple of 16 in size (ex. Vector2d, Matrix4d), they cannot be directly placed into a std::vector<> if compiling with c++11 or 14, since they need to be placed at fixed sizes in memory. See utils.hpp or mesh.hpp for how to deal with this (i.e., allocator).

This does NOT affect Vector3d, Matrix3d, or dynamic sized (ex. MatrixXd) objects, though. In addition, c++17 handles this implicitly, so no need to handle if using c++17. However, it's good to put the allocators in for fixed-size 16 data types anyways for reliability if you happen to be below c++17.

**Big TODOs**:

- *Code Restructuring*: Some major restructures would be nice, but best saved for later. For one, templating the mesh class would make the cut-mesh class easier to interface with and avoid storing a bunch of unused data in the mesh class. Potentially add cutVertex and cutHE to the struct list to accommodate this. The same could be done for the dCurvenet class, where instead of only storing deformation info, it could be used to store a bunch of other info. Combining struct info would make the code much cleaner. Look also into where we can do parallelization. Ex. during runtime, intermediate steps can be pretty cleanly parallelized to assemble all matrices.

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
