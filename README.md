# ProfileMover
Implementation of Pixar's Profile Mover in C++/Python

You will need to clone Polyscope into a folder called deps. I also haven't tested this setup on anything but Apple Silicon so the Cmakelists.txt might need some additional lines to support Windows.
Polyscope and libigl deps planned to be removed in the future.


**Update 6/14**

Major restructuring of the Mesh class to combine the half edge data structure. This is to accommodate the iterative insertion of curvenetwork vertices into the mesh during cut-mesh creation. Created functions insertVertex(), splitEdge(), and insertEdge(), which are black boxes for inserting a new vertex or edge(s) into a mesh. This should make the projection + geodesic computation much, much easier.

As a result, cut-mesh class, straightest geodesic custom files, and projected discrete curvenet were completely removed. These are/will be folded into the mesh class.

**Big TODO's**:

- *Discrete Curvenet*: Re-do discrete curvenet class. Treat as a halfedge data structure, which makes indices into the dCN on the cut-mesh easier.

- *Straightest Geodesic*: Re-do straightest geodesics. Include a fast (same-face check) version, and a slower (visibility-check) version. The first should be fairly straightforward to implement. The second will require extra Utils functions. Maybe look at other implementations for inspo? With face normals, should be fairly straightforward to just apply explicit rotations to a unit direction vector.

- *Polygon DEC*: Check polygon operators and try to find better shortcuts for constructing them (check Appendix?)

- *Profile Mover*: Integrate all components into the profilemover class and write function for actual runtime computation. 

- *Debugging and Cleaning*: Try to avoid constructing 2D basis + redundancy.

**Small Steps**:

- Clean up curvenet class again. Reduce to structs and strip to bare minimum functions. Most important part of curvenet class is fast sampling, discretization.

- Re-implement discrete curvenet as a "half-edge" type data structure where each "edge" has two directional halfedges and each halfedge can store its own scaled frame. Verts, half edges, edges, and splines can all be structs in this case, where each points to the others in a principled way. Provide functions for iterators (vertAdjHEs, vertAdjSplines), basic spline insertion, vertPairToHE map.

- Add IO function to visualize mesh edits in polyscope and debug mesh class.
