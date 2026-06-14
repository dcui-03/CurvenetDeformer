# ProfileMover
Implementation of Pixar's Profile Mover in C++/Python

You will need to clone Polyscope into a folder called deps. I also haven't tested this setup on anything but Apple Silicon so the Cmakelists.txt might need some additional lines to support Windows.
Polyscope and libigl deps planned to be removed in the future.


**Update 6/14**

Major restructuring of the Mesh class to combine the half edge data structure. This is to accommodate the iterative insertion of curvenetwork vertices into the mesh during cut-mesh creation.

As a result, removed cut-mesh class, straightest geodesic custom files, and projected discrete curvenet. These will be folded into the mesh class.

**Next Big Steps**:

- *Discrete Curvenet*: Look back over curvenet representation, frame computation, 

- *Straightest Geodesic*: Re-do straightest geodesics. Include a fast (same-face check) version, and a slower (visibility-check) version. The first should be fairly straightforward to implement. The second will require extra Utils functions.

- *Polygon DEC*: Check polygon operators and try to find better shortcuts to constructing them (check Appendix?)

- *Profile Mover*: Integrate all components into the profilemover class and write function for actual runtime computation.

- *Debugging and Cleaning*: Try to avoid constructing 2D basis + redundancy.

**Small Steps**:

- Straightest geodesic computation. Optimize using new HE structure; avoid building local bases unless necessary.

- Curvenet and discrete curvenet need to be refreshed. Rotation-minimizing frames and scaled local frames needs to be looked over again. Maybe also use more structs and fewer custom subclasses.

- Add IO function to visualize mesh edits in polyscope and debug mesh class.
