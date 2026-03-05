# ProfileFormer
Implementation of Pixar's Profile Former in C++/Python

You will need to clone libigl and Polyscope into a folder called deps. I also haven't tested this setup on anything but Apple Silicon so the Cmakelists.txt might need some additional lines to support Windows.


**Update 3/4**

I haven't compiled this yet so it's probably buggy and doesn't run...
Created template headers for all data structures (ProfileFormer, Curvenet, Cutmesh, DCurvenet, Cutmesh), such that it should be relatively straightforward to start filling in some functions. Each of these big classes has its own folder and namespace (in parentheses):

**Profile Former** (ProfileFormer): The main class the algorithm will be operating under.

**Curvenet** (Curvenet): Has its own class (*curvenet*), plus a set of components that each has its own class (*control, tangent, spline*).

**Discrete Curvenet** (DCurvenet): Has its own class (*dcurvenet*), plus a set of components that each has its own class (*dvert, dsegment, dspline*). It also has a derived class *pdcurvenet* (projected discrete curvenet) which inherits most of the structure of the *dcurvenet*, with some additional functionality for the projection.

**Cut-Mesh** (CutMesh): Has its own class (*cutmesh*).

There is also a utils folders for miscellaneous algorithms (DEC, straightest Geodesic, other) each with their own namespaces. (I'm a bit crazy about categorizing things, but better to be to organized than not enough ig...)
