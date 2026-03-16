# ProfileMover
Implementation of Pixar's Profile Mover in C++/Python

You will need to clone Polyscope into `deps/polyscope`. The current CMake setup is configured for macOS Apple Silicon and Linux/WSL.


**Update 3/12**

Added a Half-edge data structure called polyHE_t (namespace polyHE) which is modified from (https://github.com/yig/halfedge) and allows for edge querying. This comprises the topological component of the mesh class; vertices and faces are also stored separately in the mesh class for ease of access. Added some utility functions to utils.hpp, which are used by the mesh class. DECUtils should also be essentially complete. Curvenet class should also be complete, with some minor tweaks and cleaning up.

Note: there's a function for projecting a vertex onto the mesh in mesh.hpp. There's also a function for projecting a point onto a tangent plane in utils.hpp, which requires that you define a normal and a center (ex. face normal and face barycenter), as well as one for projecting a vector onto a tangent plane. Since we may have non-planar faces,I'm using vector area as the normal for faces, which is already precomputed for each face using the mesh class's fNormals list.

**Next Big Steps**:

- *Discrete Curvenet*: Integrate Leo's frame computation into the dcurvenet class

- *Projected Discrete Curvenet*: Figure out what extra labeling is needed

There is also a utils folders for miscellaneous algorithms (DEC, straightest Geodesic, other) each with their own namespaces. (I'm a bit crazy about categorizing things, but better to be to organized than not enough ig...)

## Build and Run (Current)

### 1. Dependencies

You need:
- CMake (>= 3.20)
- A C++17 compiler
- Eigen3
- nlohmann/json
- X11/OpenGL development libs (for Polyscope/GLFW on Linux/WSL)

On Ubuntu/WSL, for example:

```bash
sudo apt update
sudo apt install -y \
  build-essential \
  cmake \
  pkg-config \
  libeigen3-dev \
  nlohmann-json3-dev \
  libgl1-mesa-dev \
  libglu1-mesa-dev \
  xorg-dev
```

### 2. Dependency folder

`deps/polyscope` must exist in this repo.

```bash
git clone --recursive https://github.com/nmwsharp/polyscope.git deps/polyscope
```

### 3. Configure + build

From repo root:

```bash
CC=/usr/bin/gcc CXX=/usr/bin/g++ cmake -S . -B build-ui -DPROFILEMOVER_WITH_POLYSCOPE=ON
cmake --build build -j
```

Notes:
- For Conda environments, explicitly setting `CC/CXX` to system compilers avoids X11 header/toolchain conflicts.

### 4. Run

```bash
./build/profile_mover data/sphere.obj data/sphere-curves.json
```

Arguments:
- `argv[1]`: mesh `.obj` path (required)
- `argv[2]`: curve `.json` path (optional; defaults to `../curvenet/data/sphere-curves.json`)
- *Cut-mesh*: Figure out initialization via cutting using the mesh and projected dcurvenet as input

- *Straightest Geodesic*: Plan out straightest geodesic with consideration for non-planar faces (i.e., some sort of recursion with the *to* vector reset at each iteration)

- *Profile Mover*: Integrate all components into the profilemover class


**Small Steps**:

- Figure out projection onto non-planar faces. I'm guessing a good approximation is Newell plane projection plus mean value coordinates to compute the height field to lift back into 3D.

- Figure out straightest geodesics for non-planar faces. Maybe also a rotation to match previous normal (this would work for vertex/edge normals too) plus Newell plane projection? Needs a lot more thinking through...

- Labeling needs some thinking through in order to make code clean but also lightweight. Does dControl need to be a child class of dVert?
