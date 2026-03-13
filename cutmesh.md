# Mesh Cutting – Explanation

This section explains how a **curvenet** (a network of surface curves) is embedded into a character mesh by **cutting the mesh along the curves**.

The goal is to build a **cut-mesh** where the curvenet becomes part of the mesh topology, so that quantities can later be treated differently on each side of a curve.

---

## High-Level Idea

We start with:

- a **surface mesh** of the character
- a **curvenet** made of piecewise-linear curve segments lying on the surface

The algorithm constructs a new **cut-mesh** where:

- original mesh vertices are preserved
- curve samples become mesh vertices
- curve segments become mesh edges

That means the mesh is explicitly split wherever the curves pass through it.

---

## Step 1 — Initialize the Cut-Mesh

The process begins by making a copy of the input mesh.

So at first:

`cut-mesh = original mesh`

Then the algorithm progressively inserts new vertices and edges to encode the curvenet inside the mesh connectivity.

---

## Step 2 — Project Curvenet Samples onto the Mesh

Each curve is represented by a sequence of **samples**.

For each sample, the neutral-space position is projected to the **closest point on the surface mesh**.

After projection, each sample is classified as one of three types:

### 1. Vertex Sample
The projected sample lands on an existing mesh vertex.

Example:

    curve sample
         *
         |
    mesh vertex *

In this case, no new geometric point is needed; the existing mesh vertex is used.

### 2. Edge Sample
The sample lands somewhere along a mesh edge.

Example:

    v1 ---------*--------- v2
                 sample

In this case, the edge must be subdivided by inserting a new cut-vertex at the sample location.

### 3. Face Sample
The sample lands in the interior of a mesh face.

Example:

    +-------------+
    |      *      |
    |    sample   |
    +-------------+

In this case, a new isolated cut-vertex is created inside the face.

---

## Step 3 — Insert Cut Vertices

The action taken depends on the sample type:

- **vertex sample**: tag the existing mesh vertex
- **edge sample**: split the edge and insert a new cut-vertex
- **face sample**: create a new cut-vertex inside the face

Example of edge splitting:

Before:

    v1 ----------- v2

After:

    v1 -----*----- v2
            sample

The paper also stores a **bitmask** on each cut-vertex that says whether it corresponds to:

- an original mesh vertex
- a curvenet sample
- a segment-edge intersection

This helps track what kind of point each cut-vertex represents.

---

## Step 4 — Insert Curvenet Segments

After the sample vertices are in place, the algorithm inserts the **segments connecting them**.

This is where the actual mesh cutting happens.

The paper describes several cases.

### Case 1 — Segment Already Lies on an Existing Edge

If both samples lie on the same mesh edge, the segment may simply coincide with that edge.

Example:

    v1 ---- sample ---- sample ---- v2

Then the cut-edge is just subdivided as needed. This is the simplest case.

### Case 2 — Segment Cuts Across a Face

If the segment connects points that lie on the boundary of a face or inside the same face, then a new cut-edge is inserted through that face.

Example:

Before:

    +---------+
    |         |
    |         |
    +---------+

After:

    +----/----+
    |   /     |
    |  /      |
    +---------+

That new edge will later be used to split the face into smaller faces.

### Case 3 — Crack Inside a Face

A segment may lie entirely inside a face, with both endpoints also inside that face.

Example:

    +---------+
    |  *---*  |
    |         |
    +---------+

This creates a kind of **internal crack** in the polygon.

The paper allows this by using a customized halfedge structure where an internal cut-edge can branch inside a face.

### Case 4 — Segment Crosses Multiple Faces

Sometimes the two endpoints of a segment are not in the same face.

In that case, the segment must be traced across the surface mesh.

The paper uses the method of **Polthier and Schmies (1998)** to compute a **straightest path** over the surface mesh between the projected sample points.

As that path crosses mesh edges, new cut-vertices are inserted at the intersections.

So the single curvenet segment becomes a chain of cut-edges across multiple faces.

Conceptually:

    face1 -> face2 -> face3 -> ...

with new vertices placed where the path crosses mesh edges.

---

## Step 5 — Update Mesh Connectivity

At this stage, the algorithm has inserted the necessary vertices and edges, but the face connectivity is not yet fully rebuilt.

Now it must determine:

- which edges bound which faces
- how each original face has been split into smaller cut-faces

To do this, the paper computes a local **tangent space** at each cut-vertex.

Incident halfedges are projected into that tangent plane, then sorted counter-clockwise.

This gives a consistent circular ordering of edges around each cut-vertex, which is needed to reconstruct the new face loops.

---

## Step 6 — Compute a Tangent Space at Each Cut-Vertex

Because the cut-mesh lies on a curved surface, "sorting edges around a vertex" is not trivial in 3D.

So the paper defines a local 2D tangent space for each cut-vertex.

How that tangent space is chosen depends on where the cut-vertex lies.

### Cut-Vertex Inside a Mesh Face

If the cut-vertex lies inside a face, its tangent space is the plane orthogonal to the face normal.

So locally, the face is treated as a flat patch.

### Cut-Vertex on a Mesh Edge

If the cut-vertex lies on an edge shared by two faces, the paper unfolds those two faces into a common plane.

Example:

    face A
       \
        edge
       /
    face B

This allows the incident edges to be treated consistently in 2D.

### Cut-Vertex at a Mesh Vertex

If the cut-vertex lies exactly at an original mesh vertex, the paper flattens the entire **one-ring** of faces around that vertex.

This produces a local planar parameterization of the neighborhood, again so incident halfedges can be ordered.

---

## Step 7 — Sort Halfedges and Rebuild Faces

Once the incident halfedges at each cut-vertex are projected into the tangent plane, they are sorted counter-clockwise.

Then the algorithm walks from one halfedge to the next around the mesh to identify loops.

Each loop corresponds to a new cut-face.

Conceptually:

    edge -> next edge -> next edge -> ... -> back to start

When a loop closes, a new face is created.

This is the step that turns the partially cut structure into a fully connected cut-mesh.

---

## Step 8 — Handle Curvenet Islands Inside a Single Face

There is a special case where a cluster of curvenet segments lies entirely inside one mesh face.

Example:

    +---------+
    |  o---o  |
    |  |   |  |
    |  o---o  |
    +---------+

These are called **curvenet islands** in the explanation.

They represent details finer than the mesh resolution. In other words, the curvenet is trying to carve structure into a single polygon at a scale the mesh cannot meaningfully represent.

The paper detects these islands as connected components of the cut-mesh made only of **face-samples attached to the same mesh face**.

In their implementation, they simply remove these isolated components.

So if a tiny cluster of cuts lives entirely inside one face and does not meaningfully connect to the larger mesh structure, it gets discarded.

---

## Data Structure Details

The implementation uses a customized **halfedge data structure**.

Important additions include:

- each halfedge stores which **oriented curvenet segment** it corresponds to
- each cut-vertex stores a **bitmask** indicating whether it is:
  - a mesh vertex
  - a curvenet sample
  - a segment-edge intersection

The structure also allows special internal cut-edges whose opposite halfedges both point to the same cut-face, which is how the implementation represents cracks inside polygons.

---

## What the Final Cut-Mesh Looks Like

After all steps are complete, the cut-mesh contains:

- original mesh vertices
- inserted sample vertices
- inserted intersection vertices
- cut-edges following the curvenet
- new faces split along those cut-edges

So the curvenet is no longer just "drawn on top" of the mesh. It is now part of the mesh topology.

That is the whole purpose of the method.

---

## Why This Is Useful

Later parts of the algorithm want to evaluate or interpolate quantities differently on the two sides of a profile curve.

That only works robustly if the curve has become an actual topological boundary inside the mesh representation.

By cutting the mesh this way, the algorithm can:

- distinguish the two sides of a curve
- attach values to cut regions consistently
- propagate articulation information over the mesh in a structured way

---

## Intuition in One Sentence

The method takes a network of surface curves and turns it into actual mesh topology by projecting curve samples to the mesh, inserting vertices and edges where needed, and rebuilding faces so the mesh is explicitly split along the curves.

---

## Super Short Summary

You can think of the whole section as:

1. project curve points onto the mesh
2. classify them as lying on vertices, edges, or faces
3. insert them as cut-vertices
4. trace each curve segment across the mesh
5. add cut-edges along those traced paths
6. rebuild the mesh faces using halfedge connectivity
7. discard tiny isolated cut structures that are smaller than the mesh resolution

The result is a **cut-mesh** where the curvenet is embedded directly into the mesh.