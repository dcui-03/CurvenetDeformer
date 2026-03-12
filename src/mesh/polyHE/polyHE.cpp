#include "polyHE.hpp"

// needed for implementation
#include <cassert>
#include <set>
#include <iostream>

// TODO: What is this for?
namespace {
typedef std::map<std::pair< int, int >, int> pairToIdx;

int directed_edge2face_index( const pairToIdx& de2fi, int vertex_i, int vertex_j ) {
    assert( !de2fi.empty() );
    
    pairToIdx::const_iterator it = de2fi.find( std::make_pair( vertex_i, vertex_j ) );
    
    // If no such directed edge exists, then there's no such face in the mesh.
    // The edge must be a boundary edge.
    // In this case, the reverse orientation edge must have a face.
    if( it == de2fi.end() ) {
        assert( de2fi.find( std::make_pair( vertex_j, vertex_i ) ) != de2fi.end() );
        return -1;
    }
    
    return it->second;
}
}



namespace polyHE {

// build a polyHE_t object using the number of vertices and a face list
void polyHE_t::build_from_face_list(
    const int num_vertices,
    const std::vector<std::vector<int>>& fList
) {
    std::vector<face_t> faces;
    faces.reserve(fList.size());

    for (std::size_t fi = 0; fi < fList.size(); ++fi)
    {
        const std::vector<int>& verts = fList[fi];

        // Basic validation
        assert(verts.size() >= 3);

        // Optional: check all indices are in range
        for (std::size_t k = 0; k < verts.size(); ++k)
        {
            assert(verts[k] >= 0);
            assert(static_cast<int>(verts[k]) < num_vertices);
        }

        faces.emplace_back(verts);
    }

    std::vector<edge_t> edges;
    unordered_edges_from_faces(
        static_cast<int>(faces.size()),
        faces.data(),
        edges
    );

    build(
        num_vertices,
        static_cast<int>(faces.size()),
        faces.data(),
        static_cast<int>(edges.size()),
        edges.data()
    );
}

// build a polyHE_t object using internal types
void polyHE_t::build( const int num_vertices, const int num_faces, const face_t* faces, const int num_edges, const edge_t* edges ) {
    /*
    Generates all half edge data structures for the mesh given by its vertices 'self.vs'
    and faces 'self.faces'.
    
    Python version used heavily
    */
    
    assert( faces );
    assert( edges );
    
    pairToIdx de2fi;
    for( int fi = 0; fi < num_faces; ++fi ) {
        const face_t& face = faces[fi];
        const int n = static_cast<int>(face.v.size());
        assert(n >= 3);

        for (int r = 0; r < n; ++r)
        {
            const int a = face.v[r];
            const int b = face.v[(r + 1) % n];

            // For a manifold oriented mesh, each directed edge should belong
            // to at most one face.
            assert(de2fi.find(std::make_pair(a, b)) == de2fi.end());
            de2fi[std::make_pair(a, b)] = fi;
        }
    }
    
    clear();
    m_vertex_halfedges.resize( num_vertices, -1 );
    m_face_halfedges.resize( num_faces, -1 );
    m_edge_halfedges.resize( num_edges, -1 );
    m_halfedges.reserve( num_edges*2 );
    
    for( int ei = 0; ei < num_edges; ++ei )
    {
        const edge_t& edge = edges[ei];
        
        // Add the halfedge_t structures to the end of the list.
        const int he0index = m_halfedges.size();
        m_halfedges.push_back( halfedge_t() );
        halfedge_t& he0 = m_halfedges.back();
        
        const int he1index = m_halfedges.size();
        m_halfedges.push_back( halfedge_t() );
        halfedge_t& he1 = m_halfedges.back();
        
        // The face will be -1 if it is a boundary half-edge.
        he0.face = directed_edge2face_index( de2fi, edge.v[0], edge.v[1] );
        he0.dest = edge.v[1];
        he0.edge = ei;
        
        // The face will be -1 if it is a boundary half-edge.
        he1.face = directed_edge2face_index( de2fi, edge.v[1], edge.v[0] );
        he1.dest = edge.v[0];
        he1.edge = ei;
        
        // Store the opposite half-edge index.
        he0.twin = he1index;
        he1.twin = he0index;
        
        // Also store the index in our vertPairToHE map.
        assert( vertPairToHE.find( std::make_pair( edge.v[0], edge.v[1] ) ) == vertPairToHE.end() );
        assert( vertPairToHE.find( std::make_pair( edge.v[1], edge.v[0] ) ) == vertPairToHE.end() );
        vertPairToHE[ std::make_pair( edge.v[0], edge.v[1] ) ] = he0index;
        vertPairToHE[ std::make_pair( edge.v[1], edge.v[0] ) ] = he1index;
        
        // If the vertex pointed to by a half-edge doesn't yet have an out-going
        // halfedge, store the opposite halfedge.
        // Also, if the vertex is a boundary vertex, make sure its
        // out-going halfedge is a boundary halfedge.
        // NOTE: Halfedge data structure can't properly handle butterfly vertices.
        //       If the mesh has butterfly vertices, there will be multiple outgoing
        //       boundary halfedges.  Because we have to pick one as the vertex's outgoing
        //       halfedge, we can't iterate over all neighbors, only a single wing of the
        //       butterfly.
        if( m_vertex_halfedges[ he0.dest ] == -1 || -1 == he1.face )
        {
            m_vertex_halfedges[ he0.dest ] = he0.twin;
        }
        if( m_vertex_halfedges[ he1.dest ] == -1 || -1 == he0.face )
        {
            m_vertex_halfedges[ he1.dest ] = he1.twin;
        }
        
        // If the face pointed to by a half-edge doesn't yet have a
        // halfedge pointing to it, store the halfedge.
        if( -1 != he0.face && m_face_halfedges[ he0.face ] == -1 )
        {
            m_face_halfedges[ he0.face ] = he0index;
        }
        if( -1 != he1.face && m_face_halfedges[ he1.face ] == -1 )
        {
            m_face_halfedges[ he1.face ] = he1index;
        }
        
        // Store one of the half-edges for the edge.
        assert( m_edge_halfedges[ ei ] == -1 );
        m_edge_halfedges[ ei ] = he0index;
    }
    
    // Now that all the half-edges are created, set the remaining next field.
    // We can't yet handle boundary halfedges, so store them for later.
    std::vector< int > boundary_heis;
    for( int hei = 0; hei < m_halfedges.size(); ++hei )
    {
        halfedge_t& he = m_halfedges.at( hei );
        // Store boundary halfedges for later.
        if( -1 == he.face )
        {
            boundary_heis.push_back( hei );
            continue;
        }
        
        const face_t& face = faces[ he.face ];
        const int n = static_cast<int>(face.v.size());
        const int i = he.dest;

        int next_pos = -1;
        for (int r = 0; r < n; ++r)
        {
            if (face.v[r] == i)
            {
                next_pos = (r + 1) % n;
                break;
            }
        }
        assert(next_pos != -1);

        const int j = face.v[next_pos];
        he.next = vertPairToHE[ std::make_pair(i, j) ];
    }
    
    // Make a map from vertices to boundary halfedges (indices) originating from them.
    // NOTE: There will only be multiple originating boundary halfedges at butterfly vertices.
    std::map< int, std::set< int > > vertex2outgoing_boundary_hei;
    for( std::vector< int >::const_iterator hei = boundary_heis.begin(); hei != boundary_heis.end(); ++hei )
    {
        const int originating_vertex = m_halfedges[ m_halfedges[ *hei ].twin ].dest;
        vertex2outgoing_boundary_hei[ originating_vertex ].insert( *hei );
        if( vertex2outgoing_boundary_hei[ originating_vertex ].size() > 1 )
        {
            std::cerr << "Butterfly vertex encountered.\n";
        }
    }
    
    // For each boundary halfedge, make its next one of the boundary halfedges
    // originating at its dest.
    for( std::vector< int >::const_iterator hei = boundary_heis.begin(); hei != boundary_heis.end(); ++hei )
    {
        halfedge_t& he = m_halfedges[ *hei ];
        
        std::set< int >& outgoing = vertex2outgoing_boundary_hei[ he.dest ];
        if( !outgoing.empty() )
        {
            std::set< int >::iterator outgoing_hei = outgoing.begin();
            he.next = *outgoing_hei;
            
            outgoing.erase( outgoing_hei );
        }
    }
    
#ifndef NDEBUG
    for( std::map< int, std::set< int > >::const_iterator it = vertex2outgoing_boundary_hei.begin(); it != vertex2outgoing_boundary_hei.end(); ++it )
    {
        assert( it->second.empty() );
    }
#endif
}


// GETTERS
// Get number of edges
int polyHE_t::numEdges() {
    return static_cast<int>(m_edge_halfedges.size());
}

// Get one halfedge of the edge
int polyHE_t::edgeHalfedge(int ei) {
    return m_edge_halfedges.at(ei);
}

// Get the pair of vertices for the edge
std::pair<int,int> polyHE_t::edgeVertices(int ei) {
    int hei = m_edge_halfedges.at(ei);
    const halfedge_t& he = m_halfedges.at(hei);
    const halfedge_t& ht = m_halfedges.at(he.twin);
    return std::make_pair(ht.dest, he.dest);
}

// Get the pair of faces associated with the edge
std::pair<int,int> polyHE_t::edgeFaces(int ei) {
    int hei = m_edge_halfedges.at(ei);
    const halfedge_t& he = m_halfedges.at(hei);
    const halfedge_t& ht = m_halfedges.at(he.twin);
    return std::make_pair(he.face, ht.face);
}


// Get edge index from a vertex pair
int polyHE_t::edgeIdxFromVerts(int vi, int vj) {
    pairToIdx::const_iterator it = vertPairToHE.find(std::make_pair(vi, vj));

    if (it == vertPairToHE.end()) {
        return -1;
    }

    int hei = it->second;
    return m_halfedges.at(hei).edge;
}

// Finds all boundary vertices (by index)
std::vector<int> polyHE_t::boundary_vertices() const {
    /*
    Returns a list of the vertex indices on the boundary.
    
    untested
    */
    
    std::set<int> result;
    for( int hei = 0; hei < m_halfedges.size(); ++hei )
    {
        const halfedge_t& he = m_halfedges[ hei ];
        
        if(-1 == he.face)
        {
            // result.extend( self.he_index2directed_edge( hei ) )
            result.insert(he.dest);
            result.insert(m_halfedges[he.twin].dest);
        }
    }
    
    return std::vector<int>( result.begin(), result.end() );
}

// Returns the boundary edges as pairs of vertices
std::vector<std::pair< int, int>> polyHE_t::boundary_edges() const {
    /*
    Returns a list of undirected boundary edges (i,j).  If (i,j) is in the result, (j,i) will not be.
    
    untested
    */
    
    std::vector<std::pair<int, int>> result;
    for( int hei = 0; hei < m_halfedges.size(); ++hei )
    {
        const halfedge_t& he = m_halfedges[ hei ];
        
        if( -1 == he.face )
        {
            result.push_back( he_index2directed_edge( hei ) );
        }
    }
    
    return result;
}

// Returns ???? TODO
void unordered_edges_from_faces( const int num_faces, const face_t* faces, std::vector<edge_t>& edges_out ) {
    typedef std::set<std::pair<int, int>> edge_set_t;
    edge_set_t edges;
    for( int t = 0; t < num_faces; ++t ) {   
        int n = faces[t].v.size();
        for (int r = 0; r < n; r++) {
            edges.insert( std::make_pair( std::min( faces[t].vi(r), faces[t].vi((r+1)%n) ), std::max( faces[t].vi(r), faces[t].vi((r+1)%n) ) ) );
        }
    }
    
    edges_out.resize(edges.size());
    int e = 0;
    for( edge_set_t::const_iterator it = edges.begin(); it != edges.end(); ++it, ++e ) {
        edges_out[e].v[0] = it->first;
        edges_out[e].v[1] = it->second;
    }
}

}   // namespace polyHE
