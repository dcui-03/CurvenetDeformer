// polyHE.h
#pragma once

#include "polyHE_types.hpp" // face_t, edge_t
#include <vector>
#include <map>

namespace polyHE {

// polyHE_t::build() needs the unordered edges of the mesh.  If you don't have them, call this first.
void unordered_edges_from_faces( const int num_faces, const polyHE::face_t* faces, std::vector< polyHE::edge_t >& edges_out );

class polyHE_t
{
public:
    struct halfedge_t
    {
        // Store dest only; recover source using twin
        int dest;

        // Parent face in array.
        int face;
        // Parent edge in array.
        int edge;
        // Twin halfedge in array
        int twin;
        // Next halfedge in array.
        int next;
        
        halfedge_t() :
            dest(-1),
            face(-1),
            edge(-1),
            twin(-1),
            next(-1)
            {};
    };

    void build_from_face_list(const int num_vertices, const std::vector<std::vector<int>>& fList);
    
    // Builds the half-edge data structures from the given faces and edges.
    // NOTE: 'edges' can be derived from 'faces' by calling
    //       unordered_edges_from_faces(), above.  build() does nots
    //       but could do this for callers who do not already have edges.
    // NOTE: 'faces' and 'edges' are not needed after the call to build()
    //       completes and may be destroyed.
    void build( const int num_vertices, const int num_faces, const polyHE::face_t* faces, const int num_edges, const polyHE::edge_t* edges );
    
    void clear()
    {
        m_halfedges.clear();
        m_vertex_halfedges.clear();
        m_face_halfedges.clear();
        m_edge_halfedges.clear();
        vertPairToHE.clear();
    }
    
    const halfedge_t& halfedge( const int i ) const { return m_halfedges.at( i ); }
    
    std::pair< int, int > he_index2directed_edge( const int he_index ) const
    {
        /*
        Given the index of a halfedge_t, returns the corresponding directed edge (i,j).
        
        untested
        */
        
        const halfedge_t& he = m_halfedges[ he_index ];
        return std::make_pair( m_halfedges[ he.twin ].dest, he.dest );
    }
    
    int directed_edge2he_index( const int i, const int j ) const
    {
        /*
        Given a directed edge (i,j), returns the index of the 'halfedge_t' in
        halfedges().
        
        untested
        */
        
        /// This isn't const, and doesn't handle the case where (i,j) isn't known:
        // return vertPairToHE[ std::make_pair( i,j ) ];
        
        pairToIdx::const_iterator result = vertPairToHE.find( std::make_pair( i,j ) );
        if( result == vertPairToHE.end() ) return -1;
        
        return result->second;
    }
    
    void vertex_vertex_neighbors( const int vertex_index, std::vector< int >& result ) const {
        /*
        Returns in 'result' the vertex neighbors (as indices) of the vertex 'vertex_index'.
        
        untested
        */
        
        result.clear();
        
        const int start_hei = m_vertex_halfedges[ vertex_index ];
        int hei = start_hei;
        while( true )
        {
            const halfedge_t& he = m_halfedges[ hei ];
            result.push_back( he.dest );
            
            hei = m_halfedges[ he.twin ].next;
            if( hei == start_hei ) break;
        }
    }
    std::vector< int > vertex_vertex_neighbors( const int vertex_index ) const {
        std::vector< int > result;
        vertex_vertex_neighbors( vertex_index, result );
        return result;
    }
    
    int vertex_valence( const int vertex_index ) const {
        /*
        Returns the valence (number of vertex neighbors) of vertex with index 'vertex_index'.
        
        untested
        */
        
        std::vector< int > neighbors;
        vertex_vertex_neighbors( vertex_index, neighbors );
        return neighbors.size();
    }
    
    void vertex_face_neighbors( const int vertex_index, std::vector< int >& result ) const {
        /*
        Returns in 'result' the face neighbors (as indices) of the vertex 'vertex_index'.
        
        untested
        */
        
        result.clear();
        
        const int start_hei = m_vertex_halfedges[ vertex_index ];
        int hei = start_hei;
        while( true )
        {
            const halfedge_t& he = m_halfedges[ hei ];
            if( -1 != he.face ) result.push_back( he.face );
            
            hei = m_halfedges[ he.twin ].next;
            if( hei == start_hei ) break;
        }
    }

    std::vector< int > vertex_face_neighbors( const int vertex_index ) const {
        std::vector< int > result;
        vertex_face_neighbors( vertex_index, result );
        return result;
    }
    
    bool vertex_is_boundary( const int vertex_index ) const
    {
        /*
        Returns whether the vertex with given index is on the boundary.
        
        untested
        */
        
        return -1 == m_halfedges[ m_vertex_halfedges[ vertex_index ] ].face;
    }
    

    // GETTERS
    // Get number of edges
    int numEdges();
    // Get one halfedge of the edge
    int edgeHalfedge(int ei);
    // Get the pair of vertices for the edge
    std::pair<int,int> edgeVertices(int ei);
    // Get the pair of faces associated with the edge
    std::pair<int,int> edgeFaces(int ei);

    // Get edge index from a pair of vertices
    int edgeIdxFromVerts(int vi, int vj);



    std::vector<int> boundary_vertices() const;
    
    std::vector<std::pair<int, int>> boundary_edges() const;
    
private:
    std::vector< halfedge_t > m_halfedges;
    // Offsets into the 'halfedges' sequence, one per vertex.
    std::vector< int > m_vertex_halfedges;
    // Offset into the 'halfedges' sequence, one per face.
    std::vector< int > m_face_halfedges;
    // Offset into the 'halfedges' sequence, one per edge (unordered pair of vertex indices).
    std::vector< int > m_edge_halfedges;

    // A map from an ordered edge (an std::pair of int's) to an offset into the 'halfedge' sequence.
    typedef std::map< std::pair< int, int >, int > pairToIdx;
    pairToIdx vertPairToHE;
};

} // namespace polyHE
