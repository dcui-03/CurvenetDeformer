// polyHE_types.h
#pragma once

#include <vector>
#include <array>

namespace polyHE {
    
    // Edge type
    struct edge_t
    {   
        // Store both vertices, halfedges, plus the two adjacent faces
        std::array<int, 2> v;
        std::array<int, 2> f;
        
        int& v0() { return v[0]; }
        const int v0() const { return v[0]; }
        
        int& v1() { return v[1]; }
        const int v1() const { return v[1]; }

        int& f0() { return f[0]; }
        const int f0() const { return f[0]; }

        int& f1() { return f[1]; }
        const int& f1() const { return f[1]; }
        
        edge_t()
        {
            v[0] = v[1] = -1;
            f[0] = f[1] = -1;
        }
    };
    
    struct face_t {
        std::vector<int> v;
        
        int vi(int i) {
            if (i < 0 || i >= static_cast<int>(v.size())) {
                return -1;
            }
            return v[i];
        }
        const int vi(int i) const {
            if (i < 0 || i >= static_cast<int>(v.size())) {
                return -1;
            }
            return v[i];
        }
        
        face_t(std::vector<int> verts) {
            for (int i = 0; i < verts.size(); i++) {
                v.push_back(verts[i]);
            }
        }
    };
}
