#include "dvert.hpp"

#include <Eigen/Core>
#include <vector>

namespace DCurvenet {
    dvert::dvert() {}
    dvert::dvert(Eigen::Vector3d position): pos(position) {}
    dvert::dvert(Eigen::Vector3d position, int origin, int element_type, int element_index):
        pos(position), DC_origin(origin), projection_element(element_type), element_idx(element_index) {}

    int dvert::setPosition(Eigen::Vector3d position) {
        pos = position;
        return 1;
    }

    Eigen::Vector3d dvert::getPosition() {
        return pos;
    }

    // Returns the projection element and element idx
    std::pair<int, int> dvert::getProjection() {
        return std::pair<int, int>({projection_element, element_idx});
    }
    // Set the projection element
    void dvert::setProjection(int proj_element, int el_idx) {
        projection_element = proj_element;
        element_idx = el_idx;
    }
}   // namespace DCurvenet