#pragma once

#include "components/dsegment.hpp"
#include "components/dspline.hpp"
#include "components/dvert.hpp"
#include "curvenet/curvenet.hpp"

#include <Eigen/Core>

#include <vector>

namespace DCurvenet {

class dcurvenet {
public:
    dcurvenet() = default;
    dcurvenet(const Curvenet::curvenet& CN,
              double meanEdgeLength,
              int samplesPerMeanEdge = 5,
              int uniformRefineSamples = 64);

    const std::vector<dvert>& verts() const { return dVerts; }
    const std::vector<dsegment>& segments() const { return dSegments; }
    const std::vector<dspline>& splines() const { return dSplines; }

    double getMeanEdgeLength() const { return mean_edge_length; }
    int getSamplesPerMeanEdge() const { return samples_per_mean_edge; }
    int getUniformRefineSamples() const { return uniform_refine_samples; }

private:
    static int computeInteriorSampleCount(const Curvenet::spline& s,
                                          double meanEdgeLength,
                                          int samplesPerMeanEdge);

    static std::vector<Eigen::Vector3d> resamplePolylineEvenArcLength(const std::vector<Eigen::Vector3d>& polyline,
                                                                       int outputCount);

    bool outgoingAtControl(int controlDvertIdx,
                           int dsplineIdx,
                           int& segIdx,
                           Eigen::Vector3d& tOut,
                           double& segLen) const;
    const SegmentData* findControlSegmentData(int controlDvertIdx, int segmentIdx) const;
    void computeControlRibbonData(const std::vector<Curvenet::control>& controls,
                                  const std::vector<int>& controlToDvert);
    void interpolateRibbonAlongDSplines();

    void initializeFromCurvenet(const Curvenet::curvenet& CN);

    std::vector<dvert> dVerts;
    std::vector<dsegment> dSegments;
    std::vector<dspline> dSplines;

    double mean_edge_length = 1.0;
    int samples_per_mean_edge = 5;
    int uniform_refine_samples = 64;
};

}  // namespace DCurvenet
