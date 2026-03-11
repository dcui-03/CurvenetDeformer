// jsonUtils.cpp

#include "jsonUtils.hpp"
#include "utils.hpp"

#include <nlohmann/json.hpp>

#include <cmath>
#include <cstdint>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>

namespace JSONUtils {

using json = nlohmann::json;

namespace {

double asNumber(const json& v) {
    if (v.is_number_float()) return v.get<double>();
    if (v.is_number_integer()) return static_cast<double>(v.get<std::int64_t>());
    if (v.is_number_unsigned()) return static_cast<double>(v.get<std::uint64_t>());
    throw std::runtime_error("Expected numeric JSON value");
}

Eigen::Vector3d parseVec3(const json& v) {
    if (!v.is_array() || v.size() != 3) {
        throw std::runtime_error("Expected 3-element JSON array");
    }
    return {asNumber(v[0]), asNumber(v[1]), asNumber(v[2])};
}

bool nearlyEqual(const Eigen::Vector3d& a, const Eigen::Vector3d& b, double eps) {
    return (a - b).cwiseAbs().maxCoeff() <= eps;
}

struct QuantizedKey {
    std::int64_t x = 0;
    std::int64_t y = 0;
    std::int64_t z = 0;

    bool operator==(const QuantizedKey& rhs) const { return x == rhs.x && y == rhs.y && z == rhs.z; }
};

struct QuantizedKeyHash {
    std::size_t operator()(const QuantizedKey& k) const {
        std::size_t h1 = std::hash<std::int64_t>{}(k.x);
        std::size_t h2 = std::hash<std::int64_t>{}(k.y);
        std::size_t h3 = std::hash<std::int64_t>{}(k.z);
        return h1 ^ (h2 << 1U) ^ (h3 << 2U);
    }
};

QuantizedKey toKey(const Eigen::Vector3d& p, double eps) {
    return {
        static_cast<std::int64_t>(std::llround(p.x() / eps)),
        static_cast<std::int64_t>(std::llround(p.y() / eps)),
        static_cast<std::int64_t>(std::llround(p.z() / eps)),
    };
}

}  // namespace

CurvenetInput loadBezierCurvenetInput(const std::string& jsonPath, double dedupTolerance) {
    if (dedupTolerance <= 0.0) {
        throw std::invalid_argument("dedupTolerance must be > 0");
    }

    std::ifstream in(jsonPath);
    if (!in) {
        throw std::runtime_error("Failed to open curve JSON: " + jsonPath);
    }

    std::string text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    json root = json::parse(text);

    CurvenetInput out;
    std::unordered_map<QuantizedKey, int, QuantizedKeyHash> endpointMap;
    endpointMap.reserve(1024);

    auto addPoint = [&](const Eigen::Vector3d& p, const Eigen::Vector3d& n) -> int {
        out.controlP.push_back(p);
        out.surfaceN.push_back(n.normalized());
        return static_cast<int>(out.controlP.size() - 1);
    };

    auto addEndpoint = [&](const Eigen::Vector3d& p, const Eigen::Vector3d& n) -> int {
        const QuantizedKey key = toKey(p, dedupTolerance);
        auto it = endpointMap.find(key);
        if (it != endpointMap.end()) {
            const int idx = it->second;
            if (!nearlyEqual(out.controlP[static_cast<std::size_t>(idx)], p, dedupTolerance)) {
                throw std::runtime_error("Endpoint quantization collision while loading JSON");
            }
            return idx;
        }

        const int idx = addPoint(p, n);
        endpointMap.emplace(key, idx);
        return idx;
    };

    for (const json& curveObj : root.at("curves")) {
        for (const json& splineObj : curveObj.at("splines")) {
            if (splineObj.at("type").get<std::string>() != "BEZIER") {
                continue;
            }

            const bool cyclic = splineObj.at("is_cyclic").get<bool>();
            const json& cps = splineObj.at("control_points");
            if (!cps.is_array() || cps.size() < 2) {
                continue;
            }

            const std::size_t segCount = cyclic ? cps.size() : (cps.size() - 1);
            for (std::size_t i = 0; i < segCount; ++i) {
                const std::size_t j = (i + 1) % cps.size();

                const json& a = cps[i];
                const json& b = cps[j];

                const Eigen::Vector3d p0 = parseVec3(a.at("co"));
                const Eigen::Vector3d p1 = parseVec3(b.at("co"));
                const Eigen::Vector3d t0 = parseVec3(a.at("handle_right"));
                const Eigen::Vector3d t1 = parseVec3(b.at("handle_left"));
                const Eigen::Vector3d n0 = parseVec3(a.at("mesh_binding").at("normal_world"));
                const Eigen::Vector3d n1 = parseVec3(b.at("mesh_binding").at("normal_world"));

                const int c0 = addEndpoint(p0, n0);
                const int c1 = addEndpoint(p1, n1);
                const int h0 = addPoint(t0, n0);
                const int h1 = addPoint(t1, n1);

                out.curveC.push_back({c0, h0, h1, c1});
            }
        }
    }

    return out;
}

}  // namespace JSONUtils
