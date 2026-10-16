// SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
// SPDX-License-Identifier: Apache-2.0

// USD-free, intentionally small Python binding for the local modeler.
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <SmApiBrep.h>
#include <SmApiCurves.h>
#include <SmApiFillets.h>
#include <SmApiGeneral.h>
#include <SmApiPrimitives.h>
#include <SmApiQueries.h>
#include <SmBrep.h>
#include <SmFace.h>
#include <SmLine.h>
#include <SmPoly.h>
#include <SmVector3d.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace py = pybind11;

namespace {

void check(SmApiStatus status, const char* operation)
{
    if (status != SM_SUCCESS)
        throw std::runtime_error(std::string(operation) + " returned status " + std::to_string(status));
}

SmVector3d vector3(const py::tuple& value)
{
    if (value.size() != 3)
        throw py::value_error("Expected a 3-tuple");
    const double x = value[0].cast<double>();
    const double y = value[1].cast<double>();
    const double z = value[2].cast<double>();
    if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z))
        throw py::value_error("Coordinates must be finite");
    return SmVector3d(x, y, z);
}

void positive(double value, const char* name)
{
    if (!std::isfinite(value) || value <= 0.0)
        throw py::value_error(std::string(name) + " must be finite and positive");
}

void nonnegative(double value, const char* name)
{
    if (!std::isfinite(value) || value < 0.0)
        throw py::value_error(std::string(name) + " must be finite and non-negative");
}

void angular_sweep(double start, double end)
{
    if (!std::isfinite(start) || !std::isfinite(end) || end <= start || end - start > 360.0)
        throw py::value_error("end_angle_deg must be greater than start_angle_deg by at most 360 degrees");
}

struct BrepHandle
{
    explicit BrepHandle(SmBrep* raw) : value(raw)
    {
        if (!raw)
            throw std::runtime_error("SM_API returned a null Brep");
    }

    SmBrep* get() const
    {
        if (!value)
            throw py::value_error("Brep was consumed by a Boolean operation");
        return value.get();
    }

    std::unique_ptr<SmBrep> value;
};

struct PolyHandle
{
    explicit PolyHandle(SmPolyBrep* raw) : value(raw)
    {
        if (!raw)
            throw std::runtime_error("Tessellation returned a null PolyBrep");
    }

    std::unique_ptr<SmPolyBrep> value;
};

using Brep = std::shared_ptr<BrepHandle>;
using PolyBrep = std::shared_ptr<PolyHandle>;

struct Profile
{
    std::vector<std::array<double, 2>> points;
};

Profile make_profile(const std::vector<std::array<double, 2>>& input)
{
    Profile profile{input};
    auto& points = profile.points;
    if (points.size() > 1 && points.front() == points.back())
        points.pop_back();
    if (points.size() < 3 || points.size() > 1024)
        throw py::value_error("profile requires 3 to 1024 vertices");
    double scale = 0.0;
    for (const auto& p : points)
        for (int j = 0; j < 2; ++j)
        {
            if (!std::isfinite(p[j]))
                throw py::value_error("profile coordinates must be finite");
            scale = std::max(scale, std::abs(p[j] - points.front()[j]));
        }
    if (!std::isfinite(scale) || scale == 0.0)
        throw py::value_error("profile must have nonzero area");
    auto normalized = points;
    for (auto& p : normalized)
        for (int j = 0; j < 2; ++j)
            p[j] = (p[j] - points.front()[j]) / scale;
    auto cross = [&](size_t a, size_t b, size_t c) {
        return (normalized[b][0] - normalized[a][0]) * (normalized[c][1] - normalized[a][1])
             - (normalized[b][1] - normalized[a][1]) * (normalized[c][0] - normalized[a][0]);
    };
    constexpr double eps = 1e-14;
    const size_t n = points.size();
    double area = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        const size_t next = (i + 1) % n;
        if (std::hypot(normalized[i][0] - normalized[next][0],
                       normalized[i][1] - normalized[next][1]) <= eps)
            throw py::value_error("profile contains a zero-length edge");
        area += normalized[i][0] * normalized[next][1] - normalized[next][0] * normalized[i][1];
        const size_t prev = (i + n - 1) % n;
        if (std::abs(cross(prev, i, next)) <= eps
            && (normalized[prev][0] - normalized[i][0]) * (normalized[next][0] - normalized[i][0])
             + (normalized[prev][1] - normalized[i][1]) * (normalized[next][1] - normalized[i][1]) > eps)
            throw py::value_error("profile contains overlapping adjacent edges");
    }
    auto on_segment = [&](size_t a, size_t b, size_t p) {
        return std::abs(cross(a, b, p)) <= eps
            && normalized[p][0] >= std::min(normalized[a][0], normalized[b][0]) - eps
            && normalized[p][0] <= std::max(normalized[a][0], normalized[b][0]) + eps
            && normalized[p][1] >= std::min(normalized[a][1], normalized[b][1]) - eps
            && normalized[p][1] <= std::max(normalized[a][1], normalized[b][1]) + eps;
    };
    for (size_t a = 0; a < n; ++a)
        for (size_t c = a + 1; c < n; ++c)
        {
            const size_t b = (a + 1) % n, d = (c + 1) % n;
            if (b == c || d == a)
                continue;
            const double ab_c = cross(a, b, c), ab_d = cross(a, b, d);
            const double cd_a = cross(c, d, a), cd_b = cross(c, d, b);
            if (((ab_c > eps && ab_d < -eps) || (ab_c < -eps && ab_d > eps))
                && ((cd_a > eps && cd_b < -eps) || (cd_a < -eps && cd_b > eps)))
                throw py::value_error("profile must not self-intersect");
            if (on_segment(a, b, c) || on_segment(a, b, d) || on_segment(c, d, a) || on_segment(c, d, b))
                throw py::value_error("profile must not self-touch or overlap");
        }
    if (std::abs(area) <= eps)
        throw py::value_error("profile must have nonzero area");
    if (area < 0.0)
        std::reverse(points.begin(), points.end());
    return profile;
}

SmVector3d unit_vector(const SmVector3d& value, const char* name)
{
    const double scale = std::max({std::abs(value.x), std::abs(value.y), std::abs(value.z)});
    if (scale == 0.0)
        throw py::value_error(std::string(name) + " must be nonzero");
    SmVector3d result(value.x / scale, value.y / scale, value.z / scale);
    result *= 1.0 / result.Length();
    return result;
}

struct ProfileFrame
{
    SmVector3d origin, x, y, z;

    SmPoint3d point(double u, double v, double w) const
    {
        const SmVector3d p = origin + x * u + y * v + z * w;
        if (!std::isfinite(p.x) || !std::isfinite(p.y) || !std::isfinite(p.z))
            throw py::value_error("placed profile coordinates must be finite");
        return SmPoint3d(p.x, p.y, p.z);
    }
};

ProfileFrame profile_frame(py::tuple origin, const py::object& axis, const py::object& xdir)
{
    ProfileFrame f;
    f.origin = vector3(origin);
    f.z = axis.is_none() ? SmVector3d(0, 0, 1) : unit_vector(vector3(axis.cast<py::tuple>()), "axis");
    if (!xdir.is_none())
    {
        const SmVector3d direction = unit_vector(vector3(xdir.cast<py::tuple>()), "xdir");
        f.x = direction - f.z * direction.Dot(f.z);
        if (f.x.Length() <= 1e-12)
            throw py::value_error("xdir must not be parallel to axis");
        f.x = unit_vector(f.x, "xdir");
    }
    else
    {
        // Match orient_axis: shortest +Z rotation; -Z rotates about +X.
        const double lateral = std::hypot(f.z.x, f.z.y);
        f.x = SmVector3d(1, 0, 0);
        if (lateral > 1e-12)
        {
            const double ux = -f.z.y / lateral, uy = f.z.x / lateral;
            f.x = SmVector3d(f.z.z + ux * ux * (1 - f.z.z), ux * uy * (1 - f.z.z), -uy * lateral);
        }
    }
    f.y = SmVector3d(f.z.y * f.x.z - f.z.z * f.x.y,
                     f.z.z * f.x.x - f.z.x * f.x.z,
                     f.z.x * f.x.y - f.z.y * f.x.x);
    return f;
}

Brep make_profile_solid(const Profile& profile, py::tuple origin, double start, double end,
                        const py::object& axis, const py::object& xdir, bool revolve)
{
    if (revolve)
        angular_sweep(start, end);
    else if (!std::isfinite(start) || !std::isfinite(end) || !std::isfinite(end - start) || end <= start)
        throw py::value_error("extrude requires finite start < end");
    const ProfileFrame f = profile_frame(origin, axis, xdir);
    std::vector<SmPoint3d> points;
    const double angle = std::remainder(start, 360.0) * (std::acos(-1.0) / 180.0);
    for (const auto& p : profile.points)
    {
        if (revolve && p[0] < 0.0)
            throw py::value_error("revolve profile radius must be non-negative");
        points.push_back(revolve ? f.point(p[0] * std::cos(angle), p[0] * std::sin(angle), p[1])
                                 : f.point(p[0], p[1], start));
    }
    std::vector<std::unique_ptr<SmLine>> owned;
    SmTArray<SmCurve*> curves;
    for (size_t i = 0; i < points.size(); ++i)
    {
        SmLine* line = nullptr;
        check(SmApiCreateLineSegment(points[i], points[(i + 1) % points.size()], line), "SmApiCreateLineSegment");
        owned.emplace_back(line);
        curves.Add(line);
    }
    // SM_API copies the curves; owned releases our originals on every exit.
    SmBrep* result = nullptr;
    const SmApiStatus status = revolve
        ? SmApiCreateRotationalSweep(curves, f.origin, f.z, end - start, TRUE, result)
        : SmApiCreateLinearSweep(curves, f.z, end - start, TRUE, result);
    std::unique_ptr<SmBrep> cleanup(result);
    check(status, revolve ? "SmApiCreateRotationalSweep" : "SmApiCreateLinearSweep");
    if (!result)
        throw std::runtime_error("profile sweep returned a null Brep");
    SmBoolean manifold = FALSE;
    check(SmApiBrepIsManifoldSolid(result, manifold), "SmApiBrepIsManifoldSolid");
    if (!manifold)
        throw std::runtime_error("profile sweep did not produce a manifold solid");
    return std::make_shared<BrepHandle>(cleanup.release());
}

Brep orient_axis(SmBrep* raw, const SmVector3d& pivot, const py::object& axis)
{
    std::unique_ptr<SmBrep> result(raw);
    if (!axis.is_none())
    {
        const SmVector3d direction = vector3(axis.cast<py::tuple>());
        const double scale = std::max({std::abs(direction.x), std::abs(direction.y), std::abs(direction.z)});
        if (scale == 0.0)
            throw py::value_error("axis must be nonzero");
        const double x = direction.x / scale;
        const double y = direction.y / scale;
        const double z = direction.z / scale;
        const double lateral = std::hypot(x, y);
        const double length = std::hypot(lateral, z);
        // The shortest rotation from local +Z fixes a reproducible local X/Y
        // orientation for sector shapes. -Z uses +X as the rotation axis.
        if (lateral > 1e-12 * length)
        {
            const SmVector3d rotation_axis(-y / lateral, x / lateral, 0.0);
            const double angle = std::atan2(lateral, z) * (180.0 / std::acos(-1.0));
            check(SmApiRotate(result.get(), pivot, rotation_axis, angle), "SmApiRotate");
        }
        else if (z < 0.0)
            check(SmApiRotate(result.get(), pivot, SmVector3d(1.0, 0.0, 0.0), 180.0), "SmApiRotate");
    }
    return std::make_shared<BrepHandle>(result.release());
}

Brep make_box(py::tuple origin, double length, double width, double height)
{
    positive(length, "length");
    positive(width, "width");
    positive(height, "height");
    SmBrep* result = nullptr;
    check(SmApiCreateBox(vector3(origin), length, width, height, result), "SmApiCreateBox");
    return std::make_shared<BrepHandle>(result);
}

Brep make_sphere(py::tuple center, double radius)
{
    positive(radius, "radius");
    SmBrep* result = nullptr;
    check(SmApiCreateSphere(vector3(center), radius, result), "SmApiCreateSphere");
    return std::make_shared<BrepHandle>(result);
}

Brep make_cylinder(py::tuple base, double radius, double height, py::object axis)
{
    positive(radius, "radius");
    positive(height, "height");
    const SmVector3d pivot = vector3(base);
    SmBrep* result = nullptr;
    check(SmApiCreateCylinder(pivot, radius, height, result), "SmApiCreateCylinder");
    return orient_axis(result, pivot, axis);
}

Brep make_cone(py::tuple base, double radius_base, double radius_top, double height, py::object axis)
{
    nonnegative(radius_base, "radius_base");
    nonnegative(radius_top, "radius_top");
    if (radius_base == 0.0 && radius_top == 0.0)
        throw py::value_error("at least one cone radius must be positive");
    positive(height, "height");
    const SmVector3d pivot = vector3(base);
    SmBrep* result = nullptr;
    check(SmApiCreateCone(pivot, radius_base, radius_top, height, result), "SmApiCreateCone");
    return orient_axis(result, pivot, axis);
}

Brep make_torus(py::tuple center, double radius_major, double radius_minor, py::object axis)
{
    positive(radius_major, "radius_major");
    positive(radius_minor, "radius_minor");
    if (radius_major <= radius_minor)
        throw py::value_error("radius_major must be greater than radius_minor");
    const SmVector3d pivot = vector3(center);
    SmBrep* result = nullptr;
    check(SmApiCreateTorus(pivot, radius_major, radius_minor, result), "SmApiCreateTorus");
    return orient_axis(result, pivot, axis);
}

Brep make_pyramid(py::tuple base, double length, double height, py::object axis)
{
    positive(length, "length");
    positive(height, "height");
    const SmVector3d pivot = vector3(base);
    SmBrep* result = nullptr;
    check(SmApiCreatePyramid(pivot, length, height, result), "SmApiCreatePyramid");
    return orient_axis(result, pivot, axis);
}

Brep make_cylindrical_box(py::tuple base, double length, double inside_radius,
                          double outside_radius, double start_angle, double end_angle, py::object axis)
{
    positive(length, "length");
    positive(inside_radius, "inside_radius");
    positive(outside_radius, "outside_radius");
    if (outside_radius <= inside_radius)
        throw py::value_error("outside_radius must be greater than inside_radius");
    angular_sweep(start_angle, end_angle);
    const SmVector3d pivot = vector3(base);
    SmBrep* result = nullptr;
    check(SmApiCreateCylindricalBox(pivot, length, inside_radius, outside_radius,
                                    start_angle, end_angle, result), "SmApiCreateCylindricalBox");
    return orient_axis(result, pivot, axis);
}

Brep make_partial_sphere(py::tuple center, double radius, double start_angle, double end_angle, py::object axis)
{
    positive(radius, "radius");
    angular_sweep(start_angle, end_angle);
    const SmVector3d pivot = vector3(center);
    SmBrep* result = nullptr;
    check(SmApiCreatePartialSphere(pivot, radius, start_angle, end_angle, result),
          "SmApiCreatePartialSphere");
    return orient_axis(result, pivot, axis);
}

Brep make_partial_cone(py::tuple base, double radius_base, double radius_top, double height,
                       double start_angle, double end_angle, py::object axis)
{
    nonnegative(radius_base, "radius_base");
    nonnegative(radius_top, "radius_top");
    if (radius_base == 0.0 && radius_top == 0.0)
        throw py::value_error("at least one cone radius must be positive");
    positive(height, "height");
    angular_sweep(start_angle, end_angle);
    const SmVector3d pivot = vector3(base);
    SmBrep* result = nullptr;
    check(SmApiCreatePartialCone(pivot, radius_base, radius_top, height,
                                 start_angle, end_angle, result), "SmApiCreatePartialCone");
    return orient_axis(result, pivot, axis);
}

Brep make_partial_cylinder(py::tuple base, double radius, double height,
                           double start_angle, double end_angle, py::object axis)
{
    positive(radius, "radius");
    positive(height, "height");
    angular_sweep(start_angle, end_angle);
    const SmVector3d pivot = vector3(base);
    SmBrep* result = nullptr;
    check(SmApiCreatePartialCylinder(pivot, radius, height, start_angle, end_angle, result),
          "SmApiCreatePartialCylinder");
    return orient_axis(result, pivot, axis);
}

Brep make_partial_torus(py::tuple center, double radius_major, double radius_minor,
                        double start_angle, double end_angle, py::object axis)
{
    positive(radius_major, "radius_major");
    positive(radius_minor, "radius_minor");
    if (radius_major <= radius_minor)
        throw py::value_error("radius_major must be greater than radius_minor");
    angular_sweep(start_angle, end_angle);
    const SmVector3d pivot = vector3(center);
    SmBrep* result = nullptr;
    check(SmApiCreatePartialTorus(pivot, radius_major, radius_minor,
                                  start_angle, end_angle, result), "SmApiCreatePartialTorus");
    return orient_axis(result, pivot, axis);
}

Brep make_ellipsoid(py::tuple center, double radius_x, double radius_y, double radius_z)
{
    positive(radius_x, "radius_x");
    positive(radius_y, "radius_y");
    positive(radius_z, "radius_z");
    const SmVector3d origin = vector3(center);
    SmBrep* raw = nullptr;
    check(SmApiCreateSphere(origin, 1.0, raw), "SmApiCreateSphere");
    std::unique_ptr<SmBrep> result(raw);
    check(SmApiScaleByPt(result.get(), origin, SmVector3d(radius_x, radius_y, radius_z)),
          "SmApiScaleByPt");
    return std::make_shared<BrepHandle>(result.release());
}

Brep make_sphere_no_pole(py::tuple center, double radius)
{
    positive(radius, "radius");
    SmBrep* result = nullptr;
    check(SmApiCreateSphereNoPole(vector3(center), radius, result), "SmApiCreateSphereNoPole");
    return std::make_shared<BrepHandle>(result);
}

Brep make_plane(py::tuple origin, double length, double width)
{
    positive(length, "length");
    positive(width, "width");
    SmBrep* result = nullptr;
    check(SmApiCreatePlane(vector3(origin), length, width, result), "SmApiCreatePlane");
    return std::make_shared<BrepHandle>(result);
}

Brep make_planar_circle(py::tuple origin, double radius)
{
    positive(radius, "radius");
    SmBrep* result = nullptr;
    check(SmApiCreatePlanarCircle(vector3(origin), radius, result), "SmApiCreatePlanarCircle");
    return std::make_shared<BrepHandle>(result);
}

using BooleanFunction = SmApiStatus (*)(SmBrep*, SmBrep*, SmBrep*&);

Brep boolean_op(const Brep& a, const Brep& b, BooleanFunction operation, const char* name)
{
    if (a == b)
        throw py::value_error("Boolean operands must differ");
    SmBrep* primary = a->get();
    SmBrep* secondary = b->get();
    SmBrep* result = nullptr;
    const SmApiStatus status = operation(primary, secondary, result);
    if (status == SM_SUCCESS)
    {
        // The SM_API contract deletes B on success and returns A in place.
        b->value.release();
        if (result != primary)
            throw std::runtime_error(std::string(name) + " returned an unexpected result pointer");
    }
    check(status, name);
    return a;
}

Brep chamfer_fillet(const Brep& brep, double width)
{
    positive(width, "width");
    check(SmApiChamferFillet(brep->get(), width), "SmApiChamferFillet");
    return brep;
}

Brep circular_fillet(const Brep& brep, double radius)
{
    positive(radius, "radius");
    check(SmApiCircularFillet(brep->get(), radius), "SmApiCircularFillet");
    return brep;
}

Brep transform_brep(const Brep& brep, py::tuple scale, py::tuple euler, py::tuple position)
{
    const SmVector3d factors = vector3(scale);
    const SmVector3d angles = vector3(euler);
    const SmVector3d offset = vector3(position);
    positive(factors.x, "scale x");
    positive(factors.y, "scale y");
    positive(factors.z, "scale z");
    SmBrep* solid = brep->get();
    if (factors.x != 1.0 || factors.y != 1.0 || factors.z != 1.0)
        check(SmApiScale(solid, factors), "SmApiScale");
    const SmVector3d origin(0.0, 0.0, 0.0);
    if (angles.x != 0.0)
        check(SmApiRotate(solid, origin, SmVector3d(1.0, 0.0, 0.0), angles.x), "SmApiRotate X");
    if (angles.y != 0.0)
        check(SmApiRotate(solid, origin, SmVector3d(0.0, 1.0, 0.0), angles.y), "SmApiRotate Y");
    if (angles.z != 0.0)
        check(SmApiRotate(solid, origin, SmVector3d(0.0, 0.0, 1.0), angles.z), "SmApiRotate Z");
    if (offset.x != 0.0 || offset.y != 0.0 || offset.z != 0.0)
        check(SmApiTranslate(solid, offset), "SmApiTranslate");
    return brep;
}

PolyBrep tessellate(const Brep& brep, double chord, double curve_angle, double surface_angle)
{
    if (!std::isfinite(chord) || chord < 0 || !std::isfinite(curve_angle) || curve_angle < 0 ||
        !std::isfinite(surface_angle) || surface_angle < 0)
        throw py::value_error("Tessellation tolerances must be finite and non-negative");
    SmTessellationParams params;
    params.dChordHeightTol = chord;
    params.dCurveAngleTolDeg = curve_angle;
    params.dSurfaceAngleTolDeg = surface_angle;
    SmPolyBrep* result = nullptr;
    check(SmApiTessellate(brep->get(), result, params), "SmApiTessellate");
    if (!result)
        return nullptr;
    return std::make_shared<PolyHandle>(result);
}

py::dict mesh_arrays(const PolyHandle& handle)
{
    if (!handle.value)
        throw py::value_error("Invalid PolyBrep");
    SmTArray<SmPolyVertex*> vertices;
    SmTArray<SmPolyFace*> faces;
    handle.value->GetPolyVertices(vertices);
    handle.value->GetPolyFaces(faces);
    std::unordered_map<SmPolyVertex*, size_t> vertex_index;
    py::list points, face_stream, face_normals, corner_normals;
    for (ULONG i = 0; i < vertices.GetSize(); ++i)
    {
        SmPolyVertex* vertex = vertices[i];
        const SmPoint3d& point = vertex->GetPoint();
        vertex_index[vertex] = i;
        points.append(py::make_tuple(point.x, point.y, point.z));
    }
    for (ULONG i = 0; i < faces.GetSize(); ++i)
    {
        SmPolyFace* face = faces[i];
        SmPolyLoop* loop = face->GetOuterPolyLoop();
        if (!loop)
            continue;
        SmTArray<SmPolyEdge*> edges;
        loop->GetPolyEdges(edges);
        if (edges.GetSize() < 3)
            continue;
        py::list polygon, polygon_normals;
        bool complete = true;
        const SmVector3d flat_normal = face->GetNormal(TRUE, FALSE);
        SmFace* original_face = face->GetOriginalFace();
        for (ULONG j = 0; j < edges.GetSize(); ++j)
        {
            SmPolyVertex* vertex = edges[j]->GetStartPolyVertex();
            auto it = vertex_index.find(vertex);
            if (it == vertex_index.end())
            {
                complete = false;
                break;
            }
            polygon.append(it->second);
            // Tessellation stores a separate surface normal for each original
            // BRep face touching this vertex. Keep that association so sharp
            // edges remain sharp while a curved face shades smoothly.
            SmVector3d normal = flat_normal;
            SmTArray<SmTopology*>& linked_faces = vertex->GetFacesRef();
            SmTArray<SmVector3d>& linked_normals = vertex->GetNormalsRef();
            for (ULONG k = 0; k < linked_faces.GetSize() && k < linked_normals.GetSize(); ++k)
            {
                if (linked_faces[k] == face || (original_face && linked_faces[k] == original_face))
                {
                    normal = linked_normals[k];
                    break;
                }
            }
            const double length = normal.Length();
            if (length > 1e-12)
                polygon_normals.append(py::make_tuple(normal.x / length, normal.y / length, normal.z / length));
            else
                polygon_normals.append(py::make_tuple(flat_normal.x, flat_normal.y, flat_normal.z));
        }
        if (!complete)
            continue;
        face_stream.append(edges.GetSize());
        for (auto item : polygon)
            face_stream.append(item);
        face_normals.append(py::make_tuple(flat_normal.x, flat_normal.y, flat_normal.z));
        corner_normals.append(polygon_normals);
    }
    py::dict result;
    result["points"] = points;
    result["faces"] = face_stream;
    result["face_normals"] = face_normals;
    result["corner_normals"] = corner_normals;
    return result;
}

} // namespace

PYBIND11_MODULE(_smlib_core, module)
{
    module.doc() = "USD-free SMLib solid modeling subset for Python 3.12.";
    if (SmApiGetOrCreateContext() == nullptr)
        SmApiCreateContext();

    py::class_<BrepHandle, Brep>(module, "Brep")
        .def_property_readonly("valid", [](const BrepHandle& value) { return bool(value.value); });
    py::class_<PolyHandle, PolyBrep>(module, "PolyBrep")
        .def("to_numpy_arrays", &mesh_arrays,
             "Return points, faces, oriented face normals, and per-corner surface normals.");

    py::class_<Profile>(module, "Profile")
        .def_property_readonly("points", [](const Profile& p) { return p.points; },
                               "Return a copy of the normalized counter-clockwise 2D boundary.");
    module.def("create_profile", &make_profile, py::arg("points"),
               "Pure: create a reusable closed polygon Profile; input points are unchanged.\n\n"
               "Accept 3 to 1024 finite 2D vertices without holes or self-intersections.\n"
               "The closing edge is implicit; a repeated last point is accepted. Winding is normalized.\n"
               "Validation is performed in this binding; no kernel geometry is allocated.");
    module.def("create_extrude", [](const Profile& p, py::tuple origin, double start, double end,
                                    py::object axis, py::object xdir) {
        return make_profile_solid(p, origin, start, end, axis, xdir, false);
    }, py::arg("profile"), py::arg("origin"), py::arg("start"), py::arg("end"),
       py::arg("axis") = py::none(), py::arg("xdir") = py::none(),
       "Pure: returns new capped Brep; Profile unchanged.\n\n"
       "Extrude local XY profile from start to end along axis (default +Z), start < end.\n"
       "xdir sets local +X after projection perpendicular to axis; default is shortest +Z rotation.\n\n"
       "Wraps: SmApiCreateLineSegment, SmApiCreateLinearSweep");
    module.def("create_revolve", [](const Profile& p, py::tuple origin, double start, double end,
                                    py::object axis, py::object xdir) {
        return make_profile_solid(p, origin, start, end, axis, xdir, true);
    }, py::arg("profile"), py::arg("origin"), py::arg("start_angle_deg"), py::arg("end_angle_deg"),
       py::arg("axis") = py::none(), py::arg("xdir") = py::none(),
       "Pure: returns new capped Brep; Profile unchanged.\n\n"
       "Revolve (radius, height) profile by 0 < end-start <= 360 degrees about axis (default +Z).\n"
       "Radii must be non-negative; axis contact is allowed. xdir defines the zero-angle ray.\n"
       "Positive angles follow the right-hand rule. Partial rotations are capped.\n\n"
       "Wraps: SmApiCreateLineSegment, SmApiCreateRotationalSweep");

    module.def("create_box", &make_box, py::arg("origin"), py::arg("length"),
               py::arg("width"), py::arg("height"),
               "Create a solid box.\n\nWraps: SmApiCreateBox");
    module.def("create_sphere", &make_sphere, py::arg("center"), py::arg("radius"),
               "Create a solid sphere.\n\nWraps: SmApiCreateSphere");
    module.def("create_cylinder", &make_cylinder, py::arg("base"), py::arg("radius"),
               py::arg("height"), py::arg("axis") = py::none(),
               "Create a capped cylinder along axis (default +Z).\n\nWraps: SmApiCreateCylinder, SmApiRotate");
    module.def("create_cone", &make_cone, py::arg("base"), py::arg("radius_base"),
               py::arg("radius_top"), py::arg("height"), py::arg("axis") = py::none(),
               "Create a cone or frustum along axis (default +Z).\n\nWraps: SmApiCreateCone, SmApiRotate");
    module.def("create_torus", &make_torus, py::arg("center"), py::arg("radius_major"),
               py::arg("radius_minor"), py::arg("axis") = py::none(),
               "Create a torus about axis (default +Z).\n\nWraps: SmApiCreateTorus, SmApiRotate");
    module.def("create_pyramid", &make_pyramid, py::arg("base"), py::arg("length"),
               py::arg("height"), py::arg("axis") = py::none(),
               "Create a square pyramid along axis (default +Z).\n\nWraps: SmApiCreatePyramid, SmApiRotate");
    module.def("create_cylindrical_box", &make_cylindrical_box, py::arg("base"), py::arg("length"),
               py::arg("inside_radius"), py::arg("outside_radius"), py::arg("start_angle_deg"),
               py::arg("end_angle_deg"), py::arg("axis") = py::none(),
               "Create a hollow cylinder sector along axis (default +Z).\n\nWraps: SmApiCreateCylindricalBox, SmApiRotate");
    module.def("create_partial_sphere", &make_partial_sphere, py::arg("center"), py::arg("radius"),
               py::arg("start_angle_deg"), py::arg("end_angle_deg"), py::arg("axis") = py::none(),
               "Create a spherical azimuth sector about axis (default +Z).\n\nWraps: SmApiCreatePartialSphere, SmApiRotate");
    module.def("create_partial_cone", &make_partial_cone, py::arg("base"), py::arg("radius_base"),
               py::arg("radius_top"), py::arg("height"), py::arg("start_angle_deg"),
               py::arg("end_angle_deg"), py::arg("axis") = py::none(),
               "Create a cone sector along axis (default +Z).\n\nWraps: SmApiCreatePartialCone, SmApiRotate");
    module.def("create_partial_cylinder", &make_partial_cylinder, py::arg("base"), py::arg("radius"),
               py::arg("height"), py::arg("start_angle_deg"), py::arg("end_angle_deg"),
               py::arg("axis") = py::none(),
               "Create a cylinder sector along axis (default +Z).\n\nWraps: SmApiCreatePartialCylinder, SmApiRotate");
    module.def("create_partial_torus", &make_partial_torus, py::arg("center"), py::arg("radius_major"),
               py::arg("radius_minor"), py::arg("start_angle_deg"), py::arg("end_angle_deg"),
               py::arg("axis") = py::none(),
               "Create a torus sector about axis (default +Z).\n\nWraps: SmApiCreatePartialTorus, SmApiRotate");
    module.def("create_ellipsoid", &make_ellipsoid, py::arg("center"), py::arg("radius_x"),
               py::arg("radius_y"), py::arg("radius_z"),
               "Create an ellipsoid by scaling an exact sphere.\n\nWraps: SmApiCreateSphere, SmApiScaleByPt");
    module.def("create_sphere_no_pole", &make_sphere_no_pole, py::arg("center"), py::arg("radius"),
               "Create a sphere without pole singularities.\n\nWraps: SmApiCreateSphereNoPole");
    module.def("create_plane", &make_plane, py::arg("origin"), py::arg("length"), py::arg("width"),
               "Create a planar face (not a solid).\n\nWraps: SmApiCreatePlane");
    module.def("create_planar_circle", &make_planar_circle, py::arg("origin"), py::arg("radius"),
               "Create a circular planar face (not a solid).\n\nWraps: SmApiCreatePlanarCircle");
    module.def("boolean_union", [](const Brep& a, const Brep& b) {
        return boolean_op(a, b, SmApiBooleanUnion, "SmApiBooleanUnion");
    }, py::arg("a"), py::arg("b"), "Consume B and union it into A.\n\nWraps: SmApiBooleanUnion");
    module.def("boolean_difference", [](const Brep& a, const Brep& b) {
        return boolean_op(a, b, SmApiBooleanDifference, "SmApiBooleanDifference");
    }, py::arg("a"), py::arg("b"), "Consume B and subtract it from A.\n\nWraps: SmApiBooleanDifference");
    module.def("boolean_intersection", [](const Brep& a, const Brep& b) {
        return boolean_op(a, b, SmApiBooleanIntersection, "SmApiBooleanIntersection");
    }, py::arg("a"), py::arg("b"), "Consume B and intersect it with A.\n\nWraps: SmApiBooleanIntersection");
    module.def("chamfer_fillet", &chamfer_fillet, py::arg("brep"), py::arg("radius"),
               "Chamfer all eligible edges in place and return the Brep. Width is measured across the chamfer face.\n\nWraps: SmApiChamferFillet");
    module.def("circular_fillet", &circular_fillet, py::arg("brep"), py::arg("radius"),
               "Round all eligible edges in place and return the Brep.\n\nWraps: SmApiCircularFillet");
    module.def("transform_brep", &transform_brep, py::arg("brep"), py::arg("scale"),
               py::arg("euler"), py::arg("position"),
               "Mutate a Brep using local scale, X/Y/Z Euler degrees, then parent-space translation.\n\n"
               "Wraps: SmApiScale, SmApiRotate, SmApiTranslate");
    module.def("tessellate", &tessellate, py::arg("brep"),
               py::arg("chord_height_tolerance") = 0.0,
               py::arg("curve_angle_tolerance_deg") = 25.0,
               py::arg("surface_angle_tolerance_deg") = 25.0,
               "Tessellate a Brep without changing it.\n\nWraps: SmApiTessellate");
}
