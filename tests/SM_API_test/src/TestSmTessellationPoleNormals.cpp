// SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef SM_TESS_POLE_TEST_STANDALONE
#include "StdAfx.h"
#endif
#include <SmApiBrep.h>
#include "SM_API_test.h"
#include <SmApiCurves.h>
#include <SmApiGeneral.h>
#include <SmApiQueries.h>
#include <SmBrep.h>
#include <SmGraphicsExtern.h>
#include <SmLine.h>
#include <SmPoly.h>
#include <SmTess.h>

#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
void Require(bool ok, const char* message)
{
    if (!ok) throw std::runtime_error(message);
}

void Check(SmStatus status, const char* message)
{
    if (status != SM_SUCCESS)
        throw std::runtime_error(std::string(message) + ": " + std::to_string(status));
}

class PoleProbe : public SmPolygonOutputCallback
{
    bool m_xAxis;

    SmPoint3d Local(const SmPoint3d& p) const
    {
        return m_xAxis ? SmPoint3d(p.y, p.z + 5, p.x) : SmPoint3d(p.x, p.y, p.z + 5);
    }

    SmVector3d LocalNormal(const SmVector3d& n) const
    {
        return m_xAxis ? SmVector3d(n.y, n.z, n.x) : n;
    }

public:
    unsigned inner = 0, outer = 0, meshCalls = 0;

    PoleProbe(SmPolygonOutputType mode, bool xAxis) : SmPolygonOutputCallback(mode), m_xAxis(xAxis) {}

    void Verify(ULONG count, const SmPoint3d* input, const SmVector3d* normals)
    {
        bool hasApex = false;
        std::vector<SmPoint3d> points;
        for (ULONG i = 0; i < count; ++i)
        {
            points.push_back(Local(input[i]));
            hasApex |= std::hypot(points.back().x, points.back().y) < 1e-8;
            Require(std::abs(normals[i].Length() - 1.0) < 1e-6, "Non-unit corner normal");
        }
        if (!hasApex) return;
        for (int side = 0; side != 2; ++side)
        {
            const double slope = side == 0 ? 5.0 / 14.0 : 9.0 / 4.0;
            const double sign = side == 0 ? -1.0 : 1.0;
            bool onCone = true;
            SmVector3d ray(0, 0, 0);
            for (const auto& p : points)
            {
                onCone &= std::abs(std::hypot(p.x, p.y) - (p.z + 2) * slope) < 1e-6;
                ray += SmVector3d(p.x, p.y, 0);
            }
            if (!onCone) continue;
            Require(ray.Length() > 0.0, "Degenerate cone polygon");
            ray *= 1.0 / ray.Length();
            for (ULONG i = 0; i < count; ++i)
            {
                const auto& p = points[i];
                const double radius = std::hypot(p.x, p.y);
                SmVector3d radial = radius < 1e-8 ? ray : SmVector3d(p.x / radius, p.y / radius, 0);
                SmVector3d expected(sign * radial.x, sign * radial.y, -sign * slope);
                expected *= 1.0 / expected.Length();
                Require(expected.Dot(LocalNormal(normals[i])) > 0.99999, "Incorrect directional pole normal");
                if (radius < 1e-8)
                {
                    if (side == 0) ++inner;
                    else ++outer;
                }
            }
        }
    }

    SmStatus OutputPolygon(ULONG type, ULONG count, SmPoint3d* points, SmVector3d* normals,
                           SmPoint2d*, SmSurface*, SmFace*, SmGfxArraySet*) override
    {
        // Non-pole surfaces, including planar cut faces through the axis,
        // can still use strips/fans. Check their constituent triangles too.
        if (type <= 1) Verify(count, points, normals);
        else
            for (ULONG i = 0; i + 2 < count; ++i)
            {
                const ULONG a = type == 3 ? 0 : i;
                SmPoint3d triangle[] = {points[a], points[i + 1], points[i + 2]};
                SmVector3d cornerNormals[] = {normals[a], normals[i + 1], normals[i + 2]};
                Verify(3, triangle, cornerNormals);
            }
        return SM_SUCCESS;
    }

    SmStatus OutputMesh(SmTArray<ULONG>& counts, SmTArray<ULONG>& indices, SmTArray<SmPoint3d>& points,
                        SmTArray<SmVector3d>& normals, SmTArray<SmPoint2d>&,
                        SmSurface*, SmFace*, SmGfxArraySet*) override
    {
        ++meshCalls;
        ULONG offset = 0;
        for (ULONG i = 0; i < counts.GetSize(); ++i)
        {
            std::vector<SmPoint3d> polygon;
            std::vector<SmVector3d> cornerNormals;
            for (ULONG j = 0; j < counts[i]; ++j)
            {
                const ULONG index = indices[offset++];
                polygon.push_back(points[index]);
                cornerNormals.push_back(normals[index]);
            }
            Verify(counts[i], polygon.data(), cornerNormals.data());
        }
        return SM_SUCCESS;
    }
};

std::unique_ptr<SmBrep> Fixture(bool xAxis, double start, double end)
{
    const double rz[][2] = {{0, -2}, {9, 2}, {9, 10}, {5, 12}};
    const double angle = start * std::acos(-1.0) / 180.0;
    SmPoint3d points[4];
    for (int i = 0; i < 4; ++i)
    {
        const double x = rz[i][0] * std::cos(angle), y = rz[i][0] * std::sin(angle);
        points[i] = xAxis ? SmPoint3d(rz[i][1], x, y - 5) : SmPoint3d(x, y, rz[i][1] - 5);
    }
    SmTArray<SmCurve*> curves;
    std::vector<std::unique_ptr<SmLine>> owned;
    for (int i = 0; i < 4; ++i)
    {
        SmLine* line = nullptr;
        Check(SmApiCreateLineSegment(points[i], points[(i + 1) % 4], line), "Create line");
        owned.emplace_back(line);
        curves.Add(line);
    }
    SmBrep* raw = nullptr;
    Check(SmApiCreateRotationalSweep(curves, SmPoint3d(0, 0, -5),
                                     xAxis ? SmVector3d(1, 0, 0) : SmVector3d(0, 0, 1),
                                     end - start, TRUE, raw), "Revolve");
    return std::unique_ptr<SmBrep>(raw);
}

void RunChecks()
{
    const SmContext* context = SmApiGetOrCreateContext();
    Require(context != nullptr, "Missing context");
    for (bool xAxis : {false, true})
        for (bool partial : {false, true})
        {
            const double start = partial ? 30 : 0, end = partial ? 270 : 360;
            auto body = Fixture(xAxis, start, end);
            SmBoolean manifold = FALSE;
            Check(SmApiBrepIsManifoldSolid(body.get(), manifold), "Manifold query");
            Require(manifold == TRUE, "Invalid fixture");
            double beforeVolume = 0;
            Check(SmApiBrepComputeVolume(body.get(), 1e-7, beforeVolume), "Volume query");
            const double expected = 740 * std::acos(-1.0) * (end - start) / 360;
            Require(std::abs(beforeVolume - expected) < expected * 1e-5, "Incorrect fixture volume");

            SmPolyBrep* raw = nullptr;
            SmTessellationParams params;
            params.dChordHeightTol = 0.08;
            params.dCurveAngleTolDeg = params.dSurfaceAngleTolDeg = 20;
            Check(SmApiTessellate(body.get(), raw, params), "API tessellation");
            std::unique_ptr<SmPolyBrep> mesh(raw);
            SmTArray<SmPolyFace*> faces;
            mesh->GetPolyFaces(faces);
            PoleProbe probe(SM_PO_TRIANGLES, xAxis);
            for (ULONG i = 0; i < faces.GetSize(); ++i)
            {
                SmTArray<SmPolyEdge*> edges;
                faces[i]->GetOuterPolyLoop()->GetPolyEdges(edges);
                std::vector<SmPoint3d> points;
                std::vector<SmVector3d> normals;
                for (ULONG j = 0; j < edges.GetSize(); ++j)
                {
                    points.push_back(edges[j]->GetStartPolyVertex()->GetPoint());
                    Require(edges[j]->GetAuxData() != nullptr, "Missing corner data");
                    normals.push_back(edges[j]->GetAuxData()->m_vNormal);
                }
                probe.Verify(edges.GetSize(), points.data(), normals.data());
            }
            Require(probe.inner > 0 && probe.outer > 0, "No PolyBrep pole corners checked");

            for (auto mode : {SM_PO_TRIANGLES, SM_PO_QUADRALATERALS, SM_PO_TRIANGLE_MESH,
                              SM_PO_TRIANGLE_STRIPS, SM_PO_TRIANGLE_FANS, SM_PO_STRIPS_OR_FANS})
            {
                SmCurveTessDriver curve(0, 0.08, 20);
                SmSurfaceTessDriver surface(0.08, 20);
                SmTess tess(*context, curve, surface);
                SmBoolean failed = FALSE;
                // SmTess owns and may modify this private working copy.
                Check(tess.DoTessellation(new (*context) SmBrep(*body), failed), "Kernel tessellation");
                Require(failed == FALSE, "A fixture face failed");
                PoleProbe callback(mode, xAxis);
                Check(tess.OutputPolygons(callback), "Callback output");
                Require(callback.inner > 0 && callback.outer > 0, "No callback pole corners checked");
                if (mode == SM_PO_TRIANGLE_MESH)
                    Require(callback.meshCalls > 0, "Indexed mesh callback was bypassed");
            }
            double afterVolume = 0;
            Check(SmApiBrepComputeVolume(body.get(), 1e-7, afterVolume), "Source volume after tessellation");
            Require(std::abs(beforeVolume - afterVolume) < 1e-8, "Source geometry changed");
        }
}
} // namespace

SmStatus TestSmTessellationPoleNormals()
{
    try
    {
        RunChecks();
        return SM_SUCCESS;
    }
    catch (const std::exception& error)
    {
        std::cerr << "Tessellation pole normals: " << error.what() << '\n';
        return SM_ERR;
    }
}

#ifdef SM_TESS_POLE_TEST_STANDALONE
int main()
{
    const SmStatus status = TestSmTessellationPoleNormals();
    if (status == SM_SUCCESS)
        std::cout << "Pole normals: PolyBrep and six callback modes passed for full/partial revolutions and Z/X axes\n";
    return status == SM_SUCCESS ? 0 : 1;
}
#endif
