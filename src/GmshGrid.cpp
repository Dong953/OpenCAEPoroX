/*! \file    GmshGrid.cpp
 *  \brief   GmshGrid class declaration
 *  \author  Shizhe Li
 *  \date    Sep/10/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */


#include "../config/config.hpp"

#ifdef OCP_USE_GMSH


#include "GmshGrid.hpp"
#include <map>
OCP_Polygon::OCP_Polygon(const vector<OCP_ULL>& pIndex, const OCP_ULL& tag_in, const string& phyinfo, const OCP_USI& index)
{
    p        = pIndex;
    tag      = tag_in;
    physical = phyinfo;
    phyIndex = index;
}


void OCP_Polygon::CalCenter(const vector<OCP_DBL>& points)
{
    const USI np = p.size();
    center.Reset();
    for (USI i = 0; i < np; i++) {
        center.x += points[3 * p[i] + 0];
        center.y += points[3 * p[i] + 1];
        center.z += points[3 * p[i] + 2];
    }
    center *= 1.0 / np;
}


void OCP_Polygon::CalArea(const vector<OCP_DBL>& points)
{
    if (p.size() == 3) {
        const OCP_DBL* p0 = &points[3 * p[0]];
        const OCP_DBL* p1 = &points[3 * p[1]];
        const OCP_DBL* p2 = &points[3 * p[2]];
        const USI      x = 0;
        const USI      y = 1;
        area = 0.5 * fabs((p2[x] - p1[x]) * (p0[y] - p1[y]) - (p2[y] - p1[y]) * (p0[x] - p1[x]));
    }
    else {
        const OCP_DBL* p0 = &points[3 * p[0]];
        const OCP_DBL* p1 = &points[3 * p[1]];
        const OCP_DBL* p2 = &points[3 * p[2]];
        const OCP_DBL* p3 = &points[3 * p[3]];
        const USI      x = 0;
        const USI      y = 1;
        area = 0.5 * (fabs((p2[x] - p1[x]) * (p0[y] - p1[y]) - (p2[y] - p1[y]) * (p0[x] - p1[x]))
                      + fabs((p0[x] - p3[x]) * (p2[y] - p3[y]) - (p0[y] - p3[y]) * (p2[x] - p3[x])));
    }
}


OCP_BOOL OCP_Polygon::IfPointInElement(const Point3D& objP, const vector<OCP_DBL>& points) {
    OCP_DBL tmpArea = 0;
    const USI numP = p.size();
    if (numP >= 5) {
        auto TetraVol = [](const Point3D& a, const Point3D& b, const Point3D& c, const Point3D& d) {
            return (CrossProduct(b - a , c - a) * (d - a) / 6.0);
        };

        auto InTetra = [&](const Point3D& a, const Point3D& b, const Point3D& c, const Point3D& d) -> bool {
            OCP_DBL v0 = TetraVol(a, b, c, d);
            if (fabs(v0) < TINY) return false;

            OCP_DBL v1 = TetraVol(objP, b, c, d);
            OCP_DBL v2 = TetraVol(a, objP, c, d);
            OCP_DBL v3 = TetraVol(a, b, objP, d);
            OCP_DBL v4 = TetraVol(a, b, c, objP);

            OCP_DBL sumV = v1 + v2 + v3 + v4;
            return (fabs(sumV - v0) < TINY) && (v1 * v2 > 0 && v1 * v3 > 0 && v1 * v4 > 0);
        };

        // 提取 8 个点
        Point3D tmpPoint[8];
        for (int i = 0; i < numP; ++i) {
            tmpPoint[i] = Point3D(&points[3 * p[i]]);
        }

        // 划分为 5 个四面体
        vector<array<int, 4>> tets = {
                {0, 1, 3, 4},
                {1, 2, 3, 6},
                {1, 5, 4, 6},
                {3, 4, 7, 6},
                {1, 3, 4, 6}
        };

        for (const auto& tet : tets) {
            if (InTetra(tmpPoint[tet[0]], tmpPoint[tet[1]], tmpPoint[tet[2]], tmpPoint[tet[3]])) {
                return OCP_TRUE;
            }
        }
        return OCP_FALSE;
//        return OCP_TRUE;

    }
    else {
        for (USI i = 0; i < numP; i++) {
            const Point3D tmpPoint = CrossProduct(objP - Point3D(&points[3 * p[i % numP]]),
                                                  objP - Point3D(&points[3 * p[(i + 1) % numP]]));
            tmpArea += sqrt(tmpPoint * tmpPoint);
        }
        tmpArea *= 0.5;

        if (fabs(area - tmpArea) < TINY) {
            return OCP_TRUE;
        } else {
            return OCP_FALSE;
        }
    }
}

void OCP_Polygon::CalArea3D(const vector<OCP_DBL>& points)
{
    if (p.size() == 4) {
        const OCP_DBL* p0 = &points[3 * p[0]];
        const OCP_DBL* p1 = &points[3 * p[1]];
        const OCP_DBL* p2 = &points[3 * p[2]];
        const OCP_DBL* p3 = &points[3 * p[3]];
        const USI      x = 0;
        const USI      y = 1;
        const USI      z = 2;
        area = fabs(
                (p1[0] - p0[0]) * ((p2[1] - p0[1]) * (p3[2] - p0[2]) - (p2[2] - p0[2]) * (p3[1] - p0[1])) +
                (p1[1] - p0[1]) * ((p2[2] - p0[2]) * (p3[0] - p0[0]) - (p2[0] - p0[0]) * (p3[2] - p0[2])) +
                (p1[2] - p0[2]) * ((p2[0] - p0[0]) * (p3[1] - p0[1]) - (p2[1] - p0[1]) * (p3[0] - p0[0]))
        ) / 6.0;
    }
    else {
        const OCP_DBL* p0 = &points[3 * p[0]];
        const OCP_DBL* p1 = &points[3 * p[1]];
        const OCP_DBL* p2 = &points[3 * p[2]];
        const OCP_DBL* p3 = &points[3 * p[3]];
        const OCP_DBL* p4 = &points[3 * p[4]];
        const OCP_DBL* p5 = &points[3 * p[5]];
        const OCP_DBL* p6 = &points[3 * p[6]];
        const OCP_DBL* p7 = &points[3 * p[7]];
        auto TetraVol = [](const OCP_DBL* a, const OCP_DBL* b, const OCP_DBL* c, const OCP_DBL* d) {
            OCP_DBL ab[3] = { b[0] - a[0], b[1] - a[1], b[2] - a[2] };
            OCP_DBL ac[3] = { c[0] - a[0], c[1] - a[1], c[2] - a[2] };
            OCP_DBL ad[3] = { d[0] - a[0], d[1] - a[1], d[2] - a[2] };

            OCP_DBL cx = ac[1]*ad[2] - ac[2]*ad[1];
            OCP_DBL cy = ac[2]*ad[0] - ac[0]*ad[2];
            OCP_DBL cz = ac[0]*ad[1] - ac[1]*ad[0];

            return fabs(ab[0]*cx + ab[1]*cy + ab[2]*cz) / 6.0;
        };
        area =
                TetraVol(p0, p1, p3, p4) +   // T0
                TetraVol(p1, p2, p3, p6) +   // T1
                TetraVol(p1, p4, p5, p6) +   // T2
                TetraVol(p3, p4, p6, p7) +   // T3
                TetraVol(p1, p3, p4, p6);    // T4
    }
}

void GMSHGrid::InputGrid(const string& file)
{
    ifUse = OCP_TRUE;

    gmsh::initialize();
    gmsh::open(file);

    // Print the model name and dimension:
    std::string name;
    gmsh::model::getCurrent(name);
    std::cout << "Model " << name << " (" << gmsh::model::getDimension() << "D)\n";
    dimen = gmsh::model::getDimension();

    if (dimen == 2) {
        InputGrid2D(file);
    }
    if (dimen == 3) {
        InputGrid3D(file);
    }
    gmsh::finalize();

    Setup();
}

void GMSHGrid::InputGrid2D(const string& file)
{

    // Get all mesh nodes
    std::vector<std::size_t> nodeTags;
    std::vector<double> nodeParams;
    vector<OCP_DBL>     pointsTmp;

#if OCPFLOATTYPEWIDTH == 128

    vector<double>        pointsTmp;
	gmsh::model::mesh::getNodes(nodeTags, pointsTmp, nodeParams, -1, -1);

	points.resize(pointsTmp.size());
	for (OCP_ULL i = 0; i < pointsTmp.size(); i++) {
		points[i] = static_cast<OCP_DBL>(pointsTmp[i]);
	}
	vector<double>().swap(pointsTmp);

#else
    gmsh::model::mesh::getNodes(nodeTags, pointsTmp, nodeParams, -1, -1);
#endif

    // re-order points
    points = pointsTmp;
    for (OCP_ULL i = 0; i < nodeTags.size(); i++) {
        copy(&pointsTmp[i * 3], &pointsTmp[i * 3] + 3, &points[(nodeTags[i] - 1) * 3]);
    }
    vector<OCP_DBL>().swap(pointsTmp);


    // from m to cm for spe11a
    for (auto& p : points) {
        p *= 100;
    }

    // Get all the elementary entities in the model, as a vector of (dimension, tag) pairs:
    gmsh::vectorpair entities;
    gmsh::model::getEntities(entities);

    OCP_ULL lineTag   = 1;
    OCP_ULL faceIndex = 0;

    /// find all physical name, ordered as the ones of definitions in *.geo
    physicalNameSet.resize(dimen + 1);
    for (USI d = 0; d <= dimen; d++) {
        gmsh::vectorpair dimTags;
        string           tmp;
        gmsh::model::getPhysicalGroups(dimTags, d);
        for (const auto& p : dimTags) {
            gmsh::model::getPhysicalName(p.first, p.second, tmp);
            physicalNameSet[d].push_back(tmp);
        }
    }
    /// Allocate facies
    for (const auto& p2 : physicalNameSet[2]) {
        facies.push_back(Facies(p2));
    }

    for (const auto& e : entities) {
        // Dimension and tag of the entity:
        const int dim = e.first, tag = e.second;


        std::vector<int> physicalTags;
        gmsh::model::getPhysicalGroupsForEntity(dim, tag, physicalTags);

        if (physicalTags.empty()) {
            continue;
        }

        string physicalName;
        gmsh::model::getPhysicalName(dim, physicalTags[0], physicalName);

        // Get the mesh elements for the entity (dim, tag):
        std::vector<int> elemTypes;
        std::vector<std::vector<std::size_t> > elemTags, elemNodeTags;
        gmsh::model::mesh::getElements(elemTypes, elemTags, elemNodeTags, dim, tag);


        INT physicalIndex = -1;
        for (INT p = 0; p < physicalNameSet[dim].size(); p++) {
            if (physicalNameSet[dim][p] == physicalName) {
                physicalIndex = p;
                break;
            }
        }
        if (physicalIndex == -1) {
            OCP_ABORT("No Matched physical Body!");
        }


        if (dim == 1) {
            // for boundary lines
            for (USI t = 0; t < elemTypes.size(); t++) {
                for (OCP_ULL l = 0; l < elemTags[t].size(); l++) {
                    const OCP_ULL bId = elemNodeTags[t][2 * l];
                    const OCP_ULL eId = elemNodeTags[t][2 * l + 1];
                    edges.insert(Edge(bId - 1, eId - 1, elemTags[t][l], physicalName, physicalIndex));
                    lineTag++;
                }
            }
        }
        else if (dim == 2) {
            // for triangle and quadrangle
            USI np = 0;
            for (USI t = 0; t < elemTypes.size(); t++) {
                if (elemTypes[t] == 2)  np = 3;  // for triangle
                else                    np = 4;  // for quadrangle

                vector<OCP_ULL> indexPoints(np);
                for (OCP_ULL l = 0; l < elemTags[t].size(); l++) {
                    indexPoints.clear();

                    for (USI i = 0; i < np; i++) {
                        const OCP_ULL bId = elemNodeTags[t][np * l + (i % np)];
                        const OCP_ULL eId = elemNodeTags[t][np * l + (i + 1) % np];

                        auto iter = edges.find(Edge(bId - 1, eId - 1));
                        if (iter == edges.end()) {
                            edges.insert(Edge(bId - 1, eId - 1, lineTag++, faceIndex, i));
                        }
                        else {
                            iter->faceIndex.push_back(faceIndex);
                            iter->faceIndex.push_back(i);
                        }

                        indexPoints.push_back(elemNodeTags[t][np * l + i] - 1);
                    }
                    faceIndex++;
                    elements.push_back(OCP_Polygon(indexPoints, elemTags[t][l], physicalName, physicalIndex));

                    //cout << "------------------------------------" << endl;
                    //for (USI i = 0; i < np; i++) {
                    //	for (USI j = 0; j < 3; j++) {
                    //		cout << points[indexPoints[i] * 3 + j] << "   ";
                    //	}
                    //	cout << endl;
                    //}
                }
            }
        }
    }
}

void GMSHGrid::InputGrid3D(const string& file)
{
    // Get all mesh nodes
    std::vector<std::size_t> nodeTags;
    std::vector<double> nodeParams;
    vector<OCP_DBL>     pointsTmp;

#if OCPFLOATTYPEWIDTH == 128
    vector<double>        pointsTmp;
    gmsh::model::mesh::getNodes(nodeTags, pointsTmp, nodeParams, -1, -1);

    points.resize(pointsTmp.size());
    for (OCP_ULL i = 0; i < pointsTmp.size(); i++) {
        points[i] = static_cast<OCP_DBL>(pointsTmp[i]);
    }
    vector<double>().swap(pointsTmp);
#else
    gmsh::model::mesh::getNodes(nodeTags, pointsTmp, nodeParams, -1, -1);
#endif

    // re-order points
    points = pointsTmp;
    for (OCP_ULL i = 0; i < nodeTags.size(); i++) {
        copy(&pointsTmp[i * 3], &pointsTmp[i * 3] + 3, &points[(nodeTags[i] - 1) * 3]);
    }
    vector<OCP_DBL>().swap(pointsTmp);

    // from m to cm for spe11a
    for (auto& p : points) {
        p *= 100;
    }

    // Get all the elementary entities in the model, as a vector of (dimension, tag) pairs:
    gmsh::vectorpair entities;
    gmsh::model::getEntities(entities);

    OCP_ULL lineTag   = 1;
    OCP_ULL faceIndex = 0;

    /// find all physical name, ordered as the ones of definitions in *.geo
    physicalNameSet.resize(dimen + 1);
    for (USI d = 0; d <= dimen; d++) {
        gmsh::vectorpair dimTags;
        string           tmp;
        gmsh::model::getPhysicalGroups(dimTags, d);
        for (const auto& p : dimTags) {
            gmsh::model::getPhysicalName(p.first, p.second, tmp);
            physicalNameSet[d].push_back(tmp);
        }
    }
    /// Allocate facies
    for (const auto& p2 : physicalNameSet[3]) {
        facies.push_back(Facies(p2));
    }

    for (const auto& e : entities) {
        // Dimension and tag of the entity:
        const int dim = e.first, tag = e.second;

        std::vector<int> physicalTags;
        gmsh::model::getPhysicalGroupsForEntity(dim, tag, physicalTags);

        if (physicalTags.empty()) {
            continue;
        }

        string physicalName;
        gmsh::model::getPhysicalName(dim, physicalTags[0], physicalName);

        // Get the mesh elements for the entity (dim, tag):
        std::vector<int> elemTypes;
        std::vector<std::vector<std::size_t> > elemTags, elemNodeTags;
        gmsh::model::mesh::getElements(elemTypes, elemTags, elemNodeTags, dim, tag);

        INT physicalIndex = -1;
        for (INT p = 0; p < physicalNameSet[dim].size(); p++) {
            if (physicalNameSet[dim][p] == physicalName) {
                physicalIndex = p;
                break;
            }
        }
        if (physicalIndex == -1) {
            OCP_ABORT("No Matched physical Body!");
        }

        if (dim == 1) {
            // for boundary lines
            for (USI t = 0; t < elemTypes.size(); t++) {
                for (OCP_ULL l = 0; l < elemTags[t].size(); l++) {
                    const OCP_ULL bId = elemNodeTags[t][2 * l];
                    const OCP_ULL eId = elemNodeTags[t][2 * l + 1];
                    edges.insert(Edge(bId - 1, eId - 1, elemTags[t][l], physicalName, physicalIndex));
                    lineTag++;
                }
            }
        }
        else if (dim == 2) {
            // for boundary faces (triangles or quads)
            for (USI t = 0; t < elemTypes.size(); t++) {
                USI np = (elemTypes[t] == 2 ? 3 : 4); // triangle or quad
                vector<OCP_ULL> indexPoints(np);
                for (OCP_ULL l = 0; l < elemTags[t].size(); l++) {
                    indexPoints.clear();
                    for (USI i = 0; i < np; i++) {
                        indexPoints.push_back(elemNodeTags[t][np * l + i] - 1);
                    }
                    faces.insert(Face(indexPoints, physicalName, physicalIndex, elemTags[t][l]));
                }
            }
        }
        else if (dim == 3) {
            // for tetrahedron and hexahedron
            USI np = 0;
            for (USI t = 0; t < elemTypes.size(); t++) {
                if (elemTypes[t] == 4) np = 4;      // tetrahedron
                else if (elemTypes[t] == 5) np = 8; // hexahedron
                else continue; // unsupported

                vector<OCP_ULL> indexPoints(np);
                for (OCP_ULL l = 0; l < elemTags[t].size(); l++) {
                    indexPoints.clear();
                    for (USI i = 0; i < np; i++) {
                        indexPoints.push_back(elemNodeTags[t][np * l + i] - 1);
                    }
                    faceIndex++;
                    elements.push_back(OCP_Polygon(indexPoints, elemTags[t][l], physicalName, physicalIndex));
                }
            }
        }
    }
    if(faces.empty()) {
        // key：排序后的节点列表 → value：所有拥有此面的体元索引列表
        std::map<std::vector<OCP_ULL>, std::vector<OCP_ULL>> faceElems;

        // 枚举所有体元，按四面体/六面体提取面
        for(OCP_ULL ei = 0; ei < elements.size(); ++ei) {
            const auto& p = elements[ei].p;
            if(p.size() == 4) {
                // 四面体：4 个三角面
                static const int tetFaces[4][3] = {
                        {0,1,2},{0,1,3},{1,2,3},{2,0,3}
                };
                for(int f = 0; f < 4; ++f) {
                    std::vector<OCP_ULL> ids = {
                            p[tetFaces[f][0]],
                            p[tetFaces[f][1]],
                            p[tetFaces[f][2]]
                    };
                    std::sort(ids.begin(), ids.end());
                    faceElems[ids].push_back(ei);
                }
            }
            else if(p.size() == 8) {
                // 六面体：6 个四边形面
                static const int hexFaces[6][4] = {
                        {0,1,2,3},{4,5,6,7},
                        {0,1,5,4},{1,2,6,5},
                        {2,3,7,6},{3,0,4,7}
                };
                for(int f = 0; f < 6; ++f) {
                    std::vector<OCP_ULL> ids(4);
                    for(int k = 0; k < 4; ++k)
                        ids[k] = p[hexFaces[f][k]];
                    std::sort(ids.begin(), ids.end());
                    faceElems[ids].push_back(ei);
                }
            }
        }

        // 将所有面（不论是内部面还是边界面）插入到 faces 中
        for(auto& kv : faceElems) {
            const auto& nodes      = kv.first;   // 该面的节点索引列表
            const auto& elemList   = kv.second;  // 持有此面的体元索引列表

            // 以第一个体元的物理属性作为此面的物理属性
            const auto& firstElem = elements[elemList[0]];

            Face f(nodes,
                   firstElem.physical,  // 面的物理名称
                   firstElem.phyIndex,  // 物理索引（OCP_Polygon 的索引）
                    /*tag*/ elemList[0]  // 我们这里用第一个体元的索引作为 tag
            );
            // 把所有共享此面的体元索引都记下来
            for(auto ei : elemList)
                f.faceIndex.push_back(ei);

            faces.insert(f);
        }
    }
}


void GMSHGrid::Setup()
{
    if (dimen == 2) {
        CalAreaCenter2D();
        SetupConnAreaAndBoundary2D();
    }

    if (dimen == 3) {
        CalAreaCenter3D();
        SetupConnAreaAndBoundary3D();
    }
}


void GMSHGrid::CalAreaCenter2D()
{
    for (auto& e : elements) {
        e.CalCenter(points);
        e.CalArea(points);
    }
}

void GMSHGrid::CalAreaCenter3D()
{
    for (auto& e : elements) {
        e.CalCenter(points);
        e.CalArea3D(points);
    }
}

void GMSHGrid::SetupConnAreaAndBoundary2D()
{
    for (const auto& e : edges) {

        if (e.faceIndex.size() == 2) {
            // boundary
            OCP_Polygon& element = elements[e.faceIndex[0]];
            element.boundary   += (e.physical) + " & ";
            element.boundIndex += e.phyIndex;
            // Calculate effective area
            const OCP_ULL   bId       = element.p[e.faceIndex[1]];
            const OCP_ULL   eId       = element.p[(e.faceIndex[1] + 1) % (element.p.size())];
            const Point3D&& edgeNode0 = Point3D(&points[3 * bId]);
            const Point3D&& edgeNode1 = Point3D(&points[3 * eId]);
            const Point3D&& edgeNormal{ -(edgeNode0 - edgeNode1).y, (edgeNode0 - edgeNode1).x, 0 };
            const Point3D&& center2edge = 0.5 * (edgeNode0 + edgeNode1) - element.center;
            e.area.push_back(fabs((center2edge * edgeNormal) / sqrt(center2edge * center2edge)));

            element.boundArea = e.area[0];
        }
        else {
            // internal edge
            for (USI i = 0; i < 2; i++) {
                const OCP_Polygon& element = elements[e.faceIndex[2 * i]];
                // Calculate effective area
                const OCP_ULL   bId       = element.p[e.faceIndex[2 * i + 1]];
                const OCP_ULL   eId       = element.p[(e.faceIndex[2 * i + 1] + 1) % (element.p.size())];
                const Point3D&& edgeNode0 = Point3D(&points[3 * bId]);
                const Point3D&& edgeNode1 = Point3D(&points[3 * eId]);
                const Point3D&& edgeNormal{ -(edgeNode0 - edgeNode1).y, (edgeNode0 - edgeNode1).x, 0 };
                const Point3D&& center2edge = 0.5 * (edgeNode0 + edgeNode1) - element.center;
                e.area.push_back(fabs((center2edge * edgeNormal) / sqrt(center2edge * center2edge)));
            }
        }
    }
}

void GMSHGrid::SetupConnAreaAndBoundary3D()
{
    for (const auto& face : faces) {
        // Find the corresponding element (cell) this face belongs to
        // face.tag is the element tag (cell) which this boundary face belongs to
        OCP_Polygon& element = elements[face.tag];

        element.boundary += face.physical + " & ";
        element.boundIndex += face.phyIndex;

        Point3D faceCenter;
        for (auto pid : face.nodes)
            faceCenter += Point3D(&points[3 * pid]);
        faceCenter *= 1.0 / face.nodes.size();

        Point3D center2face = faceCenter - element.center;

        // Dummy normal (should use true face normal if available)
        Point3D faceNormal = center2face;

        OCP_DBL projected = fabs((center2face * faceNormal) / sqrt(center2face * center2face));
        element.boundArea = projected;
    }
}


void GMSHGrid::PrintElementPoint(const OCP_ULL& n) const
{
    for (const auto& p : elements[n].p) {
        for (USI i = 0; i < 3; i++) {
            cout << points[p * 3 + i] << "   ";
        }
        cout << endl;
    }
}


void GMSHGrid::InputProperty(ifstream& ifs)
{
    if (elements.empty()) {
        OCP_ABORT("INPUT KEYWORD GMSH FIRST!");
    }

    // Note that the order of physical must be consistent with the ones in *.geo
    USI i = 0;
    while (!ifs.eof()) {

        string fbuf;
        getline(ifs, fbuf);

        if (fbuf == "GMSHPROEND") break;

        vector<string> vbuf;

        if (i < facies.size()) {
            if (fbuf == facies[i].name) {
                while (true) {
                    ReadLine(ifs, vbuf);
                    if (vbuf[0] == "END") break;

                    if (vbuf[0] == "*PORO")      facies[i].poro = stod(vbuf[1]);
                    else if (vbuf[0] == "*PERM") {
                        facies[i].kx = stod(vbuf[1]);
                        facies[i].ky = stod(vbuf[1]);
                        facies[i].kz = stod(vbuf[1]);
                    }
                    else if (vbuf[0] == "*PERMX") {
                        facies[i].kx = stod(vbuf[1]);
                    }
                    else if (vbuf[0] == "*PERMY") {
                        facies[i].ky = stod(vbuf[1]);
                    }
                    else if (vbuf[0] == "*PERMZ") {
                        facies[i].kz = stod(vbuf[1]);
                    }
                }
                i++;
            }
        }

        // other params
        vbuf.clear();
        istringstream tmp(fbuf);
        while (tmp >> fbuf)  vbuf.push_back(fbuf);

        if (dimen == 2) {
            if (!vbuf.empty() && vbuf[0] == "THICKNESS") {
                thickness = stod(vbuf[1]);
            }
        }
    }


    if (i != facies.size()) {
        OCP_ABORT("Order of Physical body must be consistent with the ones in *.geo!");
    }
}


#endif

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Sep/10/2023      Create file                          */
/*----------------------------------------------------------------------------*/