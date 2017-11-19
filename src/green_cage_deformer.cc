#include <maya/MVector.h>
#include <maya/MFnTypedAttribute.h>
#include <maya/MFnNumericAttribute.h>
#include <maya/MFnPlugin.h>
#include <maya/MPointArray.h>
#include <maya/MFnMesh.h>
#include <maya/MMatrix.h>
#include <maya/MItGeometry.h>
#include <maya/MDagPath.h>

#include <numeric>
#include <cmath>
#include <vector>
#include <array>

#include "tbb/tick_count.h"
#include "tbb/parallel_for.h"
#include "tbb/task_scheduler_init.h"

#include <Eigen/Dense>
#include <Eigen/StdVector>

#include "green_cage_deformer.h"


const MTypeId GreenCageDeformer::s_node_id(0xBBBB9255);
const MString GreenCageDeformer::s_node_name("green_cage_deformer");
MObject GreenCageDeformer::s_cage_mesh;


namespace {
	const Eigen::Vector3d NULL_VECTOR(0.0, 0.0, 0);
	constexpr auto TOLERANCE = 1e-6;
}


namespace MathConstants {
	constexpr auto PI = 3.14159265358979323846;
	constexpr auto ONE_OVER_FOUR_PI = 0.07957747154594767;
	constexpr auto SQRT8 = 2.828427124746190097603;
};


const double GCTriInt(const Eigen::Vector3d& p, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2)
{
	const auto v2_v1 = v2 - v1;
	const auto p_v1 = p - v1;
	const auto v1_p = v1 - p;
	const auto v2_p = v2 - p;

	const auto alpha = std::acos(std::min(std::max(((v2_v1).dot((p_v1))) / ((v2_v1).norm() * (p_v1).norm()), -1.0), 1.0));
	if (abs(alpha - MathConstants::PI) < TOLERANCE || abs(alpha) < TOLERANCE) {
		return 0.0;
	}

	const auto beta = std::acos(std::min(std::max((((v1_p).dot((v2_p)))) / ((v1_p).norm() * (v2_p).norm()), -1.0), 1.0));
	const auto lambda = p_v1.squaredNorm() * std::sin(alpha) * std::sin(alpha);
	const auto c = p.squaredNorm();

	const auto sqrt_c = sqrt(c);
	const auto sqrt_lambda = sqrt(lambda);

	const std::array<double, 2> theta{ MathConstants::PI - alpha, MathConstants::PI - alpha - beta };
	std::array<double, 2> I;

	for (size_t i = 0; i < 2; ++i)
	{
		const auto S = std::sin(theta[i]);
		const auto C = std::cos(theta[i]);
		const auto sign = S < 0 ? -1.0 : 1.0;

		const auto SS = S*S;
		const auto half_sign = (-sign * 0.5);
		const auto tan_part = (2 * sqrt_c * std::atan2((sqrt_c * C), sqrt(lambda + (SS * c))));
		const auto log_part = log(((2 * sqrt_lambda * SS) / std::pow(1.0 - C, 2)) * (1.0 - ((2 * c * C) / ((c*(1 + C) + lambda + sqrt((lambda * lambda) + (lambda * c * SS)))))));

		I[i] = half_sign * (tan_part + (sqrt_lambda * log_part));
	}

	return -MathConstants::ONE_OVER_FOUR_PI * abs(I[0] - I[1] - sqrt_c * beta);
}


void GCDeformData::setup_cage() {
	MStatus stat;
	MFnMesh fn_mesh(m_cage_object);
	MIntArray m_tri_counts;
	fn_mesh.getTriangles(m_tri_counts, m_tri_verts);

	m_nr_of_tris = 0;
	for (unsigned int i = 0; i < m_tri_counts.length(); ++i)
	{
		m_nr_of_tris += m_tri_counts[i];
	}


	MPointArray cage_points;
	fn_mesh.getPoints(cage_points);
	m_nr_of_verts = fn_mesh.numVertices();
	m_vertices.resize(m_nr_of_verts);

	for (unsigned int i = 0; i < cage_points.length(); ++i) {
		const auto point = cage_points[i];
		m_vertices[i] = Eigen::Vector3d(point.x, point.y, point.z);
	}

	m_scale_factor.resize(m_nr_of_tris);
	m_tri_areas.resize(m_nr_of_tris);
	m_tri_normals.resize(m_nr_of_tris);
	m_tri_points.resize(m_nr_of_tris);
	m_tri_edges.resize(m_nr_of_tris);

	for (unsigned int j = 0; j < m_nr_of_tris; ++j) {
		const auto tri_idx = j * 3;
		const auto v0 = m_vertices[m_tri_verts[tri_idx + 0]];
		const auto v1 = m_vertices[m_tri_verts[tri_idx + 1]];
		const auto v2 = m_vertices[m_tri_verts[tri_idx + 2]];
		m_tri_points[j][0] = v0;
		m_tri_points[j][1] = v1;
		m_tri_points[j][2] = v2;
		// ccw order
		const Eigen::Vector3d edge1{ v1 - v0 };
		const Eigen::Vector3d edge2{ v2 - v1 };
		const Eigen::Vector3d cross = edge1.cross(edge2);
		m_tri_areas[j] = cross.norm()*0.5;
		m_tri_normals[j] = cross.normalized();
		m_tri_edges[j] = std::make_pair(edge1, edge2);
	}

}

// variable names based on pseudocode
void GCDeformData::calc_green_coords(MItGeometry& iter, const MMatrix& mat)
{
	MPointArray target_vertices;
	iter.allPositions(target_vertices);
	m_nr_of_targetverts = iter.exactCount();
	m_psi.resize(m_nr_of_tris*m_nr_of_targetverts, 0.0);
	m_phi.resize(m_nr_of_verts*m_nr_of_targetverts, 0.0);

	tbb::parallel_for(tbb::blocked_range<unsigned int>(0, m_nr_of_targetverts), [&](const tbb::blocked_range<unsigned int>& r) {
		for (unsigned int idx = r.begin(); idx < r.end(); ++idx) {
			const auto maya_vec = target_vertices[idx] * mat;
			const Eigen::Vector3d pvec(maya_vec.x, maya_vec.y, maya_vec.z);
			Eigen::Vector3d s;
			Eigen::Vector3d I;
			std::array<double, 3> II;
			std::array<Eigen::Vector3d, 3> N;
			for (unsigned int i = 0; i < m_nr_of_tris; ++i) {
				const auto nrm = m_tri_normals[i];
				std::array<Eigen::Vector3d, 3> vj;
				for (size_t l = 0; l < 3; ++l) {
					vj[l] = m_tri_points[i][l] - pvec;
				}
				const auto p = nrm * (vj[0].dot(nrm));
				for (size_t k = 0; k < 3; ++k) {
					const auto v0 = vj[k];
					const auto v1 = vj[(k + 1) % 3];
					const auto lg = ((v0 - p).cross((v1 - p))).dot(nrm);
					s[k] = lg < 0 ? -1.0 : 1.0;
					// eta always zero(?)
					I[k] = GCTriInt(p, v0, v1);
					II[k] = GCTriInt(NULL_VECTOR, v1, v0);
					N[k] = (v1.cross(v0)).normalized();
				}

				const auto I_ = -abs(s.dot(I));

				m_psi[(idx*m_nr_of_tris) + i] = -I_;

				Eigen::Vector3d w = I_ * nrm;
				for (int k = 0; k < 3; k++)
					w += (II[k] * N[k]);

				if (w.norm() > DBL_EPSILON) {
					for (unsigned int l = 0; l < 3; l++)
					{
						const auto l1 = (l + 1) % 3;
						const auto a1 = N[l1];
						const auto a2 = w;
						const auto a3 = vj[l];
						const auto val = ((a1.dot(a2)) / (a1.dot(a3)));
						m_phi[m_nr_of_verts*idx + m_tri_verts[(i * 3) + l]] += val;
					}
				}
			}
		}
	});
}


const MPoint GCDeformData::calc_pos(unsigned int idx) const
{
	Eigen::Vector3d pnt(0.0, 0.0, 0.0);
	const auto current_tri = idx*m_nr_of_tris;
	const auto current_vert = idx*m_nr_of_verts;

	for (unsigned int i = 0; i < m_nr_of_verts; ++i) {
		pnt += m_phi[current_vert + i] * m_vertices[i];
	}

	for (unsigned int i = 0; i < m_nr_of_tris; ++i) {

		pnt += m_psi[current_tri + i] * m_scale_factor[i] * m_tri_normals[i];
	}

	return MPoint(pnt.x(), pnt.y(), pnt.z());
}


void GCDeformData::calc_scale_factor()
{

	for (unsigned int i = 0; i < m_nr_of_tris; ++i) {
		const auto rest_edges = m_original_edges[i];
		const auto current_edges = m_tri_edges[i];

		const auto u0 = rest_edges.first;
		const auto v0 = rest_edges.second;

		const auto u1 = current_edges.first;
		const auto v1 = current_edges.second;

		m_scale_factor[i] = sqrt((u1.squaredNorm()) * (v0.squaredNorm()) - 2.0 * (u1.dot(v1)) * (u0.dot(v0)) + (v1.squaredNorm()) * (u0.squaredNorm())) / (MathConstants::SQRT8 * m_original_areas[i]);

	}
}


void GCDeformData::save_original_cage_data()
{
	m_original_areas.reserve(m_nr_of_tris);
	m_original_areas = m_tri_areas;

	m_original_edges.reserve(m_nr_of_tris);
	m_original_edges = m_tri_edges;
}


MStatus GreenCageDeformer::deform(MDataBlock& block, MItGeometry& iter, const MMatrix& mat, unsigned int multiIndex)
{
	MStatus stat;

	// not yet supported
	if (multiIndex > 0)
		return stat;

	if (!m_cage_connected)
		return stat;

	const auto envelope_val = block.inputValue(envelope).asDouble();
	if (envelope_val <= 0.0)
		return stat;

	const auto h_cage = block.inputValue(s_cage_mesh);

	// since we are getting the worldMesh lets not fiddle around with local to world transformations
	const auto cage_object = h_cage.asMeshTransformed();

	p_cage_data->m_cage_object = cage_object;
	p_cage_data->setup_cage();

	if (!m_green_calculated) {
		MArrayDataHandle h_input_array = block.outputArrayValue(input, &stat);
		// atm multiIndex will be only 0
		stat = h_input_array.jumpToElement(multiIndex);
		p_cage_data->calc_green_coords(iter, mat);
		p_cage_data->save_original_cage_data();
		m_green_calculated = true;
		// TODO: write out green coords w/ original cage to be able to reload scene without calculating(also correctly) everytime we open it
	}

	p_cage_data->calc_scale_factor();
	MPointArray new_arr;
	new_arr.setLength(p_cage_data->m_nr_of_targetverts);

	const auto inverse_mat = mat.inverse();
	tbb::parallel_for(tbb::blocked_range<unsigned int>(0, p_cage_data->m_nr_of_targetverts), [&](const tbb::blocked_range<unsigned int>& r) {
		for (unsigned int idx = r.begin(); idx < r.end(); ++idx) {
			new_arr[idx] = p_cage_data->calc_pos(idx) * inverse_mat;
		}
	});

	iter.setAllPositions(new_arr);
	return stat;
}

void GreenCageDeformer::postConstructor()
{
	m_cage_connected = false;
	m_green_calculated = false;
}

MStatus GreenCageDeformer::connectionMade(const MPlug & plug, const MPlug & otherPlug, bool asSrc)
{
	if (!asSrc && plug == s_cage_mesh) {		
		p_cage_data.reset(new GCDeformData());
		m_cage_connected = true;
	}

	return MPxDeformerNode::connectionMade(plug, otherPlug, asSrc);
}

MStatus GreenCageDeformer::connectionBroken(const MPlug & plug, const MPlug & otherPlug, bool asSrc)
{
	if (!asSrc && plug == s_cage_mesh) {
		p_cage_data.reset();
		m_cage_connected = false;
		m_green_calculated = false;
	}

	return MPxDeformerNode::connectionBroken(plug, otherPlug, asSrc);
}


void* GreenCageDeformer::creator()
{
	return new GreenCageDeformer();
}


MStatus GreenCageDeformer::initialize()
{
	MStatus stat;
	MFnTypedAttribute tAttr;

	s_cage_mesh = tAttr.create("cageMesh", "cm", MFnData::kMesh, &stat);
	CHECK_MSTATUS_AND_RETURN_IT(stat);
	stat = addAttribute(s_cage_mesh);
	CHECK_MSTATUS_AND_RETURN_IT(stat);
	stat = attributeAffects(s_cage_mesh, outputGeom);

	return stat;
}


MStatus initializePlugin(MObject obj)
{
	MStatus stat;

	MFnPlugin plugin(obj, "Balazs Pataki", "1.0", "Any");
	stat = plugin.registerNode(GreenCageDeformer::s_node_name, GreenCageDeformer::s_node_id, &GreenCageDeformer::creator, &GreenCageDeformer::initialize, MPxNode::kDeformerNode);
	CHECK_MSTATUS_AND_RETURN_IT(stat);

	return stat;
}


MStatus uninitializePlugin(MObject obj)
{
	MStatus stat;

	MFnPlugin plugin(obj);
	stat = plugin.deregisterNode(GreenCageDeformer::s_node_id);
	CHECK_MSTATUS_AND_RETURN_IT(stat);

	return stat;
}
