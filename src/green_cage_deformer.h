#pragma once

#include <maya/MPxDeformerNode.h>
#include <memory>
#include <Eigen/Dense>


typedef std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> Vec3Array;


struct GCDeformData
{
	GCDeformData() = default;
	~GCDeformData() = default;

	void setup_cage();
	void calc_green_coords(MItGeometry& iter, const MMatrix& mat);
	const MPoint calc_pos(unsigned int idx) const;
	void calc_scale_factor();
	void save_original_cage_data();

	MObject m_cage_object;
	MIntArray m_tri_verts;

	unsigned int m_nr_of_verts;
	unsigned int m_nr_of_tris;
	unsigned int m_nr_of_targetverts;

	std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> m_tri_edges;
	std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> m_original_edges;
	std::vector<std::array<Eigen::Vector3d, 3>> m_tri_points;
	Vec3Array m_vertices;
	Vec3Array m_tri_normals;
	std::vector<int> m_triangle_indices;
	std::vector<double> m_tri_areas;
	std::vector<double> m_original_areas;
	std::vector<double> m_psi;
	std::vector<double> m_phi;
	std::vector<double> m_scale_factor;

};


class GreenCageDeformer : MPxDeformerNode
{
public:
	static void* creator();
	GreenCageDeformer() = default;
	~GreenCageDeformer() override = default;

	MStatus deform(MDataBlock& block, MItGeometry& iter, const MMatrix& mat, unsigned int multiIndex) override;
	void postConstructor();
	MStatus connectionMade(const MPlug &plug, const MPlug &otherPlug, bool asSrc) override;
	MStatus connectionBroken(const MPlug &plug, const MPlug &otherPlug, bool asSrc) override;
	static MStatus initialize();

	static const MString s_node_name;
	static const MTypeId s_node_id;

	static MObject s_cage_mesh;
	// laziness to cleanup atm :)
	std::unique_ptr<GCDeformData> p_cage_data;
	bool m_cage_connected;
	bool m_green_calculated;
};
