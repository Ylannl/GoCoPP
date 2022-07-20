#pragma once
#define CGAL_EIGEN3_ENABLED
#include "defs.h"
#include "defs_cgal_ui.h"


#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Point_with_normal_3.h>
#include <CGAL/Shape_detection_3.h>
#include <CGAL/property_map.h>
#include <CGAL/Surface_mesh.h>
#include <fstream>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>


#include <CGAL/Polygon_mesh_processing/stitch_borders.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/compute_normal.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/compute_average_spacing.h>
#include <CGAL/pca_estimate_normals.h>
#include <CGAL/mst_orient_normals.h>
#include <CGAL/property_map.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>

#ifdef CGAL_LINKED_WITH_TBB
typedef CGAL::Parallel_tag Concurrency_tag;
#else
typedef CGAL::Sequential_tag Concurrency_tag;
#endif



typedef std::pair<Inexact_Point_3, Inexact_Vector_3> Point_with_normal;

typedef std::vector<Point_with_normal> Pwn_vector;



typedef CGAL::First_of_pair_property_map<Point_with_normal> Point_map;
typedef CGAL::Second_of_pair_property_map<Point_with_normal> Normal_map;

typedef CGAL::Shape_detection_3::Shape_detection_traits
<CGAL::Exact_predicates_inexact_constructions_kernel, Pwn_vector, Point_map, Normal_map> Traits;
typedef CGAL::Exact_predicates_inexact_constructions_kernel::Point_3 Inexact_Point_3;




class Shape_Detector
{

	struct Add_Comparator {
		bool operator()(std::pair<int, double> p1, std::pair<int, double> p2) {

			return p1.second > p2.second;
		}
	};
	struct Add_Comparator_normal {
		bool operator()(std::pair<int, double> p1, std::pair<int, double> p2) {

			return p1.second < p2.second;
		}
	};

	struct Remove_Comparator {
		bool operator()(std::pair<int, double> p1, std::pair<int, double> p2) {

			return p1.second < p2.second;
		}
	};
	struct Remove_Comparator_normal {
		bool operator()(std::pair<int, double> p1, std::pair<int, double> p2) {

			return p1.second > p2.second;
		}
	};

	struct Weight_Comparator_with_energy {
		bool operator()(std::pair<std::pair<int, std::vector<int>>, std::vector<double>> p1, std::pair<std::pair<int, std::vector<int>>, std::vector<double>> p2) {

			return p1.second[0] > p2.second[0];
		}
	};


public:
	Shape_Detector(const std::string & _filename);

	~Shape_Detector();

	bool load_points();
	

protected:
	bool load_ply();




protected:
	void set_extrema();

public:



	void set_lambda_r(double db);
	void set_weight_m(int wm);
	void set_max_steps(int wm);
	void set_lambda_c(double db);
	void set_detection_parameters(int _min_points,  int _knn, double _normal_threshold);
	void set_epsilon(double _rg_epsilon);
	void detect_shapes();

protected:
	void compute_average_spacing_and_k_nearest_neighbors();
	void do_region_growing();



public:



	void test_connected_primitives();


	void planar_shape_detection_L1();
	void planar_shape_detection_hybrid();





	void save_infor_to_txt_right_folder(std::string name);

	void planar_shape_detection_l2();



	void transfer_operator_normal();
	void transfer_operator_normal_for_hybrid();

	void transfer_operator_l2();

	double energy_changed(double dis, double numb);
	double get_final_energy(std::string meseurment, int mm);
	double get_original_energy(std::string meseurment, int mm);

	double energy_changed_merge(double dis, double numb, int nmoves);
	double energy_changed_normal(double dis, double numb);
	double energy_changed_normal_merge(double dis, double numb, int nmoves);
	double energy_changed_second_normal(double dis, double numb);

	double energy_changed_second(double dis, double numb);

	bool separate_two_out(int id, std::vector<int> & max_list, std::vector<int> & min_list, double & dif);
	bool separate_two_out_normal(int id, std::vector<int> & max_list, std::vector<int> & min_list, double & dif);
	bool convert_form_split_normal(int i, std::vector<std::vector<int>>& max_list_vector, std::vector<std::vector<int>>& min_list_vector, std::pair<std::pair<int, std::vector<int>>, std::vector<double>> & one_element);
	bool convert_form_insert_normal(int i, std::pair<std::pair<int, std::vector<int>>, std::vector<double>>& one_element);
	bool convert_form_exclude_normal(int i, std::pair<std::pair<int, std::vector<int>>, std::vector<double>>& one_element);
	bool update_good_points_normal(int id_shape, std::pair<std::pair<int, std::vector<int>>, std::vector<double>>& one_element);
	bool update_bad_points_normal(int id_shape, std::pair<std::pair<int, std::vector<int>>, std::vector<double>>& one_element);
	bool convert_form_merge_normal(int i, int j, std::pair<std::pair<int, std::vector<int>>, std::vector<double>>& one_element);
	bool update_good_points(int id_shape, std::pair<std::pair<int, std::vector<int>>, std::vector<double>>& one_element);
	bool update_bad_points(int id_shape, std::pair<std::pair<int, std::vector<int>>, std::vector<double>>& one_element);
	bool test_if_connected(int i, int j);
	void local_operators();
	void local_operators_normal();

	bool convert_form_merge_l2(int i, int j, std::pair<std::pair<int, std::vector<int>>, std::vector<double>>& one_element);
	bool convert_form_split_l2(int i, std::vector<std::vector<int>>& max_list_vector, std::vector<std::vector<int>>& min_list_vector, std::pair<std::pair<int, std::vector<int>>, std::vector<double>> & one_element);
	bool convert_form_exclude_l2(int i, std::pair<std::pair<int, std::vector<int>>, std::vector<double>>& one_element);
	bool convert_form_insert_l2(int i, std::pair<std::pair<int, std::vector<int>>, std::vector<double>>& one_element);

	double merge_distance_changed_normal_with_epsilon(int i, int j, double & dif, std::vector<int> & move_ids);


	double merge_distance_changed_with_epsilon(int i, int j, double & dif, std::vector<int> & move_ids);


protected:
	void detect_planes();

	bool inliers_arent_aligned(const std::vector<int> & indices);

	void get_coverage_and_mean_error();
	void get_coverage_and_mean_error_pure();

public:

	void get_last_information();
	void show_result(double t);
	void set_primitives_simple();


public:

	void to_vg();
	void save_alpha_shapes();
	void save_convex_hull();

	void set_constraint(bool cc);

	void get_distance_diviation();


	void get_distance_diviation_show_normal_info(double t);
	void get_distance_diviation_show_merge_info(double t);

	double add_changed_error(int id, std::vector<int> point_ids);
	double add_changed_error_normal(int id, std::vector<int> point_ids);

	void get_bad_points();
	void get_good_points();
	void get_bad_points_normal();
	void get_good_points_normal();
	double remove_changed_error(int id, std::vector<int> point_ids);
	double remove_changed_error_normal(int id, std::vector<int> point_ids);

	double get_bbox_diagonal() const;
	double get_average_spacing() const;


	double get_epsilon() const;
	int get_knn() const;

	double get_normal_threshold() const;

	int get_neighborhood_algorithm() const;
	double get_neighborhood_distance() const;



protected:
	std::string path_point_cloud;
	std::string path_point_cloud_basename;
	std::string path_point_cloud_extension;

	Pwn_vector points;


protected:
	std::vector<std::vector<int> > spherical_neighborhood;
	std::vector<Inexact_Point_3> planes_centroids;
	
protected:
	std::string refinement_type;
	int ori_primitives_number;
	double ori_coverage;
	double ori_mean_error;
	int ori_inliers_number;
	int last_primitives_number;
	double last_normal_deviation;
	double last_coverage;
	double last_mean_error;
	int primitives_number;
	double coverage;
	double x_min, x_max, y_min, y_max, z_min, z_max;
	double bbox_diagonal;
	double average_spacing;
	bool spacing_is_known;
	
	int t_l;
	double epsilon;
	int knn;
	bool should_compute_knn;
	double normal_threshold;
	int neighborhood_algorithm;
	double neighborhood_distance;
	bool should_compute_neighborhood;
	std::vector<std::set<int> > neighbors;

	bool if_constraint;


public:


	int number_iterations;



	std::vector<int> bad_points;
	std::vector<std::pair<int,int>> good_points;
	std::vector<std::pair<std::pair<int, std::vector<int>>,std::vector<double>>> good_points_shape;
	std::vector<std::pair<std::pair<int, std::vector<int>>, std::vector<double>>> bad_points_shape;


	std::vector<int> points_if_added;
	std::vector<int> points_changed;

	int min_points;

	std::vector<Inexact_Plane> planes_0;
	std::vector<Inexact_Plane> planes_1; 
	std::vector<Inexact_Plane> planes_2;
	std::vector<Inexact_Plane> planes_3;
	std::vector<std::vector<int>> primitive_connection;

	std::vector<CGAL::Color> planes_to_colors;

	std::vector<std::vector<int> > planes_to_inliers;
	std::vector<int> region_type;
	std::vector<int> region_type_m;


	std::vector<int> inliers_to_planes;



	int non_coplanar_planes;
	double mean_error;
	double all_error;
	double ori_all_error;

	double mean_normal_diviation;
	double mean_distance_diaviation;

	double interval_all;
	double all_distance_diaviation;
	double all_normal_diaviation;

	size_t number_of_assigned_points;
	double size_current_primitives;
	double mean_distance_current;
	double mean_normal_current;

	double ori_mean_normal_diviation;
	double number_inlier_before_opers;

	int t_m;
	int t_merge;
	int t_split;
	int t_insert;
	int t_exlude;
	int all_t_transfer;
	int all_t_merge;
	int all_t_split;
	int all_t_insert;
	int all_t_exlude;

	std::vector<std::vector<Inexact_Point_3> > alpha_shapes_pts;
	std::vector<CGAL::Color> alpha_shapes_colors;


	std::vector<std::vector<Inexact_Point_3> > convex_hulls_pts;
	std::vector<CGAL::Color> convex_hulls_colors;

	std::vector<std::pair<int, int> > convex_hulls_seeds;


	double lambda_r;
	int weight_mode;
	double lambda_c;

	int old_size_current_primitives;

	double old_mean_distance_diaviation;

	double old_mean_normal_diaviation;

	double old_coverage;
	int number_of_insert_exclude;

	std::vector<int> if_insert_candidate;
	double mean_distance_diaviation_before_transfer;
	double mean_normal_diviation_befor_transfer;
	int stop_iterations;

};

