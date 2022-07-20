#include "shape_detector.h"

#include <CGAL/IO/read_ply_points.h>


#include <CGAL/convex_hull_2.h>
#include <CGAL/compute_average_spacing.h>

#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <random>

Shape_Detector::Shape_Detector(const std::string & _filename)
{
	number_of_insert_exclude = 10;
	number_iterations = 0;
	
	refinement_type = "unknow";
	interval_all = 0;
	path_point_cloud = _filename;
	path_point_cloud_basename = boost::filesystem::path(path_point_cloud).stem().string();
	path_point_cloud_extension = boost::filesystem::extension(path_point_cloud);

	lambda_r = 1.0;
	weight_mode = 0;
	lambda_c = 1.0;
	knn = -1;
	should_compute_knn = true;
	should_compute_neighborhood = true;






	mean_error = 0;
	mean_normal_diviation = 0;
	//number_distance = 0;
	//number_normal = 0;
	mean_distance_diaviation = 0;

	t_exlude = 0;
	t_insert = 0;
	t_merge = 0;
	t_insert = 0;
	t_split = 0;
	all_t_transfer = 0;
	all_t_exlude = 0;
	all_t_insert = 0;
	all_t_split = 0;
	all_t_merge = 0;

	if_constraint = false;
	stop_iterations = 7;
	
}


Shape_Detector::~Shape_Detector()
{
}

void Shape_Detector::set_lambda_r(double db) {
	lambda_r = db;

}


void Shape_Detector::set_lambda_c(double db) {
	lambda_c = db;

}
void Shape_Detector::set_weight_m(int wm) {
	weight_mode = wm;
	

}
void Shape_Detector::set_max_steps(int wm) {
	stop_iterations = wm;


}

void Shape_Detector::to_vg()
{
	int X = 1, Y = 1, Z = 1;
	Pwn_vector T = points;
	// Step 1.
	// For each point, we find in which subdivision it is located.

	double x_inf = FLT_MAX, x_sup = -FLT_MAX;
	double y_inf = FLT_MAX, y_sup = -FLT_MAX;
	double z_inf = FLT_MAX, z_sup = -FLT_MAX;

	for (size_t i = 0; i < T.size(); ++i) {
		const Inexact_Point_3 & pt = T[i].first;
		const double &x = pt.x(), &y = pt.y(), &z = pt.z();

		if (x < x_inf) x_inf = x;
		if (x > x_sup) x_sup = x;
		if (y < y_inf) y_inf = y;
		if (y > y_sup) y_sup = y;
		if (z < z_inf) z_inf = z;
		if (z > z_sup) z_sup = z;
	}

	double dx = x_sup - x_inf;
	double dy = y_sup - y_inf;
	double dz = z_sup - z_inf;

	typedef std::tuple<int, int, int> Triplet;
	std::vector<Triplet> vertices_locations;
	vertices_locations.reserve(T.size());

	for (size_t i = 0; i < T.size(); ++i) {
		const Inexact_Point_3 & pt = T[i].first;

		int a = jclamp(0, int(X * (pt.x() - x_inf) / dx), X - 1);
		int b = jclamp(0, int(Y * (pt.y() - y_inf) / dy), Y - 1);
		int c = jclamp(0, int(Z * (pt.z() - z_inf) / dz), Z - 1);
		vertices_locations.push_back(Triplet(a, b, c));
	}

	// Step 2.
	// For each planar shape we need to know in which subdivisions it can be found

	std::vector< std::set<Triplet> > shapes_locations;
	shapes_locations.reserve(planes_to_inliers.size());

	for (size_t i = 0; i < planes_to_inliers.size(); ++i) {
		std::set<Triplet> L;
		for (size_t j = 0; j < planes_to_inliers[i].size(); ++j) {
			int ind_ij = planes_to_inliers[i][j];
			L.insert(vertices_locations[ind_ij]);
		}
		shapes_locations.push_back(L);
	}

	// Step 3.
	// Prints a .vg file for each subdivision

	for (int i = 0; i < X; ++i) {
		for (int j = 0; j < Y; ++j) {
			for (int k = 0; k < Z; ++k) {
				Triplet tr(i, j, k);

				// Step 3.1
				// Collects points in T_curr

				std::vector<int> local_indices(T.size(), -1);

				Pwn_vector T_curr;
				for (size_t p = 0; p < T.size(); ++p) {
					if (vertices_locations[p] == tr) {
						local_indices[p] = T_curr.size();
						T_curr.push_back(T[p]);
					}
				}

				// Step 3.2
				// Gets planes associated to the current triplet

				std::vector<Inexact_Plane> S_eq_curr;
				std::vector<std::vector<int> > S_inl_curr;

				for (size_t p = 0; p < planes_to_inliers.size(); ++p) {
					if (shapes_locations[p].find(tr) == shapes_locations[p].end()) continue;

					S_eq_curr.push_back(planes_2[p]);

					std::vector<int> I_p;
					for (size_t q = 0; q < planes_to_inliers[p].size(); ++q) {
						int ind = planes_to_inliers[p][q];
						if (local_indices[ind] != -1) {
							I_p.push_back(local_indices[ind]);
						}
					}
					S_inl_curr.push_back(I_p);
				}

				std::cout << "[" << i << ", " << j << ", " << k << "] : " << S_eq_curr.size() << std::endl;

				// Step 3.3
				// Prints a vg file

			
				std::string filename = "./" + path_point_cloud_basename + "_planar_primitives_detection.vg";
				std::ofstream stream(filename, std::ios::out);

				// Step 3.3.1
				// Provides a description of the point cloud

				size_t N = T_curr.size();
				stream << "num_points: " << N << std::endl;
				for (size_t i = 0; i < N; ++i) {
					const Inexact_Point_3 & pt = T_curr[i].first;
					stream << pt.x() << " " << pt.y() << " " << pt.z() << std::endl;
				}
				stream << "num_colors: " << N << std::endl;
				for (size_t i = 0; i < N; ++i) {
					stream << "128 128 128" << std::endl;
				}
				stream << "num_normals: " << N << std::endl;
				for (size_t i = 0; i < N; ++i) {
					const Inexact_Vector_3 & n = T_curr[i].second;
					stream << n.x() << " " << n.y() << " " << n.z() << std::endl;
				}

				// Step 3.3.2
				// Provides a description of the detected planes.

				size_t M = S_eq_curr.size();
				stream << "num_groups: " << M << std::endl;

				for (size_t i = 0; i < S_eq_curr.size(); ++i) {
					stream << "group_type: 0" << std::endl;
					stream << "num_group_parameters: 4" << std::endl;

					const Inexact_Plane & H = S_eq_curr[i];
					stream << "group_parameters: " << H.a() << " " << H.b() << " " << H.c() << " " << H.d() << std::endl;

					CGAL::Color C = planes_to_colors[i];
					const std::vector<int> & shape_assigned_pts = S_inl_curr[i];

					stream << "group_label: unknown" << std::endl;
					stream << "group_color: " << int(C.red()) << " " << int(C.green()) << " " << int(C.blue()) << std::endl;
					stream << "group_num_points: " << shape_assigned_pts.size() << std::endl;

					for (size_t j = 0; j < shape_assigned_pts.size(); ++j) {
						stream << shape_assigned_pts[j] << (j == shape_assigned_pts.size() - 1 ? '\n' : ' ');
					}
					stream << "num_children: 0" << std::endl;
				}

				stream.close();
			}
		}
	}
}




bool Shape_Detector::load_points()
{
	points.clear();

	if (path_point_cloud_extension == ".ply") {
		load_ply();
	} 


	spacing_is_known = false;
	set_extrema();

	return true;
}


bool Shape_Detector::load_ply()
{
	
	std::ifstream streamb(path_point_cloud);
	int line = 0;
	std::string s;
	while (streamb&&line<2) {
		if (!getline(streamb, s)) break;
		line++;
	}


	if(s=="format ascii 1.0"){
		
		std::ifstream stream(path_point_cloud);
		if (!stream || !CGAL::read_ply_points(stream, std::back_inserter(points), CGAL::parameters::point_map(Point_map()).normal_map(Normal_map()))) {
			return false;
		}
	}
	else {
		std::ifstream stream(path_point_cloud, std::ios_base::binary);
		
		
		if (!stream || !CGAL::read_ply_points(stream, std::back_inserter(points), CGAL::parameters::point_map(Point_map()).normal_map(Normal_map()))) {
			return false;
		}
	}
	if (points[0].second == Inexact_Vector_3(0, 0, 0)) {
		std::cout << std::endl;
		std::cout << "The Input PLY file does not contain 'Normals' or contains invalid 'Normals', The 'Normals' will be estimated."<<std::endl;

		CGAL::pca_estimate_normals<Concurrency_tag>
			(points,
				knn, // when using a neighborhood radius, K=0 means no limit on the number of neighbors returns
				CGAL::parameters::point_map(Point_map())
				.normal_map(Normal_map()));
	}
	
	return true;
}


void Shape_Detector::compute_average_spacing_and_k_nearest_neighbors()
{
	spherical_neighborhood.clear();
	spherical_neighborhood.reserve(points.size());

	std::map<Inexact_Point_3, int> map_indice_point;
	std::list<Inexact_Point_3> list_points;

	for (int i = 0; i < points.size(); i++) {
		const Inexact_Point_3 & pt = points[i].first;
		map_indice_point[pt] = i;
		list_points.push_back(pt);
	}

	Tree tree(list_points.begin(), list_points.end());
	if (!spacing_is_known) average_spacing = 0;

	std::cout << "Computing K nearest neighbors ";

	int ten_percent_of_points = points.size() / 10;

	for (int i = 0; i < points.size(); i++) {
		if (i % ten_percent_of_points == 0) std::cout << ". " << std::flush;

		Inexact_Point_3 query = points[i].first;
		Neighbor_search search(tree, query, knn + 1);

		double dist_1st_neighbor = FLT_MAX;
		std::vector<int> index_of_neighbors;
		index_of_neighbors.reserve(knn);
		
		int terms = 0;

		for (Neighbor_search::iterator it = search.begin(); it != search.end(); ++it) {
			std::map<Inexact_Point_3, int>::iterator iter = map_indice_point.begin();
			iter = map_indice_point.find(it->first);
			if (iter != map_indice_point.end() && iter->second != i) {
				index_of_neighbors.push_back(iter->second);
				
				if (!spacing_is_known) {
					double d = sqrt((query - iter->first).squared_length());
					if (d < dist_1st_neighbor) dist_1st_neighbor = d;
				}

				++terms;
			}
		}

		spherical_neighborhood.push_back(index_of_neighbors);
		if (!spacing_is_known && terms > 0) average_spacing += dist_1st_neighbor;
	}

	std::cout << std::endl;

	if (!spacing_is_known) {
		average_spacing /= double(points.size());
		std::cout << "Average spacing : " << average_spacing << std::endl;
		spacing_is_known = true;
	}

	should_compute_knn = false;
}



void Shape_Detector::set_extrema()
{
	for (size_t i = 0 ; i < points.size() ; ++i) {
		const Inexact_Vector_3 & n = points[i].second;
		const double &x = n.x(), &y = n.y(), &z = n.z();
		double l = sqrt(x * x + y * y + z * z);
		points[i].second = Inexact_Vector_3 (x / l, y / l, z / l);
	}

	x_min = y_min = z_min = FLT_MAX;
	x_max = y_max = z_max = -FLT_MAX;
	for (size_t i = 0 ; i < points.size() ; ++i) {
		const Inexact_Point_3 & pt = points[i].first;
		const double &x = pt.x(), &y = pt.y(), &z = pt.z();
		if (x < x_min) x_min = x;
		if (x > x_max) x_max = x;
		if (y < y_min) y_min = y;
		if (y > y_max) y_max = y;
		if (z < z_min) z_min = z;
		if (z > z_max) z_max = z;
	}
	
	std::cout << "X : [" << x_min << ", " << x_max << "]" << std::endl;
	std::cout << "Y : [" << y_min << ", " << y_max << "]" << std::endl;
	std::cout << "Z : [" << z_min << ", " << z_max << "]" << std::endl;

	double dx = x_max - x_min, dy = y_max - y_min, dz = z_max - z_min;
	bbox_diagonal = sqrt(dx * dx + dy * dy + dz * dz);

	std::cout << "Bounding box diagonal : " << bbox_diagonal << std::endl;
}



void Shape_Detector::set_detection_parameters(int _rg_min_points,  int _knn, double _rg_normal_threshold)
{
	if (knn != _knn) {
		knn = _knn;
		should_compute_knn = true;
	}

	min_points = _rg_min_points;
	
	normal_threshold = _rg_normal_threshold;
	
	


}

void Shape_Detector::set_epsilon(double _rg_epsilon) {

	epsilon = _rg_epsilon;
	
}




void Shape_Detector::detect_shapes()
{
	
	detect_planes();

	get_coverage_and_mean_error();
	ori_all_error = all_error;
	ori_coverage = coverage;
	ori_mean_error = mean_error;
	ori_inliers_number = number_of_assigned_points;
	ori_primitives_number = primitives_number;
	ori_mean_normal_diviation = mean_normal_diviation;

	should_compute_neighborhood = true;
	get_last_information();



	
}

void Shape_Detector::save_convex_hull() {
	std::vector<std::vector<Inexact_Point_3> > as_points;
	std::vector<CGAL::Color> as_colors;
	for (size_t p = 0; p < convex_hulls_pts.size(); ++p) {
		if (convex_hulls_pts[p].size() < 3) continue;

		const CGAL::Color & C = convex_hulls_colors[p];


		const std::vector<Inexact_Point_3> & H = convex_hulls_pts[p];
		as_points.push_back(H);
		as_colors.push_back(C);
		
	}

	std::string filename = "./" + path_point_cloud_basename + "_convex_hull.ply";
	std::ofstream stream(filename, std::ios::out);
	if (!stream.is_open()) {
		throw std::ios_base::failure("Error : cannot write into an output file");
	}



	stream.precision(12);

	size_t nb_pts = 0;
	for (size_t i = 0; i < as_points.size(); i++) nb_pts += as_points[i].size();

	stream << "ply" << std::endl;
	stream << "format ascii 1.0" << std::endl;
	stream << "element vertex " << nb_pts << std::endl;
	stream << "property float x" << std::endl;
	stream << "property float y" << std::endl;
	stream << "property float z" << std::endl;
	stream << "property uchar red" << std::endl;
	stream << "property uchar green" << std::endl;
	stream << "property uchar blue" << std::endl;
	stream << "element face " << as_points.size() << std::endl;
	stream << "property list uchar int vertex_index" << std::endl;
	stream << "end_header" << std::endl;

	// Vertices

	for (size_t i = 0; i < as_points.size(); ++i) {
		for (size_t j = 0; j < as_points[i].size(); ++j) {
			stream << as_points[i][j].x() << " " << as_points[i][j].y() << " " << as_points[i][j].z()
				<< " " << (int)as_colors[i].r() << " " << (int)as_colors[i].g() << " " << (int)as_colors[i].b() << std::endl;
		}
	}

	// Facets

	size_t cont = 0;
	for (size_t i = 0; i < as_points.size(); i++) {
		stream << as_points[i].size() << " ";
		for (size_t k = cont; k < cont + as_points[i].size(); k++) {
			stream << k << " ";
		}
		stream << std::endl;
		cont = cont + as_points[i].size();
	}

	stream << std::endl;
	stream.close();
}


void Shape_Detector::save_alpha_shapes() {
	std::vector<std::vector<Inexact_Point_3> > as_points;
	std::vector<CGAL::Color> as_colors;
	for (size_t p = 0; p < convex_hulls_pts.size(); ++p) {
		if (convex_hulls_pts[p].size()<3) continue;

		const CGAL::Color & C = convex_hulls_colors[p];

		
			// Alpha-shapes
			std::pair<int, int> tr_range = convex_hulls_seeds[p];
			for (int tr_id = tr_range.first; tr_id <= tr_range.second; ++tr_id) {
				const std::vector<Inexact_Point_3> & T = alpha_shapes_pts[tr_id];
				if (T.size() == 3) {
					as_points.push_back(T);
					as_colors.push_back(C);
				}
			}
		}

	std::string filename = "./" + path_point_cloud_basename + "_alpha_shape.ply";
	std::ofstream stream(filename, std::ios::out);
	if (!stream.is_open()) {
		throw std::ios_base::failure("Error : cannot write into an output file");
	}

	

	stream.precision(12);

	size_t nb_pts = 0;
	for (size_t i = 0; i < as_points.size(); i++) nb_pts += as_points[i].size();

	stream << "ply" << std::endl;
	stream << "format ascii 1.0" << std::endl;
	stream << "element vertex " << nb_pts << std::endl;
	stream << "property float x" << std::endl;
	stream << "property float y" << std::endl;
	stream << "property float z" << std::endl;
	stream << "property uchar red" << std::endl;
	stream << "property uchar green" << std::endl;
	stream << "property uchar blue" << std::endl;
	stream << "element face " << as_points.size() << std::endl;
	stream << "property list uchar int vertex_index" << std::endl;
	stream << "end_header" << std::endl;

	// Vertices

	for (size_t i = 0; i < as_points.size(); ++i) {
		for (size_t j = 0; j < as_points[i].size(); ++j) {
			stream << as_points[i][j].x() << " " << as_points[i][j].y() << " " << as_points[i][j].z()
				<< " " << (int)as_colors[i].r() << " " << (int)as_colors[i].g() << " " << (int)as_colors[i].b() << std::endl;
		}
	}

	// Facets

	size_t cont = 0;
	for (size_t i = 0; i < as_points.size(); i++) {
		stream << as_points[i].size() << " ";
		for (size_t k = cont; k < cont + as_points[i].size(); k++) {
			stream << k << " ";
		}
		stream << std::endl;
		cont = cont + as_points[i].size();
	}

	stream << std::endl;
	stream.close();
}


void Shape_Detector::save_infor_to_txt_right_folder(std::string name) {
	/*QDir subdirectory(QDir::currentPath() + QString("\\")
		+ QString::fromStdString(path_point_cloud_basename));
	subdirectory.mkpath(".");*/
	//std::string fileName = subdirectory.path().toStdString() + "\\" + name + ".txt";
	std::string fileName = name;
	std::ofstream stream(fileName.c_str(), std::ios::out);

	if (!stream.is_open()) {
		throw std::ios_base::failure("Error : cannot write into an output file");
	}
	stream << "***All operators: " << interval_all << " s." << std::endl;
	stream << "**Transfer operator: " << all_t_transfer << " times" << std::endl;
	stream << "**Merge operator: " << all_t_merge << " times" << std::endl;

	stream << "**Split operator: " << all_t_split << " times" << std::endl;

	stream << "**Insert operator: " << all_t_insert << " times" << std::endl;
	stream << "**Exclude operator: " << all_t_exlude << " times" << std::endl;
	stream << "**number of iterations: " << number_iterations << std::endl;
	stream << "Primitives : " << (primitives_number) << std::endl;
	stream << "Coverage  : " << (coverage) << std::endl;
	stream << "Mean error : " << (mean_error) << std::endl;
	stream << "Mean normal deviation : " << (mean_normal_diviation) << std::endl;
	stream << std::endl;
	if (refinement_type != "unknow") {
		stream << "***Energy" << std::endl
			<< "Initialized energy: " << get_original_energy(refinement_type, weight_mode) << std::endl
			<< "Final energy: " << get_final_energy(refinement_type, weight_mode) << std::endl;
	}


	stream.close();
}





void Shape_Detector::planar_shape_detection_L1() {

	int timett = 0;
	bool propagation2;

	get_last_information();

	region_type = std::vector<int>(planes_to_inliers.size(), 0);
	region_type_m = std::vector<int>(planes_to_inliers.size(), 0);
	points_if_added = std::vector<int>(points.size(), -1);

	old_size_current_primitives = last_primitives_number;
	old_mean_normal_diaviation = last_normal_deviation;
	old_mean_distance_diaviation = last_mean_error;
	old_coverage = last_coverage;
	all_t_exlude = 0;
	all_t_insert = 0;
	all_t_transfer = 0;
	all_t_merge = 0;
	all_t_split = 0;
	clock_t t_start_all = clock();
	do {
		timett++;

		number_of_insert_exclude *= ((timett) / 7 + 1);
		propagation2 = true;



		test_connected_primitives();


		clock_t t_start1 = clock();
		get_distance_diviation();
		number_inlier_before_opers = number_of_assigned_points;
		good_points_shape.clear();
		bad_points_shape.clear();
		if (lambda_c > 0) {
			get_good_points_normal();
		}
		//the first cycle, we do not exclude the point, since the initialized configuration is bad, there are plenty of points would be excluded befor transfer.
		//And the excluded points will be inserted later. 
		if (timett > 1) {
			get_bad_points_normal();
		}



		local_operators_normal();

		clock_t t_end1 = clock();
		double local_time = double(t_end1 - t_start1) / CLOCKS_PER_SEC;

		get_distance_diviation_show_merge_info(local_time);



		clock_t t_start = clock();


		//get_distance_diviation();
		mean_normal_diviation_befor_transfer = mean_normal_current;
		transfer_operator_normal();



		clock_t t_end = clock();
		double trsfer_time = double(t_end - t_start) / CLOCKS_PER_SEC;

		all_t_transfer++;
		get_distance_diviation_show_normal_info(trsfer_time);
		number_iterations++;
	
	

		if ((t_m == 0 && timett > 1) || timett > 800 || (timett > 500 && t_split == t_merge)|| timett>= stop_iterations) propagation2 = false;
	} while (propagation2);


	clock_t t_end_all = clock();
	interval_all = double(t_end_all - t_start_all) / CLOCKS_PER_SEC;

	get_distance_diviation();


	get_coverage_and_mean_error_pure();
	planes_1 = planes_2;
	planes_0 = planes_2;


}

void Shape_Detector::planar_shape_detection_hybrid() {

	//test_curvature();
	int timett = 0;
	bool propagation2;
	//SD->set_primitives(true);
	get_last_information();


	region_type = std::vector<int>(planes_to_inliers.size(), 0);
	region_type_m = std::vector<int>(planes_to_inliers.size(), 0);
	points_if_added = std::vector<int>(points.size(), -1);


	old_size_current_primitives = last_primitives_number;
	old_mean_normal_diaviation = last_normal_deviation;
	old_mean_distance_diaviation = last_mean_error;
	old_coverage = last_coverage;
	all_t_exlude = 0;
	all_t_insert = 0;
	all_t_transfer = 0;
	all_t_merge = 0;
	all_t_split = 0;
	clock_t t_start_all = clock();

	do {
		timett++;
		propagation2 = true;

		number_of_insert_exclude *= ((timett) / 7 + 1);

		test_connected_primitives();



		clock_t t_start1 = clock();
		good_points_shape.clear();
		bad_points_shape.clear();
		get_distance_diviation();
		number_inlier_before_opers = number_of_assigned_points;

		if (lambda_c > 0) {
			get_good_points();
		}
		if (timett > 1) {
			get_bad_points();

		}


		local_operators();


		clock_t t_end1 = clock();
		double local_time = double(t_end1 - t_start1) / CLOCKS_PER_SEC;

		get_distance_diviation_show_merge_info(local_time);
		clock_t t_start = clock();

		mean_distance_diaviation_before_transfer = mean_distance_diaviation;

		transfer_operator_normal_for_hybrid();



		clock_t t_end = clock();
		double trsfer_time = double(t_end - t_start) / CLOCKS_PER_SEC;

		all_t_transfer++;
		get_distance_diviation_show_normal_info(trsfer_time);
		number_iterations++;
	



		if ((t_m == 0 && timett > 1) || timett > 800 || (timett > 500 && t_split == t_merge)||timett >= stop_iterations) propagation2 = false;
	} while (propagation2);


	get_distance_diviation();
	clock_t t_end_all = clock();
	interval_all = double(t_end_all - t_start_all) / CLOCKS_PER_SEC;

	get_coverage_and_mean_error_pure();
	planes_1 = planes_2;
	planes_0 = planes_2;

}

void Shape_Detector::planar_shape_detection_l2() {


	int timett = 0;
	bool propagation2;

	get_last_information();


	region_type = std::vector<int>(planes_to_inliers.size(), 0);
	region_type_m = std::vector<int>(planes_to_inliers.size(), 0);
	points_if_added = std::vector<int>(points.size(), -1);

	old_size_current_primitives = last_primitives_number;

	old_mean_distance_diaviation = last_mean_error;
	old_coverage = last_coverage;
	all_t_exlude = 0;
	all_t_insert = 0;
	all_t_transfer = 0;
	all_t_merge = 0;
	all_t_split = 0;
	clock_t t_start_all = clock();

	do {
		timett++;
		propagation2 = true;
		number_of_insert_exclude *= ((timett) / 7 + 1);

		test_connected_primitives();
		good_points_shape.clear();
		bad_points_shape.clear();

		clock_t t_start1 = clock();
		get_distance_diviation();
		number_inlier_before_opers = number_of_assigned_points;
		if (lambda_c > 0) {
			get_good_points();
		}
		if (timett > 1) {
			get_bad_points();
		}

		local_operators();



		clock_t t_end1 = clock();
		double local_time = double(t_end1 - t_start1) / CLOCKS_PER_SEC;


		get_distance_diviation_show_merge_info(local_time);

		clock_t t_start = clock();
		mean_distance_diaviation_before_transfer = mean_distance_diaviation;


		transfer_operator_l2();

		clock_t t_end = clock();
		double trsfer_time = double(t_end - t_start) / CLOCKS_PER_SEC;

		all_t_transfer++;
		get_distance_diviation_show_normal_info(trsfer_time);

		number_iterations++;

		if ((t_m == 0 && timett > 1) || timett > 800 || (timett > 500 && t_split == t_merge) || timett >= stop_iterations) propagation2 = false;
	} while (propagation2);








	get_distance_diviation();
	clock_t t_end_all = clock();
	interval_all = double(t_end_all - t_start_all) / CLOCKS_PER_SEC;

	get_coverage_and_mean_error_pure();
	planes_1 = planes_2;
	planes_0 = planes_2;






}

bool Shape_Detector::separate_two_out(int id, std::vector<int> & max_list, std::vector<int> & min_list, double & dif) {
	std::vector<int> this_region = planes_to_inliers[id];
	Inexact_Plane this_plane = planes_2[id];

	double max_value = FLT_MIN;
	int max_index = -1;
	double min_value = FLT_MIN;
	int min_index = -1;
	double dis_bef = 0;
	//divide into two regions according to the two farthest points.
	double this_mean_normal_diviation = 0;
	Inexact_Vector_3 plane_direction = this_plane.orthogonal_vector();

	for (int ii = 0; ii < this_region.size(); ++ii) {
		double this_dis = sqrt(CGAL::squared_distance(this_plane, points[this_region[ii]].first));
		double this_normal_d = std::abs(points[this_region[ii]].second*plane_direction);
		this_mean_normal_diviation += this_normal_d;
		dis_bef += this_dis;
		if (this_plane.has_on_positive_side(points[this_region[ii]].first)) {
			if (this_dis > max_value) {
				max_value = this_dis;
				max_index = this_region[ii];

			}
		}
		else {
			if (this_dis > min_value) {
				min_value = this_dis;
				min_index = this_region[ii];
			}

		}
	}
	this_mean_normal_diviation /= double(this_region.size());
	//we do not split the primitive that has smaller normal diviation than ori mean.



	if (this_mean_normal_diviation > ori_mean_normal_diviation) {
		min_list.clear();
		max_list.clear();
		dif = 0;
		return false;
	}


	std::map<int, int> label_points;
	Inexact_Point_3 max_point = points[max_index].first;
	Inexact_Point_3 min_point = points[min_index].first;
	std::vector<Inexact_Point_3> max_point_list;
	std::vector<Inexact_Point_3> min_point_list;
	Inexact_Point_2 max_point_2d = this_plane.to_2d(max_point);
	Inexact_Point_2 min_point_2d = this_plane.to_2d(min_point);


	for (int j = 0; j < this_region.size(); ++j) {
		if ((max_point_2d - this_plane.to_2d(points[this_region[j]].first)).squared_length() < (min_point_2d - this_plane.to_2d(points[this_region[j]].first)).squared_length()) {
			max_list.push_back(this_region[j]);
			label_points[this_region[j]] = 1;
			max_point_list.push_back(points[this_region[j]].first);
		}
		else {
			label_points[this_region[j]] = -1;
			min_point_list.push_back(points[this_region[j]].first);
			min_list.push_back(this_region[j]);
		}
	}

	//transfer between max and min lists.
	Inexact_Plane plane_max;
	Inexact_Plane plane_min;

	if (max_point_list.size() < 3 || min_point_list.size() < 3) {
		return false;
	}


	linear_least_squares_fitting_3(max_point_list.begin(), max_point_list.end(), plane_max, CGAL::Dimension_tag<0>());
	linear_least_squares_fitting_3(min_point_list.begin(), min_point_list.end(), plane_min, CGAL::Dimension_tag<0>());
	bool propo = false;
	int Nb_neigh = 10;
	int refine_time = 0;

	do {
		refine_time++;
		int moving_n = 0;
		propo = true;
		std::vector<int> if_moved = std::vector<int>(this_region.size(), -1);
		for (int f = 0; f < this_region.size(); ++f) {

			if (if_moved[f] > 2) continue;
			int this_label = label_points[this_region[f]];
			if ((int)spherical_neighborhood[this_region[f]].size() < Nb_neigh) {
				Nb_neigh = (int)spherical_neighborhood[this_region[f]].size();
			}
			bool bb = false;
			for (int it = 0; it < Nb_neigh; it++) {

				int neighbor_index = spherical_neighborhood[this_region[f]][it];
				if (inliers_to_planes[neighbor_index] != id) continue; //the neighbor point should be in the splitted region
				if (label_points[neighbor_index] != this_label) {
					bb = true;
					break;
				}


			}
			if (bb == false) continue;
			Inexact_Point_3 this_point = points[this_region[f]].first;



			if (this_label == -1) {

				if (abs(points[this_region[f]].second * plane_max.orthogonal_vector()) > normal_threshold) {
					if (CGAL::squared_distance(plane_max, this_point) < CGAL::squared_distance(plane_min, this_point)) {																																														 //if ((1 - that_cos)*sqrt((this_point - planes_2[neight_id].projection(this_point)).squared_length()) < min_sin) {
						label_points[this_region[f]] = -this_label;
						moving_n++;
						if_moved[f] ++;
					}
				}
			}
			else {

				if (abs(points[this_region[f]].second * plane_min.orthogonal_vector()) > normal_threshold) {
					if (CGAL::squared_distance(plane_min, this_point) < CGAL::squared_distance(plane_max, this_point)) {																																													 //if ((1 - that_cos)*sqrt((this_point - planes_2[neight_id].projection(this_point)).squared_length()) < min_sin) {
						label_points[this_region[f]] = -this_label;
						moving_n++;
						if_moved[f] ++;
					}
				}
			}

		}

		//update the two regions
		std::map<int, int>::iterator iter;
		iter = label_points.begin();
		max_list.clear();
		min_list.clear();
		max_point_list.clear();
		min_point_list.clear();
		while (iter != label_points.end()) {
			if (iter->second == 1) {
				max_list.push_back(iter->first);
				max_point_list.push_back(points[iter->first].first);

			}
			else {
				min_list.push_back(iter->first);
				min_point_list.push_back(points[iter->first].first);
			}
			iter++;
		}
		if (max_point_list.size() < min_points || min_point_list.size() < min_points) {

			return false;
		}
		linear_least_squares_fitting_3(max_point_list.begin(), max_point_list.end(), plane_max, CGAL::Dimension_tag<0>());
		linear_least_squares_fitting_3(min_point_list.begin(), min_point_list.end(), plane_min, CGAL::Dimension_tag<0>());
		if (moving_n < 5 || refine_time>25) {
			propo = false;

		}

	} while (propo);

	double dis_after = 0;

	for (int i = 0; i < max_list.size(); ++i) {

		dis_after += sqrt(CGAL::squared_distance(plane_max, points[max_list[i]].first));

	}

	for (int i = 0; i < min_list.size(); ++i) {

		dis_after += sqrt(CGAL::squared_distance(plane_min, points[min_list[i]].first));

	}

	dif = dis_after - dis_bef;

	return true;


}
bool Shape_Detector::separate_two_out_normal(int id, std::vector<int> & max_list, std::vector<int> & min_list, double & dif) {
	std::vector<int> this_region = planes_to_inliers[id];
	Inexact_Plane this_plane = planes_2[id];

	double max_value = FLT_MIN;
	int max_index = -1;
	double min_value = FLT_MIN;
	int min_index = -1;
	double dis_bef = 0;
	double mean_distance = 0;
	//divide into two regions

	Inexact_Vector_3 plane_direction = this_plane.orthogonal_vector();

	for (int ii = 0; ii < this_region.size(); ++ii) {

		double this_dis = sqrt(CGAL::squared_distance(this_plane, points[this_region[ii]].first));
		double this_normal_d = std::abs(points[this_region[ii]].second*plane_direction);
		mean_distance += this_dis;
		dis_bef += this_normal_d;
		if (this_plane.has_on_positive_side(points[this_region[ii]].first)) {
			if (this_dis > max_value) {
				max_value = this_dis;
				max_index = this_region[ii];

			}
		}
		else {
			if (this_dis > min_value) {
				min_value = this_dis;
				min_index = this_region[ii];
			}

		}
	}
	mean_distance /= double(this_region.size());
	//we do not split the primitive that has smaller distance error than the ori mean.
	if (mean_distance < ori_mean_error) {
		min_list.clear();
		max_list.clear();
		dif = 0;
		return false;
	}
	std::map<int, int> label_points;
	Inexact_Point_3 max_point = points[max_index].first;
	Inexact_Point_3 min_point = points[min_index].first;
	std::vector<Inexact_Point_3> max_point_list;
	std::vector<Inexact_Point_3> min_point_list;
	Inexact_Point_2 max_point_2d = this_plane.to_2d(max_point);
	Inexact_Point_2 min_point_2d = this_plane.to_2d(min_point);


	for (int j = 0; j < this_region.size(); ++j) {
		if ((max_point_2d - this_plane.to_2d(points[this_region[j]].first)).squared_length() < (min_point_2d - this_plane.to_2d(points[this_region[j]].first)).squared_length()) {
			max_list.push_back(this_region[j]);
			label_points[this_region[j]] = 1;
			max_point_list.push_back(points[this_region[j]].first);
		}
		else {
			label_points[this_region[j]] = -1;
			min_point_list.push_back(points[this_region[j]].first);
			min_list.push_back(this_region[j]);
		}
	}

	//transformation
	Inexact_Plane plane_max;
	Inexact_Plane plane_min;
	if (max_point_list.size() < 3 || min_point_list.size() < 3) {

		return false;
	}
	linear_least_squares_fitting_3(max_point_list.begin(), max_point_list.end(), plane_max, CGAL::Dimension_tag<0>());
	linear_least_squares_fitting_3(min_point_list.begin(), min_point_list.end(), plane_min, CGAL::Dimension_tag<0>());
	bool propo = false;
	int Nb_neigh = 10;
	int refine_time = 0;
	do {
		refine_time++;
		int moving_n = 0;
		propo = true;
		std::vector<int> if_moved = std::vector<int>(this_region.size(), -1);
		for (int f = 0; f < this_region.size(); ++f) {
			if (if_moved[f] > 2) continue;
			int this_label = label_points[this_region[f]];
			if ((int)spherical_neighborhood[this_region[f]].size() < Nb_neigh) {
				Nb_neigh = (int)spherical_neighborhood[this_region[f]].size();
			}
			bool bb = false;
			for (int it = 0; it < Nb_neigh; it++) {

				int neighbor_index = spherical_neighborhood[this_region[f]][it];
				if (inliers_to_planes[neighbor_index] != id) continue; //the neighbor point should be in the splitted region
				if (label_points[neighbor_index] != this_label) {
					bb = true;
					break;
				}


			}
			if (bb == false) continue;
			Inexact_Point_3 this_point = points[this_region[f]].first;



			if (this_label == -1) {
				if (CGAL::squared_distance(plane_max, this_point) < epsilon) {
					if (abs(points[this_region[f]].second * plane_max.orthogonal_vector()) > abs(points[this_region[f]].second * plane_min.orthogonal_vector())) {
						//if ((1 - that_cos)*sqrt((this_point - planes_2[neight_id].projection(this_point)).squared_length()) < min_sin) {
						label_points[this_region[f]] = -this_label;
						moving_n++;
						if_moved[f]++;
					}
				}
			}
			else {
				if (CGAL::squared_distance(plane_min, this_point) < epsilon) {


					if (abs(points[this_region[f]].second * plane_max.orthogonal_vector()) < abs(points[this_region[f]].second * plane_min.orthogonal_vector())) {//need to changed for smooth																																														 //if ((1 - that_cos)*sqrt((this_point - planes_2[neight_id].projection(this_point)).squared_length()) < min_sin) {
						label_points[this_region[f]] = -this_label;
						moving_n++;
						if_moved[f]++;
					}
				}
			}







		}
		std::map<int, int>::iterator iter;
		iter = label_points.begin();
		max_list.clear();
		min_list.clear();
		max_point_list.clear();
		min_point_list.clear();
		while (iter != label_points.end()) {
			if (iter->second == 1) {
				max_list.push_back(iter->first);
				max_point_list.push_back(points[iter->first].first);

			}
			else {
				min_list.push_back(iter->first);
				min_point_list.push_back(points[iter->first].first);
			}
			iter++;
		}
		if (max_point_list.size() < min_points || min_point_list.size() < min_points) {

			return false;
		}
		linear_least_squares_fitting_3(max_point_list.begin(), max_point_list.end(), plane_max, CGAL::Dimension_tag<0>());
		linear_least_squares_fitting_3(min_point_list.begin(), min_point_list.end(), plane_min, CGAL::Dimension_tag<0>());
		if (moving_n < 5 || refine_time>25) {
			propo = false;
		}

	} while (propo);

	double dis_after = 0;

	for (int i = 0; i < max_list.size(); ++i) {
		dis_after += abs(points[max_list[i]].second * plane_max.orthogonal_vector());

	}

	for (int i = 0; i < min_list.size(); ++i) {
		dis_after += abs(points[min_list[i]].second * plane_min.orthogonal_vector());
	}

	dif = dis_after - dis_bef;

	return true;


}





double Shape_Detector::energy_changed(double dis, double numb) {
	


	double dem_number_shape = double(points.size()) / double(min_points);


	if (weight_mode == 2) {
		double term1 = double(3 - lambda_c - lambda_r)*(dis / double(number_of_assigned_points)) / (mean_distance_current);
		double term2 = lambda_r * numb / (double(size_current_primitives));

		return(term1 + term2);
	}
	else if (weight_mode == 1) {
		double term1 = double(3 - lambda_c - lambda_r)*(dis / double(number_of_assigned_points)) / (ori_mean_error);
		double term2 = lambda_r * numb / (double(ori_primitives_number));

		return(term1 + term2);

	}
	else {
		double term1 = double(3 - lambda_c - lambda_r)*(dis / double(number_of_assigned_points)) / (epsilon);
		double term2 = lambda_r * numb / (double(dem_number_shape));

		return(term1 + term2);


	}


}

double Shape_Detector::get_final_energy(std::string meseurment, int mm) {
	

	double dem_number_shape = double(points.size()) / double(min_points);

	if (meseurment == "L2" || meseurment == "L_hybrid") {
		if (mm == 2) {
			double weight1 = (1.0 / double(old_mean_distance_diaviation)) / ((1.0 / double(old_mean_distance_diaviation)) + (1.0 / double(old_coverage)) + (1.0 / double(old_size_current_primitives)));
			double weight2 = (1.0 / double(old_size_current_primitives)) / ((1.0 / double(old_mean_distance_diaviation)) + (1.0 / double(old_coverage)) + (1.0 / double(old_size_current_primitives)));
			double weight3 = (1.0 / double(old_coverage)) / ((1.0 / double(old_mean_distance_diaviation)) + (1.0 / double(old_coverage)) + (1.0 / double(old_size_current_primitives)));

			double term1 = weight1 * double(3 - lambda_c - lambda_r)*(mean_distance_current) / epsilon;
			double term2 = weight2 * lambda_r * size_current_primitives / (double(dem_number_shape));
			double term3 = weight3 * lambda_c * (1 - number_of_assigned_points / double(points.size()));

			return(term1 / 3.0 + term2 / 3.0 + term3 / 3.0);
		}
		else if (mm == 1) {
			double weight1 = (1.0 / double(ori_mean_error)) / ((1.0 / double(ori_mean_error)) + (1.0 / double(ori_coverage)) + (1.0 / double(ori_primitives_number)));
			double weight2 = (1.0 / double(ori_primitives_number)) / ((1.0 / double(ori_mean_error)) + (1.0 / double(ori_coverage)) + (1.0 / double(ori_primitives_number)));
			double weight3 = (1.0 / double(ori_coverage)) / ((1.0 / double(ori_mean_error)) + (1.0 / double(ori_coverage)) + (1.0 / double(ori_primitives_number)));

			double term1 = weight1 * double(3 - lambda_c - lambda_r)*(mean_distance_current) / epsilon;
			double term2 = weight2 * lambda_r * size_current_primitives / (double(dem_number_shape));
			double term3 = weight3 * lambda_c * (1 - number_of_assigned_points / double(points.size()));

			return(term1 / 3.0 + term2 / 3.0 + term3 / 3.0);

		}
		else {

			double term1 = double(3 - lambda_c - lambda_r)*(mean_distance_current) / epsilon;
			double term2 = lambda_r * size_current_primitives / (double(dem_number_shape));
			double term3 = lambda_c * (1 - number_of_assigned_points / double(points.size()));

			return(term1 / 3.0 + term2 / 3.0 + term3 / 3.0);


		}
	}
	else if (meseurment == "L1.1") {
		if (mm == 2) {
			double weight1 = (1.0 / double(old_mean_normal_diaviation)) / ((1.0 / double(old_mean_normal_diaviation)) + (1.0 / double(old_coverage)) + (1.0 / double(old_size_current_primitives)));
			double weight2 = (1.0 / double(old_size_current_primitives)) / ((1.0 / double(old_mean_normal_diaviation)) + (1.0 / double(old_coverage)) + (1.0 / double(old_size_current_primitives)));
			double weight3 = (1.0 / double(old_coverage)) / ((1.0 / double(old_mean_normal_diaviation)) + (1.0 / double(old_coverage)) + (1.0 / double(old_size_current_primitives)));

			double term1 = weight1 * double(3 - lambda_c - lambda_r)*(1 - mean_normal_current / normal_threshold);
			double term2 = weight2 * lambda_r * size_current_primitives / (double(dem_number_shape));
			double term3 = weight3 * lambda_c * (1 - number_of_assigned_points / double(points.size()));

			return(term1 / 3.0 + term2 / 3.0 + term3 / 3.0);
		}
		else if (mm == 1) {
			double weight1 = (1.0 / double(ori_mean_normal_diviation)) / ((1.0 / double(ori_mean_normal_diviation)) + (1.0 / double(ori_coverage)) + (1.0 / double(ori_primitives_number)));
			double weight2 = (1.0 / double(ori_primitives_number)) / ((1.0 / double(ori_mean_normal_diviation)) + (1.0 / double(ori_coverage)) + (1.0 / double(ori_primitives_number)));
			double weight3 = (1.0 / double(ori_coverage)) / ((1.0 / double(ori_mean_normal_diviation)) + (1.0 / double(ori_coverage)) + (1.0 / double(ori_primitives_number)));

			double term1 = weight1 * double(3 - lambda_c - lambda_r)*(1 - mean_normal_current / normal_threshold);

			double term2 = weight2 * lambda_r * size_current_primitives / (double(dem_number_shape));
			double term3 = weight3 * lambda_c * (1 - number_of_assigned_points / double(points.size()));

			return(term1 / 3.0 + term2 / 3.0 + term3 / 3.0);

		}
		else {

			double term1 = double(3 - lambda_c - lambda_r)*((1 - mean_normal_current) / (1 - normal_threshold));

			double term2 = lambda_r * size_current_primitives / (double(dem_number_shape));
			double term3 = lambda_c * (1 - number_of_assigned_points / double(points.size()));

			return(term1 / 3.0 + term2 / 3.0 + term3 / 3.0);


		}


	}
	else {
		throw "Using indefinded metric";
	}
}

double Shape_Detector::get_original_energy(std::string meseurment, int mm) {
	
	double dem_number_shape = double(points.size()) / double(min_points);

	if (meseurment == "L2" || meseurment == "L_hybrid") {
		if (mm == 1 || mm == 2) {
			double weight1 = (1.0 / double(ori_mean_error)) / ((1.0 / double(ori_mean_error)) + (1.0 / double(ori_coverage)) + (1.0 / double(ori_primitives_number)));
			double weight2 = (1.0 / double(ori_primitives_number)) / ((1.0 / double(ori_mean_error)) + (1.0 / double(ori_coverage)) + (1.0 / double(ori_primitives_number)));
			double weight3 = (1.0 / double(ori_coverage)) / ((1.0 / double(ori_mean_error)) + (1.0 / double(ori_coverage)) + (1.0 / double(ori_primitives_number)));

			double term1 = weight1 * double(3 - lambda_c - lambda_r)*(ori_mean_error) / epsilon;
			double term2 = weight2 * lambda_r * ori_primitives_number / (double(dem_number_shape));
			double term3 = weight3 * lambda_c * (1 - ori_inliers_number / double(points.size()));

			return(term1 / 3.0 + term2 / 3.0 + term3 / 3.0);

		}
		else {

			double term1 = double(3 - lambda_c - lambda_r)*(ori_mean_error) / epsilon;
			double term2 = lambda_r * ori_primitives_number / (double(dem_number_shape));
			double term3 = lambda_c * (1 - ori_inliers_number / double(points.size()));

			return(term1 / 3.0 + term2 / 3.0 + term3 / 3.0);


		}
	}
	else if (meseurment == "L1.1") {
		if (mm == 1 || mm == 2) {
			double weight1 = (1.0 / double(ori_mean_normal_diviation)) / ((1.0 / double(ori_mean_normal_diviation)) + (1.0 / double(ori_coverage)) + (1.0 / double(ori_primitives_number)));
			double weight2 = (1.0 / double(ori_primitives_number)) / ((1.0 / double(ori_mean_normal_diviation)) + (1.0 / double(ori_coverage)) + (1.0 / double(ori_primitives_number)));
			double weight3 = (1.0 / double(ori_coverage)) / ((1.0 / double(ori_mean_normal_diviation)) + (1.0 / double(ori_coverage)) + (1.0 / double(ori_primitives_number)));

			double term1 = (weight1) * double(3 - lambda_c - lambda_r)*(1 - ori_mean_normal_diviation / normal_threshold);
			double term2 = (weight2)* lambda_r * ori_primitives_number / (double(dem_number_shape));
			double term3 = (weight3)* lambda_c * (1 - ori_inliers_number / double(points.size()));

			return(term1 / 3.0 + term2 / 3.0 + term3 / 3.0);

		}
		else {

			double term1 = double(3 - lambda_c - lambda_r)*((1 - ori_mean_normal_diviation) / (1 - normal_threshold));
			double term2 = lambda_r * ori_primitives_number / (double(dem_number_shape));
			double term3 = lambda_c * (1 - ori_inliers_number / double(points.size()));

			return(term1 / 3.0 + term2 / 3.0 + term3 / 3.0);


		}


	}
	else {
		throw "Using indefinded metric";
	}
}

double Shape_Detector::energy_changed_merge(double dis, double numb, int nmoves) {
	
	double dem_number_shape = double(points.size()) / double(min_points);
	
	if (weight_mode == 2) {
		double term1 = double(3 - lambda_c - lambda_r)*((all_distance_diaviation + dis) / (double(number_of_assigned_points) - nmoves) - mean_distance_current) / (mean_distance_current);

		
		double term2 = lambda_r * numb / (double(size_current_primitives));
		double term3 = double(lambda_c)*nmoves / (double(number_inlier_before_opers));

		return(term1 + term2 + term3);
	}
	else if (weight_mode == 1) {
		double term1 = double(3 - lambda_c - lambda_r)*((all_distance_diaviation + dis) / (double(number_of_assigned_points) - nmoves) - mean_distance_current) / (ori_mean_error);

	
		double term2 = lambda_r * numb / (double(ori_primitives_number));
		double term3 = double(lambda_c)*nmoves / (double(ori_inliers_number));

		return(term1 + term2 + term3);

	}
	else {
		double term1 = double(3 - lambda_c - lambda_r)*((all_distance_diaviation + dis) / (double(number_of_assigned_points) - nmoves) - mean_distance_current) / (epsilon);

		
		double term2 = lambda_r * numb / (double(dem_number_shape));
		double term3 = double(lambda_c)*nmoves / (double(points.size()));

		return(term1 + term2 + term3);


	}


}
double Shape_Detector::energy_changed_normal_merge(double dis, double numb, int nmoves) {

	double fix_number_prim = double(points.size()) / double(min_points);

	if (weight_mode == 2) {
		double term1 = double(lambda_c + lambda_r - 3)*((all_normal_diaviation + dis) / (double(number_of_assigned_points) - nmoves) - mean_normal_current) / (mean_normal_current);

		
		double term2 = lambda_r * numb / (double(size_current_primitives));
		double term3 = double(lambda_c)*nmoves / (double(number_inlier_before_opers));
		return(term1 + term2 + term3);
	}
	else if (weight_mode == 1) {
		double term1 = double(lambda_c + lambda_r - 3)*((all_normal_diaviation + dis) / (double(number_of_assigned_points) - nmoves) - mean_normal_current) / (ori_mean_normal_diviation);

		double term2 = lambda_r * numb / (double(ori_primitives_number));
		double term3 = double(lambda_c)*nmoves / (double(ori_inliers_number));

		return(term1 + term2 + term3);
	}
	else {
		double term1 = double(lambda_c + lambda_r - 3)*((all_normal_diaviation + dis) / (double(number_of_assigned_points) - nmoves) - mean_normal_current) / normal_threshold;

		
		double term2 = lambda_r * numb / double(fix_number_prim);
		double term3 = double(lambda_c)*nmoves / (double(points.size()));

		return(term1 + term2 + term3);
	}

}
double Shape_Detector::energy_changed_normal(double dis, double numb) {
	
	double fix_number_prim = double(points.size()) / double(min_points);
	
	if (weight_mode == 2) {
		double term1 = double(lambda_c + lambda_r - 3)*(dis / double(number_of_assigned_points)) / (mean_normal_current);
		double term2 = lambda_r * numb / (double(size_current_primitives));
		return(term1 + term2);
	}
	else if (weight_mode == 1) {
		double term1 = double(lambda_c + lambda_r - 3)*(dis / double(number_of_assigned_points)) / (ori_mean_normal_diviation);
		double term2 = lambda_r * numb / (double(ori_primitives_number));
		return(term1 + term2);
	}
	else {
		double term1 = double(lambda_c + lambda_r - 3)*(dis / double(number_of_assigned_points)) / (normal_threshold);
		double term2 = lambda_r * numb / double(fix_number_prim);
		return(term1 + term2);
	}

}
double Shape_Detector::energy_changed_second_normal(double dis, double numb) {

	if (weight_mode == 2) {
		double change_fedilite = double(lambda_c + lambda_r - 3)*((all_normal_diaviation + dis) / (double(number_of_assigned_points) - numb) - mean_normal_current) / (mean_normal_current);
		double change_completness = double(lambda_c)*numb / (double(number_inlier_before_opers));

		return (change_fedilite + change_completness);
	}
	else if (weight_mode == 1) {
		double change_fedilite = double(lambda_c + lambda_r - 3)*((all_normal_diaviation + dis) / (double(number_of_assigned_points) - numb) - mean_normal_current) / (ori_mean_normal_diviation);
		double change_completness = double(lambda_c)*numb / (double(ori_inliers_number));

		return (change_fedilite + change_completness);
	}
	else {
		double change_fedilite = double(lambda_c + lambda_r - 3)*((all_normal_diaviation + dis) / (double(number_of_assigned_points) - numb) - mean_normal_current) / (normal_threshold);
		double change_completness = double(lambda_c)*numb / (double(points.size()));

		return (change_fedilite + change_completness);

	}


}


double Shape_Detector::energy_changed_second(double dis, double numb) {


	if (weight_mode == 2) {
		double change_fedilite = double(3 - lambda_c - lambda_r)*((all_distance_diaviation + dis) / (double(number_of_assigned_points) - numb) - mean_distance_current) / (mean_distance_current);
		double change_completness = double(lambda_c)*(numb / (double(number_inlier_before_opers)));
		return (change_fedilite + change_completness);
	}
	else if (weight_mode == 1) {
		double change_fedilite = double(3 - lambda_c - lambda_r)*((all_distance_diaviation + dis) / (double(number_of_assigned_points) - numb) - mean_distance_current) / (ori_mean_error);
		double change_completness = double(lambda_c)*(numb / (double(ori_inliers_number)));
		return (change_fedilite + change_completness);
	}
	else {
		double change_fedilite = double(3 - lambda_c - lambda_r)*((all_distance_diaviation + dis) / (double(number_of_assigned_points) - numb) - mean_distance_current) / (epsilon);
		double change_completness = double(lambda_c)*(numb / (double(ori_inliers_number)));
		return (change_fedilite + change_completness);
	}


}


void Shape_Detector::transfer_operator_normal() {

	bool propagation1;
	int times = 0;
	int number_moving;

	t_l = 0;


	std::vector<int> if_moved;
	if_moved = std::vector<int>(inliers_to_planes.size(), -1);


	std::vector<int> update_primitive_size;
	for (int mm = 0; mm < planes_to_inliers.size(); ++mm) {
		update_primitive_size.push_back((int)planes_to_inliers[mm].size());
	}



	do {
		std::vector<std::vector<int> > planes_to_inliers_last = planes_to_inliers;



		std::vector<Inexact_Plane> planes_2_last = planes_2;
		std::vector<int> inliers_to_planes_last = inliers_to_planes;

		int Nb_neigh = knn;




		times++;
		number_moving = 0;
		propagation1 = true;




		for (size_t i = 0; i < points.size(); ++i) {
			if (inliers_to_planes[i] == -1) continue;
			if (if_moved[i] >= 5) continue;
			if (update_primitive_size[inliers_to_planes[i]] <= min_points) continue;
			if ((int)spherical_neighborhood[i].size() < Nb_neigh) {
				Nb_neigh = (int)spherical_neighborhood[i].size();
			}

			std::set<int> one_neight_id;//connected primitives for the point
			//check if this point is on the boundary of primitives.
			for (int it = 0; it < Nb_neigh; it++) {

				int neighbor_index = spherical_neighborhood[i][it];


				if (inliers_to_planes[neighbor_index] != -1 && inliers_to_planes[i] != inliers_to_planes[neighbor_index]) {

					one_neight_id.insert(inliers_to_planes[neighbor_index]);

				}


			}

			if (one_neight_id.empty()) { continue; }

			Inexact_Point_3 this_point = points[i].first;


			//detect the best primitive
			double max_cos = abs(points[i].second * planes_2[inliers_to_planes[i]].orthogonal_vector());

			int changed_plane_id = -1;
			for (int neight_id : one_neight_id) {


				double that_cos = abs(points[i].second * planes_2[neight_id].orthogonal_vector());
				if (that_cos > max_cos) {
					if (sqrt(CGAL::squared_distance(planes_2[neight_id], this_point)) <= epsilon) {

						max_cos = that_cos;
						changed_plane_id = neight_id;
					}
				}

			}

			if (changed_plane_id != -1) {
				update_primitive_size[inliers_to_planes[i]]--;
				inliers_to_planes[i] = changed_plane_id;
				update_primitive_size[changed_plane_id]++;
				if_moved[i]++;
				number_moving++;
			}

		}
		//update the configuration
		int plane_size = planes_to_inliers.size();

		planes_to_inliers.clear();
		planes_to_inliers.resize(plane_size);
		for (int k = 0; k < points.size(); k++) {
			if (inliers_to_planes[k] != -1) {
				planes_to_inliers[inliers_to_planes[k]].push_back(k);
			}

		}

		planes_2.clear();




		for (int ii = 0; ii < planes_to_inliers.size(); ++ii) {


			std::vector<Inexact_Point_3> inliers_i;
			inliers_i.reserve(planes_to_inliers[ii].size());

			for (int jj = 0; jj < planes_to_inliers[ii].size(); ++jj) {
				const Inexact_Point_3 & pt = points[planes_to_inliers[ii][jj]].first;
				inliers_i.push_back(pt);

			}


			if (inliers_i.size() < 3) {

				getchar();
			}
			Inexact_Plane plane;
			linear_least_squares_fitting_3(inliers_i.begin(), inliers_i.end(), plane, CGAL::Dimension_tag<0>());


			planes_2.push_back(plane);

		}

		get_distance_diviation();
		//make sure that the transfer operation reducing the energy.
		if (mean_normal_diviation < mean_normal_diviation_befor_transfer) {
			propagation1 = false;
			planes_2 = planes_2_last;
			planes_to_inliers = planes_to_inliers_last;

			inliers_to_planes = inliers_to_planes_last;


		}

		if (if_constraint) {

			if (mean_normal_diviation < ori_mean_normal_diviation) {
				propagation1 = false;
				planes_2 = planes_2_last;
				planes_to_inliers = planes_to_inliers_last;

				inliers_to_planes = inliers_to_planes_last;


			}
		}


		t_l++;

		if (number_moving == 0 || times > 100) { propagation1 = false; }
	} while (propagation1);







}

void Shape_Detector::transfer_operator_normal_for_hybrid() {

	bool propagation1;
	int times = 0;
	int number_moving;

	t_l = 0;


	std::vector<int> if_moved;
	if_moved = std::vector<int>(inliers_to_planes.size(), -1);


	std::vector<int> update_primitive_size;
	for (int mm = 0; mm < planes_to_inliers.size(); ++mm) {
		update_primitive_size.push_back((int)planes_to_inliers[mm].size());
	}



	do {
		std::vector<std::vector<int> > planes_to_inliers_last = planes_to_inliers;



		std::vector<Inexact_Plane> planes_2_last = planes_2;
		std::vector<int> inliers_to_planes_last = inliers_to_planes;

		int Nb_neigh = knn;




		times++;
		number_moving = 0;
		propagation1 = true;




		for (size_t i = 0; i < points.size(); ++i) {
			if (inliers_to_planes[i] == -1) continue;
			if (if_moved[i] >= 5) continue;
			if (update_primitive_size[inliers_to_planes[i]] <= min_points) continue;
			if ((int)spherical_neighborhood[i].size() < Nb_neigh) {
				Nb_neigh = (int)spherical_neighborhood[i].size();
			}

			std::set<int> one_neight_id;//connected primitives for the point

			for (int it = 0; it < Nb_neigh; it++) {

				int neighbor_index = spherical_neighborhood[i][it];


				if (inliers_to_planes[neighbor_index] != -1 && inliers_to_planes[i] != inliers_to_planes[neighbor_index]) {

					one_neight_id.insert(inliers_to_planes[neighbor_index]);

				}


			}

			if (one_neight_id.empty()) { continue; }

			Inexact_Point_3 this_point = points[i].first;


			//detect the best primitive
			double max_cos = abs(points[i].second * planes_2[inliers_to_planes[i]].orthogonal_vector());

			int changed_plane_id = -1;
			for (int neight_id : one_neight_id) {


				double that_cos = abs(points[i].second * planes_2[neight_id].orthogonal_vector());
				if (that_cos > max_cos) {
					if (sqrt(CGAL::squared_distance(planes_2[neight_id], this_point)) <= epsilon) {

						max_cos = that_cos;
						changed_plane_id = neight_id;
					}
				}

			}

			if (changed_plane_id != -1) {
				update_primitive_size[inliers_to_planes[i]]--;
				inliers_to_planes[i] = changed_plane_id;
				update_primitive_size[changed_plane_id]++;
				if_moved[i]++;
				number_moving++;
			}

		}
		//update the configuration
		int plane_size = planes_to_inliers.size();

		planes_to_inliers.clear();
		planes_to_inliers.resize(plane_size);
		for (int k = 0; k < points.size(); k++) {
			if (inliers_to_planes[k] != -1) {
				planes_to_inliers[inliers_to_planes[k]].push_back(k);
			}

		}

		planes_2.clear();




		for (int ii = 0; ii < planes_to_inliers.size(); ++ii) {


			std::vector<Inexact_Point_3> inliers_i;
			inliers_i.reserve(planes_to_inliers[ii].size());

			for (int jj = 0; jj < planes_to_inliers[ii].size(); ++jj) {
				const Inexact_Point_3 & pt = points[planes_to_inliers[ii][jj]].first;
				inliers_i.push_back(pt);

			}


			if (inliers_i.size() < 3) {

		
				getchar();
			}
			Inexact_Plane plane;
			linear_least_squares_fitting_3(inliers_i.begin(), inliers_i.end(), plane, CGAL::Dimension_tag<0>());


			planes_2.push_back(plane);

		}

		get_distance_diviation();
		//make sure that the transfer operation reducing the energy.

		if (mean_distance_diaviation > mean_distance_diaviation_before_transfer) {
			propagation1 = false;
			planes_2 = planes_2_last;
			planes_to_inliers = planes_to_inliers_last;

			inliers_to_planes = inliers_to_planes_last;


		}

		if (if_constraint) {

			if (mean_distance_diaviation < ori_mean_error) {
				propagation1 = false;
				planes_2 = planes_2_last;
				planes_to_inliers = planes_to_inliers_last;

				inliers_to_planes = inliers_to_planes_last;


			}
		}


		t_l++;

		if (number_moving == 0 || times > 100) { propagation1 = false; }
	} while (propagation1);







}


void Shape_Detector::transfer_operator_l2() {

	bool propagation1;
	int times = 0;
	int number_moving;



	t_l = 0;

	std::vector<int> if_moved;
	if_moved = std::vector<int>(inliers_to_planes.size(), -1);
	std::vector<int> update_primitive_size;
	for (int mm = 0; mm < planes_to_inliers.size(); ++mm) {
		update_primitive_size.push_back(planes_to_inliers[mm].size());
	}
	int Nb_neigh = knn;
	do {
		std::vector<std::vector<int> > planes_to_inliers_last = planes_to_inliers;
		std::vector<int> region_type_last = region_type;


		std::vector<Inexact_Plane> planes_2_last = planes_2;
		std::vector<int> inliers_to_planes_last = inliers_to_planes;




		times++;
		number_moving = 0;
		propagation1 = true;


		for (size_t i = 0; i < points.size(); ++i) {
			if (inliers_to_planes[i] == -1) continue;

			if (if_moved[i] >= 5) continue;
			if (update_primitive_size[inliers_to_planes[i]] <= min_points) continue;


			if ((int)spherical_neighborhood[i].size() < Nb_neigh) {
				Nb_neigh = (int)spherical_neighborhood[i].size();
			}

			std::set<int> one_neight_id;




			//check if this point is on the boundary between primitives
			for (int it = 0; it < Nb_neigh; it++) {

				int neighbor_index = spherical_neighborhood[i][it];


				if (inliers_to_planes[neighbor_index] != -1 && inliers_to_planes[i] != inliers_to_planes[neighbor_index]) {

					one_neight_id.insert(inliers_to_planes[neighbor_index]);

				}


			}

			if (one_neight_id.empty()) { continue; }




			Inexact_Point_3 this_point = points[i].first;

			double min_distance = sqrt(CGAL::squared_distance(planes_2[inliers_to_planes[i]], points[i].first));

			int changed_plane_id = -1;
			//transfer the point to its best primitive
			for (int neight_id : one_neight_id) {



				double neight_distance = sqrt(CGAL::squared_distance(planes_2[neight_id], points[i].first));

				if (neight_distance < min_distance) {

					if (abs(points[i].second * planes_2[neight_id].orthogonal_vector()) >= normal_threshold) {
						min_distance = neight_distance;
						changed_plane_id = neight_id;
					}
				}
			}

			if (changed_plane_id != -1) {
				update_primitive_size[inliers_to_planes[i]]--;
				inliers_to_planes[i] = changed_plane_id;
				update_primitive_size[changed_plane_id]++;
				if_moved[i] ++;
				number_moving++;
			}




		}

		int plane_size = planes_to_inliers.size();
		planes_to_inliers.clear();
		planes_to_inliers.resize(plane_size);
		for (int k = 0; k < points.size(); k++) {
			if (inliers_to_planes[k] != -1) {
				planes_to_inliers[inliers_to_planes[k]].push_back(k);
			}

		}




		planes_2.clear();



		for (int ii = 0; ii < planes_to_inliers.size(); ++ii) {


			std::vector<Inexact_Point_3> inliers_i;
			inliers_i.reserve(planes_to_inliers[ii].size());

			for (int jj = 0; jj < planes_to_inliers[ii].size(); ++jj) {
				const Inexact_Point_3 & pt = points[planes_to_inliers[ii][jj]].first;
				inliers_i.push_back(pt);

			}



			if (inliers_i.size() < 3) {

		
				getchar();
			}
			Inexact_Plane plane;
			linear_least_squares_fitting_3(inliers_i.begin(), inliers_i.end(), plane, CGAL::Dimension_tag<0>());

			planes_2.push_back(plane);

		}


		get_distance_diviation();
		//make sure that the transfer operation reducing the energy.
		if (mean_distance_diaviation > mean_distance_diaviation_before_transfer) {
			propagation1 = false;
			planes_2 = planes_2_last;
			planes_to_inliers = planes_to_inliers_last;
			region_type = region_type_last;
			inliers_to_planes = inliers_to_planes_last;


		}



		if (if_constraint) {

			if (mean_distance_diaviation > ori_mean_error) {
				propagation1 = false;
				planes_2 = planes_2_last;
				planes_to_inliers = planes_to_inliers_last;
				region_type = region_type_last;
				inliers_to_planes = inliers_to_planes_last;


			}
		}
		t_l++;

		if (number_moving == 0 || times > 100) { propagation1 = false; }
	} while (propagation1);







}




void Shape_Detector::test_connected_primitives() {
	std::vector<int> one_connect;
	one_connect = std::vector<int>(planes_to_inliers.size(), 0);
	primitive_connection = std::vector<std::vector<int>>(planes_to_inliers.size(), one_connect);
	int Nb_neigh = 10;
	for (size_t i = 0; i < points.size(); i++) {
		if ((int)spherical_neighborhood[i].size() < Nb_neigh) {
			Nb_neigh = (int)spherical_neighborhood[i].size();
		}

		std::set<int> one_neight_id;




		for (int it = 0; it < Nb_neigh; it++) {

			int neighbor_index = spherical_neighborhood[i][it];

			if (inliers_to_planes[neighbor_index] != -1 && inliers_to_planes[i] != -1 && inliers_to_planes[i] != inliers_to_planes[neighbor_index]) {
				primitive_connection[inliers_to_planes[neighbor_index]][inliers_to_planes[i]] ++;
				primitive_connection[inliers_to_planes[i]][inliers_to_planes[neighbor_index]] ++;

			}


		}

	}

	for (int i = 0; i < primitive_connection.size(); ++i) {
		std::vector<int> this_raw = primitive_connection[i];
		for (int j = i; j < this_raw.size(); ++j) {
			if (this_raw[j] >= 10) {
				primitive_connection[i][j] = 1;
				primitive_connection[j][i] = 1;

			}
			else {
				primitive_connection[i][j] = 0;
				primitive_connection[j][i] = 0;


			}

		}

	}
}



double Shape_Detector::merge_distance_changed_with_epsilon(int i, int j, double & dif, std::vector<int> & move_ids) {
	Inexact_Plane p_i = planes_2[i];
	Inexact_Plane p_j = planes_2[j];
	double dis_before = 0;

	std::vector<Inexact_Point_3> assigned_pts;
	std::vector<int> assigned_pts_id;

	assigned_pts.reserve(planes_to_inliers[i].size() + planes_to_inliers[j].size());
	for (int pt_index : planes_to_inliers[i]) {
		assigned_pts.push_back(points[pt_index].first);
		assigned_pts_id.push_back(pt_index);

		dis_before += sqrt(CGAL::squared_distance(p_i, points[pt_index].first));
	}

	for (int pt_index : planes_to_inliers[j]) {
		assigned_pts_id.push_back(pt_index);

		assigned_pts.push_back(points[pt_index].first);
		dis_before += sqrt(CGAL::squared_distance(p_j, points[pt_index].first));

	}





	size_t n_poly = assigned_pts.size();



	double dP_norm = 0;

	Inexact_Plane test_plane;

	linear_least_squares_fitting_3(assigned_pts.begin(), assigned_pts.end(), test_plane, CGAL::Dimension_tag<0>());



	Inexact_Plane new_plane;
	std::vector<Inexact_Point_3> new_assigned_pts;
	std::vector<int> new_assigned_pts_ids;
	for (int pt : assigned_pts_id) {
		if (sqrt(CGAL::squared_distance(test_plane, points[pt].first)) > epsilon) {
			move_ids.push_back(pt);
		}
		else {
			new_assigned_pts.push_back(points[pt].first);
			new_assigned_pts_ids.push_back(pt);

		}
	}

	if (new_assigned_pts_ids.size() < min_points) {
		dif = FLT_MAX;
		return false;
	}
	linear_least_squares_fitting_3(new_assigned_pts.begin(), new_assigned_pts.end(), new_plane, CGAL::Dimension_tag<0>());
	for (int p : new_assigned_pts_ids) {
		dP_norm += sqrt(CGAL::squared_distance(new_plane, points[p].first));

	}

	//dP_norm /= double(n_poly);
	//dP_norm = sqrt(dP_norm);
	dif = dP_norm - dis_before;
	//	return (dif)/double(number_of_assigned_points);
		//return (dP_norm - dis_before) / double(ori_inliers_number);
	return dif;




}


double Shape_Detector::merge_distance_changed_normal_with_epsilon(int i, int j, double & dif, std::vector<int> & move_ids) {
	Inexact_Plane p_i = planes_2[i];
	Inexact_Plane p_j = planes_2[j];
	double dis_before = 0;

	std::vector<Inexact_Point_3> assigned_pts;
	std::vector<int> assigned_pts_id;

	assigned_pts.reserve(planes_to_inliers[i].size() + planes_to_inliers[j].size());
	for (int pt_index : planes_to_inliers[i]) {
		assigned_pts.push_back(points[pt_index].first);
		assigned_pts_id.push_back(pt_index);
		dis_before += abs(points[pt_index].second * p_i.orthogonal_vector());
	}

	for (int pt_index : planes_to_inliers[j]) {
		assigned_pts.push_back(points[pt_index].first);
		assigned_pts_id.push_back(pt_index);

		dis_before += abs(points[pt_index].second * p_j.orthogonal_vector());


	}





	size_t n_poly = assigned_pts.size();
	//dis_before /= double(n_poly);


	double dP_norm = 0;

	Inexact_Plane test_plane;
	linear_least_squares_fitting_3(assigned_pts.begin(), assigned_pts.end(), test_plane, CGAL::Dimension_tag<0>());
	/*for (int id_in : assigned_pts_id) {

		dP_norm += abs(points[id_in].second * test_plane.orthogonal_vector());

		move_ids.clear();
	}*/


	Inexact_Plane new_plane;
	std::vector<Inexact_Point_3> new_assigned_pts;
	std::vector<int> new_assigned_pts_ids;
	for (int pt : assigned_pts_id) {
		if (sqrt(CGAL::squared_distance(test_plane, points[pt].first)) > epsilon) {
			move_ids.push_back(pt);
		}
		else {
			new_assigned_pts.push_back(points[pt].first);
			new_assigned_pts_ids.push_back(pt);

		}
	}
	if (new_assigned_pts_ids.size() < min_points) {
		dif = -FLT_MAX;

		return false;
	}
	linear_least_squares_fitting_3(new_assigned_pts.begin(), new_assigned_pts.end(), new_plane, CGAL::Dimension_tag<0>());
	for (int p : new_assigned_pts_ids) {

		dP_norm += abs(points[p].second * new_plane.orthogonal_vector());

	}


	//dP_norm /= double(n_poly);
	//dP_norm = sqrt(dP_norm);

	dif = dP_norm - dis_before;
	//	return (dif)/double(number_of_assigned_points);
		//return (dP_norm - dis_before) / double(ori_inliers_number);
	return dif;




}




bool Shape_Detector::convert_form_merge_l2(int i, int j, std::pair<std::pair<int, std::vector<int>>, std::vector<double>>& one_element) {

	std::vector<int> o;
	o.push_back(i);
	o.push_back(j);
	double dif_m = 0;
	std::vector<int> merge_moved_ids;
	//we suppose that the pair of primitives can not merged,when their normal cosine smaller than normal_threshold.
	if (abs(planes_2[i].orthogonal_vector()*planes_2[j].orthogonal_vector()) < normal_threshold) {
		return false;
	}
	std::vector<double> s_d;
	double all_changed_distances = merge_distance_changed_with_epsilon(i, j, dif_m, merge_moved_ids);

	/*if (merge_moved_ids.size() > points.size()*0.01|| merge_moved_ids.size()>planes_to_inliers[i].size()/4|| merge_moved_ids.size() > planes_to_inliers[j].size() / 4) {
		return false;
	}*/

	if ((planes_to_inliers[i].size() + planes_to_inliers[j].size() - merge_moved_ids.size()) < min_points) {
		return false;
	}
	double change_ccc = energy_changed_merge(all_changed_distances, -1.0, merge_moved_ids.size());
	if (change_ccc >= 0) {
		return false;
	}

	s_d.push_back(change_ccc);
	o.insert(o.end(), merge_moved_ids.begin(), merge_moved_ids.end());

	s_d.push_back(dif_m);
	one_element = std::make_pair(std::make_pair(1, o), s_d);
	return true;
}

bool Shape_Detector::convert_form_split_l2(int i, std::vector<std::vector<int>>& max_list_vector, std::vector<std::vector<int>>& min_list_vector, std::pair<std::pair<int, std::vector<int>>, std::vector<double>> & one_element) {
	std::vector<int> max_list;
	std::vector<int> min_list;
	double dif_s = 0;
	bool if_splite = separate_two_out(i, max_list, min_list, dif_s);

	if (energy_changed(dif_s, 1.0) >= 0 || max_list.size() < min_points || min_list.size() < min_points || !if_splite) {
		return false;
	}
	max_list_vector.push_back(max_list);
	min_list_vector.push_back(min_list);

	std::vector<int> o;
	o.push_back(i);
	o.push_back(max_list_vector.size() - 1);
	std::vector<double> s_d;
	s_d.push_back(energy_changed(dif_s, 1.0));
	s_d.push_back(dif_s);
	one_element = std::make_pair(std::make_pair(2, o), s_d);
	return true;
}

bool Shape_Detector::convert_form_exclude_l2(int i, std::pair<std::pair<int, std::vector<int>>, std::vector<double>>& one_element) {

	std::vector<int> o;
	o.push_back(bad_points_shape[i].first.first);

	o.insert(o.end(), bad_points_shape[i].first.second.begin(), bad_points_shape[i].first.second.end());

	one_element = std::make_pair(std::make_pair(3, o), bad_points_shape[i].second);

	return true;
}

bool Shape_Detector::convert_form_insert_l2(int i, std::pair<std::pair<int, std::vector<int>>, std::vector<double>>& one_element) {
	std::vector<int> o;
	o.push_back(good_points_shape[i].first.first);
	o.insert(o.end(), good_points_shape[i].first.second.begin(), good_points_shape[i].first.second.end());


	one_element = std::make_pair(std::make_pair(4, o), good_points_shape[i].second);
	return true;
}
bool Shape_Detector::update_bad_points(int id_shape, std::pair<std::pair<int, std::vector<int>>, std::vector<double>>& one_element) {



	std::vector<std::pair<int, double>> ids_bad_inliers;
	double bigest_normal = mean_normal_diviation;
	double bigest_dis = 1.5 * mean_distance_diaviation;


	const Inexact_Plane & H = planes_2[id_shape];


	for (int j : planes_to_inliers[id_shape]) {
		const Point_with_normal & pt = points[j];
		if (points_if_added[j] > 2) continue;

		if (sqrt(CGAL::squared_distance(H, pt.first)) > bigest_dis) {
			ids_bad_inliers.push_back(std::make_pair(j, sqrt(CGAL::squared_distance(H, pt.first))));
		}

	}

	std::priority_queue<std::pair<int, double>, std::vector<std::pair<int, double>>, Remove_Comparator> p_q;
	for (std::pair<int, double> ppp : ids_bad_inliers) {
		p_q.push(ppp);

	}
	std::vector<int> removed_points_ids;
	double en = 0;
	double dif_e = 0;

	double n_b = 0;
	while (!p_q.empty() && n_b < number_of_insert_exclude) {
		removed_points_ids.push_back(p_q.top().first);
		p_q.pop();
		n_b++;
	}
	dif_e = remove_changed_error(id_shape, removed_points_ids);
	en = energy_changed_second(dif_e, n_b);

	if (removed_points_ids.size() != 0 && (planes_to_inliers[id_shape].size() - removed_points_ids.size()) >= min_points && en < 0) {
		std::vector<double> dd;
		dd.push_back(en);
		dd.push_back(dif_e);
		std::vector<int> o;
		o.push_back(id_shape);
		o.insert(o.end(), removed_points_ids.begin(), removed_points_ids.end());


		one_element = std::make_pair(std::make_pair(3, o), dd);
		return true;
	}
	else { return false; }



}
bool Shape_Detector::update_good_points(int id_shape, std::pair<std::pair<int, std::vector<int>>, std::vector<double>>& one_element) {

	int Nb_neigh = 10;
	std::vector<std::pair<int, double>> ids_good_outliers;
	for (int i = 0; i < points.size(); ++i) {
		if (if_insert_candidate[i] == 0) continue;
		if (inliers_to_planes[i] != -1) continue;
		if (points_if_added[i] > 2) continue;//avoide infinite loop
		if ((int)spherical_neighborhood[i].size() < Nb_neigh) {
			Nb_neigh = (int)spherical_neighborhood[i].size();
		}

		std::set<int> one_neight_id;

		bool if_correspond = false;
		for (int it = 0; it < Nb_neigh; ++it) {

			int neighbor_index = spherical_neighborhood[i][it];

			if (inliers_to_planes[neighbor_index] != -1) {

				one_neight_id.insert(inliers_to_planes[neighbor_index]);
				if (inliers_to_planes[neighbor_index] == id_shape) {
					if_correspond = true;
				}
			}


		}

		if (one_neight_id.empty() || !if_correspond) { continue; }

		double mm_dis = epsilon;

		int changed_plane_id = -1;
		Point_with_normal this_p = points[i];

		for (int neight_id : one_neight_id) {
			if (abs(this_p.second * planes_2[neight_id].orthogonal_vector()) > normal_threshold) {
				if (sqrt(CGAL::squared_distance(planes_2[neight_id], this_p.first)) < mm_dis) {

					mm_dis = sqrt(CGAL::squared_distance(planes_2[neight_id], this_p.first));
					changed_plane_id = neight_id;

				}
			}
		}


		if (changed_plane_id == id_shape) {
			ids_good_outliers.push_back(std::make_pair(i, mm_dis));//point_id, mean distance.

		}

	}


	std::priority_queue<std::pair<int, double>, std::vector<std::pair<int, double>>, Add_Comparator> p_q;

	for (std::pair<int, double> ppp : ids_good_outliers) {
		p_q.push(ppp);
	}

	std::vector<int> added_points_ids;
	double en = 0;
	double dif_e = 0;
	double n_a = 0;

	while (!p_q.empty() && n_a < number_of_insert_exclude) {
		added_points_ids.push_back(p_q.top().first);
		p_q.pop();
		n_a++;
	}
	dif_e = add_changed_error(id_shape, added_points_ids);
	en = energy_changed_second(dif_e, -n_a);

	if (added_points_ids.size() != 0 && en < 0) {
		for (int idd : added_points_ids) {
			if_insert_candidate[idd] = 0;
		}
		std::vector<double> dd;
		dd.push_back(en);
		dd.push_back(dif_e);

		std::vector<int> o;
		o.push_back(id_shape);
		o.insert(o.end(), added_points_ids.begin(), added_points_ids.end());


		one_element = std::make_pair(std::make_pair(4, o), dd);


		return true;

	}
	else {
		return false;
	}


}



/*
Please refer to local_operators_normal() for comment.
*/
void Shape_Detector::local_operators() {
	//points_changed = std::vector<int>(points.size(), 0);
	bool bb = false;
	all_t_merge += t_merge;
	all_t_split += t_split;
	all_t_insert += t_insert;
	all_t_exlude += t_exlude;
	t_m = 0;
	t_merge = 0;
	t_split = 0;
	t_exlude = 0;
	t_insert = 0;





	std::priority_queue<std::pair<std::pair<int, std::vector<int>>, std::vector<double>>, std::vector<std::pair<std::pair<int, std::vector<int>>, std::vector<double>>>, Weight_Comparator_with_energy> p_q;

	std::vector<std::pair<std::pair<int, std::vector<int>>, std::vector<double>>> v_p_q;

	std::vector<std::vector<int>> max_list_vector;
	std::vector<std::vector<int>> min_list_vector;


	for (int i = 0; i < primitive_connection.size(); ++i) {

		std::vector<int> one = primitive_connection[i];
		for (int j = i + 1; j < one.size(); ++j) {
			if (one[j] >= 1) {


				std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;
				bool if_satisfied = convert_form_merge_l2(i, j, one_element);
				if (!if_satisfied) continue;

				p_q.push(one_element);
				v_p_q.push_back(one_element);

			}
		}
	}


	for (int i = 0; i < planes_to_inliers.size(); ++i) {

		if (region_type[i] > 5) continue;
		std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;
		bool if_satisfied = convert_form_split_l2(i, max_list_vector, min_list_vector, one_element);
		if (!if_satisfied) continue;
		p_q.push(one_element);
		v_p_q.push_back(one_element);




	}


	for (int i = 0; i < bad_points_shape.size(); ++i) {




		std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;
		bool if_satisfied = convert_form_exclude_l2(i, one_element);
		if (!if_satisfied) continue;
		p_q.push(one_element);
		v_p_q.push_back(one_element);

	}


	for (int i = 0; i < good_points_shape.size(); ++i) {

		


		std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;
		bool if_satisfied = convert_form_insert_l2(i, one_element);
		if (!if_satisfied) continue;
	
		p_q.push(one_element);
		v_p_q.push_back(one_element);

	}








	do {

		bb = true;

		std::vector<int> if_merged;
		if_merged = std::vector<int>(planes_to_inliers.size(), -1);



		std::vector<int> one_merged;

		std::vector<int> max_list;
		std::vector<int> min_list;
		std::vector<int> id_of_moved_merge;


		int one_splite = -1;
		std::vector<int> one_remove;
		int one_remove_plane = -1;
		std::vector<int> one_add;
		int one_add_plane = -1;

		while (!p_q.empty()) {


			if (p_q.top().second[0] > 0) break;
			if (p_q.top().first.first == 1) {

				std::vector<int> this_pair = p_q.top().first.second;

				int f = this_pair[0];
				int s = this_pair[1];

				if ((planes_to_inliers[f].size() + planes_to_inliers[s].size() - 2 - this_pair.size()) < min_points) {
					p_q.pop();
					continue;
				}



				double after_operation_distance_visual = all_distance_diaviation + p_q.top().second[1];

				if (if_constraint) {
					/*if (abs(planes_2[f].orthogonal_vector() * planes_2[s].orthogonal_vector()) < normal_threshold) {

						p_q.pop();
						continue;

					}*/
					if (number_of_assigned_points + 2 - this_pair.size() <= ori_inliers_number) {
						p_q.pop();
						continue;
					}

					if (after_operation_distance_visual / double(number_of_assigned_points) > ori_mean_error) {

						p_q.pop();
						continue;
					}
					else {
						p_q.pop();
						id_of_moved_merge.clear();
						if (this_pair.size() > 2) {
							for (int k = 2; k < this_pair.size(); ++k) {
								id_of_moved_merge.push_back(this_pair[k]);
							}
						}
						number_of_assigned_points = number_of_assigned_points - id_of_moved_merge.size();
						if_merged[f] = 0;
						if_merged[s] = 0;
						one_merged.push_back(f);
						one_merged.push_back(s);
						//after_operation_distance = after_operation_distance_visual;
						all_distance_diaviation = after_operation_distance_visual;

						t_m++;
						t_merge++;
					

						break;
					}
				}
				else {
					/*if (abs(planes_2[f].orthogonal_vector() * planes_2[s].orthogonal_vector()) < normal_threshold) {

						p_q.pop();
						continue;

					}*/
					if ((planes_to_inliers[f].size() + planes_to_inliers[s].size() + 2 - this_pair.size()) < min_points) {//when the new primitives have less point than threshold.
						p_q.pop();
						continue;
					}
					else {
						p_q.pop();
						id_of_moved_merge.clear();
						if (this_pair.size() > 2) {
							for (int k = 2; k < this_pair.size(); ++k) {
								id_of_moved_merge.push_back(this_pair[k]);
							}
						}
						number_of_assigned_points = number_of_assigned_points - id_of_moved_merge.size();

						if_merged[f] = 0;
						if_merged[s] = 0;
						one_merged.push_back(f);
						one_merged.push_back(s);

						all_distance_diaviation = after_operation_distance_visual;

						t_m++;
						t_merge++;
				

						break;
					}
				}

			}
			else if (p_q.top().first.first == 2) {



				int splited_index = p_q.top().first.second[0];
				int info_index = p_q.top().first.second[1];


				double after_splited_distance = all_distance_diaviation + p_q.top().second[1];



				max_list = max_list_vector[info_index];
				min_list = min_list_vector[info_index];

				if (max_list.size() < min_points || min_list.size() < min_points) {
					p_q.pop();
					continue;
				}
				if (if_constraint) {

					if (after_splited_distance / double(number_of_assigned_points) > ori_mean_error || (planes_2.size() + 1) > ori_primitives_number) {

						p_q.pop();
						continue;

					}
					else {
						p_q.pop();
						if_merged[splited_index] = 0;

						one_splite = splited_index;

						all_distance_diaviation = after_splited_distance;

						t_m++;
						t_split++;
						
						break;
					}
				}
				else {

					p_q.pop();
					if_merged[splited_index] = 0;

					one_splite = splited_index;

					all_distance_diaviation = after_splited_distance;

					t_m++;
					t_split++;
				
					break;
				}

			}
			else if (p_q.top().first.first == 3) {
				

				int plane_id = p_q.top().first.second[0];

				//double after_delete_distance = after_operation_distance + p_q.top().second[1];
				double after_delete_distance = all_distance_diaviation + p_q.top().second[1];
				if (if_constraint) {
					if (number_of_assigned_points + 1 - p_q.top().first.second.size() <= ori_inliers_number || after_delete_distance / double(number_of_assigned_points + 1 - p_q.top().first.second.size()) >= ori_mean_error) {

						//if (after_operation_inlier_number + 1 - p_q.top().first.second.size() <= ori_inliers_number || after_delete_distance / double(after_operation_inlier_number + 1 - p_q.top().first.second.size()) >= ori_mean_error) {

						p_q.pop();
						continue;
					}
				}


				if (planes_to_inliers[plane_id].size() + 1 - p_q.top().first.second.size() < min_points) {

					p_q.pop();
					continue;
				}
				else {

					if_merged[plane_id] = 0;
					one_remove_plane = plane_id;
					std::vector<int> points_ids = p_q.top().first.second;
					points_ids.erase(points_ids.begin());
					one_remove = points_ids;
					//after_operation_distance = after_delete_distance;
					all_distance_diaviation = after_delete_distance;

					number_of_assigned_points = number_of_assigned_points - one_remove.size();
					//after_operation_inlier_number = after_operation_inlier_number - one_remove.size();
					p_q.pop();
					t_m++;
					t_exlude++;
				
					break;
				}

			}
			else if (p_q.top().first.first == 4) {
				

				int plane_id = p_q.top().first.second[0];

				//double after_add_distance = after_operation_distance + p_q.top().second[1];
				double after_add_distance = all_distance_diaviation + p_q.top().second[1];
				if (if_constraint) {
					if (after_add_distance / double(number_of_assigned_points - 1 + p_q.top().first.second.size()) >= ori_mean_error) {
						//if (after_add_distance / double(after_operation_inlier_number - 1 + p_q.top().first.second.size()) >= ori_mean_error) {

						p_q.pop();
						continue;
					}
				}
				if (planes_to_inliers[plane_id].size() - 1 + p_q.top().first.second.size() < min_points) {

					p_q.pop();
				}
				else {

					if_merged[plane_id] = 0;
					one_add_plane = plane_id;
					std::vector<int> points_ids = p_q.top().first.second;
					points_ids.erase(points_ids.begin());
					one_add = points_ids;
					//after_operation_distance = after_add_distance;
					all_distance_diaviation = after_add_distance;

					number_of_assigned_points = number_of_assigned_points + one_add.size();
					//after_operation_inlier_number = after_operation_inlier_number + one_add.size();
					p_q.pop();
					t_m++;
					t_insert++;
				
					break;
				}

			}
		}

		if (one_merged.size() != 0 && one_splite == -1 && one_remove.size() == 0 && one_add.size() == 0) {


			int last_size = planes_2.size();
			std::vector<int> respective_planes = std::vector<int>(last_size, -1);


			int id_1 = one_merged[0];
			int id_2 = one_merged[1];
			std::vector<int> one_merge_points_real;
			std::vector<int> one_merge_points = planes_to_inliers[id_1];
			one_merge_points.insert(one_merge_points.end(), planes_to_inliers[id_2].begin(), planes_to_inliers[id_2].end());
			for (int id : one_merge_points) {
				if (std::find(id_of_moved_merge.begin(), id_of_moved_merge.end(), id) == id_of_moved_merge.end()) {
					one_merge_points_real.push_back(id);
				}
			}

			one_merge_points = one_merge_points_real;
			for (int m_id : id_of_moved_merge) {
				inliers_to_planes[m_id] = -1;
			}


			if (one_merge_points.size() >= min_points) {
				if (id_2 < id_1) {
					for (int t = 0; t < last_size; ++t) {
						if (t < id_2) {
							respective_planes[t] = t;
						}
						else if (t == id_2) {
							respective_planes[t] = last_size - 2;
						}
						else if (t < id_1) {
							respective_planes[t] = t - 1;

						}
						else if (t == id_1) {
							respective_planes[t] = last_size - 2;

						}
						else {
							respective_planes[t] = t - 2;

						}


					}
					region_type.push_back(region_type[id_1] + region_type[id_2] + 1);
					region_type.erase(region_type.begin() + id_1);
					region_type.erase(region_type.begin() + id_2);



					planes_to_inliers.erase(planes_to_inliers.begin() + id_1);
					planes_to_inliers.erase(planes_to_inliers.begin() + id_2);

					planes_2.erase(planes_2.begin() + id_1);
					planes_2.erase(planes_2.begin() + id_2);


				}
				else {
					for (int t = 0; t < last_size; ++t) {
						if (t < id_1) {
							respective_planes[t] = t;
						}
						else if (t == id_1) {
							respective_planes[t] = last_size - 2;
						}
						else if (t < id_2) {
							respective_planes[t] = t - 1;

						}
						else if (t == id_2) {
							respective_planes[t] = last_size - 2;

						}
						else {
							respective_planes[t] = t - 2;

						}


					}
					region_type.push_back(region_type[id_1] + region_type[id_2] + 1);
					region_type.erase(region_type.begin() + id_2);
					region_type.erase(region_type.begin() + id_1);

					planes_to_inliers.erase(planes_to_inliers.begin() + id_2);
					planes_to_inliers.erase(planes_to_inliers.begin() + id_1);

					planes_2.erase(planes_2.begin() + id_2);
					planes_2.erase(planes_2.begin() + id_1);


				}


				planes_to_inliers.push_back(one_merge_points);


				std::vector<Inexact_Point_3> inliers;
				inliers.reserve(one_merge_points.size());

				for (int jj = 0; jj < one_merge_points.size(); ++jj) {

					inliers.push_back(points[one_merge_points[jj]].first);

				}
				Inexact_Plane plane;
				if (inliers.size() < 3) {

					continue;
				}
				linear_least_squares_fitting_3(inliers.begin(), inliers.end(), plane, CGAL::Dimension_tag<0>());

				planes_2.push_back(plane);

				std::vector<int> inliers_to_planes_merge_local;
				inliers_to_planes_merge_local = std::vector<int>(points.size(), -1);
				for (int m = 0; m < planes_to_inliers.size(); ++m) {
					for (int k = 0; k < planes_to_inliers[m].size(); ++k) {
						inliers_to_planes_merge_local[planes_to_inliers[m][k]] = m;
					}


				}
				inliers_to_planes.clear();
				inliers_to_planes = inliers_to_planes_merge_local;



				std::vector<std::pair<std::pair<int, std::vector<int>>, std::vector<double>>> local_v_p_q;
				while (!p_q.empty()) p_q.pop();
				for (std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_p_q : v_p_q) {
					if (one_p_q.first.first == 1) {//when the operator is merge.
						std::vector<int> this_pair = one_p_q.first.second;


						if (if_merged[this_pair[0]] == 0 && if_merged[this_pair[1]] == 0)//when that is the merged pair.
						{

							continue;
						}
						else if ((if_merged[this_pair[0]] == 0 && if_merged[this_pair[1]] != 0) || (if_merged[this_pair[0]] != 0 && if_merged[this_pair[1]] == 0)) {//when only one primitive is merged.



							std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;
							bool if_satisfied = convert_form_merge_l2(respective_planes[this_pair[1]], respective_planes[this_pair[0]], one_element);
							if (!if_satisfied) continue;
							p_q.push(one_element);
							local_v_p_q.push_back(one_element);



						}

						else {



							std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;

							std::vector<int> o = this_pair;
							o[0] = respective_planes[this_pair[0]];
							o[1] = respective_planes[this_pair[1]];
							one_element = std::make_pair(std::make_pair(1, o), one_p_q.second);
							p_q.push(one_element);
							local_v_p_q.push_back(one_element);




						}
					}
					else if (one_p_q.first.first == 2) {
						std::vector<int> this_reg = one_p_q.first.second;
						if (if_merged[this_reg[0]] == 0)//when the primitive has been merged.
						{

							continue;
						}

						else {
							std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;

							std::vector<int> o;
							o.push_back(respective_planes[this_reg[0]]);
							o.push_back(this_reg[1]);
							one_element = std::make_pair(std::make_pair(2, o), one_p_q.second);
							p_q.push(one_element);
							local_v_p_q.push_back(one_element);


						}

					}
					else if (one_p_q.first.first == 3) {
						std::vector<int> this_p_p = one_p_q.first.second;
						if (if_merged[this_p_p[0]] == 0) {
							continue;
						}
						else {
							std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;

							std::vector<int> o;
							o = this_p_p;
							o[0] = respective_planes[o[0]];

							one_element = std::make_pair(std::make_pair(3, o), one_p_q.second);
							p_q.push(one_element);
							local_v_p_q.push_back(one_element);


						}

					}
					else if (one_p_q.first.first == 4) {
						std::vector<int> this_p_p = one_p_q.first.second;
						if (if_merged[this_p_p[0]] == 0) {
							continue;
						}
						else {
							std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;

							std::vector<int> o;
							o = this_p_p;
							o[0] = respective_planes[o[0]];

							one_element = std::make_pair(std::make_pair(4, o), one_p_q.second);
							p_q.push(one_element);
							local_v_p_q.push_back(one_element);


						}

					}
				}

				//add new operator, split the new merged primitive
				std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element_1;
				bool if_satisfied = convert_form_split_l2(last_size - 2, max_list_vector, min_list_vector, one_element_1);
				if (if_satisfied) {
					p_q.push(one_element_1);
					local_v_p_q.push_back(one_element_1);
				}
				//add new operat, exlude and inerte for the new merged primitive
				if (all_t_transfer > 0) {
					std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element_b;
					bool if_satisfied_b = update_bad_points(last_size - 2, one_element_b);
					if (if_satisfied_b) {
						p_q.push(one_element_b);
						local_v_p_q.push_back(one_element_b);
					}
				}

				std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element_g;
				bool if_satisfied_g = update_good_points(last_size - 2, one_element_g);
				if (if_satisfied_g) {
					p_q.push(one_element_g);
					local_v_p_q.push_back(one_element_g);
				}

				v_p_q.clear();
				v_p_q = local_v_p_q;


			}

		}
		else if (one_merged.size() == 0 && one_splite != -1 && one_remove.size() == 0 && one_add.size() == 0) {
			int last_size = planes_2.size();

			std::vector<int> respective_planes = std::vector<int>(last_size, -1);


			for (int t = 0; t < last_size; ++t) {
				if (t < one_splite) {
					respective_planes[t] = t;
				}
				else if (t == one_splite) {
					respective_planes[t] = last_size;
				}

				else {
					respective_planes[t] = t - 1;

				}


			}
			planes_to_inliers.erase(planes_to_inliers.begin() + one_splite);
			region_type.push_back(region_type[one_splite]);
			region_type.push_back(region_type[one_splite]);
			region_type.erase(region_type.begin() + one_splite);




			planes_2.erase(planes_2.begin() + one_splite);







			planes_to_inliers.push_back(max_list);
			planes_to_inliers.push_back(min_list);




			std::vector<Inexact_Point_3> max_inliers;
			max_inliers.reserve(max_list.size());

			for (int jj = 0; jj < max_list.size(); ++jj) {
				//points_changed[one_merge_points[jj]] = -1;
				//const Inexact_Point_3 & pt = points[max_list[jj]].first;
				max_inliers.push_back(points[max_list[jj]].first);

			}
			Inexact_Plane max_plane;
			if (max_inliers.size() < 3) {
				
				getchar();
			}
			linear_least_squares_fitting_3(max_inliers.begin(), max_inliers.end(), max_plane, CGAL::Dimension_tag<0>());

			planes_2.push_back(max_plane);


			std::vector<Inexact_Point_3> min_inliers;
			min_inliers.reserve(min_list.size());

			for (int jj = 0; jj < min_list.size(); ++jj) {
				//points_changed[one_merge_points[jj]] = -1;
				//const Inexact_Point_3 & pt = points[min_list[jj]].first;
				min_inliers.push_back(points[min_list[jj]].first);

			}
			Inexact_Plane min_plane;
			if (min_inliers.size() < 3) {
			
				getchar();
			}
			linear_least_squares_fitting_3(min_inliers.begin(), min_inliers.end(), min_plane, CGAL::Dimension_tag<0>());

			planes_2.push_back(min_plane);

			std::vector<int> inliers_to_planes_merge_local;
			inliers_to_planes_merge_local = std::vector<int>(points.size(), -1);
			for (int m = 0; m < planes_to_inliers.size(); ++m) {
				for (int k = 0; k < planes_to_inliers[m].size(); ++k) {
					inliers_to_planes_merge_local[planes_to_inliers[m][k]] = m;
				}


			}
			inliers_to_planes.clear();
			inliers_to_planes = inliers_to_planes_merge_local;

			std::vector<std::pair<std::pair<int, std::vector<int>>, std::vector<double>>> local_v_p_q;
			while (!p_q.empty()) p_q.pop();
			for (std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_p_q : v_p_q) {
				if (one_p_q.first.first == 1) {
					std::vector<int> this_pair = one_p_q.first.second;
					//double this_weight = p_q_c.top().second;
					if (if_merged[this_pair[0]] == 0 && if_merged[this_pair[1]] == 0)//not possible case
					{
						continue;
					}
					else if (if_merged[this_pair[0]] == 0 && if_merged[this_pair[1]] != 0) {//one of the primitive is splited
						if (test_if_connected(respective_planes[this_pair[1]], planes_to_inliers.size() - 1)) {
							std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;
							bool if_satisfied = convert_form_merge_l2(planes_to_inliers.size() - 1, respective_planes[this_pair[1]], one_element);
							if (!if_satisfied) continue;

							p_q.push(one_element);
							local_v_p_q.push_back(one_element);




						}
						if (test_if_connected(respective_planes[this_pair[1]], planes_to_inliers.size() - 2)) {
							std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;
							bool if_satisfied = convert_form_merge_l2(planes_to_inliers.size() - 2, respective_planes[this_pair[1]], one_element);
							if (!if_satisfied) continue;


							p_q.push(one_element);
							local_v_p_q.push_back(one_element);


						}



					}
					else if (if_merged[this_pair[0]] != 0 && if_merged[this_pair[1]] == 0) {
						if (test_if_connected(respective_planes[this_pair[0]], planes_to_inliers.size() - 1)) {
							std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;
							bool if_satisfied = convert_form_merge_l2(planes_to_inliers.size() - 1, respective_planes[this_pair[0]], one_element);
							if (!if_satisfied) continue;

							p_q.push(one_element);
							local_v_p_q.push_back(one_element);


						}
						if (test_if_connected(respective_planes[this_pair[0]], planes_to_inliers.size() - 2)) {
							std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;
							bool if_satisfied = convert_form_merge_l2(planes_to_inliers.size() - 2, respective_planes[this_pair[0]], one_element);
							if (!if_satisfied) continue;

							p_q.push(one_element);
							local_v_p_q.push_back(one_element);


						}


					}
					else {
						std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;

						std::vector<int> o = this_pair;
						o[0] = respective_planes[this_pair[0]];
						o[1] = respective_planes[this_pair[1]];
						one_element = std::make_pair(std::make_pair(1, o), one_p_q.second);
						p_q.push(one_element);
						local_v_p_q.push_back(one_element);

					}
				}
				else if (one_p_q.first.first == 2) {
					std::vector<int> this_reg = one_p_q.first.second;
					//double this_weight = p_q_c.top().second;
					if (if_merged[this_reg[0]] == 0)
					{

						continue;
					}

					else {
						std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;

						std::vector<int> o;
						o.push_back(respective_planes[this_reg[0]]);
						o.push_back(this_reg[1]);
						one_element = std::make_pair(std::make_pair(2, o), one_p_q.second);
						p_q.push(one_element);
						local_v_p_q.push_back(one_element);


					}

				}
				else if (one_p_q.first.first == 3) {
					std::vector<int> this_p_p = one_p_q.first.second;
					if (if_merged[this_p_p[0]] == 0) {
						continue;
					}
					else {
						std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;

						std::vector<int> o;
						o = this_p_p;
						o[0] = respective_planes[o[0]];

						one_element = std::make_pair(std::make_pair(3, o), one_p_q.second);
						p_q.push(one_element);
						local_v_p_q.push_back(one_element);


					}

				}
				else if (one_p_q.first.first == 4) {
					std::vector<int> this_p_p = one_p_q.first.second;
					if (if_merged[this_p_p[0]] == 0) {
						continue;
					}
					else {
						std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;

						std::vector<int> o;
						o = this_p_p;
						o[0] = respective_planes[o[0]];

						one_element = std::make_pair(std::make_pair(4, o), one_p_q.second);
						p_q.push(one_element);
						local_v_p_q.push_back(one_element);


					}

				}
			}

			//add new operators, splite the new two primitives.
			std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element_1;
			bool if_satisfied = convert_form_split_l2(planes_to_inliers.size() - 1, max_list_vector, min_list_vector, one_element_1);
			if (if_satisfied) {
				p_q.push(one_element_1);
				local_v_p_q.push_back(one_element_1);
			}

			std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element_2;
			bool if_satisfied_2 = convert_form_split_l2(planes_to_inliers.size() - 2, max_list_vector, min_list_vector, one_element_2);
			if (if_satisfied_2) {
				p_q.push(one_element_2);
				local_v_p_q.push_back(one_element_2);
			}
			//add new operators, exlude and insert for two new primitives.
			if (all_t_transfer > 0) {
				std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element_b;
				bool if_satisfied_b = update_bad_points(planes_to_inliers.size() - 2, one_element_b);
				if (if_satisfied_b) {
					p_q.push(one_element_b);
					local_v_p_q.push_back(one_element_b);
				}
			}
			std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element_g;
			bool if_satisfied_g = update_good_points(planes_to_inliers.size() - 2, one_element_g);
			if (if_satisfied_g) {
				p_q.push(one_element_g);
				local_v_p_q.push_back(one_element_g);
			}
			if (all_t_transfer > 0) {
				std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element_b_2;
				bool if_satisfied_b_2 = update_bad_points(planes_to_inliers.size() - 1, one_element_b_2);
				if (if_satisfied_b_2) {
					p_q.push(one_element_b_2);
					local_v_p_q.push_back(one_element_b_2);
				}
			}
			std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element_g_2;
			bool if_satisfied_g_2 = update_good_points(planes_to_inliers.size() - 1, one_element_g_2);
			if (if_satisfied_g_2) {
				p_q.push(one_element_g_2);
				local_v_p_q.push_back(one_element_g_2);
			}

			v_p_q.clear();
			v_p_q = local_v_p_q;



		}
		else if (one_merged.size() == 0 && one_splite == -1 && one_remove.size() != 0 && one_add.size() == 0) {
			for (int iii : one_remove) {
				points_if_added[iii]++;
			}
			std::vector<int> chaged_in_p;
			for (int id : planes_to_inliers[one_remove_plane]) {

				if (std::find(one_remove.begin(), one_remove.end(), id) == one_remove.end()) {
					chaged_in_p.push_back(id);
				}
			}


			planes_to_inliers[one_remove_plane] = chaged_in_p;




			std::vector<Inexact_Point_3> chaged_in_p_inliers;
			chaged_in_p_inliers.reserve(chaged_in_p.size());

			for (int jj = 0; jj < chaged_in_p.size(); ++jj) {
				//points_changed[one_merge_points[jj]] = -1;
				//const Inexact_Point_3 & pt = points[max_list[jj]].first;
				chaged_in_p_inliers.push_back(points[chaged_in_p[jj]].first);

			}
			if (chaged_in_p_inliers.size() > min_points) {

				Inexact_Plane changed_plane;
				if (chaged_in_p_inliers.size() < 3) {
				
					getchar();
				}
				linear_least_squares_fitting_3(chaged_in_p_inliers.begin(), chaged_in_p_inliers.end(), changed_plane, CGAL::Dimension_tag<0>());

				planes_2[one_remove_plane] = changed_plane;
				std::vector<int> inliers_to_planes_merge_local;
				inliers_to_planes_merge_local = std::vector<int>(points.size(), -1);
				for (int m = 0; m < planes_to_inliers.size(); ++m) {
					for (int k = 0; k < planes_to_inliers[m].size(); ++k) {
						inliers_to_planes_merge_local[planes_to_inliers[m][k]] = m;
					}


				}
				inliers_to_planes.clear();
				inliers_to_planes = inliers_to_planes_merge_local;
				std::vector<std::pair<std::pair<int, std::vector<int>>, std::vector<double>>> local_v_p_q;
				while (!p_q.empty()) p_q.pop();
				for (std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_p_q : v_p_q) {
					if (one_p_q.first.first == 1) {
						std::vector<int> this_pair = one_p_q.first.second;
						//double this_weight = p_q_c.top().second;

						if (if_merged[this_pair[0]] == 0 || if_merged[this_pair[1]] == 0) {

							std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;
							bool if_satisfied = convert_form_merge_l2(this_pair[0], this_pair[1], one_element);
							if (!if_satisfied) continue;


							p_q.push(one_element);
							local_v_p_q.push_back(one_element);




						}

						else {


							p_q.push(one_p_q);
							local_v_p_q.push_back(one_p_q);





						}
					}
					else if (one_p_q.first.first == 2) {
						std::vector<int> this_reg = one_p_q.first.second;
						//double this_weight = p_q_c.top().second;
						if (if_merged[this_reg[0]] == 0)
						{
							std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element_1;
							bool if_satisfied = convert_form_split_l2(this_reg[0], max_list_vector, min_list_vector, one_element_1);
							if (!if_satisfied) continue;

							p_q.push(one_element_1);
							local_v_p_q.push_back(one_element_1);


						}

						else {

							p_q.push(one_p_q);
							local_v_p_q.push_back(one_p_q);


						}

					}
					else if (one_p_q.first.first == 3) {
						std::vector<int> this_p_p = one_p_q.first.second;
						//double this_weight = p_q_c.top().second;
						if (if_merged[this_p_p[0]] == 0)
						{
							std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element_b;
							bool if_satisfied_b = update_bad_points(this_p_p[0], one_element_b);
							if (if_satisfied_b) {
								p_q.push(one_element_b);
								local_v_p_q.push_back(one_element_b);
							}




						}


						else {

							p_q.push(one_p_q);
							local_v_p_q.push_back(one_p_q);


						}

					}
					else if (one_p_q.first.first == 4) {
						std::vector<int> this_p_p = one_p_q.first.second;


						if (if_merged[this_p_p[0]] == 0) {
							std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element_g;
							bool if_satisfied_g = update_good_points(this_p_p[0], one_element_g);
							if (if_satisfied_g) {
								p_q.push(one_element_g);
								local_v_p_q.push_back(one_element_g);
							}

						}

						else {

							p_q.push(one_p_q);
							local_v_p_q.push_back(one_p_q);


						}

					}

				}


				v_p_q.clear();
				v_p_q = local_v_p_q;

			}

		}
		else if (one_merged.size() == 0 && one_splite == -1 && one_remove.size() == 0 && one_add.size() != 0) {

			planes_to_inliers[one_add_plane].insert(planes_to_inliers[one_add_plane].end(), one_add.begin(), one_add.end());
			for (int iii : one_add) {
				points_if_added[iii]++;
			}



			std::vector<Inexact_Point_3> chaged_in_p_inliers;
			chaged_in_p_inliers.reserve(planes_to_inliers[one_add_plane].size());

			for (int jj = 0; jj < planes_to_inliers[one_add_plane].size(); ++jj) {
				//points_changed[one_merge_points[jj]] = -1;
				//const Inexact_Point_3 & pt = points[max_list[jj]].first;
				chaged_in_p_inliers.push_back(points[planes_to_inliers[one_add_plane][jj]].first);

			}
			if (chaged_in_p_inliers.size() > min_points) {
				Inexact_Plane changed_plane;
				if (chaged_in_p_inliers.size() < 3) {
				
				
					getchar();
				}
				linear_least_squares_fitting_3(chaged_in_p_inliers.begin(), chaged_in_p_inliers.end(), changed_plane, CGAL::Dimension_tag<0>());

				planes_2[one_add_plane] = changed_plane;


				std::vector<int> inliers_to_planes_merge_local;
				inliers_to_planes_merge_local = std::vector<int>(points.size(), -1);
				for (int m = 0; m < planes_to_inliers.size(); ++m) {
					for (int k = 0; k < planes_to_inliers[m].size(); ++k) {
						inliers_to_planes_merge_local[planes_to_inliers[m][k]] = m;
					}


				}
				inliers_to_planes.clear();
				inliers_to_planes = inliers_to_planes_merge_local;

				std::vector<std::pair<std::pair<int, std::vector<int>>, std::vector<double>>> local_v_p_q;
				while (!p_q.empty()) p_q.pop();
				for (std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_p_q : v_p_q) {
					if (one_p_q.first.first == 1) {
						

						std::vector<int> this_pair = one_p_q.first.second;
					

						if (if_merged[this_pair[0]] == 0 || if_merged[this_pair[1]] == 0) {

							std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;
							bool if_satisfied = convert_form_merge_l2(this_pair[0], this_pair[1], one_element);
							if (!if_satisfied) continue;
							p_q.push(one_element);
							local_v_p_q.push_back(one_element);




						}

						else {
							std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;



							p_q.push(one_p_q);
							local_v_p_q.push_back(one_p_q);

						}
					}
					else if (one_p_q.first.first == 2) {


						std::vector<int> this_reg = one_p_q.first.second;
						//double this_weight = p_q_c.top().second;
						if (if_merged[this_reg[0]] == 0)
						{
							std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element_1;
							bool if_satisfied = convert_form_split_l2(this_reg[0], max_list_vector, min_list_vector, one_element_1);
							if (!if_satisfied) continue;

							p_q.push(one_element_1);
							local_v_p_q.push_back(one_element_1);


						}

						else {

							p_q.push(one_p_q);
							local_v_p_q.push_back(one_p_q);


						}

					}
					else if (one_p_q.first.first == 3) {

						std::vector<int> this_p_p = one_p_q.first.second;

						if (if_merged[this_p_p[0]] == 0) {

							std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element_b;
							bool if_satisfied_b = update_bad_points(this_p_p[0], one_element_b);
							if (if_satisfied_b) {
								p_q.push(one_element_b);
								local_v_p_q.push_back(one_element_b);
							}
						}
						else {

							p_q.push(one_p_q);
							local_v_p_q.push_back(one_p_q);


						}

					}
					else if (one_p_q.first.first == 4) {

						std::vector<int> this_p_p = one_p_q.first.second;

						if (if_merged[this_p_p[0]] == 0) {
							std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element_g;
							bool if_satisfied_g = update_good_points(this_p_p[0], one_element_g);
							if (if_satisfied_g) {
								p_q.push(one_element_g);
								local_v_p_q.push_back(one_element_g);
							}

						}


						else {

							p_q.push(one_p_q);
							local_v_p_q.push_back(one_p_q);

						}

					}

				}


				v_p_q.clear();
				v_p_q = local_v_p_q;

			}

		}
		else {
			std::vector<int> inliers_to_planes_merge_local;
			inliers_to_planes_merge_local = std::vector<int>(points.size(), -1);
			for (int m = 0; m < planes_to_inliers.size(); ++m) {
				for (int k = 0; k < planes_to_inliers[m].size(); ++k) {
					inliers_to_planes_merge_local[planes_to_inliers[m][k]] = m;
				}


			}
			inliers_to_planes.clear();
			inliers_to_planes = inliers_to_planes_merge_local;
			bb = false;
		}
		if (t_m >= 1000) bb = false;
	
	} while (bb == true);

}

bool Shape_Detector::test_if_connected(int i, int j) {



	int Nb_neigh = 10;
	for (int f = 0; f < planes_to_inliers[i].size(); ++f) {

		if ((int)spherical_neighborhood[planes_to_inliers[i][f]].size() < Nb_neigh) {
			Nb_neigh = (int)spherical_neighborhood[planes_to_inliers[i][f]].size();
		}

		for (int it = 0; it < Nb_neigh; it++) {

			int neighbor_index = spherical_neighborhood[planes_to_inliers[i][f]][it];

			for (int id : planes_to_inliers[j]) {
				if (id == neighbor_index) {

					return true;
				}

			}

		}

	}
	return false;
}



bool Shape_Detector::convert_form_merge_normal(int i, int j, std::pair<std::pair<int, std::vector<int>>, std::vector<double>>& one_element) {

	std::vector<int> o;
	o.push_back(i);
	o.push_back(j);
	double dif_m = 0;
	std::vector<double> s_d;
	std::vector<int> merge_moved_ids;
	//we suppose that the pair of primitives can not merged,when their normal cosine smaller than normal_threshold.
	if (abs(planes_2[i].orthogonal_vector()*planes_2[j].orthogonal_vector()) < normal_threshold) {
		return false;
	}
	double all_changed_distances = merge_distance_changed_normal_with_epsilon(i, j, dif_m, merge_moved_ids);

	/*if (merge_moved_ids.size() > points.size()*0.01|| merge_moved_ids.size() > planes_to_inliers[i].size()/4|| merge_moved_ids.size() > planes_to_inliers[j].size()/4) {
		return false;
	}*/
	if ((planes_to_inliers[i].size() + planes_to_inliers[j].size() - merge_moved_ids.size()) < min_points) {
		return false;
	}
	double energy_ccc = energy_changed_normal_merge(all_changed_distances, -1, merge_moved_ids.size());
	if (energy_ccc > 0) {
		return false;
	}
	s_d.push_back(energy_ccc);
	o.insert(o.end(), merge_moved_ids.begin(), merge_moved_ids.end());

	s_d.push_back(dif_m);
	one_element = std::make_pair(std::make_pair(1, o), s_d);
	return true;
}

bool Shape_Detector::convert_form_split_normal(int i, std::vector<std::vector<int>>& max_list_vector, std::vector<std::vector<int>>& min_list_vector, std::pair<std::pair<int, std::vector<int>>, std::vector<double>> & one_element) {
	std::vector<int> max_list;
	std::vector<int> min_list;
	double dif_s = 0;
	bool if_split = separate_two_out_normal(i, max_list, min_list, dif_s);

	double energy_ccc = energy_changed_normal(dif_s, 1.0);

	if (energy_ccc >= 0 || max_list.size() < min_points || min_list.size() < min_points || !if_split) {
		return false;
	}
	max_list_vector.push_back(max_list);
	min_list_vector.push_back(min_list);
	std::vector<int> o;
	o.push_back(i);
	o.push_back(max_list_vector.size() - 1);
	std::vector<double> s_d;
	s_d.push_back(energy_ccc);

	s_d.push_back(dif_s);
	one_element = std::make_pair(std::make_pair(2, o), s_d);
	return true;
}

bool Shape_Detector::convert_form_exclude_normal(int i, std::pair<std::pair<int, std::vector<int>>, std::vector<double>>& one_element) {

	std::vector<int> o;
	o.push_back(bad_points_shape[i].first.first);

	o.insert(o.end(), bad_points_shape[i].first.second.begin(), bad_points_shape[i].first.second.end());

	one_element = std::make_pair(std::make_pair(3, o), bad_points_shape[i].second);
	return true;
}

bool Shape_Detector::convert_form_insert_normal(int i, std::pair<std::pair<int, std::vector<int>>, std::vector<double>>& one_element) {
	std::vector<int> o;
	o.push_back(good_points_shape[i].first.first);
	o.insert(o.end(), good_points_shape[i].first.second.begin(), good_points_shape[i].first.second.end());


	one_element = std::make_pair(std::make_pair(4, o), good_points_shape[i].second);
	return true;
}


bool Shape_Detector::update_good_points_normal(int id_shape, std::pair<std::pair<int, std::vector<int>>, std::vector<double>>& one_element) {

	int Nb_neigh = 10;
	std::vector<std::pair<int, double>> ids_good_outliers;
	for (int i = 0; i < points.size(); ++i) {
		if (if_insert_candidate[i] == 0) continue;
		if (inliers_to_planes[i] != -1) continue;
		if (points_if_added[i] > 2) continue;//avoide infinite loop

		if ((int)spherical_neighborhood[i].size() < Nb_neigh) {
			Nb_neigh = (int)spherical_neighborhood[i].size();
		}

		std::set<int> one_neight_id;

		bool if_correspond = false;
		for (int it = 0; it < Nb_neigh; ++it) {

			int neighbor_index = spherical_neighborhood[i][it];

			if (inliers_to_planes[neighbor_index] != -1) {

				one_neight_id.insert(inliers_to_planes[neighbor_index]);
				if (inliers_to_planes[neighbor_index] == id_shape) {

					if_correspond = true;
				}
			}


		}

		if (one_neight_id.empty() || !if_correspond) { continue; }



		double mm_normal = normal_threshold;

		int changed_plane_id = -1;
		Point_with_normal this_p = points[i];

		for (int neight_id : one_neight_id) {
			if (abs(this_p.second * planes_2[neight_id].orthogonal_vector()) > mm_normal) {
				if (sqrt(CGAL::squared_distance(planes_2[neight_id], this_p.first)) < epsilon) {

					mm_normal = abs(this_p.second * planes_2[neight_id].orthogonal_vector());
					changed_plane_id = neight_id;

				}

			}
		}


		if (changed_plane_id == id_shape) {

			ids_good_outliers.push_back(std::make_pair(i, mm_normal));//plane_id, point_id, normal_deviation.

		}

	}


	std::priority_queue<std::pair<int, double>, std::vector<std::pair<int, double>>, Add_Comparator_normal> p_q;

	for (std::pair<int, double> ppp : ids_good_outliers) {
		p_q.push(ppp);
	}

	std::vector<int> added_points_ids;
	double en = 0;
	double dif_e = 0;
	double n_a = 0;

	while (!p_q.empty() && n_a < number_of_insert_exclude) {
		for (int idd : added_points_ids) {
			if_insert_candidate[idd] = 0;
		}
		added_points_ids.push_back(p_q.top().first);
		p_q.pop();
		n_a++;
	}
	dif_e = add_changed_error_normal(id_shape, added_points_ids);
	en = energy_changed_second_normal(dif_e, -n_a);

	if (added_points_ids.size() != 0 && en < 0) {
		for (int idd : added_points_ids) {
			if_insert_candidate[idd] = 0;
		}
		std::vector<double> dd;
		dd.push_back(en);
		dd.push_back(dif_e);

		std::vector<int> o;
		o.push_back(id_shape);
		o.insert(o.end(), added_points_ids.begin(), added_points_ids.end());


		one_element = std::make_pair(std::make_pair(4, o), dd);


		return true;

	}
	else {
		return false;
	}


}


bool Shape_Detector::update_bad_points_normal(int id_shape, std::pair<std::pair<int, std::vector<int>>, std::vector<double>>& one_element) {



	std::vector<std::pair<int, double>> ids_bad_inliers;
	double bigest_normal = mean_normal_diviation;
	//double bigest_normal = 0.01;



	const Inexact_Plane & H = planes_2[id_shape];


	for (int j : planes_to_inliers[id_shape]) {
		const Point_with_normal & pt = points[j];
		if (points_if_added[j] > 2) continue;

		if (abs(H.orthogonal_vector()*pt.second) < bigest_normal) {
			ids_bad_inliers.push_back(std::make_pair(j, abs(H.orthogonal_vector()*pt.second)));
		}

	}





	std::priority_queue<std::pair<int, double>, std::vector<std::pair<int, double>>, Remove_Comparator_normal> p_q;
	for (std::pair<int, double> ppp : ids_bad_inliers) {
		p_q.push(ppp);

	}
	std::vector<int> removed_points_ids;
	double en = 0;
	double dif_e = 0;

	double n_b = 0;
	while (!p_q.empty() && n_b < number_of_insert_exclude) {
		removed_points_ids.push_back(p_q.top().first);
		p_q.pop();
		n_b++;
	}
	dif_e = remove_changed_error_normal(id_shape, removed_points_ids);
	en = energy_changed_second_normal(dif_e, n_b);

	if (removed_points_ids.size() != 0 && (planes_to_inliers[id_shape].size() - removed_points_ids.size()) >= min_points && en < 0) {
		std::vector<double> dd;
		dd.push_back(en);
		dd.push_back(dif_e);
		std::vector<int> o;
		o.push_back(id_shape);
		o.insert(o.end(), removed_points_ids.begin(), removed_points_ids.end());


		one_element = std::make_pair(std::make_pair(3, o), dd);
		return true;
	}
	else { return false; }



}

void Shape_Detector::local_operators_normal() {

	bool bb = false;
	all_t_merge += t_merge;
	all_t_split += t_split;
	all_t_insert += t_insert;
	all_t_exlude += t_exlude;
	t_m = 0;
	t_merge = 0;
	t_split = 0;
	t_exlude = 0;
	t_insert = 0;



	std::priority_queue<std::pair<std::pair<int, std::vector<int>>, std::vector<double>>, std::vector<std::pair<std::pair<int, std::vector<int>>, std::vector<double>>>, Weight_Comparator_with_energy> p_q;
	std::vector<std::pair<std::pair<int, std::vector<int>>, std::vector<double>>> v_p_q;
	std::vector<std::vector<int>> max_list_vector;
	std::vector<std::vector<int>> min_list_vector;

	//add the merging operation to priority queue
	for (int i = 0; i < primitive_connection.size(); ++i) {


		std::vector<int> one = primitive_connection[i];
		for (int j = i + 1; j < one.size(); ++j) {
			if (one[j] >= 1) {

				std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;
				bool if_satisfied = convert_form_merge_normal(i, j, one_element);
				if (!if_satisfied)continue;

				p_q.push(one_element);
				v_p_q.push_back(one_element);

			}
		}
	}

	//add the spliting operation to priority queue
	for (int i = 0; i < planes_to_inliers.size(); ++i) {
		if (region_type[i] > 5) continue;
		std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;
		bool if_satisfied = convert_form_split_normal(i, max_list_vector, min_list_vector, one_element);
		if (!if_satisfied) continue;
		p_q.push(one_element);
		v_p_q.push_back(one_element);


	}
	//add the exclusion operation to priority queue
	for (int i = 0; i < bad_points_shape.size(); ++i) {

		std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;

		bool if_satisfied = convert_form_exclude_normal(i, one_element);
		if (!if_satisfied) continue;
		p_q.push(one_element);
		v_p_q.push_back(one_element);

	}
	//add the insertion operation to priority queue
	for (int i = 0; i < good_points_shape.size(); ++i) {




		std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;
		bool if_satisfied = convert_form_insert_normal(i, one_element);
		if (!if_satisfied) continue;
		p_q.push(one_element);
		v_p_q.push_back(one_element);

	}



	do {


		bb = true;

		std::vector<int> if_merged;
		if_merged = std::vector<int>(planes_to_inliers.size(), -1);



		std::vector<int> one_merged;//ids of merged primitives.
		std::vector<int> id_of_moved_merge;//ids of the ecluded points because of merging.
		std::vector<int> max_list;
		std::vector<int> min_list;


		int one_splite = -1;//id of split primitive.
		std::vector<int> one_remove;//ids of excluded inliers.
		int one_remove_plane = -1;//id of the primitive, from which the inliers are excluded.
		std::vector<int> one_add;//ids of inserted outliers.
		int one_add_plane = -1;//id of the primitive, to which the outliers are inserted.


		//test if the top operation is available, if it's valable, the operation is used to update the configuration and the priority queue is
		//also updated, otherwise, the top operation is poped and the next operation is tested until the top opeartion can not decrease the energy.

		while (!p_q.empty()) {


			if (p_q.top().second[0] > 0) break;//when the top operation can not decrease the energy, stop.
			if (p_q.top().first.first == 1) {//merge operation in the top of priority queue.

				std::vector<int> this_pair = p_q.top().first.second;//ids of the pair of merged primitives + ids of excluded points.

				int f = this_pair[0];
				int s = this_pair[1];
				if ((planes_to_inliers[f].size() + planes_to_inliers[s].size() + 2 - this_pair.size()) < min_points) {//when the new primitives have less point than threshold, the operation will be ignoed.
					p_q.pop();
					continue;
				}



				double after_operation_distance_visual = all_normal_diaviation + p_q.top().second[1];
				if (if_constraint) {

					if (number_of_assigned_points + 2 - this_pair.size() < ori_inliers_number) {//when the current inlier number is smaller than the original one.

						p_q.pop();
						continue;
					}

					if (after_operation_distance_visual / double(number_of_assigned_points) < ori_mean_normal_diviation) {//when the current normal deviation is bigger than orig

						p_q.pop();
						continue;


					}

				}

				p_q.pop();
				id_of_moved_merge.clear();
				//saving the points that far away from the new plane.
				if (this_pair.size() > 2) {
					for (int k = 2; k < this_pair.size(); ++k) {
						id_of_moved_merge.push_back(this_pair[k]);
					}
				}
				number_of_assigned_points = number_of_assigned_points - id_of_moved_merge.size();

				if_merged[f] = 0;
				if_merged[s] = 0;
				one_merged.push_back(f);
				one_merged.push_back(s);

				all_normal_diaviation = after_operation_distance_visual;

				t_m++;
				t_merge++;

				
				break;


			}
			else if (p_q.top().first.first == 2) {//split operator in the top of priority queue.


				int splited_index = p_q.top().first.second[0];
				int info_index = p_q.top().first.second[1];

				double after_splited_distance = all_normal_diaviation + p_q.top().second[1];



				max_list = max_list_vector[info_index];
				min_list = min_list_vector[info_index];

				if (if_constraint) {
					//when current normal deviation is bigger than original one or current number of primitives is bigger than original one.
					if (after_splited_distance / double(number_of_assigned_points) < ori_mean_normal_diviation || (planes_2.size() + 1) > ori_primitives_number) {

						p_q.pop();
						continue;
					}


				}

				p_q.pop();
				if_merged[splited_index] = 0;

				one_splite = splited_index;

				all_normal_diaviation = after_splited_distance;

				t_m++;
				t_split++;
			

				break;

			}
			else if (p_q.top().first.first == 3) {//exclusion operator in the top of priority queue.

				int plane_id = p_q.top().first.second[0];

				double after_delete_distance = all_normal_diaviation + p_q.top().second[1];

				if (if_constraint) {
					//when the current number of inliers is smaller than original one or the current normal devidation is bigger than original one.
					if (number_of_assigned_points + 1 - p_q.top().first.second.size() < ori_inliers_number || after_delete_distance / double(number_of_assigned_points + 1 - p_q.top().first.second.size()) < ori_mean_normal_diviation) {

						p_q.pop();
						continue;
					}
				}
				//when the inliers number of the primitive is smaller than threhsold.
				if (planes_to_inliers[plane_id].size() + 1 - p_q.top().first.second.size() < min_points) {

					p_q.pop();
					continue;
				}


				if_merged[plane_id] = 0;
				one_remove_plane = plane_id;
				std::vector<int> points_ids = p_q.top().first.second;
				points_ids.erase(points_ids.begin());
				one_remove = points_ids;

				all_normal_diaviation = after_delete_distance;

				number_of_assigned_points = number_of_assigned_points - one_remove.size();

				p_q.pop();
				t_m++;
				t_exlude++;
			

				break;

			}
			else if (p_q.top().first.first == 4) {//insert operator in the top of priority queue.

				int plane_id = p_q.top().first.second[0];


				double after_add_distance = all_normal_diaviation + p_q.top().second[1];
				if (if_constraint) {
					//when the current normal deviation is bigger than original one.
					if (after_add_distance / double(number_of_assigned_points - 1 + p_q.top().first.second.size()) < ori_mean_normal_diviation) {
						p_q.pop();
						continue;
					}
				}
				if_merged[plane_id] = 0;
				one_add_plane = plane_id;
				std::vector<int> points_ids = p_q.top().first.second;
				points_ids.erase(points_ids.begin());
				one_add = points_ids;

				all_normal_diaviation = after_add_distance;

				number_of_assigned_points = number_of_assigned_points + one_add.size();

				p_q.pop();
				t_m++;
				t_insert++;

			

				break;


			}
		}



		//implement the operation and update the priority queue.
		if (one_merged.size() != 0 && one_splite == -1 && one_remove.size() == 0 && one_add.size() == 0) {//merge operator.
			int last_size = planes_2.size();
			std::vector<int> respective_planes = std::vector<int>(last_size, -1);


			int id_1 = one_merged[0];
			int id_2 = one_merged[1];
			std::vector<int> one_merge_points_real;
			std::vector<int> one_merge_points = planes_to_inliers[id_1];
			one_merge_points.insert(one_merge_points.end(), planes_to_inliers[id_2].begin(), planes_to_inliers[id_2].end());
			for (int id : one_merge_points) {
				if (std::find(id_of_moved_merge.begin(), id_of_moved_merge.end(), id) == id_of_moved_merge.end()) {
					one_merge_points_real.push_back(id);
				}
			}

			one_merge_points = one_merge_points_real;
			for (int m_id : id_of_moved_merge) {//exclude the far away inliers.
				inliers_to_planes[m_id] = -1;

			}
			if (one_merge_points.size() >= min_points) {
				if (id_2 < id_1) {
					for (int t = 0; t < last_size; ++t) {
						if (t < id_2) {
							respective_planes[t] = t;
						}
						else if (t == id_2) {
							respective_planes[t] = last_size - 2;
						}
						else if (t < id_1) {
							respective_planes[t] = t - 1;

						}
						else if (t == id_1) {
							respective_planes[t] = last_size - 2;

						}
						else {
							respective_planes[t] = t - 2;

						}


					}
					region_type.push_back(region_type[id_1] + region_type[id_2] + 1);
					region_type.erase(region_type.begin() + id_1);
					region_type.erase(region_type.begin() + id_2);



					planes_to_inliers.erase(planes_to_inliers.begin() + id_1);
					planes_to_inliers.erase(planes_to_inliers.begin() + id_2);
					planes_2.erase(planes_2.begin() + id_1);
					planes_2.erase(planes_2.begin() + id_2);


				}
				else {
					for (int t = 0; t < last_size; ++t) {
						if (t < id_1) {
							respective_planes[t] = t;
						}
						else if (t == id_1) {
							respective_planes[t] = last_size - 2;
						}
						else if (t < id_2) {
							respective_planes[t] = t - 1;

						}
						else if (t == id_2) {
							respective_planes[t] = last_size - 2;

						}
						else {
							respective_planes[t] = t - 2;

						}


					}
					region_type.push_back(region_type[id_1] + region_type[id_2] + 1);
					region_type.erase(region_type.begin() + id_2);
					region_type.erase(region_type.begin() + id_1);

					planes_to_inliers.erase(planes_to_inliers.begin() + id_2);
					planes_to_inliers.erase(planes_to_inliers.begin() + id_1);
					planes_2.erase(planes_2.begin() + id_2);
					planes_2.erase(planes_2.begin() + id_1);


				}


				planes_to_inliers.push_back(one_merge_points);




				std::vector<Inexact_Point_3> inliers;
				inliers.reserve(one_merge_points.size());

				for (int jj = 0; jj < one_merge_points.size(); ++jj) {

					inliers.push_back(points[one_merge_points[jj]].first);

				}
				Inexact_Plane plane;
				if (inliers.size() < 3) {

					continue;
				}
				linear_least_squares_fitting_3(inliers.begin(), inliers.end(), plane, CGAL::Dimension_tag<0>());

				planes_2.push_back(plane);

				//update inliers_to_planes.
				std::vector<int> inliers_to_planes_merge_local;
				inliers_to_planes_merge_local = std::vector<int>(points.size(), -1);
				for (int m = 0; m < planes_to_inliers.size(); ++m) {
					for (int k = 0; k < planes_to_inliers[m].size(); ++k) {
						inliers_to_planes_merge_local[planes_to_inliers[m][k]] = m;
					}


				}
				inliers_to_planes.clear();
				inliers_to_planes = inliers_to_planes_merge_local;

				//update priority queue
				std::vector<std::pair<std::pair<int, std::vector<int>>, std::vector<double>>> local_v_p_q;
				while (!p_q.empty()) p_q.pop();//clean p_q.
				//update the existing operations
				for (std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_p_q : v_p_q) {
					if (one_p_q.first.first == 1) {//when the operator is merge.
						std::vector<int> this_pair = one_p_q.first.second;
						if (if_merged[this_pair[0]] == 0 && if_merged[this_pair[1]] == 0)//when that is the merged pair.
						{


							continue;
						}
						else if ((if_merged[this_pair[0]] == 0 && if_merged[this_pair[1]] != 0) || (if_merged[this_pair[0]] != 0 && if_merged[this_pair[1]] == 0)) {//when only one primitive is merged.


							std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;
							bool if_satisfied = convert_form_merge_normal(respective_planes[this_pair[1]], respective_planes[this_pair[0]], one_element);
							if (!if_satisfied) continue;
							p_q.push(one_element);
							local_v_p_q.push_back(one_element);




						}
						else {//when the operation is not affected.



							std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;

							std::vector<int> o = this_pair;
							o[0] = respective_planes[this_pair[0]];
							o[1] = respective_planes[this_pair[1]];

							one_element = std::make_pair(std::make_pair(1, o), one_p_q.second);
							p_q.push(one_element);
							local_v_p_q.push_back(one_element);


						}
					}
					else if (one_p_q.first.first == 2) {//when the operator is split.
						std::vector<int> this_reg = one_p_q.first.second;
						if (if_merged[this_reg[0]] == 0)//when the primitive has been merged.
						{
							continue;
						}

						else {//when the operation is not affected.
							std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;

							std::vector<int> o;
							o.push_back(respective_planes[this_reg[0]]);
							o.push_back(this_reg[1]);
							one_element = std::make_pair(std::make_pair(2, o), one_p_q.second);
							p_q.push(one_element);
							local_v_p_q.push_back(one_element);

						}

					}
					else if (one_p_q.first.first == 3) {//when the operator is exclude.
						std::vector<int> this_p_p = one_p_q.first.second;
						if (if_merged[this_p_p[0]] == 0) {
							continue;
						}
						else {
							std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;

							std::vector<int> o;
							o = this_p_p;
							o[0] = respective_planes[o[0]];

							one_element = std::make_pair(std::make_pair(3, o), one_p_q.second);
							p_q.push(one_element);
							local_v_p_q.push_back(one_element);


						}

					}
					else if (one_p_q.first.first == 4) {//when the operator is insert.
						std::vector<int> this_p_p = one_p_q.first.second;
						if (if_merged[this_p_p[0]] == 0) {
							continue;
						}
						else {
							std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;
							std::vector<int> o;
							o = this_p_p;
							o[0] = respective_planes[o[0]];
							one_element = std::make_pair(std::make_pair(4, o), one_p_q.second);
							p_q.push(one_element);
							local_v_p_q.push_back(one_element);


						}

					}
				}


				//add new operator, split the new merged primitive
				std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element_1;
				bool if_satisfied = convert_form_split_normal(last_size - 2, max_list_vector, min_list_vector, one_element_1);
				if (if_satisfied) {
					p_q.push(one_element_1);
					local_v_p_q.push_back(one_element_1);
				}
				//add new operat, exlude and inerte for the new merged primitive
				if (all_t_transfer > 0) {
					std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element_b;
					bool if_satisfied_b = update_bad_points_normal(last_size - 2, one_element_b);
					if (if_satisfied_b) {
						p_q.push(one_element_b);
						local_v_p_q.push_back(one_element_b);
					}
				}
				std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element_g;
				bool if_satisfied_g = update_good_points_normal(last_size - 2, one_element_g);
				if (if_satisfied_g) {
					p_q.push(one_element_g);
					local_v_p_q.push_back(one_element_g);
				}

				v_p_q.clear();
				v_p_q = local_v_p_q;

			}

		}
		else if (one_merged.size() == 0 && one_splite != -1 && one_remove.size() == 0 && one_add.size() == 0) {//split operator.
			int last_size = planes_2.size();


			std::vector<int> respective_planes = std::vector<int>(last_size, -1);


			for (int t = 0; t < last_size; ++t) {
				if (t < one_splite) {
					respective_planes[t] = t;
				}
				else if (t == one_splite) {
					respective_planes[t] = last_size;
				}

				else {
					respective_planes[t] = t - 1;

				}


			}
			planes_to_inliers.erase(planes_to_inliers.begin() + one_splite);
			region_type.push_back(region_type[one_splite]);
			region_type.push_back(region_type[one_splite]);
			region_type.erase(region_type.begin() + one_splite);



			planes_2.erase(planes_2.begin() + one_splite);


			planes_to_inliers.push_back(max_list);
			planes_to_inliers.push_back(min_list);


			std::vector<Inexact_Point_3> max_inliers;
			max_inliers.reserve(max_list.size());

			for (int jj = 0; jj < max_list.size(); ++jj) {

				max_inliers.push_back(points[max_list[jj]].first);

			}
			Inexact_Plane max_plane;
			if (max_inliers.size() < 3) {
				getchar();
			}
			linear_least_squares_fitting_3(max_inliers.begin(), max_inliers.end(), max_plane, CGAL::Dimension_tag<0>());

			planes_2.push_back(max_plane);


			std::vector<Inexact_Point_3> min_inliers;
			min_inliers.reserve(min_list.size());

			for (int jj = 0; jj < min_list.size(); ++jj) {

				min_inliers.push_back(points[min_list[jj]].first);

			}
			Inexact_Plane min_plane;
			if (min_inliers.size() < 5) {
				getchar();
			}
			linear_least_squares_fitting_3(min_inliers.begin(), min_inliers.end(), min_plane, CGAL::Dimension_tag<0>());

			planes_2.push_back(min_plane);
			//update inliers to plane
			std::vector<int> inliers_to_planes_merge_local;
			inliers_to_planes_merge_local = std::vector<int>(points.size(), -1);
			for (int m = 0; m < planes_to_inliers.size(); ++m) {
				for (int k = 0; k < planes_to_inliers[m].size(); ++k) {
					inliers_to_planes_merge_local[planes_to_inliers[m][k]] = m;
				}


			}
			inliers_to_planes.clear();
			inliers_to_planes = inliers_to_planes_merge_local;
			//update priority queue
			std::vector<std::pair<std::pair<int, std::vector<int>>, std::vector<double>>> local_v_p_q;
			while (!p_q.empty()) p_q.pop();
			for (std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_p_q : v_p_q) {
				if (one_p_q.first.first == 1) {
					std::vector<int> this_pair = one_p_q.first.second;

					if (if_merged[this_pair[0]] == 0 && if_merged[this_pair[1]] == 0)//not possible case
					{
						continue;
					}
					else if (if_merged[this_pair[0]] == 0 && if_merged[this_pair[1]] != 0) {//one of the primitive is splited
						if (test_if_connected(respective_planes[this_pair[1]], planes_to_inliers.size() - 1)) {
							std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;
							bool if_satisfied = convert_form_merge_normal(planes_to_inliers.size() - 1, respective_planes[this_pair[1]], one_element);
							if (!if_satisfied) continue;
							p_q.push(one_element);
							local_v_p_q.push_back(one_element);

						}
						if (test_if_connected(respective_planes[this_pair[1]], planes_to_inliers.size() - 2)) {
							std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;

							bool if_satisfied = convert_form_merge_normal(planes_to_inliers.size() - 2, respective_planes[this_pair[1]], one_element);
							if (!if_satisfied) continue;
							p_q.push(one_element);
							local_v_p_q.push_back(one_element);

						}



					}
					else if (if_merged[this_pair[0]] != 0 && if_merged[this_pair[1]] == 0) {//one of the primitive is splited
						if (test_if_connected(respective_planes[this_pair[0]], planes_to_inliers.size() - 1)) {
							std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;
							bool if_satisfied = convert_form_merge_normal(planes_to_inliers.size() - 1, respective_planes[this_pair[0]], one_element);
							if (!if_satisfied) continue;

							p_q.push(one_element);
							local_v_p_q.push_back(one_element);


						}
						if (test_if_connected(respective_planes[this_pair[0]], planes_to_inliers.size() - 2)) {
							std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;
							bool if_satisfied = convert_form_merge_normal(planes_to_inliers.size() - 2, respective_planes[this_pair[0]], one_element);
							if (!if_satisfied) continue;

							p_q.push(one_element);
							local_v_p_q.push_back(one_element);

						}


					}
					else {// is not affected
						std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;

						std::vector<int> o = this_pair;
						o[0] = respective_planes[this_pair[0]];
						o[1] = respective_planes[this_pair[1]];


						one_element = std::make_pair(std::make_pair(1, o), one_p_q.second);
						p_q.push(one_element);
						local_v_p_q.push_back(one_element);




					}
				}
				else if (one_p_q.first.first == 2) {
					std::vector<int> this_reg = one_p_q.first.second;

					if (if_merged[this_reg[0]] == 0)
					{
						continue;

					}

					else {
						std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;

						std::vector<int> o;
						o.push_back(respective_planes[this_reg[0]]);
						o.push_back(this_reg[1]);
						one_element = std::make_pair(std::make_pair(2, o), one_p_q.second);
						p_q.push(one_element);
						local_v_p_q.push_back(one_element);


					}

				}
				else if (one_p_q.first.first == 3) {
					std::vector<int> this_p_p = one_p_q.first.second;
					if (if_merged[this_p_p[0]] == 0) {
						continue;
					}
					else {
						std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;

						std::vector<int> o;
						o = this_p_p;
						o[0] = respective_planes[o[0]];

						one_element = std::make_pair(std::make_pair(3, o), one_p_q.second);
						p_q.push(one_element);
						local_v_p_q.push_back(one_element);


					}

				}
				else if (one_p_q.first.first == 4) {
					std::vector<int> this_p_p = one_p_q.first.second;
					if (if_merged[this_p_p[0]] == 0) {
						continue;
					}
					else {
						std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;

						std::vector<int> o;
						o = this_p_p;
						o[0] = respective_planes[o[0]];

						one_element = std::make_pair(std::make_pair(4, o), one_p_q.second);
						p_q.push(one_element);
						local_v_p_q.push_back(one_element);


					}

				}
			}

			//add new operators, splite the new two primitives.
			std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element_1;
			bool if_satisfied = convert_form_split_normal(planes_to_inliers.size() - 1, max_list_vector, min_list_vector, one_element_1);
			if (if_satisfied) {
				p_q.push(one_element_1);
				local_v_p_q.push_back(one_element_1);
			}

			std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element_2;
			bool if_satisfied_2 = convert_form_split_normal(planes_to_inliers.size() - 2, max_list_vector, min_list_vector, one_element_2);
			if (if_satisfied_2) {
				p_q.push(one_element_2);
				local_v_p_q.push_back(one_element_2);
			}
			//add new operators, exlude and insert for two new primitives.
			if (all_t_transfer > 0) {
				std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element_b;
				bool if_satisfied_b = update_bad_points_normal(planes_to_inliers.size() - 2, one_element_b);
				if (if_satisfied_b) {
					p_q.push(one_element_b);
					local_v_p_q.push_back(one_element_b);
				}
			}

			std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element_g;
			bool if_satisfied_g = update_good_points_normal(planes_to_inliers.size() - 2, one_element_g);
			if (if_satisfied_g) {
				p_q.push(one_element_g);
				local_v_p_q.push_back(one_element_g);
			}
			if (all_t_transfer > 0) {
				std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element_b_2;
				bool if_satisfied_b_2 = update_bad_points_normal(planes_to_inliers.size() - 1, one_element_b_2);
				if (if_satisfied_b_2) {
					p_q.push(one_element_b_2);
					local_v_p_q.push_back(one_element_b_2);
				}
			}

			std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element_g_2;
			bool if_satisfied_g_2 = update_good_points_normal(planes_to_inliers.size() - 1, one_element_g_2);
			if (if_satisfied_g_2) {
				p_q.push(one_element_g_2);
				local_v_p_q.push_back(one_element_g_2);
			}

			v_p_q.clear();
			v_p_q = local_v_p_q;



		}
		else if (one_merged.size() == 0 && one_splite == -1 && one_remove.size() != 0 && one_add.size() == 0) {
			for (int iii : one_remove) {
				points_if_added[iii]++;
			}
			std::vector<int> chaged_in_p;
			for (int id : planes_to_inliers[one_remove_plane]) {

				if (std::find(one_remove.begin(), one_remove.end(), id) == one_remove.end()) {
					chaged_in_p.push_back(id);
				}
			}


			planes_to_inliers[one_remove_plane] = chaged_in_p;




			std::vector<Inexact_Point_3> chaged_in_p_inliers;
			chaged_in_p_inliers.reserve(chaged_in_p.size());

			for (int jj = 0; jj < chaged_in_p.size(); ++jj) {

				chaged_in_p_inliers.push_back(points[chaged_in_p[jj]].first);

			}
			if (chaged_in_p_inliers.size() > min_points) {
				Inexact_Plane changed_plane;
				if (chaged_in_p_inliers.size() < 3) {

					getchar();
				}
				linear_least_squares_fitting_3(chaged_in_p_inliers.begin(), chaged_in_p_inliers.end(), changed_plane, CGAL::Dimension_tag<0>());

				planes_2[one_remove_plane] = changed_plane;
				//update inliers to plane
				std::vector<int> inliers_to_planes_merge_local;
				inliers_to_planes_merge_local = std::vector<int>(points.size(), -1);
				for (int m = 0; m < planes_to_inliers.size(); ++m) {
					for (int k = 0; k < planes_to_inliers[m].size(); ++k) {
						inliers_to_planes_merge_local[planes_to_inliers[m][k]] = m;
					}


				}
				inliers_to_planes.clear();
				inliers_to_planes = inliers_to_planes_merge_local;

				//update priority queue
				std::vector<std::pair<std::pair<int, std::vector<int>>, std::vector<double>>> local_v_p_q;
				while (!p_q.empty()) p_q.pop();
				for (std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_p_q : v_p_q) {
					if (one_p_q.first.first == 1) {
						std::vector<int> this_pair = one_p_q.first.second;


						if (if_merged[this_pair[0]] == 0 || if_merged[this_pair[1]] == 0) {

							std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;
							bool if_satisfied = convert_form_merge_normal(this_pair[0], this_pair[1], one_element);
							if (!if_satisfied) continue;

							p_q.push(one_element);
							local_v_p_q.push_back(one_element);


						}

						else {

							p_q.push(one_p_q);
							local_v_p_q.push_back(one_p_q);


						}
					}
					else if (one_p_q.first.first == 2) {
						std::vector<int> this_reg = one_p_q.first.second;

						if (if_merged[this_reg[0]] == 0)
						{
							std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element_1;
							bool if_satisfied = convert_form_split_normal(this_reg[0], max_list_vector, min_list_vector, one_element_1);
							if (!if_satisfied) continue;

							p_q.push(one_element_1);
							local_v_p_q.push_back(one_element_1);


						}

						else {

							p_q.push(one_p_q);
							local_v_p_q.push_back(one_p_q);


						}

					}
					else if (one_p_q.first.first == 3) {
						std::vector<int> this_p_p = one_p_q.first.second;

						if (if_merged[this_p_p[0]] == 0)
						{
							std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element_b;
							bool if_satisfied_b = update_bad_points_normal(this_p_p[0], one_element_b);
							if (if_satisfied_b) {

								p_q.push(one_element_b);
								local_v_p_q.push_back(one_element_b);
							}


						}

						else {

							p_q.push(one_p_q);
							local_v_p_q.push_back(one_p_q);


						}

					}
					else if (one_p_q.first.first == 4) {
						std::vector<int> this_p_p = one_p_q.first.second;


						if (if_merged[this_p_p[0]] == 0) {
							std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element_g;
							bool if_satisfied_g = update_good_points_normal(this_p_p[0], one_element_g);
							if (if_satisfied_g) {
								p_q.push(one_element_g);
								local_v_p_q.push_back(one_element_g);
							}

						}

						else {

							p_q.push(one_p_q);
							local_v_p_q.push_back(one_p_q);


						}

					}

				}


				v_p_q.clear();
				v_p_q = local_v_p_q;


			}
		}
		else if (one_merged.size() == 0 && one_splite == -1 && one_remove.size() == 0 && one_add.size() != 0) {

			planes_to_inliers[one_add_plane].insert(planes_to_inliers[one_add_plane].end(), one_add.begin(), one_add.end());
			for (int iii : one_add) {
				points_if_added[iii]++;
			}



			std::vector<Inexact_Point_3> chaged_in_p_inliers;
			chaged_in_p_inliers.reserve(planes_to_inliers[one_add_plane].size());

			for (int jj = 0; jj < planes_to_inliers[one_add_plane].size(); ++jj) {

				chaged_in_p_inliers.push_back(points[planes_to_inliers[one_add_plane][jj]].first);

			}
			if (chaged_in_p_inliers.size() > min_points) {
				Inexact_Plane changed_plane;

				linear_least_squares_fitting_3(chaged_in_p_inliers.begin(), chaged_in_p_inliers.end(), changed_plane, CGAL::Dimension_tag<0>());

				planes_2[one_add_plane] = changed_plane;




				//update inlier to plane
				std::vector<int> inliers_to_planes_merge_local;
				inliers_to_planes_merge_local = std::vector<int>(points.size(), -1);
				for (int m = 0; m < planes_to_inliers.size(); ++m) {
					for (int k = 0; k < planes_to_inliers[m].size(); ++k) {
						inliers_to_planes_merge_local[planes_to_inliers[m][k]] = m;
					}


				}
				inliers_to_planes.clear();
				inliers_to_planes = inliers_to_planes_merge_local;
				//update priority queue
				std::vector<std::pair<std::pair<int, std::vector<int>>, std::vector<double>>> local_v_p_q;
				while (!p_q.empty()) p_q.pop();
				for (std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_p_q : v_p_q) {
					if (one_p_q.first.first == 1) {


						std::vector<int> this_pair = one_p_q.first.second;


						if (if_merged[this_pair[0]] == 0 || if_merged[this_pair[1]] == 0) {

							std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element;
							bool if_satisfied = convert_form_merge_normal(this_pair[0], this_pair[1], one_element);
							if (!if_satisfied) continue;

							p_q.push(one_element);
							local_v_p_q.push_back(one_element);

						}

						else {

							p_q.push(one_p_q);
							local_v_p_q.push_back(one_p_q);


						}
					}
					else if (one_p_q.first.first == 2) {


						std::vector<int> this_reg = one_p_q.first.second;

						if (if_merged[this_reg[0]] == 0)
						{
							std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element_1;
							bool if_satisfied = convert_form_split_normal(this_reg[0], max_list_vector, min_list_vector, one_element_1);
							if (!if_satisfied) continue;

							p_q.push(one_element_1);
							local_v_p_q.push_back(one_element_1);


						}

						else {

							p_q.push(one_p_q);
							local_v_p_q.push_back(one_p_q);


						}

					}
					else if (one_p_q.first.first == 3) {


						std::vector<int> this_p_p = one_p_q.first.second;

						if (if_merged[this_p_p[0]] == 0) {

							std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element_b;
							bool if_satisfied_b = update_bad_points_normal(this_p_p[0], one_element_b);
							if (if_satisfied_b) {
								p_q.push(one_element_b);
								local_v_p_q.push_back(one_element_b);
							}
						}
						else {

							p_q.push(one_p_q);
							local_v_p_q.push_back(one_p_q);


						}

					}
					else if (one_p_q.first.first == 4) {


						std::vector<int> this_p_p = one_p_q.first.second;

						if (if_merged[this_p_p[0]] == 0) {

							std::pair<std::pair<int, std::vector<int>>, std::vector<double>> one_element_g;
							bool if_satisfied_g = update_good_points_normal(this_p_p[0], one_element_g);
							if (if_satisfied_g) {

								p_q.push(one_element_g);
								local_v_p_q.push_back(one_element_g);
							}

						}

						else {

							p_q.push(one_p_q);
							local_v_p_q.push_back(one_p_q);

						}

					}

				}


				v_p_q.clear();
				v_p_q = local_v_p_q;

			}

		}
		else {


			std::vector<int> inliers_to_planes_merge_local;
			inliers_to_planes_merge_local = std::vector<int>(points.size(), -1);
			for (int m = 0; m < planes_to_inliers.size(); ++m) {
				for (int k = 0; k < planes_to_inliers[m].size(); ++k) {
					inliers_to_planes_merge_local[planes_to_inliers[m][k]] = m;
				}


			}
			inliers_to_planes.clear();
			inliers_to_planes = inliers_to_planes_merge_local;
			bb = false;
		}
		if (t_m >= 1000) bb = false;
	
	} while (bb == true);

}




void Shape_Detector::show_result(double t) {
	std::cout << std::endl;
	std::cout << "***All operators: " << t << " s." << std::endl;
	std::cout << "**Transfer operator: " << all_t_transfer << " times" << std::endl;
	std::cout << "**Merge operator: " << all_t_merge << " times" << std::endl;

	std::cout << "**Split operator: " << all_t_split << " times" << std::endl;

	std::cout << "**Insert operator: " << all_t_insert << " times" << std::endl;
	std::cout << "**Exclude operator: " << all_t_exlude << " times" << std::endl;
	std::cout << "**number of iterations: " << number_iterations << std::endl;
	std::cout << "Primitives : " << (primitives_number) << std::endl;
	std::cout << "Coverage  : " << (coverage) << std::endl;
	std::cout << "Mean error : " << (mean_error) << std::endl;
	std::cout << "Mean normal deviation : " << (mean_normal_diviation) << std::endl;

	std::cout << "Primitives reducing : " << (last_primitives_number - primitives_number) << std::endl;
	std::cout << "Coverage adding   : " << (coverage - last_coverage) << std::endl;
	std::cout << "Mean error reducing : " << (last_mean_error - mean_error) << std::endl;
	std::cout << "Mean normal deviation adding : " << (mean_normal_diviation - last_normal_deviation) << std::endl;







}



void Shape_Detector::get_last_information() {
	last_coverage = coverage;
	last_mean_error = mean_error;
	last_primitives_number = primitives_number;
	last_normal_deviation = mean_normal_diviation;


}


void Shape_Detector::do_region_growing()
{
	int class_index = -1;
	int ten_percent_of_points = points.size() / 10;

	std::cout << "Region growing "<<points.size()<<std::endl;

	for (size_t i = 0; i < points.size(); i++) {
		if (i % ten_percent_of_points == 0) std::cout << ". " << std::flush;

		if (inliers_to_planes[i] == -1) {

			// Updates the index of primitive
			//inliers_to_planes[i] = ++class_index;
			++class_index;
			inliers_to_planes[i] = planes_to_inliers.size();

			int conti = 0; 	//for accelerate least_square fitting 

							// Characteristics of the seed
			Inexact_Point_3 pt_seed = points[i].first;
			Inexact_Vector_3 normal_seed = points[i].second;
			Inexact_Plane optimal_plane(pt_seed, normal_seed);

			// Initialization containers
			std::vector<int> index_container; index_container.push_back(i);
			std::vector<int> index_container_former_ring; index_container_former_ring.push_back(i);
			std::list<int> index_container_current_ring;

			// Propagation
			bool propagation = true;
			do {
				propagation = false;

				for (int k = 0; k < (int)index_container_former_ring.size(); k++) {

					int point_index = index_container_former_ring[k];

					int Nb_neigh = (int)spherical_neighborhood[point_index].size();

					if ((int)spherical_neighborhood[point_index].size() > knn) {
						Nb_neigh = knn;
					}
					for (int it = 0; it < Nb_neigh; it++) {

						int neighbor_index = spherical_neighborhood[point_index][it];

						if (inliers_to_planes[neighbor_index] == -1) {

							Inexact_Point_3 neighbor = points[neighbor_index].first;
							Inexact_Point_3 neighbor_projection = optimal_plane.projection(neighbor);
							double distance = sqrt((neighbor - neighbor_projection).squared_length());

							if (distance < epsilon && abs(points[neighbor_index].second * optimal_plane.orthogonal_vector()) > normal_threshold) {

								//inliers_to_planes[neighbor_index] = class_index;
								inliers_to_planes[neighbor_index] = planes_to_inliers.size();

								propagation = true;
								index_container_current_ring.push_back(neighbor_index);
								conti++;

								if ((conti < 50 && conti % 10 == 0) || (conti > 50 && conti % 500 == 0)) {
									std::list<Inexact_Point_3> listp;
									for (int pm = 0; pm < (int)index_container.size(); pm++) {
										Inexact_Point_3 ptza = points[index_container[pm]].first;
										listp.push_back(ptza);
									}

									if (listp.size() >= 3) {
										Inexact_Plane reajusted_plane;
										linear_least_squares_fitting_3(listp.begin(), listp.end(), reajusted_plane, CGAL::Dimension_tag<0>());
										optimal_plane = reajusted_plane;
									}
								}
							}
							/*else if (distance > epsilon && abs(points[neighbor_index].second * optimal_plane.orthogonal_vector()) > normal_threshold) {
								number_distance++;


							}
							else if (distance < epsilon && abs(points[neighbor_index].second * optimal_plane.orthogonal_vector()) < normal_threshold) {
								number_normal++;


							}
							else {
								number_distance++;
								number_normal++;


							}*/
						}
					}
				}

				// Updates containers
				index_container_former_ring.clear();
				for (std::list<int>::iterator it = index_container_current_ring.begin(); it != index_container_current_ring.end(); ++it) {
					index_container_former_ring.push_back(*it);
					index_container.push_back(*it);
				}
				index_container_current_ring.clear();

			} while (propagation);

			// Tests the number of inliers
			if (index_container.size() >= min_points && inliers_arent_aligned(index_container)) {
				planes_to_inliers.push_back(index_container);
			}

			else {
				--class_index;
				inliers_to_planes[i] = -1;
				for (int k = 0; k < (int)index_container.size(); k++) inliers_to_planes[index_container[k]] = -1;
			}
		}
	}
}


void Shape_Detector::detect_planes()
{
	clock_t t_detect_start = clock();

	// Part 1.
	// Performs a region-growing algorithm.

	// Initializes structures

	planes_0.clear();
	planes_1.clear();
	planes_2.clear();
	non_coplanar_planes = 0;
	planes_centroids.clear();
	planes_to_colors.clear();

	if (path_point_cloud_extension == ".ply"|| path_point_cloud_extension == ".off") {
		planes_to_inliers.clear();
		inliers_to_planes = std::vector<int>(points.size(), -1);

		
			if (should_compute_knn) compute_average_spacing_and_k_nearest_neighbors();
			do_region_growing();

		
		
	}
	else {
		compute_average_spacing_and_k_nearest_neighbors();
	}

	// Part 2.
	// Extracts planes

	
	std::default_random_engine generator;
	std::uniform_int_distribution<int> uniform_distribution(100, 225);
	unsigned char r = 0, g = 0, b = 0;
	for (size_t i = 0; i < planes_to_inliers.size(); ++i) {

		//Inexact_Point_3 centroid = CGAL::ORIGIN;
	

		std::vector<Inexact_Point_3> inliers_i;
		inliers_i.reserve(planes_to_inliers[i].size());

		for (size_t j = 0; j < planes_to_inliers[i].size(); ++j) {
			const Inexact_Point_3 & pt = points[planes_to_inliers[i][j]].first;
			inliers_i.push_back(pt);
			//centroid = CGAL::barycenter(centroid, j, pt, 1);
		
		}
		
	

		Inexact_Plane plane;
		linear_least_squares_fitting_3(inliers_i.begin(), inliers_i.end(), plane, CGAL::Dimension_tag<0>());
	

		planes_0.push_back(plane);
		planes_1.push_back(plane);
		planes_2.push_back(plane);

		
		r = unsigned char(uniform_distribution(generator));
		g = unsigned char(uniform_distribution(generator));
		b = unsigned char(uniform_distribution(generator));
		planes_to_colors.push_back(CGAL::Color(r, g, b));

	}

	// export_region_growing_results();

	clock_t t_detect_end = clock();
	std::cout << std::endl << "** Plane detection : done in " << double(t_detect_end - t_detect_start) / CLOCKS_PER_SEC << " s." << std::endl;
}

bool Shape_Detector::inliers_arent_aligned(const std::vector<int> & inds)
{
	size_t n = inds.size();

	for (size_t j = 0; j < n / 3; j++) {
		size_t id_0 = j, id_1 = n / 3 + j, id_2 = 2 * n / 3 + j;
		const Inexact_Point_3 & P_0 = points[id_0].first, &P_1 = points[id_1].first, &P_2 = points[id_2].first;

		Inexact_Vector_3 N = CGAL::cross_product(P_1 - P_0, P_2 - P_0);
		if (N != CGAL::NULL_VECTOR) {
			return true;
		}
	}

	return false;
}


void Shape_Detector::get_coverage_and_mean_error()
{
	number_of_assigned_points = 0;

	mean_error = 0;
	mean_normal_diviation = 0;

	for (size_t i = 0; i < planes_to_inliers.size(); ++i) {
		number_of_assigned_points += planes_to_inliers[i].size();
		const Inexact_Plane & H = planes_2[i];
		for (int j : planes_to_inliers[i]) {
			const Point_with_normal & pt = points[j];
			mean_error += sqrt(CGAL::squared_distance(H, pt.first));
			mean_normal_diviation += abs(pt.second * H.orthogonal_vector());

		}
	}
	all_error = mean_error;
	mean_error /= double(number_of_assigned_points);
	mean_normal_diviation /= double(number_of_assigned_points);


	coverage = double(number_of_assigned_points) / double(points.size());
	primitives_number = planes_to_inliers.size();
	std::cout << "Primitives : " << primitives_number << std::endl;
	std::cout << "Coverage   : " << coverage << std::endl;
	std::cout << "Mean error : " << mean_error << std::endl;
	std::cout << "Mean normal deviation     : " << mean_normal_diviation << std::endl;

}

void Shape_Detector::get_coverage_and_mean_error_pure()
{
	number_of_assigned_points = 0;

	mean_error = 0;
	mean_normal_diviation = 0;
	for (size_t i = 0; i < planes_to_inliers.size(); ++i) {
		number_of_assigned_points += planes_to_inliers[i].size();
		const Inexact_Plane & H = planes_2[i];
		for (int j : planes_to_inliers[i]) {
			const Point_with_normal & pt = points[j];
			mean_error += sqrt(CGAL::squared_distance(H, pt.first));
			mean_normal_diviation += abs(H.orthogonal_vector()*pt.second);
		}
	}
	mean_error /= double(number_of_assigned_points);
	mean_normal_diviation /= double(number_of_assigned_points);
	coverage = double(number_of_assigned_points) / double(points.size());
	primitives_number = planes_to_inliers.size();

}


void Shape_Detector::get_distance_diviation() {
	number_of_assigned_points = 0;
	mean_distance_diaviation = 0;
	size_current_primitives = 0;
	mean_distance_current = 0;
	mean_normal_current = 0;
	mean_normal_diviation = 0;
	for (size_t i = 0; i < planes_to_inliers.size(); ++i) {
		number_of_assigned_points += planes_to_inliers[i].size();
		const Inexact_Plane & H = planes_2[i];


		for (int j : planes_to_inliers[i]) {
			const Point_with_normal & pt = points[j];

			mean_distance_diaviation += sqrt(CGAL::squared_distance(H, pt.first));
			mean_normal_diviation += abs(pt.second * H.orthogonal_vector());
		}
	}
	size_current_primitives = planes_to_inliers.size();
	all_distance_diaviation = mean_distance_diaviation;
	all_normal_diaviation = mean_normal_diviation;
	mean_distance_diaviation /= double(number_of_assigned_points);
	mean_normal_diviation /= double(number_of_assigned_points);
	mean_distance_current = mean_distance_diaviation;
	mean_normal_current = mean_normal_diviation;

}

void Shape_Detector::get_distance_diviation_show_normal_info(double t) {


	number_of_assigned_points = 0;
	mean_distance_diaviation = 0;
	size_current_primitives = 0;
	mean_distance_current = 0;
	mean_normal_current = 0;
	mean_normal_diviation = 0;




	for (size_t i = 0; i < planes_to_inliers.size(); ++i) {
		number_of_assigned_points += planes_to_inliers[i].size();
		const Inexact_Plane & H = planes_2[i];


		for (int j : planes_to_inliers[i]) {
			const Point_with_normal & pt = points[j];

			mean_distance_diaviation += sqrt(CGAL::squared_distance(H, pt.first));
			mean_normal_diviation += abs(pt.second * H.orthogonal_vector());

		}
	}
	size_current_primitives = planes_to_inliers.size();
	all_distance_diaviation = mean_distance_diaviation;
	all_normal_diaviation = mean_normal_diviation;


	mean_distance_diaviation /= double(number_of_assigned_points);
	mean_normal_diviation /= double(number_of_assigned_points);
	mean_normal_current = mean_normal_diviation;

	mean_distance_current = mean_distance_diaviation;
	double new_coverage = double(number_of_assigned_points) / double(points.size());



	std::cout << std::endl;


	std::cout << "** Transfer operator: " << t << " s." << std::endl;
	std::cout << "Primitives : " << (size_current_primitives) << std::endl;

	std::cout << "Coverage  : " << (new_coverage) << std::endl;
	std::cout << "Mean error : " << (mean_distance_diaviation) << std::endl;
	std::cout << "Mean normal deviation : " << (mean_normal_diviation) << std::endl;

	std::cout << "Primitives reducing : " << (old_size_current_primitives - size_current_primitives) << std::endl;

	std::cout << "Coverage adding   : " << new_coverage - old_coverage << std::endl;
	std::cout << "Mean error reducing : " << (old_mean_distance_diaviation - mean_distance_diaviation) << std::endl;
	std::cout << "Mean normal deviation adding : " << (mean_normal_diviation - old_mean_normal_diaviation) << std::endl;

	old_size_current_primitives = size_current_primitives;
	old_mean_distance_diaviation = mean_distance_diaviation;
	old_coverage = new_coverage;
	old_mean_normal_diaviation = mean_normal_diviation;




}

void Shape_Detector::get_distance_diviation_show_merge_info(double t) {


	number_of_assigned_points = 0;
	mean_distance_diaviation = 0;
	size_current_primitives = 0;
	mean_distance_current = 0;
	mean_normal_current = 0;
	mean_normal_diviation = 0;

	for (size_t i = 0; i < planes_to_inliers.size(); ++i) {
		number_of_assigned_points += planes_to_inliers[i].size();
		const Inexact_Plane & H = planes_2[i];


		for (int j : planes_to_inliers[i]) {
			const Point_with_normal & pt = points[j];

			mean_distance_diaviation += sqrt(CGAL::squared_distance(H, pt.first));
			mean_normal_diviation += abs(pt.second * H.orthogonal_vector());

		}
	}
	size_current_primitives = planes_to_inliers.size();
	all_distance_diaviation = mean_distance_diaviation;
	all_normal_diaviation = mean_normal_diviation;

	mean_distance_diaviation /= double(number_of_assigned_points);
	mean_normal_diviation /= double(number_of_assigned_points);
	mean_normal_current = mean_normal_diviation;

	mean_distance_current = mean_distance_diaviation;
	double new_coverage = double(number_of_assigned_points) / double(points.size());

	std::cout << std::endl;
	std::cout << "** Local operators: " << t << " s." << std::endl;
	std::cout << "* Merge operators: " << t_merge << " times." << std::endl;
	std::cout << "* Split operators: " << t_split << " times." << std::endl;
	std::cout << "* Insert operators: " << t_insert << " times." << std::endl;
	std::cout << "* Exclude operators: " << t_exlude << " times." << std::endl;

	std::cout << "Primitives : " << (size_current_primitives) << std::endl;
	
	std::cout << "Coverage  : " << (new_coverage) << std::endl;
	std::cout << "Mean error : " << (mean_distance_diaviation) << std::endl;
	std::cout << "Mean normal deviation : " << (mean_normal_diviation) << std::endl;

	std::cout << "Primitives reducing : " << (old_size_current_primitives - size_current_primitives) << std::endl;

	std::cout << "Coverage adding   : " << new_coverage - old_coverage << std::endl;
	std::cout << "Mean error reducing : " << (old_mean_distance_diaviation - mean_distance_diaviation) << std::endl;
	std::cout << "Mean normal deviation adding : " << (mean_normal_diviation - old_mean_normal_diaviation) << std::endl;

	old_size_current_primitives = size_current_primitives;
	old_mean_distance_diaviation = mean_distance_diaviation;
	old_coverage = new_coverage;
	old_mean_normal_diaviation = mean_normal_diviation;



}


double Shape_Detector::add_changed_error(int id, std::vector<int> point_ids) {
	std::vector<int> befor_ids = planes_to_inliers[id];
	Inexact_Plane befor_p = planes_2[id];

	double dif_befor = 0;
	for (int u : befor_ids) {
		dif_befor += sqrt(CGAL::squared_distance(befor_p, points[u].first));
	}
	befor_ids.insert(befor_ids.end(), point_ids.begin(), point_ids.end());
	std::vector<Inexact_Point_3> after_points;
	for (int n : befor_ids) {
		after_points.push_back(points[n].first);
	}
	Inexact_Plane after_p;

	linear_least_squares_fitting_3(after_points.begin(), after_points.end(), after_p, CGAL::Dimension_tag<0>());

	double dif_after = 0;
	for (Inexact_Point_3 u : after_points) {
		dif_after += sqrt(CGAL::squared_distance(after_p, u));
	}

	return (dif_after - dif_befor);


}
void Shape_Detector::get_good_points() {
	if_insert_candidate = std::vector<int>(points.size(), -1);
	good_points_shape.clear();
	int Nb_neigh = 10;
	std::map<int, std::vector<std::pair<int, double>>> plane_outliers;
	for (int i = 0; i < points.size(); ++i) {
		if (inliers_to_planes[i] != -1) continue;
		if (points_if_added[i] > 2) continue;

		if ((int)spherical_neighborhood[i].size() < Nb_neigh) {
			Nb_neigh = (int)spherical_neighborhood[i].size();
		}

		std::set<int> one_neight_id;



		for (int it = 0; it < Nb_neigh; ++it) {

			int neighbor_index = spherical_neighborhood[i][it];

			if (inliers_to_planes[neighbor_index] != -1) {

				one_neight_id.insert(inliers_to_planes[neighbor_index]);

			}


		}

		if (one_neight_id.empty()) { continue; }

		double mm_dis = epsilon;

		int changed_plane_id = -1;
		Point_with_normal this_p = points[i];

		for (int neight_id : one_neight_id) {

			if (abs(this_p.second * planes_2[neight_id].orthogonal_vector()) > normal_threshold) {
				if (sqrt(CGAL::squared_distance(planes_2[neight_id], this_p.first)) < mm_dis) {

					mm_dis = sqrt(CGAL::squared_distance(planes_2[neight_id], this_p.first));
					changed_plane_id = neight_id;

				}
			}

		}

		if (changed_plane_id != -1) {
			plane_outliers[changed_plane_id].push_back(std::make_pair(i, mm_dis));//plane_id, point_id, distance.

		}

	}
	std::map<int, std::vector<std::pair<int, double>>>::iterator g_i;
	g_i = plane_outliers.begin();


	while (g_i != plane_outliers.end()) {
		std::priority_queue<std::pair<int, double>, std::vector<std::pair<int, double>>, Add_Comparator> p_q;

		for (std::pair<int, double> ppp : g_i->second) {
			p_q.push(ppp);
		}

		std::vector<int> added_points_ids;
		double en = 0;
		double dif_e = 0;
		double n_a = 0;
		while (!p_q.empty() && n_a < number_of_insert_exclude) {
			added_points_ids.push_back(p_q.top().first);
			p_q.pop();
			n_a++;
		}
		dif_e = add_changed_error(g_i->first, added_points_ids);
		en = energy_changed_second(dif_e, -n_a);
		if (added_points_ids.size() != 0 && en < 0) {
			for (int idd : added_points_ids) {
				if_insert_candidate[idd] = 0;
			}
			std::vector<double> dd;
			dd.push_back(en);
			dd.push_back(dif_e);
			good_points_shape.push_back(std::make_pair(std::make_pair(g_i->first, added_points_ids), dd));
		}
		g_i++;
	}

}
double Shape_Detector::add_changed_error_normal(int id, std::vector<int> point_ids) {
	std::vector<int> befor_ids = planes_to_inliers[id];
	Inexact_Plane befor_p = planes_2[id];

	double dif_befor = 0;
	for (int u : befor_ids) {

		dif_befor += abs(befor_p.orthogonal_vector()*points[u].second);

	}
	befor_ids.insert(befor_ids.end(), point_ids.begin(), point_ids.end());
	std::vector<Inexact_Point_3> after_points;

	for (int n : befor_ids) {
		after_points.push_back(points[n].first);

	}
	Inexact_Plane after_p;

	linear_least_squares_fitting_3(after_points.begin(), after_points.end(), after_p, CGAL::Dimension_tag<0>());

	double dif_after = 0;
	for (int u : befor_ids) {
		dif_after += abs(after_p.orthogonal_vector()*points[u].second);
	}

	return (dif_after - dif_befor);


}
void Shape_Detector::get_good_points_normal() {
	if_insert_candidate = std::vector<int>(points.size(), -1);
	good_points_shape.clear();
	int Nb_neigh = 10;
	std::map<int, std::vector<std::pair<int, double>>> plane_outliers;
	for (int i = 0; i < points.size(); ++i) {
		if (inliers_to_planes[i] != -1) continue;
		if (points_if_added[i] > 2) continue;
		if ((int)spherical_neighborhood[i].size() < Nb_neigh) {
			Nb_neigh = (int)spherical_neighborhood[i].size();
		}

		std::set<int> one_neight_id;


		for (int it = 0; it < Nb_neigh; ++it) {

			int neighbor_index = spherical_neighborhood[i][it];

			if (inliers_to_planes[neighbor_index] != -1) {

				one_neight_id.insert(inliers_to_planes[neighbor_index]);

			}


		}

		if (one_neight_id.empty()) { continue; }

		double mm_normal = normal_threshold;

		int changed_plane_id = -1;
		Point_with_normal this_p = points[i];

		for (int neight_id : one_neight_id) {
			if (abs(this_p.second * planes_2[neight_id].orthogonal_vector()) > mm_normal) {
				if (sqrt(CGAL::squared_distance(planes_2[neight_id], this_p.first)) < epsilon) {

					mm_normal = abs(this_p.second * planes_2[neight_id].orthogonal_vector());
					changed_plane_id = neight_id;

				}

			}
		}


		if (changed_plane_id != -1) {
			plane_outliers[changed_plane_id].push_back(std::make_pair(i, mm_normal));//plane_id, point_id, normal_deviation.

		}

	}
	std::map<int, std::vector<std::pair<int, double>>>::iterator g_i;
	g_i = plane_outliers.begin();
	while (g_i != plane_outliers.end()) {
		std::priority_queue<std::pair<int, double>, std::vector<std::pair<int, double>>, Add_Comparator_normal> p_q;

		for (std::pair<int, double> ppp : g_i->second) {
			p_q.push(ppp);
		}

		std::vector<int> added_points_ids;
		double en = 0;
		double dif_e = 0;
		double n_a = 0;

		while (!p_q.empty() && n_a < number_of_insert_exclude) {
			added_points_ids.push_back(p_q.top().first);
			p_q.pop();
			n_a++;
		}
		dif_e = add_changed_error_normal(g_i->first, added_points_ids);
		en = energy_changed_second_normal(dif_e, -n_a);

		if (added_points_ids.size() != 0 && en < 0) {
			for (int idd : added_points_ids) {
				if_insert_candidate[idd] = 0;
			}
			std::vector<double> dd;
			dd.push_back(en);
			dd.push_back(dif_e);
			good_points_shape.push_back(std::make_pair(std::make_pair(g_i->first, added_points_ids), dd));//plane_id,point_ids,energy_change,distance
		}
		g_i++;
	}

}


double Shape_Detector::remove_changed_error(int id, std::vector<int> point_ids) {
	std::vector<int> befor_ids = planes_to_inliers[id];
	Inexact_Plane befor_p = planes_2[id];

	double dif_befor = 0;
	std::vector<int> after_ids;
	for (int u : befor_ids) {
		dif_befor += sqrt(CGAL::squared_distance(befor_p, points[u].first));
		/*bool sma = false;
		for (int t : point_ids) {
			if (t == u) {
				sma = true;
				break;
			}
		}
		if (!sma)after_ids.push_back(u);*/
		if (std::find(point_ids.begin(), point_ids.end(), u) == point_ids.end()) {
			after_ids.push_back(u);

		}
	}





	//befor_ids.insert(befor_ids.end(), point_ids.begin(), point_ids.end());







	std::vector<Inexact_Point_3> after_points;

	for (int n : after_ids) {
		after_points.push_back(points[n].first);
	}
	Inexact_Plane after_p;

	linear_least_squares_fitting_3(after_points.begin(), after_points.end(), after_p, CGAL::Dimension_tag<0>());

	double dif_after = 0;
	for (Inexact_Point_3 u : after_points) {
		dif_after += sqrt(CGAL::squared_distance(after_p, u));
	}

	return (dif_after - dif_befor);


}

void Shape_Detector::get_bad_points() {

	bad_points_shape.clear();
	if (old_coverage > ori_coverage) {
		std::map<int, std::vector<std::pair<int, double>>> plane_inliers;
		double bigest_dis = 1.5 * mean_distance_diaviation;

		for (size_t i = 0; i < planes_to_inliers.size(); ++i) {

			const Inexact_Plane & H = planes_2[i];


			for (int j : planes_to_inliers[i]) {
				const Point_with_normal & pt = points[j];
				if (points_if_added[j] > 2) continue;
				if (sqrt(CGAL::squared_distance(H, pt.first)) > bigest_dis) {
					plane_inliers[i].push_back(std::make_pair(j, sqrt(CGAL::squared_distance(H, pt.first))));
				}

			}
		}

		std::map<int, std::vector<std::pair<int, double>>>::iterator b_o;
		b_o = plane_inliers.begin();

		while (b_o != plane_inliers.end()) {


			std::priority_queue<std::pair<int, double>, std::vector<std::pair<int, double>>, Remove_Comparator> p_q;
			for (std::pair<int, double> ppp : b_o->second) {
				p_q.push(ppp);

			}
			std::vector<int> removed_points_ids;
			double en = 0;
			double dif_e = 0;

			double n_b = 0;
			while (!p_q.empty() && n_b < number_of_insert_exclude) {
				removed_points_ids.push_back(p_q.top().first);


				p_q.pop();
				n_b++;
			}
			dif_e = remove_changed_error(b_o->first, removed_points_ids);
			en = energy_changed_second(dif_e, n_b);

			if (removed_points_ids.size() != 0 && (planes_to_inliers[b_o->first].size() - removed_points_ids.size()) >= min_points && en < 0) {
				std::vector<double> dd;
				dd.push_back(en);
				dd.push_back(dif_e);
				bad_points_shape.push_back(std::make_pair(std::make_pair(b_o->first, removed_points_ids), dd));
			}
			b_o++;
		}

	}
}

double Shape_Detector::remove_changed_error_normal(int id, std::vector<int> point_ids) {
	std::vector<int> befor_ids = planes_to_inliers[id];
	Inexact_Plane befor_p = planes_2[id];

	double dif_befor = 0;
	std::vector<int> after_ids;
	for (int u : befor_ids) {
		//dif_befor += sqrt(CGAL::squared_distance(befor_p, points[u].first));
		dif_befor += abs(points[u].second * befor_p.orthogonal_vector());

		if (std::find(point_ids.begin(), point_ids.end(), u) == point_ids.end()) {
			after_ids.push_back(u);

		}
	}




	std::vector<Inexact_Point_3> after_points;
	std::vector<int> after_points_ids;
	for (int n : after_ids) {
		after_points.push_back(points[n].first);
		after_points_ids.push_back(n);
	}
	Inexact_Plane after_p;

	linear_least_squares_fitting_3(after_points.begin(), after_points.end(), after_p, CGAL::Dimension_tag<0>());

	double dif_after = 0;
	for (int u : after_points_ids) {
		//dif_after += sqrt(CGAL::squared_distance(after_p, u));
		dif_after += abs(points[u].second * after_p.orthogonal_vector());
	}

	return (dif_after - dif_befor);


}

void Shape_Detector::get_bad_points_normal() {

	bad_points_shape.clear();
	if (old_coverage > ori_coverage) {
		std::map<int, std::vector<std::pair<int, double>>> plane_inliers;
		double bigest_normal = mean_normal_diviation;
		//double bigest_normal = 0.01;

		for (size_t i = 0; i < planes_to_inliers.size(); ++i) {

			const Inexact_Plane & H = planes_2[i];


			for (int j : planes_to_inliers[i]) {
				const Point_with_normal & pt = points[j];
				if (points_if_added[j] > 2) continue;

				if (abs(H.orthogonal_vector()*pt.second) < bigest_normal) {



					plane_inliers[i].push_back(std::make_pair(j, abs(H.orthogonal_vector()*pt.second)));
				}

			}
		}

		std::map<int, std::vector<std::pair<int, double>>>::iterator b_o;
		b_o = plane_inliers.begin();

		while (b_o != plane_inliers.end()) {


			std::priority_queue<std::pair<int, double>, std::vector<std::pair<int, double>>, Remove_Comparator_normal> p_q;
			for (std::pair<int, double> ppp : b_o->second) {
				p_q.push(ppp);

			}
			std::vector<int> removed_points_ids;
			double en = 0;
			double dif_e = 0;

			double n_b = 0;
			while (!p_q.empty() && n_b < number_of_insert_exclude) {
				removed_points_ids.push_back(p_q.top().first);
				p_q.pop();
				n_b++;
			}
			dif_e = remove_changed_error_normal(b_o->first, removed_points_ids);
			en = energy_changed_second_normal(dif_e, n_b);

			if (removed_points_ids.size() != 0 && (planes_to_inliers[b_o->first].size() - removed_points_ids.size()) >= min_points && en < 0) {
				std::vector<double> dd;
				dd.push_back(en);
				dd.push_back(dif_e);
				bad_points_shape.push_back(std::make_pair(std::make_pair(b_o->first, removed_points_ids), dd));
			}
			b_o++;
		}

	}
}


void Shape_Detector::set_primitives_simple()
{
	

	alpha_shapes_pts.clear();
	
	convex_hulls_pts.clear();

	convex_hulls_seeds.clear();
	alpha_shapes_colors.clear();
	convex_hulls_colors.clear();
	
	clock_t t_start = clock();

	std::vector<Inexact_Plane> planes_3;
	std::default_random_engine generator;
	std::uniform_int_distribution<int> uniform_distribution(100, 225);
	unsigned char r = 0, g = 0, b = 0;

	
	for (size_t u = 0; u < planes_2.size(); ++u) {



		Inexact_Plane H = planes_2[u];
		r = unsigned char(uniform_distribution(generator));
		g = unsigned char(uniform_distribution(generator));
		b = unsigned char(uniform_distribution(generator));

		const CGAL::Color & col = (CGAL::Color(r, g, b));

		std::vector<Inexact_Point_3> assigned_pts;
		assigned_pts.reserve(planes_to_inliers[u].size());
		for (int pt_index : planes_to_inliers[u]) {
			assigned_pts.push_back(points[pt_index].first);
		}

		// Computes alpha-shapes.

		Inexact_Vector_3 n = H.orthogonal_vector(), b1 = H.base1(), b2 = H.base2();
		n = n * 1.0 / (sqrt(n*n));
		b1 = b1 * 1.0 / (sqrt(b1*b1));
		b2 = b2 * 1.0 / (sqrt(b2*b2));

		std::vector<Inexact_Point_2> p2d;
		std::map<Inexact_Point_2, int> vMap, vMap2;
		std::vector<Inexact_Point_3> p3d;

		auto it = assigned_pts.begin();
		Inexact_Point_3 O = H.projection(*it);
		p2d.push_back(Inexact_Point_2(0, 0));
		vMap[p2d.back()] = p2d.size() - 1;
		it++;
		while (it != assigned_pts.end()) {
			Inexact_Vector_3 p = (*it) - O;
			p2d.push_back(Inexact_Point_2(p * b1, p * b2));
			vMap[p2d.back()] = p2d.size() - 1;
			it++;
		}

		Alpha_Shape as(p2d.begin(), p2d.end());
		as.set_alpha(double(0.005 * bbox_diagonal));

		int ind_min = alpha_shapes_pts.size();

		std::list<Inexact_Point_3> las3d;

		CDT2 cdt;
		std::map<Inexact_Point_2, CDT2::Vertex_handle> pVMap;

		auto vit = as.alpha_shape_vertices_begin();
		while (vit != as.alpha_shape_vertices_end()) {
			CDT2::Vertex_handle v = cdt.insert((*vit)->point());
			pVMap[(*vit)->point()] = v;

			vit++;
		}


		auto eit = as.alpha_shape_edges_begin();
		while (eit != as.alpha_shape_edges_end()) {
			switch (as.classify(*eit)) {
			case Alpha_Shape::SINGULAR:
				eit++;
				continue;
			default:
				break;
			}

			CDT2::Vertex_handle v1 = pVMap[as.segment(*eit).source()];
			CDT2::Vertex_handle v2 = pVMap[as.segment(*eit).target()];
			if (!v1->is_valid() || !v2->is_valid())
				std::cout << "invalid!" << std::endl;

			eit++;
		}

		int outside = 0;

		auto fit = cdt.finite_faces_begin();
		while (fit != cdt.finite_faces_end()) {
			int res;
			Inexact_Point_2 p1 = fit->vertex(0)->point();
			Inexact_Point_2 p2 = fit->vertex(1)->point();
			Inexact_Point_2 p3 = fit->vertex(2)->point();
			if (CGAL::collinear(p1, p2, p3))
				continue;
			Inexact_Point_2 center = Inexact_Point_2(0, 0) + (((p1 - Inexact_Point_2(0, 0)) + (p2 - Inexact_Point_2(0, 0)) + (p3 - Inexact_Point_2(0, 0))) * 1.0 / 3.0);
			res = as.classify(center);
			if (res == Alpha_Shape::INTERIOR) {
			
				for (int i = 0; i < 3; i++) {
					Inexact_Point_2 p = fit->vertex(i)->point();

					// if point is used in cdt, copy it to 
					if (vMap2.count(p) == 0) {
						p3d.push_back(H.projection(assigned_pts[vMap[p]]));

					

						vMap2[p] = vMap[p];
						vMap[p] = p3d.size() - 1;
					}
				}

			}
			else outside++;
			fit++;
		}

		fit = cdt.finite_faces_begin();
		while (fit != cdt.finite_faces_end()) {
			int res;
			//res = as.classify(fit);
			Inexact_Point_2 p1 = fit->vertex(0)->point();
			Inexact_Point_2 p2 = fit->vertex(1)->point();
			Inexact_Point_2 p3 = fit->vertex(2)->point();
			Inexact_Point_2 center = Inexact_Point_2(0, 0) + (((p1 - Inexact_Point_2(0, 0)) + (p2 - Inexact_Point_2(0, 0)) + (p3 - Inexact_Point_2(0, 0))) * 1.0 / 3.0);
			res = as.classify(center);
			if (res == Alpha_Shape::INTERIOR) {
				if (vMap.count(p1) == 0) {
					continue;
				}
				if (vMap.count(p2) == 0) {
					continue;
				}
				if (vMap.count(p3) == 0) {
					continue;
				}

				std::vector<Inexact_Point_3> tr(3);
				tr[0] = p3d[vMap[p1]];
				tr[1] = p3d[vMap[p2]];
				tr[2] = p3d[vMap[p3]];
				alpha_shapes_pts.push_back(tr);
				alpha_shapes_colors.push_back(col);

				las3d.push_back(tr[0]);
				las3d.push_back(tr[1]);
				las3d.push_back(tr[2]);

			}
			fit++;
		}

		int ind_max = alpha_shapes_pts.size() - 1;

		// Computes convex hulls.

		std::vector<Inexact_Point_3> as3d(las3d.begin(), las3d.end());
		std::vector<Inexact_Point_2> as2d(as3d.size());

		for (size_t j = 0; j < as3d.size(); ++j) {
			as2d[j] = H.to_2d(as3d[j]);


		}

		std::vector<Inexact_Point_2> ch2d;
		CGAL::ch_graham_andrew(as2d.begin(), as2d.end(), std::back_inserter(ch2d));

		std::vector<Inexact_Point_3> ch3d;
		ch3d.reserve(ch2d.size());
		for (size_t j = 0; j < ch2d.size(); ++j) {
			ch3d.push_back(H.to_3d(ch2d[j]));
		}
		convex_hulls_seeds.push_back(std::make_pair(ind_min, ind_max));

		convex_hulls_pts.push_back(ch3d);

		convex_hulls_colors.push_back(col);

	}
}



void Shape_Detector::set_constraint(bool cc) {
	
	if_constraint = cc;
}



double Shape_Detector::get_bbox_diagonal() const
{
	return bbox_diagonal;
}


double Shape_Detector::get_average_spacing() const
{
	return average_spacing;
}



double Shape_Detector::get_epsilon() const
{
	return epsilon;
}


int Shape_Detector::get_knn() const
{
	return knn;
}


double Shape_Detector::get_normal_threshold() const
{
	return normal_threshold;
}


int Shape_Detector::get_neighborhood_algorithm() const
{
	return neighborhood_algorithm;
}


double Shape_Detector::get_neighborhood_distance() const
{
	return neighborhood_distance;
}











