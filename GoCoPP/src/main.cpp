
#include "shape_detector.h"

#include <boost/filesystem.hpp>


int main(int argc, char *argv[])
{
	
	double pd_epsilon = 0;
	char parameter[1024];
	int pd_sigma = 50;
	int stop_iterations = 7;
	double pd_normal_deviation = 0.85;
	int pd_nn = 20;
	int rr = 2;
	double pd_c = 1;
	double pd_s = 1;
	std::string pd_norm = "hybrid";
	std::string path_point_cloud = argv[1];
	bool pd_out_vg = false;
	bool pd_out_alpha_shape = true;
	bool pd_out_convex_hull = false;
	
	
	//std::string path_point_cloud = "D:/planar_shape_detection_command/build/bin/Release/torus-100K.ply";

	bool pd_if_constraint = false;
	int pd_weight_mode = 0;
	while (rr < argc) {
		if (!strcmp(argv[rr], "--epsilon") && rr + 1 < argc) {
			pd_epsilon = atof(argv[rr + 1]);
			

			rr += 2;
		}
		else if (!strcmp(argv[rr], "--sigma") && rr + 1 < argc) {
			pd_sigma = atoi(argv[rr + 1]);


			rr += 2;
		}
		else if (!strcmp(argv[rr], "--nn") && rr + 1 < argc) {
			pd_nn = atoi(argv[rr + 1]);


			rr += 2;
		}
		else if (!strcmp(argv[rr], "--normal_deviation") && rr + 1 < argc) {
			pd_normal_deviation = atof(argv[rr + 1]);


			rr += 2;
		}
		else if (!strcmp(argv[rr], "--s") && rr + 1 < argc) {
			pd_s = atof(argv[rr + 1]);


			rr += 2;
		}
		else if (!strcmp(argv[rr], "--c") && rr + 1 < argc) {
			pd_c = atof(argv[rr + 1]);


			rr += 2;
		}
		
		else if (!strcmp(argv[rr], "--norm") && rr + 1 < argc) {
			pd_norm = argv[rr + 1];

			
	
			rr += 2;
		}
		else if (!strcmp(argv[rr], "--weight_mode") && rr + 1 < argc) {
			pd_weight_mode = atoi(argv[rr + 1]);


			rr += 2;
		}
		else if (!strcmp(argv[rr], "--max_steps") && rr + 1 < argc) {
			stop_iterations = atoi(argv[rr + 1]);


			rr += 2;
		}
		else if (!strcmp(argv[rr], "--constraint")) {
			pd_if_constraint = true;


			rr += 1;
		}
		else if (!strcmp(argv[rr], "--vg")) {
			pd_out_vg = true;


			rr += 1;
		}
		else if (!strcmp(argv[rr], "--alpha")) {
			pd_out_alpha_shape = true;


			rr += 1;
		}
		else if (!strcmp(argv[rr], "--hull")) {
			pd_out_convex_hull = true;


			rr += 1;
		}
		
		else {
			rr++;
		}
		
	}
	
	Shape_Detector* CS;

	
	

	std::string path_point_cloud_basename = boost::filesystem::path(path_point_cloud).stem().string();
	if (!boost::filesystem::exists(path_point_cloud)) {
		std::cout << "Can not find the file." << std::endl << std::endl;
			
		return 0;
	}

	CS = new Shape_Detector(path_point_cloud);
	CS->set_detection_parameters(pd_sigma,  pd_nn, pd_normal_deviation);
	
	CS->set_max_steps(stop_iterations);
	CS->set_weight_m(pd_weight_mode);
	CS->set_constraint(pd_if_constraint);
	
	CS->set_lambda_r(pd_s);
	CS->set_lambda_c(pd_c);



	if (!CS->load_points()) {
		std::stringstream error_message;
		error_message << "The selected point cloud couldn't be processed." << std::endl << std::endl
			<< "Please make that the file doesn't have any unhandled property. The expected ones are : x, y, z, nx, ny, nz." << std::endl;

		delete CS;
		return 0;
	}



	if (pd_epsilon == 0) {
		double pd_epsilon_0 = 0.004 * CS->get_bbox_diagonal();
		pd_epsilon = (round(1000 * pd_epsilon_0)) / 1000;
	}

	CS->set_epsilon(pd_epsilon);
	std::cout << "**Epsilon: " << pd_epsilon << std::endl;
	std::cout << "**Sigma: " << pd_sigma << std::endl;
	std::cout << "**KNN: " << pd_nn << std::endl;
	std::cout << "**Normal threshold: " << pd_normal_deviation << std::endl;
	std::cout << "**Norm: " << pd_norm << std::endl;
	std::cout << "**If Constraint: " << pd_if_constraint << std::endl;
	

	
	CS->detect_shapes();
	CS->set_primitives_simple();


	
	
	if (pd_norm == "normal") {
		clock_t t_start = clock();

		
		CS->planar_shape_detection_L1();


		clock_t t_end = clock();

		double t_all = double(t_end - t_start) / CLOCKS_PER_SEC;
		CS->show_result(t_all);
		CS->set_primitives_simple();
		if (pd_out_convex_hull) {
			CS->save_convex_hull();
		}
		if (pd_out_alpha_shape) {
			CS->save_alpha_shapes();
		}
		if (pd_out_vg) {
			CS->to_vg();
		}
	}
	else if (pd_norm == "L2") {

		

		clock_t t_start = clock();


		CS->planar_shape_detection_l2();



		clock_t t_end = clock();

		double t_all = double(t_end - t_start) / CLOCKS_PER_SEC;
		CS->show_result(t_all);
		CS->set_primitives_simple();
		if (pd_out_convex_hull) {
			CS->save_convex_hull();
		}
		if (pd_out_alpha_shape) {
			CS->save_alpha_shapes();
		}
		if (pd_out_vg) {
			CS->to_vg();
		}

	}
	else if (pd_norm == "hybrid") {
		

		clock_t t_start = clock();


		CS->planar_shape_detection_hybrid();



		clock_t t_end = clock();

		double t_all = double(t_end - t_start) / CLOCKS_PER_SEC;
		CS->show_result(t_all);

		CS->set_primitives_simple();
		if (pd_out_convex_hull) {
			CS->save_convex_hull();
		}
		if (pd_out_alpha_shape) {
			CS->save_alpha_shapes();
		}
		if (pd_out_vg) {
			CS->to_vg();
		}

	}
	return 1;
	
	
	
}
