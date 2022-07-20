# [Finding Good Configurations of Planar Primitives in Unorganized Point Clouds](https://hal.inria.fr/hal-03621896)

[Mulin Yu](http://www-sop.inria.fr/members/Mulin.Yu/), [Florent Lafarge](http://www-sop.inria.fr/members/Florent.Lafarge/) 

Inria - Université Côte d'Azur

firstname.lastname@inria.fr

## Citation

```
@INPROCEEDINGS{{Yu_cvpr22,
  Author = {Yu, Mulin and Lafarge, Florent},
  Title = {Finding Good Configurations of Planar Primitives in Unorganized Point Clouds},
  booktitle = {Proc. of the IEEE conference on Computer Vision and Pattern Recognition (CVPR)},
  Year = {2022},
  address = {New Orleans, US},
}
```

## Introduction

We present an algorithm for detecting planar primitives from unorganized 3D point clouds. Departing from an initial configuration, the algorithm refines both the continuous plane parameters and the discrete assignment of input points to them by seeking high fidelity, high simplicity and high completeness. Our key contribution relies upon the design of an exploration mechanism guided by a multiobjective energy function. The transitions within the large solution space are handled by five geometric operators that create, remove and modify primitives. We demonstrate the potential of our method on a variety of scenes, from organic shapes to man-made objects, and sensors, from multiview stereo to laser. We show its efficacy with respect to existing primitive fitting approaches and illustrate its applicative interest in compact mesh reconstruction, when combined with a plane assembly method.


## Requirement

CGAL 5.2.2

Eigen

Boost_1_76_0

CMake

Visual Studio 2017

## Compile Guideline

Open sources/CMakeLists.txt in cmake-gui and fill in all the relevant fields. This will generate a Visual Studio solution in folder “build”. Open this solution and compile the projects `GoCoPP`. The output executable files are created in subfolder “build/bin”.


## Input

The input data is a [PLY](https://en.wikipedia.org/wiki/PLY_(file_format)) file. It should contain the coordinates of input points (x y z). The normals of input points (nx ny nz) are also encouraged to be included. If not included, they will be estimated by [PCA](https://doc.cgal.org/5.2.4/Point_set_processing_3/index.html) just after loading the file.

## Parameters



- The parameter `--epsilon` is the fitting tolerance that specifies the maximal distance of an inlier to its supporting plane (default `0.4% * bounding box diagonal`).

- The parameter `--sigma` is the minimal primitive size that allows primitives with a too low number of inliers to be discarded (default `50`).

- The parameter `--nn` is the k of the k-nearest neighbor graph and also used to estimate the points' normals if they are not provided (default `20`).

- The parameter `--normal_deviation` is (i) the minimum cosine of the angle between an inlier and its supporting plane, which is used for region growing, and (ii) the minimum cosine of the angle between two primitives that can be merged (default `0.85`).

- The parameter `--s` and `--c` are the weight for simplicity and completeness, note that `0 < s + c < 3` (default `1` and `1`).

- The parameter `--norm` decide the fidelity metric, which can be `normal`, `L2` or `hybrid` (default `hybrid`). `L2` represents the Euclidean distance between the inlier points and the associated supporting planes, `normal` presents the deviation between inliers' normals and the associated supporting planes' normals and `hybrid` mode uses `L2` in the priority queue and `normal` in the transfer operator.

- The parameter `--max_steps` is the maximum iterations of our exploration mechanism (default `7`).

- The parameter `--constraint` is an option to make sure that the final configuration does not degrade the fidelity, simplicity and completeness of the initial configuration (default `False`). 

- The parameter `--vg` is an option to output the results in 'vg' form([Vertex Group](https://github.com/LiangliangNan/PolyFit/blob/main/ReadMe-data.md)) (default `False`). 

- The parameter `--alpha` is an option to output the primitives that are represented as alpha shapes (default `True`). 

- The parameter `--hull` is an option to output the primitives that are represented as convex hulls (default `False`). 

## Usage

- An example:
```
.\GoCoPP.exe path_to_your_file.ply --constraint --norm L2
```


## License

Inria



