#pragma once

#ifndef Q_MOC_RUN
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>
#endif

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

// R-tree 2D
typedef bg::model::point<double, 2, bg::cs::cartesian> Boost_Point_2;
typedef bg::model::box<Boost_Point_2> Boost_Box_2;
typedef std::pair<Boost_Box_2, unsigned int> Boost_Value_2;
typedef bgi::rtree<Boost_Value_2, bgi::quadratic<16> > Boost_RTree_2;

// R-tree 3D
typedef bg::model::point<double, 3, bg::cs::cartesian> Boost_Point_3;
typedef bg::model::box<Boost_Point_3> Boost_Box_3;
typedef std::pair<Boost_Box_3, unsigned int> Boost_Value_3;
typedef bgi::rtree<Boost_Value_3, bgi::quadratic<16> > Boost_RTree_3;

template <typename T>
bool jin(const T & a, const T & x, const T & b) {
	return (a <= x && x <= b);
}

template <typename T>
bool jinr(const T & a, const T & x, const T & b) {
	return (a < x && x < b);
}

template <typename T>
const T & jmin(const T & a, const T & b) {
	return (a < b ? a : b);
}

template <typename T>
const T & jmax(const T & a, const T & b) {
	return (a > b ? a : b);
}

template <typename T>
const T & jclamp(const T & a, const T & x, const T & b) {
	return (x < a ? a : (x > b ? b : x));
}

#ifndef My_PI//(mulin)
#define My_PI 3.141592653589783238462643383279
#endif