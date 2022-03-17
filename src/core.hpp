#ifndef CORE_HPP
#define CORE_HPP

#ifndef PROFILER
#define PROFILER 0
#endif

#if PROFILER
#include <easy/profiler.h>
#endif

#include <array>
#include <assert.h>
#include <bitset>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <unordered_map>
#include <vector>
#include <chrono>

typedef bool QSBool;
typedef char QSChar;
typedef int QSInt;
typedef unsigned int QSSize;
typedef long long QSLInt;
typedef unsigned long long int QSLSize;
typedef float QSFloat;
typedef double QSDouble;
typedef std::string QSString;
typedef std::stringstream QSStringStream;
typedef std::fstream QSFileStream;

template <class Key, class T>
using QSHashMap = std::unordered_map<Key, T>;

template <class Key, class T>
using QSMap = std::map<Key, T>;

template <class T>
using QSVector = std::vector<T>;

template <class T>
using QSArray = QSVector<QSVector<T>>;

/** @brief Class representing 2D Point.

		Detailed description.
		@author Konrad Jabłoński
		@date December 2018
		*/
struct Point {
	Point(QSDouble x_, QSDouble y_)
	{
		x = x_;
		y = y_;
	}
	QSDouble x;
	QSDouble y;
};

enum QTSide { ROOT,
	NW,
	NE,
	SW,
	SE };

enum LBNodeType {
	INTERIOR,
	WALLBB,
	NORTHVELOCITY,
	SOUTHVELOCITY,
	EASTVELOCITY,
	WESTVELOCITY,
	NORTHPRESSURE,
	SOUTHPRESSURE,
	EASTPRESSURE,
	WESTPRESSURE,
	GHOST
};

#endif // CORE_HPP
