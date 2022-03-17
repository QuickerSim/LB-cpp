#include "mpiblock.hpp"
#include "quadtree.hpp"

int main()
{
    QuadTree root(0.0, 1.0, 0.0, 1.0);

    std::vector<Point> points;

//		points.push_back(Point(0.49, 0.51));
//		points.push_back(Point(0.51, 0.51));
//		points.push_back(Point(0.49, 0.49));
//		points.push_back(Point(0.51, 0.49));
//		points.push_back(Point(0, 1));
//		points.push_back(Point(0, 0.1));
//		points.push_back(Point(0, 0.2));
//		points.push_back(Point(0, 0.3));
//		points.push_back(Point(0, 0.4));
//		points.push_back(Point(0, 0.5));
//		points.push_back(Point(0, 0.6));
//		points.push_back(Point(0, 0.7));
//		points.push_back(Point(0, 0.8));
//		points.push_back(Point(0, 0.9));
//		points.push_back(Point(0, 1.0));

//		points.push_back(Point(1, 1));
//		points.push_back(Point(1, 0.1));
//		points.push_back(Point(1, 0.2));
//		points.push_back(Point(1, 0.3));
//		points.push_back(Point(1, 0.4));
//		points.push_back(Point(1, 0.5));
//		points.push_back(Point(1, 0.6));
//		points.push_back(Point(1, 0.7));
//		points.push_back(Point(1, 0.8));
//		points.push_back(Point(1, 0.9));
//		points.push_back(Point(1, 1.0));

		points.push_back(Point(0.35, 0.56));
		points.push_back(Point(0.25, 0.56));

		points.push_back(Point(0.25, 0.46));
		points.push_back(Point(0.35, 0.46));

//		points.push_back(Point(0.251, 1));
//		points.push_back(Point(0.249, 1));
//		points.push_back(Point(0.751, 1));
//		points.push_back(Point(0.749, 1));

		// root.CreateTreeByPointCount(points);
		root.CreateTreeAtPoints(points, 5);

    root.BalanceTree(1);

    MPIBlock dummyblock;

    root.PopulateMpiBlock(dummyblock);

		QSDouble omega = 1.0 / 0.51;

		dummyblock.InitNodes(32, 32, 2, omega, 0.1, 0.0);

		QSSize maxtimestep = 10000;
		QSSize sizeinterval = 1000;

    dummyblock.Simulate(maxtimestep, sizeinterval);

    return 0;
}
