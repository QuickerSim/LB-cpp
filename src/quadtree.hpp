#ifndef QUADTREE_H
#define QUADTREE_H

#include "core.hpp"
#include "lbblock.hpp"
#include "treeid.hpp"

class MPIBlock;

/** @brief Class representing QuadTree.

		Detailed description.
		@author Konrad Jabłoński
		@date January 2019
		*/
class QuadTree {
public:
    ~QuadTree();
    QuadTree(QSDouble, QSDouble, QSDouble, QSDouble);

    TreeID& rTreeID() { return _treeid; }
    QSDouble& rXmin() { return _xmin; }
    QSDouble& rXmax() { return _xmax; }
    QSDouble& rYmin() { return _ymin; }
    QSDouble& rYmax() { return _ymax; }
    QTSide& rSide() { return _side; }
    QuadTree* rNw() { return _nw; }
    QuadTree* rNe() { return _ne; }
    QuadTree* rSw() { return _sw; }
    QuadTree* rSe() { return _se; }
    LBBlock* rNode() { return _node; }

    const TreeID& cTreeID() const { return _treeid; }
    const QSDouble& cXmin() const { return _xmin; }
    const QSDouble& cXmax() const { return _xmax; }
    const QSDouble& cYmin() const { return _ymin; }
    const QSDouble& cYmax() const { return _ymax; }
    const QTSide& cSide() const { return _side; }
    QuadTree* const& cNw() const { return _nw; }
    QuadTree* const& cNe() const { return _ne; }
    QuadTree* const& cSw() const { return _sw; }
    QuadTree* const& cSe() const { return _se; }
    LBBlock* const& cNode() const { return _node; }

    void CreateTreeAtPoints(QSVector<Point>, QSSize);
    void CreateTreeByPointCount(QSVector<Point>);
    void CreateTreeByPointCount(QSVector<Point>, QSSize);
    void BalanceTree(QSSize);
    void PopulateMpiBlock(MPIBlock& mpiblock);

    QuadTree* NorthNeighbor() const;
    QuadTree* SouthNeighbor() const;
    QuadTree* EastNeighbor() const;
    QuadTree* WestNeighbor() const;

    void InitNode(QuadTree* lbsgowner) { _node = new LBBlock(lbsgowner); }

private:
    QuadTree(TreeID, QSDouble, QSDouble, QSDouble, QSDouble, QuadTree*, QTSide);

    void CreateChildrens();
    QSBool RefineAtPoints(QSVector<Point>);
    QSBool RefineByPointCount(QSVector<Point>);
    QSBool RefineByPointCount(QSVector<Point>, QSSize);
    void AddLeafs(QSVector<QuadTree*>&);

    TreeID _treeid;

    QuadTree* _parent;

    QuadTree* _nw;
    QuadTree* _ne;
    QuadTree* _sw;
    QuadTree* _se;

    QTSide _side;

    QSDouble _xmin;
    QSDouble _xmax;
    QSDouble _ymin;
    QSDouble _ymax;

    LBBlock* _node;
};

#endif // QUADTREE_H
