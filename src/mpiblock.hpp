#ifndef MPIBLOCK_HPP
#define MPIBLOCK_HPP

#include "core.hpp"
#include "lbblock.hpp"
#include "quadtree.hpp"

/** @brief Class representing OpenMPI Block.

		Detailed description.
		@author Konrad Jabłoński
		@date January 2019
		*/
class MPIBlock {
public:
    MPIBlock() {}

    QSHashMap<QSSize, LBBlock*>& rBlocksVec() { return _blocksmap; }

    const QSHashMap<QSSize, LBBlock*>& cBlocksVec() const { return _blocksmap; }

    void AddLBBlock(LBBlock* block);

    void InitNodes(QSSize nx, QSSize ny, QSSize ng, QSDouble omega, QSDouble uini,
        QSDouble vini);
    void Simulate(QSSize maxtimesteps, QSSize savetimestep);

    void WriteStep(QSSize timestep, QSFileStream& file);
    void WriteHeader(QSFileStream& file);
    void WriteFooter(QSFileStream& file);

private:
    QSHashMap<QSSize, LBBlock*> _blocksmap;
    QSMap<QSSize, QSVector<QSSize>> _depthmap;

    QSSize _mindepth;
    QSSize _maxdepth;

    void FindNeighbors();
    void CreateDepthMap();
    void DoStepOnDepth(QSSize);

    void FindNorth(LBBlock& b);
    void FindEast(LBBlock& b);
    void FindSouth(LBBlock& b);
    void FindWest(LBBlock& b);
    void FindNorthWest(LBBlock& b);
    void FindNorthEast(LBBlock& b);
    void FindSouthWest(LBBlock& b);
    void FindSouthEast(LBBlock& b);
};

#endif // MPIBLOCK_HPP
