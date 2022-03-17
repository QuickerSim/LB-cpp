#ifndef LBNODE_H
#define LBNODE_H

#include "core.hpp"

class QuadTree;

/** @brief Class representing Lattice Boltzmann Block.

		Detailed description.
		@author Konrad Jabłoński
		@date January 2019
		*/
class LBBlock {
public:
    LBBlock() {}
    LBBlock(QuadTree* lbblockowner) { _lbblockowner = lbblockowner; }
    ~LBBlock() {}

    // non-const getters
    QSSize& rNx() { return _nx; }
    QSSize& rNy() { return _ny; }

    QSVector<LBBlock*>& rNorth() { return _north; }
    QSVector<LBBlock*>& rEast() { return _east; }
    QSVector<LBBlock*>& rSouth() { return _south; }
    QSVector<LBBlock*>& rWest() { return _west; }
    QSVector<LBBlock*>& rNorthWest() { return _northwest; }
    QSVector<LBBlock*>& rNorthEast() { return _northeast; }
    QSVector<LBBlock*>& rSouthWest() { return _southwest; }
    QSVector<LBBlock*>& rSouthEast() { return _southeast; }

    QuadTree* rLBBlockOwner() { return _lbblockowner; }

    // const getters
    const QSSize& cNx() const { return _nx; }
    const QSSize& cNy() const { return _ny; }

    const QSVector<LBBlock*>& cNorth() const { return _north; }
    const QSVector<LBBlock*>& cEast() const { return _east; }
    const QSVector<LBBlock*>& cSouth() const { return _south; }
    const QSVector<LBBlock*>& cWest() const { return _west; }
    const QSVector<LBBlock*>& cNorthWest() const { return _northwest; }
    const QSVector<LBBlock*>& cNorthEast() const { return _northeast; }
    const QSVector<LBBlock*>& cSouthWest() const { return _southwest; }
    const QSVector<LBBlock*>& cSouthEast() const { return _southeast; }

    QuadTree* const& cLBBlockOwner() const { return _lbblockowner; }

    // methods
		// init methods
    void Init(QSSize, QSSize, QSSize, QSDouble, QSDouble, QSDouble);

		// 'step' methods
    void CollisionStep();
    void StreamingStep();
    void SameLevelCommunicationStep();
    void InterpolateCoarseToFine();
    void InterpolateFineToCoarse();

		// collision methods
		void Collision      (const QSSize i, const QSSize j);
		void EmptyCollision (const QSSize i, const QSSize j);

		// boundary conditions methods
		void NorthVelocity  (const QSSize i, const QSSize j, const QSDouble u0, const QSDouble v0);
		void SouthVelocity  (const QSSize i, const QSSize j, const QSDouble u0, const QSDouble v0);
		void EastVelocity   (const QSSize i, const QSSize j, const QSDouble u0, const QSDouble v0);
		void WestVelocity   (const QSSize i, const QSSize j, QSDouble u0, QSDouble v0);
		void NorthPressure  (const QSSize i, const QSSize j, const QSDouble p0);
		void SouthPressure  (const QSSize i, const QSSize j, const QSDouble p0);
		void EastPressure   (const QSSize i, const QSSize j, const QSDouble p0);
		void WestPressure   (const QSSize i, const QSSize j, const QSDouble p0);
		void BounceBack     (const QSSize i, const QSSize j);

		// streaming methods
		void Streaming      (const QSSize i, const QSSize j);
		void Streamingf0    (const QSSize i, const QSSize j);
		void Streamingf1    (const QSSize i, const QSSize j);
		void Streamingf2    (const QSSize i, const QSSize j);
		void Streamingf3    (const QSSize i, const QSSize j);
		void Streamingf4    (const QSSize i, const QSSize j);
		void Streamingf5    (const QSSize i, const QSSize j);
		void Streamingf6    (const QSSize i, const QSSize j);
		void Streamingf7    (const QSSize i, const QSSize j);
		void Streamingf8    (const QSSize i, const QSSize j);

    void GhostInfo();

		void ComputeDensity(const QSSize i, const QSSize j);

		void ComputeVelocityX(const QSSize i, const QSSize j);
		void ComputeVelocityY(const QSSize i, const QSSize j);

		QSDouble CoordinateX(const QSSize i, const QSSize j);
		QSDouble CoordinateY(const QSSize i, const QSSize j);

    void WriteToVTI(QSSize lp, QSSize time);

		// TESTING METHODS
		void InitZero(QSSize, QSSize, QSSize, QSDouble, QSDouble, QSDouble);
		void InitWave();
		void InitWaveZero();
		void EmptyCollisionStep();

private:
    QuadTree* _lbblockowner;

    QSSize _nx;
    QSSize _ny;
    QSSize _ng;

    QSVector<LBBlock*> _north;
    QSVector<LBBlock*> _east;
    QSVector<LBBlock*> _south;
    QSVector<LBBlock*> _west;
    QSVector<LBBlock*> _northwest;
    QSVector<LBBlock*> _northeast;
    QSVector<LBBlock*> _southwest;
    QSVector<LBBlock*> _southeast;

    QSDouble _omega;

    QSDouble _u_ini;
    QSDouble _v_ini;

    QSArray<QSVector<QSDouble>> _f;

    QSArray<QSVector<QSDouble>> _f_buf;

    QSArray<LBNodeType> _flag;

    QSArray<QSDouble> _u;
    QSArray<QSDouble> _v;
    QSArray<QSDouble> _rho;
};

#endif // LBNODE_H
