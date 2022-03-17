#ifndef TREEID_HPP
#define TREEID_HPP

#include "core.hpp"

const QSSize qtelemsum[16] = {
    1, 5, 21, 85, 341, 1365, 5461, 21845,
    87381, 349525, 1398101, 5592405, 22369621, 89478485, 357913941, 1431655765
};
// 5726623061, 22906492245, 91625968981, 366503875925, 1466015503701 };

/** @brief Class representing implicit QuadTree numerating.

    Detailed description.
    @author Konrad Jabłoński
    @date December 2018
    */
class TreeID {

private:
    QSSize _id;
    QSSize _depth;

public:
    /** Default constructor */
    TreeID() {}

    /** Constructor for initializing with given id and depth
    @param id
    @param depth */
    TreeID(QSSize id, QSSize depth)
        : _id{ id }
        , _depth{ depth }
    {
    }

    /** Setter for ID */
    QSSize& rID() { return _id; }
    /** Setter for depth */
    QSSize& rDepth() { return _depth; }

    /** Getter for ID */
    const QSSize& cID() const { return _id; }
    /** Getter for ID */
    const QSSize& cDepth() const { return _depth; }

    /** Find TreeID for north neighbor */
    TreeID moveN() const;

    /** Find TreeID for east neighbor */
    TreeID moveE() const;

    /** Find TreeID for south neighbor */
    TreeID moveS() const;

    /** Find TreeID for west neighbor */
    TreeID moveW() const;

    /** Find TreeID for north-west neighbor */
    TreeID moveNW() const;

    /** Find TreeID for north-east neighbor */
    TreeID moveNE() const;

    /** Find TreeID for south-east neighbor */
    TreeID moveSE() const;

    /** Find TreeID for south-west neighbor */
    TreeID moveSW() const;

    /** Find TreeID of parent */
    TreeID Parent();

    /** Find TreeID's of TreeID children
    @return std::array of 4 TreeID children */
    std::array<TreeID, 4> Childs();
};

#endif // TREEID_HPP
