#include "quadtree.hpp"
#include "mpiblock.hpp"

QuadTree::QuadTree(TreeID treeid, QSDouble left, QSDouble right,
                   QSDouble bottom, QSDouble top, QuadTree *par, QTSide si) {
  _treeid = treeid;

  _xmin = left;
  _xmax = right;
  _ymin = bottom;
  _ymax = top;
  _parent = par;
  _side = si;
  _node = nullptr;

  _nw = nullptr;
  _ne = nullptr;
  _sw = nullptr;
  _se = nullptr;
}

QuadTree::QuadTree(QSDouble left, QSDouble right, QSDouble bottom,
                   QSDouble top) {
  _treeid.rID() = 0;
  _treeid.rDepth() = 0;

  _xmin = left;
  _xmax = right;
  _ymin = bottom;
  _ymax = top;
  _parent = nullptr;
  _side = QTSide::ROOT;
  _node = nullptr;

  _nw = nullptr;
  _ne = nullptr;
  _sw = nullptr;
  _se = nullptr;
}

QuadTree::~QuadTree() {
  delete _nw;
  delete _ne;
  delete _sw;
  delete _se;

  if (_node != nullptr) {
    delete _node;
  }
}

void QuadTree::CreateChildrens() {
  double yhalf = (_ymax - _ymin) / 2;
  double xhalf = (_xmax - _xmin) / 2;

  _sw = new QuadTree(TreeID(4 * _treeid.cID() + 1, _treeid.cDepth() + 1), _xmin,
                     _xmin + xhalf, _ymin, _ymin + yhalf, this, QTSide::SW);

  _se = new QuadTree(TreeID(4 * _treeid.cID() + 2, _treeid.cDepth() + 1),
                     _xmin + xhalf, _xmax, _ymin, _ymin + yhalf, this,
                     QTSide::SE);

  _nw = new QuadTree(TreeID(4 * _treeid.cID() + 3, _treeid.cDepth() + 1), _xmin,
                     _xmin + xhalf, _ymin + yhalf, _ymax, this, QTSide::NW);

  _ne = new QuadTree(TreeID(4 * _treeid.cID() + 4, _treeid.cDepth() + 1),
                     _xmin + xhalf, _xmax, _ymin + yhalf, _ymax, this,
                     QTSide::NE);
}

QSBool QuadTree::RefineAtPoints(QSVector<Point> points) {
  for (auto p : points) {
    if (p.x >= _xmin && p.x <= _xmax) {
      if (p.y >= _ymin && p.y <= _ymax) {
        return true;
      }
    }
  }

  return false;
}

QSBool QuadTree::RefineByPointCount(QSVector<Point> points) {
  QSSize numberofpoints = 0;

  for (auto p : points) {
    if (p.x >= _xmin && p.x <= _xmax) {
      if (p.y >= _ymin && p.y <= _ymax) {
        numberofpoints++;
      }
    }
  }

  if (numberofpoints > 1) {
    return true;
  } else {
    return false;
  }
}

QSBool QuadTree::RefineByPointCount(QSVector<Point> points, QSSize npoints) {
  QSSize numberofpoints = 0;

  for (auto p : points) {
    if (p.x >= _xmin && p.x <= _xmax) {
      if (p.y >= _ymin && p.y <= _ymax) {
        numberofpoints++;
      }
    }
  }

  if (numberofpoints > npoints) {
    return true;
  } else {
    return false;
  }
}

void QuadTree::CreateTreeAtPoints(QSVector<Point> points, QSSize max_depth) {
  if (_treeid.cDepth() >= max_depth) {
    return;
  }

  if (RefineAtPoints(points) == true) {
    CreateChildrens();

    _nw->CreateTreeAtPoints(points, max_depth);
    _ne->CreateTreeAtPoints(points, max_depth);
    _se->CreateTreeAtPoints(points, max_depth);
    _sw->CreateTreeAtPoints(points, max_depth);
  }
}

void QuadTree::CreateTreeByPointCount(QSVector<Point> points) {
  if (RefineByPointCount(points) == true) {
    CreateChildrens();

    _nw->CreateTreeByPointCount(points);
    _ne->CreateTreeByPointCount(points);
    _se->CreateTreeByPointCount(points);
    _sw->CreateTreeByPointCount(points);
  }
}

void QuadTree::CreateTreeByPointCount(QSVector<Point> points, QSSize npoints) {
  if (RefineByPointCount(points, npoints) == true) {
    CreateChildrens();

    _nw->CreateTreeByPointCount(points, npoints);
    _ne->CreateTreeByPointCount(points, npoints);
    _se->CreateTreeByPointCount(points, npoints);
    _sw->CreateTreeByPointCount(points, npoints);
  }
}

void QuadTree::AddLeafs(QSVector<QuadTree *> &leafs) {
  if (_nw != nullptr) {
    _nw->AddLeafs(leafs);
    _ne->AddLeafs(leafs);
    _sw->AddLeafs(leafs);
    _se->AddLeafs(leafs);
  } else
    leafs.push_back(this);
}

void QuadTree::PopulateMpiBlock(MPIBlock &mpiblock) {
  if (_nw != nullptr) {
    _nw->PopulateMpiBlock(mpiblock);
    _ne->PopulateMpiBlock(mpiblock);
    _sw->PopulateMpiBlock(mpiblock);
    _se->PopulateMpiBlock(mpiblock);
  } else {
    this->InitNode(this);
    mpiblock.AddLBBlock(this->cNode());
  }
}

void QuadTree::BalanceTree(QSSize depth_diff) {
  if (_side != QTSide::ROOT)
    throw;

  QSVector<QuadTree *> leafs;

  AddLeafs(leafs);

  while (!leafs.empty()) {
    if (leafs.back()->_side == QTSide::NE) { // NE
      if (leafs.back()->NorthNeighbor() != nullptr) {
        if (leafs.back()->cTreeID().cDepth() -
                leafs.back()->NorthNeighbor()->cTreeID().cDepth() >
            depth_diff) {
          QuadTree *tmp = leafs.back()->NorthNeighbor();

          tmp->CreateChildrens();

          leafs.pop_back();
          leafs.push_back(tmp->_nw);
          leafs.push_back(tmp->_ne);
          leafs.push_back(tmp->_sw);
          leafs.push_back(tmp->_se);
          continue;
        }
      }

      if (leafs.back()->EastNeighbor() != nullptr) {
        if (leafs.back()->cTreeID().cDepth() -
                leafs.back()->EastNeighbor()->cTreeID().cDepth() >
            depth_diff) {
          QuadTree *tmp = leafs.back()->EastNeighbor();

          tmp->CreateChildrens();

          leafs.pop_back();
          leafs.push_back(tmp->_nw);
          leafs.push_back(tmp->_ne);
          leafs.push_back(tmp->_sw);
          leafs.push_back(tmp->_se);
          continue;
        }
      }

      leafs.pop_back();

    } else if (leafs.back()->_side == QTSide::NW) { // NW
      if (leafs.back()->NorthNeighbor() != nullptr) {
        if (leafs.back()->cTreeID().cDepth() -
                leafs.back()->NorthNeighbor()->cTreeID().cDepth() >
            depth_diff) {
          QuadTree *tmp = leafs.back()->NorthNeighbor();

          tmp->CreateChildrens();

          leafs.pop_back();
          leafs.push_back(tmp->_nw);
          leafs.push_back(tmp->_ne);
          leafs.push_back(tmp->_sw);
          leafs.push_back(tmp->_se);
          continue;
        }
      }

      if (leafs.back()->WestNeighbor() != nullptr) {
        if (leafs.back()->cTreeID().cDepth() -
                leafs.back()->WestNeighbor()->cTreeID().cDepth() >
            depth_diff) {
          QuadTree *tmp = leafs.back()->WestNeighbor();

          tmp->CreateChildrens();

          leafs.pop_back();
          leafs.push_back(tmp->_nw);
          leafs.push_back(tmp->_ne);
          leafs.push_back(tmp->_sw);
          leafs.push_back(tmp->_se);
          continue;
        }
      }

      leafs.pop_back();

    } else if (leafs.back()->_side == QTSide::SW) { // SW
      if (leafs.back()->SouthNeighbor() != nullptr) {
        if (leafs.back()->cTreeID().cDepth() -
                leafs.back()->SouthNeighbor()->cTreeID().cDepth() >
            depth_diff) {
          QuadTree *tmp = leafs.back()->SouthNeighbor();

          tmp->CreateChildrens();

          leafs.pop_back();
          leafs.push_back(tmp->_nw);
          leafs.push_back(tmp->_ne);
          leafs.push_back(tmp->_sw);
          leafs.push_back(tmp->_se);
          continue;
        }
      }

      if (leafs.back()->WestNeighbor() != nullptr) {
        if (leafs.back()->cTreeID().cDepth() -
                leafs.back()->WestNeighbor()->cTreeID().cDepth() >
            depth_diff) {
          QuadTree *tmp = leafs.back()->WestNeighbor();

          tmp->CreateChildrens();

          leafs.pop_back();
          leafs.push_back(tmp->_nw);
          leafs.push_back(tmp->_ne);
          leafs.push_back(tmp->_sw);
          leafs.push_back(tmp->_se);
          continue;
        }
      }

      leafs.pop_back();

    } else if (leafs.back()->_side == QTSide::SE) { // SE
      if (leafs.back()->SouthNeighbor() != nullptr) {
        if (leafs.back()->cTreeID().cDepth() -
                leafs.back()->SouthNeighbor()->cTreeID().cDepth() >
            depth_diff) {
          QuadTree *tmp = leafs.back()->SouthNeighbor();

          tmp->CreateChildrens();

          leafs.pop_back();
          leafs.push_back(tmp->_nw);
          leafs.push_back(tmp->_ne);
          leafs.push_back(tmp->_sw);
          leafs.push_back(tmp->_se);
          continue;
        }
      }

      if (leafs.back()->EastNeighbor() != nullptr) {
        if (leafs.back()->cTreeID().cDepth() -
                leafs.back()->EastNeighbor()->cTreeID().cDepth() >
            depth_diff) {
          QuadTree *tmp = leafs.back()->EastNeighbor();

          tmp->CreateChildrens();

          leafs.pop_back();
          leafs.push_back(tmp->_nw);
          leafs.push_back(tmp->_ne);
          leafs.push_back(tmp->_sw);
          leafs.push_back(tmp->_se);
          continue;
        }
      }

      leafs.pop_back();
    } else if (leafs.back()->_side == QTSide::ROOT) {
      leafs.pop_back();
    }
  }
}

QuadTree *QuadTree::NorthNeighbor() const {
  if (_side == QTSide::ROOT) {
    return nullptr;
  } else {
    if (_side == 3)
      return _parent->_nw;
    else if (_side == 4)
      return _parent->_ne;
  }

  QuadTree *mu = nullptr;
  mu = _parent->NorthNeighbor();

  if (mu == nullptr) {
    return mu;
  } else {
    if (mu->_nw == nullptr)
      return mu;

    if (_side == 1)
      return mu->_sw;
    else if (_side == 2)
      return mu->_se;
  }

  std::cout << "North neighbor finding throw error" << std::endl;
  throw;
}

QuadTree *QuadTree::SouthNeighbor() const {
  if (_side == 0) {
    return nullptr;
  } else {
    if (_side == 1)
      return _parent->_sw;
    else if (_side == 2)
      return _parent->_se;
  }

  QuadTree *mu = nullptr;
  mu = _parent->SouthNeighbor();

  if (mu == nullptr) {
    return mu;
  } else {
    if (mu->_nw == nullptr)
      return mu;

    if (_side == 3)
      return mu->_nw;
    else if (_side == 4)
      return mu->_ne;
  }

  std::cout << "South neighbor finding throw error" << std::endl;
  throw;
}

QuadTree *QuadTree::EastNeighbor() const {
  if (_side == 0) {
    return nullptr;
  } else {
    if (_side == 1)
      return _parent->_ne;
    else if (_side == 3)
      return _parent->_se;
  }

  QuadTree *mu = nullptr;
  mu = _parent->EastNeighbor();

  if (mu == nullptr) {
    return mu;
  } else {
    if (mu->_nw == nullptr)
      return mu;

    if (_side == 2)
      return mu->_nw;
    else if (_side == 4)
      return mu->_sw;
  }

  std::cout << "East neighbor finding throw error" << std::endl;
  throw;
}

QuadTree *QuadTree::WestNeighbor() const {
  if (_side == 0) {
    return nullptr;
  } else {
    if (_side == 2)
      return _parent->_nw;
    else if (_side == 4)
      return _parent->_sw;
  }

  QuadTree *mu = nullptr;
  mu = _parent->WestNeighbor();

  if (mu == nullptr) {
    return mu;
  } else {
    if (mu->_nw == nullptr)
      return mu;

    if (_side == 1)
      return mu->_ne;
    else if (_side == 3)
      return mu->_se;
  }

  std::cout << "West neighbor finding throw error" << std::endl;
  throw;
}
