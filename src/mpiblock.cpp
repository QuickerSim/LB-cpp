#include "mpiblock.hpp"

void MPIBlock::AddLBBlock(LBBlock* block) {
	_blocksmap.insert( std::make_pair(block->cLBBlockOwner()->cTreeID().cID(), block) );
}

void MPIBlock::FindNeighbors() {
  for (auto& b : _blocksmap) {
    FindNorth(*b.second);

    FindEast(*b.second);

    FindSouth(*b.second);

    FindWest(*b.second);

    FindNorthWest(*b.second);

    FindNorthEast(*b.second);

    FindSouthWest(*b.second);

    FindSouthEast(*b.second);
  }
}

void MPIBlock::CreateDepthMap() {
  for (auto& i : _blocksmap) {
    _depthmap[i.second->cLBBlockOwner()->cTreeID().cDepth()].push_back(i.first);
  }
  _mindepth = _depthmap.begin()->first;
  _maxdepth = _depthmap.rbegin()->first;
}

void MPIBlock::InitNodes(QSSize nx, QSSize ny, QSSize ng, QSDouble omega, QSDouble uini, QSDouble vini) {
  for (auto& b : _blocksmap) {
		b.second->Init(nx, ny, ng, omega, uini, vini);
  }

  FindNeighbors();
  CreateDepthMap();
}

void MPIBlock::Simulate(QSSize maxtimesteps, QSSize savetimestep) {
  QSStringStream filename;
  filename << "sym.pvd";

  QSFileStream file(filename.str(), std::ios::out);

  WriteHeader(file);
  WriteStep(0, file);
	QSVector<std::chrono::duration<QSDouble>> times;


  for (QSSize t = 0; t <= maxtimesteps; ++t) {
		auto start = std::chrono::system_clock::now();
    DoStepOnDepth(_mindepth);

		auto end = std::chrono::system_clock::now();
		times.push_back(end-start);

    if (t % savetimestep == 0) {
			std::cout << "elapsed time: " << times.rbegin()[0].count() << " s\n";
      WriteStep(t + 1, file);
    }
  }

  WriteFooter(file);
  file.close();
}

void MPIBlock::DoStepOnDepth(QSSize depth) {
  for (auto& block : _depthmap[depth]) {
		_blocksmap.at(block)->CollisionStep();
  }

  if (depth != _maxdepth) {
    DoStepOnDepth(depth + 1);
  }

  if (depth != _mindepth) {
    for (auto& block : _depthmap[depth]) {
      _blocksmap.at(block)->InterpolateCoarseToFine();
    }
  }

  for (auto& block : _depthmap[depth]) {
    _blocksmap.at(block)->SameLevelCommunicationStep();
  }

  for (auto& block : _depthmap[depth]) {
    _blocksmap.at(block)->StreamingStep();
  }

  if (depth != _maxdepth) {
    for (auto& block : _depthmap[depth]) {
      _blocksmap.at(block)->InterpolateFineToCoarse();
    }
  }

  if (depth == _mindepth) {
    return;
  }

  for (auto& block : _depthmap[depth]) {
		_blocksmap.at(block)->CollisionStep();
  }

  if (depth != _maxdepth) {
    DoStepOnDepth(depth + 1);
  }

  for (auto& block : _depthmap[depth]) {
    _blocksmap.at(block)->SameLevelCommunicationStep();
  }

  for (auto& block : _depthmap[depth]) {
    _blocksmap.at(block)->StreamingStep();
  }

  if (depth != _maxdepth) {
    for (auto& block : _depthmap[depth]) {
      _blocksmap.at(block)->InterpolateFineToCoarse();
    }
  }

}

void MPIBlock::WriteHeader(QSFileStream& file) {
  file << "<?xml version=\"1.0\"?>" << std::endl;
  file << "<VTKFile type=\"Collection\" version=\"0.1\" "
          "byte_order=\"LittleEndian\">"
       << std::endl;
  file << "<Collection>" << std::endl;
}

void MPIBlock::WriteFooter(QSFileStream& file) {
  file << "</Collection>" << std::endl;
  file << "</VTKFile>" << std::endl;
}

void MPIBlock::WriteStep(QSSize timestep, QSFileStream& file) {
  QSSize counter = 0;

  std::cout << "time: " << timestep << std::endl;

  for (auto& b : _blocksmap) {
    b.second->WriteToVTI(counter, timestep);
    file << "<DataSet timestep=\"" << timestep
         << "\" group=\"\" part=\"0\" file=\"t_" << timestep << "/sym_"
         << counter << ".vti\"/>" << std::endl;
    ++counter;
  }
}

void MPIBlock::FindNorth(LBBlock& b) {
  TreeID tidnorth = b.cLBBlockOwner()->cTreeID().moveN();

  if (tidnorth.cID() == 0) {
    return;
  }

  if (_blocksmap.find(tidnorth.cID()) != _blocksmap.end()) {
    b.rNorth().push_back(_blocksmap.at(tidnorth.cID()));
  } else {
    if (_blocksmap.find(tidnorth.Parent().cID()) != _blocksmap.end()) {
      b.rNorth().push_back(_blocksmap.at(tidnorth.Parent().cID()));
    } else {
      std::array<TreeID, 4> nchilds = tidnorth.Childs();
      if (_blocksmap.find(nchilds[0].cID()) != _blocksmap.end()) {
        b.rNorth().push_back(_blocksmap.at(nchilds[0].cID()));
      }
      if (_blocksmap.find(nchilds[1].cID()) != _blocksmap.end()) {
        b.rNorth().push_back(_blocksmap.at(nchilds[1].cID()));
      }
    }
  }
}

void MPIBlock::FindEast(LBBlock& b) {
  TreeID tideast = b.cLBBlockOwner()->cTreeID().moveE();

  if (tideast.cID() == 0) {
    return;
  }

  if (_blocksmap.find(tideast.cID()) != _blocksmap.end()) {
    b.rEast().push_back(_blocksmap.at(tideast.cID()));
  } else {
    if (_blocksmap.find(tideast.Parent().cID()) != _blocksmap.end()) {
      b.rEast().push_back(_blocksmap.at(tideast.Parent().cID()));
    } else {
      std::array<TreeID, 4> nchilds = tideast.Childs();
      if (_blocksmap.find(nchilds[0].cID()) != _blocksmap.end()) {
        b.rEast().push_back(_blocksmap.at(nchilds[0].cID()));
      }
      if (_blocksmap.find(nchilds[2].cID()) != _blocksmap.end()) {
        b.rEast().push_back(_blocksmap.at(nchilds[2].cID()));
      }
    }
  }
}

void MPIBlock::FindSouth(LBBlock& b) {
  TreeID tidsouth = b.cLBBlockOwner()->cTreeID().moveS();

  if (tidsouth.cID() == 0) {
    return;
  }

  if (_blocksmap.find(tidsouth.cID()) != _blocksmap.end()) {
    b.rSouth().push_back(_blocksmap.at(tidsouth.cID()));
  } else {
    if (_blocksmap.find(tidsouth.Parent().cID()) != _blocksmap.end()) {
      b.rSouth().push_back(_blocksmap.at(tidsouth.Parent().cID()));
    } else {
      std::array<TreeID, 4> nchilds = tidsouth.Childs();
      if (_blocksmap.find(nchilds[2].cID()) != _blocksmap.end()) {
        b.rSouth().push_back(_blocksmap.at(nchilds[2].cID()));
      }
      if (_blocksmap.find(nchilds[3].cID()) != _blocksmap.end()) {
        b.rSouth().push_back(_blocksmap.at(nchilds[3].cID()));
      }
    }
  }
}

void MPIBlock::FindWest(LBBlock& b) {
  TreeID tidwest = b.cLBBlockOwner()->cTreeID().moveW();

  if (tidwest.cID() == 0) {
    return;
  }

  if (_blocksmap.find(tidwest.cID()) != _blocksmap.end()) {
    b.rWest().push_back(_blocksmap.at(tidwest.cID()));
  } else {
    if (_blocksmap.find(tidwest.Parent().cID()) != _blocksmap.end()) {
      b.rWest().push_back(_blocksmap.at(tidwest.Parent().cID()));
    } else {
      std::array<TreeID, 4> nchilds = tidwest.Childs();
      if (_blocksmap.find(nchilds[1].cID()) != _blocksmap.end()) {
        b.rWest().push_back(_blocksmap.at(nchilds[1].cID()));
      }
      if (_blocksmap.find(nchilds[3].cID()) != _blocksmap.end()) {
        b.rWest().push_back(_blocksmap.at(nchilds[3].cID()));
      }
    }
  }
}

void MPIBlock::FindNorthWest(LBBlock& b) {
  TreeID tidnorthwest = b.cLBBlockOwner()->cTreeID().moveNW();

  if (tidnorthwest.cID() == 0) {
    return;
  }

  if (_blocksmap.find(tidnorthwest.cID()) != _blocksmap.end()) {
    b.rNorthWest().push_back(_blocksmap.at(tidnorthwest.cID()));
  } else {
    if (_blocksmap.find(tidnorthwest.Parent().cID()) != _blocksmap.end()) {
      b.rNorthWest().push_back(_blocksmap.at(tidnorthwest.Parent().cID()));
    } else {
      std::array<TreeID, 4> nchilds = tidnorthwest.Childs();
			if (_blocksmap.find(nchilds[1].cID()) != _blocksmap.end()) {
				b.rNorthWest().push_back(_blocksmap.at(nchilds[1].cID()));
      }
    }
  }
}

void MPIBlock::FindNorthEast(LBBlock& b) {
  TreeID tidnortheast = b.cLBBlockOwner()->cTreeID().moveNE();

  if (tidnortheast.cID() == 0) {
    return;
  }

  if (_blocksmap.find(tidnortheast.cID()) != _blocksmap.end()) {
    b.rNorthEast().push_back(_blocksmap.at(tidnortheast.cID()));
  } else {
    if (_blocksmap.find(tidnortheast.Parent().cID()) != _blocksmap.end()) {
      b.rNorthEast().push_back(_blocksmap.at(tidnortheast.Parent().cID()));
    } else {
      std::array<TreeID, 4> nchilds = tidnortheast.Childs();
			if (_blocksmap.find(nchilds[0].cID()) != _blocksmap.end()) {
				b.rNorthEast().push_back(_blocksmap.at(nchilds[0].cID()));
      }
    }
  }
}

void MPIBlock::FindSouthWest(LBBlock& b) {
  TreeID tidsouthwest = b.cLBBlockOwner()->cTreeID().moveSW();

  if (tidsouthwest.cID() == 0) {
    return;
  }

  if (_blocksmap.find(tidsouthwest.cID()) != _blocksmap.end()) {
    b.rSouthWest().push_back(_blocksmap.at(tidsouthwest.cID()));
  } else {
    if (_blocksmap.find(tidsouthwest.Parent().cID()) != _blocksmap.end()) {
      b.rSouthWest().push_back(_blocksmap.at(tidsouthwest.Parent().cID()));
    } else {
      std::array<TreeID, 4> nchilds = tidsouthwest.Childs();
      if (_blocksmap.find(nchilds[3].cID()) != _blocksmap.end()) {
        b.rSouthWest().push_back(_blocksmap.at(nchilds[3].cID()));
      }
    }
  }
}

void MPIBlock::FindSouthEast(LBBlock& b) {
  TreeID tidsoutheast = b.cLBBlockOwner()->cTreeID().moveSE();

  if (tidsoutheast.cID() == 0) {
    return;
  }

  if (_blocksmap.find(tidsoutheast.cID()) != _blocksmap.end()) {
    b.rSouthEast().push_back(_blocksmap.at(tidsoutheast.cID()));
  } else {
    if (_blocksmap.find(tidsoutheast.Parent().cID()) != _blocksmap.end()) {
      b.rSouthEast().push_back(_blocksmap.at(tidsoutheast.Parent().cID()));
    } else {
      std::array<TreeID, 4> nchilds = tidsoutheast.Childs();
      if (_blocksmap.find(nchilds[2].cID()) != _blocksmap.end()) {
        b.rSouthEast().push_back(_blocksmap.at(nchilds[2].cID()));
      }
    }
  }
}
