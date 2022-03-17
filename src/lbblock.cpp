#include "lbblock.hpp"
#include "quadtree.hpp"

// ############################################################################################################# //

void LBBlock::Init(QSSize nx, QSSize ny, QSSize ng, QSDouble omega, QSDouble uini, QSDouble vini)
{
	QSDouble W1 = 4.0 / 9.0;
	QSDouble W2 = 1.0 / 9.0;
	QSDouble W3 = 1.0 / 36.0;

	std::vector<QSDouble> weigths;

	weigths.resize(9);

	weigths[0] = W1;
	weigths[1] = W2;
	weigths[2] = W2;
	weigths[3] = W2;
	weigths[4] = W2;
	weigths[5] = W3;
	weigths[6] = W3;
	weigths[7] = W3;
	weigths[8] = W3;

	_nx = nx;
	_ny = ny;
	_ng = ng;
	_omega = omega;

	_u_ini = uini;
	_v_ini = vini;

	_f.resize(9);

	for (QSSize i = 0; i < 9; ++i) {
		_f[i].resize(_nx + 2 * _ng);
		for (auto& j : _f[i]) {
			j.resize(_ny + 2 * _ng);
    }
	}

	_f_buf.resize(9);

	for (QSSize i = 0; i < 9; ++i) {
		_f_buf[i].resize(_nx + 2 * _ng);
		for (auto& j : _f_buf[i]) {
			j.resize(_ny + 2 * _ng);
    }
	}

	for (QSSize i = 0; i < _nx + 2 * _ng; ++i) {
		for (QSSize j = 0; j < _ny + 2 * _ng; ++j) {
			for (QSSize k = 0; k < 9; ++k) {
				_f[k][i][j] = weigths[k];
				_f_buf[k][i][j] = weigths[k];
			}
    }
	}

	_flag.resize(_nx + 2 * _ng);

	for (auto& i : _flag) {
		i.resize(_ny + 2 * _ng);
	}

//	for (QSSize j = 0; j < _ny + 2 * _ng; ++j) {
//		for (QSSize i = 0; i < _nx + 2 * _ng; ++i) {
//			if (i < _ng || i >= _nx + _ng || j < _ng || j >= _ny + _ng) {
//				_flag[i][j] = LBNodeType::GHOST;
//			} else if (fabs(CoordinateY(i, j)) < 1e-2 && j == _ng) { // bottom
//				_flag[i][j] = LBNodeType::WALLBB;
//			} else if (fabs(CoordinateY(i, j) - 1) < 1e-2 && j == _ny + _ng - 1) { // top
//				_flag[i][j] = LBNodeType::WALLBB;
//			} else if (fabs(CoordinateX(i, j)) < 1e-2 && i == _ng) { // left
//				_flag[i][j] = LBNodeType::WESTVELOCITY;
//			} else if (fabs(CoordinateX(i, j) - 1) < 1e-2 && i == _nx + _ng - 1) { // right
//				_flag[i][j] = LBNodeType::EASTPRESSURE;
//			} else
//				_flag[i][j] = INTERIOR;
//    }
//	}

	for (QSSize j = 0; j < _ny + 2 * _ng; ++j) {
		for (QSSize i = 0; i < _nx + 2 * _ng; ++i) {
			if (i < _ng || i >= _nx + _ng || j < _ng || j >= _ny + _ng) {
				_flag[i][j] = LBNodeType::GHOST;
			} else if ( (CoordinateX(i, j) - 0.3)*(CoordinateX(i, j) - 0.3) + (CoordinateY(i, j) - 0.51)*(CoordinateY(i, j) - 0.51) < 0.05*0.05) {
				_flag[i][j] = LBNodeType::WALLBB;
			} else if (fabs(CoordinateY(i, j)) < 1e-4) { // bottom
				_flag[i][j] = LBNodeType::WALLBB;
			} else if (fabs(CoordinateY(i, j) - 1) < 1e-4) { // top
				_flag[i][j] = LBNodeType::WALLBB;
			} else if (fabs(CoordinateX(i, j)) < 1e-4) { // left
				_flag[i][j] = LBNodeType::WESTVELOCITY;
			} else if (fabs(CoordinateX(i, j) - 1) < 1e-4) { // right
				_flag[i][j] = LBNodeType::EASTPRESSURE;
			} else
				_flag[i][j] = INTERIOR;
		}
	}

	_u.resize(_nx + 2 * _ng);
	for (auto& i : _u) {
		i.resize(_ny + 2 * _ng);
	}

	_v.resize(_nx + 2 * _ng);
	for (auto& i : _v) {
		i.resize(_ny + 2 * _ng);
	}

	_rho.resize(_nx + 2 * _ng);
	for (auto& i : _rho) {
		i.resize(_ny + 2 * _ng);
	}
}

// ############################################################################################################# //

void LBBlock::CollisionStep()
{
	#pragma omp parallel for
	for (QSSize i = 0; i < _nx + 2 * _ng; ++i) {
		for (QSSize j = 0; j < _ny + 2 * _ng; ++j) {
			switch (_flag[i][j]) {
			case LBNodeType::INTERIOR: // fluid
				ComputeDensity(i, j);
				ComputeVelocityX(i, j);
				ComputeVelocityY(i, j);
				Collision(i, j);
				// KForce(i, j);
				break;

			case LBNodeType::NORTHVELOCITY: // ZuoHe on north wall with given velocity
				NorthVelocity(i, j, _u_ini, _v_ini);
				ComputeDensity(i, j);
				ComputeVelocityX(i, j);
				ComputeVelocityY(i, j);
				Collision(i, j);
				break;

			case LBNodeType::SOUTHVELOCITY: // ZuoHe on south wall with given velocity
				SouthVelocity(i, j, _u_ini, _v_ini);
				ComputeDensity(i, j);
				ComputeVelocityX(i, j);
				ComputeVelocityY(i, j);
				Collision(i, j);
				break;

			case LBNodeType::EASTVELOCITY: // ZuoHe on east wall with given velocity
				EastVelocity(i, j, _u_ini, _v_ini);
				ComputeDensity(i, j);
				ComputeVelocityX(i, j);
				ComputeVelocityY(i, j);
				Collision(i, j);
				break;

			case LBNodeType::WESTVELOCITY: // ZuoHe on west wall with given velocity
				WestVelocity(i, j, _u_ini, _v_ini);
				ComputeDensity(i, j);
				ComputeVelocityX(i, j);
				ComputeVelocityY(i, j);
				Collision(i, j);
				break;

			case LBNodeType::NORTHPRESSURE: // ZuoHe on north wall with given pressure
				NorthPressure(i, j, 1);
				ComputeDensity(i, j);
				ComputeVelocityX(i, j);
				ComputeVelocityY(i, j);
				Collision(i, j);
				break;

			case LBNodeType::SOUTHPRESSURE: // ZuoHe on south wall with given pressure
				SouthPressure(i, j, 1);
				ComputeDensity(i, j);
				ComputeVelocityX(i, j);
				ComputeVelocityY(i, j);
				Collision(i, j);
				break;

			case LBNodeType::EASTPRESSURE: // ZuoHe on east wall with given pressure
				EastPressure(i, j, 1);
				ComputeDensity(i, j);
				ComputeVelocityX(i, j);
				ComputeVelocityY(i, j);
				Collision(i, j);
				break;

			case LBNodeType::WESTPRESSURE: // ZuoHe on west wall with given pressure
				WestPressure(i, j, 1);
				ComputeDensity(i, j);
				ComputeVelocityX(i, j);
				ComputeVelocityY(i, j);
				Collision(i, j);
				break;

			case LBNodeType::WALLBB: // BounceBack
				BounceBack(i, j);
				break;

			case LBNodeType::GHOST:
				EmptyCollision(i, j);
				break;

			default:
				throw;
			}
    }
	}
}

void LBBlock::StreamingStep()
{
	#pragma omp parallel for
	for (QSSize i = 1; i < _nx + 2 * _ng - 1; ++i) {
		for (QSSize j = 1; j < _ny + 2 * _ng - 1; ++j) {
			Streaming(i, j);
    }
	}

	for (QSSize j = 1; j < _ny + 2 * _ng - 1; ++j) {
		Streamingf0(0, j);
		Streamingf2(0, j);
		Streamingf3(0, j);
		Streamingf4(0, j);
		Streamingf6(0, j);
		Streamingf7(0, j);
	}

	for (QSSize j = 1; j < _ny + 2 * _ng - 1; ++j) {
		Streamingf0(_nx + 2 * _ng - 1, j);
		Streamingf1(_nx + 2 * _ng - 1, j);
		Streamingf2(_nx + 2 * _ng - 1, j);
		Streamingf4(_nx + 2 * _ng - 1, j);
		Streamingf5(_nx + 2 * _ng - 1, j);
		Streamingf8(_nx + 2 * _ng - 1, j);
	}

	for (QSSize i = 1; i < _nx + 2 * _ng - 1; ++i) {
		Streamingf0(i, 0);
		Streamingf1(i, 0);
		Streamingf3(i, 0);
		Streamingf4(i, 0);
		Streamingf7(i, 0);
		Streamingf8(i, 0);
	}

	for (QSSize i = 1; i < _nx + 2 * _ng - 1; ++i) {
		Streamingf0(i, _ny + 2 * _ng - 1);
		Streamingf1(i, _ny + 2 * _ng - 1);
		Streamingf2(i, _ny + 2 * _ng - 1);
		Streamingf3(i, _ny + 2 * _ng - 1);
		Streamingf5(i, _ny + 2 * _ng - 1);
		Streamingf6(i, _ny + 2 * _ng - 1);
	}

	Streamingf0(0, 0);
	Streamingf3(0, 0);
	Streamingf4(0, 0);
	Streamingf7(0, 0);

	Streamingf0(_nx + 2 * _ng - 1, 0);
	Streamingf1(_nx + 2 * _ng - 1, 0);
	Streamingf4(_nx + 2 * _ng - 1, 0);
	Streamingf8(_nx + 2 * _ng - 1, 0);

	Streamingf0(0, _ny + 2 * _ng - 1);
	Streamingf2(0, _ny + 2 * _ng - 1);
	Streamingf3(0, _ny + 2 * _ng - 1);
	Streamingf6(0, _ny + 2 * _ng - 1);

	Streamingf0(_nx + 2 * _ng - 1, _ny + 2 * _ng - 1);
	Streamingf1(_nx + 2 * _ng - 1, _ny + 2 * _ng - 1);
	Streamingf2(_nx + 2 * _ng - 1, _ny + 2 * _ng - 1);
	Streamingf5(_nx + 2 * _ng - 1, _ny + 2 * _ng - 1);
}

void LBBlock::SameLevelCommunicationStep()
{
	GhostInfo();
}

void LBBlock::InterpolateCoarseToFine()
{
	// NORTH
	if (_north.size() == 1) {
		if (_lbblockowner->cTreeID().cDepth() > _north[0]->_lbblockowner->cTreeID().cDepth()) {
			if (_lbblockowner->cSide() == QTSide::NW) {
				for (QSSize i = _ng; i < _nx + _ng; ++i) {
					for (QSSize j = _ny + _ng; j < _ny + 2 * _ng; ++j) {
						for (QSSize k = 0; k < 9; ++k) {
							_f_buf[k][i][j] = _north[0]->_f_buf[k][i / 2 + 1][(j - _ny) / 2 + 1];
						}
					}
				}
			} else if (_lbblockowner->cSide() == QTSide::NE) {
			for (QSSize k = 0; k < 9; ++k) {
				for (QSSize i = _ng; i < _nx + _ng; ++i) {
					for (QSSize j = _ny + _ng; j < _ny + 2 * _ng; ++j) {
							_f_buf[k][i][j] = _north[0]->_f_buf[k][(_nx + 2 * _ng) / 2 + i / 2 - 1][(j - _ny) / 2 + 1];
						}
					}
				}
			}
		}
	}

	// SOUTH
	if (_south.size() == 1) {
		if (_lbblockowner->cTreeID().cDepth() > _south[0]->_lbblockowner->cTreeID().cDepth()) {
			if (_lbblockowner->cSide() == QTSide::SW) {
				for (QSSize i = _ng; i < _nx + _ng; ++i) {
					for (QSSize j = 0; j < _ng; ++j) {
						for (QSSize k = 0; k < 9; ++k) {
							_f_buf[k][i][j] = _south[0]->_f_buf[k][i / 2 + 1][j / 2 + _ny + 1];
						}
					}
				}
			} else if (_lbblockowner->cSide() == QTSide::SE) {
				for (QSSize i = _ng; i < _nx + _ng; ++i) {
					for (QSSize j = 0; j < _ng; ++j) {
						for (QSSize k = 0; k < 9; ++k) {
							_f_buf[k][i][j] = _south[0]->_f_buf[k][(_nx + 2 * _ng) / 2 + i / 2 - 1][j / 2 + _ny + 1];
						}
					}
				}
			}
		}
	}

	// EAST
	if (_east.size() == 1) {
		if (_lbblockowner->cTreeID().cDepth() > _east[0]->_lbblockowner->cTreeID().cDepth()) {
			if (_lbblockowner->cSide() == QTSide::SE) {
				for (QSSize i = _nx + _ng; i < _nx + 2 * _ng; ++i) {
					for (QSSize j = _ng; j < _ny + _ng; ++j) {
						for (QSSize k = 0; k < 9; ++k) {
							_f_buf[k][i][j] = _east[0]->_f_buf[k][(i - _nx) / 2 + 1][j / 2 + 1];
						}
					}
				}
			} else if (_lbblockowner->cSide() == QTSide::NE) {
				for (QSSize i = _nx + _ng; i < _nx + 2 * _ng; ++i) {
					for (QSSize j = _ng; j < _ny + _ng; ++j) {
						for (QSSize k = 0; k < 9; ++k) {
							_f_buf[k][i][j] = _east[0]->_f_buf[k][(i - _nx) / 2 + 1][(_ny + 2 * _ng) / 2 + j / 2 - 1];
						}
					}
				}
			}
		}
	}

	// WEST
	if (_west.size() == 1) {
		if (_lbblockowner->cTreeID().cDepth() > _west[0]->_lbblockowner->cTreeID().cDepth()) {
			if (_lbblockowner->cSide() == QTSide::SW) {
				for (QSSize i = 0; i < _ng; ++i) {
					for (QSSize j = _ng; j < _ny + _ng; ++j) {
						for (QSSize k = 0; k < 9; ++k) {
							_f_buf[k][i][j] = _west[0]->_f_buf[k][i / 2 + _nx + 1][j / 2 + 1];
						}
					}
				}
			} else if (_lbblockowner->cSide() == QTSide::NW) {
				for (QSSize i = 0; i < _ng; ++i) {
					for (QSSize j = _ng; j < _ny + _ng; ++j) {
						for (QSSize k = 0; k < 9; ++k) {
							_f_buf[k][i][j] = _west[0]->_f_buf[k][i / 2 + _nx + 1][(_ny + 2 * _ng) / 2 + j / 2 - 1];
						}
					}
				}
			}
		}
	}

	// NORTHWEST
	if (_northwest.size() == 1) {
		if (_lbblockowner->cTreeID().cDepth() > _northwest[0]->_lbblockowner->cTreeID().cDepth()) {
			if (_lbblockowner->cSide() == QTSide::NW) {
				for (QSSize i = 0; i < _ng; ++i) {
					for (QSSize j = _ny + _ng; j < _ny + 2 * _ng; ++j) {
						for (QSSize k = 0; k < 9; ++k) {
							_f_buf[k][i][j] = _northwest[0]->_f_buf[k][i / 2 + _nx + 1][(j - _ny) / 2 + 1];
						}
					}
				}
			} else if (_lbblockowner->cSide() == QTSide::NE) {
				for (QSSize i = 0; i < _ng; ++i) {
					for (QSSize j = _ny + _ng; j < _ny + 2 * _ng; ++j) {
						for (QSSize k = 0; k < 9; ++k) {
							_f_buf[k][i][j] = _northwest[0]->_f_buf[k][i / 2 + _nx / 2 + 1][(j - _ny) / 2 + 1];
						}
					}
				}
			} else if (_lbblockowner->cSide() == QTSide::SW) {
				for (QSSize i = 0; i < _ng; ++i) {
					for (QSSize j = _ny + _ng; j < _ny + 2 * _ng; ++j) {
						for (QSSize k = 0; k < 9; ++k) {
							_f_buf[k][i][j] = _northwest[0]->_f_buf[k][i / 2 + _nx + 1][j / 2];
						}
					}
				}
			}
		}
	}

	// NORTHEAST
	if (_northeast.size() == 1) {
		if (_lbblockowner->cTreeID().cDepth() >	_northeast[0]->_lbblockowner->cTreeID().cDepth()) {
			if (_lbblockowner->cSide() == QTSide::NE) {
				for (QSSize i = _nx + _ng; i < _nx + 2 * _ng; ++i) {
					for (QSSize j = _ny + _ng; j < _ny + 2 * _ng; ++j) {
						for (QSSize k = 0; k < 9; ++k) {
							_f_buf[k][i][j] = _northeast[0]->_f_buf[k][(i - _nx) / 2 + 1][(j - _ny) / 2 + 1];
						}
					}
				}
			} else if (_lbblockowner->cSide() == QTSide::NW) {
				for (QSSize i = _nx + _ng; i < _nx + 2 * _ng; ++i) {
					for (QSSize j = _ny + _ng; j < _ny + 2 * _ng; ++j) {
						for (QSSize k = 0; k < 9; ++k) {
							_f_buf[k][i][j] = _northeast[0]->_f_buf[k][i / 2][(j - _ny) / 2 + 1];
						}
					}
				}
			} else if (_lbblockowner->cSide() == QTSide::SE) {
				for (QSSize i = _nx + _ng; i < _nx + 2 * _ng; ++i) {
					for (QSSize j = _ny + _ng; j < _ny + 2 * _ng; ++j) {
						for (QSSize k = 0; k < 9; ++k) {
							_f_buf[k][i][j] = _northeast[0]->_f_buf[k][(i - _nx) / 2 + 1][j / 2];
						}
					}
				}
			}
		}
	}

	// SOUTHWEST
	if (_southwest.size() == 1) {
		if (_lbblockowner->cTreeID().cDepth() > _southwest[0]->_lbblockowner->cTreeID().cDepth()) {
			if (_lbblockowner->cSide() == QTSide::SW) {
				for (QSSize i = 0; i < _ng; ++i) {
					for (QSSize j = 0; j < _ng; ++j) {
						for (QSSize k = 0; k < 9; ++k) {
							_f_buf[k][i][j] = _southwest[0]->_f_buf[k][i / 2 + _nx + 1][j / 2 + _ny + 1];
						}
					}
				}
			} else if (_lbblockowner->cSide() == QTSide::NW) {
				for (QSSize i = 0; i < _ng; ++i) {
					for (QSSize j = 0; j < _ng; ++j) {
						for (QSSize k = 0; k < 9; ++k) {
							_f_buf[k][i][j] = _southwest[0]->_f_buf[k][i / 2 + _nx + 1][j / 2 + _ny / 2 + 1];
						}
					}
				}
			} else if (_lbblockowner->cSide() == QTSide::SE) {
				for (QSSize i = 0; i < _ng; ++i) {
					for (QSSize j = 0; j < _ng; ++j) {
						for (QSSize k = 0; k < 9; ++k) {
							_f_buf[k][i][j] = _southwest[0]->_f_buf[k][i / 2 + _nx / 2 + 1][j / 2 + _ny + 1];
						}
					}
				}
			}
		}
	}

	// SOUTHEAST
	if (_southeast.size() == 1) {
		if (_lbblockowner->cTreeID().cDepth() > _southeast[0]->_lbblockowner->cTreeID().cDepth()) {
			if (_lbblockowner->cSide() == QTSide::SE) {
				for (QSSize i = _nx + _ng; i < _nx + 2 * _ng; ++i) {
					for (QSSize j = 0; j < _ng; ++j) {
						for (QSSize k = 0; k < 9; ++k) {
							_f_buf[k][i][j] = _southeast[0]->_f_buf[k][(i - _nx) / 2 + 1][j / 2 + _ny + 1];
						}
					}
				}
			} else if (_lbblockowner->cSide() == QTSide::NE) {
				for (QSSize i = _nx + _ng; i < _nx + 2 * _ng; ++i) {
					for (QSSize j = 0; j < _ng; ++j) {
						for (QSSize k = 0; k < 9; ++k) {
							_f_buf[k][i][j] = _southeast[0]->_f_buf[k][(i - _nx) / 2 + 1][j / 2 + _ny / 2 + 1];
						}
					}
				}
			} else if (_lbblockowner->cSide() == QTSide::SW) {
				for (QSSize i = _nx + _ng; i < _nx + 2 * _ng; ++i) {
					for (QSSize j = 0; j < _ng; ++j) {
						for (QSSize k = 0; k < 9; ++k) {
							_f_buf[k][i][j] = _southeast[0]->_f_buf[k][i / 2][j / 2 + _ny + 1];
						}
					}
				}
			}
		}
	}
}

void LBBlock::InterpolateFineToCoarse()
{
	// NORTH
	if (_north.size() == 2) {
		if (_lbblockowner->cTreeID().cDepth() < _north[0]->_lbblockowner->cTreeID().cDepth()) {
			for (QSSize i = _ng; i < _nx / 2 + _ng; ++i) {
				for (QSSize j = _ny + _ng - _ng / 2; j < _ny + _ng; ++j) {
					_f[4][i][j] = (_north[0]->_f[4][2 * i - _ng][2 * (j - _ny -  _ng / 2)		 ] + _north[0]->_f[4][2 * i - _ng + 1][2 * (j - _ny -  _ng / 2)		 ] +
												 _north[0]->_f[4][2 * i - _ng][2 * (j - _ny -  _ng / 2) + 1] + _north[0]->_f[4][2 * i - _ng + 1][2 * (j - _ny -  _ng / 2) + 1]) / 4;
					_f[7][i][j] = (_north[0]->_f[7][2 * i - _ng][2 * (j - _ny -  _ng / 2)		 ] + _north[0]->_f[7][2 * i - _ng + 1][2 * (j - _ny -  _ng / 2)		 ] +
												 _north[0]->_f[7][2 * i - _ng][2 * (j - _ny -  _ng / 2) + 1] + _north[0]->_f[7][2 * i - _ng + 1][2 * (j - _ny -  _ng / 2) + 1]) / 4;
					_f[8][i][j] = (_north[0]->_f[8][2 * i - _ng][2 * (j - _ny -  _ng / 2)		 ] + _north[0]->_f[8][2 * i - _ng + 1][2 * (j - _ny -  _ng / 2)		 ] +
												 _north[0]->_f[8][2 * i - _ng][2 * (j - _ny -  _ng / 2) + 1] + _north[0]->_f[8][2 * i - _ng + 1][2 * (j - _ny -  _ng / 2) + 1]) / 4;
				}
			}

			for (QSSize i = _nx / 2 + _ng; i < _nx + _ng; ++i) {
				for (QSSize j = _ny + _ng - _ng / 2; j < _ny + _ng; ++j) {
					_f[4][i][j] = (_north[1]->_f[4][2 * (i - (_nx + _ng) / 2)][2 * (j - _ny -  _ng / 2)		 ] + _north[1]->_f[4][2 * (i - (_nx + _ng) / 2) + 1][2 * (j - _ny -  _ng / 2)    ]
											 + _north[1]->_f[4][2 * (i - (_nx + _ng) / 2)][2 * (j - _ny -  _ng / 2) + 1] + _north[1]->_f[4][2 * (i - (_nx + _ng) / 2) + 1][2 * (j - _ny -  _ng / 2) + 1]) / 4;
					_f[7][i][j] = (_north[1]->_f[7][2 * (i - (_nx + _ng) / 2)][2 * (j - _ny -  _ng / 2)		 ] + _north[1]->_f[7][2 * (i - (_nx + _ng) / 2) + 1][2 * (j - _ny -  _ng / 2)    ]
											 + _north[1]->_f[7][2 * (i - (_nx + _ng) / 2)][2 * (j - _ny -  _ng / 2) + 1] + _north[1]->_f[7][2 * (i - (_nx + _ng) / 2) + 1][2 * (j - _ny -  _ng / 2) + 1]) / 4;
					_f[8][i][j] = (_north[1]->_f[8][2 * (i - (_nx + _ng) / 2)][2 * (j - _ny -  _ng / 2)		 ] + _north[1]->_f[8][2 * (i - (_nx + _ng) / 2) + 1][2 * (j - _ny -  _ng / 2)    ]
											 + _north[1]->_f[8][2 * (i - (_nx + _ng) / 2)][2 * (j - _ny -  _ng / 2) + 1] + _north[1]->_f[8][2 * (i - (_nx + _ng) / 2) + 1][2 * (j - _ny -  _ng / 2) + 1]) / 4;
				}
			}
		}
	}

	// SOUTH
	if (_south.size() == 2) {
		if (_lbblockowner->cTreeID().cDepth() < _south[0]->_lbblockowner->cTreeID().cDepth()) {
			for (QSSize i = _ng; i < (_nx + 2 * _ng) / 2; ++i) {
				for (QSSize j = _ng; j < _ng + _ng / 2; ++j) {
					_f[2][i][j] = (_south[0]->_f[2][2 * i - _ng][2 * j + _ny - _ng		] + _south[0]->_f[2][2 * i - _ng + 1][2 * j + _ny - _ng    ]
											 + _south[0]->_f[2][2 * i - _ng][2 * j + _ny - _ng + 1] + _south[0]->_f[2][2 * i - _ng + 1][2 * j + _ny - _ng + 1]) / 4;
					_f[5][i][j] = (_south[0]->_f[5][2 * i - _ng][2 * j + _ny - _ng		] + _south[0]->_f[5][2 * i - _ng + 1][2 * j + _ny - _ng    ]
											 + _south[0]->_f[5][2 * i - _ng][2 * j + _ny - _ng + 1] + _south[0]->_f[5][2 * i - _ng + 1][2 * j + _ny - _ng + 1]) / 4;
					_f[6][i][j] = (_south[0]->_f[6][2 * i - _ng][2 * j + _ny - _ng		] + _south[0]->_f[6][2 * i - _ng + 1][2 * j + _ny - _ng    ]
											 + _south[0]->_f[6][2 * i - _ng][2 * j + _ny - _ng + 1] + _south[0]->_f[6][2 * i - _ng + 1][2 * j + _ny - _ng + 1]) / 4;
				}
			}

			for (QSSize i = (_nx + 2 * _ng) / 2; i < _nx + _ng; ++i) {
				for (QSSize j = _ng; j < _ng + _ng / 2; ++j) {
					_f[2][i][j] = (_south[1]->_f[2][2 * (i - (_nx + _ng) / 2)][2 * j + _ny - _ng    ] + _south[1]->_f[2][2 * (i - (_nx + _ng) / 2) + 1][2 * j + _ny - _ng		 ]
											 + _south[1]->_f[2][2 * (i - (_nx + _ng) / 2)][2 * j + _ny - _ng + 1] + _south[1]->_f[2][2 * (i - (_nx + _ng) / 2) + 1][2 * j + _ny - _ng + 1]) / 4;
					_f[5][i][j] = (_south[1]->_f[5][2 * (i - (_nx + _ng) / 2)][2 * j + _ny - _ng		] + _south[1]->_f[5][2 * (i - (_nx + _ng) / 2) + 1][2 * j + _ny - _ng		 ]
											 + _south[1]->_f[5][2 * (i - (_nx + _ng) / 2)][2 * j + _ny - _ng + 1] + _south[1]->_f[5][2 * (i - (_nx + _ng) / 2) + 1][2 * j + _ny - _ng + 1]) / 4;
					_f[6][i][j] = (_south[1]->_f[6][2 * (i - (_nx + _ng) / 2)][2 * j + _ny - _ng		] + _south[1]->_f[6][2 * (i - (_nx + _ng) / 2) + 1][2 * j + _ny - _ng		 ]
											 + _south[1]->_f[6][2 * (i - (_nx + _ng) / 2)][2 * j + _ny - _ng + 1] + _south[1]->_f[6][2 * (i - (_nx + _ng) / 2) + 1][2 * j + _ny - _ng + 1]) / 4;
				}
			}
		}
	}

	// EAST
	if (_east.size() == 2) {
		if (_lbblockowner->cTreeID().cDepth() < _east[0]->_lbblockowner->cTreeID().cDepth()) {
			for (QSSize i = _nx + _ng - _ng / 2; i < _nx + _ng; ++i) {
				for (QSSize j = _ng; j < _ny / 2 + _ng; ++j) {
					_f[3][i][j] = (_east[0]->_f[3][2 * (i - _nx -  _ng / 2)    ][2 * j - _ng] + _east[0]->_f[3][2 * (i - _nx -  _ng / 2)    ][2 * j - _ng + 1] +
												 _east[0]->_f[3][2 * (i - _nx -  _ng / 2) + 1][2 * j - _ng] + _east[0]->_f[3][2 * (i - _nx -  _ng / 2) + 1][2 * j - _ng + 1]) / 4;
					_f[6][i][j] = (_east[0]->_f[6][2 * (i - _nx -  _ng / 2)    ][2 * j - _ng] + _east[0]->_f[6][2 * (i - _nx -  _ng / 2)    ][2 * j - _ng + 1] +
												 _east[0]->_f[6][2 * (i - _nx -  _ng / 2) + 1][2 * j - _ng] + _east[0]->_f[6][2 * (i - _nx -  _ng / 2) + 1][2 * j - _ng + 1]) / 4;
					_f[7][i][j] = (_east[0]->_f[7][2 * (i - _nx -  _ng / 2)    ][2 * j - _ng] + _east[0]->_f[7][2 * (i - _nx -  _ng / 2)    ][2 * j - _ng + 1] +
												 _east[0]->_f[7][2 * (i - _nx -  _ng / 2) + 1][2 * j - _ng] + _east[0]->_f[7][2 * (i - _nx -  _ng / 2) + 1][2 * j - _ng + 1]) / 4;
				}
			}
			for (QSSize i = _nx + _ng - _ng / 2; i < _nx + _ng; ++i) {
				for (QSSize j = _ny / 2 + _ng; j < _ny + _ng; ++j) {
					_f[3][i][j] = (_east[1]->_f[3][2 * (i - _nx -  _ng / 2)    ][2 * (j - (_ny + _ng) / 2)] + _east[1]->_f[3][2 * (i - _nx -  _ng / 2)    ][2 * (j - (_ny + _ng) / 2) + 1] +
												 _east[1]->_f[3][2 * (i - _nx -  _ng / 2) + 1][2 * (j - (_ny + _ng) / 2)] + _east[1]->_f[3][2 * (i - _nx -  _ng / 2) + 1][2 * (j - (_ny + _ng) / 2) + 1]) / 4;
					_f[6][i][j] = (_east[1]->_f[6][2 * (i - _nx -  _ng / 2)    ][2 * (j - (_ny + _ng) / 2)] + _east[1]->_f[6][2 * (i - _nx -  _ng / 2)    ][2 * (j - (_ny + _ng) / 2) + 1] +
												 _east[1]->_f[6][2 * (i - _nx -  _ng / 2) + 1][2 * (j - (_ny + _ng) / 2)] + _east[1]->_f[6][2 * (i - _nx -  _ng / 2) + 1][2 * (j - (_ny + _ng) / 2) + 1]) / 4;
					_f[7][i][j] = (_east[1]->_f[7][2 * (i - _nx -  _ng / 2)    ][2 * (j - (_ny + _ng) / 2)] + _east[1]->_f[7][2 * (i - _nx -  _ng / 2)    ][2 * (j - (_ny + _ng) / 2) + 1] +
												 _east[1]->_f[7][2 * (i - _nx -  _ng / 2) + 1][2 * (j - (_ny + _ng) / 2)] + _east[1]->_f[7][2 * (i - _nx -  _ng / 2) + 1][2 * (j - (_ny + _ng) / 2) + 1]) / 4;
				}
			}
		}
	}

	// WEST
	if (_west.size() == 2) {
		if (_lbblockowner->cTreeID().cDepth() < _west[0]->_lbblockowner->cTreeID().cDepth()) {
			for (QSSize i = _ng; i < _ng + _ng / 2; ++i) {
				for (QSSize j = _ng; j < (_ny + 2 * _ng) / 2; ++j) {
					_f[1][i][j] = (_west[0]->_f[1][2 * i + _nx - _ng    ][2 * j - _ng] + _west[0]->_f[1][2 * i + _nx - _ng    ][2 * j - _ng + 1]
											 + _west[0]->_f[1][2 * i + _nx - _ng + 1][2 * j - _ng] + _west[0]->_f[1][2 * i + _nx - _ng + 1][2 * j - _ng + 1]) / 4;
					_f[5][i][j] = (_west[0]->_f[5][2 * i + _nx - _ng    ][2 * j - _ng] + _west[0]->_f[5][2 * i + _nx - _ng    ][2 * j - _ng + 1]
											 + _west[0]->_f[5][2 * i + _nx - _ng + 1][2 * j - _ng] + _west[0]->_f[5][2 * i + _nx - _ng + 1][2 * j - _ng + 1]) / 4;
					_f[8][i][j] = (_west[0]->_f[8][2 * i + _nx - _ng    ][2 * j - _ng] + _west[0]->_f[8][2 * i + _nx - _ng    ][2 * j - _ng + 1]
											 + _west[0]->_f[8][2 * i + _nx - _ng + 1][2 * j - _ng] + _west[0]->_f[8][2 * i + _nx - _ng + 1][2 * j - _ng + 1]) / 4;
				}
			}
			for (QSSize i = _ng; i < _ng + _ng / 2; ++i) {
				for (QSSize j = (_ny + 2 * _ng) / 2; j < _ny + _ng; ++j) {
					_f[1][i][j] = (_west[1]->_f[1][2 * i + _nx - _ng    ][2 * (j - (_ny + _ng) / 2)] + _west[1]->_f[1][2 * i + _nx - _ng    ][2 * (j - (_ny + _ng) / 2) + 1]
											 + _west[1]->_f[1][2 * i + _nx - _ng + 1][2 * (j - (_ny + _ng) / 2)] + _west[1]->_f[1][2 * i + _nx - _ng + 1][2 * (j - (_ny + _ng) / 2) + 1]) / 4;
					_f[5][i][j] = (_west[1]->_f[5][2 * i + _nx - _ng    ][2 * (j - (_ny + _ng) / 2)] + _west[1]->_f[5][2 * i + _nx - _ng    ][2 * (j - (_ny + _ng) / 2) + 1]
											 + _west[1]->_f[5][2 * i + _nx - _ng + 1][2 * (j - (_ny + _ng) / 2)] + _west[1]->_f[5][2 * i + _nx - _ng + 1][2 * (j - (_ny + _ng) / 2) + 1]) / 4;
					_f[8][i][j] = (_west[1]->_f[8][2 * i + _nx - _ng    ][2 * (j - (_ny + _ng) / 2)] + _west[1]->_f[8][2 * i + _nx - _ng    ][2 * (j - (_ny + _ng) / 2) + 1]
											 + _west[1]->_f[8][2 * i + _nx - _ng + 1][2 * (j - (_ny + _ng) / 2)] + _west[1]->_f[8][2 * i + _nx - _ng + 1][2 * (j - (_ny + _ng) / 2) + 1]) / 4;
				}
			}
		}
	}

	// NORTHWEST
	if (_northwest.size() == 1) {
		if (_lbblockowner->cTreeID().cDepth() < _northwest[0]->_lbblockowner->cTreeID().cDepth()) {
			for (QSSize i = _ng; i < _ng + _ng / 2; ++i) {
				for (QSSize j = _ny + _ng - _ng / 2; j < _ny + _ng; ++j) {
					_f[8][i][j] = (_northwest[0]->_f[8][2 * i + _nx - _ng    ][2 * (j - _ny -  _ng / 2)] + _northwest[0]->_f[8][2 * i + _nx - _ng	   ][2 * (j - _ny -  _ng / 2) + 1]
											 + _northwest[0]->_f[8][2 * i + _nx - _ng + 1][2 * (j - _ny -  _ng / 2)] + _northwest[0]->_f[8][2 * i + _nx - _ng + 1][2 * (j - _ny -  _ng / 2) + 1] ) / 4;
				}
			}
		}
	}

	// NORTHEAST
	if (_northeast.size() == 1) {
		if (_lbblockowner->cTreeID().cDepth() < _northeast[0]->_lbblockowner->cTreeID().cDepth()) {
			for (QSSize i = _nx + _ng - _ng / 2; i < _nx + _ng; ++i) {
				for (QSSize j = _ny + _ng - _ng / 2; j < _ny + _ng; ++j) {
					_f[7][i][j] = (_northeast[0]->_f[7][i - _nx    ][j - _ny] + _northeast[0]->_f[7][i - _nx    ][j - _ny - 1]
											 + _northeast[0]->_f[7][i - _nx - 1][j - _ny] + _northeast[0]->_f[7][i - _nx - 1][j - _ny - 1] ) / 4;
				}
			}
		}
	}

	// SOUTHWEST
	if (_southwest.size() == 1) {
		if (_lbblockowner->cTreeID().cDepth() < _southwest[0]->_lbblockowner->cTreeID().cDepth()) {
			for (QSSize i = _ng; i < _ng + _ng / 2; ++i) {
				for (QSSize j = _ng; j < _ng + _ng / 2; ++j) {
					_f[5][i][j] = (_southwest[0]->_f[5][2 * i + _nx - _ng    ][2 * j + _ny - _ng] + _southwest[0]->_f[5][2 * i + _nx - _ng    ][2 * j + _ny - _ng + 1]
											 + _southwest[0]->_f[5][2 * i + _nx - _ng + 1][2 * j + _ny - _ng] + _southwest[0]->_f[5][2 * i + _nx - _ng + 1][2 * j + _ny - _ng + 1] ) / 4;
				}
			}
		}
	}

	// SOUTHEAST
	if (_southeast.size() == 1) {
		if (_lbblockowner->cTreeID().cDepth() < _southeast[0]->_lbblockowner->cTreeID().cDepth()) {
			for (QSSize i = _nx + _ng - _ng / 2; i < _nx + _ng; ++i) {
				for (QSSize j = _ng; j < _ng + _ng / 2; ++j) {
					_f[6][i][j] = (_southeast[0]->_f[6][i - _nx    ][j + _ny] + _southeast[0]->_f[6][i - _nx    ][j + _ny + 1]
											 + _southeast[0]->_f[6][i - _nx - 1][j + _ny] + _southeast[0]->_f[6][i - _nx - 1][j + _ny + 1]) / 4;
				}
			}
		}
	}
}

// ############################################################################################################# //

void LBBlock::Collision(const QSSize i, const QSSize j)
{
	QSDouble W1 = 4.0 / 9.0;
	QSDouble W2 = 1.0 / 9.0;
	QSDouble W3 = 1.0 / 36.0;

	double f_eq0;
	double f_eq1;
	double f_eq2;
	double f_eq3;
	double f_eq4;
	double f_eq5;
	double f_eq6;
	double f_eq7;
	double f_eq8;

	f_eq0 = W1 * _rho[i][j] * (1.0 - 1.5 * (_u[i][j] * _u[i][j] + _v[i][j] * _v[i][j]));
	f_eq1 = W2 * _rho[i][j] * (1.0 + 3.0 * _u[i][j] + 4.5 * _u[i][j] * _u[i][j] - 1.5 * (_u[i][j] * _u[i][j] + _v[i][j] * _v[i][j]));
	f_eq2 = W2 * _rho[i][j] * (1.0 + 3.0 * _v[i][j] + 4.5 * _v[i][j] * _v[i][j] - 1.5 * (_u[i][j] * _u[i][j] + _v[i][j] * _v[i][j]));
	f_eq3 = W2 * _rho[i][j] * (1.0 - 3.0 * _u[i][j] + 4.5 * _u[i][j] * _u[i][j] - 1.5 * (_u[i][j] * _u[i][j] + _v[i][j] * _v[i][j]));
	f_eq4 = W2 * _rho[i][j] * (1.0 - 3.0 * _v[i][j] + 4.5 * _v[i][j] * _v[i][j] - 1.5 * (_u[i][j] * _u[i][j] + _v[i][j] * _v[i][j]));
	f_eq5 = W3 * _rho[i][j] * (1.0 + 3.0 * (_u[i][j] + _v[i][j]) + 4.5 * (_u[i][j] + _v[i][j]) * (_u[i][j] + _v[i][j]) - 1.5 * (_u[i][j] * _u[i][j] + _v[i][j] * _v[i][j]));
	f_eq6 = W3 * _rho[i][j] * (1.0 + 3.0 * (-_u[i][j] + _v[i][j]) + 4.5 * (-_u[i][j] + _v[i][j]) * (-_u[i][j] + _v[i][j]) - 1.5 * (_u[i][j] * _u[i][j] + _v[i][j] * _v[i][j]));
	f_eq7 = W3 * _rho[i][j] * (1.0 - 3.0 * (_u[i][j] + _v[i][j]) + 4.5 * (_u[i][j] + _v[i][j]) * (_u[i][j] + _v[i][j]) - 1.5 * (_u[i][j] * _u[i][j] + _v[i][j] * _v[i][j]));
	f_eq8 = W3 * _rho[i][j] * (1.0 + 3.0 * (_u[i][j] - _v[i][j]) + 4.5 * (_u[i][j] - _v[i][j]) * (_u[i][j] - _v[i][j]) - 1.5 * (_u[i][j] * _u[i][j] + _v[i][j] * _v[i][j]));

	QSSize L = this->_lbblockowner->cTreeID().cDepth();

	double omega = 2 * _omega / ( pow(2, L + 1) + (1 - pow(2, L)) * _omega );

	_f_buf[0][i][j] = omega * f_eq0 + (1 - omega) * _f[0][i][j];
	_f_buf[1][i][j] = omega * f_eq1 + (1 - omega) * _f[1][i][j];
	_f_buf[2][i][j] = omega * f_eq2 + (1 - omega) * _f[2][i][j];
	_f_buf[3][i][j] = omega * f_eq3 + (1 - omega) * _f[3][i][j];
	_f_buf[4][i][j] = omega * f_eq4 + (1 - omega) * _f[4][i][j];
	_f_buf[5][i][j] = omega * f_eq5 + (1 - omega) * _f[5][i][j];
	_f_buf[6][i][j] = omega * f_eq6 + (1 - omega) * _f[6][i][j];
	_f_buf[7][i][j] = omega * f_eq7 + (1 - omega) * _f[7][i][j];
	_f_buf[8][i][j] = omega * f_eq8 + (1 - omega) * _f[8][i][j];
}

void LBBlock::EmptyCollision(QSSize i, QSSize j)
{
	_f_buf[0][i][j] = _f[0][i][j];
	_f_buf[1][i][j] = _f[1][i][j];
	_f_buf[2][i][j] = _f[2][i][j];
	_f_buf[3][i][j] = _f[3][i][j];
	_f_buf[4][i][j] = _f[4][i][j];
	_f_buf[5][i][j] = _f[5][i][j];
	_f_buf[6][i][j] = _f[6][i][j];
	_f_buf[7][i][j] = _f[7][i][j];
	_f_buf[8][i][j] = _f[8][i][j];
}

// ############################################################################################################# //

void LBBlock::ComputeDensity(const QSSize i, const QSSize j)
{
	_rho[i][j] = _f[0][i][j] + _f[1][i][j] + _f[2][i][j] + _f[3][i][j] + _f[4][i][j] + _f[5][i][j] + _f[6][i][j] + _f[7][i][j] + _f[8][i][j];
}

void LBBlock::ComputeVelocityX(const QSSize i, const QSSize j)
{
	_u[i][j] = (_f[1][i][j] + _f[5][i][j] + _f[8][i][j] - (_f[3][i][j] + _f[6][i][j] + _f[7][i][j])) / _rho[i][j];
}

void LBBlock::ComputeVelocityY(const QSSize i, const QSSize j)
{
	_v[i][j] = (_f[5][i][j] + _f[2][i][j] + _f[6][i][j] - (_f[7][i][j] + _f[4][i][j] + _f[8][i][j])) / _rho[i][j];
}

QSDouble LBBlock::CoordinateX(const QSSize i, const QSSize j){
	return _lbblockowner->cXmin() + (i-_ng) * (_lbblockowner->cXmax() - _lbblockowner->cXmin())/(_nx - 1);
}

QSDouble LBBlock::CoordinateY(const QSSize i, const QSSize j){
	return _lbblockowner->cYmin() + (j-_ng) * (_lbblockowner->cYmax() - _lbblockowner->cYmin())/(_ny - 1);
}

// ############################################################################################################# //

void LBBlock::NorthVelocity(const QSSize i, const QSSize j, const QSDouble u0, const QSDouble v0)
{
//	QSDouble densityN = _f[0][i][j] + _f[2][i][j] + _f[4][i][j] + 2 * (_f[5][i][j] + _f[2][i][j] + _f[6][i][j]);
	QSDouble densityN = _f[0][i][j] + _f[1][i][j] + _f[2][i][j] + _f[3][i][j] + _f[4][i][j] + _f[5][i][j] + _f[6][i][j] + _f[7][i][j] + _f[8][i][j];

	_f[4][i][j] = _f[2][i][j] - 2.0 / 3.0 * densityN * v0;
	_f[7][i][j] = _f[5][i][j] - 1.0 / 6.0 * densityN * v0 - 0.5 * u0 * densityN + 0.5 * (_f[1][i][j] - _f[3][i][j]);
	_f[8][i][j] = _f[6][i][j] - 1.0 / 6.0 * densityN * v0 + 0.5 * u0 * densityN + 0.5 * (_f[3][i][j] - _f[1][i][j]);
}

void LBBlock::SouthVelocity(const QSSize i, const QSSize j, const QSDouble u0, const QSDouble v0)
{
//	double densityN = _f[0][i][j] + _f[2][i][j] + _f[4][i][j] + 2 * (_f[7][i][j] + _f[4][i][j] + _f[8][i][j]);
	QSDouble densityN = _f[0][i][j] + _f[1][i][j] + _f[2][i][j] + _f[3][i][j] + _f[4][i][j] + _f[5][i][j] + _f[6][i][j] + _f[7][i][j] + _f[8][i][j];

	_f[2][i][j] = _f[4][i][j] + 2.0 / 3.0 * densityN * v0;
	_f[5][i][j] = _f[7][i][j] + 1.0 / 6.0 * densityN * v0 + 0.5 * u0 * densityN + 0.5 * (_f[3][i][j] - _f[1][i][j]);
	_f[6][i][j] = _f[8][i][j] + 1.0 / 6.0 * densityN * v0 - 0.5 * u0 * densityN + 0.5 * (_f[1][i][j] - _f[3][i][j]);
}

void LBBlock::EastVelocity(const QSSize i, const QSSize j, const QSDouble u0, const QSDouble v0)
{
//	double densityN = _f[0][i][j] + _f[1][i][j] + _f[3][i][j] + 2 * (_f[5][i][j] + _f[1][i][j] + _f[8][i][j]);
	QSDouble densityN = _f[0][i][j] + _f[1][i][j] + _f[2][i][j] + _f[3][i][j] + _f[4][i][j] + _f[5][i][j] + _f[6][i][j] + _f[7][i][j] + _f[8][i][j];

	_f[3][i][j] = _f[1][i][j] - 2.0 / 3.0 * densityN * u0;
	_f[6][i][j] = _f[8][i][j] - 1.0 / 6.0 * densityN * u0 + 0.5 * v0 * densityN + 0.5 * (_f[4][i][j] - _f[2][i][j]);
	_f[7][i][j] = _f[5][i][j] - 1.0 / 6.0 * densityN * u0 - 0.5 * v0 * densityN + 0.5 * (_f[2][i][j] - _f[4][i][j]);
}

void LBBlock::WestVelocity(const QSSize i, const QSSize j, QSDouble u0, QSDouble v0)
{
	QSDouble densityN = _f[0][i][j] + _f[1][i][j] + _f[2][i][j] + _f[3][i][j] + _f[4][i][j] + _f[5][i][j] + _f[6][i][j] + _f[7][i][j] + _f[8][i][j];

	QSDouble y = CoordinateY(i, j);

	u0 = (-y * (y - 1) ) / 4;

	_f[1][i][j] = _f[3][i][j] + 2.0 / 3.0 * densityN * u0;
	_f[5][i][j] = _f[7][i][j] + 1.0 / 6.0 * densityN * u0 + 0.5 * v0 * densityN + 0.5 * (_f[4][i][j] - _f[2][i][j]);
	_f[8][i][j] = _f[6][i][j] + 1.0 / 6.0 * densityN * u0 - 0.5 * v0 * densityN + 0.5 * (_f[2][i][j] - _f[4][i][j]);
}

// ############################################################################################################# //

void LBBlock::NorthPressure(const QSSize i, const QSSize j, const QSDouble p0)
{
	QSDouble uN = (_f[1][i][j] + _f[5][i][j] + _f[8][i][j] - (_f[3][i][j] + _f[6][i][j] + _f[7][i][j])) / p0;
	QSDouble vN = (_f[5][i][j] + _f[2][i][j] + _f[6][i][j] - (_f[7][i][j] + _f[4][i][j] + _f[8][i][j])) / p0;

	_f[4][i][j] = _f[2][i][j] - 2.0 / 3.0 * p0 * vN;
	_f[7][i][j] = _f[5][i][j] - 1.0 / 6.0 * p0 * vN - 0.5 * uN * p0 + 0.5 * (_f[1][i][j] - _f[3][i][j]);
	_f[8][i][j] = _f[6][i][j] - 1.0 / 6.0 * p0 * vN + 0.5 * uN * p0 + 0.5 * (_f[3][i][j] - _f[1][i][j]);
}

void LBBlock::SouthPressure(const QSSize i, const QSSize j, const QSDouble p0)
{
	QSDouble uN = (_f[1][i][j] + _f[5][i][j] + _f[8][i][j] - (_f[3][i][j] + _f[6][i][j] + _f[7][i][j])) / p0;
	QSDouble vN = (_f[5][i][j] + _f[2][i][j] + _f[6][i][j] - (_f[7][i][j] + _f[4][i][j] + _f[8][i][j])) / p0;

	_f[2][i][j] = _f[4][i][j] + 2.0 / 3.0 * p0 * vN;
	_f[5][i][j] = _f[7][i][j] + 1.0 / 6.0 * p0 * vN + 0.5 * uN * p0 + 0.5 * (_f[3][i][j] - _f[1][i][j]);
	_f[6][i][j] = _f[8][i][j] + 1.0 / 6.0 * p0 * vN - 0.5 * uN * p0 + 0.5 * (_f[1][i][j] - _f[3][i][j]);
}

void LBBlock::EastPressure(const QSSize i, const QSSize j, const QSDouble p0)
{
	QSDouble uN = (_f[1][i][j] + _f[5][i][j] + _f[8][i][j] - (_f[3][i][j] + _f[6][i][j] + _f[7][i][j])) / p0;
	QSDouble vN = (_f[5][i][j] + _f[2][i][j] + _f[6][i][j] - (_f[7][i][j] + _f[4][i][j] + _f[8][i][j])) / p0;

	_f[3][i][j] = _f[1][i][j] - 2.0 / 3.0 * p0 * uN;
	_f[6][i][j] = _f[8][i][j] - 1.0 / 6.0 * p0 * uN + 0.5 * vN * p0 + 0.5 * (_f[4][i][j] - _f[2][i][j]);
	_f[7][i][j] = _f[5][i][j] - 1.0 / 6.0 * p0 * uN - 0.5 * vN * p0 + 0.5 * (_f[2][i][j] - _f[4][i][j]);
}

void LBBlock::WestPressure(const QSSize i, const QSSize j, const QSDouble p0)
{
	QSDouble uN = (_f[1][i][j] + _f[5][i][j] + _f[8][i][j] - (_f[3][i][j] + _f[6][i][j] + _f[7][i][j])) / p0;
	QSDouble vN = (_f[5][i][j] + _f[2][i][j] + _f[6][i][j] - (_f[7][i][j] + _f[4][i][j] + _f[8][i][j])) / p0;

	_f[1][i][j] = _f[3][i][j] + 2.0 / 3.0 * p0 * uN;
	_f[5][i][j] = _f[7][i][j] + 1.0 / 6.0 * p0 * uN + 0.5 * vN * p0 + 0.5 * (_f[4][i][j] - _f[2][i][j]);
	_f[8][i][j] = _f[6][i][j] + 1.0 / 6.0 * p0 * uN - 0.5 * vN * p0 + 0.5 * (_f[2][i][j] - _f[4][i][j]);
}

// ############################################################################################################# //

void LBBlock::BounceBack(const QSSize i, const QSSize j)
{
	_f_buf[0][i][j] = _f[0][i][j];
	_f_buf[1][i][j] = _f[3][i][j];
	_f_buf[2][i][j] = _f[4][i][j];
	_f_buf[3][i][j] = _f[1][i][j];
	_f_buf[4][i][j] = _f[2][i][j];
	_f_buf[5][i][j] = _f[7][i][j];
	_f_buf[6][i][j] = _f[8][i][j];
	_f_buf[7][i][j] = _f[5][i][j];
	_f_buf[8][i][j] = _f[6][i][j];
}

// ############################################################################################################# //

void LBBlock::Streaming(const QSSize i, const QSSize j)
{
	_f[0][i][j] = _f_buf[0][i][j];
	_f[1][i][j] = _f_buf[1][i - 1][j];
	_f[2][i][j] = _f_buf[2][i][j - 1];
	_f[3][i][j] = _f_buf[3][i + 1][j];
	_f[4][i][j] = _f_buf[4][i][j + 1];
	_f[5][i][j] = _f_buf[5][i - 1][j - 1];
	_f[6][i][j] = _f_buf[6][i + 1][j - 1];
	_f[7][i][j] = _f_buf[7][i + 1][j + 1];
	_f[8][i][j] = _f_buf[8][i - 1][j + 1];
}

void LBBlock::Streamingf0(const QSSize i, const QSSize j)
{
	_f[0][i][j] = _f_buf[0][i][j];
}

void LBBlock::Streamingf1(const QSSize i, const QSSize j)
{
	_f[1][i][j] = _f_buf[1][i - 1][j];
}

void LBBlock::Streamingf2(const QSSize i, const QSSize j)
{
	_f[2][i][j] = _f_buf[2][i][j - 1];
}

void LBBlock::Streamingf3(const QSSize i, const QSSize j)
{
	_f[3][i][j] = _f_buf[3][i + 1][j];
}

void LBBlock::Streamingf4(const QSSize i, const QSSize j)
{
	_f[4][i][j] = _f_buf[4][i][j + 1];
}

void LBBlock::Streamingf5(const QSSize i, const QSSize j)
{
	_f[5][i][j] = _f_buf[5][i - 1][j - 1];
}

void LBBlock::Streamingf6(const QSSize i, const QSSize j)
{
	_f[6][i][j] = _f_buf[6][i + 1][j - 1];
}

void LBBlock::Streamingf7(const QSSize i, const QSSize j)
{
	_f[7][i][j] = _f_buf[7][i + 1][j + 1];
}

void LBBlock::Streamingf8(const QSSize i, const QSSize j)
{
	_f[8][i][j] = _f_buf[8][i - 1][j + 1];
}

// ############################################################################################################# //

void LBBlock::GhostInfo()
{
	// NORTH
	if (_north.size() == 1) {
		if (_lbblockowner->cTreeID().cDepth() == _north[0]->_lbblockowner->rTreeID().cDepth()) {
			for (QSSize i = _ng; i < _nx + _ng; ++i) {
				for (QSSize j = _nx + _ng; j < _nx + 2 * _ng; ++j) {
					_f_buf[0][i][j] = _north[0]->_f_buf[0][i][j - _ny];
					_f_buf[1][i][j] = _north[0]->_f_buf[1][i][j - _ny];
					_f_buf[2][i][j] = _north[0]->_f_buf[2][i][j - _ny];
					_f_buf[3][i][j] = _north[0]->_f_buf[3][i][j - _ny];
					_f_buf[4][i][j] = _north[0]->_f_buf[4][i][j - _ny];
					_f_buf[5][i][j] = _north[0]->_f_buf[5][i][j - _ny];
					_f_buf[6][i][j] = _north[0]->_f_buf[6][i][j - _ny];
					_f_buf[7][i][j] = _north[0]->_f_buf[7][i][j - _ny];
					_f_buf[8][i][j] = _north[0]->_f_buf[8][i][j - _ny];
        }
			}
		}
	}

	// EAST
	if (_east.size() == 1) {
		if (_lbblockowner->cTreeID().cDepth() == _east[0]->_lbblockowner->cTreeID().cDepth()) {
			for (QSSize i = _nx + _ng; i < _nx + 2 * _ng; ++i) {
				for (QSSize j = _ng; j < _ny + _ng; ++j) {
					_f_buf[0][i][j] = _east[0]->_f_buf[0][i - _nx][j];
					_f_buf[1][i][j] = _east[0]->_f_buf[1][i - _nx][j];
					_f_buf[2][i][j] = _east[0]->_f_buf[2][i - _nx][j];
					_f_buf[3][i][j] = _east[0]->_f_buf[3][i - _nx][j];
					_f_buf[4][i][j] = _east[0]->_f_buf[4][i - _nx][j];
					_f_buf[5][i][j] = _east[0]->_f_buf[5][i - _nx][j];
					_f_buf[6][i][j] = _east[0]->_f_buf[6][i - _nx][j];
					_f_buf[7][i][j] = _east[0]->_f_buf[7][i - _nx][j];
					_f_buf[8][i][j] = _east[0]->_f_buf[8][i - _nx][j];
        }
			}
		}
	}

	// SOUTH
	if (_south.size() == 1) {
		if (_lbblockowner->cTreeID().cDepth() == _south[0]->_lbblockowner->cTreeID().cDepth()) {
			for (QSSize i = _ng; i < _nx + _ng; ++i) {
				for (QSSize j = 0; j < _ng; ++j) {
					_f_buf[0][i][j] = _south[0]->_f_buf[0][i][j + _ny];
					_f_buf[1][i][j] = _south[0]->_f_buf[1][i][j + _ny];
					_f_buf[2][i][j] = _south[0]->_f_buf[2][i][j + _ny];
					_f_buf[3][i][j] = _south[0]->_f_buf[3][i][j + _ny];
					_f_buf[4][i][j] = _south[0]->_f_buf[4][i][j + _ny];
					_f_buf[5][i][j] = _south[0]->_f_buf[5][i][j + _ny];
					_f_buf[6][i][j] = _south[0]->_f_buf[6][i][j + _ny];
					_f_buf[7][i][j] = _south[0]->_f_buf[7][i][j + _ny];
					_f_buf[8][i][j] = _south[0]->_f_buf[8][i][j + _ny];
        }
			}
		}
	}

	// WEST
	if (_west.size() == 1) {
		if (_lbblockowner->cTreeID().cDepth() == _west[0]->_lbblockowner->cTreeID().cDepth()) {
			for (QSSize i = 0; i < _ng; ++i) {
				for (QSSize j = _ng; j < _ny + _ng; ++j) {
					_f_buf[0][i][j] = _west[0]->_f_buf[0][i + _nx][j];
					_f_buf[1][i][j] = _west[0]->_f_buf[1][i + _nx][j];
					_f_buf[2][i][j] = _west[0]->_f_buf[2][i + _nx][j];
					_f_buf[3][i][j] = _west[0]->_f_buf[3][i + _nx][j];
					_f_buf[4][i][j] = _west[0]->_f_buf[4][i + _nx][j];
					_f_buf[5][i][j] = _west[0]->_f_buf[5][i + _nx][j];
					_f_buf[6][i][j] = _west[0]->_f_buf[6][i + _nx][j];
					_f_buf[7][i][j] = _west[0]->_f_buf[7][i + _nx][j];
					_f_buf[8][i][j] = _west[0]->_f_buf[8][i + _nx][j];
        }
			}
		}
	}

	// NORTHWEST
	if (_northwest.size() == 1) {
		if (_lbblockowner->cTreeID().cDepth() == _northwest[0]->_lbblockowner->cTreeID().cDepth()) {
			for (QSSize i = 0; i < _ng; ++i) {
				for (QSSize j = _ny + _ng; j < _ny + 2 * _ng; ++j) {
					_f_buf[0][i][j] = _northwest[0]->_f_buf[0][i + _nx][j - _ny];
					_f_buf[1][i][j] = _northwest[0]->_f_buf[1][i + _nx][j - _ny];
					_f_buf[2][i][j] = _northwest[0]->_f_buf[2][i + _nx][j - _ny];
					_f_buf[3][i][j] = _northwest[0]->_f_buf[3][i + _nx][j - _ny];
					_f_buf[4][i][j] = _northwest[0]->_f_buf[4][i + _nx][j - _ny];
					_f_buf[5][i][j] = _northwest[0]->_f_buf[5][i + _nx][j - _ny];
					_f_buf[6][i][j] = _northwest[0]->_f_buf[6][i + _nx][j - _ny];
					_f_buf[7][i][j] = _northwest[0]->_f_buf[7][i + _nx][j - _ny];
					_f_buf[8][i][j] = _northwest[0]->_f_buf[8][i + _nx][j - _ny];
        }
			}
		}
	}

	// NORTHEAST
	if (_northeast.size() == 1) {
		if (_lbblockowner->cTreeID().cDepth() == _northeast[0]->_lbblockowner->cTreeID().cDepth()) {
			for (QSSize i = _nx + _ng; i < _nx + 2 * _ng; ++i) {
				for (QSSize j = _ny + _ng; j < _ny + 2 * _ng; ++j) {
					_f_buf[0][i][j] = _northeast[0]->_f_buf[0][i - _nx][j - _ny];
					_f_buf[1][i][j] = _northeast[0]->_f_buf[1][i - _nx][j - _ny];
					_f_buf[2][i][j] = _northeast[0]->_f_buf[2][i - _nx][j - _ny];
					_f_buf[3][i][j] = _northeast[0]->_f_buf[3][i - _nx][j - _ny];
					_f_buf[4][i][j] = _northeast[0]->_f_buf[4][i - _nx][j - _ny];
					_f_buf[5][i][j] = _northeast[0]->_f_buf[5][i - _nx][j - _ny];
					_f_buf[6][i][j] = _northeast[0]->_f_buf[6][i - _nx][j - _ny];
					_f_buf[7][i][j] = _northeast[0]->_f_buf[7][i - _nx][j - _ny];
					_f_buf[8][i][j] = _northeast[0]->_f_buf[8][i - _nx][j - _ny];
        }
			}
		}
	}

	// SOUTHWEST
	if (_southwest.size() == 1) {
		if (_lbblockowner->cTreeID().cDepth() == _southwest[0]->_lbblockowner->cTreeID().cDepth()) {
			for (QSSize i = 0; i < _ng; ++i) {
				for (QSSize j = 0; j < _ng; ++j) {
					_f_buf[0][i][j] = _southwest[0]->_f_buf[0][i + _nx][j + _ny];
					_f_buf[1][i][j] = _southwest[0]->_f_buf[1][i + _nx][j + _ny];
					_f_buf[2][i][j] = _southwest[0]->_f_buf[2][i + _nx][j + _ny];
					_f_buf[3][i][j] = _southwest[0]->_f_buf[3][i + _nx][j + _ny];
					_f_buf[4][i][j] = _southwest[0]->_f_buf[4][i + _nx][j + _ny];
					_f_buf[5][i][j] = _southwest[0]->_f_buf[5][i + _nx][j + _ny];
					_f_buf[6][i][j] = _southwest[0]->_f_buf[6][i + _nx][j + _ny];
					_f_buf[7][i][j] = _southwest[0]->_f_buf[7][i + _nx][j + _ny];
					_f_buf[8][i][j] = _southwest[0]->_f_buf[8][i + _nx][j + _ny];
        }
			}
		}
	}

	// SOUTHEAST
	if (_southeast.size() == 1) {
		if (_lbblockowner->cTreeID().cDepth() == _southeast[0]->_lbblockowner->cTreeID().cDepth()) {
			for (QSSize i = _nx + _ng; i < _nx + 2 * _ng; ++i) {
				for (QSSize j = 0; j < _ng; ++j) {
					_f_buf[0][i][j] = _southeast[0]->_f_buf[0][i - _nx][j + _ny];
					_f_buf[1][i][j] = _southeast[0]->_f_buf[1][i - _nx][j + _ny];
					_f_buf[2][i][j] = _southeast[0]->_f_buf[2][i - _nx][j + _ny];
					_f_buf[3][i][j] = _southeast[0]->_f_buf[3][i - _nx][j + _ny];
					_f_buf[4][i][j] = _southeast[0]->_f_buf[4][i - _nx][j + _ny];
					_f_buf[5][i][j] = _southeast[0]->_f_buf[5][i - _nx][j + _ny];
					_f_buf[6][i][j] = _southeast[0]->_f_buf[6][i - _nx][j + _ny];
					_f_buf[7][i][j] = _southeast[0]->_f_buf[7][i - _nx][j + _ny];
					_f_buf[8][i][j] = _southeast[0]->_f_buf[8][i - _nx][j + _ny];
        }
			}
		}
	}
}

void LBBlock::WriteToVTI(QSSize lp, QSSize time)
{
	const QSDouble dx = (_lbblockowner->cXmax() - _lbblockowner->cXmin()) / _nx;
	const QSDouble dy = (_lbblockowner->cYmax() - _lbblockowner->cYmin()) / _ny;

	std::stringstream filename;
	std::stringstream foldername;

	foldername << "t_" << time;
	filename << "t_" << time << "/sym_" << lp << ".vti";

	mkdir(foldername.str().c_str(), 0700);
	QSFileStream file(filename.str(), std::ios::out | std::ios::binary);

	if (file.good()) {
		file << "<VTKFile type=\"ImageData\" version=\"0.1\" " "byte_order=\"LittleEndian\">" << std::endl;
		file << "<ImageData WholeExtent=\"0 " << _nx << " 0 " << _ny << " 0 0\" Origin=\"" << _lbblockowner->cXmin() << " " << _lbblockowner->cYmin() << " 0\" Spacing=\"" << dx << " " << dy
				 << " " << "0.01\">" << std::endl;
		file << "<CellData Scalars=\"Pressure\" Vectors=\"Velocity\">" << std::endl;

		file << "<DataArray type=\"Float64\" Name=\"Pressure\" format=\"ascii\">" << std::endl;

		for (QSSize j = _ng; j < _ny + _ng; j++) {
			for (QSSize i = _ng; i < _nx + _ng; i++) {
				ComputeDensity(i, j);
				file << _rho[i][j] << std::endl;
			}
		}
		file << std::endl << "</DataArray>" << std::endl;

		file << "<DataArray type=\"Float64\" Name=\"Velocity\" format=\"ascii\" " "NumberOfComponents=\"3\">" << std::endl;

		for (QSSize j = _ng; j < _ny + _ng; j++) {
			for (QSSize i = _ng; i < _nx + _ng; i++) {
				ComputeVelocityX(i, j);
				ComputeVelocityY(i, j);
				file << _u[i][j] << " " << _v[i][j] << " 0" << std::endl;
			}
		}

		file << std::endl << "</DataArray>" << std::endl;

		file << "<DataArray type=\"Int32\" Name=\"Flag\" format=\"ascii\">" << std::endl;
		for (QSSize j = _ng; j < _ny + _ng; j++) {
			for (QSSize i = _ng; i < _nx + _ng; i++) {
				file << _flag[i][j] << std::endl;
			}
    }
		file << std::endl << "</DataArray>" << std::endl;

		file << "</CellData>" << std::endl;
		file << "</ImageData>" << std::endl;
		file << "</VTKFile>" << std::endl;
	} else {
		throw;
	}

	file.close();
}

// ############################################################################################################# //
// TESTING METHODS

void LBBlock::InitZero(QSSize nx, QSSize ny, QSSize ng, QSDouble omega, QSDouble uini, QSDouble vini)
{
	QSDouble W1 = 4.0 / 9.0;
	QSDouble W2 = 1.0 / 9.0;
	QSDouble W3 = 1.0 / 36.0;

	std::vector<QSDouble> weigths;

	weigths.resize(9);

	weigths[0] = W1;
	weigths[1] = W2;
	weigths[2] = W2;
	weigths[3] = W2;
	weigths[4] = W2;
	weigths[5] = W3;
	weigths[6] = W3;
	weigths[7] = W3;
	weigths[8] = W3;

	_nx = nx;
	_ny = ny;
	_ng = ng;
	_omega = omega;

	_u_ini = uini;
	_v_ini = vini;

	_f.resize(9);

	for (QSSize i = 0; i < 9; ++i) {
		_f[i].resize(_nx + 2 * _ng);
		for (auto& j : _f[i]) {
			j.resize(_ny + 2 * _ng);
    }
	}

	_f_buf.resize(9);

	for (QSSize i = 0; i < 9; ++i) {
		_f_buf[i].resize(_nx + 2 * _ng);
		for (auto& j : _f_buf[i]) {
			j.resize(_ny + 2 * _ng);
    }
	}

	for (QSSize i = 0; i < _nx + 2 * _ng; ++i) {
		for (QSSize j = 0; j < _ny + 2 * _ng; ++j) {
			for (QSSize k = 0; k < 9; ++k) {
				_f[k][i][j] = 1;
				_f_buf[k][i][j] = 1;
			}
    }
	}

	_flag.resize(_nx + 2 * _ng);

	for (auto& i : _flag) {
		i.resize(_ny + 2 * _ng);
	}

	for (QSSize j = 0; j < _ny + 2 * _ng; ++j) {
		for (QSSize i = 0; i < _nx + 2 * _ng; ++i) {
			if (i < _ng || i >= _nx + _ng || j < _ng || j >= _ny + _ng) {
				_flag[i][j] = LBNodeType::GHOST;
			} else if (fabs(_lbblockowner->cYmax() - 1) < 1e-6 && j == _ny + _ng - 1) { // top
				_flag[i][j] = LBNodeType::NORTHVELOCITY;
			} else if (fabs(_lbblockowner->cXmin()) < 1e-6 && i == _ng) { // left
				_flag[i][j] = LBNodeType::WALLBB;
			} else if (fabs(_lbblockowner->cXmax() - 1) < 1e-6 && i == _nx + _ng - 1) { // right
				_flag[i][j] = LBNodeType::WALLBB;
			} else if (fabs(_lbblockowner->cYmin()) < 1e-6 && j == _ng) { // bottom
				_flag[i][j] = LBNodeType::WALLBB;
			} else
				_flag[i][j] = INTERIOR;
    }
	}

	_u.resize(_nx + 2 * _ng);
	for (auto& i : _u) {
		i.resize(_ny + 2 * _ng);
	}

	_v.resize(_nx + 2 * _ng);
	for (auto& i : _v) {
		i.resize(_ny + 2 * _ng);
	}

	_rho.resize(_nx + 2 * _ng);
	for (auto& i : _rho) {
		i.resize(_ny + 2 * _ng);
	}
}

void LBBlock::InitWave()
{
	QSDouble W1 = 4.0 / 9.0;
	QSDouble W2 = 1.0 / 9.0;
	QSDouble W3 = 1.0 / 36.0;

	_f[0][(_nx + 2 * _ng) / 2][(_ny + 2 * _ng) / 2] = W1;
	_f[1][(_nx + 2 * _ng) / 2][(_ny + 2 * _ng) / 2] = W2;
	_f[2][(_nx + 2 * _ng) / 2][(_ny + 2 * _ng) / 2] = W2;
	_f[3][(_nx + 2 * _ng) / 2][(_ny + 2 * _ng) / 2] = W2;
	_f[4][(_nx + 2 * _ng) / 2][(_ny + 2 * _ng) / 2] = W2;
	_f[5][(_nx + 2 * _ng) / 2][(_ny + 2 * _ng) / 2] = W3;
	_f[6][(_nx + 2 * _ng) / 2][(_ny + 2 * _ng) / 2] = W3;
	_f[7][(_nx + 2 * _ng) / 2][(_ny + 2 * _ng) / 2] = W3;
	_f[8][(_nx + 2 * _ng) / 2][(_ny + 2 * _ng) / 2] = W3;
}

void LBBlock::InitWaveZero()
{
	QSDouble W2 = 1.0 / 9.0;
	QSDouble W3 = 1.0 / 36.0;

	_f[0][(_nx + 2 * _ng) / 2][(_ny + 2 * _ng) / 2] = 1;
	_f[1][(_nx + 2 * _ng) / 2][(_ny + 2 * _ng) / 2] = W2;
	_f[2][(_nx + 2 * _ng) / 2][(_ny + 2 * _ng) / 2] = W2;
	_f[3][(_nx + 2 * _ng) / 2][(_ny + 2 * _ng) / 2] = W2;
	_f[4][(_nx + 2 * _ng) / 2][(_ny + 2 * _ng) / 2] = W2;
	_f[5][(_nx + 2 * _ng) / 2][(_ny + 2 * _ng) / 2] = W3;
	_f[6][(_nx + 2 * _ng) / 2][(_ny + 2 * _ng) / 2] = W3;
	_f[7][(_nx + 2 * _ng) / 2][(_ny + 2 * _ng) / 2] = W3;
	_f[8][(_nx + 2 * _ng) / 2][(_ny + 2 * _ng) / 2] = W3;
}

void LBBlock::EmptyCollisionStep()
{
	for (QSSize i = 0; i < _nx + 2 * _ng; ++i) {
		for (QSSize j = 0; j < _ny + 2 * _ng; ++j) {
			switch (_flag[i][j]) {
			case LBNodeType::INTERIOR: // fluid
				EmptyCollision(i, j);
				break;

			case LBNodeType::NORTHVELOCITY: // ZuoHe
				BounceBack(i, j);
				break;

			case LBNodeType::WALLBB: // BounceBack
				BounceBack(i, j);
				break;

			case LBNodeType::GHOST:
				EmptyCollision(i, j);
				break;

			default:
				throw;
			}
		}
	}
}
