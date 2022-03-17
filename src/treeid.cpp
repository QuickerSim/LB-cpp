#include "treeid.hpp"

// TODO - rework of whole file - too many implicit conversions
TreeID TreeID::moveN() const
{
  std::bitset<sizeof(QSSize) * 8> idbits(_id - qtelemsum[_depth - 1]);
  bool carry = true;
  for (int i = 1; i < sizeof(QSSize) * 8; i = i + 2) {
    if (carry == true) {
      if (idbits[i] == true) {
        idbits.set(i, false);
      } else {
        idbits.set(i, true);
        carry = false;
      }
    } else {
      break;
    }
  }
	if (idbits.to_ulong() >= qtelemsum[_depth] - qtelemsum[_depth - 1]) {
    return TreeID(0, 0);
  }
  return TreeID(idbits.to_ulong() + qtelemsum[_depth - 1], _depth);
}

TreeID TreeID::moveE() const
{
  std::bitset<sizeof(QSSize) * 8> idbits(_id - qtelemsum[_depth - 1]);

  bool carry = true;
  for (int i = 0; i < sizeof(QSSize) * 8; i = i + 2) {
    if (carry == true) {
      if (idbits[i] == true) {
        idbits.set(i, false);
      } else {
        idbits.set(i, true);
        carry = false;
      }
    } else {
      break;
    }
  }
	if (idbits.to_ulong() >= qtelemsum[_depth] - qtelemsum[_depth - 1]) {
    return TreeID(0, 0);
  }
  return TreeID(idbits.to_ulong() + qtelemsum[_depth - 1], _depth);
}

TreeID TreeID::moveS() const
{
  std::bitset<sizeof(QSSize) * 8> idbits(_id - qtelemsum[_depth - 1]);

  bool carry = true;
  for (int i = 1; i < 64; i = i + 2) {
    if (carry == true) {
      if (idbits[i] == true) {
        if (i == 1) {
          idbits.set(i, false);
          carry = false;
          break;
        }

        idbits.set(i, false);
        for (int j = i - 2; j >= 0; j = j - 2) {
          if (idbits[j] == false) {
            idbits.set(j, true);
          } else {
            idbits.set(j, false);
            carry = false;
            break;
          }

          if (j <= 1) {
            carry = false;
          }
        }
      }
    } else {
      break;
    }
  }

  if (carry == true) {
    return TreeID(0, 0);
  }
  return TreeID(idbits.to_ulong() + qtelemsum[_depth - 1], _depth);
}

TreeID TreeID::moveW() const
{
  std::bitset<sizeof(QSSize) * 8> idbits(_id - qtelemsum[_depth - 1]);

  bool carry = true;
  for (int i = 0; i < 64; i = i + 2) {
    if (carry == true) {
      if (idbits[i] == true) {
        if (i == 0) {
          idbits.set(i, false);
          carry = false;
          break;
        }

        idbits.set(i, false);
        for (int j = i - 2; j >= 0; j = j - 2) {
          if (idbits[j] == false) {
            idbits.set(j, true);
          } else {
            idbits.set(j, false);
            carry = false;
            break;
          }
          if (j <= 1) {
            carry = false;
          }
        }
      }
    } else {
      break;
    }
  }
  if (carry == true) {
    return TreeID(0, 0);
  }

  return TreeID(idbits.to_ulong() + qtelemsum[_depth - 1], _depth);
}

TreeID TreeID::moveNW() const
{
  TreeID tmptree = this->moveN();

  if (tmptree.cDepth() == 0) {
    return TreeID(0, 0);
  }

  return tmptree.moveW();
}

TreeID TreeID::moveNE() const
{
  TreeID tmptree = this->moveN();

  if (tmptree.cDepth() == 0) {
    return TreeID(0, 0);
  }

  return tmptree.moveE();
}

TreeID TreeID::moveSW() const
{
  TreeID tmptree = this->moveS();

  if (tmptree.cDepth() == 0) {
    return TreeID(0, 0);
  }

  return tmptree.moveW();
}

TreeID TreeID::moveSE() const
{
  TreeID tmptree = this->moveS();

  if (tmptree.cDepth() == 0) {
    return TreeID(0, 0);
  }

  return tmptree.moveE();
}

TreeID TreeID::Parent()
{
  return TreeID(
      static_cast<QSSize>(floor((static_cast<QSDouble>(_id) - 1.0) / 4.0)),
      _depth + 1);
}

std::array<TreeID, 4> TreeID::Childs()
{
  std::array<TreeID, 4> childs;

  for (QSSize i = 0; i < 4; ++i) {
    childs[i] = TreeID(4 * _id + i + 1, _depth + 1);
  }

  return childs;
}
