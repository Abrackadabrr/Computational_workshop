#ifndef MESHUTILS_HPP
#define MESHUTILS_HPP

#include "inmost.h"

namespace Mesh {

enum CellType {
  cut = 1,
  endocard = 2,
  epicard = 3,
  central_line = -1,
};

}

#endif // MESHUTILS_HPP
