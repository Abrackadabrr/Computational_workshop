#include "integration/Gaussian.hpp"
#include "integration/Quadrature.hpp"
#include "mesh/MeshUtils.hpp"
#include "types/BasicTypes.hpp"

#include <inmost.h>
#include <iostream>

Types::scalar function(const Types::array_s<3> &x) { return 1; }

int main() {
  INMOST::Mesh mesh{};
  Mesh::CellType cell_type = Mesh::CellType::central_line;
  mesh.LoadMSH("/home/evgen/Education/MasterDegree/5_level/Chapter_1/INM_RAS/"
               "Computing_workshop/Workshop/Lab1/ventrical.msh");
  Types::scalar res = Math::Integration::integrate<
      Math::Integration::QuadratureTypes::GaussianPoints<1>>(function, mesh,
                                                             cell_type);
  std::cout << res << std::endl;

  // проверка на работоспособность
  int amount = 0;
  for (auto ielem = mesh.BeginElement(INMOST::CELL | INMOST::FACE);
       ielem != mesh.EndElement(); ++ielem) {
    if (cell_type == ielem->Integer(mesh.GetTag("GMSH_TAGS")))
      amount++;
  };
  std::cout << amount << std::endl;
}
