/*---------------------------------------------------------------------------*\
  Compatibility for the OpenFOAM Windows runtime.

  GCC 13 MinGW <iostream> emits an undefined reference to
  std::ios_base_library_init. The win64Mingw OpenFOAM package ships the
  GCC 7 libstdc++, which does not export that symbol, and LoadLibrary then
  fails. The definition below keeps the reference inside this DLL.
  Iostream objects are still constructed by the libstdc++ that is loaded
  with the OpenFOAM libraries.
\*---------------------------------------------------------------------------*/

namespace std
{

void ios_base_library_init()
{}

}

// ************************************************************************* //
