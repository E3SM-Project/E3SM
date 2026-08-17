//===-- drivers/coupled/MoabInterface.cpp - MOAB coupling bridge
//-----------===//
//
//===---------------------------------------------------------------------===//

#include "MoabInterface.h"

#include "Decomp.h"
#include "GlobalConstants.h"
#include "HorzMesh.h"
#include "Logging.h"
#include "MachEnv.h"
#include "OmegaKokkos.h"

#include "moab/iMOAB.h"

#include <cstdio>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <vector>

namespace OMEGA {

namespace {

int Pid = -1;

// Full CIME x2o/o2x field-name lists (colon-separated, null-terminated),
// captured by moabDefineTagStorage and reused by moabImportTagStorage/
// moabExportTagStorage.
std::string X2oTagNames;
std::string O2xTagNames;

void checkMoabErr(ErrCode Err, const std::string &Msg) {
   if (Err != moab::MB_SUCCESS)
      LOG_CRITICAL("MOAB error in {}: code {}", Msg, static_cast<int>(Err));
}

// iMOAB's tag storage API is double-only. These are templates (not plain
// `if constexpr` in a non-template function) specifically so the untaken
// branch is genuinely discarded-and-not-instantiated for the concrete T,
// rather than needing to type-check for every T in every build: when T is
// double (the common case) the buffer is read/written directly with no
// intermediate copy; only a build where Real != double needs the scratch
// conversion.
template <typename T>
void getDoubleTagStorage(const std::string &TagNames, T *Buffer, int Len) {
   int EntType = 1; // elements
   if constexpr (std::is_same_v<T, double>) {
      ErrCode Err = iMOAB_GetDoubleTagStorage(&Pid, TagNames.c_str(), &Len,
                                              &EntType, Buffer);
      checkMoabErr(Err, "iMOAB_GetDoubleTagStorage");
   } else {
      std::vector<double> Scratch(Len);
      ErrCode Err = iMOAB_GetDoubleTagStorage(&Pid, TagNames.c_str(), &Len,
                                              &EntType, Scratch.data());
      checkMoabErr(Err, "iMOAB_GetDoubleTagStorage");
      for (int I = 0; I < Len; ++I)
         Buffer[I] = static_cast<T>(Scratch[I]);
   }
}

template <typename T>
void setDoubleTagStorage(const std::string &TagNames, const T *Buffer,
                         int Len) {
   int EntType = 1; // elements
   if constexpr (std::is_same_v<T, double>) {
      // iMOAB_SetDoubleTagStorage takes a non-const double* even though it
      // only reads from it; the const_cast is safe, not a real mutation.
      ErrCode Err = iMOAB_SetDoubleTagStorage(&Pid, TagNames.c_str(), &Len,
                                              &EntType,
                                              const_cast<double *>(Buffer));
      checkMoabErr(Err, "iMOAB_SetDoubleTagStorage");
   } else {
      std::vector<double> Scratch(Len);
      for (int I = 0; I < Len; ++I)
         Scratch[I] = static_cast<double>(Buffer[I]);
      ErrCode Err = iMOAB_SetDoubleTagStorage(&Pid, TagNames.c_str(), &Len,
                                              &EntType, Scratch.data());
      checkMoabErr(Err, "iMOAB_SetDoubleTagStorage");
   }
}

// Registers this ocean instance as a MOAB application and returns its pid
int registerApplication(MPI_Comm Comm, int OcnID) {

   ErrCode Err;
   int LocalPid;
   int CompID = OcnID;

   char OcnIDStr[8];
   std::snprintf(OcnIDStr, sizeof(OcnIDStr), "%02d", OcnID);
   std::string AppName = std::string("OMEGA_MB_") + OcnIDStr;

   Err = iMOAB_RegisterApplication(AppName.c_str(), &Comm, &CompID, &LocalPid);
   checkMoabErr(Err, "iMOAB_RegisterApplication");

   return LocalPid;
}

// Defines a dense double tag with 1 component per entity on elements and
// sets it from Data (size NCellsOwned).
void defineAndSetElemDoubleTag(int LocalPid, const std::string &TagName,
                                int NCellsOwned, const double *Data) {
   ErrCode Err;
   int TagType = DENSE_DOUBLE;
   int NumCo   = 1;
   int TagIndex;

   Err = iMOAB_DefineTagStorage(&LocalPid, TagName.c_str(), &TagType, &NumCo,
                                &TagIndex);
   checkMoabErr(Err, "iMOAB_DefineTagStorage(" + TagName + ")");

   int EntType = 1; // elements
   int Len     = NCellsOwned;
   Err = iMOAB_SetDoubleTagStorage(&LocalPid, TagName.c_str(), &Len, &EntType,
                                   const_cast<double *>(Data));
   checkMoabErr(Err, "iMOAB_SetDoubleTagStorage(" + TagName + ")");
}

// Builds the MOAB mesh for this ocean instance from the default
// Decomp/HorzMesh, scoped to owned cells only (matches MPAS-Ocean's
// init_moab_mpas, which never registers halo elements or ghost layers for
// its component-side MOAB application).
void createMesh(int LocalPid, int OcnID) {

   Decomp *DefDecomp   = Decomp::getDefault();
   HorzMesh *DefMesh   = HorzMesh::getDefault();
   const I4 NCellsOwned = DefDecomp->NCellsOwned;
   const I4 MaxEdges    = DefDecomp->MaxEdges;

   deepCopy(DefDecomp->CellIDH, DefDecomp->CellID);
   deepCopy(DefDecomp->VertexIDH, DefDecomp->VertexID);
   deepCopy(DefDecomp->VerticesOnCellH, DefDecomp->VerticesOnCell);
   deepCopy(DefDecomp->NEdgesOnCellH, DefDecomp->NEdgesOnCell);
   deepCopy(DefMesh->XVertexH, DefMesh->XVertex);
   deepCopy(DefMesh->YVertexH, DefMesh->YVertex);
   deepCopy(DefMesh->ZVertexH, DefMesh->ZVertex);

   // Single-pass compaction: remap Decomp's local vertex indices (which
   // range over owned+halo vertices) to a compact, first-encountered
   // numbering covering only the vertices actually referenced by owned
   // cells, and build the (NCellsOwned, MaxEdges) connectivity array in
   // that compact numbering. Cells with fewer than MaxEdges sides are
   // padded by repeating the cell's last real vertex, matching
   // mpas_moabmesh.F's convention.
   std::unordered_map<I4, I4> CompactIndex;
   std::vector<I4> CompactToLocal;
   std::vector<int> Connectivity(static_cast<size_t>(NCellsOwned) * MaxEdges);

   I4 Offset = 0;
   for (I4 Cell = 0; Cell < NCellsOwned; ++Cell) {
      const I4 NEdges = DefDecomp->NEdgesOnCellH(Cell);
      I4 LastCompact  = -1;
      for (I4 Edge = 0; Edge < NEdges; ++Edge) {
         const I4 LocalVertex = DefDecomp->VerticesOnCellH(Cell, Edge);
         auto It              = CompactIndex.find(LocalVertex);
         I4 Compact;
         if (It == CompactIndex.end()) {
            Compact = static_cast<I4>(CompactToLocal.size());
            CompactIndex.emplace(LocalVertex, Compact);
            CompactToLocal.push_back(LocalVertex);
         } else {
            Compact = It->second;
         }
         // iMOAB_CreateElements indexes into the just-created vertex range
         // as connectivity[j] + firstVertex - 1, i.e. it expects 1-based
         // indices; Compact itself stays 0-based everywhere else (it also
         // indexes Coords/CompactToLocal directly).
         Connectivity[Offset + Edge] = Compact + 1;
         LastCompact                 = Compact;
      }
      for (I4 Edge = NEdges; Edge < MaxEdges; ++Edge)
         Connectivity[Offset + Edge] = LastCompact + 1;
      Offset += MaxEdges;
   }

   int NCompactVerts = static_cast<int>(CompactToLocal.size());

   // Vertex coordinates, unit-sphere Cartesian, in compact order.
   std::vector<double> Coords(3 * static_cast<size_t>(NCompactVerts));
   for (int V = 0; V < NCompactVerts; ++V) {
      const I4 LocalVertex = CompactToLocal[V];
      Coords[3 * V + 0] = DefMesh->XVertexH[LocalVertex] / DefMesh->SphereRadius;
      Coords[3 * V + 1] = DefMesh->YVertexH[LocalVertex] / DefMesh->SphereRadius;
      Coords[3 * V + 2] = DefMesh->ZVertexH[LocalVertex] / DefMesh->SphereRadius;
   }

   ErrCode Err;
   int CoordsLen = 3 * NCompactVerts;
   int Dim       = 3;
   Err = iMOAB_CreateVertices(&LocalPid, &CoordsLen, &Dim, Coords.data());
   checkMoabErr(Err, "iMOAB_CreateVertices");

   int NElem          = static_cast<int>(NCellsOwned);
   int MBType         = 4; // MBPOLYGON
   int NNodesPerElem  = static_cast<int>(MaxEdges);
   int BlockID = 100 * OcnID + MachEnv::getDefault()->getMyTask();
   Err = iMOAB_CreateElements(&LocalPid, &NElem, &MBType, &NNodesPerElem,
                              Connectivity.data(), &BlockID);
   checkMoabErr(Err, "iMOAB_CreateElements");

   // GLOBAL_ID tags on vertices and elements, needed both for
   // ResolveSharedEntities and so the coupler-side offline weight file's
   // row/column numbering lines up with this mesh.
   std::string TagName = "GLOBAL_ID";
   int TagType          = DENSE_INTEGER;
   int NumCo            = 1;
   int TagIndex;
   Err = iMOAB_DefineTagStorage(&LocalPid, TagName.c_str(), &TagType, &NumCo,
                                &TagIndex);
   checkMoabErr(Err, "iMOAB_DefineTagStorage(GLOBAL_ID)");

   std::vector<int> VertexGlobalIDs(NCompactVerts);
   for (int V = 0; V < NCompactVerts; ++V)
      VertexGlobalIDs[V] = DefDecomp->VertexIDH[CompactToLocal[V]];

   int VertEntType = 0; // vertices
   Err = iMOAB_SetIntTagStorage(&LocalPid, TagName.c_str(), &NCompactVerts,
                                &VertEntType, VertexGlobalIDs.data());
   checkMoabErr(Err, "iMOAB_SetIntTagStorage(GLOBAL_ID, vertices)");

   std::vector<int> CellGlobalIDs(NElem);
   for (I4 Cell = 0; Cell < NCellsOwned; ++Cell)
      CellGlobalIDs[Cell] = DefDecomp->CellIDH[Cell];

   int ElemEntType = 1; // elements
   Err = iMOAB_SetIntTagStorage(&LocalPid, TagName.c_str(), &NElem,
                                &ElemEntType, CellGlobalIDs.data());
   checkMoabErr(Err, "iMOAB_SetIntTagStorage(GLOBAL_ID, elements)");

   Err = iMOAB_ResolveSharedEntities(&LocalPid, &NCompactVerts,
                                     VertexGlobalIDs.data());
   checkMoabErr(Err, "iMOAB_ResolveSharedEntities");

   Err = iMOAB_UpdateMeshInfo(&LocalPid);
   checkMoabErr(Err, "iMOAB_UpdateMeshInfo");
}

// Defines and sets the domain tags (lon/lat/area/aream/mask/frac) that the
// coupler-side mapper reads for the offline weight file's area-weighted
// interpolation. Values are computed the same way as the MCT path's
// ocn_set_domain_mct (via omega_get_lonlat_cell/omega_get_area_cell in
// omega_cxx2f_interface.cpp) so both drivers see identical numbers. aream
// is given a placeholder, same as ocn_set_domain_mct does for MCT: the
// coupler computes the real value from the mapping file and pushes it back
// down onto this tag (component_init_areacor_moab's 'x2c' exchange), but
// the tag must already exist here for that exchange to have somewhere to
// write it, and for its combined "area:aream:mask" read to succeed.
void setDomainTags(int LocalPid) {

   HorzMesh *DefMesh    = HorzMesh::getDefault();
   const I4 NCellsOwned = DefMesh->NCellsOwned;

   deepCopy(DefMesh->LonCellH, DefMesh->LonCell);
   deepCopy(DefMesh->LatCellH, DefMesh->LatCell);
   deepCopy(DefMesh->AreaCellH, DefMesh->AreaCell);

   const Real SphereRadius2 = DefMesh->SphereRadius * DefMesh->SphereRadius;

   std::vector<double> Lon(NCellsOwned), Lat(NCellsOwned), Area(NCellsOwned);
   std::vector<double> MaskFrac(NCellsOwned, 1.0);
   std::vector<double> AreaM(NCellsOwned, -9999.0);
   for (I4 Cell = 0; Cell < NCellsOwned; ++Cell) {
      Lon[Cell]  = static_cast<double>(DefMesh->LonCellH[Cell] * Rad2Deg);
      Lat[Cell]  = static_cast<double>(DefMesh->LatCellH[Cell] * Rad2Deg);
      Area[Cell] = static_cast<double>(DefMesh->AreaCellH[Cell] / SphereRadius2);
   }

   defineAndSetElemDoubleTag(LocalPid, "lon", NCellsOwned, Lon.data());
   defineAndSetElemDoubleTag(LocalPid, "lat", NCellsOwned, Lat.data());
   defineAndSetElemDoubleTag(LocalPid, "area", NCellsOwned, Area.data());
   defineAndSetElemDoubleTag(LocalPid, "aream", NCellsOwned, AreaM.data());
   defineAndSetElemDoubleTag(LocalPid, "mask", NCellsOwned, MaskFrac.data());
   defineAndSetElemDoubleTag(LocalPid, "frac", NCellsOwned, MaskFrac.data());
}

} // namespace

int moabInit(MPI_Comm Comm, int OcnID) {
   Pid = registerApplication(Comm, OcnID);
   createMesh(Pid, OcnID);
   setDomainTags(Pid);

   return Pid;
}

void moabDefineTagStorage(int LocalPid, const std::string &Cpl2OcnFieldNames,
                           const std::string &Ocn2CplFieldNames) {

   X2oTagNames = Cpl2OcnFieldNames;
   O2xTagNames = Ocn2CplFieldNames;

   ErrCode Err;
   int TagType = DENSE_DOUBLE;
   int NumCo   = 1;
   int TagIndex;

   Err = iMOAB_DefineTagStorage(&LocalPid, X2oTagNames.c_str(), &TagType,
                                &NumCo, &TagIndex);
   checkMoabErr(Err, "iMOAB_DefineTagStorage(x2o fields)");

   Err = iMOAB_DefineTagStorage(&LocalPid, O2xTagNames.c_str(), &TagType,
                                &NumCo, &TagIndex);
   checkMoabErr(Err, "iMOAB_DefineTagStorage(o2x fields)");
}

void moabImportTagStorage(Real *Buffer, int Len) {
   getDoubleTagStorage(X2oTagNames, Buffer, Len);
}

void moabExportTagStorage(const Real *Buffer, int Len) {
   setDoubleTagStorage(O2xTagNames, Buffer, Len);
}

} // namespace OMEGA
