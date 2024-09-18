/*!
 * @file        grid_polyhedra.h
 * @brief       inline functions of class PolyGrid
 *
 * @author      Wisces 〔wwwangqs17@163.com〕
 * @version     1.0
 * @date        2023-05-21
 * @copyright   Copyright (c) 2023  PERSONAL
 *
 * @par 修改日志:
 * <table>
 * <tr><th> Date        <th> Version  <th> Author  <th> Description
 * <tr><td> 2023-05-21  <td> 1.0      <td> Wisces  <td> 〔内容〕
 * </table>
 */

#ifdef TDTREE
inline void PolyGrid::SetTDTree(TDTree *in)
{
    sSet(TDTreeRoot, in);
}

inline TDTree *PolyGrid::GetTDTree() const
{
    return TDTreeRoot;
}
#endif
#ifdef GASPI

inline pair<int,int> **PolyGrid::GetC2B()
{
    return C2B;
}

inline pair<int,int> **PolyGrid::GetN2B()
{
    return N2B;
}

inline IntType *PolyGrid::GetnC2B()
{
    return nC2B;
}

inline IntType *PolyGrid::GetnN2B()
{
    return nN2B;
}

/*inline vector<pair<int,int>> *PolyGrid::Getf2b()
{
    return f2b;
}
*/
inline IntType PolyGrid::GetblockSize()
{
    return GASPI_blockSize;
}

inline IntType **PolyGrid::GetGhost()
{
    return Ghost;
}
#endif
//!< flow variables, dt_timestep, residuals, DQ - gField
//!< add by Wisces
inline void PolyGrid::SetLUSGSLayer(IntType *layer)
{
    sSet(LUSGSLayer, layer);
}

inline void PolyGrid::SetLUSGSCellOrder(IntType *luorder)
{
    sSet(LUSGSCellOrder, luorder);
}

inline void PolyGrid::SetLUSGScellsPerlayer(IntType *cellsPerlayer)
{
    sSet(LUSGScellsPerlayer, cellsPerlayer);
}

inline void PolyGrid::SetRho(RealFlow *q)
{
    sSet(rho, q);
}

inline void PolyGrid::SetU(RealFlow *q)
{
    sSet(u, q);
}

inline void PolyGrid::SetV(RealFlow *q)
{
    sSet(v, q);
}

inline void PolyGrid::SetW(RealFlow *q)
{
    sSet(w, q);
}

inline void PolyGrid::SetP(RealFlow *q)
{
    sSet(p, q);
}

inline void PolyGrid::SetDqdx(RealFlow *dq)
{
    sSet(dqdx, dq);
}

inline void PolyGrid::SetDqdy(RealFlow *dq)
{
    sSet(dqdy, dq);
}

inline void PolyGrid::SetDqdz(RealFlow *dq)
{
    sSet(dqdz, dq);
}

inline void PolyGrid::SetDt(RealFlow *q)
{
    sSet(dt, q);
}

inline void PolyGrid::SetRes(RealFlow *q)
{
    sSet(res, q);
}

inline void PolyGrid::SetDQ(RealFlow *q)
{
    sSet(DQ, q);
}

inline IntType *PolyGrid::GetLUSGSLayer() const
{
    return LUSGSLayer;
}

inline IntType *PolyGrid::GetLUSGSCellOrder() const
{
    return LUSGSCellOrder;
}

inline IntType *PolyGrid::GetLUSGScellsPerlayer() const
{
    return LUSGScellsPerlayer;
}

inline RealFlow *PolyGrid::GetRho() const
{
    return rho;
}

inline RealFlow *PolyGrid::GetU() const
{
    return u;
}

inline RealFlow *PolyGrid::GetV() const
{
    return v;
}

inline RealFlow *PolyGrid::GetW() const
{
    return w;
}

inline RealFlow *PolyGrid::GetP() const
{
    return p;
}

inline RealFlow *PolyGrid::GetDqdx() const
{
    return dqdx;
}

inline RealFlow *PolyGrid::GetDqdy() const
{
    return dqdy;
}

inline RealFlow *PolyGrid::GetDqdz() const
{
    return dqdz;
}

inline RealFlow *PolyGrid::GetDt() const
{
    return dt;
}

inline RealFlow *PolyGrid::GetRes() const
{
    return res;
}

inline RealFlow *PolyGrid::GetDQ() const
{
    return DQ;
}

//!< Some parameters for the Gradient computations
//!< add by Wisces
inline void PolyGrid::SetNodeSymm(IntType *symm)
{
    sSet(nodesymm, symm);
}

inline void PolyGrid::SetXfnNSymm(RealGeom *symm)
{
    sSet(xfn_n_symm, symm);
}

inline void PolyGrid::SetYfnNSymm(RealGeom *symm)
{
    sSet(yfn_n_symm, symm);
}

inline void PolyGrid::SetZfnNSymm(RealGeom *symm)
{
    sSet(zfn_n_symm, symm);
}

inline IntType *PolyGrid::GetNodeSymm() const
{
    return nodesymm;
}

inline RealGeom *PolyGrid::GetXfnNSymm() const
{
    return xfn_n_symm;
}

inline RealGeom *PolyGrid::GetYfnNSymm() const
{
    return yfn_n_symm;
}

inline RealGeom *PolyGrid::GetZfnNSymm() const
{
    return zfn_n_symm;
}

inline void PolyGrid::SetNTFace(const IntType in)
{
    nTFace = in;
}

inline void PolyGrid::SetNTCell(const IntType in)
{
    nTCell = in;
}

inline void PolyGrid::SetNBFace(const IntType in)
{
    nBFace = in;
}

inline void PolyGrid::SetNIFace(const IntType in)
{
    nIFace = in;
}

inline void PolyGrid::SetNINode(const IntType in)
{
    nINode = in;
}

inline IntType PolyGrid::GetNTFace() const
{
    return nTFace;
}

inline IntType PolyGrid::GetNTCell() const
{
    return nTCell;
}

inline IntType PolyGrid::GetNBFace() const
{
    return nBFace;
}

inline IntType PolyGrid::GetNIFace() const
{
    return nIFace;
}

inline IntType PolyGrid::GetNINode() const
{
    return nINode;
}

inline void PolyGrid::SetnNPF(IntType *in)
{
    sSet(nNPF, in);
}

inline void PolyGrid::SetnFPC(IntType *in)
{
    sSet(nFPC, in);
}

inline void PolyGrid::SetnNPC(IntType *in)
{
    sSet(nNPC, in);
}

inline void PolyGrid::Setf2n(IntType *in)
{
    sSet(f2n, in);
}

inline void PolyGrid::Setf2c(IntType *in)
{
    sSet(f2c, in);
}

inline RealGeom *PolyGrid::GetXcc() const
{
    return xcc;
}

inline RealGeom *PolyGrid::GetYcc() const
{
    return ycc;
}

inline RealGeom *PolyGrid::GetZcc() const
{
    return zcc;
}

// Get and set cell volume
inline RealGeom *PolyGrid::GetCellVol() const
{
    return vol;
}

inline RealGeom *PolyGrid::GetXfc() const
{
    return xfc;
}

inline RealGeom *PolyGrid::GetYfc() const
{
    return yfc;
}

inline RealGeom *PolyGrid::GetZfc() const
{
    return zfc;
}

inline RealGeom *PolyGrid::GetXfn() const
{
    return xfn;
}

inline RealGeom *PolyGrid::GetYfn() const
{
    return yfn;
}

inline RealGeom *PolyGrid::GetZfn() const
{
    return zfn;
}

// Get and set face aera
inline RealGeom *PolyGrid::GetFaceArea() const
{
    return area;
}

// inline void PolyGrid::SetcGrid(PolyGrid *grid)
// {
//     cGrid = grid;
// }

// inline void PolyGrid::SetfGrid(PolyGrid *grid)
// {
//     fGrid = grid;
// }

// inline Grid *PolyGrid::GetcGrid() const
// {
//     return cGrid;
// }

// inline Grid *PolyGrid::GetfGrid() const
// {
//     return fGrid;
// }

inline void PolyGrid::Setc2cc(IntType *in)
{
    sSet(c2cc, in);
}

inline IntType *PolyGrid::Getc2cc() const
{
    return c2cc;
}

inline void PolyGrid::SetnCPC(IntType *in)
{
    sSet(nCPC, in);
}

inline void PolyGrid::Setc2c(IntType **in)
{
    sSet(c2c, in);
}

inline void PolyGrid::SetC2N(IntType **in)
{
    sSet(C2N, in);
}

inline void PolyGrid::SetC2F(IntType **in)
{
    sSet(C2F, in);
}

//!< F2N is special and just a reference to f2n, so use 1D array operator, tangj
inline void PolyGrid::SetF2N(IntType **in)
{
    if ((in != NULL) && (F2N != in))
    {
        sdel_array_1D(F2N);
    }
    F2N = in;
}

inline IntType *PolyGrid::GetnNPF() const
{
    return nNPF;
}

inline IntType *PolyGrid::GetnFPC() const
{
    return nFPC;
}

inline IntType *PolyGrid::GetnNPC() const
{
    return nNPC;
}

inline IntType *PolyGrid::Getf2n() const
{
    return f2n;
}

inline IntType *PolyGrid::Getf2c() const
{
    return f2c;
}

inline IntType *PolyGrid::GetnCPC() const
{
    return nCPC;
}

inline IntType **PolyGrid::Getc2c() const
{
    return c2c;
}

inline IntType **PolyGrid::GetC2N() const
{
    return C2N;
}

inline IntType **PolyGrid::GetC2F() const
{
    return C2F;
}

inline IntType **PolyGrid::GetF2N() const
{
    return F2N;
}

// inline void PolyGrid::SetnCCPN(IntType *in)
// {
//     sSet(nCCPN, in);
// }

// inline IntType *PolyGrid::GetnCCPN() const
// {
//     return nCCPN;
// }

// inline void PolyGrid::SetN2CC(IntType **in)
// {
//     sSet(N2CC, in);
// }

// inline IntType **PolyGrid::GetN2CC() const
// {
//     return N2CC;
// }

// inline void PolyGrid::SetnVPN(IntType *in)
// {
//     sSet(nVPN, in);
// }

// inline IntType *PolyGrid::GetnVPN() const
// {
//     return nVPN;
// }

// inline void PolyGrid::SetN2V(IntType **in)
// {
//     sSet(N2V, in);
// }

// inline IntType **PolyGrid::GetN2V() const
// {
//     return N2V;
// }

// inline void PolyGrid::SetWeightNodeProl(RealGeom *weight_prol)
// {
//     sSet(WeightNodeProl, weight_prol);
// }

// inline RealGeom *PolyGrid::GetWeightNodeProl(void) const
// {
//     return WeightNodeProl;
// }

inline void PolyGrid::SetNumberOfFaceNeighbors(const IntType n)
{
    nNeighbor = n;
}

inline IntType PolyGrid::GetNumberOfFaceNeighbors() const
{
    return nNeighbor;
}

inline void PolyGrid::SetNumberOfNodeNeighbors(const IntType n)
{
    nNeighborN = n;
}

inline IntType PolyGrid::GetNumberOfNodeNeighbors() const
{
    return nNeighborN;
}

inline void PolyGrid::SetFaceNeighborZones(IntType *fnz)
{
    sSet(nb, fnz);
}

inline IntType *PolyGrid::GetFaceNeighborZones() const
{
    return nb;
}

inline void PolyGrid::SetNodeNeighborZones(IntType *nnz)
{
    sSet(nbN, nnz);
}

// inline IntType *PolyGrid::GetNodeNeighborZones() const
// {
//     return nbN;
// }

// inline void PolyGrid::SetNeighborGrids(PolyGrid **grids)
// {
//     //!< This is an object array, we only delete the array and do not delete the objects passed in the array.
//     if ((nbg != NULL) && (nbg != grids))
//     {
//         sdel_array_1D(nbg);
//     }
//     nbg = grids;
// }

inline void PolyGrid::SetnbZ(IntType *in)
{
    sSet(nbZ, in);
}

inline IntType *PolyGrid::GetnbZ() const
{
    return nbZ;
}

inline void PolyGrid::SetnbBF(IntType *in)
{
    sSet(nbBF, in);
}

inline IntType *PolyGrid::GetnbBF() const
{
    return nbBF;
}

inline void PolyGrid::SetnbSN(IntType *in)
{
    sSet(nbSN, in);
}

inline IntType *PolyGrid::GetnbSN() const
{
    return nbSN;
}

inline void PolyGrid::SetnbZN(IntType *in)
{
    sSet(nbZN, in);
}

inline IntType *PolyGrid::GetnbZN() const
{
    return nbZN;
}

inline void PolyGrid::SetnbRN(IntType *in)
{
    sSet(nbRN, in);
}

inline IntType *PolyGrid::GetnbRN() const
{
    return nbRN;
}

// inline IntType PolyGrid::GetLevel() const
// {
//     return level;
// }

// inline void PolyGrid::SetLevel(IntType lin)
// {
//     level = lin;
// }

inline void PolyGrid::Setbcr(BCRecord **bcrs)
{
    //!< This is an object array, we only delete the array and do not delete the objects passed in the array.
    if ((bcr != NULL) && (bcr != bcrs))
    {
        sdel_array_1D(bcr);
    }
    bcr = bcrs;
}

inline BCRecord **PolyGrid::Getbcr() const
{
    return bcr;
}

// inline void PolyGrid::SetVolAvg(RealGeom in)
// {
//     VolAvg = in;
// }

// inline RealGeom PolyGrid::GetVolAvg() const
// {
//     return VolAvg;
// }

// //---------------------------------------------------------
// inline RealGeom *PolyGrid::GetGridQualityFaceCentroidSkewness(void) const
// {
//     return facecentroidskewness;
// }

// inline IntType *PolyGrid::GetGridQualityCellWallNumber(void) const
// {
//     return cellwallnumber;
// }

//---------------------------------------------------------
// Flow reconstruction from cell to node
//---------------------------------------------------------
inline void PolyGrid::SetNodeType(IntType *node_type)
{
    sSet(Nmark, node_type);
}

inline IntType *PolyGrid::GetNodeType(void) const
{
    return Nmark;
}

inline void PolyGrid::SetWeightNodeDist(RealGeom *weight)
{
    sSet(WeightNodeDist, weight);
}

inline RealGeom *PolyGrid::GetWeightNodeDist(void) const
{
    return WeightNodeDist;
}

inline void PolyGrid::SetWeightNodeC2N(RealGeom **weight_c2n)
{
    sSet(WeightNodeC2N, weight_c2n);
}

inline RealGeom **PolyGrid::GetWeightNodeC2N(void) const
{
    return WeightNodeC2N;
}

inline void PolyGrid::SetWeightNodeN2C(RealGeom **weight_n2c)
{
    sSet(WeightNodeN2C, weight_n2c);
}

inline RealGeom **PolyGrid::GetWeightNodeN2C(void) const
{
    return WeightNodeN2C;
}

inline void PolyGrid::SetWeightNodeBFace2C(RealGeom **WeightNodebFace2c)
{
    sSet(WeightNodeBFace2C, WeightNodebFace2c);
}

inline RealGeom **PolyGrid::GetWeightNodeBFace2C(void) const
{
    return WeightNodeBFace2C;
}

//---------------------------------------------------------
// moving grid
//---------------------------------------------------------
inline RealGeom *PolyGrid::GetFaceNormalVelocity(void) const
{
    return vgn;
}

inline RealGeom *PolyGrid::GetBoundaryFaceVelocityX(void) const
{
    return BFacevgx;
}

inline RealGeom *PolyGrid::GetBoundaryFaceVelocityY(void) const
{
    return BFacevgy;
}

inline RealGeom *PolyGrid::GetBoundaryFaceVelocityZ(void) const
{
    return BFacevgz;
}

// // for line implicit
// inline void PolyGrid::SetLI_nLine(const IntType in)
// {
//     nLinesForLI = in;
// }

// inline IntType PolyGrid::GetLI_nLine() const
// {
//     return nLinesForLI;
// }

// inline void PolyGrid::SetLI_nCellsInLine(IntType *in)
// {
//     nCellsInLineForLI = in;
// }

// inline IntType *PolyGrid::GetLI_nCellsInLine() const
// {
//     return nCellsInLineForLI;
// }

// inline void PolyGrid::SetLI_CellOfLine(IntType **in)
// {
//     CellsOfLinesForLI = in;
// }

// inline IntType **PolyGrid::GetLI_CellOfLine() const
// {
//     return CellsOfLinesForLI;
// }

// inline void PolyGrid::SetLI_FaceOfLine(IntType **in)
// {
//     FacesOfLinesForLI = in;
// }

// inline IntType **PolyGrid::GetLI_FaceOfLine() const
// {
//     return FacesOfLinesForLI;
// }

// inline void PolyGrid::SetLI_LineOfCell(IntType *in)
// {
//     LineIndexOfCellsForLI = in;
// }

// inline IntType *PolyGrid::GetLI_LineOfCell() const
// {
//     return LineIndexOfCellsForLI;
// }

// inline void PolyGrid::SetLI_nLine_wall_start(const IntType in)
// {
//     n_lines_wall_start = in;
// }

// inline IntType PolyGrid::GetLI_nLine_wall_start() const
// {
//     return n_lines_wall_start;
// }

// inline void PolyGrid::SetLI_nLine_wall(const IntType in)
// {
//     n_lines_wall = in;
// }

// inline IntType PolyGrid::GetLI_nLine_wall() const
// {
//     return n_lines_wall;
// }

// inline void PolyGrid::SetLI_nCellinLine_wall(IntType *in)
// {
//     n_cells_in_line_wall = in;
// }

// inline IntType *PolyGrid::GetLI_nCellinLine_wall() const
// {
//     return n_cells_in_line_wall;
// }

// inline RealGeom &PolyGrid::span_length_of_2d()
// {
//     return span_len_of_2d;
// }

// inline RealGeom PolyGrid::span_length_of_2d() const
// {
//     return span_len_of_2d;
// }

// node to cell relationships:
inline IntType **PolyGrid::GetN2C() const
{
    return N2C;
}

inline IntType *PolyGrid::GetnCPN() const
{
    return nCPN;
}

inline void PolyGrid::SetN2C(IntType **in)
{
    sSet(N2C, in);
}

inline void PolyGrid::SetnCPN(IntType *in)
{
    sSet(nCPN, in);
}
