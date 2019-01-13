//
//  UniCellBKDTGridSBSolver.cpp
//  274F16NearestSB
//
//  Created by nick on 12/30/18.
//  Copyright © 2018 Yunlong Liu. All rights reserved.
//

#include "UniCellBKDTGridSBSolver.hpp"

//#include <omp.h>




template <template <size_t, typename, typename Point<3>::DistType> class Tree>
UniCellBKDTGridSBSolver<Tree>::
UniCellBKDTGridSBSolver(double alpc) : BKDTSBSolver<Tree>(), AVE_LOC_PER_CELL(alpc) {}

template <template <size_t, typename, typename Point<3>::DistType> class Tree>
void UniCellBKDTGridSBSolver<Tree>::printSolverInfo() const {
    this->locKdt.printTreeInfo();
}


template <template <size_t, typename elem, typename Point<3>::DistType> class Tree>
void UniCellBKDTGridSBSolver<Tree>::fillGridCache(size_t rowSize, double sideLen) {
    gridCache.reserve(this->locKdt.size()*1.2/AVE_LOC_PER_CELL);
    thisRowStartIdx.reserve(rowSize);
    std::vector<std::pair<Point<3>, const SBLoc*>> ptLocPairs(this->locKdt.size());
    size_t totalTreeSize = 0, singleLocs = 0;
    //#pragma omp parallel for num_threads(std::thread::hardware_concurrency())\
    //default(none) schedule(guided) shared(diff) firstprivate(ptLocPairs) \
    //reduction(+:totalTreeSize, singleLocs) collapse(2)
    double thisLat = -0.5 * M_PI, thisCtrLat = thisLat + 0.5 * latInc;
    for (size_t r = 0, idx = 0; r < rowSize;
         ++r, thisCtrLat += latInc, thisLat += latInc) {
        size_t thisColSize = static_cast<size_t>(2*M_PI * SBLoc::EARTH_RADIUS *
                             cos(thisLat>0?thisLat-latInc:thisLat)/sideLen) + 2,
               thisEndIdx = idx + thisColSize;
        double thisLngInc = 2*M_PI/thisColSize +
                            2*M_PI/(thisColSize*thisColSize*0xFFFF),
               thisDiff = SBLoc::xyzDistFromLngLat(r*latInc - 0.5*M_PI,
                          (r+1)*latInc - 0.5*M_PI, thisLngInc),
               thisCtrLng = 0.5 * thisLngInc - M_PI;
        thisRowStartIdx.emplace_back(idx, thisLngInc);
        for (; idx < thisEndIdx; ++idx, thisCtrLng += thisLngInc) {
            auto locsEnd = this->locKdt.rangeDiffKNNPairs
                           (SBLoc::latLngToCart3DPt(thisCtrLng, thisCtrLat),
                           thisDiff, ptLocPairs.begin());
            size_t locsSize = locsEnd - ptLocPairs.begin();
            if (locsSize > 1) {
                gridCache.emplace_back(std::piecewise_construct,
                                       std::forward_as_tuple
                                       (ptLocPairs.begin(), locsEnd),
                                       std::forward_as_tuple(nullptr));
            } else {
                gridCache.emplace_back(std::piecewise_construct,
                                       std::forward_as_tuple(),
                                       std::forward_as_tuple(ptLocPairs[0].second));
                singleLocs++;
            }
            totalTreeSize += locsSize;
        }
    }
    std::cout << "Total tree nodes: " << totalTreeSize
              << "\nAve tree size: " << totalTreeSize/gridCache.size()
              << "\nAve tree height: "
              << static_cast<size_t>(log2(totalTreeSize/gridCache.size()+1)) + 1
              << "\nRatio of tree nodes over num locs: "
              << totalTreeSize/this->locKdt.size()
              << "\nSingle loc cells: " << singleLocs << "\nMulti-loc cells:"
              << gridCache.size() - singleLocs << std::endl;
}

template <template <size_t, typename elem, typename Point<3>::DistType> class Tree>
double UniCellBKDTGridSBSolver<Tree>::calcSideLenFromAlpc() {
    double surfaceArea = 4*M_PI*SBLoc::EARTH_RADIUS*SBLoc::EARTH_RADIUS;
    double numCells = this->locKdt.size()/AVE_LOC_PER_CELL;
    return sqrt(surfaceArea/numCells);
}

template <template <size_t, typename elem, typename Point<3>::DistType> class Tree>
void UniCellBKDTGridSBSolver<Tree>::
build(const std::shared_ptr<std::vector<SBLoc>> &locData) {
    BKDTSBSolver<Tree>::generateKDT(locData);
    double sideLen = calcSideLenFromAlpc();
    latInc = std::fabs(SBLoc::latFromHavDist(sideLen, 0));
    size_t rowSize = std::ceil(M_PI/(latInc - latInc*latInc/(M_PI*0xFFFF)));
    fillGridCache(rowSize, sideLen);
}

template <template <size_t, typename elem, typename Point<3>::DistType> class Tree>
const SBLoc* UniCellBKDTGridSBSolver<Tree>::
findNearest(double lng, double lat) const {
    auto[startIdx, thisLngInc] = thisRowStartIdx[(lat+0.5*M_PI)/latInc];
    const auto& [cacheTree, singleLoc] =
        gridCache[startIdx + static_cast<size_t>((lng+M_PI)/thisLngInc)];
    return singleLoc ? singleLoc :
           cacheTree.kNNValue(SBLoc::latLngToCart3DPt(lng, lat), 1);
}





template class UniCellBKDTGridSBSolver<KDTree>;
template class UniCellBKDTGridSBSolver<KDTreeCusMem>;
