//
//  BKDT3DCartGridSBSolver.cpp
//  274F16NearestSB
//
//  Created by Yunlong Liu on 1/27/17.
//  Copyright © 2017 Yunlong Liu. All rights reserved.
//

#include "BKDT3DCartGridSBSolver.hpp"

BKDT3DCartGridSBSolver::BKDT3DCartGridSBSolver(double aveLocsPerCell) :
AVE_LOC_PER_CELL(aveLocsPerCell),
minX(DOUBLE_MAX), minY(DOUBLE_MAX), minZ(DOUBLE_MAX),
maxX(DOUBLE_MIN), maxY(DOUBLE_MIN), maxZ(DOUBLE_MIN) {}

size_t BKDT3DCartGridSBSolver::
getIdx(size_t gridX, size_t gridY, size_t gridZ) const {
    return gridX*gridYSize*gridZSize + gridY*gridZSize + gridZ;
}

size_t BKDT3DCartGridSBSolver::getIdx(double x, double y, double z) const {
    return getIdx(static_cast<size_t>((x-midX)/sideLen + gridXSize/2),
                  static_cast<size_t>((y-midY)/sideLen + gridYSize/2),
                  static_cast<size_t>((z-midZ)/sideLen + gridZSize/2));
}


void BKDT3DCartGridSBSolver::build() {
    std::vector<std::pair<Point<3>, const SBLoc*>> kdtData;
    for (auto const &l : *sbData) {
        auto pt = SBLoc::latLngToCart3DXYZ(l.lng, l.lat);
        minX = std::min(minX, pt[0]);
        minY = std::min(minY, pt[1]);
        minZ = std::min(minZ, pt[2]);
        maxX = std::max(maxX, pt[0]);
        maxY = std::max(maxY, pt[1]);
        maxZ = std::max(maxZ, pt[2]);
        kdtData.emplace_back(pt, &l);
    }
    

    sbKdt = KDTree<3, const SBLoc*, DistType::EUC>(kdtData.begin(), kdtData.end());
    std::cout << "Tree height is " << sbKdt.height() << std::endl;
    
    
    gridXSize = pow(sbData->size(), 1.0/3) / AVE_LOC_PER_CELL;
    sideLen = (maxX - minX)/gridXSize;
    gridXSize *= 2;
    gridYSize = (maxY - minY)/sideLen  +1;
    gridZSize = (maxZ - minZ)/sideLen + 1;
    gridCache = std::vector<KDTree<3, const SBLoc*, DistType::EUC>>(gridXSize*gridYSize*gridZSize);
    midX = (minX + maxX)/2, midY = (minY + maxY)/2, midZ = (minZ + maxZ)/2;
    size_t totalTreeSize = 0;
    std::vector<std::pair<Point<3>, const SBLoc*>> cacheCellLocs(sbData->size());
    
    for (size_t gridX = 0; gridX < gridXSize; ++gridX) {
        for (size_t gridY = 0; gridY < gridYSize; ++gridY) {
            for (size_t gridZ = 0; gridZ < gridZSize; ++gridZ) {
                Point<3> cellPt;
                cellPt[0] = midX + ((int)gridX-(int)gridXSize/2+0.5)*sideLen;
                cellPt[1] = midY + ((int)gridY-(int)gridYSize/2+0.5)*sideLen;
                cellPt[2] = midZ + ((int)gridZ-(int)gridZSize/2+0.5)*sideLen;
                auto end = sbKdt.rangeDiffKNNPairs(cellPt, sideLen*sqrt(3)/2,
                           cacheCellLocs.begin());
                totalTreeSize += end - cacheCellLocs.begin();
                gridCache[getIdx(gridX, gridY, gridZ)] = KDTree<3, const SBLoc*, DistType::EUC>(cacheCellLocs.begin(), end);
            }
        }
    }
    std::cout << "Average tree size is: "
              << totalTreeSize/(gridXSize*gridYSize*gridZSize) << std::endl;
}

const SBLoc* BKDT3DCartGridSBSolver::findNearest(double lng, double lat) const {
    auto pt = SBLoc::latLngToCart3DXYZ(lng, lat);
    return gridCache[getIdx(pt[0], pt[1], pt[2])].kNNValue(pt, 1);
}
