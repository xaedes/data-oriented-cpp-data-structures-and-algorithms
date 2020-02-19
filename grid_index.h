#pragma once

// MIT License

// Copyright (c) 2020 xaedes

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include <functional> 
#include <algorithm> 
#include <vector> 
#include <cmath> 

// Will use simple C++ only implementation of Vec or opencv if you wish so.
#include "vec.h"

/**
 * @brief      This class describes a grid by its bounds in the world and
 *             cell size.
 *
 *             It provides functions to convert coordinates between world
 *             coordinates, grid cell coordinates and one dimensional cell
 *             indices.
 */
template<typename TValueType, unsigned int TNumDimensions>
class GridMeta
{
public:
    typedef Vec<TValueType, TNumDimensions> ValueVec;
    typedef Vec<int, TNumDimensions> IntVec;
    typedef Vec<bool, TNumDimensions> BoolVec;
    GridMeta(ValueVec lowerBound, ValueVec upperBound, ValueVec cellSize, BoolVec wrap)
    {
        set(lowerBound, upperBound, cellSize, wrap);
    }

    void set(ValueVec lowerBound, ValueVec upperBound, ValueVec cellSize, BoolVec wrap)
    {
        this->m_lowerBound = lowerBound;
        this->m_upperBound = upperBound;
        this->m_cellSize = cellSize;
        this->wrap = wrap;
        update();
    }

    void setCellSize(ValueVec cellSize)
    {
        this->m_cellSize = cellSize;
    }

    IntVec gridCoord(int cellIndex)
    {
        IntVec result;
        result[0] = (int)(cellIndex);
        for (int i = 0; i < TNumDimensions; ++i)
        {
            result[i] = (cellIndex / m_steps[i]) % m_gridSize[i];
            cellIndex -= result[i] * m_steps[i];
        }
        assert(cellIndex == 0);
        return result;
    }

    IntVec gridCoord(ValueVec worldCoord)
    {
        // auto scaled_coord = (worldCoord - m_lowerBound) / m_cellSize;
        IntVec result;
        for (int i = 0; i < TNumDimensions; ++i)
        {

            result[i] = (int)((worldCoord[i] - m_lowerBound[i]) / m_cellSize[i]);
            result[i] = result[i] % m_gridSize[i];
            // % can give negative result (up to -(m_gridSize[i]-1))
            if (result[i] < 0) result[i] += m_gridSize[i];
        }
        return result;
    }

    int cellIndex(ValueVec worldCoord)
    {

        return cellIndex(gridCoord(worldCoord));
    }

    int cellIndex(IntVec gridCoord)
    {
        int index = 0;
        for (int i = 0; i < TNumDimensions; ++i)
        {
            index += m_steps[i] * gridCoord[i];
        }
        return index;
    }

    ValueVec worldCoord(int cellIndex)
    {

        return worldCoord(gridCoord(cellIndex));
    }
    ValueVec worldCoord(IntVec gridCoord)
    {
        ValueVec result;
        for (int i = 0; i < TNumDimensions; ++i)
            result[i] = gridCoord[i] * m_cellSize[i] + m_lowerBound[i];
        return result;
    }

    ValueVec lowerBound() { return m_lowerBound; }
    ValueVec upperBound() { return m_upperBound; }
    ValueVec cellSize() { return m_cellSize; }
    ValueVec gridDimensions() { return m_gridDimensions; }
    IntVec gridSize() { return m_gridSize; }
    IntVec steps() { return m_steps; }
    int numCells() { return m_numCells; }

    BoolVec wrap;

    bool bounds_check(int dimension, int& index)
    {
        if (wrap[dimension])
        {
            index = index % m_gridSize[dimension];
            if (index < 0) index += m_gridSize[dimension];
            return true;
        }
        else
        {
            if (index < 0) return false;
            if (index >= m_gridSize[dimension]) return false;
            return true;
        }
    }

    bool bounds_check(IntVec& gridCoord)
    {
        bool inBounds = true;
        for (int i = 0; i < TNumDimensions; ++i)
        {
            inBounds &= bounds_check(i, gridCoord[i]);
            if (!inBounds) break;
        }
        return inBounds;
    }

protected:
    void update()
    {
        for (int i = 0; i < TNumDimensions; ++i)
        {
            float lower = std::min(m_lowerBound[i], m_upperBound[i]);
            float upper = std::max(m_lowerBound[i], m_upperBound[i]);
            m_lowerBound[i] = lower;
            m_upperBound[i] = upper;
        }
        m_gridDimensions = m_upperBound - m_lowerBound;
        m_numCells = 1;
        for (int i = 0; i < TNumDimensions; ++i)
        {
            m_gridSize[i] = (int)ceil(m_gridDimensions[i] / m_cellSize[i]);
            if (m_gridSize[i] < 1) m_gridSize[i] = 1;
            m_steps[i] = m_numCells; // stride for dimension 0 is always 1
            m_numCells *= m_gridSize[i];
        }
    }


    ValueVec m_lowerBound;
    ValueVec m_upperBound;
    ValueVec m_cellSize;

    ValueVec m_gridDimensions;
    IntVec m_gridSize;
    IntVec m_steps; // steps between elements in kth dimension (m_steps[0] == 1)
    int m_numCells;
};

struct GridCell
{
    // start index into cell contents array
    int start;
    // cell capacity, i.e. how many items will be contained
    int capacity;
    // cell size, i.e. how many items are contained so far
    int size;
    // contains current seq number, is used to check the age of cell to
    // discard not freshly updated ones (or treat them as zero)
    int seqNum;

    GridCell() : start(0), capacity(0), size(0), seqNum(0) {}
    GridCell(int start, int capacity, int size, int seqNum)
        : start(start), capacity(capacity), size(size), seqNum(seqNum)
    {}
};

// https://stackoverflow.com/a/27270738/798588
template <int A, int B>
struct get_power
{
    static const int value = A * get_power<A, B - 1>::value;
};
template <int A>
struct get_power<A, 0>
{
    static const int value = 1;
};

template<typename TValueType, unsigned int TNumDimensions>
class GridIndex
{
public:
    typedef Vec<TValueType, TNumDimensions> ValueVec;
    typedef Vec<int, TNumDimensions> IntVec;
    typedef Vec<bool, TNumDimensions> BoolVec;

    const int NumNeighborhoodCells = get_power<3, TNumDimensions>::value;

    GridIndex(GridMeta<TValueType, TNumDimensions> gridMeta)
        : gridMeta(gridMeta), numItems(0), indexedCollection()
    {
    }
    GridIndex(GridMeta<TValueType, TNumDimensions> gridMeta, const std::vector<ValueVec>& items)
        : gridMeta(gridMeta), numItems(0), indexedCollection()
    {
        setItems(items);
    }

    void setItems(const std::vector<ValueVec>& items)
    {
        seqNum++;
        indexedCollection = items; // copy assignment

        numItems = indexedCollection.size();
        memoryAllocation();
        initializeRelevantCells();
        computeCellCapacities();
        insertItems();
    }

    std::vector<int> collectNeighbors(int item_index)
    {
        std::vector<int> result;
        int numNeighborCells = cellNeighborhoodSizes[item_index];
        int offset = item_index * NumNeighborhoodCells;
        for (int i = 0; i < numNeighborCells; ++i)
        {
            int neighborCellIndex = cellNeighborhoodContents[offset + i];
            GridCell cell = cells[neighborCellIndex];
            if (cell.seqNum < seqNum)
            {
                cell = GridCell();
            }
            for (int j = 0; j < cell.size; ++j)
            {
                int k = cellContents[cell.start + j];
                result.push_back(k);
            }
        }
        return result;
    }

    std::vector<int> getVectorCollectNeighbors(int item_index)
    {
        std::vector<int> results;
        collectNeighbors(item_index, results);
        return results;
    }
    std::vector<int> getVectorCollectNeighborsMaxDistance(int item_index, TValueType maxDistance)
    {
        std::vector<int> results;
        collectNeighborsMaxDistance(item_index, maxDistance, results);
        return results;
    }

    void collectNeighbors(int item_index, std::vector<int>& results)
    {
        int numNeighborCells = cellNeighborhoodSizes[item_index];
        int offset = item_index * NumNeighborhoodCells;
        for (int i = 0; i < numNeighborCells; ++i)
        {
            int neighborCellIndex = cellNeighborhoodContents[offset + i];
            GridCell cell = cells[neighborCellIndex];
            if (cell.seqNum < seqNum)
            {
                cell = GridCell();
            }
            for (int j = 0; j < cell.size; ++j)
            {
                int k = cellContents[cell.start + j];
                results.push_back(k);
            }
        }
    }

    void collectNeighborsMaxDistance(int item_index, TValueType maxDistance, std::vector<int>& results)
    {
        ValueVec item = indexedCollection[item_index];
        int numNeighborCells = cellNeighborhoodSizes[item_index];
        int offset = item_index * NumNeighborhoodCells;
        for (int i = 0; i < numNeighborCells; ++i)
        {
            int neighborCellIndex = cellNeighborhoodContents[offset + i];
            GridCell cell = cells[neighborCellIndex];
            if (cell.seqNum < seqNum)
            {
                cell = GridCell();
            }
            for (int j = 0; j < cell.size; ++j)
            {
                int k = cellContents[cell.start + j];
                ValueVec diff = item - indexedCollection[k];
                TValueType distance = norm(diff);
                if (distance <= maxDistance)
                {
                    results.push_back(k);
                }
            }
        }
    }

    void iterateNeighbors(int item_index, std::function<void(int)> callback)
    {
        int numNeighborCells = cellNeighborhoodSizes[item_index];
        int offset = item_index * NumNeighborhoodCells;
        for (int i = 0; i < numNeighborCells; ++i)
        {
            int neighborCellIndex = cellNeighborhoodContents[offset + i];
            GridCell cell = cells[neighborCellIndex];
            if (cell.seqNum < seqNum)
            {
                cell = GridCell();
            }
            for (int j = 0; j < cell.size; ++j)
            {
                int k = cellContents[cell.start + j];
                callback(k);
            }
        }
    }

    void iterateNeighborsMaxDistance(int item_index, TValueType maxDistance, std::function<void(int)> callback)
    {
        ValueVec item = indexedCollection[item_index];
        std::function<void(int)> guardedCallback = [this, &item, maxDistance, &callback](int other_index)
        {
            ValueVec diff = item - this->indexedCollection[other_index];
            TValueType distance = norm(diff);
            if (distance <= maxDistance)
            {
                callback(other_index);
            }
        };
        iterateNeighbors(item_index, guardedCallback);
    }

protected:
    void initializeRelevantCells()
    {
        for (int i = 0; i < numItems; ++i)
        {
            ValueVec item = indexedCollection[i];
            int cellIndex = gridMeta.cellIndex(item);
            cellIndices[i] = cellIndex;
            initializeCellNeighborhood(i);
        }
    }

    bool nextNeighbor(IntVec& delta)
    {
        // iterate all coordinates generated from [-1,0,+1]^TNumDimensions 

        // do trinary counting with digits out of [-1,0,+1] by incrementing
        // the lowest dimension and carry over to higher dimensions 

        // advance on first dimension
        delta[0] += 1;

        // handle carry over to higher dimensions
        int dim_index = 0;
        while (delta[dim_index] > 1)
        {
            // overflow happened on current dimension
            // set current dimension to lowest digit (that is -1)
            delta[dim_index] = -1;
            // and carry over to next dimension
            dim_index++;
            if (dim_index < TNumDimensions)
            {
                delta[dim_index] += 1;
            }
            else
            {
                // carry overflow
                // i.e. all neighbors where iterated
                return false;
            }
        }

        // there are more neighbors to iterate over
        return true;
    }

    void initializeCellNeighborhood(int item_index)
    {
        // whole [+1,0,-1]**NumDims neighborhood could be accessed, 
        // initialize these cells
        // also insert neighborhoods for later use
        cellNeighborhoodSizes[item_index] = 0;

        IntVec gridCoord = gridMeta.gridCoord(cellIndices[item_index]);

        int offset = item_index * NumNeighborhoodCells;

        // list all neighbors by trinary counting with TNumDimensions digits
        // out of [-1,0,+1]
        IntVec delta = -1;
        do
        {
            IntVec currentCoord = gridCoord + delta;
            if (gridMeta.bounds_check(currentCoord))
            {
                int currentCellIndex = gridMeta.cellIndex(currentCoord);
                // insert cell into neighborhood
                cellNeighborhoodContents[offset + cellNeighborhoodSizes[item_index]] = currentCellIndex;
                ++cellNeighborhoodSizes[item_index];
                // zero the cell and update its age to the current seqNum
                // when others find cells with old seqNum they should assume
                // it was not zeroed due too lazyness (improve performance)
                // and treat it as zero
                GridCell& cell = cells[currentCellIndex];
                if (cell.seqNum < seqNum)
                {
                    cell.start = 0;
                    cell.capacity = 0;
                    cell.size = 0;
                    cell.seqNum = seqNum;
                }
            }

        } while (nextNeighbor(delta));

    }

    void computeCellCapacities()
    {
        // all relevant cells are zeroed out
        for (int i = 0; i < numItems; ++i)
        {
            cells[cellIndices[i]].capacity += 1;
        }
    }

    void insertItems()
    {
        // all relevant cells are zeroed out and capacity is already
        // computed

        int nextContentIndex = 0;
        // insert item indices into cells. memory for cells is "allocated" as
        // memory blocks in cellContents on the fly
        for (int i = 0; i < numItems; ++i)
        {
            int k = cellIndices[i];
            GridCell& cell = cells[k];
            if (cell.size == 0)
            {
                // no item inserted in this cell so far, "allocate" memory for it
                cell.start = nextContentIndex;
                nextContentIndex += cell.capacity;
            }
            int insertIndex = cell.start + cell.size;
            cellContents[insertIndex] = i;
            ++cell.size;
        }
    }

    void memoryAllocation()
    {
        int minCapacity = 50;
        int numCells = gridMeta.numCells();
        bool cellArraysTooSmall = cells.size() < std::max(minCapacity, numCells);
        bool cellArraysTooBig = cells.size() > numCells;
        // just too be sure, but it isn't nice because it could result in infrequent periodical stutter
        // bool tooOld = seqNum > 1000;
        if (cellArraysTooSmall)
        {
            cells.resize(std::max(minCapacity, numCells));
            int sizeBeforeEnlarge = cells.size();
            std::fill(cells.begin()+ sizeBeforeEnlarge, cells.end(), GridCell());
        }
        else if (cellArraysTooBig)
        {
            cells.resize(std::min(minCapacity, numCells));
        }
        if (cellContents.size() < numItems)
        {
            cellContents.resize(numItems);
        }
        if (cellIndices.size() < numItems)
        {
            cellIndices.resize(numItems);
        }
        if (cellNeighborhoodContents.size() < numItems * NumNeighborhoodCells)
        {
            cellNeighborhoodContents.resize(numItems * NumNeighborhoodCells);
        }
        if (cellNeighborhoodSizes.size() < numItems)
        {
            cellNeighborhoodSizes.resize(numItems);
        }
    }

public:
    std::vector<ValueVec> indexedCollection;

    // contains all item indices grouped by cell in contigous blocks
    // the offsets of each cell can be retrieved via `cells`
    std::vector<int> cellContents;
    // contains cell numbers for each item of indexed collection
    std::vector<int> cellIndices;
    // stores for each cell the offsets into `cellContents`, the cell size &
    // capacity and its age
    std::vector<GridCell> cells;
    // contains neighborhood cell numbers grouped by item of indexed
    // collection in contigous blocks of 3**TNumDimensions
    std::vector<int> cellNeighborhoodContents;
    // contains actual number of neighborhood (including self) cells for
    // each item. due to world bounds this can be less than
    // 3**TNumDimensions
    std::vector<int> cellNeighborhoodSizes;

    GridMeta<TValueType, TNumDimensions> gridMeta;
    int numItems;
    int seqNum;
    friend class NeighborVisitor;
};


