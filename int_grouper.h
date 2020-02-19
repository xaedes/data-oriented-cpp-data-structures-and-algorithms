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

#include <vector>

// Group a list of integers by its values in O(N) time and O(N+M) memory.
//
// Memory usage is linear to number of items plus maximum integer value / aka
// group number. The resulting item order is stored in groupedIndices, which
// contains indices into the given list of group numbers. The position of the
// groups inside this array is stored in groupStarts and groupSizes.
// groupPointers is only internally used and you can ignore it.
class IntGrouper
{
public:
    IntGrouper()
    {
    }
    void group(
        const std::vector<int>& groupNumbers
    )
    {
        int numItems = groupNumbers.size();
        if (numItems == 0) 
        {
            group(groupNumbers, numItems, 0);
            return;
        }
        int numGroups = groupNumbers[0];
        for (int i = 1; i < numItems; ++i)
            if (groupNumbers[i] > numGroups) 
                numGroups = groupNumbers[i];

        group(groupNumbers, numItems, numGroups);

    }
    void group(
        const std::vector<int>& groupNumbers,
        int numGroups
    )
    {
        int numItems = groupNumbers.size();
        group(groupNumbers, numItems, numGroups);
    }
    void group(
        const std::vector<int>& groupNumbers,
        int numItems,
        int numGroups
    )
    {
        // initialize group info
        groupStarts.resize(numGroups);
        groupSizes.resize(numGroups);
        groupPointers.resize(numGroups);
        
        // when number of groups is not big, we can zero them all another
        // implementation strategy would be to only zero relevant groups. to
        // ensure consistency with remaining group data, each group is assigned
        // a sequence number, which can be used to check whether a group is
        // fresh or old and should therefore be treated as zero
        for (int k=0; k<numGroups; ++k)
        {
            groupStarts[k] = 0;
            groupSizes[k] = 0;
            groupPointers[k] = 0;
        }
        // determine group sizes
        for (int i=0; i<numItems; ++i)
        {
            int groupNumber = groupNumbers[i];
            ++groupSizes[groupNumber];
        }
        // determine group starts
        int nextGroupStart = 0;
        for (int k=0; k<numGroups; ++k)
        {
            groupStarts[k] = nextGroupStart;
            nextGroupStart += groupSizes[k];
        }
        assert(nextGroupStart == numItems);
        // group items
        groupedIndices.resize(numItems);
        for (int i = 0; i < numItems; ++i)
        {
            int groupNumber = groupNumbers[i];
            int insertIndex = groupStarts[groupNumber] + groupPointers[groupNumber];
            groupedIndices[insertIndex] = i;
            ++groupPointers[groupNumber];
        }
    }

    std::vector<int> groupedIndices;
    std::vector<int> groupStarts;
    std::vector<int> groupSizes;
    std::vector<int> groupPointers;
};
