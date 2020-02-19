#pragma once
#include <vector>

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
