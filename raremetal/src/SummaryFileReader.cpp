#include "SummaryFileReader.h"
#include "Error.h"
#include "Meta.h"

int SummaryFileReader::counter = 0;

//SummaryFileReader::SummaryFileReader() {}

bool SummaryFileReader::ReadTabix(String filename)
{
    myFilePtr = ifopen(filename, "r");
    old_pos = 0;
    old_chr = "";
    if (myFilePtr == NULL)
    {
        return (false);
    }
    String tabixFilename;
    tabixFilename = filename + ".tbi";
    StatGenStatus::Status tabixStat = myTabix.readIndex(tabixFilename.c_str());
    if (tabixStat != StatGenStatus::SUCCESS)
    {
        return (false);
    }
    return (true);
}

bool SummaryFileReader::ReadRecord(String chr, int pos)
{
    //buffer = "";
    myFilePtr->disableBuffering();
    if (old_chr == chr && old_pos == pos)
    {
        return (true);
    }
    //if current file pointer is close enough to the destination
    StringArray tmp;
    if (old_chr == chr && old_pos > pos - 16000 && old_pos < pos)
    {
        while (!ifeof(myFilePtr))
        {
            buffer.ReadLine(myFilePtr);
            if (buffer.Length() == 0)
            {
                return (false);
            }
            tmp.Clear();
//printf("old_chr=%s, chr=%s, old_pos=%d, pos=%d, buffer=%s\n",old_chr.c_str(), chr.c_str(), old_pos, pos, buffer.c_str());
            tmp.AddTokens(buffer, "\t");
            counter++;
            markerPosHash.SetInteger(tmp[0] + tmp[1], counter);
            markerNearby.Push(tmp[marker_col]);
            markerNearbyCov.Push(tmp[cov_col]);
            old_pos = tmp[1].AsInteger();
            old_chr = tmp[0];
            if (old_chr == chr && old_pos == pos)
            {
                marker_nearby = tmp[marker_col];
                marker_cov = tmp[cov_col];
                return (true);
            }
            if (old_pos > pos)
            {
                return (false);
            }
        }
    }

    uint64_t startPos = 0;
    // Find where this section starts in the file.
    if (!myTabix.getStartPos(chr.c_str(), pos, startPos))
    {
        // Didn't find the position.
        return (false);
    }
    if (startPos != (uint64_t) iftell(myFilePtr))
    {
        // Seek to the start position.
        if (ifseek(myFilePtr, startPos, SEEK_SET) != true)
        {
            return (false);
        }
    }

    if (ifeof(myFilePtr))
    {
        return (false);
    }

    while (!ifeof(myFilePtr))
    {
        buffer.ReadLine(myFilePtr);
        if (buffer.Length() == 0)
        {
            return (false);
        }
        tmp.Clear();
        tmp.AddTokens(buffer, "\t");
        counter++;
        markerPosHash.SetInteger(tmp[0] + tmp[1], counter);
        markerNearby.Push(tmp[marker_col]);
        markerNearbyCov.Push(tmp[cov_col]);
        old_pos = tmp[1].AsInteger();
        old_chr = tmp[0];
        if (old_chr == chr && old_pos == pos)
        {
            marker_nearby = tmp[marker_col];
            marker_cov = tmp[cov_col];
            return (true);
        }
        if (old_pos > pos)
        {
            return (false);
        }
    }

    if (ifeof(myFilePtr))
    {
        return (false);
    }
    return (true);
}
