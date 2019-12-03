using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace _20190618_GlycoTools_V2
{
    class ExtractFeatures
    {
        public static void GetApexPeak(LFPeptide targetPeptide, bool isLibrary = true)
        {
            List<RTPeak> XICPeaks = new List<RTPeak>();
            int startIndex = 0;
            if (isLibrary)
            {
                XICPeaks.AddRange(targetPeptide.SmoothLibrary);
                startIndex = XICPeaks.BinarySearch(new RTPeak(0, 0, targetPeptide.parentMS1Time));
                if (startIndex < 0)
                {
                    startIndex = ~startIndex;
                }
            }
            else
            {
                XICPeaks.AddRange(targetPeptide.SmoothExperiment);
                startIndex = XICPeaks.BinarySearch(targetPeptide.ExperimentLookupPeak);
                if (startIndex < 0)
                {
                    startIndex = ~startIndex;
                }
            }
            int leftIndex = 0;
            int apexIndex = startIndex;
            int rightIndex = 0;

            double sumLeft = 0;
            double sumRight = 0;

            for (int i = 1; i <= 5; i++)
            {
                if ((apexIndex - i) >= 0)
                {
                    sumLeft += XICPeaks[apexIndex - i].Intensity;
                }
                if ((apexIndex + i) < XICPeaks.Count)
                {
                    sumRight += XICPeaks[apexIndex + i].Intensity;
                }
            }
            if (sumLeft == 0 && sumRight == 0 && isLibrary)
            {
                targetPeptide.apexPeakLibrary = new RTPeak(targetPeptide.AdjustedMZ_Library, 0, targetPeptide.parentMS1Time);
                targetPeptide.rightPeakLibrary = new RTPeak(targetPeptide.AdjustedMZ_Library, 0, targetPeptide.parentMS1Time);
                targetPeptide.leftPeakLibrary = new RTPeak(targetPeptide.AdjustedMZ_Library, 0, targetPeptide.parentMS1Time);
                return;
            }
            if (sumLeft == 0 && sumRight == 0 && !isLibrary)
            {
                targetPeptide.apexPeakExperiment = new RTPeak(targetPeptide.AdjustedMZ_Experiment, 0, targetPeptide.parentMS1Time);
                targetPeptide.rightPeakExperiment = new RTPeak(targetPeptide.AdjustedMZ_Experiment, 0, targetPeptide.parentMS1Time);
                targetPeptide.leftPeakExperiment = new RTPeak(targetPeptide.AdjustedMZ_Experiment, 0, targetPeptide.parentMS1Time);
                return;
            }
            try
            {
                if (sumRight > sumLeft) //assume that the max lies to the right
                {
                    bool isAscending = true;
                    bool isDescending = false;
                    RTPeak currentMax = XICPeaks[startIndex];
                    double threshold = currentMax.Intensity * 0.01;
                    RTPeak currentMinRight = null;
                    RTPeak currentMinLeft = null;
                    int currentIndex = 1 + startIndex;
                    while (isAscending)
                    {
                        if (currentIndex < XICPeaks.Count)
                        {
                            var currentPeak = XICPeaks[currentIndex];
                            if (currentPeak.Intensity > currentMax.Intensity) //Reset max and threshold, continue on
                            {
                                currentMax = currentPeak;
                                threshold = currentMax.Intensity * 0.01;
                                currentIndex++;
                                continue;
                            }
                            else //Could be that I was at the apex, or I'm in a dip
                            {
                                if (currentPeak.Intensity < (currentMax.Intensity * (1.0 / 1.8)))
                                //I was at the max, so I want to keep this point
                                {
                                    //and start looking for the right minimum
                                    isAscending = false;
                                    isDescending = true;
                                    break;
                                }
                                else
                                //could be the start of the descent or some weird dip. Look ahead to the next peaks to find out.
                                {
                                    int AdditionalCounter = 1;
                                    while (true) //keep in this loop until one of the conditions is satisfied.
                                    {
                                        int aheadIndex = currentIndex + AdditionalCounter;
                                        if ((aheadIndex) < XICPeaks.Count)
                                        {

                                            var aheadPeak = XICPeaks[aheadIndex];
                                            if (aheadPeak.Intensity < (currentMax.Intensity * (1.0 / 1.8)))
                                            {
                                                isAscending = false;
                                                isDescending = true;
                                                break;
                                            }
                                            if (aheadPeak.Intensity > currentMax.Intensity)
                                            {
                                                currentMax = aheadPeak;
                                                currentIndex += AdditionalCounter;
                                                //currentIndex++;
                                                break;
                                            }
                                            AdditionalCounter++;
                                        }
                                        else
                                        {
                                            isAscending = false;
                                            isDescending = true;
                                            break;
                                        }
                                    }
                                }
                            }
                        }
                    }
                    //you have located the max. assume that the peak to right is the local minima
                    currentIndex++;
                    if (currentIndex < XICPeaks.Count)
                    {
                        currentMinRight = XICPeaks[currentIndex];
                    }
                    currentIndex++;
                    while (isDescending)
                    {
                        if (currentIndex < XICPeaks.Count)
                        {
                            var currentPeak = XICPeaks[currentIndex];
                            if (currentPeak.Intensity <= currentMinRight.Intensity)
                            {
                                currentMinRight = currentPeak;
                                if (currentPeak.Intensity < threshold)
                                {
                                    isDescending = false;
                                    break;
                                }
                            }
                            if (currentPeak.Intensity > currentMinRight.Intensity)
                            //could be the local min or some weird dip
                            {
                                int AdditionalCounter = 1;
                                while (true)
                                {
                                    int aheadIndex = currentIndex + AdditionalCounter;
                                    if ((aheadIndex) < XICPeaks.Count)
                                    {

                                        var aheadPeak = XICPeaks[aheadIndex];
                                        if (aheadPeak.Intensity > (1.8 * currentMinRight.Intensity))
                                        //There is a new peak coming up. You have found the local min. Break out.
                                        {
                                            isDescending = false;
                                            break;
                                        }
                                        if (aheadPeak.Intensity < currentMinRight.Intensity || aheadPeak.Intensity == 0)
                                        {
                                            currentMinRight = aheadPeak;
                                            currentIndex += AdditionalCounter;
                                            //currentIndex++;
                                            break;
                                        }
                                        AdditionalCounter++;
                                    }
                                    else
                                    {
                                        isDescending = false;
                                        break;
                                    }
                                }
                            }
                        }
                        else
                        {
                            isDescending = false;
                            break;
                        }
                        currentIndex++;
                    }
                    //you now have the max, and the right minima. Assume that the peak to the left of the start point is the left minima.
                    currentIndex = startIndex - 1;
                    if (currentIndex >= 0)
                    {
                        currentMinLeft = XICPeaks[currentIndex];
                    }
                    currentIndex--;
                    isDescending = true;
                    while (isDescending)
                    {
                        if (currentIndex >= 0)
                        {
                            var currentPeak = XICPeaks[currentIndex];
                            if (currentPeak.Intensity <= currentMinLeft.Intensity)
                            {
                                currentMinLeft = currentPeak;
                                if (currentPeak.Intensity < threshold)
                                {
                                    isDescending = false;
                                    break;
                                }
                            }
                            if (currentPeak.Intensity > currentMinLeft.Intensity)
                            //could be the local min or some weird dip
                            {
                                int AdditionalCounter = 1;
                                while (true)
                                {
                                    int aheadIndex = currentIndex - AdditionalCounter;
                                    if ((aheadIndex) >= 0)
                                    {

                                        var aheadPeak = XICPeaks[aheadIndex];
                                        if (aheadPeak.Intensity > (1.8 * currentMinLeft.Intensity))
                                        //There is a new peak coming up. You have found the local min. Break out.
                                        {
                                            isDescending = false;
                                            break;
                                        }
                                        if (aheadPeak.Intensity < currentMinLeft.Intensity || aheadPeak.Intensity == 0)
                                        {
                                            currentMinLeft = aheadPeak;
                                            currentIndex -= AdditionalCounter;
                                            //currentIndex--;
                                            break;
                                        }
                                        AdditionalCounter++;
                                    }
                                    else
                                    {
                                        isDescending = false;
                                        break;
                                    }
                                }
                            }
                            currentIndex--;
                        }
                        else
                        {
                            isDescending = false;
                        }

                    }
                    if (isLibrary)
                    {
                        targetPeptide.apexPeakLibrary = currentMax;
                        targetPeptide.leftPeakLibrary = currentMinLeft;
                        targetPeptide.rightPeakLibrary = currentMinRight;
                    }
                    else
                    {
                        targetPeptide.apexPeakExperiment = currentMax;
                        targetPeptide.leftPeakExperiment = currentMinLeft;
                        targetPeptide.rightPeakExperiment = currentMinRight;
                    }
                }
                else //assume the max lies to the left;
                {
                    bool isAscending = true;
                    bool isDescending = false;
                    RTPeak currentMax = XICPeaks[startIndex];
                    double threshold = currentMax.Intensity * 0.01;
                    RTPeak currentMinRight = null;
                    RTPeak currentMinLeft = null;
                    int currentIndex = startIndex - 1;
                    while (isAscending)
                    {
                        if (currentIndex >= 0)
                        {
                            var currentPeak = XICPeaks[currentIndex];
                            if (currentPeak.Intensity > currentMax.Intensity) //Reset max and threshold, continue on
                            {
                                currentMax = currentPeak;
                                threshold = currentMax.Intensity * 0.01;
                                currentIndex--;
                                continue;
                            }
                            else //Could be that I was at the apex, or I'm in a dip
                            {
                                if (currentPeak.Intensity < (currentMax.Intensity * (1.0 / 1.8)))
                                //I was at the max, so I want to keep this point
                                {
                                    //and start looking for the right minimum
                                    isAscending = false;
                                    isDescending = true;
                                    break;
                                }
                                else
                                //could be the start of the descent or some weird dip. Look ahead to the next peaks to find out.
                                {
                                    int AdditionalCounter = 1;
                                    while (true) //keep in this loop until one of the conditions is satisfied.
                                    {
                                        int aheadIndex = currentIndex - AdditionalCounter;
                                        if ((aheadIndex) >= 0)
                                        {

                                            var aheadPeak = XICPeaks[aheadIndex];
                                            if (aheadPeak.Intensity < (currentMax.Intensity * (1.0 / 1.8)))
                                            {
                                                isAscending = false;
                                                isDescending = true;
                                                break;
                                            }
                                            if (aheadPeak.Intensity > currentMax.Intensity)
                                            {
                                                currentMax = aheadPeak;
                                                currentIndex -= AdditionalCounter;
                                                //currentIndex--;
                                                break;
                                            }
                                            AdditionalCounter++;
                                        }
                                        else
                                        {
                                            isAscending = false;
                                            isDescending = true;
                                            break;
                                        }
                                    }
                                }
                            }
                        }
                    }
                    currentIndex--;
                    if (currentIndex >= 0)
                    {
                        currentMinLeft = XICPeaks[currentIndex];
                    }
                    currentIndex--;
                    isDescending = true;
                    while (isDescending)
                    {
                        if (currentIndex >= 0)
                        {
                            var currentPeak = XICPeaks[currentIndex];
                            if (currentPeak.Intensity <= currentMinLeft.Intensity)
                            {
                                currentMinLeft = currentPeak;
                                if (currentPeak.Intensity < threshold)
                                {
                                    isDescending = false;
                                    break;
                                }
                            }
                            if (currentPeak.Intensity > currentMinLeft.Intensity)
                            //could be the local min or some weird dip
                            {
                                int AdditionalCounter = 1;
                                while (true)
                                {
                                    int aheadIndex = currentIndex - AdditionalCounter;
                                    if ((aheadIndex) >= 0)
                                    {

                                        var aheadPeak = XICPeaks[aheadIndex];
                                        if (aheadPeak.Intensity > (1.8 * currentMinLeft.Intensity))
                                        //There is a new peak coming up. You have found the local min. Break out.
                                        {
                                            isDescending = false;
                                            break;
                                        }
                                        if (aheadPeak.Intensity < currentMinLeft.Intensity || aheadPeak.Intensity == 0)
                                        {
                                            currentMinLeft = aheadPeak;
                                            currentIndex -= AdditionalCounter;
                                            //currentIndex--;
                                            break;
                                        }
                                        AdditionalCounter++;
                                    }
                                    else
                                    {
                                        isDescending = false;
                                        break;
                                    }
                                }
                            }
                            currentIndex--;
                        }
                        else
                        {
                            isDescending = false;
                        }

                    }
                    isDescending = true;
                    currentIndex = startIndex + 1;
                    if (currentIndex < XICPeaks.Count)
                    {
                        currentMinRight = XICPeaks[currentIndex];
                    }
                    currentIndex++;
                    while (isDescending)
                    {
                        if (currentIndex < XICPeaks.Count)
                        {
                            var currentPeak = XICPeaks[currentIndex];
                            if (currentPeak.Intensity <= currentMinRight.Intensity)
                            {
                                currentMinRight = currentPeak;
                                if (currentPeak.Intensity < threshold)
                                {
                                    isDescending = false;
                                    break;
                                }
                            }
                            if (currentPeak.Intensity > currentMinRight.Intensity)
                            //could be the local min or some weird dip
                            {
                                int AdditionalCounter = 1;
                                while (true)
                                {
                                    int aheadIndex = currentIndex + AdditionalCounter;
                                    if ((aheadIndex) < XICPeaks.Count)
                                    {

                                        var aheadPeak = XICPeaks[aheadIndex];
                                        if (aheadPeak.Intensity > (1.8 * currentMinRight.Intensity))
                                        //There is a new peak coming up. You have found the local min. Break out.
                                        {
                                            isDescending = false;
                                            break;
                                        }
                                        if (aheadPeak.Intensity < currentMinRight.Intensity || aheadPeak.Intensity == 0)
                                        {
                                            currentMinRight = aheadPeak;
                                            currentIndex += AdditionalCounter;
                                            //currentIndex++;
                                            break;
                                        }
                                        AdditionalCounter++;
                                    }
                                    else
                                    {
                                        isDescending = false;
                                        break;
                                    }
                                }
                            }
                            currentIndex++;
                        }
                        else
                        {
                            isDescending = false;
                        }

                    }
                    if (isLibrary)
                    {
                        targetPeptide.apexPeakLibrary = currentMax;
                        targetPeptide.leftPeakLibrary = currentMinLeft;
                        targetPeptide.rightPeakLibrary = currentMinRight;
                    }
                    else
                    {
                        targetPeptide.apexPeakExperiment = currentMax;
                        targetPeptide.leftPeakExperiment = currentMinLeft;
                        targetPeptide.rightPeakExperiment = currentMinRight;
                    }
                }
            }
            catch (Exception e)
            {
                if (isLibrary)
                {
                    targetPeptide.apexPeakLibrary = new RTPeak(targetPeptide.AdjustedMZ_Library, 0, targetPeptide.parentMS1Time);
                    targetPeptide.leftPeakLibrary = new RTPeak(targetPeptide.AdjustedMZ_Library, 0, targetPeptide.parentMS1Time);
                    targetPeptide.rightPeakLibrary = new RTPeak(targetPeptide.AdjustedMZ_Library, 0, targetPeptide.parentMS1Time);
                }
                else
                {
                    targetPeptide.apexPeakExperiment = new RTPeak(targetPeptide.AdjustedMZ_Library, 0, targetPeptide.parentMS1Time);
                    targetPeptide.leftPeakExperiment = new RTPeak(targetPeptide.AdjustedMZ_Library, 0, targetPeptide.parentMS1Time);
                    targetPeptide.rightPeakExperiment = new RTPeak(targetPeptide.AdjustedMZ_Library, 0, targetPeptide.parentMS1Time);
                }
                return;
            }
        }

        public static void GetFWHMWindow(LFPeptide peptide)
        {
            var apexPeak = peptide.apexPeakLibrary;
            var threshold = apexPeak.Intensity * .5;
            var startIndex = peptide.SmoothLibrary.BinarySearch(apexPeak);
            if (startIndex < 0)
            {
                startIndex = ~startIndex;
            }

            RTPeak LeftFWHM = null;
            RTPeak RightFWHM = null;
            for (int i = startIndex; i < peptide.SmoothLibrary.Count; i++)
            {
                var currPeak = peptide.SmoothLibrary[i];
                if (currPeak.Intensity < threshold)
                {
                    break;
                }
                RightFWHM = currPeak;
            }
            for (int i = startIndex; i >= 0; i--)
            {
                var currPeak = peptide.SmoothLibrary[i];
                if (currPeak.Intensity < threshold)
                {
                    break;
                }
                LeftFWHM = currPeak;
            }

            peptide.LeftFWHM = LeftFWHM;
            peptide.RightFWHM = RightFWHM;
        }

        public static ChromaPeak GetChromaPeak(RTPeak startPeak, List<RTPeak> XICPeaks)
        {
            int startIndex = XICPeaks.BinarySearch(startPeak);
            if (startIndex < 0)
            {
                startIndex = ~startIndex;
            }
            int leftIndex = 0;
            int apexIndex = startIndex;
            int rightIndex = 0;

            double sumLeft = 0;
            double sumRight = 0;

            for (int i = 1; i <= 5; i++)
            {
                if ((apexIndex - i) >= 0)
                {
                    sumLeft += XICPeaks[apexIndex - i].Intensity;
                }
                if ((apexIndex + i) < XICPeaks.Count)
                {
                    sumRight += XICPeaks[apexIndex + i].Intensity;
                }
            }
            if (sumLeft == 0 && sumRight == 0)
            {
                return null;
            }
            var returnPeak = new ChromaPeak();
            returnPeak.initialPeak = startPeak;
            try
            {
                if (sumRight > sumLeft) //assume that the max lies to the right
                {
                    bool isAscending = true;
                    bool isDescending = false;
                    RTPeak currentMax = XICPeaks[startIndex];
                    double threshold = currentMax.Intensity * 0.01;
                    RTPeak currentMinRight = null;
                    RTPeak currentMinLeft = null;
                    int currentIndex = 1 + startIndex;
                    while (isAscending)
                    {
                        if (currentIndex < XICPeaks.Count)
                        {
                            var currentPeak = XICPeaks[currentIndex];
                            if (currentPeak.Intensity > currentMax.Intensity) //Reset max and threshold, continue on
                            {
                                currentMax = currentPeak;
                                threshold = currentMax.Intensity * 0.01;
                                currentIndex++;
                                continue;
                            }
                            else //Could be that I was at the apex, or I'm in a dip
                            {
                                if (currentPeak.Intensity < (currentMax.Intensity * (1.0 / 1.8)))
                                //I was at the max, so I want to keep this point
                                {
                                    //and start looking for the right minimum
                                    isAscending = false;
                                    isDescending = true;
                                    break;
                                }
                                else
                                //could be the start of the descent or some weird dip. Look ahead to the next peaks to find out.
                                {
                                    int AdditionalCounter = 1;
                                    while (true) //keep in this loop until one of the conditions is satisfied.
                                    {
                                        int aheadIndex = currentIndex + AdditionalCounter;
                                        if ((aheadIndex) < XICPeaks.Count)
                                        {

                                            var aheadPeak = XICPeaks[aheadIndex];
                                            if (aheadPeak.Intensity < (currentMax.Intensity * (1.0 / 1.8)))
                                            {
                                                isAscending = false;
                                                isDescending = true;
                                                break;
                                            }
                                            if (aheadPeak.Intensity > currentMax.Intensity)
                                            {
                                                currentMax = aheadPeak;
                                                currentIndex += AdditionalCounter;
                                                //currentIndex++;
                                                break;
                                            }
                                            AdditionalCounter++;
                                        }
                                        else
                                        {
                                            isAscending = false;
                                            isDescending = true;
                                            break;
                                        }
                                    }
                                }
                            }
                        }
                    }
                    //you have located the max. assume that the peak to right is the local minima
                    currentIndex++;
                    if (currentIndex < XICPeaks.Count)
                    {
                        currentMinRight = XICPeaks[currentIndex];
                    }
                    currentIndex++;
                    while (isDescending)
                    {
                        if (currentIndex < XICPeaks.Count)
                        {
                            var currentPeak = XICPeaks[currentIndex];
                            if (currentPeak.Intensity <= currentMinRight.Intensity)
                            {
                                currentMinRight = currentPeak;
                                if (currentPeak.Intensity < threshold)
                                {
                                    isDescending = false;
                                    break;
                                }
                            }
                            if (currentPeak.Intensity > currentMinRight.Intensity)
                            //could be the local min or some weird dip
                            {
                                int AdditionalCounter = 1;
                                while (true)
                                {
                                    int aheadIndex = currentIndex + AdditionalCounter;
                                    if ((aheadIndex) < XICPeaks.Count)
                                    {

                                        var aheadPeak = XICPeaks[aheadIndex];
                                        if (aheadPeak.Intensity > (1.8 * currentMinRight.Intensity))
                                        //There is a new peak coming up. You have found the local min. Break out.
                                        {
                                            isDescending = false;
                                            break;
                                        }
                                        if (aheadPeak.Intensity < currentMinRight.Intensity || aheadPeak.Intensity == 0)
                                        {
                                            currentMinRight = aheadPeak;
                                            currentIndex += AdditionalCounter;
                                            //currentIndex++;
                                            break;
                                        }
                                        AdditionalCounter++;
                                    }
                                    else
                                    {
                                        isDescending = false;
                                        break;
                                    }
                                }
                            }
                            currentIndex++;
                        }
                        else
                        {
                            isDescending = false;
                        }
                    }
                    //you now have the max, and the right minima. Assume that the peak to the left of the start point is the left minima.
                    currentIndex = startIndex - 1;
                    if (currentIndex >= 0)
                    {
                        currentMinLeft = XICPeaks[currentIndex];
                    }
                    currentIndex--;
                    isDescending = true;
                    while (isDescending)
                    {
                        if (currentIndex >= 0)
                        {
                            var currentPeak = XICPeaks[currentIndex];
                            if (currentPeak.Intensity <= currentMinLeft.Intensity)
                            {
                                currentMinLeft = currentPeak;
                                if (currentPeak.Intensity < threshold)
                                {
                                    isDescending = false;
                                    break;
                                }
                            }
                            if (currentPeak.Intensity > currentMinLeft.Intensity)
                            //could be the local min or some weird dip
                            {
                                int AdditionalCounter = 1;
                                while (true)
                                {
                                    int aheadIndex = currentIndex - AdditionalCounter;
                                    if ((aheadIndex) >= 0)
                                    {

                                        var aheadPeak = XICPeaks[aheadIndex];
                                        if (aheadPeak.Intensity > (1.8 * currentMinLeft.Intensity))
                                        //There is a new peak coming up. You have found the local min. Break out.
                                        {
                                            isDescending = false;
                                            break;
                                        }
                                        if (aheadPeak.Intensity < currentMinLeft.Intensity || aheadPeak.Intensity == 0)
                                        {
                                            currentMinLeft = aheadPeak;
                                            currentIndex -= AdditionalCounter;
                                            //currentIndex--;
                                            break;
                                        }
                                        AdditionalCounter++;
                                    }
                                    else
                                    {
                                        isDescending = false;
                                        break;
                                    }
                                }
                            }
                            currentIndex--;
                        }
                        else
                        {
                            isDescending = false;
                        }

                    }

                    returnPeak.apexPeak = currentMax;
                    returnPeak.leftPeak = currentMinLeft;
                    returnPeak.rightPeak = currentMinRight;
                    return returnPeak;
                }
                else //assume the max lies to the left;
                {
                    bool isAscending = true;
                    bool isDescending = false;
                    RTPeak currentMax = XICPeaks[startIndex];
                    double threshold = currentMax.Intensity * 0.01;
                    RTPeak currentMinRight = null;
                    RTPeak currentMinLeft = null;
                    int currentIndex = startIndex - 1;
                    while (isAscending)
                    {
                        if (currentIndex >= 0)
                        {
                            var currentPeak = XICPeaks[currentIndex];
                            if (currentPeak.Intensity > currentMax.Intensity) //Reset max and threshold, continue on
                            {
                                currentMax = currentPeak;
                                threshold = currentMax.Intensity * 0.01;
                                currentIndex--;
                                continue;
                            }
                            else //Could be that I was at the apex, or I'm in a dip
                            {
                                if (currentPeak.Intensity < (currentMax.Intensity * (1.0 / 1.8)))
                                //I was at the max, so I want to keep this point
                                {
                                    //and start looking for the right minimum
                                    isAscending = false;
                                    isDescending = true;
                                    break;
                                }
                                else
                                //could be the start of the descent or some weird dip. Look ahead to the next peaks to find out.
                                {
                                    int AdditionalCounter = 1;
                                    while (true) //keep in this loop until one of the conditions is satisfied.
                                    {
                                        int aheadIndex = currentIndex - AdditionalCounter;
                                        if ((aheadIndex) >= 0)
                                        {

                                            var aheadPeak = XICPeaks[aheadIndex];
                                            if (aheadPeak.Intensity < (currentMax.Intensity * (1.0 / 1.8)))
                                            {
                                                isAscending = false;
                                                isDescending = true;
                                                break;
                                            }
                                            if (aheadPeak.Intensity > currentMax.Intensity)
                                            {
                                                currentMax = aheadPeak;
                                                currentIndex -= AdditionalCounter;
                                                //currentIndex--;
                                                break;
                                            }
                                            AdditionalCounter++;
                                        }
                                        else
                                        {
                                            isAscending = false;
                                            isDescending = true;
                                            break;
                                        }
                                    }
                                }
                            }
                        }
                        else
                        {
                            isAscending = false;
                            break;
                        }
                    }
                    currentIndex--;
                    if (currentIndex >= 0)
                    {
                        currentMinLeft = XICPeaks[currentIndex];
                    }
                    currentIndex--;
                    isDescending = true;
                    while (isDescending)
                    {
                        if (currentIndex >= 0)
                        {
                            var currentPeak = XICPeaks[currentIndex];
                            if (currentPeak.Intensity <= currentMinLeft.Intensity)
                            {
                                currentMinLeft = currentPeak;
                                if (currentPeak.Intensity < threshold)
                                {
                                    isDescending = false;
                                    break;
                                }
                            }
                            if (currentPeak.Intensity > currentMinLeft.Intensity)
                            //could be the local min or some weird dip
                            {
                                int AdditionalCounter = 1;
                                while (true)
                                {
                                    int aheadIndex = currentIndex - AdditionalCounter;
                                    if ((aheadIndex) >= 0)
                                    {

                                        var aheadPeak = XICPeaks[aheadIndex];
                                        if (aheadPeak.Intensity > (1.8 * currentMinLeft.Intensity))
                                        //There is a new peak coming up. You have found the local min. Break out.
                                        {
                                            isDescending = false;
                                            break;
                                        }
                                        if (aheadPeak.Intensity < currentMinLeft.Intensity || aheadPeak.Intensity == 0)
                                        {
                                            currentMinLeft = aheadPeak;
                                            currentIndex -= AdditionalCounter;
                                            //currentIndex--;
                                            break;
                                        }
                                        AdditionalCounter++;
                                    }
                                    else
                                    {
                                        isDescending = false;
                                        break;
                                    }
                                }
                            }
                            currentIndex--;
                        }
                        else
                        {
                            isDescending = false;
                        }

                    }
                    isDescending = true;
                    currentIndex = startIndex + 1;
                    if (currentIndex < XICPeaks.Count)
                    {
                        currentMinRight = XICPeaks[currentIndex];
                    }
                    currentIndex++;
                    while (isDescending)
                    {
                        if (currentIndex < XICPeaks.Count)
                        {
                            var currentPeak = XICPeaks[currentIndex];
                            if (currentPeak.Intensity <= currentMinRight.Intensity)
                            {
                                currentMinRight = currentPeak;
                                if (currentPeak.Intensity < threshold)
                                {
                                    isDescending = false;
                                    break;
                                }
                            }
                            if (currentPeak.Intensity > currentMinRight.Intensity)
                            //could be the local min or some weird dip
                            {
                                int AdditionalCounter = 1;
                                while (true)
                                {
                                    int aheadIndex = currentIndex + AdditionalCounter;
                                    if ((aheadIndex) < XICPeaks.Count)
                                    {

                                        var aheadPeak = XICPeaks[aheadIndex];
                                        if (aheadPeak.Intensity > (1.8 * currentMinRight.Intensity))
                                        //There is a new peak coming up. You have found the local min. Break out.
                                        {
                                            isDescending = false;
                                            break;
                                        }
                                        if (aheadPeak.Intensity < currentMinRight.Intensity || aheadPeak.Intensity == 0)
                                        {
                                            currentMinRight = aheadPeak;
                                            currentIndex += AdditionalCounter;
                                            //currentIndex++;
                                            break;
                                        }
                                        AdditionalCounter++;
                                    }
                                    else
                                    {
                                        isDescending = false;
                                        break;
                                    }
                                }
                            }
                            currentIndex++;
                        }
                        else
                        {
                            isDescending = false;
                        }
                    }

                    returnPeak.apexPeak = currentMax;
                    returnPeak.leftPeak = currentMinLeft;
                    returnPeak.rightPeak = currentMinRight;
                    return returnPeak;
                }
            }
            catch (Exception e)
            {
                return null;
            }
        }
    }
}
