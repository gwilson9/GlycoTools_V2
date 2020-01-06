using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
//using CenterSpace.NMath.Core;
//using CenterSpace.NMath.Stats;

namespace _20190618_GlycoTools_V2
{
    class Smoothing
    {
        /**
        public static List<RTPeak> SavitskyGolaySmooth(LFPeptide targetPeptide, int left, int right, int polyN)
        {
            var copyPeaks = new List<RTPeak>();
            foreach(var peak in targetPeptide.XICLibrary)
            {
                copyPeaks.Add(peak);
            }

            copyPeaks = FillInXICGaps(copyPeaks);

            var c = MovingWindowFilter.SavitzkyGolayCoefficients(left, right, polyN);

            var filter = new MovingWindowFilter(left, right, c);

            var x = new DoubleVector(copyPeaks.Select(a => a.RT).ToList().ToArray());

            var y = new DoubleVector(copyPeaks.Select(a => a.Intensity).ToList().ToArray());

            var yFiltered = filter.Filter(y, MovingWindowFilter.BoundaryOption.PadWithZeros);

            var returnList = new List<RTPeak>();

            for (int i = 0; i < x.Count(); i++)
            {
                var peak = new RTPeak(targetPeptide.UserMZ, yFiltered[i], x[i]);

                returnList.Add(peak);
            }

            return returnList;
        }
    **/

        public static List<RTPeak> NoSmoothing(LFPeptide targetPeptide)
        {
            var returnList = new List<RTPeak>();

            foreach(var peak in targetPeptide.XICLibrary)
            {
                var newPeak = new RTPeak(peak.MZ, peak.Intensity, peak.RT);
                returnList.Add(newPeak);
            }

            return returnList;
        }

        public static List<RTPeak> GetRollingAveragePeaks(LFPeptide targetPeptide, int period, bool isLibrary = true)
        {    
           
            RTPeak start = new RTPeak(0, 0, 0);
            RTPeak end = new RTPeak(0, 0, 100000);
            List<RTPeak> copyPeaks = new List<RTPeak>();
            if (isLibrary)
            {
                copyPeaks.Add(start);
                copyPeaks.Add(start);
                copyPeaks.Add(start);
                copyPeaks.Add(start);
                foreach (var peak in targetPeptide.XICLibrary)
                {
                    copyPeaks.Add(peak);
                }
                copyPeaks.Add(end);
                copyPeaks.Add(end);
                copyPeaks.Add(end);
                copyPeaks.Add(end);
                copyPeaks = FillInXICGaps(copyPeaks);
                return GetRollingAveragePeaks(copyPeaks, period);
            }
            else
            {
                copyPeaks.Add(start);
                copyPeaks.Add(start);
                copyPeaks.Add(start);
                copyPeaks.Add(start);
                foreach (var peak in targetPeptide.XICExperiment)
                {
                    copyPeaks.Add(peak);
                }
                copyPeaks.Add(end);
                copyPeaks.Add(end);
                copyPeaks.Add(end);
                copyPeaks.Add(end);
                copyPeaks = FillInXICGaps(copyPeaks);
                return GetRollingAveragePeaks(copyPeaks, period);
            }
        }

        private static List<RTPeak> FillInXICGaps(List<RTPeak> peaks)
        {
            List<RTPeak> returnList = new List<RTPeak>();
            returnList.Add(peaks[0]);
            for (int i = 1; i < peaks.Count - 1; i++)
            {
                if (peaks[i - 1].Intensity != 0 && peaks[i].Intensity == 0 && peaks[i + 1].Intensity != 0)
                {
                    //this is a drop out which confuses analysis so fill it in
                    double intermediateIntensity = peaks[i - 1].Intensity + peaks[i + 1].Intensity;
                    intermediateIntensity = intermediateIntensity / 2;
                    double intermediateSN = peaks[i - 1].Sn + peaks[i + 1].Sn;
                    intermediateSN = intermediateSN / 2;
                    RTPeak newPeak = new RTPeak(peaks[i].MZ, intermediateIntensity, intermediateSN, peaks[i].RT);
                    returnList.Add(newPeak);
                }
                else
                {
                    returnList.Add(peaks[i]);
                }
            }
            return returnList;
        }

        public static List<RTPeak> GetRollingAveragePeaks(List<RTPeak> peaks, int period)
        {
            List<RTPeak> outPeaks = new List<RTPeak>();
            if (peaks.Count > 0)
            {
                for (int i = 0; i < peaks.Count; i++)
                {
                    if (i >= period)
                    {
                        double total = 0;
                        double totalIntensity = 0;
                        for (int x = i; x > (i - period); x--)
                        {
                            total += peaks[x].SN;
                            totalIntensity += peaks[x].Intensity;
                        }
                        double average = total / (double)period;
                        double averageInt = totalIntensity / (double)period;
                        RTPeak newPeak = new RTPeak(peaks[i].MZ, averageInt, average, peaks[i - (period / 2)].RT);
                        outPeaks.Add(newPeak);
                    }
                }
            }
            return outPeaks;
        }

        public static List<RTPeak> GetBoxcarAvgPeaks(List<RTPeak> peaks, int period)
        {
            var tmpList = new List<RTPeak>();
            int peaksToAdd = (period / 2) + 1;
            var startPeak = new RTPeak(0, 0, 0);
            var stopPeak = new RTPeak(0, 0, 100000);
            for (int i = 0; i < peaksToAdd; i++)
            {
                tmpList.Add(startPeak);
            }
            foreach (var peak in peaks)
            {
                tmpList.Add(peak);
            }
            for (int i = 0; i < peaksToAdd - 1; i++)
            {
                tmpList.Add(stopPeak);
            }
            List<RTPeak> outPeaks = new List<RTPeak>();
            if (tmpList.Count > 0)
            {
                for (int i = 0; i < tmpList.Count; i++)
                {
                    if (i >= period)
                    {
                        double total = 0;
                        double totalIntensity = 0;
                        for (int x = i; x > (i - period); x--)
                        {
                            //total += tmpList[x].SN;
                            totalIntensity += tmpList[x].Intensity;
                        }
                        // double average = total / (double)period;
                        double averageInt = totalIntensity / (double)period;
                        RTPeak newPeak = new RTPeak(tmpList[i].MZ, averageInt, 0, tmpList[i - (period / 2)].RT);
                        outPeaks.Add(newPeak);
                    }
                }
            }
            return outPeaks;
        }
    }
}
