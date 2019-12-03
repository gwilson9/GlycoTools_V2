using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using CSMSL;
using CSMSL.Spectral;
using CSMSL.IO.Thermo;
using CSMSL.Proteomics;

namespace _20190618_GlycoTools_V2
{
    public class LFPeptide : ICloneable
    {
        public string sequence;
        public int charge;
        public int ms2ScanNumber;
        public double TheoreticalMZ;
        public double AdjustedMZ_Library;
        public double AdjustedMZ_Experiment;
        public List<RTPeak> XICLibrary;
        public List<RTPeak> XICUsedPeaks;
        public List<RTPeak> XICExperiment;
        public List<RTPeak> SmoothExperiment;
        public List<RTPeak> SmoothLibrary;
        public int parentMS1;
        public double parentMS1Time;
        public int startIndex;
        public double UserMZ;
        public double Area;
        public double startLookupTime;
        public double stopLookupTime;
        public double IdentificationRT;
        public double ExperimentRT;
        public bool include;
        public int FirstScan;
        public int LastScan;
        public RTPeak apexPeakLibrary;
        public RTPeak leftPeakLibrary;
        public RTPeak rightPeakLibrary;
        public RTPeak apexPeakExperiment;
        public RTPeak leftPeakExperiment;
        public RTPeak rightPeakExperiment;
        public int RunScanNumber;
        public MZPeak RunPeak;
        public double PPMError_LibraryToTheoretical;
        public double PPMError_ExperimentToLibrary;
        public RTPeak ExperimentLookupPeak;
        public List<RTPeak> DotProductLibPeaks;
        public List<RTPeak> OffsetExpPeaks;
        public List<RTPeak> OffsetExpPeaksSmooth;
        public DoubleRange lookupRange;
        public List<ChromaPeak> ExpChromaPeaks;
        public bool doneBuildingXIC = false;
        public bool extractedXIC = false;
        public double Mass;
        public double ppmError;

        public RTPeak LeftFWHM;
        public RTPeak RightFWHM;

        public LFPeptide()
        {
            XICLibrary = new List<RTPeak>();
            XICUsedPeaks = new List<RTPeak>();
            SmoothLibrary = new List<RTPeak>();
            XICExperiment = new List<RTPeak>();
            SmoothExperiment = new List<RTPeak>();
            OffsetExpPeaks = new List<RTPeak>();
            OffsetExpPeaksSmooth = new List<RTPeak>();
            ExpChromaPeaks = new List<ChromaPeak>();
        }

        public Object Clone()
        {
            return this.MemberwiseClone();
        }

        public void AddLibraryPeak(MZPeak peak, double RT)
        {
            var newRTPeak = new RTPeak(peak, RT);
            XICLibrary.Add(newRTPeak);
        }

        public void AddLibraryPeak(RTPeak peak)
        {
            XICLibrary.Add(peak);
        }

        public void AddExperimentPeak(MZPeak peak, double RT)
        {
            var newRTPeak = new RTPeak(peak, RT);
            XICExperiment.Add(newRTPeak);
        }

        public void AddExperimentPeak(RTPeak peak)
        {
            XICExperiment.Add(peak);
        }

        public override string ToString()
        {
            return sequence;
        }
    }
}
