using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using CSMSL;
using CSMSL.Spectral;
using CSMSL.Proteomics;
using CSMSL.IO.Thermo;
using LumenWorks.Framework.IO.Csv;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearRegression;
using MathNet.Numerics;

namespace _20190618_GlycoTools_V2
{
    class LFQProcessor
    {

        // Assign targetPeptides to the list of PSMs in LFQProcessor 

        public string rawpath;
        public string peptidesPath;
        public string appendedPath;
        public static List<LFPeptide> targetPeptides;

        public LFQProcessor(string rawpath)
        {
            this.rawpath = rawpath;
            //this.peptidesPath = peptidesPath;
            //this.appendedPath = Path.GetDirectoryName(peptidesPath) + "\\" + Path.GetFileName(peptidesPath) + "_quant.csv";
            targetPeptides = new List<LFPeptide>();
        }


        public List<PSM> crunch(List<PSM> psms)
        {
            var rawFile = new ThermoRawFile(rawpath);
            rawFile.Open();

            List<int> msOneScanNumbers = new List<int>();
            foreach (var scan in rawFile)
            {
                if (scan.MsnOrder == 1)
                {
                    msOneScanNumbers.Add(scan.SpectrumNumber);
                }
            }
            //ReadInTargets(peptidesPath, rawFile);

            foreach (PSM psm in psms)
            {
                LFPeptide pep = psmToLFPeptide(psm, rawFile);
                targetPeptides.Add(pep);
            }

            // Question: how do we sync the list of PSMs with the intensity of the LFPeptide
            // They will be in same order

            List<LFPeptide> pepsWithSmooth = new List<LFPeptide>();

            if (rawFile.GetMzAnalyzer(msOneScanNumbers[0]).Equals(MZAnalyzerType.Orbitrap))
            {
                var firstSpectrumRT = rawFile.GetRetentionTime(rawFile.FirstSpectrumNumber);
                var lastSpectrumRT = rawFile.GetRetentionTime(rawFile.LastSpectrumNumber);

                foreach (var pep in targetPeptides)
                {
                    var startTime = Math.Max(pep.parentMS1Time - 2, firstSpectrumRT);
                    var stopTime = Math.Min(pep.parentMS1Time + 2, lastSpectrumRT);
                    pep.startLookupTime = startTime;
                    pep.stopLookupTime = stopTime;
                    pep.FirstScan = rawFile.GetSpectrumNumber(pep.startLookupTime);
                    pep.LastScan = rawFile.GetSpectrumNumber(pep.stopLookupTime);
                    pep.lookupRange = DoubleRange.FromPPM(pep.UserMZ, 50);
                }

                foreach (var ms1 in msOneScanNumbers)
                {
                    var spectrum = rawFile.GetSpectrum(ms1);
                    GetXICs(spectrum, ms1, rawFile.GetRetentionTime(ms1));
                    rawFile.ClearCachedScans();

                    List<LFPeptide> peptidesForExtract = targetPeptides.Where(x => x.doneBuildingXIC && !x.extractedXIC).ToList();
                    foreach (var pep in peptidesForExtract)
                    {
                        try
                        {
                            pep.SmoothLibrary = Smoothing.GetRollingAveragePeaks(pep, 11, true);
                            if (pep.SmoothLibrary.Count != 0)
                            {
                                ExtractFeatures.GetApexPeak(pep, true);
                                ExtractFeatures.GetFWHMWindow(pep);
                                LFPeptide pepSmooth = (LFPeptide)pep.Clone();
                                pepsWithSmooth.Add(pepSmooth);
                                pep.XICLibrary.Clear();
                            }

                        }
                        catch (Exception e)
                        {/**
                            System.Windows.Forms.MessageBox.Show("XICLibraryCount: " + pep.XICLibrary.Count 
                                                                    + "Peptide: " + pep.sequence
                                                                    + "\n" + e.ToString());
    **/
                        }
                    }
                }

                List<LFPeptide> peptidesToWrite = targetPeptides.Where(x => x.doneBuildingXIC && !x.extractedXIC).ToList();

                //Changed from psms.Count to pepswithsmooth.Count
                //for(int i = 0; i < psms.Count; i++)
                for (int i = 0; i < pepsWithSmooth.Count; i++)
                {
                    /**
                    if(psms[i].scanTime == 20.45)
                    {
                        double[] xs = new double[pepsWithSmooth[i].SmoothLibrary.Count];
                        double[] ys = new double[pepsWithSmooth[i].SmoothLibrary.Count];

                        StreamWriter writer = new StreamWriter(@"C:\Users\gwilson\Desktop\GlycoQuant - Injs\TLN[NGlycan_1216.42286271]CSGAHVK_0.5_2_PolynomialFit.csv");
                        for (int k = 0; k < pepsWithSmooth[i].SmoothLibrary.Count; k++)
                        {
                            xs[k] = pepsWithSmooth[i].SmoothLibrary[k].RT;
                            ys[k] = pepsWithSmooth[i].SmoothLibrary[k].Intensity;                            

                            //writer.WriteLine(pepsWithSmooth[i].SmoothLibrary[k].RT + "," + pepsWithSmooth[i].SmoothLibrary[k].Intensity);
                            
                            
                            if (pepsWithSmooth[i].SmoothLibrary[k].RT > targetPeptides[i].LeftFWHM.RT && pepsWithSmooth[i].SmoothLibrary[k].RT < targetPeptides[i].RightFWHM.RT)
                            {
                                writer.WriteLine(pepsWithSmooth[i].SmoothLibrary[k].RT + "," + pepsWithSmooth[i].SmoothLibrary[k].Intensity);
                            }
                            
                        }
                                                
                        double[] p = Fit.Polynomial(xs, ys, 3)

                        for(int j = 0; j < p.Length; j++)
                        {
                            //writer.WriteLine()
                        }

                        writer.Close();
                    }
                    **/


                    double intensity = 0;

                    try
                    {
                        for (int j = 0; j < pepsWithSmooth[i].SmoothLibrary.Count - 1; j++)
                        {
                            // Intensity is the sum of peak intensities from Left FWHM to Right FWHM
                            /**
                            if (pepsWithSmooth[i].SmoothLibrary[j].RT > targetPeptides[i].LeftFWHM.RT && pepsWithSmooth[i].SmoothLibrary[j].RT < targetPeptides[i].RightFWHM.RT)
                            {
                                intensity += pepsWithSmooth[i].SmoothLibrary[j].Intensity;
                            }
                            **/

                            // Retention times
                            double RT1 = pepsWithSmooth[i].SmoothLibrary[j].RT;
                            double RT2 = pepsWithSmooth[i].SmoothLibrary[j + 1].RT;

                            // Intensities
                            double int1 = pepsWithSmooth[i].SmoothLibrary[j].Intensity;
                            double int2 = pepsWithSmooth[i].SmoothLibrary[j + 1].Intensity;

                            //Rectangle area
                            double rectArea = (RT2 - RT1) * Math.Min(int1, int2);

                            //Triangle Area
                            double triArea = ((RT2 - RT1) * Math.Abs(int1 - int2)) / 2;

                            intensity += rectArea + triArea;

                        }
                    }
                    catch (Exception e)
                    {
                        Console.ReadKey();
                    }

                    psms[i].intensity = intensity;
                    psms[i].peakElution = pepsWithSmooth[i].SmoothLibrary;
                    psms[i].scanNumberofMaxElutionIntensity = rawFile.GetSpectrumNumber(pepsWithSmooth[i].SmoothLibrary.Aggregate((i1, i2) => i1.Intensity > i2.Intensity ? i1 : i2).RT);
                    //psms[i].intensity = targetPeptides[i].apexPeakLibrary.Intensity;
                }

                rawFile.Dispose();

                return psms;
            }

            rawFile.Dispose();
            return psms;
        }

        public static void GetXICs(ThermoSpectrum currentSpectrum, int specNumber, double rt)
        {
            List<LFPeptide> donePeptides = targetPeptides.Where(x => x.LastScan < specNumber).ToList();
            foreach (var pep in donePeptides)
            {
                pep.doneBuildingXIC = true;
            }

            List<LFPeptide> currPeptides = targetPeptides.Where(x => x.FirstScan <= specNumber && x.LastScan >= specNumber).ToList();
            foreach (var pep in currPeptides)
            {
                List<ThermoMzPeak> outPeaks = new List<ThermoMzPeak>();
                if (currentSpectrum.TryGetPeaks(pep.lookupRange, out outPeaks))
                {
                    var matchPeak = GetClosestPeak(outPeaks, pep.UserMZ);
                    var newRTPeak = new RTPeak(matchPeak, rt);
                    pep.XICLibrary.Add(newRTPeak);
                }
                else
                {
                    var newRTPeak = new RTPeak(pep.UserMZ, 0, rt);
                    pep.XICLibrary.Add(newRTPeak);
                }
            }
        }

        public static ThermoMzPeak GetClosestPeak(List<ThermoMzPeak> peaks, double expMZ)
        {
            double diff = double.MaxValue;
            ThermoMzPeak returnPeak = null;
            foreach (var peak in peaks)
            {
                var currDiff = Math.Abs(peak.MZ - expMZ);
                if (currDiff < diff)
                {
                    diff = currDiff;
                    returnPeak = peak;
                }
            }
            return returnPeak;
        }

        public LFPeptide psmToLFPeptide(PSM psm, ThermoRawFile rawFile)
        {
            var pep = new LFPeptide
            {
                ms2ScanNumber = psm.scanNumber,
                sequence = psm.peptide.Sequence,
                charge = psm.charge,
                TheoreticalMZ = psm.mzCalc
            };

            var parentScan = rawFile.GetParentSpectrumNumber(pep.ms2ScanNumber);
            while (true)
            {                
                if (rawFile.GetMsnOrder(parentScan) == 2)
                {
                    parentScan = rawFile.GetParentSpectrumNumber(parentScan);
                }
                else
                {
                    break;
                }
            }
            pep.parentMS1 = psm.MasterScan;
            pep.parentMS1Time = rawFile.GetRetentionTime(psm.MasterScan);
            //pep.UserMZ = psm.mzObs;
            pep.UserMZ = rawFile.GetPrecursorMz(pep.ms2ScanNumber);

            return pep;
        }


        public static void ReadInTargets(string filePath, ThermoRawFile rawFile)
        {
            targetPeptides = new List<LFPeptide>();
            using (var csv = new CsvReader(new StreamReader(filePath), true, '\t'))
            {
                int fieldCount = csv.FieldCount;
                string[] headers = csv.GetFieldHeaders();
                //"Spectrum number"
                //"Peptide"
                //"Charge"
                //"Precursor Theoretical m/z (Th)"
                while (csv.ReadNextRecord())
                {
                    var pep = new LFPeptide
                    {
                        ms2ScanNumber = int.Parse(csv["Spectrum Number"]),
                        sequence = csv["Peptide"],
                        charge = int.Parse(csv["Charge"]),
                        TheoreticalMZ = double.Parse(csv["Precursor Theoretical m/z (Th)"])
                    };
                    pep.parentMS1 = rawFile.GetParentSpectrumNumber(pep.ms2ScanNumber);
                    pep.parentMS1Time = rawFile.GetRetentionTime(pep.parentMS1);
                    pep.UserMZ = rawFile.GetPrecursorMz(pep.ms2ScanNumber);
                    targetPeptides.Add(pep);
                }
            }
        }
    }
}
