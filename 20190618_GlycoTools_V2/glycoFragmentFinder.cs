using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using CSMSL;
using CSMSL.IO.Thermo;
using CSMSL.Spectral;
using CSMSL.Util;
using LumenWorks.Framework.IO.Csv;
using CSMSL.Analysis;
using CSMSL.Proteomics;
using CSMSL.Chemistry;
using System.Data.SQLite;
using System.Diagnostics;

namespace _20190618_GlycoTools_V2
{
    class glycoFragmentFinder
    {
        public string rawfilePath { get; set; }
        public string outputPath { get; set; }
        public bool localizeGlycan { get; set; }
        public List<GlycoPSM> glycoPSMs { get; set; }
        //public StreamWriter outputEachFragment { get; set; }
        public SQLiteConnection sqlReader { get; set; }
        public Stopwatch watch { get; set; }
        public string useType;

        public glycoFragmentFinder(string rawfilePath, GlycoPSM glycoPSM)
        {
            this.rawfilePath = rawfilePath;
            this.outputPath = "";
            this.glycoPSMs = new List<GlycoPSM> { glycoPSM };
        }

        public glycoFragmentFinder(string rawfilePath, List<GlycoPSM> glycoPSMs, string outputPath, bool localize)
        {
            this.rawfilePath = rawfilePath;
            this.outputPath = outputPath;
            this.glycoPSMs = glycoPSMs;
            this.watch = new Stopwatch();
            watch.Start();
            //outputEachFragment = new StreamWriter(@outputPath + "\\" + Path.GetFileNameWithoutExtension(this.rawfilePath) + "_AnnotatedSpectra.txt");
            var sqlPath = string.Format("Data Source={0}\\MyDatabase.sqlite; Version=3;", outputPath);
            this.sqlReader = new SQLiteConnection(sqlPath);
            sqlReader.Open();

            var fragTableCreationQuery = "CREATE TABLE IF NOT EXISTS FragmentIons (File STRING, Condition STRING, Replicate STRING, ScanNumber INT, Peptide STRING, Sequence STRING, Glycan STRING, ByonicScore REAL, DeltaModScore REAL, LogProb REAL, IonType STRING, Fragment STRING, FragmentNumber INT, MZ REAL, Charge INT, Intensity REAL)";
            var fragTableCreationCommand = new SQLiteCommand(fragTableCreationQuery, sqlReader);
            fragTableCreationCommand.ExecuteNonQuery();

            var pepTableCreationQuery = "CREATE TABLE IF NOT EXISTS AnnotatedPeptides (File STRING, Condition STRING, Replicate STRING, Peptide STRING, Sequence STRING, Glycans STRING, ScanNumber INT, Fragmentation STRING, MZ REAL, Charge INT, PeptideMass REAL, GlycanMasses REAL, Score REAL, DeltaModScore REAL, FDR2D REAL, PEP2D REAL, SeqCoverage REAL, GlycanSeqCoverage REAL, PercentTICExplained REAL, NumBackboneFragments INT, NumYIons INT, NumOxoniumIons INT)";
            var pepTableCreationCommand = new SQLiteCommand(pepTableCreationQuery, sqlReader);
            pepTableCreationCommand.ExecuteNonQuery();

        }

        public void crunch()
        {
            try
            {
                var rawFile = new ThermoRawFile(rawfilePath);
                rawFile.Open();

                foreach (GlycoPSM glycoPSM in glycoPSMs)
                {
                    if (glycoPSM.glycans.Count > 0)
                    {                       

                        int glycoPosition = glycoPSM.glycanPositions[0];
                        double glycanMass = glycoPSM.glycanMasses[0];
                        string glycan = glycoPSM.glycans[0];
                        var spectrum = rawFile.GetSpectrum(glycoPSM.scanNumber);

                        List<FragmentMatch> peptideFragments = new List<FragmentMatch>();
                        List<FragmentMatch> peptideNeutralLossFragments = new List<FragmentMatch>();
                        List<FragmentMatch> peptideFragmentsMustIncludeGlycan = new List<FragmentMatch>();
                        List<FragmentMatch> peptideIntactGlycans = new List<FragmentMatch>();
                        List<FragmentMatch> oxoniumIons = new List<FragmentMatch>();

                        int hexCount = 0;
                        int hexNAcCount = 0;
                        int fucCount = 0;
                        int neuAcCount = 0;
                        int neuGcCount = 0;

                        foreach (string glycanInList in glycoPSM.glycans)
                        {
                            hexCount += GetSugarCount(glycanInList, "Hex");
                            hexNAcCount += GetSugarCount(glycanInList, "HexNAc");
                            fucCount += GetSugarCount(glycanInList, "Fuc");
                            neuAcCount += GetSugarCount(glycanInList, "NeuAc");
                            neuGcCount += GetSugarCount(glycanInList, "NeuGc");
                        }

                        //generate theoretical fragments
                        List<Fragment> theoFragments = new List<Fragment>();

                        List<Fragment> bTypeFragments = glycoPSM.peptide.Fragment(FragmentTypes.b).ToList();
                        List<Fragment> yTypeFragments = glycoPSM.peptide.Fragment(FragmentTypes.y).ToList();
                        List<Fragment> cTypeFragments = glycoPSM.peptide.Fragment(FragmentTypes.c).ToList();
                        List<Fragment> zdotTypeFragments = glycoPSM.peptide.Fragment(FragmentTypes.zdot).ToList();

                        theoFragments.AddRange(bTypeFragments);
                        theoFragments.AddRange(yTypeFragments);

                        var fragmentation = rawFile.GetDissociationType(glycoPSM.scanNumber).ToString();
                        if (fragmentation.Equals("ETD"))
                        {
                            theoFragments.AddRange(cTypeFragments);
                            theoFragments.AddRange(zdotTypeFragments);
                        }

                        bool notInterferingWithPrecursorPeaks = false;

                        //search for fragments in spectrum
                        foreach (Fragment theoFragment in theoFragments)
                        {
                            for (int i = glycoPSM.charge - 1; i > 0; i--)
                            {
                                double theoFragmentMZ = (theoFragment.MonoisotopicMass + (Constants.Hydrogen * i)) / ((double)i);

                                notInterferingWithPrecursorPeaks = NoPrecursorPeakInterference(theoFragmentMZ, glycoPSM.peptide, glycoPSM.charge, fragmentation);

                                if (notInterferingWithPrecursorPeaks)
                                {
                                    double theoFragmentIsoMZ = (theoFragmentMZ + (1 * (Constants.Hydrogen)) / ((double)i));
                                    var range = DoubleRange.FromPPM(theoFragmentMZ, 20);

                                    List<ThermoMzPeak> outPeaks = new List<ThermoMzPeak>();

                                    if (spectrum.TryGetPeaks(range, out outPeaks))
                                    {
                                        var closestPeak = GetClosestPeak(outPeaks, theoFragmentMZ);
                                        double intensity = closestPeak.Intensity;

                                        var rangeForIsotope = DoubleRange.FromPPM(theoFragmentIsoMZ, 20);
                                        List<ThermoMzPeak> outPeaksIso = new List<ThermoMzPeak>();
                                        if (spectrum.TryGetPeaks(rangeForIsotope, out outPeaksIso))
                                        {
                                            var closestPeakIso = GetClosestPeak(outPeaksIso, theoFragmentIsoMZ);
                                            intensity += closestPeakIso.Intensity;
                                        }

                                        double basePeak = spectrum.GetBasePeakIntensity();

                                        if (closestPeak.Charge == i)
                                        {
                                            FragmentMatch fragmentMatch = new FragmentMatch(theoFragment.ToString(), theoFragment.Type.ToString(), theoFragment.Number, i, closestPeak.MZ, intensity);
                                            peptideFragments.Add(fragmentMatch);
                                            peptideFragmentsMustIncludeGlycan.Add(fragmentMatch);
                                        }
                                        else if (rawFile.GetSpectrum(glycoPSM.scanNumber).TryGetPeaks(rangeForIsotope, out outPeaksIso))
                                        {
                                            FragmentMatch fragmentMatch = new FragmentMatch(theoFragment.ToString(), theoFragment.Type.ToString(), theoFragment.Number, i, closestPeak.MZ, intensity);
                                            peptideFragments.Add(fragmentMatch);
                                            peptideFragmentsMustIncludeGlycan.Add(fragmentMatch);
                                        }
                                        else if (intensity >= (0.01 * basePeak) && closestPeak.Charge == 0)
                                        {
                                            FragmentMatch fragmentMatch = new FragmentMatch(theoFragment.ToString(), theoFragment.Type.ToString(), theoFragment.Number, i, closestPeak.MZ, intensity);
                                            peptideFragments.Add(fragmentMatch);
                                            peptideFragmentsMustIncludeGlycan.Add(fragmentMatch);
                                        }
                                    }
                                }
                            }
                        }

                        //search for fragments without glycan that should have glycan
                        foreach (Fragment theoFragment in theoFragments)
                        {
                            LookForFragmentWithNoGlycan(rawFile, theoFragment, glycoPSM.charge, glycoPosition, glycanMass, peptideFragments, peptideNeutralLossFragments,
                                                        glycoPSM.scanNumber, glycoPSM.peptide, fragmentation);
                        }

                        //search for glycan peaks still attached to peptide
                        GetIntactPepGlycanFragments(glycoPSM, glycoPSM.charge, hexCount, hexNAcCount, fucCount, neuAcCount,
                                                    neuGcCount, peptideIntactGlycans, spectrum, glycan, glycanMass);
                        GetOxoniumIons(glycoPSM, hexCount, hexNAcCount, fucCount, neuAcCount, neuGcCount, glycan, glycanMass, spectrum, oxoniumIons);



                        double sequenceCoverage = GetSequenceCoverage(glycoPSM, peptideFragments);
                        double sequenceCoverageOnlyGlycoFragments = GetSequenceCoverage(glycoPSM, peptideFragmentsMustIncludeGlycan);

                        bool anyFragmentsContainGlycan = SeeIfFragmentsHaveGlycan(glycoPSM, peptideFragmentsMustIncludeGlycan, "");
                        bool ntermFragmentsContainGlycan = SeeIfFragmentsHaveGlycan(glycoPSM, peptideFragmentsMustIncludeGlycan, "Nterm");
                        bool ctermFragmentsContainGlycan = SeeIfFragmentsHaveGlycan(glycoPSM, peptideFragmentsMustIncludeGlycan, "Cterm");

                        double glycanSequenceCoverage = GetGlycanSequenceCoverage(glycoPSM, peptideIntactGlycans, hexCount, hexNAcCount, fucCount, neuAcCount, neuGcCount);

                        string fragmentList = "";
                        double fragmentIntensity = 0;

                        double TIC = rawFile.GetSpectrum(glycoPSM.scanNumber).TotalIonCurrent;

                        double numberOfNeutralLossFragments = peptideNeutralLossFragments.Count;

                        peptideFragments = peptideFragments.Distinct().ToList();
                        peptideNeutralLossFragments = peptideNeutralLossFragments.Distinct().ToList();
                        peptideIntactGlycans = peptideIntactGlycans.Distinct().ToList();
                        peptideFragmentsMustIncludeGlycan = peptideFragmentsMustIncludeGlycan.Distinct().ToList();
                        oxoniumIons = oxoniumIons.Distinct().ToList();

                        if (useType.Equals("BatchAnnotation"))
                        {
                            var pepFragmentIntensity = 0.0;
                            foreach(FragmentMatch fragment in peptideFragments)
                            {
                                var commandString = string.Format("INSERT INTO FragmentIons('File','Condition','Replicate','ScanNumber','Peptide','Sequence','Glycan','ByonicScore','DeltaModScore','LogProb','IonType','Fragment','FragmentNumber','MZ','Charge','Intensity') VALUES ({0})", FragToSQLString(glycoPSM, fragment, "Backbone"));
                                var command = new SQLiteCommand(commandString, sqlReader);
                                command.ExecuteNonQuery();

                                pepFragmentIntensity += fragment.fragmentSignal;
                            }

                            double glycanFragmentIntensity = 0.0;
                            foreach (FragmentMatch fragment in peptideIntactGlycans)
                            {
                                var commandString = string.Format("INSERT INTO FragmentIons('File','Condition','Replicate','ScanNumber','Peptide','Sequence','Glycan','ByonicScore','DeltaModScore','LogProb','IonType','Fragment','FragmentNumber','MZ','Charge','Intensity') VALUES ({0})", FragToSQLString(glycoPSM, fragment, "YIon"));
                                var command = new SQLiteCommand(commandString, sqlReader);
                                command.ExecuteNonQuery();

                                glycanFragmentIntensity += fragment.fragmentSignal;
                            }

                            var oxoniumIonIntensity = 0.0;
                            foreach(FragmentMatch fragment in oxoniumIons)
                            {
                                var commandString = string.Format("INSERT INTO FragmentIons('File','Condition','Replicate','ScanNumber','Peptide','Sequence','Glycan','ByonicScore','DeltaModScore','LogProb','IonType','Fragment','FragmentNumber','MZ','Charge','Intensity') VALUES ({0})", FragToSQLString(glycoPSM, fragment, "Oxonium"));
                                var command = new SQLiteCommand(commandString, sqlReader);
                                command.ExecuteNonQuery();

                                oxoniumIonIntensity += fragment.fragmentSignal;
                            }

                            var ticExplained = ((pepFragmentIntensity + glycanFragmentIntensity + oxoniumIonIntensity) / TIC) * 100;
                            var insertPepQuery = string.Format("INSERT INTO AnnotatedPeptides('File','Condition','Replicate','Peptide','Sequence','Glycans','ScanNumber','Fragmentation','MZ','Charge','PeptideMass','GlycanMasses','Score','DeltaModScore','FDR2D','PEP2D','SeqCoverage','GlycanSeqCoverage','PercentTICExplained','NumBackboneFragments','NumYIons','NumOxoniumIons') VALUES ({0})", peptideToSQLString(glycoPSM, sequenceCoverage, glycanSequenceCoverage, peptideFragments.Count(), peptideIntactGlycans.Count(), oxoniumIons.Count(), ticExplained));
                            var insertPepCommand = new SQLiteCommand(insertPepQuery, sqlReader);
                            insertPepCommand.ExecuteNonQuery();

                        }

                        if(useType.Equals("SingleSpectrum"))
                        {
                            var returnData = new FragmentDataReturnArgs(peptideFragments, peptideNeutralLossFragments, peptideIntactGlycans, peptideFragmentsMustIncludeGlycan, oxoniumIons, rawFile.GetSpectrum(glycoPSM.scanNumber), glycoPSM);
                            onReturnData(returnData);
                            //onReturnData(peptideFragments, peptideNeutralLossFragments, peptideIntactGlycans, peptideFragmentsMustIncludeGlycan, oxoniumIons, rawFile.GetSpectrum(glycoPSM.scanNumber), glycoPSM);
                        }                
                    }
                }

                if (useType.Equals("BatchAnnotation"))
                {
                    watch.Stop();
                    var time = Math.Round((watch.ElapsedMilliseconds / 60000.0), 2).ToString() + " minutes.";
                    //outputEachFragment.Close();
                    var updateReturnData = new FragmentDataReturnArgs(Path.GetFileName(rawfilePath), time);
                    onReturnData(updateReturnData);                    
                }

                rawFile.Dispose();
            }catch(Exception e)
            {
                var error = e.StackTrace;
            }

            
        }

        public string peptideToSQLString(GlycoPSM psm, double seqCoverage, double glySeqCoverage, int numFrags, int numYions, int numOxoniumIons, double ticExplained)
        {
            return string.Format("'{0}','{1}','{2}','{3}','{4}','{5}','{6}','{7}','{8}','{9}','{10}','{11}','{12}','{13}','{14}','{15}','{16}','{17}','{18}','{19}','{20}','{21}'",
                                    psm.file, psm.Condition, psm.Replicate, psm.peptide.ToString(), psm.peptide.Sequence, psm.glycansToString(),
                                    psm.scanNumber, psm.fragmentation, psm.mzObs, psm.charge, psm.peptideMonoMass, String.Join(";", psm.glycanMasses),
                                    psm.score, psm.deltaModScore, psm.FDR2D, psm.PEP2D, seqCoverage, glySeqCoverage, ticExplained, numFrags,
                                    numYions, numOxoniumIons);
        }

        public string FragToSQLString(GlycoPSM psm, FragmentMatch frag, string type)
        {
            return string.Format("'{0}','{1}','{2}','{3}','{4}','{5}','{6}','{7}','{8}','{9}','{10}','{11}','{12}','{13}','{14}','{15}'",
                                 psm.file, psm.Condition, psm.Replicate,psm.scanNumber, psm.peptide.ToString(), psm.peptide.Sequence, 
                                 psm.glycansToString(), psm.score, psm.deltaModScore, psm.logProb, type, frag.fragmentName, 
                                 frag.fragmentNumber, frag.fragmentMZ, frag.fragmentCharge, frag.fragmentSignal);
        }

        public static string GetModName(string mod)
        {
            string[] modParse1 = mod.Split('(');
            string[] modParse2 = modParse1[1].Split(')');
            string[] modParse3 = modParse2[0].Split('/');
            string[] modParse4 = modParse3[0].Split(' ');

            string modName = modParse4[0];

            return modName;
        }

        public static double GetModMass(string mod)
        {
            string[] modParse1 = mod.Split('(');
            string[] modParse2 = modParse1[1].Split(')');
            string[] modParse3 = modParse2[0].Split('/');

            string massSubstring = modParse3[1].Substring(1);

            double modMass = Convert.ToDouble(massSubstring);

            return modMass;
        }

        public static int GetModPosition(string mod)
        {
            string[] modParse1 = mod.Split('(');

            string modPositionString = modParse1[0].Substring(1);

            int modPosition = Convert.ToInt32(modPositionString);

            return modPosition;
        }

        public static int GetSugarCount(string glycan, string sugar)
        {
            int sugarCount = 0;

            string[] sugarsArray = glycan.Split(')');

            for (int i = 0; i < sugarsArray.Length; i++)
            {
                string[] glycanSep = sugarsArray[i].Split('(');
                if (glycanSep[0].Equals(sugar))
                {
                    sugarCount = Convert.ToInt32(glycanSep[1]);
                }
            }

            return sugarCount;
        }

        //Changed to most intense
        public static ThermoMzPeak GetClosestPeak(List<ThermoMzPeak> peaks, double theoMZ)
        {
            return peaks.OrderByDescending(x => x.Intensity).ToList()[0];
            double diff = double.MaxValue;
            ThermoMzPeak returnPeak = null;
            foreach (var peak in peaks)
            {
                var currDiff = Math.Abs(peak.MZ - theoMZ);
                if (currDiff < diff)
                {
                    diff = currDiff;
                    returnPeak = peak;
                }
            }
            return returnPeak;
        }

        public static bool NoPrecursorPeakInterference(double mzOfFragment, Peptide peptide, int precursorChargeState, string fragmentation)
        {
            bool noInterference = false;

            bool noInterferenceIntactPrecursor = false;
            bool noInterferenceCRP = false;

            double intactPrecursorMZ = peptide.ToMz(precursorChargeState);

            if (mzOfFragment > (intactPrecursorMZ + 1) || mzOfFragment < (intactPrecursorMZ - 1))
            {
                noInterferenceIntactPrecursor = true;
            }

            if (fragmentation.Equals("ETD"))
            {
                for (int i = precursorChargeState - 1; i > 0; i--)
                {
                    double crp = (peptide.MonoisotopicMass + (Constants.Hydrogen * precursorChargeState)) / ((double)i);
                    if (mzOfFragment > (crp + 0.5) || mzOfFragment < (crp - 0.5))
                    {
                        noInterferenceCRP = true;
                    }
                    else
                    {
                        noInterferenceCRP = false;
                    }
                }
            }

            if (fragmentation.Equals("HCD"))
            {
                if (noInterferenceIntactPrecursor)
                {
                    noInterference = true;
                }
            }
            else if (fragmentation.Equals("ETD"))
            {
                if (noInterferenceIntactPrecursor && noInterferenceCRP)
                {
                    noInterference = true;
                }
            }

            return noInterference;
        }

        public static void LookForFragmentWithNoGlycan(ThermoRawFile rawFile, Fragment theoFragment, int charge, int glycanPosition, double glycanMass, 
                                                       List<FragmentMatch> peptideFragments, List<FragmentMatch> peptideNeutralLossFragments, int scanNumber, 
                                                       Peptide peptide, string fragmentation)
        {
            if (theoFragment.Type.Equals(FragmentTypes.b) || theoFragment.Type.Equals(FragmentTypes.c))
            {
                if (theoFragment.Number >= glycanPosition)
                {
                    for (int i = charge - 1; i > 0; i--)
                    {
                        double theoFragmentMZnoGlycan = (theoFragment.MonoisotopicMass + (Constants.Hydrogen * i) - glycanMass) / ((double)i);

                        if (NoPrecursorPeakInterference(theoFragmentMZnoGlycan, peptide, charge, fragmentation))
                        {
                            double theoFragmentIsoMZnoGlycan = (theoFragmentMZnoGlycan + (1 * (Constants.Hydrogen)) / ((double)i));
                            var range = DoubleRange.FromPPM(theoFragmentMZnoGlycan, 20);

                            List<ThermoMzPeak> outPeaks = new List<ThermoMzPeak>();

                            if (rawFile.GetSpectrum(scanNumber).TryGetPeaks(range, out outPeaks))
                            {
                                var closestPeak = GetClosestPeak(outPeaks, theoFragmentMZnoGlycan);
                                double intensity = closestPeak.Intensity;

                                var rangeForIsotope = DoubleRange.FromPPM(theoFragmentIsoMZnoGlycan, 20);
                                List<ThermoMzPeak> outPeaksIso = new List<ThermoMzPeak>();
                                if (rawFile.GetSpectrum(scanNumber).TryGetPeaks(rangeForIsotope, out outPeaksIso))
                                {
                                    var closestPeakIso = GetClosestPeak(outPeaksIso, theoFragmentIsoMZnoGlycan);
                                    intensity += closestPeakIso.Intensity;
                                }

                                string fragmentName = "~" + theoFragment.ToString();

                                double basePeak = rawFile.GetSpectrum(scanNumber).GetBasePeakIntensity();

                                if (closestPeak.Charge == i)
                                {
                                    FragmentMatch fragmentMatch = new FragmentMatch(fragmentName, theoFragment.Type.ToString(), theoFragment.Number, i, closestPeak.MZ, intensity);
                                    peptideFragments.Add(fragmentMatch);
                                    peptideNeutralLossFragments.Add(fragmentMatch);
                                }
                                else if (rawFile.GetSpectrum(scanNumber).TryGetPeaks(rangeForIsotope, out outPeaksIso))
                                {
                                    FragmentMatch fragmentMatch = new FragmentMatch(fragmentName, theoFragment.Type.ToString(), theoFragment.Number, i, closestPeak.MZ, intensity);
                                    peptideFragments.Add(fragmentMatch);
                                    peptideNeutralLossFragments.Add(fragmentMatch);
                                }
                                else if (intensity >= (0.01 * basePeak) && closestPeak.Charge == 0)
                                {
                                    FragmentMatch fragmentMatch = new FragmentMatch(fragmentName, theoFragment.Type.ToString(), theoFragment.Number, i, closestPeak.MZ, intensity);
                                    peptideFragments.Add(fragmentMatch);
                                    peptideNeutralLossFragments.Add(fragmentMatch);
                                }
                            }
                        }

                        //Also look for 1 HexNAc still attached
                        ChemicalFormula hexNAc = new ChemicalFormula("C8O5NH13");

                        double theoFragmentMZplusHexNAc = (theoFragment.MonoisotopicMass + (Constants.Hydrogen * i) - glycanMass + hexNAc.MonoisotopicMass) / ((double)i);

                        if (NoPrecursorPeakInterference(theoFragmentMZplusHexNAc, peptide, charge, fragmentation))
                        {
                            double theoFragmentIsoMZplusHexNAc = (theoFragmentMZplusHexNAc + (1 * (Constants.Hydrogen)) / ((double)i));
                            var rangeHexNAcFragment = DoubleRange.FromPPM(theoFragmentMZplusHexNAc, 20);

                            List<ThermoMzPeak> outPeaksHexNAc = new List<ThermoMzPeak>();

                            if (rawFile.GetSpectrum(scanNumber).TryGetPeaks(rangeHexNAcFragment, out outPeaksHexNAc))
                            {
                                var closestPeak = GetClosestPeak(outPeaksHexNAc, theoFragmentMZplusHexNAc);
                                double intensity = closestPeak.Intensity;

                                var rangeForIsotopeHexNAc = DoubleRange.FromPPM(theoFragmentIsoMZplusHexNAc, 20);
                                List<ThermoMzPeak> outPeaksIsoHexNAc = new List<ThermoMzPeak>();
                                if (rawFile.GetSpectrum(scanNumber).TryGetPeaks(rangeForIsotopeHexNAc, out outPeaksIsoHexNAc))
                                {
                                    var closestPeakIso = GetClosestPeak(outPeaksIsoHexNAc, theoFragmentIsoMZplusHexNAc);
                                    intensity += closestPeakIso.Intensity;
                                }

                                string fragmentName = "~" + theoFragment.ToString() + "+203";

                                double basePeak = rawFile.GetSpectrum(scanNumber).GetBasePeakIntensity();

                                if (closestPeak.Charge == i)
                                {
                                    FragmentMatch fragmentMatch = new FragmentMatch(fragmentName, theoFragment.Type.ToString(), theoFragment.Number, i, closestPeak.MZ, intensity);
                                    peptideFragments.Add(fragmentMatch);
                                    peptideNeutralLossFragments.Add(fragmentMatch);
                                }
                                else if (rawFile.GetSpectrum(scanNumber).TryGetPeaks(rangeForIsotopeHexNAc, out outPeaksIsoHexNAc))
                                {
                                    FragmentMatch fragmentMatch = new FragmentMatch(fragmentName, theoFragment.Type.ToString(), theoFragment.Number, i, closestPeak.MZ, intensity);
                                    peptideFragments.Add(fragmentMatch);
                                    peptideNeutralLossFragments.Add(fragmentMatch);
                                }
                                else if (intensity >= (0.01 * basePeak) && closestPeak.Charge == 0)
                                {
                                    FragmentMatch fragmentMatch = new FragmentMatch(fragmentName, theoFragment.Type.ToString(), theoFragment.Number, i, closestPeak.MZ, intensity);
                                    peptideFragments.Add(fragmentMatch);
                                    peptideNeutralLossFragments.Add(fragmentMatch);
                                }
                            }
                        }
                    }
                }
            }
            else if (theoFragment.Type.Equals(FragmentTypes.y) || theoFragment.Type.Equals(FragmentTypes.zdot))
            {
                if (theoFragment.Number >= peptide.Length - glycanPosition)
                {
                    for (int i = charge - 1; i > 0; i--)
                    {
                        double theoFragmentMZnoGlycan = (theoFragment.MonoisotopicMass + (Constants.Hydrogen * i) - glycanMass) / ((double)i);

                        if (NoPrecursorPeakInterference(theoFragmentMZnoGlycan, peptide, charge, fragmentation))
                        {
                            double theoFragmentIsoMZnoGlycan = (theoFragmentMZnoGlycan + (1 * (Constants.Hydrogen)) / ((double)i));
                            var range = DoubleRange.FromPPM(theoFragmentMZnoGlycan, 20);

                            List<ThermoMzPeak> outPeaks = new List<ThermoMzPeak>();

                            if (rawFile.GetSpectrum(scanNumber).TryGetPeaks(range, out outPeaks))
                            {
                                var closestPeak = GetClosestPeak(outPeaks, theoFragmentMZnoGlycan);
                                double intensity = closestPeak.Intensity;

                                var rangeForIsotope = DoubleRange.FromPPM(theoFragmentIsoMZnoGlycan, 20);
                                List<ThermoMzPeak> outPeaksIso = new List<ThermoMzPeak>();
                                if (rawFile.GetSpectrum(scanNumber).TryGetPeaks(rangeForIsotope, out outPeaksIso))
                                {
                                    var closestPeakIso = GetClosestPeak(outPeaksIso, theoFragmentIsoMZnoGlycan);
                                    intensity += closestPeakIso.Intensity;
                                }

                                string fragmentName = "~" + theoFragment.ToString();

                                double basePeak = rawFile.GetSpectrum(scanNumber).GetBasePeakIntensity();

                                if (closestPeak.Charge == i)
                                {
                                    FragmentMatch fragmentMatch = new FragmentMatch(fragmentName, theoFragment.Type.ToString(), theoFragment.Number, i, closestPeak.MZ, intensity);
                                    peptideFragments.Add(fragmentMatch);
                                    peptideNeutralLossFragments.Add(fragmentMatch);
                                }
                                else if (rawFile.GetSpectrum(scanNumber).TryGetPeaks(rangeForIsotope, out outPeaksIso))
                                {
                                    FragmentMatch fragmentMatch = new FragmentMatch(fragmentName, theoFragment.Type.ToString(), theoFragment.Number, i, closestPeak.MZ, intensity);
                                    peptideFragments.Add(fragmentMatch);
                                    peptideNeutralLossFragments.Add(fragmentMatch);
                                }
                                else if (intensity >= (0.01 * basePeak) && closestPeak.Charge == 0)
                                {
                                    FragmentMatch fragmentMatch = new FragmentMatch(fragmentName, theoFragment.Type.ToString(), theoFragment.Number, i, closestPeak.MZ, intensity);
                                    peptideFragments.Add(fragmentMatch);
                                    peptideNeutralLossFragments.Add(fragmentMatch);
                                }
                            }
                        }

                        //Also look for 1 HexNAc still attached
                        ChemicalFormula hexNAc = new ChemicalFormula("C8O5NH13");

                        double theoFragmentMZplusHexNAc = (theoFragment.MonoisotopicMass + (Constants.Hydrogen * i) - glycanMass + hexNAc.MonoisotopicMass) / ((double)i);

                        if (NoPrecursorPeakInterference(theoFragmentMZplusHexNAc, peptide, charge, fragmentation))
                        {
                            double theoFragmentIsoMZplusHexNAc = (theoFragmentMZplusHexNAc + (1 * (Constants.Hydrogen)) / ((double)i));
                            var rangeHexNAcFragment = DoubleRange.FromPPM(theoFragmentMZplusHexNAc, 20);

                            List<ThermoMzPeak> outPeaksHexNAc = new List<ThermoMzPeak>();

                            if (rawFile.GetSpectrum(scanNumber).TryGetPeaks(rangeHexNAcFragment, out outPeaksHexNAc))
                            {
                                var closestPeak = GetClosestPeak(outPeaksHexNAc, theoFragmentMZplusHexNAc);
                                double intensity = closestPeak.Intensity;

                                var rangeForIsotopeHexNAc = DoubleRange.FromPPM(theoFragmentIsoMZplusHexNAc, 20);
                                List<ThermoMzPeak> outPeaksIsoHexNAc = new List<ThermoMzPeak>();
                                if (rawFile.GetSpectrum(scanNumber).TryGetPeaks(rangeForIsotopeHexNAc, out outPeaksIsoHexNAc))
                                {
                                    var closestPeakIso = GetClosestPeak(outPeaksIsoHexNAc, theoFragmentIsoMZplusHexNAc);
                                    intensity += closestPeakIso.Intensity;
                                }

                                string fragmentName = "~" + theoFragment.ToString() + "+203";

                                double basePeak = rawFile.GetSpectrum(scanNumber).GetBasePeakIntensity();

                                if (closestPeak.Charge == i)
                                {
                                    FragmentMatch fragmentMatch = new FragmentMatch(fragmentName, theoFragment.Type.ToString(), theoFragment.Number, i, closestPeak.MZ, intensity);
                                    peptideFragments.Add(fragmentMatch);
                                    peptideNeutralLossFragments.Add(fragmentMatch);
                                }
                                else if (rawFile.GetSpectrum(scanNumber).TryGetPeaks(rangeForIsotopeHexNAc, out outPeaksIsoHexNAc))
                                {
                                    FragmentMatch fragmentMatch = new FragmentMatch(fragmentName, theoFragment.Type.ToString(), theoFragment.Number, i, closestPeak.MZ, intensity);
                                    peptideFragments.Add(fragmentMatch);
                                    peptideNeutralLossFragments.Add(fragmentMatch);
                                }
                                else if (intensity >= (0.01 * basePeak) && closestPeak.Charge == 0)
                                {
                                    FragmentMatch fragmentMatch = new FragmentMatch(fragmentName, theoFragment.Type.ToString(), theoFragment.Number, i, closestPeak.MZ, intensity);
                                    peptideFragments.Add(fragmentMatch);
                                    peptideNeutralLossFragments.Add(fragmentMatch);
                                }
                            }
                        }

                    }
                }
            }
        }

        public static void GetIntactPepGlycanFragments(GlycoPSM glycoPSM, int charge, int hexCount, int hexNAcCount, int fucCount, int neuAcCount, int neuGcCount, List<FragmentMatch> peptideIntactGlycanFragments, ThermoSpectrum spectrum, string glycanMod, double glycanMass)
        {
            Glycan glycan = new Glycan(glycanMod);

            int totalSugarIntactGlycan = hexNAcCount + hexCount + fucCount + neuAcCount + neuGcCount;

            glycan.GenerateCombinations();

            var glycanPieces = glycan.AllFragments;

            var filteredPieces = GlycanMethods.GetValidGlycanStructures(glycanPieces);
            foreach (var piece in filteredPieces)
            {
                for (int i = charge; i > 0; i--)
                {
                    double peptidePlusGlycan_MZ = ((glycoPSM.peptide.MonoisotopicMass - glycanMass) + (Constants.Hydrogen * i) + GlycanMethods.GetGlycanMass(piece.Value)) / ((double)i);
                    if (peptidePlusGlycan_MZ < 3000)
                    {
                        double peptidePlusGlycan_MZiso1 = (peptidePlusGlycan_MZ + (1 * Constants.Hydrogen)) / ((double)i);
                        double peptidePlusGlycan_MZiso2 = (peptidePlusGlycan_MZ + (2 * Constants.Hydrogen)) / ((double)i);
                        double peptidePlusGlycan_MZiso3 = (peptidePlusGlycan_MZ + (3 * Constants.Hydrogen)) / ((double)i);

                        DoubleRange range = DoubleRange.FromPPM(peptidePlusGlycan_MZ, 40);
                        List<ThermoMzPeak> outPeaks = new List<ThermoMzPeak>();
                        if (spectrum.TryGetPeaks(range, out outPeaks))
                        {
                            var closestPeak = GetClosestPeak(outPeaks, peptidePlusGlycan_MZ);
                            double intensity = closestPeak.Intensity;

                            DoubleRange rangeIso1 = DoubleRange.FromPPM(peptidePlusGlycan_MZiso1, 40);
                            List<ThermoMzPeak> outpeaksIso1 = new List<ThermoMzPeak>();
                            if (spectrum.TryGetPeaks(rangeIso1, out outpeaksIso1))
                            {
                                var closestPeakIso1 = GetClosestPeak(outpeaksIso1, peptidePlusGlycan_MZiso1);
                                intensity += closestPeakIso1.Intensity;
                            }

                            DoubleRange rangeIso2 = DoubleRange.FromPPM(peptidePlusGlycan_MZiso2, 40);
                            List<ThermoMzPeak> outpeaksIso2 = new List<ThermoMzPeak>();
                            if (spectrum.TryGetPeaks(rangeIso2, out outpeaksIso2))
                            {
                                var closestPeakIso2 = GetClosestPeak(outpeaksIso2, peptidePlusGlycan_MZiso2);
                                intensity += closestPeakIso2.Intensity;
                            }

                            DoubleRange rangeIso3 = DoubleRange.FromPPM(peptidePlusGlycan_MZiso3, 40);
                            List<ThermoMzPeak> outpeaksIso3 = new List<ThermoMzPeak>();
                            if (spectrum.TryGetPeaks(rangeIso3, out outpeaksIso3))
                            {
                                var closestPeakIso3 = GetClosestPeak(outpeaksIso3, peptidePlusGlycan_MZiso3);
                                intensity += closestPeakIso3.Intensity;
                            }

                            double basePeak = spectrum.GetBasePeakIntensity();

                            string fragmentName = "Pep+" + piece.Key;

                            int sugarCountThisFragment = 0;
                            sugarCountThisFragment += GetSugarCount(piece.Key, "Hex");
                            sugarCountThisFragment += GetSugarCount(piece.Key, "HexNAc");
                            sugarCountThisFragment += GetSugarCount(piece.Key, "Fuc");
                            sugarCountThisFragment += GetSugarCount(piece.Key, "NeuAc");
                            sugarCountThisFragment += GetSugarCount(piece.Key, "NeuGc");

                            int sugarFragmentNumber = totalSugarIntactGlycan - sugarCountThisFragment;

                            if (closestPeak.Charge == i)
                            {
                                FragmentMatch fragmentMatch = new FragmentMatch(fragmentName, "Glyco", sugarFragmentNumber, i, peptidePlusGlycan_MZ, intensity);
                                peptideIntactGlycanFragments.Add(fragmentMatch);
                            }
                            else if (spectrum.TryGetPeaks(rangeIso1, out outpeaksIso1))
                            {
                                FragmentMatch fragmentMatch = new FragmentMatch(fragmentName, "Glyco", sugarFragmentNumber, i, peptidePlusGlycan_MZ, intensity);
                                peptideIntactGlycanFragments.Add(fragmentMatch);
                            }
                            else if (intensity > (basePeak * 0.01) && closestPeak.Charge == 0)
                            {
                                FragmentMatch fragmentMatch = new FragmentMatch(fragmentName, "Glyco", sugarFragmentNumber, i, peptidePlusGlycan_MZ, intensity);
                                peptideIntactGlycanFragments.Add(fragmentMatch);
                            }
                        }
                    }
                }
            }


        }

        public static void GetOxoniumIons(GlycoPSM glycoPSM, int hexCount, int hexNAcCount, int fucCount, 
                                          int neuAcCount, int neuGcCount, string glycanMod, double glycanMass, 
                                          ThermoSpectrum spectrum, List<FragmentMatch> oxoniumIons)
        {
            Glycan glycan = new Glycan(glycanMod);

            int totalSugarIntactGlycan = hexNAcCount + hexCount + fucCount + neuAcCount + neuGcCount;

            glycan.GenerateCombinations();

            var glycanPieces = glycan.AllFragments;

            var filteredPieces = GlycanMethods.GetValidGlycanStructures(glycanPieces);

            foreach (var piece in filteredPieces)
            {
                //Console.WriteLine(piece.Key + GlycanMethods.GetGlycanMZ(piece.Value, 1));
                FindOxoniumIon(GlycanMethods.GetGlycanMZ(piece.Value, 1), spectrum, oxoniumIons, glycoPSM.scanNumber, piece.Key);
            }

            FindOxoniumIon(366.14, spectrum, oxoniumIons, glycoPSM.scanNumber, "Hex(1)HexNAc(1)");
            FindOxoniumIon(528.19, spectrum, oxoniumIons, glycoPSM.scanNumber, "Hex(2)HexNAc(1)");
            FindOxoniumIon(690.245, spectrum, oxoniumIons, glycoPSM.scanNumber, "Hex(3)HexNAc(1)");

            if (hexCount > 0)
            {
                FindOxoniumIon(163.06, spectrum, oxoniumIons, glycoPSM.scanNumber, "Hex(1)");
            }
            if (fucCount > 0)
            {
                FindOxoniumIon(147.065, spectrum, oxoniumIons, glycoPSM.scanNumber, "Fuc(1)");
            }
            if (neuAcCount > 0)
            {
                FindOxoniumIon(274.092, spectrum, oxoniumIons, glycoPSM.scanNumber, "NeuAcWaterLoss");
                FindOxoniumIon(292.103, spectrum, oxoniumIons, glycoPSM.scanNumber, "NeuAc(1)");
                FindOxoniumIon(495.1382, spectrum, oxoniumIons, glycoPSM.scanNumber, "HexNAc(1)NeuAc(1)");
                FindOxoniumIon(657.23, spectrum, oxoniumIons, glycoPSM.scanNumber, "Hex(1)HexNAc(1)NeuAc(1)");

            }
            if (neuGcCount > 0)
            {
                FindOxoniumIon(290.087, spectrum, oxoniumIons, glycoPSM.scanNumber, "NeuGcWaterLoss");
                FindOxoniumIon(308.098, spectrum, oxoniumIons, glycoPSM.scanNumber, "NeuGc(1)");
            }

        }

        public static void FindOxoniumIon(double mz, ThermoSpectrum spectrum, List<FragmentMatch> list, int scanNumber, string name)
        {
            double mzIso = mz + Constants.Hydrogen;

            DoubleRange mzRange = DoubleRange.FromPPM(mz, 20);
            List<ThermoMzPeak> outPeaks = new List<ThermoMzPeak>();

            if (spectrum.TryGetPeaks(mzRange, out outPeaks))
            {
                var closestPeak = GetClosestPeak(outPeaks, mz);
                double intensity = closestPeak.Intensity;

                DoubleRange isoMZrange = DoubleRange.FromPPM(mzIso, 20);
                List<ThermoMzPeak> outPeaksIso = new List<ThermoMzPeak>();
                if (spectrum.TryGetPeaks(isoMZrange, out outPeaksIso))
                {
                    var closestPeakIso = GetClosestPeak(outPeaksIso, mzIso);
                    intensity += closestPeakIso.Intensity;
                }

                double basePeak = spectrum.GetBasePeakIntensity();

                if (intensity > (0.01 * basePeak))
                {
                    FragmentMatch match = new FragmentMatch(name, "Oxonium", 0, 1, mz, intensity);
                    list.Add(match);
                }
            }

        }

        public static double GetSequenceCoverage(GlycoPSM glycoPSM, List<FragmentMatch> fragments)
        {
            double sequenceCoverage = 0;

            char[] peptidePostitionArrary = new char[glycoPSM.peptide.Length - 1];

            foreach (FragmentMatch fragment in fragments)
            {
                if (fragment.fragmentType.Equals("b") || fragment.fragmentType.Equals("c"))
                {
                    int bondNumberBroken = fragment.fragmentNumber;
                    int positionInArray = bondNumberBroken - 1;
                    peptidePostitionArrary[positionInArray] = 'X';
                }
                if (fragment.fragmentType.Equals("y") || fragment.fragmentType.Equals("zdot"))
                {
                    int bondNumberBroken = glycoPSM.peptide.Length - fragment.fragmentNumber;
                    int positionInArray = bondNumberBroken - 1;
                    peptidePostitionArrary[positionInArray] = 'X';
                }
            }

            int numberOfBondsBroken = 0;
            for (int i = 0; i < peptidePostitionArrary.Length; i++)
            {
                if (peptidePostitionArrary[i].Equals('X'))
                    numberOfBondsBroken++;
            }
            sequenceCoverage = (double)numberOfBondsBroken / (double)peptidePostitionArrary.Length;

            return sequenceCoverage;
        }

        public static double GetGlycanSequenceCoverage(GlycoPSM glycoPSM, List<FragmentMatch> fragments, int hexCount, int hexNAcCount, int fucCount, int neuAcCount, int neuGcCount)
        {
            double sequenceCoverage = 0;
            Glycan glycan = new Glycan(glycoPSM.glycans[0]);

            glycan.GenerateCombinations();

            var glycanPieces = glycan.AllFragments;

            var filteredPieces = GlycanMethods.GetValidGlycanStructures(glycanPieces);

            HashSet<string> theorecticalFragments = new HashSet<string>();

            //Console.WriteLine(glycoPSM.glycans[0]);
            foreach (var piece in filteredPieces)
            {
                string fragment = "Pep+" + piece.Key;
                theorecticalFragments.Add(fragment);
            }

            HashSet<string> matches = new HashSet<string>();

            foreach (FragmentMatch match in fragments)
            {
                string[] glycanFrag = match.fragmentName.Split('_');
                matches.Add(glycanFrag[0]);
            }

            int numberOfGlycoFragmentsExplained = 0;
            int numberOfGlycoFragmentsTotal = theorecticalFragments.Count;

            foreach (string theoFrag in theorecticalFragments)
            {
                if (matches.Contains(theoFrag))
                {
                    numberOfGlycoFragmentsExplained++;
                }
            }

            sequenceCoverage = (double)numberOfGlycoFragmentsExplained / (double)numberOfGlycoFragmentsTotal;

            return sequenceCoverage;
        }

        public static bool SeeIfFragmentsHaveGlycan(GlycoPSM glycoPSM, List<FragmentMatch> fragments, string terminus)
        {
            bool returnBoolAll = false;

            bool returnBoolNterm = false;
            bool returnBoolCterm = false;



            int glycanPosition = glycoPSM.glycanPositions[0];
            //Console.WriteLine(glycoPSM.peptide);
            //Console.WriteLine(glycanPosition);

            foreach (FragmentMatch fragment in fragments)
            {
                //Console.WriteLine(fragment.fragmentName);
                if (fragment.fragmentType.Equals("b") || fragment.fragmentType.Equals("c"))
                {
                    if (fragment.fragmentNumber >= glycanPosition)
                    {
                        returnBoolAll = true;
                        returnBoolNterm = true;
                    }
                }

                if (fragment.fragmentType.Equals("y") || fragment.fragmentType.Equals("zdot"))
                {
                    if (fragment.fragmentNumber > (glycoPSM.peptide.Length - glycanPosition))
                    {
                        returnBoolAll = true;
                        returnBoolCterm = true;
                    }
                }
            }

            //Console.WriteLine(returnBoolAll);
            //Console.WriteLine(returnBoolNterm);
            //Console.WriteLine(returnBoolCterm);
            //Console.ReadKey();
            if (terminus.Equals("NTerm"))
            {
                return returnBoolNterm;
            }
            else if (terminus.Equals("Cterm"))
            {
                return returnBoolCterm;
            }
            else
            {
                return returnBoolAll;
            }
        }

        public EventHandler<FragmentDataReturnArgs> ReturnData;

        //protected virtual void onReturnData(List<FragmentMatch> peptideFragments, List<FragmentMatch> peptideNeutralLossFragments, List<FragmentMatch> YIons, List<FragmentMatch> peptideFragmentsMustIncludeGlycan, List<FragmentMatch> oxoniumIons, ThermoSpectrum spectrum, GlycoPSM glycoPSM, string file, double timeTaken, string returnType)
        protected virtual void onReturnData(FragmentDataReturnArgs args)
        {
            var handler = ReturnData;

            if(handler != null)
            {
                handler(this, args);
            }
        }        
    }
}