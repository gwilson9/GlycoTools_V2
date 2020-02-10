using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using CSMSL.IO.Thermo;
using CSMSL.Proteomics;

namespace _20190618_GlycoTools_V2
{
    public class PSM : IEquatable<PSM>
    {
        //String PID;                                                                      
        public string protRank { get; set; }
        public string PQMsID { get; set; }
        public string sequence { get; set; }
        public string sequenceNoMods { get; set; }
        public string peptidesToBeParsed { get; set; }
        public int peptideStartPosition { get; set; }
        public string modsToBeParsed { get; set; }
        public string glycansToBeParsed { get; set; }
        public double PEP2D { get; set; }
        public double PEP1D { get; set; }
        public double logProb { get; set; }
        public double score { get; set; }
        public double deltaScore { get; set; }
        public double deltaModScore { get; set; }
        public int charge { get; set; }
        public double mzObs { get; set; }
        public double mzCalc { get; set; }
        public double ppmError { get; set; }
        public double obsMH { get; set; }
        public double calcMH { get; set; }
        public string cleavage { get; set; }
        public string glycanPositions { get; set; }
        public string proteinName { get; set; }
        public int protID { get; set; }
        public int scanNumber { get; set; }
        public double intensity { get; set; }
        public string scanTime { get; set; }
        public string modsFixed { get; set; }
        public double FDR2D { get; set; }
        public double FDR1D { get; set; }
        public double FDR2Dunique { get; set; }
        public double FDR1Dunique { get; set; }
        public double qvalue2D { get; set; }
        public double qvalue1D { get; set; }
        public bool isGlycopeptide { get; set; }
        public bool modsPassedCheck { get; set; }
        public bool positionsPassedCheck { get; set; }
        public string evidenceType { get; set; }
        

        // Assigned outside of class
        public string DissociationType { get; set; }
        public int MasterScan { get; set; }
        public string File { get; set; }
        public string Condition { get; set; }
        public string Replicate { get; set; }        
        public string rawFile { get; set; }
        public Peptide peptide { get; set; }
        public ThermoMzPeak Peak126 { get; set; }
        public ThermoMzPeak Peak138 { get; set; }
        public ThermoMzPeak Peak144 { get; set; }
        public ThermoMzPeak Peak168 { get; set; }
        public ThermoMzPeak Peak186 { get; set; }
        public ThermoMzPeak Peak204 { get; set; }
        public ThermoMzPeak Peak274 { get; set; }
        public ThermoMzPeak Peak292 { get; set; }
        public ThermoMzPeak Peak366 { get; set; }
        public List<RTPeak> peakElution { get; set; }
        public List<double> localMaxima { get; set; }
        public List<Glycan> glycans { get; set; }
        public int scanNumberofMaxElutionIntensity { get; set; }
        public bool isSourceFragment { get; set; }
        public bool wasRelocalized { get; set; }


        public PSM()
        {

        }

        public PSM(string PQMsID, string sequence, string peptidesToBeParsed, int peptideStartPosition,
                    double PEP2D, double PEP1D, double score, double deltaScore, double deltaModScore,
                    int charge, double mzObs, double mzCalc, double obsMH, double calcMH, string cleavage,
                    string proteinName, int protID, int scanNumber, double FDR2D, double FDR1D, double FDR2Dunique, 
                    double FDR1Dunique, double qvalue2D, double qvalue1D, double intensity, string scanTime, string protRank)
        {
            this.PQMsID = PQMsID;
            this.sequence = sequence;
            this.peptidesToBeParsed = peptidesToBeParsed;
            this.peptideStartPosition = peptideStartPosition;
            this.PEP2D = PEP2D;
            this.PEP1D = PEP1D;
            this.score = score;
            this.deltaScore = deltaScore;
            this.deltaModScore = deltaModScore;
            this.charge = charge;
            this.mzObs = mzObs;
            this.mzCalc = mzCalc;
            this.obsMH = obsMH;
            this.calcMH = calcMH;
            this.cleavage = cleavage;
            this.proteinName = proteinName;
            this.protID = protID;
            this.scanNumber = scanNumber;
            this.FDR1D = FDR1D;
            this.FDR2D = FDR2D;
            this.FDR1D = FDR1D;
            this.FDR2Dunique = FDR2Dunique;
            this.FDR1Dunique = FDR1Dunique;
            this.qvalue2D = qvalue2D;
            this.qvalue1D = qvalue1D;
            this.intensity = intensity;
            this.scanTime = scanTime;
            this.ppmError = (mzCalc - mzObs) / mzCalc * 1000000;
            this.protRank = protRank;
            this.logProb = 0 - Math.Log10(PEP1D);
            this.evidenceType = "";
            glycans = new List<Glycan>();
            this.sequenceNoMods = peptidesToBeParsed.Split(',')[0];
            this.wasRelocalized = false;
            this.isSourceFragment = false;
        }

        public bool Equals(PSM other)
        {
            return (this.scanNumber == other.scanNumber && this.File.Equals(other.File));
        }

        public override int GetHashCode()
        {
            return this.scanNumber.GetHashCode() * this.File.GetHashCode();
        }

        public void checkMods()
        {
            if (!String.IsNullOrEmpty(modsToBeParsed))
            {
                string[] modsParsedArrary = modsToBeParsed.Split(';');
                for (int i = 0; i < modsParsedArrary.Length; i++)
                {
                    string mod = modsParsedArrary[i];
                 
                    if (mod[0].Equals(' '))
                    {
                        mod = mod.Substring(1);
                    }

                    if (mod.Contains("NGlycan") && !mod[0].Equals('N'))
                    {
                        modsPassedCheck = false;
                    }

                    if (mod.Contains("NGlycan") && mod[0].Equals('N'))
                    {
                        modsPassedCheck = true;
                    }

                    if (mod.Contains("OGlycan") && mod[0].Equals('S'))
                    {
                        modsPassedCheck = true;
                    }
                    if (mod.Contains("OGlycan") && mod[0].Equals('T'))
                    {
                        modsPassedCheck = true;
                    }
                }


                List<int> glycanPositionsList = new List<int>();
                if (!String.IsNullOrEmpty(glycanPositions))
                {
                    string[] glycansPosParsedArray = glycanPositions.Split(';');
                    for (int i = 0; i < glycansPosParsedArray.Length; i++)
                    {
                        int glycanPos = Convert.ToInt32(glycansPosParsedArray[i]);
                        glycanPositionsList.Add(glycanPos);
                    }
                }
                
                if (glycanPositionsList.Count != glycanPositionsList.Distinct().Count())
                {
                    positionsPassedCheck = false;
                }
                else
                {
                    positionsPassedCheck = true;
                }
            }  
        }

        public override string ToString()
        {
            
            var scanTimeInMin = Math.Round(double.Parse(scanTime) / 60, 2);

            var retentionTimes = "";
            var intensities = "";
            foreach(var peak in peakElution)
            {
                retentionTimes += peak.RT.ToString() + ',';
                intensities += peak.Intensity.ToString() + ',';
            }

            var absoluteMaxima = 0.0; //peakElution.Aggregate((i1, i2) => i1.Intensity > i2.Intensity ? i1 : i2).RT;
            var maxIntensity = 0.0; // peakElution.Aggregate((i1, i2) => i1.Intensity > i2.Intensity ? i1 : i2).Intensity;

            var localMaxima = new List<double>();
            for(int i = 1; i < peakElution.Count()-1; i++)
            {
                if (peakElution[i].Intensity > peakElution[i - 1].Intensity && peakElution[i].Intensity > peakElution[i + 1].Intensity && (peakElution[i].Intensity / maxIntensity) > 0.1)
                    localMaxima.Add(peakElution[i].RT);

                if (peakElution[i].Intensity > absoluteMaxima)
                {
                    absoluteMaxima = peakElution[i].RT;
                    maxIntensity = peakElution[i].Intensity;
                }
                    
            }            

            var glycanTypes = string.Join(";", glycans.Select(x => x.glycanType).ToArray());

            var protString = proteinName.Replace("'", "''");
            string returnString = string.Format("'{0}','{1}','{2}','{3}','{4}','{5}','{6}','{7}','{8}','{9}'," +
                                                "'{10}','{11}','{12}','{13}','{14}','{15}','{16}','{17}','{18}'," +
                                                "'{19}','{20}','{21}','{22}','{23}','{24}','{25}','{26}','{27}'," +
                                                "'{28}','{29}','{30}','{31}','{32}','{33}','{34}','{35}','{36}'" +
                                                ",'{37}','{38}','{39}','{40}','{41}','{42}','{43}','{44}', '{45}', '{46}', '{47}'",
                                                sequence, peptidesToBeParsed, sequenceNoMods, peptideStartPosition, modsToBeParsed, glycansToBeParsed,             // 0-5
                                                glycanTypes, PEP2D, PEP1D, logProb, score, deltaScore, deltaModScore, charge, mzObs, mzCalc,                       // 6-15
                                                ppmError, obsMH, calcMH, cleavage, glycanPositions, protString, scanTimeInMin, scanNumber,                         // 16-23
                                                modsFixed, FDR2D, FDR1D, FDR2Dunique, FDR1Dunique, qvalue2D, qvalue1D, isGlycopeptide,                             // 24-31
                                                modsPassedCheck, positionsPassedCheck, DissociationType, MasterScan, intensity,                                    // 32-36
                                                Path.GetFileNameWithoutExtension(File), Condition, Replicate, isSourceFragment.ToString(), wasRelocalized.ToString(),  // 37-41
                                                retentionTimes, intensities, evidenceType, scanNumberofMaxElutionIntensity, absoluteMaxima, string.Join(",", localMaxima));   // 42-47                                                                                                                                                                                             

            return returnString;
        }        
    }
}
