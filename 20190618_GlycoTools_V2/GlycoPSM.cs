using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using CSMSL.Proteomics;
using System.Data.SQLite;
using System.IO;

namespace _20190618_GlycoTools_V2
{
    public class GlycoPSM
    {
        public Peptide peptide { get; set; }
        public double peptideMonoMass { get; set; }
        public string sequenceWithMods { get; set; }
        public List<string> mods { get; set; }
        public List<string> glycans { get; set; }
        public List<double> glycanMasses { get; set; }
        public List<int> glycanPositions { get; set; }
        public string uniprotID { get; set; }
        public double PEP2D { get; set; }
        public double logProb { get; set; }
        public double score { get; set; }
        public double deltaScore { get; set; }
        public double deltaModScore { get; set; }
        public double mzObs { get; set; }
        public int charge { get; set; }
        public int numberOfSites { get; set; }
        public double ppmError { get; set; }
        public double obsMH { get; set; }
        public string cleavage { get; set; }
        public string proteinName { get; set; }
        public int peptideStartPosition { get; set; }
        public double scanTime { get; set; }
        public int scanNumber { get; set; }
        public double FDR2D { get; set; }
        public double FDR2Dunique { get; set; }
        public double qvalue2D { get; set; }
        public string fragmentation { get; set; }
        public bool isGlycoPeptide { get; set; }
        public bool seenWithHCD { get; set; }
        public bool seenWithETD { get; set; }
        public bool NXSmotif { get; set; }
        public bool NXTmotif { get; set; }
        public bool isLocalized { get; set; }
        public bool Nlinked { get; set; }
        public bool Olinked { get; set; }
        public bool matchedToUniprot { get; set; }
        public string evidenceType { get; set; }
        public int evidenceNumber { get; set; }
        public string file;
        public string Condition;
        public string Replicate;
        public string formattedVarMods;

        public GlycoPSM()
        {

        }

        public GlycoPSM(PSM psm, int glycanPosition)
        {
            var sequence = psm.sequence;
            var peptidesToBeParsed = psm.peptidesToBeParsed;
            var peptideStartPosition = psm.peptideStartPosition;
            var modsToBeParsed = psm.modsToBeParsed;
            var glycansToBeParsed = psm.glycansToBeParsed;
            var PEP2D = psm.PEP2D;
            var logProb = psm.logProb;
            var score = psm.score;
            var deltaScore = psm.deltaScore;
            var deltaModScore = psm.deltaModScore;
            var charge = psm.charge;
            var mzObs = psm.mzObs;
            var ppmError = psm.ppmError;
            var obsMH = psm.obsMH;
            var cleavage = psm.cleavage;
            var glycanPositions = psm.glycanPositions;
            var proteinName = psm.proteinName;
            var scanTime = psm.scanTime;
            var scanNumber = psm.scanNumber;
            var modsFixed = psm.modsFixed;
            var FDR2D = psm.FDR2D;
            var FDR2Dunique = psm.FDR2Dunique;
            var qvalue2D = psm.qvalue2D;
            var isGlycoPeptide = true;
            //var test1 = reader["isGlycoPeptide"].ToString(); // 'True' in SQLite Table reads in '0', should be '1'. Not a big issue because all psms are already filtered for glycopeptides at this point
            //var test = reader["isGlycoPeptide"]; //.Equals("True") ? true : false;
            //var fragmentation = psm.DissociationType;

            // Manually set glycan and position here
            List<int> glycanPositionsList = new List<int>() { glycanPosition };
            List<string> glycans = new List<string>() { glycansToBeParsed };

            List<string> mods = new List<string>();
            

            bool seenWithHCD = false;
            bool seenWithETD = false;
            bool NXSmotif = false;
            bool NXTmotif = false;
            bool isLocalized = false;
            bool Nlinked = false;
            bool Olinked = false;

            List<double> glycanMasses = new List<double>();

            bool matchedToUniprot = false;
            string uniprotEvidenceType = "None";
            int uniprotEvidenceNumber = 0;

            string[] parsedPeptide = peptidesToBeParsed.Split(',');
            Peptide peptide = new Peptide(parsedPeptide[0]);
            Peptide peptideNoGlycan = new Peptide(parsedPeptide[0]);
            double peptideMonoMass = peptide.MonoisotopicMass;

            char[] peptideTermini = parsedPeptide[1].ToCharArray();
            char peptideCterminusNextResidue = peptideTermini[1];

            string[] parsedSequenceWithMods = sequence.Split('.');
            string sequenceWithMods = parsedSequenceWithMods[1];

            string[] parsedProteinName = proteinName.Split('|');
            string uniprotID = parsedProteinName[1];

            //if (fragmentation.Equals("HCD"))
            //    seenWithHCD = true;

            //if (fragmentation.Equals("ETD"))
            //    seenWithETD = true;


            int numberOfSites = glycans.Count;

            if (!String.IsNullOrEmpty(modsToBeParsed))
            {
                string[] modsParsedArray = modsToBeParsed.Split(';');
                for (int i = 0; i < modsParsedArray.Length; i++)
                {
                    string mod = modsParsedArray[i];
                    if (mod[0].Equals(' '))
                    {
                        mod = mod.Substring(1);
                    }

                    mods.Add(mod);
                    string modName = GetModName(mod);
                    double modMass = GetModMass(mod);
                    int modPosition = GetModPosition(mod);
                    if (modName.Contains("Glycan"))
                    {
                        if (!peptide.Sequence[glycanPosition - 1].Equals('N'))
                        {
                            modName = "OGlycan";                            
                        }
                        else
                        {
                            modName = "NGlycan";
                        }

                        var newVarMod = peptide.Sequence[glycanPosition-1].ToString() + glycanPosition.ToString() + '(' + modName + " / " + modMass + ')';
                        formattedVarMods += ";" + newVarMod;
                        formattedVarMods = formattedVarMods.Trim(';');
                        modName = modName + "_" + modMass;
                        glycanMasses.Add(modMass);
                        modPosition = glycanPosition;
                    }
                    else
                    {
                        formattedVarMods += ";" + mod;
                        formattedVarMods = formattedVarMods.Trim(';');
                    }

                    if (modName.Contains("NGlycan"))
                        Nlinked = true;

                    if (modName.Contains("OGlycan"))
                        Olinked = true;

                    Modification modToAdd = new Modification(modMass, modName);
                    peptide.AddModification(modToAdd, modPosition);

                    if (modName.Contains("NGlycan"))
                    {
                        if ((modPosition + 2) > peptide.Length)
                        {
                            if (!peptide.GetResidue(peptide.Length - 1).Equals('P'))
                            {
                                if (peptideCterminusNextResidue.Equals('S'))
                                {
                                    NXSmotif = true;
                                }
                                if (peptideCterminusNextResidue.Equals('T'))
                                {
                                    NXTmotif = true;
                                }
                            }
                        }
                        else
                        {
                            if (!peptide.GetResidue(modPosition).Letter.Equals('P'))
                            {
                                if (peptide.GetResidue(modPosition + 1).Letter.Equals('S'))
                                {
                                    NXSmotif = true;
                                }
                                if (peptide.GetResidue(modPosition + 1).Letter.Equals('T'))
                                {
                                    NXTmotif = true;
                                }
                            }
                        }
                    }
                }
            }

            if (!String.IsNullOrEmpty(modsFixed))
            {
                string[] modsParsedArray = modsFixed.Split(';');
                for (int i = 0; i < modsParsedArray.Length; i++)
                {
                    string mod = modsParsedArray[i];
                    if (mod[0].Equals(' '))
                    {
                        mod = mod.Substring(1);
                    }
                    mods.Add(mod);
                    string modName = GetModName(mod);
                    double modMass = GetModMass(mod);
                    int modPosition = GetModPosition(mod);

                    Modification modToAdd = new Modification(modMass, modName);
                    peptide.AddModification(modToAdd, modPosition);
                }
            }

            if (deltaModScore >= 10)
                isLocalized = true;

            this.peptide = peptide;
            this.peptideMonoMass = peptideMonoMass;
            this.sequenceWithMods = peptide.SequenceWithModifications;
            this.mods = mods;
            this.glycans = glycans;
            this.glycanMasses = glycanMasses;
            this.glycanPositions = glycanPositionsList;
            this.uniprotID = uniprotID;
            this.PEP2D = PEP2D;
            this.logProb = logProb;
            this.score = score;
            this.deltaScore = deltaScore;
            this.deltaModScore = deltaModScore;
            this.mzObs = mzObs;
            this.charge = charge;
            this.numberOfSites = numberOfSites;
            this.ppmError = ppmError;
            this.obsMH = obsMH;
            this.cleavage = cleavage;
            this.proteinName = proteinName;
            this.peptideStartPosition = peptideStartPosition;
            this.scanTime = double.Parse(scanTime);
            this.scanNumber = scanNumber;
            this.FDR2D = FDR2D;
            this.FDR2Dunique = FDR2Dunique;
            this.qvalue2D = qvalue2D;
            //this.fragmentation = fragmentation;
            this.isGlycoPeptide = isGlycoPeptide;
            this.seenWithHCD = seenWithHCD;
            this.seenWithETD = seenWithETD;
            this.NXSmotif = NXSmotif;
            this.NXTmotif = NXTmotif;
            this.isLocalized = isLocalized;
            this.Nlinked = Nlinked;
            this.Olinked = Olinked;
            this.matchedToUniprot = matchedToUniprot;
            this.evidenceType = evidenceType;
            this.evidenceNumber = evidenceNumber;
            this.file = psm.File;
            this.Condition = psm.Condition;
            this.Replicate = psm.Replicate;

        }

        public static int GetModPosition(string mod)
        {
            string[] modParse1 = mod.Split('(');

            string modPositionString = modParse1[0].Substring(1);

            int modPosition = Convert.ToInt32(modPositionString);

            return modPosition;
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

        public GlycoPSM(Peptide peptide, double peptideMonoMass, string sequenceWithMods, List<string> mods, List<string> glycans, List<double> glycanMasses, List<int> glycanPositions, string uniprotID, double PEP2D, double logProb, double score,
            double deltaScore, double deltaModScore, double mzObs, int charge, int numberOfSites, double ppmError, double obsMH, string cleavage, string proteinName,
            int peptideStartPosition, double scanTime, int scanNumber, double FDR2D, double FDR2Dunique, double qvalue2D, string fragmentation, bool isGlycoPeptide,
            bool seenWithHCD, bool seenWithETD, bool NXSmotif, bool NXTmotif, bool isLocalized, bool Nlinked, bool Olinked, bool matchedToUniprot, string evidenceType, int evidenceNumber)
        {
            this.peptide = peptide;
            this.peptideMonoMass = peptideMonoMass;
            this.sequenceWithMods = peptide.SequenceWithModifications;
            this.mods = mods;
            this.glycans = glycans;
            this.glycanMasses = glycanMasses;
            this.glycanPositions = glycanPositions;
            this.uniprotID = uniprotID;
            this.PEP2D = PEP2D;
            this.logProb = logProb;
            this.score = score;
            this.deltaScore = deltaScore;
            this.deltaModScore = deltaModScore;
            this.mzObs = mzObs;
            this.charge = charge;
            this.numberOfSites = numberOfSites;
            this.ppmError = ppmError;
            this.obsMH = obsMH;
            this.cleavage = cleavage;
            this.proteinName = proteinName;
            this.peptideStartPosition = peptideStartPosition;
            this.scanTime = scanTime;
            this.scanNumber = scanNumber;
            this.FDR2D = FDR2D;
            this.FDR2Dunique = FDR2Dunique;
            this.qvalue2D = qvalue2D;
            this.fragmentation = fragmentation;
            this.isGlycoPeptide = isGlycoPeptide;
            this.seenWithHCD = seenWithHCD;
            this.seenWithETD = seenWithETD;
            this.NXSmotif = NXSmotif;
            this.NXTmotif = NXTmotif;
            this.isLocalized = isLocalized;
            this.Nlinked = Nlinked;
            this.Olinked = Olinked;
            this.matchedToUniprot = matchedToUniprot;
            this.evidenceType = evidenceType;
            this.evidenceNumber = evidenceNumber;
        }

        public string glycansToString()
        {
            string returnString = "";
            foreach (string glycan in glycans)
            {
                returnString += glycan + ';';
            }

            return returnString.Trim(';');
        }
    }
}
