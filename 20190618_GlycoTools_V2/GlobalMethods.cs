using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Data.SQLite;
using CSMSL;
using CSMSL.Proteomics;


namespace _20190618_GlycoTools_V2
{
    class GlobalMethods
    {
        public PSM PSMFromGlycoPSMsSQLTable(SQLiteDataReader reader)
        {
            string protRank = reader["ProteinRank"].ToString();
            string pqmsID = reader["pqmsID"].ToString();
            string sequence = reader["DebugText"].ToString().Substring(2, (reader["DebugText"].ToString().Length - 4));   //Trins the proximal amino acids off of the original text
            string peptidesToBeParsed = reader["PeptideParse"].ToString();
            int peptideStartPosition = int.Parse(reader["ProteinStartPosition"].ToString());
            double PEP2D = double.Parse(reader["PosteriorErrorProbability2"].ToString());
            double PEP1D = double.Parse(reader["PosteriorErrorProbability1"].ToString());
            double logProb = Math.Abs(Math.Log10(PEP1D));                                 //Ignored for now
            double score = double.Parse(reader["Score"].ToString());
            double deltaScore = double.Parse(reader["DeltaScoreSeq"].ToString());
            double deltaModScore = double.Parse(reader["DeltaScoreSeqMod"].ToString());
            int charge = int.Parse(reader["Charge"].ToString());
            double mzObs = double.Parse(reader["ObservedMz"].ToString());
            double mzCalc = double.Parse(reader["CalcMz"].ToString());
            double obsMH = double.Parse(reader["ObservedMh"].ToString());
            double calcMH = double.Parse(reader["CalcMH"].ToString());
            string cleavage = reader["Cleavage"].ToString();
            string proteinName = reader["ProteinName"].ToString();
            int protID = int.Parse(reader["Id"].ToString());
            string scanTime = reader["ScanTimeList"].ToString();

            // For scan number
            var scanNumberString = reader["ScanNumberList"].ToString();
            var stringParts = scanNumberString.Split('=');
            int scanNumber = int.Parse(stringParts[stringParts.Length - 1]);

            double FDR2D = double.Parse(reader["FalseDiscoveryRate2"].ToString());
            double FDR1D = double.Parse(reader["FalseDiscoveryRate1"].ToString());
            double FDR2Dunique = double.Parse(reader["FalseDiscoveryRateUnique2"].ToString());
            double FDR1Dunique = double.Parse(reader["FalseDiscoveryRateUnique1"].ToString());
            double qvalue2D = double.Parse(reader["PosteriorErrorProbability2_sum"].ToString());
            double qvalue1D = double.Parse(reader["PosteriorErrorProbability1_sum"].ToString());
            double intensity = double.Parse(reader["intensity"].ToString());

            PSM newPSM = new PSM(pqmsID, sequence, peptidesToBeParsed, peptideStartPosition, PEP2D, PEP1D, score,
                deltaScore, deltaModScore, charge, mzObs, mzCalc, obsMH, calcMH, cleavage, proteinName, protID,
                scanNumber, FDR2D, FDR1D, FDR2Dunique, FDR1Dunique, qvalue2D, qvalue1D, intensity, scanTime, protRank);

            return newPSM;
        }

        // Gather additional info from Byonic Results for PSM
        public void ModifyPSM(PSM psm, SQLiteConnection sqlReader)
        {
            var query = "SELECT Modifications.Classification, Modifications.Composition, PQMsPeptideToModifications.ModificationsPeptidePosition," +
                        "Modifications.AllowedSites, Modifications.MonoMassShiftTotal " +
                        "FROM PQMsPeptideToModifications " +
                        "JOIN PQMs ON PQMsPeptideToModifications.PQMsId = PQMs.Id " +
                        "JOIN Modifications ON Modifications.Id = PQMsPeptideToModifications.ModificationsId " +
                        "WHERE PQMsPeptideToModifications.PQMsId=" + psm.PQMsID;
            var command = new SQLiteCommand(query, sqlReader);
            var reader = command.ExecuteReader();

            string varMods = "";
            string fixedMods = "";
            string glycans = "";
            string glycanPos = "";

            while (reader.Read())
            {
                if (reader["Classification"].ToString().Equals("ngly") || reader["Classification"].ToString().Equals("ogly"))
                {
                    glycans += reader["Composition"] + ";";
                    glycanPos += reader["ModificationsPeptidePosition"] + ";";


                    if (reader["Classification"].ToString().Equals("ngly"))
                    {
                        varMods += "N" + reader["ModificationsPeptidePosition"] + "(" + reader["AllowedSites"] + " / "
                                    + reader["MonoMassShiftTotal"] + ");";
                    }
                    else
                    {
                        var position = Int32.Parse(reader["ModificationsPeptidePosition"].ToString());
                        varMods += psm.peptidesToBeParsed[position - 1] + reader["ModificationsPeptidePosition"].ToString() +
                                    "(" + reader["AllowedSites"] + " / " + reader["MonoMassShiftTotal"] + ");";
                    }

                }
                else
                {
                    if (reader["AllowedSites"].ToString().Equals("M"))
                    {
                        varMods += "M" + reader["ModificationsPeptidePosition"] + "(Oxidation / 15.9949);";
                    }

                    if (reader["AllowedSites"].ToString().Equals("NTerm E"))
                    {
                        varMods += "E" + reader["ModificationsPeptidePosition"] + "(Glu->pyro-Glu / -18.0106);";
                    }
                    if (reader["AllowedSites"].ToString().Equals("NTerm Q"))
                    {
                        varMods += "Q" + reader["ModificationsPeptidePosition"] + "(Gln->pyro-Glu / -17.026549);";
                    }
                    if (reader["AllowedSites"].ToString().Equals("N"))
                    {
                        varMods += "N" + reader["ModificationsPeptidePosition"] + "(Deamidated / 0.9840);";
                    }
                    if (reader["AllowedSites"].ToString().Equals("C"))
                    {
                        fixedMods += "C" + reader["ModificationsPeptidePosition"] + "(Carbamidomethyl / 57.021464);";
                    }
                }
            }
            psm.modsFixed = fixedMods.TrimEnd(';');
            psm.modsToBeParsed = varMods.TrimEnd(';');
            psm.glycansToBeParsed = glycans.TrimEnd(';');
            psm.glycanPositions = glycanPos.TrimEnd(';');
            psm.checkMods();
            psm.isGlycopeptide = psm.modsToBeParsed.Contains("Glycan") ? true : false;
            psm.peptide = new Peptide(psm.peptidesToBeParsed.Split(',')[0]);

            if (psm.isGlycopeptide)
            {
                foreach (var glycan in psm.glycansToBeParsed.Split(';'))
                {
                    psm.glycans.Add(new Glycan(glycan));
                }
            }
        }




    }
}
