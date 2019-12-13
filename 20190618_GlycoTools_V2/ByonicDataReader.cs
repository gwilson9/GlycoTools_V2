using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Data.SQLite;
using System.IO;
using CSMSL;
using CSMSL.Analysis.Identification;
using CSMSL.IO;
using CSMSL.IO.Thermo;
using CSMSL.Proteomics;
using System.Threading.Tasks;
using System.Diagnostics;
using LumenWorks.Framework.IO.Csv;


namespace _20190618_GlycoTools_V2
{
    class ByonicDataReader
    {
        //Fields
        public List<string[]> files { get; set; }
        public string outputPath { get; set; }
        public double scoreFilter { get; set; }
        public double logProbFilter { get; set; }
        public double deltaModFilter { get; set; }
        public int pepLengthFilter { get; set; }
        public int glycanCountFilter { get; set; }
        public string organism { get; set; }
        public bool performProteinInference { get; set; }
        public string AdditionalPeptidesPath { get; set; }
        public int shortestPeptide;
        public int longestPeptide;
        public string fastaFile;
        public string proteaseString;
        public int MaxMissedCleavage;
        public bool performGlycanLocalization { get; set; }
        private List<FragmentDataReturnArgs> fragData;
        public bool IdentifyForInsourceFragments;

        //Constructor
        public ByonicDataReader(List<string[]> files, string outputPath)
        {
            this.files = files;
            this.outputPath = outputPath;
        }

        public void getNewData()
        {
            using (var writer = createDB())
            {
                var returnDict = new Dictionary<string, Dictionary<string, List<string>>>();
                using (var transaction = writer.BeginTransaction())
                {
                    var filteredGlycoPSMs = new List<PSM>();

                    foreach (var file in files)
                    {
                        var clock = new Stopwatch();
                        clock.Start();

                        List<PSM> allFilePSMs = new List<PSM>();

                        var path = file[0];
                        var condition = file[2];
                        var rep = file[3];
                        var rawFile = file[1];
                        var isControl = file[4];

                        OnUpdateProgress(string.Format("Reading file: {0}", Path.GetFileName(path)));

                        var sqlReader = new SQLiteConnection(@"Data Source=" + path);
                        sqlReader.Open();

                        var reader = readData(path, sqlReader);

                        while (reader.Read())
                        {
                            var psm = readRecord(reader);
                            psm.File = path;
                            psm.Condition = condition;
                            psm.Replicate = rep;
                            psm.isControl = isControl;
                            psm.rawFile = rawFile;

                            // Adds mods to PSM
                            ModifyPSM(psm, sqlReader);

                            if (psm.score >= scoreFilter &&
                                psm.deltaModScore >= deltaModFilter &&
                                psm.sequenceNoMods.Length >= pepLengthFilter &&
                                psm.glycanPositions.Split(';').Count() <= glycanCountFilter &&
                                psm.FDR1D <= 0.01 &&
                                !psm.proteinName.Contains(">Reverse") &&
                                psm.isGlycopeptide)
                            {
                                if (performGlycanLocalization)
                                {
                                    //OnUpdateProgress("Localizing glycan modifications...");

                                    fragData = new List<FragmentDataReturnArgs>();

                                    //Compare each result in fragData by # peptide backbone fragments and remember the original ID
                                    var originalPosition = Int32.Parse(psm.glycanPositions);

                                    var possiblePositions = new List<int>();

                                    //Add original position
                                    possiblePositions.Add(originalPosition);

                                    for (int i = 0; i < psm.sequenceNoMods.Length; i++)
                                    {
                                        //Offset 0 vs 1 indexing
                                        if (i + 1 != originalPosition)
                                        {                                            
                                            if ("ST".Contains(psm.sequenceNoMods[i]))
                                                possiblePositions.Add(i + 1);

                                            if (psm.sequenceNoMods.Length - 1 > (i + 2))
                                            {
                                                if ("N".Contains(psm.sequenceNoMods[i]) && "ST".Contains(psm.sequenceNoMods[i + 2]) && !"P".Contains(psm.sequenceNoMods[i + 1]))
                                                    possiblePositions.Add(i + 1);
                                            }
                                        }    
                                    }

                                    foreach (var position in possiblePositions)
                                    {
                                        var annotator = new glycoFragmentFinder(psm.rawFile, new GlycoPSM(psm, position));
                                        annotator.useType = "SingleSpectrum";
                                        annotator.ReturnData += GlycoFragDataHandler;
                                        annotator.crunch();
                                    }

                                    OnUpdateProgress("Localized Spectra for " + psm.peptidesToBeParsed);
                                    
                                    var originalIndex = -1;
                                    for(int i = 0; i < fragData.Count(); i++)
                                    {
                                        if (fragData[i].glycoPSM.glycanPositions[0] == originalPosition)
                                            originalIndex = i;
                                    }

                                    var bestMatches = fragData[originalIndex].peptideFragmentsMustIncludeGlycan.Count();
                                    var bestIndex = -1;

                                    for(int i = 0; i < fragData.Count(); i++)
                                    {
                                        if(fragData[i].peptideFragmentsMustIncludeGlycan.Count() > bestMatches)
                                        {
                                            bestMatches = fragData[i].peptideFragmentsMustIncludeGlycan.Count();
                                            bestIndex = i;
                                        }
                                    }

                                    // If another location has more fragment ions with intact glycan
                                    if(bestIndex != -1)
                                    {
                                        //Change relevant details of the psm to the better matching peptide
                                        psm.glycanPositions = string.Join(";", fragData[bestIndex].glycoPSM.glycanPositions);                                      
                                        psm.modsToBeParsed = fragData[bestIndex].glycoPSM.formattedVarMods;
                                        psm.sequence = fragData[bestIndex].glycoPSM.peptide.ToString();

                                        var parsedPeptide = psm.peptidesToBeParsed;
                                        foreach(var position in fragData[bestIndex].glycoPSM.glycanPositions)
                                        {
                                            parsedPeptide = parsedPeptide.Replace(originalPosition + "@", position + "@");
                                        }

                                        allFilePSMs.Add(psm);
                                    }
                                    else
                                    {
                                        //Original psm is correctly localized
                                        allFilePSMs.Add(psm);
                                    }

                                }
                                else
                                {
                                    allFilePSMs.Add(psm);
                                }                                
                            }
                        }

                        //Perform all rawFile Tasks
                        allFilePSMs = GetDissociationAndMasterScans(allFilePSMs, rawFile);
                        var psmsQuant = quantifyPeptides(allFilePSMs, rawFile);
                        filteredGlycoPSMs.AddRange(psmsQuant);

                        clock.Stop();
                        OnUpdateProgress("Finished reading file in " + Math.Round((clock.ElapsedMilliseconds / 60000.0), 2).ToString() + " minutes.");

                    }

                    if (performProteinInference)
                    {
                        var protease = Protease.GetProtease(proteaseString);
                        var inferencedGlycoPSMs = InferProteins(filteredGlycoPSMs, AdditionalPeptidesPath, fastaFile, protease);

                        if (inferencedGlycoPSMs != null)
                        {
                            filteredGlycoPSMs = inferencedGlycoPSMs;
                        }
                    }

                    if (IdentifyForInsourceFragments)
                    {

                        SearchForSourceFrags(filteredGlycoPSMs);                        

                    }

                    if (string.IsNullOrEmpty(organism))
                    {
                        OnUpdateProgress("No uniprot taxanomic ID provided. Skipping step.");
                    }
                    else
                    {
                        OnUpdateProgress("Fetching glycosite annotations from Uniprot");
                        var uniprotData = getUniprotData(organism);

                        if (null == uniprotData)
                        {
                            OnUpdateProgress("Failed to get Uniprot Data");
                        }
                        else
                        {
                            var writer2 = new StreamWriter(@"C:\Users\gwilson\Desktop\Temp2\UniprotAPICall.csv");

                            foreach (var pair in uniprotData)
                            {
                                foreach (var pair2 in pair.Value)
                                {
                                    writer2.WriteLine(string.Format("{0},{1},{2}", pair.Key, pair2.Key, pair2.Value));
                                }
                            }

                            writer2.Close();

                            foreach (var psm in filteredGlycoPSMs)
                            {
                                var uniprotID = psm.proteinName.Split('|')[1];
                                if (uniprotData.ContainsKey(uniprotID))
                                {
                                    foreach (var pos in psm.glycanPositions.Split(';'))
                                    {
                                        var siteInProtein = (Int32.Parse(pos) + psm.peptideStartPosition - 1).ToString();

                                        if (uniprotData[uniprotID].ContainsKey(siteInProtein))
                                        {
                                            psm.evidenceType = uniprotData[uniprotID][siteInProtein] + ';';
                                        }
                                        else
                                        {
                                            psm.evidenceType += "None;";
                                        }
                                    }
                                }
                                else
                                {
                                    psm.evidenceType += "None;";
                                }
                                psm.evidenceType = psm.evidenceType.Trim(';');
                            }
                        }
                    }

                    createSQLTables(filteredGlycoPSMs, writer);

                    returnDict = getReturnDict(filteredGlycoPSMs);
                    transaction.Commit();
                }
                writer.Dispose();
                GC.Collect();
                onFinish(returnDict);
            }
        }

        private void SearchForSourceFrags(List<PSM> psms)
        {
            //var insourceFragsFound = 0;
            var searchResults = new Dictionary<PSM, List<PSM>>();
            foreach (var psm in psms)
            {
                OnUpdateProgress(string.Format("Working on {0}", psm.sequence));
                var searcher = new InsourceFragSearcher(psm);
                searcher.fillPossibleParents();
                searcher.getPeakElutions();

                foreach(var match in searcher.matchedParentPeaks)
                {
                    if (searchResults.ContainsKey(psm))
                    {
                        searchResults[psm].Add(match);
                    }
                    else
                    {
                        searchResults.Add(psm, new List<PSM>());
                        searchResults[psm].Add(match);
                    }                    
                }
            }


            /////TEMP FOR VIEWING POSSIBLE INSOURCE FRAGS///////
            if (true)
            {
                SQLiteConnection.CreateFile(outputPath + "\\InsourceFragView.sqlite");
                SQLiteConnection sqlWriter2 = new SQLiteConnection("Data Source=" + outputPath + "\\InsourceFragView.sqlite; Version=3;");
                sqlWriter2.Open();
                using (var transaction1 = sqlWriter2.BeginTransaction())
                {
                    var commandString = "CREATE TABLE IF NOT EXISTS Data (FILE STRING, ID STRING,ID_Glycan STRING, ID_GlycanType STRING,RT STRING,ID_MZ STRING,ID_SCANNUM STRING,ID_RTS STRING, ID_INTENSITIES STRING,PARENT STRING, PARENT_GLYCAN STRING, PARENT_GLYCANTYPE STRING, PARENT_MZ STRING, PARENT_RTS STRING, PARENT_INTENSITIES STRING)";
                    var command = new SQLiteCommand(commandString, sqlWriter2);
                    command.ExecuteNonQuery();


                    foreach (var result in searchResults)
                    {
                        foreach (var parent in result.Value)
                        {
                            var idRTs = string.Join(";", result.Key.peakElution.Select(x => x.RT).ToList());
                            var parentRTs = string.Join(";", parent.peakElution.Select(x => x.RT).ToList());
                            var idInts = string.Join(";", result.Key.peakElution.Select(x => x.Intensity).ToList());
                            var parentInts = string.Join(";", parent.peakElution.Select(x => x.Intensity).ToList());

                            var idGlycanTypes = string.Join(";", result.Key.glycans.Select(x => x.glycanType).ToArray());

                            var parentGlycan = new Glycan(parent.sequence.Split('[')[1].Split(']')[0]);

                            var insertCommandString = string.Format("INSERT into Data ('FILE', `ID`, `ID_Glycan`, 'ID_GlycanType', `RT`, `ID_MZ`, `ID_SCANNUM`,'ID_RTS','ID_INTENSITIES', `PARENT`, 'PARENT_GLYCAN', 'PARENT_GLYCANTYPE', `PARENT_MZ`, 'PARENT_RTS', 'PARENT_INTENSITIES') VALUES ('{0}', '{1}', '{2}', '{3}', '{4}', '{5}', '{6}', '{7}', '{8}', '{9}', '{10}', '{11}', '{12}', '{13}', '{14}')", result.Key.File, result.Key.sequence, result.Key.glycansToBeParsed, idGlycanTypes, (Double.Parse(result.Key.scanTime) / 60).ToString(), result.Key.mzObs, result.Key.scanNumber, idRTs, idInts, parent.sequence, parentGlycan.CoreStructure, parentGlycan.glycanType, parent.mzObs, parentRTs, parentInts);
                            var insertCommand = new SQLiteCommand(insertCommandString, sqlWriter2);
                            var reader = insertCommand.ExecuteReader();
                        }
                    }
                    transaction1.Commit();
                }
                sqlWriter2.Close();
            }

            ////////////////////////////////////////////////////////////////////
        }

        private void GlycoFragDataHandler(object sender, FragmentDataReturnArgs e)
        {
            fragData.Add(e);
        }

        public List<PSM> InferProteins(List<PSM> psms, string peptidesFile, string fastaFile, Protease protease)
        {
            var peptides = GetAllUniquePeptides(psms, peptidesFile);

            var proteins = new List<InferenceProtein>();
            proteins = GetMappedProteinsFromFasta(fastaFile, peptides, protease);

            if (proteins == null || proteins.Count == 0 )
                return null;

            var proteinGroups = GroupProteins(proteins);

            foreach(PSM psm in psms)
            {
                InferencePeptide pep;
                if(peptides.TryGetValue(psm.sequenceNoMods.Replace("I","L"), out pep))
                {
                    pep.ProteinGroups.Sort(InferenceProteinGroup.CompareIncreasing);
                    pep.BestPG = pep.ProteinGroups[0];
                    psm.proteinName = pep.BestPG.Description;
                }
                else
                {
                    OnUpdateProgress("Could not find peptide for psm: " + psm.sequence);
                }
            }

            return psms;
             
        }

        public List<InferenceProteinGroup> GroupProteins(List<InferenceProtein> proteins)
        {
            List<InferenceProteinGroup> proteinGroups = new List<InferenceProteinGroup>();

            // 1) Find Indistinguishable Proteins and group them together into Protein Groups
            // If they are not indistinguishable, then they are still converted to Protein Groups
            // but only contain one protein.
            // A 1 2 3 4
            // B 1 2 3 4
            // C 1   3 4
            // Proteins A and B are indistinguisable (have same set of peptides 1,2,3,4), and thus would become a Protein Group (PG1 [a,b])
            // C is distinguishable and would become a Protein Group (PG2 [c]).
            #region Indistinguishable

            // Loop over each protein
            int p1 = 0;
            while (p1 < proteins.Count)
            {
                // Grab the next protein and its associated peptides from the list of all proteins
                InferenceProtein protein = proteins[p1];
                HashSet<InferencePeptide> peptides = protein.Peptides;

                // Check to see if this protein has enough peptides to be considered indentified
                //if (peptides.Count < MinPeptidesPerGroup)
                //{
                //    // This protein didn't have enough peptides, so remove it from future consideration
                //    proteins.RemoveAt(p1);

                //    // Increase the counter
                //    numberRemovedForNotEnoughPeptides++;

                //    // Go to the next protein on the list
                //    continue;
                //}

                // Start off making the protein into a protein group with its associated peptides
                InferenceProteinGroup pg = new InferenceProteinGroup(protein, peptides);

                // Start looking at the next protein in the list
                int p2 = p1 + 1;

                // Loop over each other protein skipping the one you just made into the PG
                while (p2 < proteins.Count)
                {
                    // Does the next protein contain the same set of peptides as the protein group?
                    if (proteins[p2].Peptides.SetEquals(peptides))
                    {
                        // Yes they are indistinguishable (i.e. proteins A and B from above), so add this protein to the protein group
                        pg.Add(proteins[p2]);

                        // Then remove this protein from the list of all proteins as not to make it into its own PG later
                        proteins.RemoveAt(p2);
                    }
                    else
                    {
                        // Go to next protein in question
                        p2++;
                    }
                }

                // We have gone through every protein possible and thus have completed the grouping of this PG
                proteinGroups.Add(pg);
                p1++;
            }
            
            OnUpdateProgress(string.Format("{0} protein groups are left after combining indistinguishable proteins (having the exact same set of peptides)", proteinGroups.Count));

            #endregion Indistinguishable

            // 2) Find Subsumable Proteins
            // Sort proteins from worst to best to remove the worst scoring groups first (note well, lower p-values mean better scores)
            // Case Example: P-Value, Protein Group, Peptides
            // 0.1  A 1 2
            // 0.05 B 1   3
            // 0.01 C   2 3
            // These are subsumable and we remove the worst scoring protein group (in this case, Protein Group A at p-value of 0.1) first. This would leave:
            // 0.05 B 1   3
            // 0.01 C   2 3
            // Which would mean Protein Group B and C are distinct groups, but share a common peptide (3), peptides 1 and 2 would remain unshared.
            // Protein Group A is removed, as it its peptides can be explained by groups B and C.
            #region Subsumable

            // First, make sure all the peptides know which protein groups they belong too, so we can determine shared peptides
            // and thus get correct p-value for the PGs.
            //MappedPeptidesToProteinGroups(proteinGroups);

            // First update each protein's p-value
            foreach (InferenceProteinGroup proteinGroup in proteinGroups)
            {
                proteinGroup.UpdatePValue();
            }

            // Then sort the groups on decreasing p-values
            proteinGroups.Sort(InferenceProteinGroup.CompareDecreasing);

            p1 = 0;
            while (p1 < proteinGroups.Count)
            {
                // Get the peptides in the protein group
                InferenceProteinGroup proteinGroup = proteinGroups[p1];
                HashSet<InferencePeptide> referencePeptides = proteinGroup.Peptides;

                // Check if all the peptides are shared, if they are then the protein group is subsumable and should be removed
                if (referencePeptides.All(p => p.IsShared))
                {
                    // Since this protein group is being eliminated, remove its reference from all the peptides
                    foreach (InferencePeptide pep in referencePeptides)
                    {
                        pep.ProteinGroups.Remove(proteinGroup);
                    }

                    // Remove the protein group from the master list
                    proteinGroups.RemoveAt(p1);
                }
                else
                {
                    p1++;
                }
            }

            OnUpdateProgress(string.Format("{0} protein groups are left after removing subsumable groups (peptides can be explain by other groups)", proteinGroups.Count));

            #endregion Subsumable

            //// 3) Apply false discovery filtering at the protein level
            //#region FDR filtering

            //proteinGroups.Sort();
            //// Mark each protein group that passes fdr filtering
            //int count = 0;
            //foreach (InferenceProteinGroup proteinGroup in FalseDiscoveryRate<InferenceProteinGroup, double>.Filter(proteinGroups, MaxFdr / 100, true))
            //{
            //    proteinGroup.PassesFDR = true;
            //    count++;
            //}

            //#endregion FDR filtering

            return proteinGroups;
        }

        public List<InferenceProtein> GetMappedProteinsFromFasta(string fastaFile, Dictionary<string, InferencePeptide> peptides, Protease protease)
        {
            if (!File.Exists(fastaFile))
            {
                OnUpdateProgress("Fasta file does not exist, skipping inference.");
                return null;
            }

            OnUpdateProgress(string.Format("Performing {0} digestion on {1}...", protease.Name, Path.GetFileName(fastaFile)));

            int pepsMapped = 0;
            long totalBytes = new FileInfo(fastaFile).Length;

            // A hashset of all proteins that have a peptide that was in the input files
            Dictionary<InferenceProtein, InferenceProtein> proteins = new Dictionary<InferenceProtein, InferenceProtein>(1 << 13);

            int minLength = shortestPeptide - 1;
            int maxLength = longestPeptide + 1;

            using(FastaReader reader = new FastaReader(fastaFile))
            {
                foreach (Fasta fasta in reader.ReadNextFasta())
                {
                    InferenceProtein prot = new InferenceProtein(fasta.Description, fasta.Sequence);

                    foreach (string pepSeq in AminoAcidPolymer.Digest(prot.Sequence, protease, MaxMissedCleavage, minLength, maxLength, true, false))
                    {
                        InferencePeptide pep;
                        if (!peptides.TryGetValue(pepSeq.Replace('I', 'L'), out pep))
                            continue;                        

                        if (!proteins.ContainsKey(prot))
                        {
                            proteins.Add(prot, prot);                            
                        }

                        prot.AddPeptide(pep);

                        if (!pep.IsMapped)
                        {
                            pepsMapped++;
                            pep.IsMapped = true;
                        }
                    }
                }                    
            }

            if(peptides.Count > pepsMapped)
            {
                var upmapped = peptides.Where(x => !x.Value.IsMapped).ToList();

                OnUpdateProgress("Could not map all peptides in Fasta File, skipping inference.\nVerify protein database, protease, and miscleavages.");
                return null;
            }

            OnUpdateProgress(string.Format("Successfully mapped {0} peptides to at least one protein.", pepsMapped));

            return proteins.Values.ToList();
        }

        public Dictionary<string, InferencePeptide> GetAllUniquePeptides(List<PSM> psms, string peptidesFile)
        {
            Dictionary<string, InferencePeptide> peptides = new Dictionary<string, InferencePeptide>(); //Optional Capacity: (1 << 16);

            //Process psms
            foreach (var psm in psms)
            {
                var seq = psm.sequenceNoMods;
                var leuSeq = seq.Replace('I', 'L');              

                var infPSM = new InferencePSM(seq, psm.PEP1D);

                InferencePeptide realPep;
                if(peptides.TryGetValue(leuSeq, out realPep))
                {
                    realPep.PSMs.Add(infPSM);
                }
                else
                {
                    realPep = new InferencePeptide(leuSeq);
                    realPep.PSMs.Add(infPSM);

                    peptides.Add(leuSeq, realPep);

                    if (leuSeq.Length < shortestPeptide)
                        shortestPeptide = leuSeq.Length;

                    if (leuSeq.Length > longestPeptide)
                        longestPeptide = leuSeq.Length;
                }
            }

            var byonicPepCount = peptides.Count();
            OnUpdateProgress(string.Format("{0} unique peptide sequences were found from {1} PSMs loaded from Byonic result files", peptides.Count(), psms.Count()));

            if (File.Exists(peptidesFile))
            {
                using (var reader = new CsvReader(new StreamReader(peptidesFile), true))
                {
                    while (reader.ReadNextRecord())
                    {
                        var seq = reader["Peptide"].ToString().ToUpper();
                        var leuSeq = seq.Replace("I", "L");
                        var pvalue = double.Parse(reader["PEP"]);

                        var infPSM = new InferencePSM(seq, pvalue);
                        InferencePeptide realPep;
                        if(peptides.TryGetValue(leuSeq, out realPep))
                        {
                            realPep.PSMs.Add(infPSM);
                        }
                        else
                        {
                            realPep = new InferencePeptide(seq);
                            realPep.PSMs.Add(infPSM);
                            peptides.Add(leuSeq, realPep);

                            if (leuSeq.Length < shortestPeptide)
                                shortestPeptide = leuSeq.Length;

                            if (leuSeq.Length > longestPeptide)
                                longestPeptide = leuSeq.Length;

                        }
                    }
                }
            }
            else
            {
                OnUpdateProgress("No additional peptides provided for inference.");
                return peptides;
            }

            var extraPepCount = peptides.Count() - byonicPepCount;
            OnUpdateProgress(string.Format("{0} additional peptides loaded from {1} to use for inference.", extraPepCount, Path.GetFileName(AdditionalPeptidesPath)));
            return peptides;
        }        

        public void createSQLTables (List<PSM> psms, SQLiteConnection db)
        {
            //All Data Tables
            foreach(var psm in psms)
            {
                writePSMtoDB(psm, db, "AllGlycoPSMs");
            }
            WriteGlycoProteinTable(psms, db, "AllGlycoproteins");
            WriteGlycoSiteTable(psms, db, "AllGlycosites");
            WriteGlycoPeptidesTable(psms, db, "AllGlycopeptides");
            WriteGlycansTable(psms, db, "AllGlycans");

            var allConditions = psms.Select(x => x.Condition).Distinct().ToList();
            var allReplicates = psms.Select(x => x.Replicate).Distinct().ToList();
            foreach (var i in allConditions)
            {
                //_PSMs
                var query = string.Format("CREATE TABLE GlycoPSMs_Condition_{0} AS SELECT * FROM AllGlycoPSMs WHERE AllGlycoPSMs.Condition='{0}'", i);
                var command = new SQLiteCommand(query, db);
                var reader = command.ExecuteNonQuery();

                //Other tables
                var PSMsubList_Condition = psms.Where(x => x.Condition == i).ToList();
                WriteGlycoProteinTable(PSMsubList_Condition, db, string.Format("Glycoproteins_Condition_{0}", i));
                WriteGlycoSiteTable(PSMsubList_Condition, db, string.Format("Glycosites_Condition_{0}", i));
                WriteGlycoPeptidesTable(PSMsubList_Condition, db, string.Format("Glycopeptides_Condition_{0}", i));
                WriteGlycansTable(PSMsubList_Condition, db, string.Format("Glycans_Condition_{0}", i));

                foreach(var j in allReplicates)
                {
                    //_PSMs
                    var query2 = string.Format("CREATE TABLE GlycoPSMs_Condition_{0}_Replicate_{1} AS SELECT * FROM AllGlycoPSMs WHERE AllGlycoPSMs.Condition='{0}' AND AllGlycoPSMs.Replicate='{1}'", i, j);
                    var command2 = new SQLiteCommand(query2, db);
                    var reader2 = command2.ExecuteNonQuery();

                    //Other table
                    var PSMsubList_Condition_Rep = PSMsubList_Condition.Where(x => x.Replicate == j).ToList();
                    WriteGlycoProteinTable(PSMsubList_Condition_Rep, db, string.Format("Glycoproteins_Condition_{0}_Replicate_{1}", i, j));
                    WriteGlycoSiteTable(PSMsubList_Condition_Rep, db, string.Format("Glycosites_Condition_{0}_Replicate_{1}", i, j));
                    WriteGlycoPeptidesTable(PSMsubList_Condition_Rep, db, string.Format("Glycopeptides_Condition_{0}_Replicate_{1}", i, j));
                    WriteGlycansTable(PSMsubList_Condition_Rep, db, string.Format("Glycans_Condition_{0}_Replicate_{1}", i, j));

                }
            }
        }

        public Dictionary<string, Dictionary<string,string>> getUniprotData(string organism)
        {
            var API = new UniprotAPI(organism);
            return API.GetGlycosites(); 
        }

        public Dictionary<string, Dictionary<string,List<string>>> getReturnDict(List<PSM> psms)
        {
            var returnDict = new Dictionary<string, Dictionary<string, List<string>>>();        

            foreach(var psm in psms)
            {
                var file = Path.GetFileNameWithoutExtension(psm.File);
                var prot = psm.proteinName;

                if (returnDict.ContainsKey(file))
                {
                    if (returnDict[file].ContainsKey(prot))
                    {
                        if(!returnDict[file][prot].Contains(psm.sequence))
                            returnDict[file][prot].Add(psm.sequence);
                    }
                    else
                    {
                        returnDict[file].Add(prot, new List<string>());
                        returnDict[file][prot].Add(psm.sequence);
                    }
                }
                else
                {
                    returnDict.Add(file, new Dictionary<string, List<string>>());
                    returnDict[file].Add(prot, new List<string>());
                    returnDict[file][prot].Add(psm.sequence);
                }                
            }
            return returnDict;
        }

        public List<PSM> quantifyPeptides(List<PSM> psms, string rawPath)
        {
            var processor = new LFQProcessor(rawPath);
            var psmsQuant = processor.crunch(psms);
            return psmsQuant;
        }

        // Generate a new SQLite DB with a PSMs Table Initialized
        public SQLiteConnection createDB()
        {            
            SQLiteConnection.CreateFile(outputPath + "\\MyDatabase.sqlite");            
            SQLiteConnection sqlWriter = new SQLiteConnection("Data Source=" + outputPath + "\\MyDatabase.sqlite; Version=3;");            
            sqlWriter.Open();

            //GlycoPSMs
            var createPSMTableText = CreateGlycoPSMTableString("AllGlycoPSMs");
            var command = new SQLiteCommand(createPSMTableText, sqlWriter);
            command.ExecuteNonQuery();

            return sqlWriter;
        }

        // Read in Byonic Results
        public SQLiteDataReader readData(string path, SQLiteConnection sqlReader)
        {     
            var query = "SELECT ProteinsFoundPQMs.DebugText, PQMs.PeptideParse, PQMs.ProteinStartPosition, " +
                                    "PQMs.PosteriorErrorProbability2, PQMs.PosteriorErrorProbability1," +
                                    "PQMs.FalseDiscoveryRate2, PQMs.FalseDiscoveryRate1, PQMs.FalseDiscoveryRateUnique2, " +
                                    "PQMs.FalseDiscoveryRateUnique1, PQMs.PosteriorErrorProbability2_sum," +
                                    "PQMs.PosteriorErrorProbability1_sum, PQMs.CalcMz, ProteinsFoundPQMs.PQMsId, " +
                                    "PQMsQueriesSummary.Intensity, Queries.ObservedMz, Queries.ScanTimeList, " +
                                    "ProteinsFoundPQMs.ProteinRank, PQMs.Score, PQMs.DeltaScoreSeq, PQMs.DeltaScoreSeqMod," +
                                    "PQMs.Charge, PQMs.ObservedMh, PQMs.CalcMH, ProteinsFoundPQMs.Cleavage, " +
                                    "Proteins.ProteinName,Proteins.Id, Queries.ScanNumberList " +
                            "FROM PQMs JOIN ProteinsFoundPQMs ON PQMS.Id = ProteinsFoundPQMs.PQMsId " +
                            "JOIN Queries ON PQMs.QueriesId = Queries.Id " +
                            "JOIN Proteins ON Proteins.Id = PQMs.ProteinsId " +
                            "JOIN PQMsQueriesSummary ON PQMsQueriesSummary.QueriesID = PQMs.QueriesID";


            var command = new SQLiteCommand(query, sqlReader);
            var reader = command.ExecuteReader();

            return reader;
        }

        // Helper to read in Byonic Results
        public PSM readRecord(SQLiteDataReader reader)
        {
            string protRank = reader["ProteinRank"].ToString();
            string pqmsID = reader["pqmsID"].ToString();
            string sequence = reader["DebugText"].ToString().Substring(2,(reader["DebugText"].ToString().Length - 4));   //Trins the proximal amino acids off of the original text
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
                foreach(var glycan in psm.glycansToBeParsed.Split(';'))
                {
                    psm.glycans.Add(new Glycan(glycan));
                }
            }
        }

        // Write PSM to SQLiteDB
        public void writePSMtoDB(PSM psm, SQLiteConnection writer, string table)
        {
            var insertPSM = string.Format("INSERT into {0}(`Sequence`, `PeptideParseFriendly`, `Peptide`, `Position`, `ModsVar`, `Glycans`, " +
                                              "`GlycanTypes`, `PEP2D`, `PEP1D`, `logProb`, `Score`, `DeltaScore`, `DeltaModScore`, `Charge`, " +
                                              "`ObsMZ`, `CalcMZ`, `ppmerr`, `ObsMH`, `CalcMH`, `Cleavage`, `GlycansPos`, " +
                                              "`ProteinFasta`, `ScanTime`, `ScanNum`, `ModsFixed`, `FDR2D`, `FDR1D`, `FDRuniq2D`, " +
                                              "`FDRuniq1D`, `qValue2D`, `qValue1D`, `isGlycoPeptide`, `modsPassedCheck`, `positionPassedCheck`, " +
                                              "`DissociationType`, `MasterScan`, `LFQIntensity`, `File`, `Condition`, `Replicate`, `isControl`, " +
                                              "'RetentionTimes', 'Intensities', 'EvidenceType', 'MS1ScanNumOfMaxIntensity', 'RTofMaxIntensity', 'RTofLocalMaxima') " +
                                          "VALUES ({1})", table, psm.ToString());

            var command = new SQLiteCommand(insertPSM, writer);
            var reader = command.ExecuteReader();
        }

        public event EventHandler<ProgressEventArgs> UpdateProgress;

        protected virtual void OnUpdateProgress(string progress)
        {
            var handler = UpdateProgress;
            if(handler != null)
            {
                handler(this, new ProgressEventArgs(progress));
            }
        }

        public EventHandler<DataReturnArgs> finish;

        protected virtual void onFinish(Dictionary<string, Dictionary<string, List<string>>> data)
        {
            var handler = finish;

            if (handler != null)
            {
                handler(this, new DataReturnArgs(data));
            }
        }

        public List<PSM> GetDissociationAndMasterScans(List<PSM> psms, string rawPath)
        {
            var raw = new ThermoRawFile(rawPath);
            raw.Open();

            foreach(var psm in psms)
            {
                psm.DissociationType = raw.GetDissociationType(psm.scanNumber).ToString();

                //Triggered scans will have a parent scan number of another MS2...must find parent MS1
                var parentScan = raw.GetParentSpectrumNumber(psm.scanNumber);
                while (true)
                {
                    if (raw.GetMsnOrder(parentScan) == 2)
                    {
                        parentScan = raw.GetParentSpectrumNumber(parentScan);
                    }
                    else
                    {
                        break;
                    }
                }
                psm.MasterScan = parentScan;
            }

            raw.Dispose();       

            return psms;
        }       

        public void WriteGlycoProteinTable(List<PSM> psms, SQLiteConnection writer, string table)
        {
            var glycoProts = new Dictionary<string, GlycoProtein>();            

            foreach(var psm in psms)
            {
                var fasta = psm.proteinName;

                if (!glycoProts.ContainsKey(fasta))
                {
                    glycoProts.Add(fasta, new GlycoProtein(fasta));

                    glycoProts[fasta].psms.Add(psm);
                    
                    foreach(var site in psm.glycanPositions.Split(';'))
                    {
                        glycoProts[fasta].uniqueSites.Add(psm.peptideStartPosition + Int32.Parse(site) - 1 );
                    }

                    foreach(var glycan in psm.glycansToBeParsed.Split(';'))
                    {
                        glycoProts[fasta].uniqueGlycans.Add(glycan);
                    }

                    if(psm.deltaModScore > 10)
                    {                       

                        foreach (var site in psm.glycanPositions.Split(';'))
                        {
                            glycoProts[fasta].uniqueLocalizedSites.Add(psm.peptideStartPosition + Int32.Parse(site) - 1);
                        }

                        foreach (var glycan in psm.glycansToBeParsed.Split(';'))
                        {
                            glycoProts[fasta].uniqueLocalizedGlycans.Add(glycan);
                        }
                    }

                }
                else
                {
                    glycoProts[fasta].psms.Add(psm);

                    foreach (var site in psm.glycanPositions.Split(';'))
                    {
                        glycoProts[fasta].uniqueSites.Add(psm.peptideStartPosition + Int32.Parse(site) - 1);
                    }

                    foreach (var glycan in psm.glycansToBeParsed.Split(';'))
                    {
                        glycoProts[fasta].uniqueGlycans.Add(glycan);
                    }

                    if (psm.deltaModScore > 10)
                    {

                        foreach (var site in psm.glycanPositions.Split(';'))
                        {
                            glycoProts[fasta].uniqueLocalizedSites.Add(psm.peptideStartPosition + Int32.Parse(site) - 1);
                        }

                        foreach (var glycan in psm.glycansToBeParsed.Split(';'))
                        {
                            glycoProts[fasta].uniqueLocalizedGlycans.Add(glycan);
                        }
                    }
                }                
            }

            var CreateProteinTableText = CreateGlycoproteinsTableString(table);
            var command2 = new SQLiteCommand(CreateProteinTableText, writer);
            command2.ExecuteNonQuery();

            foreach (var pair in glycoProts)
            {
                var insertProt = string.Format("INSERT into {0}(`UniprotID`, `Description`, `InUniprotAsGlycoprotein`, `LocalizedPSMs`, `LocalizedGlycoSites`, " +
                                                                        "`LocalizedGlyanCount`,  `LocalizedGlycans`, `PSMs`, `GlycoSites`, `GlycanCount`, `Glycans`) " +
                                              "VALUES ({1})", table, pair.Value.ToString());

                var command = new SQLiteCommand(insertProt, writer);
                var reader = command.ExecuteReader();
            }
        }

        public void WriteGlycoSiteTable(List<PSM> psms, SQLiteConnection writer, string table)
        {
            var allSites = new Dictionary<string, GlycoSite>();

            foreach(var psm in psms)
            {
                var positions = psm.glycanPositions.Split(';');

                for(int i = 0; i < positions.Count(); i++)
                {
                    var key = psm.proteinName + "_" + (psm.peptideStartPosition + Int32.Parse(positions[i]) - 1).ToString();

                    if (allSites.ContainsKey(key))
                    {
                        allSites[key].glycans.Add(psm.glycansToBeParsed.Split(';')[i]);
                        allSites[key].peptides.Add(psm.peptidesToBeParsed.Split(',')[0]);
                        allSites[key].psms++;
                        if (!allSites[key].localized)
                            allSites[key].localized = psm.deltaModScore >= 10;
                    }
                    else
                    {
                        allSites.Add(key, new GlycoSite(psm, i));
                    }
                }
            }

            var CreateGlycositesTableText = CreateGlycoSitesTableString(table);
            var command3 = new SQLiteCommand(CreateGlycositesTableText, writer);
            command3.ExecuteNonQuery();

            foreach (var pair in allSites)
            {
                var insertSite = string.Format("INSERT into {0}(`UniprotID`, `Description`, `Site`, `Localized`, `GlycoType`, " +
                                                                        "`PSMs`,  `UniqueSequences`, `GlycanCount`, `SeenInUniprot`, `EvidenceType`, `Glycans`) " +
                                              "VALUES ({1})", table, pair.Value.ToString());

                var command = new SQLiteCommand(insertSite, writer);
                var reader = command.ExecuteReader();
            }

        }

        public void WriteGlycoPeptidesTable(List<PSM> psms, SQLiteConnection writer, string table)
        {
            var allPeps = new Dictionary<string, Glycopeptide>();

            foreach(var psm in psms)
            {
                if (allPeps.ContainsKey(psm.sequence + "_" + psm.proteinName))
                {
                    allPeps[psm.sequence + "_" + psm.proteinName].psms.Add(psm);
                    if (allPeps[psm.sequence + "_" + psm.proteinName].bestPSM.score < psm.score)
                        allPeps[psm.sequence + "_" + psm.proteinName].bestPSM = psm;
                }
                else
                {
                    allPeps.Add(psm.sequence + "_" + psm.proteinName, new Glycopeptide(psm));                    
                }
            }

            var CreateGlycopeptidesTableText = CreateGlycoPeptidesTableString(table);
            var command4 = new SQLiteCommand(CreateGlycopeptidesTableText, writer);
            command4.ExecuteNonQuery();

            foreach (var pair in allPeps)
            {
                var insertPep = string.Format("INSERT into {0}(`PeptideWithMods`, `Peptide`, `Mods`, `uniprotID`, `ProteinName`, `Localized`, " +
                                                                         "`PSMCount`, `LocalizedPSMCount`, `GlycoSiteCount`, `Glycans`, `GlycanTypes`,  " +
                                                                         "`Linkage`, `NxS`, `NxT`, `ChargeStates`, `AverageScore`, `ScoreVariance`, " +
                                                                         "`BestScore`, `FragmentationBestPSM`, `inUniprot`, `EvidenceType`, `ObservedHCD`, " +
                                                                         "`ObservedETD`, `AverageRetentionTime`, `RetentionTimeVariance`, `PrecursorMZ_BestScore`, " +
                                                                         "`BestSpectrum`, `BestSpectrumFile`, `LFQIntensity`) " +
                                              "VALUES ({1})", table, pair.Value.ToString());
                var command = new SQLiteCommand(insertPep, writer);
                var reader = command.ExecuteReader();
            }


        }

        public void WriteGlycansTable(List<PSM> psms, SQLiteConnection writer, string table)
        {
            var allGlycans = new Dictionary<string, Glycan>();

            foreach(var psm in psms)
            {
                var glycans = psm.glycansToBeParsed.Split(';');
                var positions = psm.glycanPositions.Split(';');

                for(int i = 0; i < glycans.Count(); i++)
                {
                    var glycosite = psm.proteinName + '_' + (psm.peptideStartPosition + Int32.Parse(positions[i]) - 1);
                    var uniprotId = psm.proteinName.Split('|')[1];

                    if (allGlycans.ContainsKey(glycans[i]))
                    {
                        allGlycans[glycans[i]].glycosites.Add(glycosite);
                        allGlycans[glycans[i]].uniprotIDs.Add(uniprotId);
                    }
                    else
                    {
                        var glycan = new Glycan(glycans[i]);
                        glycan.uniprotIDs.Add(uniprotId);
                        glycan.glycosites.Add(glycosite);

                        allGlycans.Add(glycans[i], glycan);
                    }
                }
            }

            var createGlycanTableText = CreateGlycansTableString(table);
            var command5 = new SQLiteCommand(createGlycanTableText, writer);
            command5.ExecuteNonQuery();

            foreach (var pair in allGlycans)
            {
                var insertGlycan = string.Format("INSERT into {0}(`Glycan`, `GlycanType`, `Mass`, `UniqueSites`, `UniqueProteins`, `UniprotIDs`) " +
                                              "VALUES ({1})", table, pair.Value.ToString());
                var command = new SQLiteCommand(insertGlycan, writer);
                var reader = command.ExecuteReader();
            }


        }

        public string CreateGlycoPSMTableString(string name)
        {
            return String.Format("CREATE TABLE IF NOT EXISTS {0} (Sequence STRING, " +
                                                                    "PeptideParseFriendly STRING, " +
                                                                    "Peptide STRING, " +
                                                                    "Position INT, " +
                                                                    "ModsVar STRING, " +
                                                                    "Glycans STRING, " +
                                                                    "GlycanTypes STRING, " +
                                                                    "PEP2D REAL, " +
                                                                    "PEP1D REAL, " +
                                                                    "logProb REAL, " +
                                                                    "Score REAL, " +
                                                                    "DeltaScore REAL, " +
                                                                    "DeltaModScore REAL, " +
                                                                    "Charge INT, " +
                                                                    "ObsMZ REAL, " +
                                                                    "CalcMZ REAL, " +
                                                                    "ppmerr REAL, " +
                                                                    "ObsMH REAL, " +
                                                                    "CalcMH REAL, " +
                                                                    "Cleavage STRING, " +
                                                                    "GlycansPos STRING, " +
                                                                    "ProteinFasta STRING, " +
                                                                    "ScanTime REAL, " +
                                                                    "ScanNum INT, " +
                                                                    "ModsFixed STRING, " +
                                                                    "FDR2D REAL, " +
                                                                    "FDR1D REAL, " +
                                                                    "FDRuniq2D REAL, " +
                                                                    "FDRuniq1D REAL, " +
                                                                    "qValue2D REAL, " +
                                                                    "qValue1D REAL, " +
                                                                    "isGlycoPeptide INT, " +
                                                                    "modsPassedCheck INT, " +
                                                                    "positionPassedCheck INT, " +
                                                                    "DissociationType STRING, " +
                                                                    "MasterScan INT, " +
                                                                    "LFQIntensity REAL, " +
                                                                    "File STRING, " +
                                                                    "Condition STRING, " +
                                                                    "Replicate STRING, " +
                                                                    "isControl INT, " +
                                                                    "RetentionTimes STRING, " +
                                                                    "Intensities STRING, " +
                                                                    "EvidenceType STRING, " +
                                                                    "MS1ScanNumOfMaxIntensity INT, " +
                                                                    "RTofMaxIntensity REAL, " +
                                                                    "RTofLocalMaxima STRING)", name);
        }

        public string CreateGlycoPeptidesTableString(string name)
        {
            return String.Format("CREATE TABLE IF NOT EXISTS {0} (PeptideWithMods STRING, " +
                                                                                         "Peptide STRING, " +
                                                                                         "Mods STRING, " +
                                                                                         "uniprotID STRING, " +
                                                                                         "ProteinName STRING, " +
                                                                                         "Localized STRING, " +
                                                                                         "PSMCount INT, " +
                                                                                         "LocalizedPSMCount INT, " +
                                                                                         "GlycoSiteCount INT, " +
                                                                                         "Glycans STRING, " +
                                                                                         "GlycanTypes STRING, " +
                                                                                         "Linkage STRING, " +
                                                                                         "NxS STRING, " +
                                                                                         "NxT STRING, " +
                                                                                         "ChargeStates STRING, " +
                                                                                         "AverageScore REAL, " +
                                                                                         "ScoreVariance REAL, " +
                                                                                         "BestScore REAL, " +
                                                                                         "FragmentationBestPSM STRING, " +
                                                                                         "inUniprot STRING, " +
                                                                                         "EvidenceType STRING, " +
                                                                                         "ObservedHCD STRING, " +
                                                                                         "ObservedETD STRING, " +
                                                                                         "AverageRetentionTime REAL, " +
                                                                                         "RetentionTimeVariance REAL, " +
                                                                                         "PrecursorMZ_BestScore REAL, " +
                                                                                         "BestSpectrum INT, " +
                                                                                         "BestSpectrumFile STRING, " +
                                                                                         "LFQIntensity REAL)", name);
        }

        public string CreateGlycoSitesTableString(string name)
        {
            return string.Format("CREATE TABLE IF NOT EXISTS {0} (UniprotID STRING, " +
                                                                        "Description STRING, " +
                                                                        "Site INT, " +
                                                                        "Localized STRING, " +
                                                                        "GlycoType STRING, " +
                                                                        "PSMs INT, " +
                                                                        "UniqueSequences INT, " +
                                                                        "GlycanCount INT, " +
                                                                        "SeenInUniprot STRING, " +
                                                                        "EvidenceType STRING, " +
                                                                        "Glycans STRING)", name);
        }

        public string CreateGlycoproteinsTableString(string name)
        {
            return string.Format("CREATE TABLE IF NOT EXISTS {0} (UniprotID STRING, " +
                                                                    "Description STRING, " +
                                                                    "InUniprotAsGlycoprotein STRING, " +
                                                                    "LocalizedPSMs INT, " +
                                                                    "LocalizedGlycoSites STRING, " +
                                                                    "LocalizedGlyanCount STRING, " +
                                                                    "LocalizedGlycans STRING, " +
                                                                    "PSMs STRING, " +
                                                                    "GlycoSites STRING, " +
                                                                    "GlycanCount STRING, " +
                                                                    "Glycans STRING)", name);
        }

        public string CreateGlycansTableString(string name)
        {
            return string.Format("CREATE TABLE IF NOT EXISTS {0} (Glycan STRING, " +
                                                                        "GlycanType STRING, " +
                                                                        "Mass REAL, " +
                                                                        "UniqueSites INT, " +
                                                                        "UniqueProteins INT, " +
                                                                        "UniprotIDs STRING)", name);
        }

    }
}
