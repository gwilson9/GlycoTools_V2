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
    class InsourceFragSearcher
    {

        public PSM idPSM;
        public List<PSM> possibleParentPeaks;
        public List<PSM> matchedParentPeaks;

        public InsourceFragSearcher(PSM psm)
        {
            this.idPSM = psm;
            possibleParentPeaks = new List<PSM>();
            matchedParentPeaks = new List<PSM>();
        }

        public void fillPossibleParents()
        {
            var psmGlycans = idPSM.glycans;
            var glycanMasses = idPSM.glycans.Select(x => x.mass).ToList();

            // <ID'd mass, List<All larger glycan masses>
            var searchSpace = new Dictionary<Glycan, List<Glycan>>();

            var glycanStrings = new GlycanConstants();

            foreach (var psmGlycan in psmGlycans)
            {
                var monomers = psmGlycan._allSugars.Select(x => x.Name).ToList().Distinct();

                foreach (var glycanString in glycanStrings.glycans)
                {
                    var glycan = new Glycan(glycanString);

                    //Check to see if glycan to be considered would be a logical parent glycan of the ID's glycan
                    var passes = true;
                    foreach (var monomer in monomers)
                    {
                        var idMonomerCount = psmGlycan._allSugars.Where(x => x.Name.Equals(monomer)).ToList().Count();
                        var glycanMonomerCount = glycan._allSugars.Where(x => x.Name.Equals(monomer)).ToList().Count();
                        if (glycanMonomerCount < idMonomerCount || glycan.mass <= psmGlycan.mass)
                        {
                            passes = false;
                            break;
                        }
                    }

                    if (!passes)
                        continue;


                    if (searchSpace.ContainsKey(psmGlycan))
                    {
                        searchSpace[psmGlycan].Add(glycan);
                    }
                    else
                    {
                        searchSpace.Add(psmGlycan, new List<Glycan>());
                        searchSpace[psmGlycan].Add(glycan);
                    }
                }
            }



            var pepsToSearch = new List<Peptide>();

            if (idPSM.glycanPositions.Count() == 1 && searchSpace.Count() > 0)
            {
                foreach (var glycan in searchSpace[idPSM.glycans[0]])
                {
                    var pep = new Peptide(idPSM.sequenceNoMods);

                    if (!String.IsNullOrEmpty(idPSM.modsToBeParsed))
                    {
                        foreach (var mod in idPSM.modsToBeParsed.Split(';'))
                        {
                            if (!String.IsNullOrEmpty(mod) && !mod.Contains("Glycan"))
                            {
                                var modMass = Double.Parse(mod.Split(' ')[2].Trim(')'));
                                var modName = mod.Split('(')[1].Split(' ')[0];
                                var newMod = new Modification(modMass, modName);
                                pep.AddModification(newMod);
                            }
                        }
                        foreach (var mod in idPSM.modsFixed.Split(';'))
                        {
                            if (!String.IsNullOrEmpty(mod) && !mod.Contains("Glycan"))
                            {
                                var massString = mod.Split(' ')[2].Trim(')');
                                var modMass = Double.Parse(massString);
                                var modName = mod.Split('(')[1].Split(' ')[0];
                                var newMod = new Modification(modMass, modName);
                                pep.AddModification(newMod);
                            }
                        }

                        var glyMod = new Modification(glycan.mass, glycan._coreStructure);
                        pep.AddModification(glyMod, Int32.Parse(idPSM.glycanPositions));

                        pepsToSearch.Add(pep);
                    }
                }
            }


            //var raw = new ThermoRawFile(psm.rawFile);
            //raw.Open();

            pepsToSearch = pepsToSearch.OrderByDescending(x => x.MonoisotopicMass).ToList();

            foreach (var pep in pepsToSearch)
            {
                //OnUpdateProgress("Still working...");
                var mz = (pep.MonoisotopicMass + (idPSM.charge * Constants.Hydrogen)) / idPSM.charge;

                var searchPSM = new PSM();
                searchPSM.scanNumber = idPSM.scanNumber;
                searchPSM.peptide = idPSM.peptide;
                searchPSM.sequence = pep.SequenceWithModifications;
                searchPSM.charge = idPSM.charge;
                searchPSM.mzCalc = idPSM.mzCalc;
                searchPSM.mzObs = mz;

                possibleParentPeaks.Add(searchPSM);
            }


                /**
                var lastSpectrum = raw.LastSpectrumNumber;
                
                var centerSpectrum = raw.GetSpectrum(psm.scanNumberofMaxElutionIntensity);

                var left = psm.scanNumberofMaxElutionIntensity;
                while (true)
                {
                    if (raw.GetMsnOrder(--left) == 1)
                        break;
                }

                var right = psm.scanNumberofMaxElutionIntensity;
                while (true)
                {
                    if (raw.GetMsnOrder(++right) == 1)
                        break;
                }

                var leftSpectrum = raw.GetSpectrum(left);
                var rightSpectrum = raw.GetSpectrum(right);


                
                var centerIntensity = 0.0;
                var leftIntensity = 0.0;
                var rightIntensity = 0.0;

                var outpeaks = new List<ThermoMzPeak>();
                if(centerSpectrum.TryGetPeaks(DoubleRange.FromPPM(mz, 20), out outpeaks))
                {
                    if(outpeaks.Count() > 0)
                    {
                        centerIntensity = outpeaks.OrderByDescending(x => x.Intensity).ToList()[0].Intensity;
                    }
                    else
                    {
                        continue;
                    }
                }                                

                var outpeaks2 = new List<ThermoMzPeak>();
                if (leftSpectrum.TryGetPeaks(DoubleRange.FromPPM(mz, 20), out outpeaks2))
                {
                    if(outpeaks2.Count() > 0)
                    {
                        leftIntensity = outpeaks2.OrderByDescending(x => x.Intensity).ToList()[0].Intensity;
                    }
                    else
                    {
                        continue;
                    }
                }                                

                var outpeaks3 = new List<ThermoMzPeak>();
                if (rightSpectrum.TryGetPeaks(DoubleRange.FromPPM(mz, 20), out outpeaks3))
                {
                    if(outpeaks3.Count() > 0)
                    {
                        rightIntensity = outpeaks3.OrderByDescending(x => x.Intensity).ToList()[0].Intensity;
                    }
                    else
                    {
                        continue;
                    }
                }
               


                var maxIntensity = psm.peakElution.Aggregate((i1, i2) => i1.Intensity > i2.Intensity ? i1 : i2).Intensity;
                if (centerIntensity > maxIntensity && centerIntensity > leftIntensity && centerIntensity > rightIntensity)
                {
                    if (psm.sequence.Equals("N[+1540.52851]ATLAEQAK") && pep.SequenceWithModifications.Equals("N[HexNAc(2)Hex(8)]ATLAEQAK"))
                    {                                        
                        var peaks = psm.peakElution.OrderByDescending(x=>x.Intensity).ToList()[0];
                    }

                    insourceFragsFound++;

                    if (searchResults.ContainsKey(psm))
                    {
                        searchResults[psm].Add(pep);
                    }
                    else
                    {
                        searchResults.Add(psm, new List<Peptide>());
                        searchResults[psm].Add(pep);
                    }
                    
                    //break;
                }
**/
            
            //raw.Dispose(); 
        }

        public void getPeakElutions()
        {
            var lfqProcessor = new LFQProcessor(idPSM.rawFile);
            var pepLFQs = lfqProcessor.crunch(possibleParentPeaks);

            foreach(var pepLFQ in pepLFQs)
            {
                var matchMax = pepLFQ.peakElution.OrderByDescending(x => x.Intensity).ToList()[0].Intensity;
                var idMax = idPSM.peakElution.OrderByDescending(x => x.Intensity).ToList()[0].Intensity;
                if (pepLFQ.scanNumberofMaxElutionIntensity == idPSM.scanNumberofMaxElutionIntensity && idMax < matchMax)
                {
                    matchedParentPeaks.Add(pepLFQ);
                }
            }            
        }
    }
}
