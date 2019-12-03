using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics.Statistics;
using CSMSL;
using System.IO;

namespace _20190618_GlycoTools_V2
{
    class Glycopeptide
    {
        public List<PSM> psms = new List<PSM>();
        public List<Glycan> glycans = new List<Glycan>();
        public PSM bestPSM;

        public Glycopeptide(PSM psm)
        {
            psms.Add(psm);
            bestPSM = psm;
            foreach(var glycan in psm.glycansToBeParsed.Split(';'))
            {
                glycans.Add(new Glycan(glycan));
            }
        }

        public override string ToString()
        {
            var averageRetentionTime = psms.Select(x => double.Parse(x.scanTime)).ToList().Average();
            var averageRTMinute = averageRetentionTime / 60;
            var stdevRetentionTime = psms.Select(x => double.Parse(x.scanTime)).ToList().StdDev();
            var stdevRTMinute = stdevRetentionTime / 60;
            var averageScore = psms.Select(x => x.score).ToList().Average();
            var stdevScore = psms.Select(x => x.score).ToList().StdDev();
            var glycansString = string.Join(";", glycans.Select(x => x._coreStructure));
            var localized = ((psms.Select(x => x.deltaModScore).ToList()).Max() > 10).ToString();
            var psmCount = psms.Count().ToString();
            var localizedCount = psms.Select(x => x.deltaModScore > 10).Count();
            var glycositeCount = glycans.Count();
            var glycanTypes = string.Join(";", glycans.Select(x => x.glycanType).ToList());
            var linkage = bestPSM.modsToBeParsed.Contains("NGlycan") ? "NLinked" : "OLinked";
            var pep = bestPSM.peptidesToBeParsed.Split(',')[0];

            var NxS = false;
            var NxT = false;

            try
            {
                foreach (var site in bestPSM.glycanPositions.Split(';'))
                {
                    if (pep[Int32.Parse(site) + 1].Equals('S'))
                        NxS = true;

                    if (pep[Int32.Parse(site) + 1].Equals('T'))
                        NxT = true;
                }
            }catch(Exception e)
            {

            }            

            var inUniprot = (bestPSM.evidenceType.Equals("None") || string.IsNullOrEmpty(bestPSM.evidenceType)) ? false : true;

            var chargeStates = string.Join(";", psms.Select(x => x.charge).Distinct().OrderBy(x => x).ToList());
            var obsHCD = psms.Select(x => x.DissociationType).ToList().Contains("HCD");
            var obsETD = psms.Select(x => x.DissociationType).ToList().Contains("ETD");

            var uniprotID = bestPSM.proteinName.Split('|')[1];

            var returnString = string.Format("'{0}','{1}','{2}','{3}','{4}','{5}','{6}','{7}','{8}','{9}'," +
                                                "'{10}','{11}','{12}','{13}','{14}','{15}','{16}','{17}','{18}'," +
                                                "'{19}','{20}','{21}','{22}','{23}','{24}','{25}','{26}','{27}'," +
                                                "'{28}'", 
                                                bestPSM.sequence, pep, bestPSM.modsToBeParsed, uniprotID, bestPSM.proteinName.Replace("'", "''"),
                                                localized, psmCount, localizedCount, glycositeCount, glycansString, glycanTypes,
                                                linkage, NxS, NxT, chargeStates, averageScore, stdevScore, bestPSM.score, bestPSM.DissociationType,
                                                inUniprot, bestPSM.evidenceType, obsHCD, obsETD, averageRTMinute, stdevRTMinute, 
                                                bestPSM.mzObs, bestPSM.scanNumber, Path.GetFileNameWithoutExtension(bestPSM.File), bestPSM.intensity);

            return returnString;
        }



    }
}
