using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace _20190618_GlycoTools_V2
{
    class GlycoSite
    {
        public HashSet<string> glycans = new HashSet<string>();
        public HashSet<string> peptides = new HashSet<string>();
        public int site;
        public string fasta;
        public string linkage;
        public int psms;
        public bool inUniprot;
        public string evidenceType;
        public bool localized;

        public GlycoSite(PSM psm, int siteIndex)
        {
            var glycanParts = psm.glycansToBeParsed.Split(';');
            glycans.Add(glycanParts[siteIndex]);

            var pepParts = psm.peptidesToBeParsed.Split(',');
            peptides.Add(pepParts[0]);

            this.fasta = psm.proteinName;
            
            this.site = psm.peptideStartPosition + Int32.Parse(psm.glycanPositions.Split(';')[siteIndex]) - 1;

            linkage = psm.modsToBeParsed.Contains("NGlycan") ? "NLinked" : "OLinked";

            psms = 1;

            inUniprot = (psm.evidenceType.Equals("None") || string.IsNullOrEmpty(psm.evidenceType)) ? false : true;

            if (!string.IsNullOrEmpty(psm.evidenceType))
            {
                evidenceType = psm.evidenceType.Split(';')[siteIndex];
            }
            else
            {
                evidenceType = "";
            }
                

            localized = psm.deltaModScore >= 10;
        }

        public override string ToString()
        {
            var uniprotID = fasta.Split('|')[1];

            var returnString = string.Format("'{0}','{1}','{2}','{3}','{4}','{5}','{6}','{7}','{8}','{9}','{10}'",
                                               uniprotID, fasta.Replace("'", "''"), site, localized.ToString(), linkage, psms, peptides.Count(),
                                               glycans.Count(), inUniprot.ToString(), evidenceType, string.Join(";", glycans));

            return returnString;
        }




    }
}
