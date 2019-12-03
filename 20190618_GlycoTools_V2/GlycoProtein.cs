using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace _20190618_GlycoTools_V2
{
    class GlycoProtein
    {
        public string fasta;
        public string uniprotID;
        public string inUniprot;
        public List<PSM> psms = new List<PSM>();
        public HashSet<string> uniqueLocalizedGlycans = new HashSet<string>();
        public HashSet<int> uniqueLocalizedSites = new HashSet<int>();
        public HashSet<string> uniqueGlycans = new HashSet<string>();
        public HashSet<int> uniqueSites = new HashSet<int>();

        public GlycoProtein(string fasta)
        {
            this.fasta = fasta;
            this.uniprotID = fasta.Split('|')[1];
        }

        public override string ToString()
        {

            var PSMCount = psms.Count();
            var localizedPSMCount = psms.Where(x => x.deltaModScore > 10).ToList().Count();
            var localizedGlycanCount = uniqueLocalizedGlycans.Count();
            var localizedSitesCount = uniqueLocalizedSites.Count();
            var glycanCount = uniqueGlycans.Count();
            var siteCount = uniqueSites.Count();

            var returnString = string.Format("'{0}','{1}','{2}','{3}','{4}','{5}','{6}','{7}','{8}','{9}','{10}'",
                                               uniprotID, fasta.Replace("'", "''"), inUniprot, localizedPSMCount, localizedSitesCount, localizedGlycanCount,
                                               string.Join(";", uniqueLocalizedGlycans), PSMCount, siteCount, glycanCount, string.Join(";", uniqueGlycans));

            return returnString;
        }



    }
}
