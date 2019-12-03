using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Net.Http;
using System.Net.Http.Headers;
using System.Net;
using System.Net.Http;
using System.Net.Http.Headers;
using System.Net.Http.Formatting;
using System.IO;

namespace _20190618_GlycoTools_V2
{
    public class DataObject
    {
        public string Name { get; set; }
    }

    class UniprotAPI
    {
        private string organism;

        public UniprotAPI(string organism)
        {
            this.organism = organism;
        }

        public Dictionary<string, Dictionary<string, string>> GetGlycosites()
        {
            string urlParameters = string.Format("?query=reviewed:yes+AND+organism:{0}&columns=id,feature(GLYCOSYLATION)&format=tab", organism);
            var numTries = 0;
            var client = new HttpClient();
            client.BaseAddress = new Uri("https://www.uniprot.org/uniprot/");
            client.DefaultRequestHeaders.Accept.Clear();
            client.DefaultRequestHeaders.Accept.Add(new MediaTypeWithQualityHeaderValue("application/txt"));

            HttpResponseMessage response = null;
            response = SendHttpRequest(client, urlParameters, numTries);

            if (response.IsSuccessStatusCode)
            {
                var data = response.Content.ReadAsStringAsync().Result.Split('\n');
                if (!data[0].Equals("Entry\tGlycosylation"))
                {
                    //OnUpdateProgress("Failed to get Uniprot Data");
                    return null;
                }

                // Key: UniprotID, Value: List of glycosites positions and evidence type
                var returnDict = new Dictionary<string, Dictionary<string,string>>();

                //var writer = new StreamWriter(@"C:\Users\gwilson\Desktop\Temp2\UniprotAPICall.csv");

                foreach (var line in data)
                {
                    var lineParts = line.Split('\t');
                    var uniprotID = lineParts[0];

                    // Skip unglycosylation proteins and header line
                    if (lineParts.Count() > 1 && !string.IsNullOrEmpty(lineParts[1]) && !lineParts[1].Equals("Glycosylation"))
                    {
                        var glycoSites = lineParts[1].Split(';');
                        foreach (var site in glycoSites)
                        {
                            //e.g."CARBOHYD 36 36 N-linked (GlcNAc...) asparagine. {ECO:0000255}."
                            var siteData = site.Trim(' ').Split(' ');

                            if (uniprotID.Equals("P10909"))
                            {
                                var x = 0;
                            }

                            if (siteData.Count() > 5)
                            {
                                var position = siteData[1];
                                var evidenceType = "";

                                if (site.Contains("ECO:0000255"))
                                    evidenceType = "Sequence Analysis";

                                if (site.Contains("ECO:0000250"))
                                    evidenceType = "By Similarity";

                                if (site.Contains("ECO:0000269"))
                                    evidenceType = "Reported in Literature";

                                if (returnDict.ContainsKey(uniprotID))
                                {
                                    if(!returnDict[uniprotID].ContainsKey(position))
                                        returnDict[uniprotID].Add(position, evidenceType);
                                }
                                else
                                {
                                    returnDict.Add(uniprotID, new Dictionary<string, string>());
                                    returnDict[uniprotID].Add(position, evidenceType);
                                }
                            }                                                        
                        }
                    }
                }
                return returnDict;                
            }
            //OnUpdateProgress("Failed to get Uniprot Data");
            return null;
        }

        public static HttpResponseMessage SendHttpRequest(HttpClient client, string urlParameters, int numTries)
        {
            try
            {
                return client.GetAsync(urlParameters).Result;
            }
            catch (AggregateException e)
            {
                //Console.WriteLine("Failed to get data. Trying again: " + numTries);
                return SendHttpRequest(client, urlParameters, numTries + 1);
            }
        }

        public event EventHandler<ProgressEventArgs> UpdateProgress;

        protected virtual void OnUpdateProgress(string progress)
        {
            var handler = UpdateProgress;
            if (handler != null)
            {
                handler(this, new ProgressEventArgs(progress));
            }
        }
    }
}
