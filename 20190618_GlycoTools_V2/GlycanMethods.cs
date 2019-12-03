using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using CSMSL;

namespace _20190618_GlycoTools_V2
{
    class GlycanMethods
    {
        public static double GetGlycanMass(List<SugarMoiety> AllSugars)
        {
            var returnVal = 0.0;
            foreach (var sugar in AllSugars)
            {
                returnVal += sugar.ChemicalFormula.MonoisotopicMass;
            }
            return returnVal;
        }

        public static double GetGlycanMZ(List<SugarMoiety> AllSugars, int Charge)
        {
            var returnVal = 0.0;
            foreach (var sugar in AllSugars)
            {
                returnVal += sugar.ChemicalFormula.MonoisotopicMass;
            }
            returnVal += (Charge * Constants.Hydrogen);
            returnVal /= Charge;
            return returnVal;
        }

        public static Dictionary<string, List<SugarMoiety>> GetValidGlycanStructures(Dictionary<string, List<SugarMoiety>> AllFragments)
        {
            Dictionary<string, List<SugarMoiety>> returnDict = new Dictionary<string, List<SugarMoiety>>();
            //RULES GO HERE
            //each fragment must have at least one HexNAc
            //if there is a Hex there must be at least two HexNAcs

            foreach (var fragment in AllFragments)
            {
                var tmpComponentDict = new Dictionary<string, int>();
                foreach (var entry in fragment.Value)
                {
                    if (!tmpComponentDict.ContainsKey(entry.Name))
                    {
                        tmpComponentDict.Add(entry.Name, 0);
                    }
                    tmpComponentDict[entry.Name]++;
                }
                //Sequentially apply rules, if a criterion is not met then go back to the top of the loop via continue.
                //If you make it through all the criteria then add this entry to the returnDict
                if (!tmpComponentDict.ContainsKey("HexNAc"))//each fragment must have at least one HexNAc
                {
                    continue;
                }
                if (tmpComponentDict.ContainsKey("Hex") && tmpComponentDict["HexNAc"] < 2)
                {
                    continue;
                }
                if (tmpComponentDict.ContainsKey("Fuc"))
                {
                    if (tmpComponentDict["Fuc"] > 1)
                    {
                        if (tmpComponentDict["HexNAc"] < 2)
                        {
                            continue;
                        }
                        if (tmpComponentDict.ContainsKey("Hex"))
                        {
                            if (tmpComponentDict["Hex"] < 3)
                            {
                                continue;
                            }
                        }
                    }
                }
                returnDict.Add(fragment.Key, fragment.Value);
            }
            return returnDict;
        }
    }
}
