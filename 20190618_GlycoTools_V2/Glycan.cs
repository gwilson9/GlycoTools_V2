using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace _20190618_GlycoTools_V2
{
    public class Glycan : IEquatable<Glycan>
    {
        public string _coreStructure;
        public string glycanType;
        public List<SugarMoiety> _allSugars = new List<SugarMoiety>();
        public SugarConstants _sugarConstants = new SugarConstants();
        public Dictionary<string, List<SugarMoiety>> _allCombos;

        public double mass;
        public HashSet<string> uniprotIDs = new HashSet<string>();
        public HashSet<string> glycosites = new HashSet<string>();

        public Glycan(string composition)
        {
            this._coreStructure = composition;
            ParseSugar();
            this.mass = _allSugars.Select(x => x.ChemicalFormula.MonoisotopicMass).Sum();  
            getType();            
        }

        public bool Equals(Glycan other)
        {
            return other._coreStructure.Equals(this._coreStructure);
        }

        public override int GetHashCode()
        {
            return this._coreStructure.GetHashCode();
        }

        public override string ToString()
        {
            var uniprotIDsString = string.Join(";", uniprotIDs);

            var returnString = string.Format("'{0}','{1}','{2}','{3}','{4}','{5}'",
                                              _coreStructure, glycanType, mass, glycosites.Count(),
                                              uniprotIDs.Count(), uniprotIDsString);

            return returnString;
        }

        public void ParseSugar()
        {
            // Ex HexNAc(4)Hex(5)NeuAc(2)
            var sugarsParts = _coreStructure.Split(')');

            foreach (var sugarPart in sugarsParts)
            {
                if (!string.IsNullOrEmpty(sugarPart))
                {
                    var sugar = sugarPart.Split('(')[0];
                    var count = sugarPart.Split('(')[1];

                    for (int i = 0; i < Int32.Parse(count); i++)
                    {
                        switch (sugar)
                        {
                            case "Hex":
                                _allSugars.Add(_sugarConstants.Hex);
                                break;

                            case "HexNAc":
                                _allSugars.Add(_sugarConstants.HexNAc);
                                break;

                            case "Fuc":
                                _allSugars.Add(_sugarConstants.Fuc);
                                break;

                            case "NeuAc":
                                _allSugars.Add(_sugarConstants.NeuAc);
                                break;

                            case "NeuGc":
                                _allSugars.Add(_sugarConstants.NeuGc);
                                break;

                            case "Phospho":
                                _allSugars.Add(_sugarConstants.Phospho);
                                break;

                            case "Pent":
                                _allSugars.Add(_sugarConstants.Pent);
                                break;
                        }
                    }
                }                
            }
        }

        public void AddHexNAc(int Count)
        {
            for (int i = 0; i < Count; i++)
            {
                _allSugars.Add(_sugarConstants.HexNAc);
            }
        }

        public void AddHex(int Count)
        {
            for (int i = 0; i < Count; i++)
            {
                _allSugars.Add(_sugarConstants.Hex);
            }
        }

        public void AddFuc(int Count)
        {
            for (int i = 0; i < Count; i++)
            {
                _allSugars.Add(_sugarConstants.Fuc);
            }
        }

        public void AddNeuAc(int Count)
        {
            for (int i = 0; i < Count; i++)
            {
                _allSugars.Add(_sugarConstants.NeuAc);
            }
        }

        public void AddNeuGc(int Count)
        {
            for (int i = 0; i < Count; i++)
            {
                _allSugars.Add(_sugarConstants.NeuGc);
            }
        }

        public void GenerateCombinations()
        {
            var allCombos = GetAllCombos(_allSugars);
            GetUniqueCombos(allCombos);
            var t = "";
        }

        public string CoreStructure
        {
            get { return _coreStructure; }
        }

        public Dictionary<string, List<SugarMoiety>> AllFragments
        {
            get { return _allCombos; }
        }

        public void GetUniqueCombos(List<List<SugarMoiety>> allCombos)
        {
            _allCombos = new Dictionary<string, List<SugarMoiety>>();
            foreach (var list in allCombos)
            {
                Dictionary<string, int> names = new Dictionary<string, int>();
                foreach (var item in list)
                {
                    if (!names.ContainsKey(item.Name))
                    {
                        names.Add(item.Name, 0);
                    }
                    names[item.Name]++;
                }
                var nameKeys = names.Keys.ToList();
                nameKeys.OrderBy(x => x).ToList();
                var key = "";
                foreach (var name in nameKeys)
                {
                    key += string.Format(name + "({0})", names[name]);
                }
                if (!_allCombos.ContainsKey(key))
                {
                    _allCombos.Add(key, list);
                }
            }
        }

        public static List<List<T>> GetAllCombos<T>(List<T> list)
        {
            List<List<T>> result = new List<List<T>>();
            // head
            result.Add(new List<T>());
            result.Last().Add(list[0]);
            if (list.Count == 1)
                return result;
            // tail
            List<List<T>> tailCombos = GetAllCombos(list.Skip(1).ToList());
            tailCombos.ForEach(combo =>
            {
                if (!result.Any(c => c.SequenceEqual(combo)))
                {
                    result.Add(new List<T>(combo));
                }

                combo.Add(list[0]);

                if (!result.Any(c => c.SequenceEqual(combo)))
                {
                    result.Add(new List<T>(combo));
                }
            });
            return result;
        }

        public void getType()
        {
            for (int i = 4; i < 13; i++)
            {
                string name = "HexNAc(2)Hex(" + i + ")";

                if (_coreStructure.Equals(name))
                {
                    glycanType = "High Mannose";
                    //glycanTypeNumber = 2;
                }
            }

            if (_coreStructure.Equals("HexNAc(1)") || _coreStructure.Equals("HexNAc(1)Fuc(1)") || _coreStructure.Equals("HexNAc(2)") || _coreStructure.Equals("HexNAc(2)Fuc(1)") ||
                _coreStructure.Equals("HexNAc(2)Hex(1)") || _coreStructure.Equals("HexNAc(2)Hex(1)Fuc(1)") || _coreStructure.Equals("HexNAc(2)Hex(2)") ||
                _coreStructure.Equals("HexNAc(2)Hex(2)Fuc(1)") || _coreStructure.Equals("HexNAc(2)Hex(3)") || _coreStructure.Equals("HexNAc(2)Hex(3)Fuc(1)"))
            {
                glycanType = "Paucimannose";
                //glycanTypeNumber = 1;
            }

            if (_coreStructure.Contains("NeuAc"))
            {
                glycanType = "Sialylated";
                //glycanTypeNumber = 6;
            }

            if (String.IsNullOrEmpty(glycanType) && _coreStructure.Contains("Fuc") && !_coreStructure.Contains("NeuAc"))
            {
                glycanType = "Fucosylated";
                //glycanTypeNumber = 5;
            }

            if (_coreStructure.Equals("HexNAc(2)Hex(6)Phospho(1)"))
            {
                glycanType = "Phosphomannose";
                //glycanTypeNumber = 3;
            }

            if (String.IsNullOrEmpty(glycanType))
            {
                glycanType = "Complex/Hybrid";
                //glycanTypeNumber = 4;
            }
        }
    }
}
