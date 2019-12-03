using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace _20190618_GlycoTools_V2
{
    public class InferenceProteinGroup : IEnumerable<InferenceProtein>, IEquatable<InferenceProteinGroup>, IComparer<InferenceProteinGroup>, IComparable<InferenceProteinGroup>
    {
        public static int GroupNumber = 1;

        public int LongestProteinLen
        {
            get
            {
                return RepresentativeProtein == null ? 0 : RepresentativeProtein.Length;
            }
        }

        public double SequenceCoverage
        {
            get
            {
                return RepresentativeProtein.CalculateSequenceCoverage(Peptides);
            }
        }

        public double SequenceRedundacy { get { return RepresentativeProtein.CalculateSequenceRedundancy(Peptides); } }

        public string Description { get { return RepresentativeProtein.Description; } }

        public InferenceProtein RepresentativeProtein { get; private set; }

        public string Name { get; private set; }

        public bool PassesFDR = false;

        public int Count
        {
            get
            {
                return _proteins.Count;
            }
        }

        private double _pScore = double.NaN;

        public double PScore
        {
            get
            {
                // Force a p-score update if one is not given already
                if (double.IsNaN(_pScore))
                {
                    CalculatePScore();
                }
                return _pScore;
            }
        }

        private int _uniquePeptides = -1;

        public int UniquePeptides
        {
            get
            {
                if (_uniquePeptides < 0)
                {
                    _uniquePeptides = NumberofUniquePeptides();
                }
                return _uniquePeptides;
            }
        }

        private readonly List<InferenceProtein> _proteins;

        private HashSet<InferencePeptide> _peptides;
        public HashSet<InferencePeptide> Peptides
        {
            get
            {
                return _peptides;
            }
            set
            {
                _peptides = value;
                _pScore = double.NaN;
                _uniquePeptides = -1;
            }
        }

        public InferenceProtein this[int index]
        {
            get
            {
                return _proteins[index];
            }
        }

        public HashSet<string> ProteinIDs { get; private set; }

        public HashSet<string> GeneNames { get; private set; }

        public InferenceProteinGroup(InferenceProtein protein, HashSet<InferencePeptide> peptides)
        {
            Name = "PG" + GroupNumber;
            GroupNumber++;
           
            _proteins = new List<InferenceProtein>(5);
            ProteinIDs = new HashSet<string>();
            GeneNames = new HashSet<string>();
            Add(protein);
            Peptides = peptides;
            foreach (InferencePeptide pep in peptides)
            {
                pep.AddProteinGroup(this);
            }

        }

        // Add a protein to this protein group because you cannot tell the proteins apart
        public void Add(InferenceProtein prot)
        {
            // A protein group is a decoy group if any of its proteins are decoy
            if (prot.IsDecoy) _decoy = true;

            //
            if (prot.Length >= LongestProteinLen)
            {
                RepresentativeProtein = prot;
            }

            _hCode ^= prot.GetHashCode();

            // Add the protein ids
            if (!string.IsNullOrEmpty(prot.ProteinID))
                ProteinIDs.Add(prot.ProteinID);

            // Add the gene names
            if (!string.IsNullOrEmpty(prot.GeneName))
                GeneNames.Add(prot.GeneName);

            // Add the protein to the internal list
            _proteins.Add(prot);
        }

        public void UpdatePValue()
        {
            _pScore = CalculatePScore();
        }

        private double CalculatePScore()
        {
            double score = 1;
            foreach (InferencePeptide pep in Peptides)
            {
                score *= pep.PSMs.LowestPvalue;
            }
            return score;
        }

        public override string ToString()
        {
            return string.Format("{0} (p-value = {1:G3}", Name, PScore);
        }

        private int NumberofUniquePeptides()
        {
            return Peptides.Sum(pep => pep.PSMs.Count);
        }




        public string ProteinIdsString()
        {
            StringBuilder sb = new StringBuilder();
            bool inId = false;
            foreach (string proteinID in ProteinIDs)
            {
                inId = true;
                sb.Append(proteinID);
                sb.Append('|');
            }
            if (inId)
                sb.Remove(sb.Length - 1, 1);
            return sb.ToString();
        }

        public string GeneNamesString()
        {
            StringBuilder sb = new StringBuilder();
            bool inId = false;
            foreach (string genename in GeneNames)
            {
                inId = true;
                sb.Append(genename);
                sb.Append('|');
            }
            if (inId)
                sb.Remove(sb.Length - 1, 1);
            return sb.ToString();
        }

        private string GetPeptideString()
        {
            StringBuilder sb = new StringBuilder();
            return sb.ToString();
        }

        private string GetProteinStrings()
        {
            StringBuilder sb = new StringBuilder();
            foreach (InferenceProtein prot in _proteins)
            {
                sb.Append(',');
                sb.Append("\"" + prot.Description + "\"");
                sb.Append(',');
                sb.Append(prot.Sequence);
                sb.Append(',');
                sb.Append(prot.Length);
                sb.Append(',');
                double seqcov = prot.CalculateSequenceCoverage(Peptides);
                sb.Append((int)(seqcov * prot.Length) / 100);
                sb.Append(',');
                sb.Append(seqcov.ToString("g4"));
                sb.Append(',');
                sb.AppendLine(prot.CalculateSequenceRedundancy(Peptides).ToString("g4"));
            }
            return sb.ToString();
        }

        public IEnumerator<InferenceProtein> GetEnumerator()
        {
            return _proteins.GetEnumerator();
        }

        System.Collections.IEnumerator System.Collections.IEnumerable.GetEnumerator()
        {
            return _proteins.GetEnumerator();
        }

        private int _hCode = 1;

        public override int GetHashCode()
        {
            return _hCode;
        }

        public bool Equals(InferenceProteinGroup other)
        {
            if (ReferenceEquals(other, null)) return false;
            if (ReferenceEquals(this, other)) return true;
            if (Count != other.Count) return false;
            return _proteins.All(prot => other._proteins.Contains(prot));
        }

        public int Compare(InferenceProteinGroup pg1, InferenceProteinGroup pg2)
        {
            return CompareIncreasing(pg1, pg2);
        }

        public static int CompareDecreasing(InferenceProteinGroup pg1, InferenceProteinGroup pg2)
        {
            return Compare(pg1, pg2, -1);
        }

        public static int CompareIncreasing(InferenceProteinGroup pg1, InferenceProteinGroup pg2)
        {
            return Compare(pg1, pg2, 1);
        }

        public static int Compare(InferenceProteinGroup pg1, InferenceProteinGroup pg2, int direction)
        {
            direction = Math.Sign(direction);
            if (pg1 == null)
            {
                // If both x and y are null they are equal, otherwise y is greater
                return (pg2 == null) ? 0 : -direction;
            }
            // If y is null x is greater, otherwise compare their pscores
            return (pg2 == null) ? direction : direction * pg1.PScore.CompareTo(pg2.PScore);
        }

        public int CompareTo(InferenceProteinGroup other)
        {
            return Compare(this, other);
        }

        private bool _decoy;

        public bool IsDecoy
        {
            get { return _decoy; }
        }

        public double FdrScoreMetric
        {
            get { return PScore; }
        }
    }
}
