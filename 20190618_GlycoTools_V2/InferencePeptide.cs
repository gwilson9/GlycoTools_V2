using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace _20190618_GlycoTools_V2
{
    public class InferencePeptide : IEquatable<InferencePeptide>
    {
        public bool IsComplete;        
        public bool IsMapped;        
        public string Sequence;
        public string LeucineSequence;
        private readonly int _hCode;
        public static int MappedCount = 0;

        public InferencePsmList PSMs;
        public List<InferenceProteinGroup> ProteinGroups { get; private set; }
        public InferenceProteinGroup BestPG;

        public InferencePeptide(string seq)
        {
            Sequence = seq;
            LeucineSequence = seq.Replace("I", "L");
            _hCode = seq.GetHashCode();
            PSMs = new InferencePsmList();

        }

        public int NumberOfSharingProteinGroups
        {
            get
            {
                return ProteinGroups != null ? ProteinGroups.Count : 0;
            }
        }

        public bool IsShared
        {
            get
            {
                return NumberOfSharingProteinGroups > 1;
            }
        }

        public int Length
        {
            get
            {
                return LeucineSequence.Length;
            }
        }

        public bool Equals(InferencePeptide other)
        {
            return ReferenceEquals(this, other) || LeucineSequence.Equals(other.LeucineSequence);
        }

        public override int GetHashCode()
        {
            return _hCode;
        }

        public void MarkAsMapped()
        {
            if (IsMapped)
                return;
            MappedCount++;
            IsMapped = true;
        }

        public override string ToString()
        {
            return LeucineSequence;
        }

        internal void AddProteinGroup(InferenceProteinGroup pg)
        {
            if (ProteinGroups == null)
                ProteinGroups = new List<InferenceProteinGroup>();

            ProteinGroups.Add(pg);
        }

    }
}
