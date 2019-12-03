using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace _20190618_GlycoTools_V2
{
    public class InferencePsmList : IEnumerable<InferencePSM>
    {
        private List<InferencePSM> Psms { get; set; }

        public int Count { get { return Psms.Count; } }

        public double LowestPvalue { get; private set; }

        public double CumulativePValue { get; private set; }

        public InferencePSM this[int index]
        {
            get
            {
                return Psms[index];
            }
        }

        public InferencePsmList()
        {
            Psms = new List<InferencePSM>();
            LowestPvalue = double.MaxValue;
            CumulativePValue = 1;
        }

        public InferencePsmList(InferencePSM psm)
            : this()
        {
            Add(psm);
        }

        public void Add(InferencePSM psm)
        {
            Psms.Add(psm);
            if (psm.PValue < LowestPvalue)
            {
                LowestPvalue = psm.PValue;
            }
            CumulativePValue *= psm.PValue;
        }

        public override string ToString()
        {
            return "Count = " + Count;
        }

        public IEnumerator<InferencePSM> GetEnumerator()
        {
            return Psms.GetEnumerator();
        }

        System.Collections.IEnumerator System.Collections.IEnumerable.GetEnumerator()
        {
            return Psms.GetEnumerator();
        }
    }
}
