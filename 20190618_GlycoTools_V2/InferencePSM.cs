using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace _20190618_GlycoTools_V2
{
    public class InferencePSM
    {
        public string Sequence { get; set; }
        public double PValue { get; set; }

        public InferencePSM(string sequence, double pvalue)
        {
            Sequence = sequence;
            PValue = pvalue;
        }

    }
}
