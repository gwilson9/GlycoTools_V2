using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace _20190618_GlycoTools_V2
{
    class SpectrumFragment
    {

        public double mz;
        public double intensity;
        public string label;

        public SpectrumFragment(double mz, double intensity, string label)
        {
            this.mz = mz;
            this.intensity = intensity;
            this.label = label;
        }


    }
}
