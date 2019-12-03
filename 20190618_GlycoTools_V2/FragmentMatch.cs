using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace _20190618_GlycoTools_V2
{
    class FragmentMatch : IEquatable<FragmentMatch>
    {
        public string fragmentName { get; set; }
        public string fragmentType { get; set; }
        public int fragmentNumber { get; set; }
        public int fragmentCharge { get; set; }
        public double fragmentMZ { get; set; }
        public double fragmentSignal { get; set; }
        public string label { get; set; }

        public FragmentMatch(string fragmentName, string fragmentType, int fragmentNumber, int fragmentCharge, double fragmentMZ, double fragmentSignal)
        {
            this.fragmentName = fragmentName;
            this.fragmentType = fragmentType;
            this.fragmentNumber = fragmentNumber;
            this.fragmentCharge = fragmentCharge;
            this.fragmentMZ = fragmentMZ;
            this.fragmentSignal = fragmentSignal;
            this.label = fragmentName + "_z" + fragmentCharge;
        }

        public override bool Equals(Object obj)
        {
            var other = obj as FragmentMatch;
            if (other == null)
                return false;

            return Equals(other);
        }

        public override int GetHashCode()
        {
            return this.fragmentName.GetHashCode() +
                   this.fragmentType.GetHashCode() +
                   this.fragmentNumber.GetHashCode() +
                   this.fragmentCharge.GetHashCode() +
                   this.fragmentMZ.GetHashCode() +
                   this.fragmentSignal.GetHashCode();
        }

        public bool Equals(FragmentMatch other)
        {
            return (this.fragmentName.Equals(other.fragmentName) &&
                    this.fragmentType.Equals(other.fragmentType) &&
                    this.fragmentNumber.Equals(other.fragmentNumber) &&
                    this.fragmentCharge.Equals(other.fragmentCharge) &&
                    this.fragmentMZ.Equals(other.fragmentMZ) &&
                    this.fragmentSignal.Equals(other.fragmentSignal));
        }

    }
}
