using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using CSMSL.Spectral;
using CSMSL.IO.Thermo;

namespace _20190618_GlycoTools_V2
{
    public class RTPeak : IComparable
    {
        public MZPeak Peak;
        private double _rt;
        private double _intensity;
        private double _mz;
        public double Sn;
        public int charge;
        public bool isValid;

        public RTPeak()
        {
        }

        public RTPeak(MZPeak peak, double RT)
        {
            this._intensity = peak.Intensity;
            this._mz = peak.MZ;
            this._rt = RT;
            this.Peak = peak;
            ThermoMzPeak labelPeak = ((ThermoMzPeak)peak);
            this.charge = labelPeak.Charge;
            this.Sn = labelPeak.GetSignalToNoise();
        }

        public RTPeak(double MZ, double Intensity, double RT)
        {
            _mz = MZ;
            this._intensity = Intensity;
            this._rt = RT;
        }

        public RTPeak(double MZ, double Intensity, double SN, double RT)
        {
            this._mz = MZ;
            this._intensity = Intensity;
            this.Sn = SN;
            this._rt = RT;
        }

        public double RT
        {
            get { return this._rt; }
            set { this._rt = value; }
        }

        public double Intensity
        {
            get { return this._intensity; }
            set { this._intensity = value; }
        }

        public double MZ
        {
            get { return this._mz; }
            set { this._mz = value; }
        }

        public double SN
        {
            get { return this.Sn; }
            set { this.Sn = value; }
        }

        public MZPeak MZPeak
        {
            get { return new MZPeak(this._mz, this._intensity); }
        }

        public int Compare(double other)
        {
            return RT.CompareTo(other);
        }

        public int CompareTo(RTPeak other)
        {
            return RT.CompareTo(other.RT);
        }

        public int CompareTo(Object other)
        {
            RTPeak otherPeak = (RTPeak)other;
            return RT.CompareTo(otherPeak.RT);
        }

        public bool Equals(RTPeak obj)
        {
            return obj is RTPeak && Equals((RTPeak)obj);
        }
    }
}
