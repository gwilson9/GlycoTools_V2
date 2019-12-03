using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using CSMSL.IO.Thermo;

namespace _20190618_GlycoTools_V2
{
    class FragmentDataReturnArgs : EventArgs
    {
        public string returnType;

        public List<FragmentMatch> peptideFragments;
        public List<FragmentMatch> peptideNeutralLossFragments;
        public List<FragmentMatch> YIons;
        public List<FragmentMatch> peptideFragmentsMustIncludeGlycan;
        public List<FragmentMatch> oxoniumIons;
        public ThermoSpectrum spectrum;
        public GlycoPSM glycoPSM;

        //For batch annotation progress update
        public string timeTaken;
        public string file;

        //For peptide level data
        double seqCoverage;
        double glySeqCoverage;
        int numFrags;
        int YIonCount;
        int oxoniumIonCount;


        //For writing out each fragment
        public FragmentDataReturnArgs(List<FragmentMatch> peptideFragments, 
                                      List<FragmentMatch> peptideNeutralLossFragments,
                                      List<FragmentMatch> YIons,
                                      List<FragmentMatch> peptideFragmentsMustIncludeGlycan,
                                      List<FragmentMatch> oxoniumIons,
                                      ThermoSpectrum spectrum,
                                      GlycoPSM glycoPSM)
        {
            this.peptideFragments = peptideFragments;
            this.peptideNeutralLossFragments = peptideNeutralLossFragments;
            this.YIons = YIons;
            this.peptideFragmentsMustIncludeGlycan = peptideFragmentsMustIncludeGlycan;
            this.oxoniumIons = oxoniumIons;
            this.spectrum = spectrum;
            this.glycoPSM = glycoPSM;
            this.returnType = "UpdateSpectrumPlot";
        }

        //For displaying progress
        public FragmentDataReturnArgs(string file, string timeTaken)
        {
            this.file = file;
            this.timeTaken = timeTaken;
            this.returnType = "BatchAnnotationUpdate";
        }

        //For writing peptide level data
        public FragmentDataReturnArgs(double seqCoverage, double glySeqCoverage, int numFrags, int YIonCount, int oxoniumIonCount, GlycoPSM psm, string file)
        {
            this.seqCoverage = seqCoverage;
            this.glySeqCoverage = glySeqCoverage;
            this.numFrags = numFrags;
            this.YIonCount = YIonCount;
            this.oxoniumIonCount = oxoniumIonCount;
            this.glycoPSM = psm;
            this.file = file;
            this.returnType = "PeptideLevel";
        }
    }
}
