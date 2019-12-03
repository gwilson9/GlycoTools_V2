using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace _20190618_GlycoTools_V2
{
    class ProgressEventArgs : EventArgs
    {
        public string progressText { get; set; }

        public ProgressEventArgs(string progressText)
        {
            this.progressText = progressText;
        }

    }
}
