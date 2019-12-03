using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace _20190618_GlycoTools_V2
{
    class DataReturnArgs : EventArgs
    {
        public Dictionary<string, Dictionary<string, List<string>>> data { get; set; }

        public DataReturnArgs(Dictionary<string, Dictionary<string, List<string>>> data)
        {
            this.data = data;
        }
    }
}
