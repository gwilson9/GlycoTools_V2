using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace _20190618_GlycoTools_V2
{
    public class SugarConstants
    {
        public readonly SugarMoiety HexNAc = new SugarMoiety("HexNAc", "Square", "C8O5NH13");
        public readonly SugarMoiety Hex = new SugarMoiety("Hex", "Circle", "C6O5H10");
        public readonly SugarMoiety Fuc = new SugarMoiety("Fuc", "Triangle", "C6O4H10");
        public readonly SugarMoiety NeuAc = new SugarMoiety("NeuAc", "Diamond (Closed)", "C11O8NH17");
        public readonly SugarMoiety NeuGc = new SugarMoiety("NeuGc", "Diamond (Open)", "C11O9NH17");
    }
}
