using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using CSMSL.Chemistry;

namespace _20190618_GlycoTools_V2
{
    public class SugarMoiety : IEquatable<SugarMoiety>
    {
        private string _name;
        private string _symbol;
        private ChemicalFormula _chemicalFormula;

        public SugarMoiety(string Name, string Symbol, string ChemicalFormula)
        {
            _name = Name;
            _symbol = Symbol;
            _chemicalFormula = new ChemicalFormula(ChemicalFormula);
        }

        public string Name
        {
            get { return _name; }
        }

        public string Symbol
        {
            get { return _symbol; }
        }

        public ChemicalFormula ChemicalFormula
        {
            get { return _chemicalFormula; }
        }

        public bool Equals(SugarMoiety other)
        {
            if (other == null)
                return false;

            if (_name == other._name && _symbol == other._symbol)
                return true;

            return false;
        }

        public override int GetHashCode()
        {
            return _name.GetHashCode() + _symbol.GetHashCode();
        }
    }
}
