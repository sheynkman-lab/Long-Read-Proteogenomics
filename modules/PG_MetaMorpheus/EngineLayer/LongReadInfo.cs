using Proteomics;
using System;
using System.Collections.Generic;
using System.Text;

namespace EngineLayer
{
    public class LongReadInfo
    {
        public readonly Protein Protein;
        public readonly double CPM;

        public LongReadInfo(Protein protein, double cpm)
        {
            this.Protein = protein;
            this.CPM = cpm;
        }

        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();

            sb.Append(CPM);

            return sb.ToString();
        }
    }
}
